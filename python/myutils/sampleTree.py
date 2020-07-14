#!/usr/bin/env python
from __future__ import print_function

import ROOT
import os
import sys
import time
import glob
import BetterConfigParser
from BranchList import BranchList
from FileLocator import FileLocator
from samplesclass import Sample
import array
import resource
import gc
from XbbConfig import XbbConfigReader, XbbConfigTools
import BetterConfigParser
from sample_parser import ParseInfo
import hashlib
import subprocess
from copytreePSI import filelist

# ------------------------------------------------------------------------------
# sample tree class
# 
# reads many .root files and chain them as a TChain and adds all the TH1D count
# histograms together. Provides an iterator which handles the updating of
# TTreeFormulas automatically when the tree number changes.
#
# create: (A)  sampleTree = SampleTree('path/to/indexfile.txt')
# create: (B)  sampleTree = SampleTree(['path/to/file1.root', 'path/to/file2.root'])
# create: (C)  sampleTree = SampleTree({'name': 'DY50to100', 'folder': ...}, config=config)
# create: (D)  sampleTree = SampleTree({'sample': <object of Sample class>, 'folder': ...}, config=config)
#
# add formula: sampleTree.addFormula('ptCut','pt>100')
#
# loop over events:
#
# for event in sampleTree: 
#     print 'pt cut:', sampleTree.evaluate('ptCut')
#     print 'pt:', event.pt
# ------------------------------------------------------------------------------
class SampleTree(object):

    def __init__(self, samples, treeName=None, limitFiles=-1, splitFilesChunkSize=-1, chunkNumber=1, countOnly=False, verbose=True, config=None, saveMemory=False, xrootdRedirector=None, fileNamesToProcess=None, fileLocator=None):

        self.new_formulas = []
        self.track_new = False

        # sequentialProcessing not much tested yet, can maybe be used as workaround for problems with ROOT version 6.12.07. (not needed in 6.10.09)
        self.sequentialProcessing = True
        self.verbose = verbose
        self.debug = 'XBBDEBUG' in os.environ
        self.debugProfiling = 'XBBPROFILING' in os.environ
        self.config = config
        self.saveMemory = saveMemory
        self.outputTreeBasketSize = None
        if self.config and self.config.has_option('Configuration', 'outputTreeBasketSize'):
            self.outputTreeBasketSize = eval(self.config.get('Configuration', 'outputTreeBasketSize'))
        self.monitorPerformance = True
        self.disableBranchesInOutput = True
        if self.config is not None:
            self.requireAllInputTrees = eval(self.config.get('Configuration', 'requireAllInputTrees')) if self.config.has_option('Configuration', 'requireAllInputTrees') else True
        else:
            self.requireAllInputTrees = False
        self.samples = samples
        self.tree = None
        self.nonEmptyTrees = []
        self.processedTrees = []
        #print('here')
        self.fileLocator = FileLocator(config=self.config, xrootdRedirector=xrootdRedirector) if fileLocator is None else fileLocator
        #print('here_again')
        self.sampleIdentifier = None
        self.numParts = -1
        self.subcut = None

        # friend trees (TODO: not fully implemented/tested)
        # {'name': 'DNN_reeval_v2', 'tree': tree, 'treeName': 'Events'}
        self.friends = []
        self.friendTreeIndex = ['run', 'event']
        # keep list of all temporary files to be deleted when processing finishes
        self.temporaryFiles = []

        # process only partial sample root file list
        self.splitFilesChunkSize = splitFilesChunkSize
        self.chunkNumber = chunkNumber

        # get list of sample root files to process
        sampleFileNamesParts = self.getSampleFileNameChunks()
        #print('sampleFileNamesParts', sampleFileNamesParts)  #[['root://t3dcachedb03.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/krgedia/VHbb/Zll/VHbbPostNano2018/sys_5may/ZZ_TuneCP5_13TeV-pythia8/tree_RunIIAutumn18NanoAOD-102X_upgr86_190416_171253_0000_10_619c0318d0c592698755987d1bd9208fea26b3e3fea0d4dd3a9a7f54.root', 'root://t3dcachedb03.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/krgedia/VHbb/Zll/VHbbPostNano2018/sys_5may/ZZ_TuneCP5_........]]
        #print('self.numParts', self.numParts)  #1
        #print('self.chunkNumber', self.chunkNumber) #1
        if self.chunkNumber > 0 and self.chunkNumber <= self.numParts:
            if len(sampleFileNamesParts) == self.numParts:
                chunkIndex = self.chunkNumber - 1  
                self.sampleFileNames = sampleFileNamesParts[chunkIndex]
                #print('self.sampleFileNames', self.sampleFileNames)  #['root://t3dcachedb03.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/krgedia/VHbb/Zll/VHbbPostNano2018/sys_5may/ZZ_TuneCP5_13TeV-pythia8/tree_RunIIAutumn18NanoAOD-102X_upgr86_190416_171253_0000_10_61....] 

                # per default add all files to the chain. By using fileNamesToProcess it is possible to only add a part of the files to the chain but at
                #  the same time retain the normalization histograms from ALL files, this is needed for some split processing modes of the samples which then
                #  can be added together (e.g. for datacards)
                self.sampleFileNamesToProcess = fileNamesToProcess if fileNamesToProcess else self.sampleFileNames
            else:
                raise Exception("InvalidNumberOfSplitParts")
        else:
            print("\x1b[31mERROR: wrong chunk number ", self.chunkNumber, "\x1b[0m")
            raise Exception("InvalidChunkNumber")
        if self.verbose:
            print ("INFO: reading part ", self.chunkNumber, " of ", self.numParts)

        self.status = 0
        if not treeName:
            if self.config and self.config.has_option('Configuration', 'treeName'):
                self.treeName = self.config.get('Configuration', 'treeName')
                #print('self.treeName', self.treeName) #Events
            else:
                # HEPPY default
                self.treeName = 'tree'
        else:
            self.treeName = treeName
        self.formulas = {}
        self.formulaDefinitions = []
        self.oldTreeNum = -1
        self.limitFiles = int(limitFiles) 
        self.timeStart = time.time()
        self.timeETA = 0
        self.eventsRead = 0
        self.outputTrees = []
        self.callbacks = {}
        self.removeBranches = []

        # e.g. for additional branches to be added
        self.newBranches = []

        # check existence of sample .txt file which contains list of .root files
        self.sampleTextFileName = ''

        # add all .root files to chain and add count histograms
        self.chainedFiles = []
        self.brokenFiles = []
        self.histograms = {}
        self.nanoTreeCounts = {}
        self.totalNanoTreeCounts = {}

        if not countOnly:
            self.tree = ROOT.TChain(self.treeName)

            # loop over all given .root files 
            for rootFileName in self.sampleFileNames:
                if self.debug:
                    print('DEBUG: next file is:', rootFileName, ", check existence")

                # check root file existence
                if self.fileLocator.exists(rootFileName, attempts=3):
                    remoteRootFileName = self.fileLocator.getRemoteFileName(rootFileName)
                    input = ROOT.TFile.Open(remoteRootFileName, 'read')
                    #print('remoteRootFileName', remoteRootFileName)

                    # check file validity
                    if input and not input.IsZombie() and input.GetNkeys() > 0 and not input.TestBit(ROOT.TFile.kRecovered):
                        if self.debug:
                            print('DEBUG: file exists and is good!')

                        # add count histograms, since they are not in the TChain
                        #print('input.GetListOfKeys()', input.GetListOfKeys())
                        for key in input.GetListOfKeys():
                            obj = key.ReadObj()
                            #print('obj.GetName()', obj.GetName())   #iterated: Runs, autoPU, Events
                            if obj.GetName() == self.treeName:
                                continue
                            histogramName = obj.GetName()
                            #print('histogramName', histogramName)   #iterated: Runs, autoPU

                            # nanoAOD: use branch of a tree instead of histogram for counting
                            if histogramName == 'Runs':
                                branchList = [x.GetName() for x in obj.GetListOfBranches()]
                                #print('branchList', branchList)   #['run', 'genEventSumw2', 'nLHEPdfSumw', 'nLHEScaleSumw', 'genEventSumw', 'genEventCount']
                                if self.debug:
                                    print ("DEBUG: nano counting tree has the following BRANCHES:", branchList)
                                for branch in branchList:
                                    if branch not in self.nanoTreeCounts:
                                        self.nanoTreeCounts[branch] = []
                                nEntries = obj.GetEntries()
                                #print('nEntries', nEntries) 1 always
                                for i in range(nEntries):
                                    obj.GetEntry(i)
                                    for branch in branchList:
                                        self.nanoTreeCounts[branch].append(getattr(obj, branch))
                                #print('self.nanoTreeCounts', self.nanoTreeCounts)  #gets append interate:  {'run': [1L], 'genEventSumw2': [132971.46875], 'nLHEPdfSumw': [0L], 'nLHEScaleSumw': [0L], 'genEventSumw': [132985.65625], 'genEventCount': [133000L]}... {'run': [1L, 1L, 1L], 'genEventSumw2': [132971.46875, 6998.5078125, 48989.41796875], 'nLHEPdfSumw': [0L, 0L, 0L], 'nLHEScaleSumw': [0L, 0L, 0L], 'genEventSumw': [132985.65625, 6999.25, 48994.6796875], 'genEventCount': [133000L, 7000L, 49000L]}

                            if histogramName in self.histograms:
                                if obj.IsA().InheritsFrom(ROOT.TTree.Class()):
                                    if self.debug:
                                        print("DEBUG: object is a tree and will be skipped:", obj.GetName())
                                else:
                                    if self.histograms[histogramName]:
                                        #print('here autoPU hists are added using hist1.Add(hist2)')
                                        self.histograms[histogramName].Add(obj)
                                    else:
                                        print ("ERROR: histogram object was None!!!")
                                        raise Exception("CountHistogramMissing")
                            else:
                                #only for first time... autoPU hist only
                                # add all TH*'s in one single histogram
                                if obj.IsA().InheritsFrom(ROOT.TH1.Class()):
                                    self.histograms[histogramName] = obj.Clone(obj.GetName())
                                    self.histograms[histogramName].SetDirectory(0)
                                else:
                                    if self.debug:
                                        print("DEBUG: omitting object ", obj, type(obj), " since it is neither TH1 or TTree!")

                        input.Close()

                        # add file to chain
                        if rootFileName in self.sampleFileNamesToProcess:
                            chainTree = '%s/%s'%(remoteRootFileName.strip(), self.treeName.strip())

                            if self.debug:
                                print ('\x1b[42mDEBUG: chaining '+chainTree,'\x1b[0m')

                            nEntriesBefore = self.tree.GetEntries()
                            statusCode = self.tree.Add(chainTree)
                            nEntriesAfter = self.tree.GetEntries()

                            if nEntriesAfter>nEntriesBefore:
                                self.nonEmptyTrees.append(rootFileName)

                            #chainTree = root://t3dcachedb03.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/krgedia/VHbb/Zll/VHbbPostNano2018/sys_5may/ZZ_TuneCP5_13TeV-pythia8/tree_RunIIAutumn18NanoAOD-102X_upgr86_190416_171253_0000_8_58d6472a6c2138455ddc891e46742bf0bf4635409bea6cbea1891808.root/Events                         
                            #print('This is where each tree is added to TChain object self.tree')
                            #print('statusCode', statusCode)   #1 always (if added)
                            if self.debug:
                                print ('\x1b[42mDEBUG: ---> %r'%statusCode,'\x1b[0m')
     
                            # check for errors in chaining the file
                            if statusCode != 1:
                                print ('ERROR: failed to chain ' + chainTree + ', returned: ' + str(statusCode), 'tree:', self.tree)
                                raise Exception("TChain method Add failure")
                            elif not self.tree:
                                print ('\x1b[31mERROR: tree died after adding %s.\x1b[0m'%rootFileName)
                            else:
                                self.treeEmpty = False
                                self.chainedFiles.append(rootFileName)
                                #print('self.chainedFiles', self.chainedFiles) #interated: ['root://t3dcachedb03.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/krgedia/VHbb/Zll/VHbbPostNano2018/sys_5may/ZZ_TuneCP5_13TeV-pythia8/tree_RunIIAutumn18NanoAOD-102X_upgr86_190416_171253_0000_10_619c0318d0c592698755987d1bd9208fea26b3e3fea0d4dd3a9a7f54.root', 'root://t3dcachedb03.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/krgedia/VHbb/Zll/VHbbPostNano2018/sys_5may/ZZ_TuneCP5_13TeV-pythia8/tree_RunIIAutumn18NanoAOD-102X_upgr86_190416_171253_0000_16_ebe7502f4b0d884060d03cadfddaf82ef6201c8a8f3c6c12c73d290d.root'............]
                                if self.limitFiles > 0 and len(self.chainedFiles) >= self.limitFiles:
                                    print ('\x1b[35mDEBUG: limit reached! no more files will be chained!!!\x1b[0m')
                                    break
                    else:
                        print ('\x1b[31mERROR: file is damaged: %s\x1b[0m'%rootFileName)
                        if input:
                            print ('DEBUG: Zombie:', input.IsZombie(), '#keys:', input.GetNkeys(), 'recovered:', input.TestBit(ROOT.TFile.kRecovered))
                        self.brokenFiles.append(rootFileName)
                else:
                    print ('\x1b[31mERROR: file is missing: %s\x1b[0m'%rootFileName)

            if self.verbose or self.debug:
                print ('INFO: # files chained: %d'%len(self.chainedFiles))
                if len(self.brokenFiles) > 0:
                    print ('INFO: # files broken : %d'%len(self.brokenFiles))
            
            if len(self.chainedFiles) < 1:
                self.tree = None

            if self.tree:
                self.tree.SetCacheSize(50*1024*1024)

            #print('self.nanoTreeCounts', self.nanoTreeCounts)

            # merge nano counting trees
            if self.nanoTreeCounts:
                # TODO: per run if possible, sum LHE weights if present

                # sum the contributions from the subtrees
                self.totalNanoTreeCounts = {key: sum(values) for key,values in self.nanoTreeCounts.iteritems() if len(values) > 0 and type(values[0]) in [int, float, long]}
                #print('self.totalNanoTreeCounts', self.totalNanoTreeCounts)  #self.totalNanoTreeCounts {'run': 19L, 'genEventSumw2': 1978563.8623046875, 'nLHEPdfSumw': 0L, 'nLHEScaleSumw': 0L, 'genEventSumw': 1978776.7690429688, 'genEventCount': 1979000L}

                # print summary table
                countBranches = self.totalNanoTreeCounts.keys()
                #print('countBranches', countBranches) #['run', 'genEventSumw2', 'nLHEPdfSumw', 'nLHEScaleSumw', 'genEventSumw', 'genEventCount']
                depth = None
                for key,values in self.nanoTreeCounts.iteritems():
                    if values and len(values)>1 and type(values[0]) in [int, float, long]:
                        depth = len(values)
                        break
                #print('depth', depth) #19 total trees chained       
                print("-"*160)
                print("tree".ljust(25), ''.join([countBranch.ljust(25) for countBranch in countBranches]))
                #print("tree".ljust(25), ''.join([countBranch.ljust(25) for countBranch in countBranches]))
                if depth:
                    for treeNum in range(depth):
                        try:
                            print(("%d"%(treeNum+1)).ljust(25),''.join([('%r'%self.nanoTreeCounts[countBranch][treeNum]).ljust(25) for countBranch in countBranches]))
                        except Exception as e:
                            print(e)

                print("\x1b[34m","sum".ljust(24), ''.join([('%r'%self.totalNanoTreeCounts[countBranch]).ljust(25) for countBranch in countBranches]),"\x1b[0m")
                print("-"*160)

                #print('self.histograms',self.histograms)  #{'autoPU': <ROOT.TH1D object ("autoPU") at 0x67ac770>}

                # fill summed tree (create new tree) 
                #keys are different branches and values are corresponding summed values of all tree printed just above. You fill them in a new TTree 'Runs' in dictonary self.histogram with same bracnhes and their summed values.
                self.histograms['Runs'] = ROOT.TTree('Runs', 'count histograms for nano')
                nanoTreeCountBuffers = {}
                for key, value in self.totalNanoTreeCounts.iteritems():
                    #print('value', value)
                    if type(value) == int:
                        # 64 bit signed int 
                        typeCode = 'L'
                    elif type(value) == long:
                        typeCode = 'L'
                    elif type(value) == float:
                        typeCode = 'f'
                    nanoTreeCountBuffers[key] = array.array(typeCode, [value])
                    self.histograms['Runs'].Branch(key, nanoTreeCountBuffers[key], '{name}/{typeCode}'.format(name=key, typeCode=typeCode))
                    
                self.histograms['Runs'].Fill()
                #print('self.histograms',self.histograms)  #{'Runs': <ROOT.TTree object ("Runs") at 0x6640b90>, 'autoPU': <ROOT.TH1D object ("autoPU") at 0x67ac770>}

    def __del__(self):
        self.delete()

    def delete(self):
        self.callbacks = None
        # close possible left open files referencing the TChain and delete output trees
        try:
            if self.tree:
                self.tree.Reset()
        except:
            pass
        self.fileLocator = None
        self.config = None
        if hasattr(self, "outputTrees"):
            for outputTree in self.outputTrees:
                del outputTree['file']
        try:
            if hasattr(self, "formulas"):
                for formulaName, formula in self.formulas.iteritems():
                    if formula:
                        del formula
                        formula = None
        except Exception as e:
            print("EXCEPTION: (destructor) ", e)
        try:
            if hasattr(self, "outputTrees"):
                for outputTree in self.outputTrees:
                    if outputTree['tree']:
                        del outputTree['tree']
                        outputTree['tree'] = None
        except Exception as e:
            print("EXCEPTION: (destructor)", e)
        try:
            if self.tree:
                del self.tree
                self.tree = None
        except Exception as e:
            print("EXCEPTION: (destructor)", e)

    def getFriendTreeFileNames(self, directory):
        # get list of files in this chain and make list of friend file names
        sampleFileNames = self.sampleFileNamesToProcess

        sampleFileNamesFriend = []
        for fullFileName in sampleFileNames:

            fileNameParts = fullFileName.split('/')
            fileName = fileNameParts[-1]
            parentDirectory = fileNameParts[-2] if len(fileNameParts) > 1 else './'

            friendFileName = directory
            if '{fileName}' in friendFileName:
                friendFileName = friendFileName.replace('{fileName}', fileName)
            if '{sampleIdentifier}' in friendFileName:
                friendFileName = friendFileName.replace('{sampleIdentifier}', self.sampleIdentifier)
            if '{parentDirectory}'  in friendFileName:
                friendFileName = friendFileName.replace('{parentDirectory}', parentDirectory)

            sampleFileNamesFriend.append(friendFileName)
        return sampleFileNamesFriend

    def addFriend(self, directory, treeName='Events_f', index=None, name='<friend>', alias=None, copy=False):
        if not index:
            index = ['run', 'event']

        friend = {
                'tree': ROOT.TChain(treeName),
                'treeName': treeName,
                'name': name
                }

        sampleFileNamesFriend = self.getFriendTreeFileNames(directory)
        for fileName in sampleFileNamesFriend:
            # set copy=True is useful in case the number of concurrent connections to storage is a bottleneck
            if copy:
                scratchDir = self.config.get('Directories','scratch') + '/friendTreeCache/'
                try:
                    os.makedirs(scratchDir)
                except:
                    pass
                fileNameHash = hashlib.sha224(fileName).hexdigest()
                tempFileName = scratchDir + '/' + fileNameHash + '_as_' + (alias if alias is not None else 'friendTree') + '.root' 
                if os.path.isfile(tempFileName):
                    self.fileLocator.rm(tempFileName)
                self.fileLocator.cp(fileName, tempFileName, force=True)
                friend['tree'].Add(tempFileName + '/' + treeName)
                self.temporaryFiles.append(tempFileName)
                if self.debug:
                    print("DEBUG: copied to temporary file:", tempFileName)
            else:
                friend['tree'].Add(fileName + '/' + treeName)

        #friend['tree'].BuildIndex()
        if alias is not None:
            self.tree.AddFriend(friend['tree'], alias)
        else:
            self.tree.AddFriend(friend['tree'])

        self.friends.append(friend)
        if self.debug:
            print("INFO: added friend tree", treeName)
            print("DEBUG: \x1b[42m-> ", sampleFileNamesFriend, "\x1b[0m")

    # ------------------------------------------------------------------------------
    # check if all files have been chained
    # ------------------------------------------------------------------------------
    def isCompleteTree(self):
        return len(self.chainedFiles) == len(self.sampleFileNamesToProcess) and len(self.sampleFileNamesToProcess) > 0

    # ------------------------------------------------------------------------------
    # return full list of sample root files 
    # ------------------------------------------------------------------------------
    def getAllSampleFileNames(self, samples=None): 

        # if no argument is given, use default files
        if not samples:
            samples = self.samples

        # given argument is list -> this is already the list of root file names
        if type(samples) == list:
            sampleFileNames = samples
        # given argument is name and folder -> glob
        elif type(samples) == dict:
            isData = False
            if 'sample' in samples:
                sampleName = samples['sample'].identifier
                if samples['sample'].subsample:
                    self.subcut = samples['sample'].subcut
                isData = isData or samples['sample'].isData()
            else:
                sampleName = samples['name']
            self.sampleIdentifier = sampleName
            sampleFolder = samples['folder']

            # this requires all the trees from the original /store input to be still present (not needed for MC) 
            if self.requireAllInputTrees or isData:
                # get list of original filenames at /store/...
                samplefiles    = self.config.get('Directories', 'samplefiles')
                sampleFileList = filelist(samplefiles, self.sampleIdentifier)

                # get list of filenames 
                sampleFileNames = ["{path}/{sample}/{fileName}".format(path=sampleFolder, sample=self.sampleIdentifier, fileName=self.fileLocator.getFilenameAfterPrep(fileName)) for fileName in sampleFileList]  

            else:
                samplesMask = self.fileLocator.getLocalFileName(sampleFolder) + '/' + sampleName + '/*.root'
                redirector = self.fileLocator.getRedirector(sampleFolder)
                if self.verbose:
                    print ("INFO: use ", samplesMask)

                sampleFileNames = self.fileLocator.glob_with_fallback(samplesMask)
                if sampleFileNames is None or len(sampleFileNames) < 1:
                    print("\x1b[31mERROR: no tree files found for this sample in",sampleFolder + "/\x1b[35m" + sampleName + "\x1b[31m !\x1b[0m")
                    raise Exception("FilesMissing")

                sampleFileNames = [self.fileLocator.addRedirector(redirector, x) for x in sampleFileNames]

            if self.verbose:
                print ("INFO: found ", len(sampleFileNames), " files.")
                if self.debug:
                    for sampleFileName in sampleFileNames:
                        print("DEBUG: >", sampleFileName)

        # given argument is a single file name -> read this .txt file 
        else:
            sampleTextFileName = samples
            
            if os.path.isfile(sampleTextFileName):
                self.sampleTextFileName = sampleTextFileName
                if self.verbose:
                    print('open samples .txt file: %s' % self.sampleTextFileName)
            else:
                print("\x1b[31mERROR: file not found: %s \x1b[0m" % sampleTextFileName)
                return

            with open(self.sampleTextFileName, 'r') as sampleTextFile:
                sampleFileNames = sampleTextFile.readlines()
        return sampleFileNames

    # ------------------------------------------------------------------------------
    # return lists of sample root files, split into chunks with certain size  
    # ------------------------------------------------------------------------------
    def getSampleFileNameChunks(self):
        sampleFileNames = self.getAllSampleFileNames()
        if self.splitFilesChunkSize > 0 and len(sampleFileNames) > self.splitFilesChunkSize:
            sampleFileNamesParts = [sampleFileNames[i:i + self.splitFilesChunkSize] for i in xrange(0, len(sampleFileNames), self.splitFilesChunkSize)]
        else:
            sampleFileNamesParts = [sampleFileNames]
        self.numParts = len(sampleFileNamesParts)
        #print('self.numParts = len(sampleFileNamesParts)', self.numParts)   #1
        #print('len(sampleFileNames)', len(sampleFileNames))  #19
        return sampleFileNamesParts

    # ------------------------------------------------------------------------------
    # return lists of sample root files for a single chunk  
    # ------------------------------------------------------------------------------
    def getSampleFileNameChunk(self, chunkNumber):
        chunks = self.getSampleFileNameChunks()
        if chunkNumber > 0 and chunkNumber <= len(chunks):
            return chunks[chunkNumber-1]
        else:
            print("\x1b[31mERROR: invalid chunk number {n} \x1b[0m".format(n=chunkNumber))

    def getNumberOfParts(self):
        return self.numParts

    # ------------------------------------------------------------------------------
    # add a TTreeFormula connected to the TChain
    # ------------------------------------------------------------------------------
    def addFormula(self, formulaName, formula=None):
        if formulaName in self.formulas:
            if self.debug:
                print('formulaName in self.formulas')

            if self.formulas[formulaName]:
                if self.debug:
                    print(self.formulas[formulaName])

                return

        if self.debug:
            print('Adding:', formulaName)

        if self.track_new:
            print('And Tracking')
            self.new_formulas.append(formulaName)

        if formula is None:
            formula = formulaName

        # there might be an undocumented limit on the length of cutstrings in ROOT...
        if len(formula) > 1023:
            print("\x1b[41m\x1b[97m------------------------------------------------------------------------------")
            print(" WARNING !!! ROOT.TTreeFormula of length %d, this might cause problems !!"%len(formula))
            print(" reduce length of formulas if problems occur, e.g. by passing lists of cut formulas!")
            print("\x1b[32m",formula,"\x1b[97m")
            print("------------------------------------------------------------------------------\x1b[0m")

        self.formulaDefinitions.append({'name': formulaName, 'formula': formula})
        self.formulas[formulaName] = ROOT.TTreeFormula(formulaName, formula, self.tree)
        #print('self.formulas', self.formulas)
        if self.formulas[formulaName].GetNdim() == 0:
            print("DEBUG: formula is:", formula)
            print("\x1b[31mERROR: adding the tree formula failed! Check branches of input tree and loaded namespaces.\x1b[0m")
            raise Exception("SampleTreeAddTTreeFormulaFailed")        
        #self.formulas for weightF formula: ... {'genWeight*puWeight*((Vtype==1&&HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) + (Vtype==0)*(4.767/41.298*HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 + 36.531/41.298*HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8))*muonSF_Iso[0]*muonSF_Id[0]*electronSF_IdIso[0]*electronSF_trigger[0]*bTagWeightDeepCSV*EWKw[0]*weightLOtoNLO*FitCorr[0]': <ROOT.TTreeFormula object ("genWeight*puWeight*((Vtype==1&&HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) + (Vtype==0)*(4.767/41.298*HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 + 36.531/41.298*HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8))*muonSF_Iso[0]*muonSF_Id[0]*electronSF_IdIso[0]*electronSF_trigger[0]*bTagWeightDeepCSV*EWKw[0]*weightLOtoNLO*FitCorr[0]") at 0x67a6f60>}

    # ------------------------------------------------------------------------------
    # return list of formulas
    # ------------------------------------------------------------------------------
    def getFormulas(self):
        return self.formulas

    # ------------------------------------------------------------------------------
    # add a new branch
    # ------------------------------------------------------------------------------
    def addOutputBranch(self, branchName, formula, branchType='f', length=1, arguments=None, leaflist=None, arrayStyle=False):
        # this is needed to overwrite the branch if it already exists!
        self.addBranchToBlacklist(branchName)

        # function
        if callable(formula):
            newBranch = {'name': branchName, 'function': formula, 'type': branchType, 'length': length}
            if arguments:
                newBranch['arguments'] = arguments
        # string which contains a TTreeFormula expression
        else:
            formulaName = 'alias:' + branchName
            self.addFormula(formulaName, formula)
            newBranch = {'name': branchName, 'formula': formulaName, 'type': branchType, 'length': length}
        if leaflist:
            newBranch['leaflist'] = leaflist
        if arrayStyle:
            newBranch['arrayStyle'] = True
        self.newBranches.append(newBranch)

    # ------------------------------------------------------------------------------
    # pass a list of dictionaries of branches to add
    # TODO: avoid detour via addOutputBranch and set dictionary directly
    # ------------------------------------------------------------------------------
    def addOutputBranches(self, branchDictList):
        for branchDict in branchDictList:
            self.addOutputBranch(
                branchName=branchDict['name'],
                formula=branchDict['formula'],
                branchType=branchDict['type'] if 'type' in branchDict else 'f',
                length=branchDict['length'] if 'length' in branchDict else 1,
                arguments=branchDict['arguments'] if 'arguments' in branchDict else None,
                leaflist=branchDict['leaflist'] if 'leaflist' in branchDict else None,
                arrayStyle=branchDict['arrayStyle'] if 'arrayStyle' in branchDict else False,
            )

    # ------------------------------------------------------------------------------
    # implement iterator for TChain, with updating TTreeFormula objects on tree
    # switching and show performance statistics during loop
    # ------------------------------------------------------------------------------
    def next(self):
        self.treeIterator.next()
        self.eventsRead += 1
        if self.debug and (self.eventsRead % 1000 == 0 or self.eventsRead == 10 or self.eventsRead == 100):
            throughput = 1.0*self.eventsRead / (1.0*(time.time() - self.timeStart) + 0.001)
            print('DEBUG: %d events read, %1.2f /s'%(self.eventsRead, throughput))
        treeNum = self.tree.GetTreeNumber()
        # TTreeFormulas have to be updated when the tree number changes in a TChain
        if treeNum != self.oldTreeNum:

            # update ETA estimates
            if treeNum == 0:
                self.timeStart = time.time()
                perfStats = '?'
            else:
                fraction = 1.0*treeNum/len(self.chainedFiles)
                passedTime = time.time() - self.timeStart
                self.timeETA = (1.0-fraction)/fraction * passedTime if fraction > 0 else 0
                perfStats = 'INPUT: {erps}/s, OUTPUT: {ewps}/s '.format(erps=self.eventsRead / passedTime if passedTime>0 else 0, ewps=sum([x['passed'] for x in self.outputTrees]) / passedTime if passedTime>0 else 0)

            # output status
            if self.verbose or self.debug:
                percentage = 100.0*treeNum/len(self.chainedFiles)
                if treeNum == 0:
                    print ('INFO: time ', time.ctime())
                if self.debug:
                    perfStats = perfStats + ' max mem used = %d'%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
                print ('INFO: switching trees --> %d (=%1.1f %%, ETA: %s min, %s)'%(treeNum, percentage, self.getETA(), perfStats))

                # in case files of the chain get inaccessible in the time between making the chain (file exsitence will be checked and error raised if file does not exist!) 
                # and the time the file is read (ROOT will silently ignore it ad there will be just Xrootd error message printed)
                # add processed numbers to array and check if all trees have been processed in the end
                self.processedTrees.append(treeNum)

                if self.debugProfiling:
                    self.tree.PrintCacheStats()
                sys.stdout.flush()
            self.oldTreeNum = treeNum
            # update TTreeFormula's
            for formulaName, treeFormula in self.formulas.iteritems():
                treeFormula.UpdateFormulaLeaves()
        return self.tree

    def __iter__(self):
        self.treeIterator = self.tree.__iter__()
        return self

    # wrapper to evaluate a formula, which has been added to the formula dictionary
    # vector valued formulas are not supported
    def evaluate(self, formulaName):  #https://root-forum.cern.ch/t/ttreeformula-evalinstance-return-0-0/16366/11
        if formulaName in self.formulas:
            if self.debug:
                if not self.formulas[formulaName]:
                    print('Evaluate to fail:', formulaName)

            if self.formulas[formulaName].GetNdata() > 0:
                return self.formulas[formulaName].EvalInstance()
            else:
                return 0
        else:
            existingFormulas = [x for x,y in self.formulas.iteritems()]
            print ("existing formulas are: ", existingFormulas)
            raise Exception("SampleTree::evaluate: formula '%s' not found!"%formulaName)

    # evaluates a vector formula and fills an array, returns the number of dimensions of the formula
    def evaluateArray(self, formulaName, destinationArray):
        if formulaName in self.formulas:
            nData = self.formulas[formulaName].GetNdata()
            for i in range(nData):
                destinationArray[i] = self.formulas[formulaName].EvalInstance(i)
            return nData
        else:
            existingFormulas = [x for x, y in self.formulas.iteritems()]
            print("existing formulas are: ", existingFormulas)
            raise Exception("SampleTree::evaluate: formula '%s' not found!" % formulaName)

    # return string of ETA in minutes
    def getETA(self):
        return '%1.1f'%(self.timeETA/60.0) if self.timeETA > 0 else '?'

    def GetListOfBranches(self):
        return self.tree.GetListOfBranches()

    # ------------------------------------------------------------------------------
    # handle 'tree-typed' cuts, passed as dictionary:
    # e.g. cut = {'OR': [{'AND': ["pt>20","eta<3"]}, "data==1"]}
    # short-circuit evaluation is handled by builtins any() and all()
    # ------------------------------------------------------------------------------
    def addCutDictRecursive(self, cutDict):
        if type(cutDict) == str:
            if cutDict not in self.formulas:
                self.addFormula(cutDict, cutDict)
        elif 'OR' in cutDict and 'AND' in cutDict:
            raise Exception("BadTreeTypeCutDict")
        elif 'OR' in cutDict:
            for subDict in cutDict['OR']:
                self.addCutDictRecursive(subDict)
        elif 'AND' in cutDict:
            for subDict in cutDict['AND']:
                self.addCutDictRecursive(subDict)
        else:
            raise Exception("BadTreeTypeCutDict")

    def evaluateCutDictRecursive(self, cutDict):
        if type(cutDict) == str:
            if self.formulaResults[cutDict] is None:
                print ("FORMULA:", cutDict)
                raise Exception("UnevaluatedFormula!!")
            return self.formulaResults[cutDict]
        elif 'OR' in cutDict and 'AND' in cutDict:
            raise Exception("BadTreeTypeCutDict")
        elif 'OR' in cutDict:
            return any([self.evaluateCutDictRecursive(subDict) for subDict in cutDict['OR']])
        elif 'AND' in cutDict:
            return all([self.evaluateCutDictRecursive(subDict) for subDict in cutDict['AND']])
        else:
            raise Exception("BadTreeTypeCutDict")
    
    # set callback function, which MUST return a boolean. To continue processing this event, the function must return True. False means skip this event!
    def setCallback(self, category, fcn):
        if category in self.callbacks:
            print("WARNING: callback function for ", category, " is overwritten!")
        self.callbacks[category] = [fcn]

    # add callback function, which MUST return a boolean. To continue processing this event, the function must return True. False means skip this event!
    def addCallback(self, category, fcn):
        if category not in self.callbacks:
             self.callbacks[category] = []
        self.callbacks[category].append(fcn)

    # ------------------------------------------------------------------------------
    # add output tree to be written during the process() function
    # ------------------------------------------------------------------------------
    def addOutputTree(self, outputFileName, cut, hash='', branches=None, callbacks=None, cutSequenceMode='AND', name='', friend=False):

        # write events which satisfy either ONE of the conditions given in the list or ALL
        if cutSequenceMode not in ['AND', 'OR', 'TREE']:
            raise Exception("InvalidCutSequenceMode")

        if len([x for x in self.outputTrees if x['fileName'] == outputFileName])>0:
            print("WARNING: skipping duplicate file ", outputFileName, "!")
            return False

        outputTree = {
            'tree': None, # will create this tree later, after it is known which branches will be enabled 
            'name': name,
            'fileName': outputFileName,
            'file': None, 
            'cut': cut,
            'cutSequence': [],
            'cutSequenceMode': cutSequenceMode,
            'hash': hash,
            'branches': branches,
            'callbacks': callbacks,
            'passed': 0,
        }

        if friend:
            outputTree['friend'] = {'treeName': self.treeName + '_f'}

        self.outputTrees.append(outputTree)

    def disableOutput(self):
        print("INFO: output disabled - no files will be written!")
        self.outputTrees = []

    def getNumberOfOutputTrees(self):
        return len(self.outputTrees)

    # these branches are ALWAYS removed (e.g. because they will be recomputed), even when they are in the 'keep_branches' list
    def addBranchToBlacklist(self, branchName):
        if branchName != '*':
            self.removeBranches.append(branchName)
        else:
            print("WARNING: can't add branch '*' to blacklist => igonre it!")

    # wrapper to enable/disable branches in the TChain
    def SetBranchStatus(self, branchName, branchStatus):
        listOfExistingBranches = self.GetListOfBranches()
        if listOfExistingBranches.FindObject(branchName) or '*' in branchName:
            self.tree.SetBranchStatus(branchName, branchStatus)

    # enables ONLY the given branches (* wildcards supported) and checks existence before enabling them to avoid warning messages during tree iteration
    def enableBranches(self, listOfBranchesToKeep):
        listOfExistingBranches = self.GetListOfBranches()   
        #print('self.GetListOfBranches()', len(self.GetListOfBranches())) #1762
        #print('self.tree.GetListOfBranches()', len(self.tree.GetListOfBranches())) #1762
        self.tree.SetBranchStatus("*", 0)
        enabledBranches = []
        for branchName in listOfBranchesToKeep:
            if listOfExistingBranches.FindObject(branchName) or '*' in branchName:
                self.tree.SetBranchStatus(branchName, 1)
                enabledBranches.append(branchName)
        print("INFO: reduced number of enabled branches from", len(listOfExistingBranches), " to", len(enabledBranches), " (branches with wildcards may not be correctly counted)")
        if self.verbose:
            print ("INFO: branches:", BranchList(enabledBranches).getShortRepresentation())

    # ------------------------------------------------------------------------------
    # run callbacks registered for output trees 
    # ------------------------------------------------------------------------------
    def _runCallbacks(self, name):
        if name in self.callbacks:
            for fcn in self.callbacks[name]:
                fcn()

        for outputTree in self.outputTrees:
            if outputTree['callbacks'] and name in outputTree['callbacks']:
                try:
                    outputTree['callbacks'][name]()
                except Exception as e:
                    print("\x1b[31mWARNING: exception during '", name, "' callback:", e, "\x1b[0m")
                    raise

    # ------------------------------------------------------------------------------
    # if self.sequentialProcessing is enabled, the input tree is traversed
    # multiple times which can lead to a slowdown for larger trees
    # ------------------------------------------------------------------------------
    def process(self):
        try: 
            if self.sequentialProcessing:
                print("\x1b[31mINFO: sequential processing of trees is enabled!\x1b[0m")

                # backup full list of trees
                allTrees = list(self.outputTrees)

                # process output files sequentially
                for singleTree in allTrees:
                    self.outputTrees = [singleTree]
                    self.process_do()

                # restore
                self.outputTrees = allTrees
            else:
                self.process_do()

        # handling of exceptions during the event loop
        except Exception as e:
            print("\x1b[41m\x1b[37mERROR: exception during event loop\x1b[0m")
            print("->", e)
            try:
                print("DEBUG: sampleTree:")
                print("DEBUG:  - GetEntries() =", self.tree.GetEntries())
                print("DEBUG:  - GetReadEntry() =", self.tree.GetReadEntry())
                for x in ['run', 'event', 'evt']:
                    if hasattr(self.tree, x):
                        print("DEBUG:  - %s ="%x, getattr(self.tree, x))
            except Exception as e2:
                print("ERROR: an exception occurred while handling an exception:", e2)
            raise


    def clear_new_formulas(self):
        for new_formula in self.new_formulas:
            if self.debug:
                print('Deleting:', new_formula)

            formula = self.formulas[new_formula]
            print(formula)
            del formula
            self.formulas.pop(new_formula)

        self.new_formulas = []

        self.track_new = True

    # ------------------------------------------------------------------------------
    # loop over all entries in the TChain and copy events to output trees, if the
    # cuts are fulfilled.
    # ------------------------------------------------------------------------------
    def process_do(self):
        self.clear_new_formulas()

        if self.debug:
            rsrc = resource.RLIMIT_DATA
            # restrict memory
            # resource.setrlimit(rsrc, (2.0*1024*1024*1024, 6*1024*1024*1024))
            soft, hard = resource.getrlimit(rsrc)
            print('DEBUG: mem limits soft/hard:', soft, hard)
            rsrc = resource.RLIMIT_AS
            # restrict memory
            # resource.setrlimit(rsrc, (2.0*1024*1024*1024, 6*1024*1024*1024))
            soft, hard = resource.getrlimit(rsrc)
            print('DEBUG: AS limits soft/hard:', soft, hard)
            rsrc = resource.RLIMIT_STACK
            soft, hard = resource.getrlimit(rsrc)
            print('DEBUG: stack limits soft/hard:', soft, hard)
            print('DEBUG: max mem used:', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        if self.verbose:
            print ('OUTPUT TREES:')
            for outputTree in self.outputTrees:
                cutString = "%r"%outputTree['cut']
                if len(cutString) > 50:
                    cutString = cutString[0:50] + '...(%s more chars)'%(len(cutString)-50)
                print (' > ', outputTree['fileName'], ' <== ', outputTree['name'] if 'name' in outputTree else outputTree['hash'], ' cut: ', cutString)
            print ('FORMULAS:')
            for formulaName, formula in self.formulas.iteritems():
                print (' > \x1b[35m', formulaName, '\x1b[0m ==> ', formula)

        # find common set of branches which needs to be enabled for cuts and desired variables in all of the output trees
        listOfBranchesToKeep = []
        for outputTree in self.outputTrees:
            if 'branches' in outputTree and outputTree['branches']:
                listOfBranchesToKeep += outputTree['branches']
            if 'cut' in outputTree and outputTree['cut']:
                listOfBranchesToKeep += BranchList(outputTree['cut']).getListOfBranches()
            for formula in self.formulaDefinitions:
                listOfBranchesToKeep += BranchList(formula['formula']).getListOfBranches()

        # keep the branches stated in config, (unless they will be recomputed)
        if self.config:
            listOfBranchesToKeep += eval(self.config.get('Branches', 'keep_branches'))
        listOfBranchesToKeep = list(set(listOfBranchesToKeep))

        # disable the branches in the input if there is no output tree which wants to have all branches
        if '*' not in listOfBranchesToKeep and len(listOfBranchesToKeep) > 0:
            self.enableBranches(listOfBranchesToKeep)
        else:
            if len(self.removeBranches) < 1:
                print("INFO: keep all branches")
            else:
                print("INFO: keep all branches but the following:")
                print("INFO:", ", ".join(self.removeBranches))

        # now disable all branches, which will be e.g. recomputed
        for branchName in self.removeBranches:
            self.SetBranchStatus(branchName, 0)
        
        # callbacks before loop
        self._runCallbacks('prepareOutput')

        # initialize the output trees, this has to be called after the calls to SetBranchStatus
        for outputTree in self.outputTrees:
            outputTree['file'] = ROOT.TFile.Open(outputTree['fileName'], 'recreate')
            if not outputTree['file'] or outputTree['file'].IsZombie():
                print ("\x1b[31mERROR: output file broken\x1b[0m")
                raise Exception("OutputFileBroken")
        
            # copy count histograms to output files
            outputTree['histograms'] = {}
            for histogramName, histogram in self.histograms.iteritems():
                    outputTree['histograms'][histogramName] = histogram.Clone(histogram.GetName())
                    outputTree['histograms'][histogramName].SetDirectory(outputTree['file'])

            # clone tree structure, but don't copy any entries
            outputTree['file'].cd()
            if 'friend' in outputTree:
                outputTree['tree'] = ROOT.TTree(outputTree['friend']['treeName'], 'friend tree' ) 
                outputTree['friend']['addr_run'] = array.array('i', [-1])
                outputTree['friend']['addr_event'] = array.array('l', [-1])
                outputTree['tree'].Branch('run', outputTree['friend']['addr_run'], 'run/i')
                outputTree['tree'].Branch('event', outputTree['friend']['addr_event'], 'event/l')
                print("INFO: create friend tree!")
            else:
                outputTree['tree'] = self.tree.CloneTree(0)
            # can be used to reduce memory consumption
            if self.outputTreeBasketSize:
                outputTree['tree'].SetBasketSize("*", self.outputTreeBasketSize)
            if not outputTree['tree']:
                print ("\x1b[31mWARNING: output tree broken. try to recover!\x1b[0m")
                # if input tree has 0 entries, don't copy 0 entries to the output tree, but ALL of them instead! (sic!)
                # (this is done by omitting the argument to CloneTree)
                outputTree['tree'] = self.tree.CloneTree()
                if not outputTree['tree']:
                    print ("\x1b[31mERROR: output tree broken, input tree: ", self.tree, " \x1b[0m")
                else:
                    print ("\x1b[32mINFO: recovered\x1b[0m")
            outputTree['tree'].SetDirectory(outputTree['file'])

        # add CUT formulas, this has to be called after the calls to SetBranchStatus
        for outputTree in self.outputTrees:
            if outputTree['cutSequenceMode'] == 'TREE' and type(outputTree['cut']) == dict:
                outputTree['cutSequence'] = outputTree['cut']
                # now recursively parse the cut-tree and add all contained cut formulas
                self.addCutDictRecursive(outputTree['cut'])
            elif type(outputTree['cut']) == dict:
                print ("HINT: use cutSequenceMode='TREE' to pass dictionaries!")
                raise Exception("InvalidCutSequenceMode")
            else:
                # cut passed as string or list of strings
                cutList = outputTree['cut'] if type(outputTree['cut']) == list else [outputTree['cut']]
                for i, cutString in enumerate(cutList):
                    formulaName = cutString.replace(' ', '')
                    if formulaName not in self.formulas:
                        self.addFormula(formulaName, cutString)
                    outputTree['cutSequence'].append(formulaName)

        # prepare memory for new branches to be written
        pyTypes = {'O': 'i'}
        for outputTree in self.outputTrees:
            outputTree['newBranchArrays'] = {}
            outputTree['newBranches'] = {}
            for branch in self.newBranches:
                # convert ROOT type convention to python array type convetion if necessary
                pyType = pyTypes[branch['type']] if branch['type'] in pyTypes else branch['type'] 
                outputTree['newBranchArrays'][branch['name']] = array.array(pyType, [0] * branch['length'])
                if 'leaflist' in branch:
                    leafList = branch['leaflist']
                else:
                    # TODO: branch type conversion could be done with dicts instead of upper to use unsigned ints...
                    leafList = '{name}{length}/{type}'.format(name=branch['name'], length='[%d]'%branch['length'] if branch['length'] > 1 else '', type=branch['type'].upper())
                outputTree['newBranches'][branch['name']] = outputTree['tree'].Branch(branch['name'], outputTree['newBranchArrays'][branch['name']], leafList)
        if len(self.newBranches) > 0:
            print("ADD NEW BRANCHES:")
            for branch in self.newBranches:
                print(" > \x1b[32m{name}\x1b[0m {formula}".format(
                    name=(branch['name']+('[%d]'%branch['length'] if branch['length'] > 1 else '')).ljust(30),
                    formula=branch['formula'] if 'formula' in branch else 'function:\x1b[33m{fct}\x1b[0m'.format(fct=branch['function'].__name__)
                ))

        # callbacks before loop
        self._runCallbacks('beforeLoop')

        print ("------------------")
        print (" start processing ")
        print ("------------------")
        # loop over all events and write to output branches
        for event in self:

            # new event callback
            if self.callbacks and 'event' in self.callbacks:
                ## old: if callbacks return false, skip event but still evaluate callbacks per event
                #callbackResults = [fcn(event) for fcn in self.callbacks['event']]
                #if any(x==False for x in callbackResults):
                #    continue

                # NEW: if callbacks return false, skip event! don't evaluate all other callbacks
                if any(fcn(event)==False for fcn in self.callbacks['event']):
                    continue

            # fill branches
            for branch in self.newBranches:
                # evaluate result either as function applied on the tree entry or as TTreeFormula
                if branch['length'] == 1 and not 'arrayStyle' in branch:
                    if 'function' in branch:
                        if 'arguments' in branch:
                            branchResult = branch['function'](event, arguments=branch['arguments'])
                        else:
                            branchResult = branch['function'](event)
                    else:
                        branchResult = self.evaluate(branch['formula'])
                    if 'type' in branch and branch['type'] == 'i':
                        branchResult = int(branchResult)
                    # fill it for all the output trees
                    for outputTree in self.outputTrees:
                        outputTree['newBranchArrays'][branch['name']][0] = branchResult
                # for arrays pass the pointer to the array to the evaluation function to save the list->array conversion
                else:
                    if 'function' in branch:
                        # todo: make it more efficient by using a shared memory block for all of the output trees'
                        # todo: branches, this would help in case one adds new branches and writes to several trees at once
                        for outputTree in self.outputTrees:
                            if 'arguments' in branch:
                                branch['function'](event, destinationArray=outputTree['newBranchArrays'][branch['name']], arguments=branch['arguments'])
                            else:
                                branch['function'](event, destinationArray=outputTree['newBranchArrays'][branch['name']])
                    else:
                        for outputTree in self.outputTrees:
                            self.evaluateArray(branch['formula'], destinationArray=outputTree['newBranchArrays'][branch['name']])

            # evaluate all formulas
            self.formulaResults = {}
            for formulaName, formula in self.formulas.iteritems():
                self.formulaResults[formulaName] = self.evaluate(formulaName)

            # evaluate cuts for all output trees
            for outputTree in self.outputTrees:

                # evaluate all cuts of the sequence and abort early if one is not satisfied
                if outputTree['cutSequenceMode'] == 'AND':
                    passedCut = True
                    for cutFormulaName in outputTree['cutSequence']:
                        passedCut = passedCut and self.formulaResults[cutFormulaName]
                        if not passedCut:
                            break
                elif outputTree['cutSequenceMode'] == 'OR':
                    passedCut = False
                    for cutFormulaName in outputTree['cutSequence']:
                        passedCut = passedCut or self.formulaResults[cutFormulaName]
                        if passedCut:
                            break
                elif outputTree['cutSequenceMode'] == 'TREE':
                    passedCut = self.evaluateCutDictRecursive(outputTree['cutSequence'])
                else:
                    raise Exception("InvalidCutSequenceMode")
                

                # fill event if it passed the selection
                if passedCut:
                    
                    # fill index variables for friend trees
                    if 'friend' in outputTree:
                        outputTree['friend']['addr_run'][0] = self.tree.run
                        outputTree['friend']['addr_event'][0] = self.tree.event

                    outputTree['tree'].Fill()
                    outputTree['passed'] += 1

        print('INFO: end of processing. time ', time.ctime())
        passedTime = time.time() - self.timeStart
        perfStats = 'INPUT: {erps:1.4f}/s, OUTPUT: {ewps:1.4f}/s '.format(erps=self.eventsRead / passedTime if passedTime>0 else 0, ewps=sum([x['passed'] for x in self.outputTrees]) / passedTime if passedTime>0 else 0)
        print('INFO: throughput:', perfStats)
        print('INFO: number of trees processed:', len(self.processedTrees))
        sys.stdout.flush()
      
        if not self.checkProcessingComplete():
            print('ERROR: some files might have gone unavailable during iteration of the TChain!')
            print('ERROR:  -> expected number of trees:', len(self.nonEmptyTrees))
            print('ERROR:  -> actual number of processed trees:', len(self.processedTrees))
            raise Exception('FileException')

        # callbacks after loop 
        self._runCallbacks('finish')

        # build indices for friend trees
        for outputTree in self.outputTrees:
            if 'friend' in outputTree:
                if self.debug:
                    print("DEBUG: building index for friend tree...")
                indexStatus = outputTree['tree'].BuildIndex(self.friendTreeIndex[0], self.friendTreeIndex[1])
                if indexStatus < 0:
                    print("\x1b[31mERROR: building the index of friend tree failed with status ", indexStatus, "\x1b[0m")
                    raise Exception("sampleTree__process_do__build_index_failed")
        
        # callbacks before write
        self._runCallbacks('beforeWrite')

        # write files
        for outputTree in self.outputTrees:
            outputTree['file'].Write()
            outputTree['file'].Close()
        print('INFO: files written')
        print('INFO: saveMemory is ', self.saveMemory)
        sys.stdout.flush()

        if self.saveMemory:
            if not self.sequentialProcessing:
                self.tree.Reset()
                self.tree = None

            for outputTree in self.outputTrees:
                outputTree['tree'] = None

            print('INFO: trees in memory destroyed!')

        # callbacks after having written file
        self._runCallbacks('afterWrite')

        print('INFO: done. time ', time.ctime(), ' events read:', self.eventsRead)
        sys.stdout.flush()

        for outputTree in self.outputTrees:
            passedSelectionFraction = 100.0*outputTree['passed']/self.eventsRead if self.eventsRead>0 else '?'
            print (' > \x1b[34m{name}\x1b[0m {passed} ({fraction}%) => {outputFile}'.format(name=outputTree['name'], passed=outputTree['passed'], fraction=passedSelectionFraction, outputFile=outputTree['fileName']))
        sys.stdout.flush()

        # delete temporary files
        for tempFileName in self.temporaryFiles:
            try:
                if os.path.isfile(tempFileName):
                    self.fileLocator.rm(tempFileName)
                    if self.debug:
                        print("DEBUG: delete ", tempFileName)
            except:
                pass

    def checkProcessingComplete(self):
        if self.debug:
            print('DEBUG: verify that all tree files were seen once during iteration') 
            print('DEBUG:  -> number of trees:', len(self.sampleFileNames))
            print('DEBUG:  -> number of non-empty trees:', len(self.nonEmptyTrees))
            print('DEBUG:  -> actual number of processed trees:', len(self.processedTrees))
        return len(self.nonEmptyTrees) == len(self.processedTrees)

    @staticmethod
    def countSampleFiles(samples):
        # get list of sample root files
        if type(samples) == list:
            return len(samples)
        else:
            sampleTextFileName = samples
            if os.path.isfile(sampleTextFileName):
                with open(sampleTextFileName, 'r') as sampleTextFile:
                    sampleFileNames = sampleTextFile.readlines()
                return len(sampleFileNames)
            else:
                print ('ERROR: sample list text file does not exist:', sampleTextFileName)
        return -1

    def getNumSampleFiles(self):
        return len(self.sampleFileNames)

    # return the total scale for the sample, calculated from all count histograms from the TChain
    def getScale(self, sample, countHistogram=None):
        try:
            sample.xsec = sample.xsec[0]
        except Exception as e:
            if self.verbose:
                print("DEBUG: XS = ", sample.xsec, sample, type(sample))

        #print('self.totalNanoTreeCounts', self.totalNanoTreeCounts)  # {'run': 19L, 'genEventSumw2': 1978563.8623046875, 'nLHEPdfSumw': 0L, 'nLHEScaleSumw': 0L, 'genEventSumw': 1978776.7690429688, 'genEventCount': 1979000L}

        if self.totalNanoTreeCounts:
            if self.config.has_option('Configuration', 'countsFromAutoPU') and eval(self.config.get('Configuration', 'countsFromAutoPU')): #countsFromAutoPU is False by default.
                count = self.histograms['autoPU'].GetEntries()
                countHistogram = "autoPU.GetEntries()"
            else:
                if not countHistogram:
                    countHistogram = self.config.get('Configuration', 'countTreeName') if self.config.has_option('Configuration', 'countTreeName') else 'genEventSumw'
                # NANOAODv6 workaround
                if countHistogram not in self.totalNanoTreeCounts and countHistogram + '_' in self.totalNanoTreeCounts:
                    countHistogram = countHistogram + '_'
                count = self.totalNanoTreeCounts[countHistogram]
                #print('count', count)  #1978776.76904
        else:
            if not countHistogram:
                try:
                    posWeight = self.histograms['CountPosWeight'].GetBinContent(1)
                    negWeight = self.histograms['CountNegWeight'].GetBinContent(1)
                    count = posWeight - negWeight
                    countHistogram = 'CountPosWeight - CountNegWeight'
                except:
                    if self.verbose:
                        print("sampleTree: no CountPosWeight/CountNegWeight: using Count instead!!!!!!!!!!!")
                    try:
                        count = self.histograms['Count'].GetBinContent(1)
                    except Exception as e:
                        print ("EXCEPTION:", e)
                        print ("ERROR: no weight histograms found in sampleTree => terminate")
                        print ("HISTOGRAMS:", self.histograms)
                        exit(0)
            else:
                count = self.histograms[countHistogram].GetBinContent(1)

        # override event counts: config needs a section 'EventCounts' with sample identifier strings as keys and the new count as value
        # [EventCounts]
        # SampleIdentifier = 12345
        try:
            if self.sampleIdentifier and self.config.has_section('EventCounts') and self.config.has_option('EventCounts', self.sampleIdentifier): #not there
                countNew = eval(self.config.get('EventCounts', self.sampleIdentifier))
                print("\x1b[97m\x1b[41mINFO: overwrite event counts with values from config!!!\n value from file:", count, "\n value from config:", countNew," <--- will be used!\x1b[0m")
                count = countNew
            #else:
            #    print("--> don't overwrite counts!", self.sampleIdentifier, self.config.has_section('EventCounts'), self.config.has_option('EventCounts', self.sampleIdentifier))
        except Exception as e:
            print("\x1b[31mException:",e," -> overwriting of event counts has been disabled\x1b[0m")


        lumi = float(sample.lumi)
        theScale = lumi * sample.xsec * sample.sf / float(count)

        if self.verbose:
            print("sampleTree.getScale(): sample: ", sample, "lumi: ", lumi, "xsec: ", sample.xsec, "sample.sf: ", sample.sf, "count (", countHistogram, "):", count, " ---> using scale: ", theScale)
        return theScale

    def getCollection(self, name):
        collection = lambda: None
        branches = [x.GetName() for x in self.GetListOfBranches()]
        numBranch = 'n' + name
        if numBranch in branches:
            collection.size = lambda: getattr(self.tree, numBranch)
            for x in branches:
                if x.startswith(name + '_'):
                    setattr(collection, x.split(name + '_')[1], getattr(self.tree, x))
        return collection

    # create a unique string representation of the total cut, e.g. used to calculate the hash for cached samples 
    # this is not required to be a 'real' cut string, used by TTreeFormula etc.
    @staticmethod
    def findMinimumCut(cutList, cutSequenceMode='AND'):
        if type(cutList) == list or type(cutList) == dict:
            cuts = cutList
        else:
            cuts = [cutList]
        if cutSequenceMode == 'TREE' or type(cutList) == dict:
            minCut = "%r"%cuts
        elif cutSequenceMode == 'AND':
            minCut = '&&'.join(['(%s)'%x.replace(' ', '') for x in sorted(cuts)])
        elif cutSequenceMode == 'OR':
            minCut = '||'.join(['(%s)'%x.replace(' ', '') for x in sorted(list(set(cuts)))])
        else:
            minCut = "%r"%cuts
        #print('\x1b[31mmin cut:',minCut,'\x1b[0m')
        return minCut

    def GetEntries(self):
        return self.tree.GetEntries()

    def Print(self):
        print("\x1b[34m\x1b[1m---- SampleTree ----")
        print("# this snippet below can be used to load this sample:")
        if len(self.sampleFileNamesToProcess) != len(self.sampleFileNames):
            print("# commented file names below mean only a part of sample was being processed in this job")
            print("# normalization is always correct for the full sample only (or if all parts are added afterwards)")
        print("import ROOT")
        print("from myutils.sampleTree import SampleTree")
        print("sampleTree = SampleTree([")
        for fileName in self.sampleFileNames:
            if fileName in self.sampleFileNamesToProcess:
                print("    '" + fileName + "',")
            else:
                print("#    '" + fileName + "',")
        print("], treeName='Events', xrootdRedirector='" + self.fileLocator.getXrootdRedirector() + "')")
        print("---- end ----\x1b[0m")

if __name__ == '__main__':

    str_sample = 'ZZ_TuneCP5_13TeV-pythia8'
    str_sample2 = 'WW_TuneCP5_13TeV-pythia8'
    config = XbbConfigReader.read('Zll2018')
    directory = config.get('Directories', 'sysOUT')
    #config = XbbConfigReader.read('Zll2018')
    
    path = 'xyz'
    sampleInfo = ParseInfo(config,path,config=config)
    sample = [x for x in sampleInfo if x.identifier==str_sample and x.subsample==False]
     

    sampleTree = SampleTree({'sample': sample[0], 'folder': directory}, config=config) 
    
    scaleXStoLumi = sampleTree.getScale(sample[0])

    weightExpression = config.get('Weights','weightF')
    signalRegionSelection = config.get('Cuts', 'ZllBDT_lowpt')
 
    sampleTree.enableBranches(BranchList([weightExpression,signalRegionSelection]).getListOfBranches())

    sampleTree.addFormula(weightExpression)

