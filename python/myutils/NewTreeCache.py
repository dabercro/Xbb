from __future__ import print_function
import glob
from Hash import Hash
import ROOT
import subprocess
import os
import sys
from samplesclass import Sample
from sampleTree import SampleTree as SampleTree
from FileLocator import FileLocator
import resource

# ------------------------------------------------------------------------------
# TreeCache
#
# cache for ROOT trees (abstracted by sampleTree class!) with specific cuts or
# branch selections
#
# usage e.g.:
#
# tc = TreeCache.TreeCache(
#     sample='ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8_ext1',
#     cutList='V_pt>100',
#     inputFolder='/scratch/p/',
#     outputFolder='/scratch/p/tmp/',
#     debug=True
# )
#
# #obtain sampleTree of already cached tree:
#
# if tc.isCached():
#    sampleTree = tc.getTree()
#
# #cache one (or several at a time!) new tree
#
# else:
#    sampleTreeBkg = SampleTree('test/bkg.txt')
#    tc.setSampleTree(sampleTreeBkg)
#    tc.cache()
#    tc2.setSampleTree(sampleTreeBkg)
#    tc2.cache()
#    tc3.setSampleTree(sampleTreeBkg)
#    tc3.cache()
#    sampleTreeBkg.process()
#
# one can add many different cuts to process at the same time with multiple
# calls of:
#    tc.setSampleTree(sampleTreeBkg)
#    tc.cache()
# with different tc objects and the same sampleTreeBkg object. In the end, one
# has to run:
#    sampleTreeBkg.process()
# to write the files
# ------------------------------------------------------------------------------

class TmpOutputName:
    def __init__(self, formatstring):
        self.outputFileNameFormat = formatstring

    def format(self, **kwargs):
        outputFileDirectory = '{outputFolder}/{hash1}/{hash2}'.format(
            outputFolder=kwargs['outputFolder'],
            hash1=kwargs['hash'][0:2],
            hash2=kwargs['hash'][2:4]
            )

        if not os.path.exists(outputFileDirectory):
            os.makedirs(outputFileDirectory)
        
        kwargs['outputFolder'] = outputFileDirectory

        return self.outputFileNameFormat.format(**kwargs)
    

class TreeCache:

    def __init__(self, sample, cutList='1', branches=None, inputFolder=None, tmpFolder=None, outputFolder=None, chunkNumber=-1, splitFilesChunks=-1, splitFilesChunkSize=-1, debug=False, fileList=None, cutSequenceMode='AND', name='', config=None, fileLocator=None):
        self.config = config
        self.fileLocator = fileLocator if fileLocator is not None else FileLocator(config=self.config)
        self.debug = debug or ('XBBDEBUG' in os.environ)

        # SAMPLE
        if isinstance(sample, Sample):
            # sample passed as Sample object
            # count number of chunks the cached data is split into
            defaultChunkSize = int(config.get('General', 'mergeCachingSize')) if config.has_option('General', 'mergeCachingSize') else 100
            splitFilesChunkSize = sample.mergeCachingSize if sample.mergeCachingSize > 0 else defaultChunkSize
            splitFilesChunks = SampleTree({'name': sample.identifier, 'folder': inputFolder}, countOnly=True, splitFilesChunkSize=splitFilesChunkSize, config=config, verbose=self.debug, fileLocator=self.fileLocator).getNumberOfParts()
            # if sample passed as object, it can be a 'subsample' and habe different name and identifier
            self.sample = sample.name
            self.sampleIdentifier = sample.identifier
            if self.debug:
                print ("INFO: use sample=", sample.name, " #parts = ", splitFilesChunks)
        else:
            # sample identifier passed as string
            self.sample = sample
            self.sampleIdentifier = sample
        self.name = name

        # CUTS
        self.cutList = cutList
        self.cutSequenceMode = cutSequenceMode
        self.minCut = SampleTree.findMinimumCut(self.cutList, cutSequenceMode=self.cutSequenceMode)

        # PATHS
        self.inputFolder = inputFolder
        self.outputFolder = (config.get('Directories', 'tmpSamples') if config else 'cache/') if outputFolder is None else outputFolder
        self.tmpFolder = (config.get('Directories', 'scratch') if config else 'tmp/') if tmpFolder is None else tmpFolder
        self.cachedFileNames = []
        self.tmpFiles = []
        self.outputFileNameFormat = TmpOutputName('{outputFolder}/tmp_{hash}_{part}of{parts}.root')

        # BRANCHES and chunk information
        self.branches = branches
        self.branchesForHash = None     # for now make hash independent of selecte branches 
        self.hash = Hash(sample=sample, minCut=self.minCut, branches=self.branchesForHash, splitFilesChunkSize=splitFilesChunkSize, debug=False, inputPath=self.inputFolder).get()
        self.chunkNumber = chunkNumber
        self.splitFilesChunks = splitFilesChunks if splitFilesChunks > 1 else 1
        self.splitFilesChunkSize = splitFilesChunkSize
        
        # identifier is just used as an arbitrary name for print-out
        cutUsedForIdentifier = (self.minCut if len(self.minCut) < 60 else self.minCut[0:50] + '...').replace(' ', '')
        self.identifier = '{sample}[{cut}]of{parts}'.format(sample=self.sample, cut=cutUsedForIdentifier, parts=self.splitFilesChunks)
        self.sampleTree = None
        self.isCachedChecked = False

        self.createFolders()

    # free memory
    def deleteSampleTree(self):
        self.sampleTree = None

    # file, where skimmed tree is written to
    def getTmpFileName(self):
        return self.outputFileNameFormat.format(
            outputFolder=self.tmpFolder,
            hash=self.hash,
            part=self.chunkNumber if self.chunkNumber > 0 else 1,
            parts='%d'%self.splitFilesChunks
        )

    # file, where skimmed tree is moved to after it has been written completely
    def getOutputFileName(self):
        return self.outputFileNameFormat.format(
            outputFolder=self.outputFolder,
            hash=self.hash,
            part=self.chunkNumber if self.chunkNumber > 0 else 1,
            parts='%d'%self.splitFilesChunks
        )

    # check existence of files with skimmed trees
    def findCachedFileNames(self, chunkNumber=-1):
        cachedFilesMaskRaw = self.outputFileNameFormat.format(
            outputFolder=self.outputFolder,
            hash=self.hash,
            part='*' if chunkNumber < 1 else '%d'%chunkNumber,
            parts=self.splitFilesChunks
        )
        cachedFilesMask = self.fileLocator.getLocalFileName(cachedFilesMaskRaw)

        # this does not work reliably on T3 worker nodes anymore 
        #self.cachedFileNames = glob.glob(cachedFilesMask)

        # workaround: use xrootd for directory listing 
        #self.cachedFileNames = self.fileLocator.glob_with_fallback(cachedFilesMaskRaw)

        # this solution uses a loop over all possible files and uses xrdfs stat instead of xrdfs ls
        self.cachedFileNames = self.fileLocator.get_numbered_file_list(cachedFilesMaskRaw, 1, self.splitFilesChunks)

        if self.debug:
            print ('DEBUG: search files:', cachedFilesMask)
            print ('\x1b[32mDEBUG: files:')
            for fileName in self.cachedFileNames:
                print (' > ', fileName)
            if len(self.cachedFileNames) < 1:
                print ('none!')
            print ('\x1b[0m(%d files found)'%len(self.cachedFileNames))

        # sort
        self.cachedFileNames = sorted(self.cachedFileNames, key=lambda x: int(x.split('_')[-1].split('of')[0]) if 'of' in x and '_' in x else -1)
        return self.cachedFileNames

    def getTotalNumberOfOutputFiles(self):
        return self.splitFilesChunks

    # check if a single part is cached, (only checks existence of the file, not validity!)
    def partIsCached(self):
        cachedFilesMaskRaw = self.outputFileNameFormat.format(
            outputFolder=self.outputFolder,
            hash=self.hash,
            part=self.chunkNumber,
            parts=self.splitFilesChunks
        )
        # this does not work reliably on T3 worker nodes anymore 
        #cachedFilesMask = self.fileLocator.getLocalFileName(cachedFilesMaskRaw)
        #return len(glob.glob(cachedFilesMask)) > 0
        return self.fileLocator.remoteFileExists(cachedFilesMaskRaw)

    # isCached == all files containing the skimmed tree found!
    def isCached(self):
        self.findCachedFileNames()
        if (len(self.cachedFileNames) != self.splitFilesChunks and self.splitFilesChunks > 1) or len(self.cachedFileNames) == 0:
            if self.debug:
                print ('\x1b[32mDEBUG: not cached:', self.identifier, '\x1b[0m')
            return False
        self.isCachedChecked = True
        return True

    # check if an existing file can be opened without errors by ROOT
    def checkFileValidity(self, rawFileName):
        xrootdFileName = self.fileLocator.getXrootdFileName(rawFileName)
        f = ROOT.TFile.Open(xrootdFileName, 'read')
        if not f or f.GetNkeys() == 0 or f.TestBit(ROOT.TFile.kRecovered) or f.IsZombie():
            print ('\x1b[31mWARNING: broken file:', rawFileName, ' => redo caching!\x1b[0m')
            if f:
                f.Close()
            self.deleteFile(rawFileName)
            return False
        if f:
            f.Close()
        return True

    # check if all cached files are valid
    def isCachedAndValid(self):
        valid = True
        if self.isCached():
            # check file integrity
            for fileName in self.cachedFileNames:
                valid = valid and self.checkFileValidity(fileName)
        else:
            valid = False
        return valid

    # set input sampleTree object
    def setSampleTree(self, sampleTree):
        self.sampleTree = sampleTree
        return self

    # this prepares the caching by telling the sampleTree object what to write during processing of the file
    # note: does not run the caching by itself! needs an additional sampleTree.process()
    def cache(self):
        if self.sampleTree:
            outputFileName = self.getTmpFileName()
            callbacks = {'afterWrite': self.moveFilesToFinalLocation}
            self.sampleTree.addOutputTree(outputFileName=outputFileName, cut=self.cutList, hash=self.hash, branches=self.branches, callbacks=callbacks, cutSequenceMode=self.cutSequenceMode, name=self.name)
            self.tmpFiles.append(outputFileName)
            if self.debug:
                print ('\x1b[32mDEBUG: output file for ', self.identifier, ' is ', outputFileName, '\x1b[0m')
        else:
            print ('\x1b[31mERROR: no sample tree connected!:', self.identifier, ' set the sampleTree first with "setSampleTree(sampleTree)" \x1b[0m')
        return self

    # return sample tree class of cached samples if all files found
    def getTree(self, chunkSize=-1, chunkNumber=-1):
        # if it has already been checked if tree is cached, then use this result dierctly
        isCached = self.isCachedChecked
        if not isCached:
            isCached = self.isCached()
        if isCached:
            if chunkSize > 0 and chunkNumber > 0:
                fileNames = self.cachedFileNames[(chunkNumber-1)*chunkSize:chunkNumber*chunkSize]
            elif chunkSize < 0 and chunkNumber < 0:
                fileNames = self.cachedFileNames
            else:
                raise Exception("InvalidParameters")
            self.sampleTree = SampleTree(self.cachedFileNames, config=self.config, fileNamesToProcess=fileNames)
            self.sampleTree.sampleIdentifier = self.sampleIdentifier

            # check if even though all files exist, they couldn't be accessed for some reason
            # and therefore the tree would be incomplete
            if not self.sampleTree.isCompleteTree():
                raise Exception("IncompleteTree")

        return self.sampleTree

    # delete file
    def deleteFile(self, rawFileName):
        if self.debug:
            print ('DELETE:', rawFileName)
        self.fileLocator.rm(rawFileName)

    # delete cached files
    def deleteCachedFiles(self, chunkNumber=-1):
        cachedFileNames = self.findCachedFileNames(chunkNumber=chunkNumber)
        for fileName in cachedFileNames:
            if self.fileLocator.fileExists(fileName):
                self.deleteFile(fileName)

    # create folders
    def createFolders(self):
        tmpfolderLocal = self.fileLocator.getLocalFileName(self.tmpFolder)
        if not os.path.isdir(tmpfolderLocal):
            print("DOES NOT EXIST:", tmpfolderLocal)
            try:
                xrootdFileName = self.fileLocator.getXrootdFileName(self.tmpFolder)
                if '://' not in xrootdFileName:
                    os.makedirs(self.tmpFolder)
                else:
                    command = 'gfal-mkdir %s' % (xrootdFileName)
                    returnCode = subprocess.call([command], shell=True)
                    if self.debug:
                        print(command, ' => ', returnCode)
                        print ()
            except:
                pass

        if not self.fileLocator.exists(self.outputFolder):
            print("INFO: output folder does not exist and will be created:", self.outputFolder)
            self.fileLocator.makedirs(self.outputFolder)

    # move files from temporary to final location
    def moveFilesToFinalLocation(self, raiseOnFailure=True):
        success = True
        # free some memory for file copy command
        if self.debug:
            print('DEBUG: max mem used A:', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        self.deleteSampleTree()
        if self.debug:
            print('DEBUG: max mem used B:', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)

        for tmpFileName in self.tmpFiles:
            outputFileName = self.outputFolder + '/' + self.tmpFolder.join(tmpFileName.split(self.tmpFolder)[1:])
            print ('copy ', tmpFileName, ' to ', outputFileName)

            if not os.path.exists(os.path.dirname(outputFileName)):
                os.makedirs(os.path.dirname(outputFileName))

            if self.fileLocator.fileExists(outputFileName):
                self.deleteFile(outputFileName)
            copySuccessful = self.fileLocator.cp(tmpFileName, outputFileName)
            if not copySuccessful:
                print("WARNING: first copy attempt failed! retry once!")
                self.fileLocator.debug = True
                copySuccessful = self.fileLocator.cp(tmpFileName, outputFileName)
                if not copySuccessful:
                    success = False
                    print('\x1b[31mERROR: copy failed for {tmpfile}->{outputfile} !\x1b[0m'.format(tmpfile=tmpFileName,
                                                                                                outputfile=outputFileName))
                    if raiseOnFailure:
                        raise Exception("CopyToFinalDestinationFailed")
            if success:
                # delete temporary file if copy was successful
                self.deleteFile(tmpFileName)
        return success
