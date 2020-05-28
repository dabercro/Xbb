#!/usr/bin/env python
import sys
import os
import ROOT
ROOT.gROOT.SetBatch(True)
from optparse import OptionParser
from myutils.FileList import FileList
from myutils import BetterConfigParser, ParseInfo, LeptonSF
from myutils.FileLocator import FileLocator
import importlib
from myutils.sampleTree import SampleTree
import resource
import uuid

# ntuple processing class (former "sys-step")

class XbbRun:

    def __init__(self, opts):

        # get file list
        self.filelist = FileList.decompress(opts.fileList) if len(opts.fileList) > 0 else None
        print "len(filelist)",len(self.filelist),
        if len(self.filelist) > 0:
            print "filelist[0]:", self.filelist[0]
        else:
            print ''

        # config
        self.debug = 'XBBDEBUG' in os.environ
        self.verifyCopy = True
        self.opts = opts
        self.config = BetterConfigParser()
        self.config.read(opts.config)
        self.channel = self.config.get('Configuration', 'channel')

        # load namespace, TODO
        VHbbNameSpace = self.config.get('VHbbNameSpace', 'library')
        ROOT.gSystem.Load(VHbbNameSpace)

        # directories
        self.pathIN = self.config.get('Directories', opts.inputDir)
        self.pathOUT = self.config.get('Directories', opts.outputDir)
        self.tmpDir = self.config.get('Directories', 'scratch')
        print 'INput samples:\t%s'%self.pathIN
        print 'OUTput samples:\t%s'%self.pathOUT

        self.fileLocator = FileLocator(config=self.config)

        # check if given sample identifier uniquely matches a samples from config
        matchingSamples = ParseInfo(samples_path=self.pathIN, config=self.config).find(identifier=opts.sampleIdentifier)
        if len(matchingSamples) != 1:
            print "ERROR: need exactly 1 sample identifier as input with -S !!"
            print matchingSamples
            exit(1)
        self.sample = matchingSamples[0]

        # collections
        self.collections = [x.strip() for x in opts.addCollections.split(',') if len(x.strip()) > 0] if len(opts.addCollections.strip())>0  else []
        if len(self.collections) < 1:
            print "\x1b[31mWARNING: no collections added! Specify the collections to add with the --addCollections option!\x1b[0m"
        print 'collections to add:', self.collections
        self.collections = self.parseCollectionList(self.collections)
        print 'after parsing:', self.collections

        # temorary folder to save the files of this job on the scratch
        temporaryName = self.sample.identifier + '/' + uuid.uuid4().hex

        # input files
        self.subJobs = []
        if opts.join:
            print("INFO: join input files! This is an experimental feature!")

            # translate naming convention of .txt file to imported files after the prep step
            inputFileNamesAfterPrep = [self.fileLocator.getFilenameAfterPrep(x) for x in self.filelist]

            self.subJobs.append({
                'inputFileNames': self.filelist,
                'localInputFileNames': ["{path}/{subfolder}/{filename}".format(path=self.pathIN, subfolder=self.sample.identifier, filename=localFileName) for localFileName in inputFileNamesAfterPrep],
                'outputFileName': "{path}/{subfolder}/{filename}".format(path=self.pathOUT, subfolder=self.sample.identifier, filename=inputFileNamesAfterPrep[0]),
                'tmpFileName': "{path}/{subfolder}/{filename}".format(path=self.tmpDir, subfolder=temporaryName, filename=inputFileNamesAfterPrep[0]),
                })

        else:
            
            # create separate subjob for all files (default!)
            for inputFileName in self.filelist:
                inputFileNamesAfterPrep = [self.fileLocator.getFilenameAfterPrep(inputFileName)]

                self.subJobs.append({
                    'inputFileNames': [inputFileName],
                    'localInputFileNames': ["{path}/{subfolder}/{filename}".format(path=self.pathIN, subfolder=self.sample.identifier, filename=localFileName) for localFileName in inputFileNamesAfterPrep],
                    'outputFileName': "{path}/{subfolder}/{filename}".format(path=self.pathOUT, subfolder=self.sample.identifier, filename=inputFileNamesAfterPrep[0]),
                    'tmpFileName': "{path}/{subfolder}/{filename}".format(path=self.tmpDir, subfolder=temporaryName, filename=inputFileNamesAfterPrep[0]),
                    })

    # lists of single modules can be given instead of a module, "--addCollections Sys.all"
    # [Sys]
    # all = ['Sys.Vtype', 'Sys.Leptons', ...]
    # TODO: make it fully recursive
    def parseCollectionList(self, collections): 
        collectionsListsReplaced = []
        for collection in collections:
            if '.' in collection:
                section = collection.split('.')[0]
                key = collection.split('.')[1]
                listExpression = self.config.get(section, key).strip()
                if listExpression.startswith('[') and listExpression.endswith(']'):
                    listParsed = eval(listExpression)
                    for i in listParsed:
                        collectionsListsReplaced.append(i)
                else:
                    collectionsListsReplaced.append(collection)
            else:
                collectionsListsReplaced.append(collection)

        if len(collectionsListReplaced) != len(set(collectionsListReplaced)):
            print('Collection list has some duplicate. Check .ini files.')
            print(collectionsListReplaced)
            raise Exception("ConfigError")

        return collectionsListsReplaced

    # run all subjobs
    def run(self):

        nFilesProcessed = 0
        nFilesFailed = 0

        for subJob in self.subJobs:

            # only process if output is non-existing/broken or --force was used
            if self.opts.force or not self.fileLocator.isValidRootFile(subJob['outputFileName']):

                # create directories
                outputFolder = '/'.join(subJob['outputFileName'].split('/')[:-1])
                tmpFolder = '/'.join(subJob['tmpFileName'].split('/')[:-1])
                self.fileLocator.makedirs(outputFolder)
                self.fileLocator.makedirs(tmpFolder)

                # load sample tree
                sampleTree = SampleTree(subJob['localInputFileNames'], config=self.config)
                if not sampleTree.tree:
                    print "trying fallback...", len(subJob['inputFileNames'])

                    if len(subJob['inputFileNames']) == 1:
                        # try original naming scheme if reading directly from Heppy/Nano ntuples (without prep)
                        fileNameOriginal = self.pathIN + '/' + subJob['inputFileNames'][0]
                        print "FO:", fileNameOriginal
                        xrootdRedirector = self.fileLocator.getRedirector(fileNameOriginal)
                        sampleTree = SampleTree([fileNameOriginal], config=self.config, xrootdRedirector=xrootdRedirector)
                        if not sampleTree.tree:
                            print "\x1b[31mERROR: file does not exist or is broken, will be SKIPPED!\x1b[0m"
                            nFilesFailed += 1
                            continue
                    else:
                        print "\x1b[31mERROR: file does not exist or is broken, will be SKIPPED! (old naming scheme not supported for joining multipel files)\x1b[0m"
                        nFilesFailed += 1
                        continue

                # to use this syntax, use "--addCollections Sys.Vtype" for a config file entry like this:
                # [Sys]
                # Vtype = VtypeCorrector.VtypeCorrector(channel='Zll')
                # (instead of passing the tree in the constructor, the setTree method can be used)
                pyModules = []
                versionTable = []
                for collection in self.collections:
                    if '.' in collection:
                        section = collection.split('.')[0]
                        key = collection.split('.')[1]
                        if self.config.has_section(section) and self.config.has_option(section, key):
                            pyCode = self.config.get(section, key)
                        elif '(' in collection and collection.endswith(')'):
                            print "WARNING: config option", collection, " not found, interpreting it as Python code!"
                            pyCode = collection 
                        else:
                            print "\x1b[31mERROR: config option not found:", collection, ". To specify Python code directly, pass a complete constructor, e.g. --addCollections 'Module.Class()'. Module has to be placed in python/myutils/ folder.\x1b[0m"
                            raise Exception("ConfigError")

                        # import module from myutils
                        moduleName = pyCode.split('(')[0].split('.')[0].strip()
                        if self.debug:
                            print "DEBUG: import module:", moduleName
                            print("\x1b[33mDEBUG: " + collection + ": run PYTHON code:\n"+pyCode+"\x1b[0m")
                        globals()[moduleName] = importlib.import_module(".{module}".format(module=moduleName), package="myutils")

                        # get object
                        wObject = eval(pyCode)

                        # pass the tree and other variables if needed to finalize initialization
                        if hasattr(wObject, "customInit") and callable(getattr(wObject, "customInit")):
                            wObject.customInit({'config': self.config,
                                                'sampleTree': sampleTree,
                                                'tree': sampleTree.tree,
                                                'sample': self.sample,
                                                'channel': self.channel,
                                                'pathIN': self.pathIN,
                                                'pathOUT': self.pathOUT,
                                                })

                        # add callbacks if the objects provides any
                        if hasattr(wObject, "processEvent") and callable(getattr(wObject, "processEvent")):
                            sampleTree.addCallback('event', wObject.processEvent)
                        for cb in ["finish", "prepareOutput"]:
                            if hasattr(wObject, cb) and callable(getattr(wObject, cb)):
                                sampleTree.addCallback(cb, getattr(wObject, cb))

                        # add branches
                        if hasattr(wObject, "getBranches") and callable(getattr(wObject, "getBranches")):
                            sampleTree.addOutputBranches(wObject.getBranches())

                        pyModules.append(wObject)

                        versionTable.append([moduleName, wObject.getVersion() if hasattr(wObject, "getVersion") else 0])
                    else:
                        print "\x1b[31mERROR: config option not found:", collection, " the format should be: [Section].[Option]\x1b[0m"
                        raise Exception("ConfigError")

                for moduleName, moduleVersion in versionTable:
                    print " > {m}:{v}".format(m=moduleName, v=moduleVersion)

                # DEPRECATED, do not use anymore ---> use BranchTools.TreeFormulas()
                if 'addbranches' in self.collections:
                    writeNewVariables = eval(self.config.get("Regression", "writeNewVariablesDict"))
                    sampleTree.addOutputBranches(writeNewVariables)
                
                # DEPRECATED, do not use anymore ---> use BranchTools.Drop()
                if 'removebranches' in self.collections:
                    bl_branch = eval(config.get('Branches', 'useless_branch'))
                    for br in bl_branch:
                        sampleTree.addBranchToBlacklist(br)
                    bl_branch = eval(config.get('Branches', 'useless_after_sys'))
                    for br in bl_branch:
                        sampleTree.addBranchToBlacklist(br)

                # define output file 
                sampleTree.addOutputTree(subJob['tmpFileName'], cut='1', branches='*', friend=self.opts.friend)

                # run processing
                for pyModule in pyModules:
                    if hasattr(pyModule, "beforeProcessing"):
                        getattr(pyModule, "beforeProcessing")()

                sampleTree.process()

                for pyModule in pyModules:
                    if hasattr(pyModule, "afterProcessing"):
                        getattr(pyModule, "afterProcessing")()

                # if output trees have been produced: copy temporary file to output folder
                if sampleTree.getNumberOfOutputTrees() > 0: 
                    try:
                        self.fileLocator.cp(subJob['tmpFileName'], subJob['outputFileName'], force=True)
                        print 'copy ', subJob['tmpFileName'], subJob['outputFileName']

                        if self.verifyCopy:
                            if not self.fileLocator.isValidRootFile(subJob['outputFileName']):
                                print 'INFO: output at final destination broken, try to copy again from scratch disk to final destination...'
                                self.fileLocator.cp(subJob['tmpFileName'], subJob['outputFileName'], force=True)
                                print 'INFO: second attempt copy done!'
                                if not self.fileLocator.isValidRootFile(subJob['outputFileName']):
                                    print '\x1b[31mERROR: output still broken!\x1b[0m'
                                    nFilesFailed += 1
                                    raise Exception("FileCopyError")
                                else:
                                    print 'INFO: file is good after second attempt!'
                    except Exception as e:
                        print e
                        print "\x1b[31mERROR: copy from scratch to final destination failed!!\x1b[0m"

                    # delete temporary file
                    try:
                        self.fileLocator.rm(subJob['tmpFileName'])
                    except Exception as e:
                        print e
                        print "WARNING: could not delete file on scratch!"


                # clean up
                if hasattr(wObject, "cleanUp") and callable(getattr(wObject, "cleanUp")):
                    getattr(wObject, "cleanUp")()

            else:
                print 'SKIP:', subJob['inputFileNames']

        if nFilesFailed > 0:
            raise Exception("ProcessingIncomplete")

argv = sys.argv
parser = OptionParser()
parser.add_option("-S", "--sampleIdentifier", dest="sampleIdentifier", default="", help="samples you want to run on")
parser.add_option("-C", "--config", dest="config", default=[], action="append", help="configuration defining the plots to make")
parser.add_option("-f", "--fileList", dest="fileList", default="", help="list of files you want to run on")
parser.add_option("-b", "--addCollections", dest="addCollections", default="", help="collections to add, e.g. Sys.all")
parser.add_option("-F", "--force", dest="force", action="store_true", help="overwrite existing files", default=False)
parser.add_option("-J", "--join", dest="join", action="store_true", help="chain all files of the sample", default=False)
parser.add_option("-d", "--friend", dest="friend", action="store_true", help="create a friend tree", default=False)
parser.add_option("-I", "--inputDir", dest="inputDir", default="SYSin", help="name of input folder in config, e.g. SYSin")
parser.add_option("-O", "--outputDir", dest="outputDir", default="SYSout", help="name of output folder in config, e.g. SYSout")


(opts, args) = parser.parse_args(argv)
if opts.config == "":
    opts.config = "config"

xr = XbbRun(opts)
xr.run()

