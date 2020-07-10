#!/usr/bin/env python
from __future__ import print_function
from optparse import OptionParser
import ROOT
ROOT.gROOT.SetBatch(True)

from myutils import NewTreeCache as TreeCache
from myutils.sampleTree import SampleTree as SampleTree
from myutils.BranchList import BranchList

import sys

class CacheTraining(object):

    def __init__(self, config, sampleIdentifier, trainingRegions, splitFilesChunks=1, chunkNumber=1, splitFilesChunkSize=-1, force=False):
        self.config = config
        self.force = force
        self.sampleIdentifier = sampleIdentifier
        self.trainingRegions = trainingRegions

        self.sampleTree = None
        if config.has_option('Directories', 'trainingSamples'):
            self.samplesPath = self.config.get('Directories', 'trainingSamples')
        else:
            self.samplesPath = self.config.get('Directories', 'MVAin')
        self.samplesInfo = ParseInfo(samples_path=self.samplesPath, config=self.config) 
        self.sampleFilesFolder = self.config.get('Directories', 'samplefiles')

        self.backgroundSampleNames = list(set(sum([eval(self.config.get(trainingRegion, 'backgrounds')) for trainingRegion in self.trainingRegions], [])))
        self.signalSampleNames = list(set(sum([eval(self.config.get(trainingRegion, 'signals')) for trainingRegion in self.trainingRegions], [])))
        # can include DATA in the .h5 files for training
        self.dataSampleNames = list(set(sum([eval(self.config.get(trainingRegion, 'data')) if self.config.has_option(trainingRegion, 'data') else []  for trainingRegion in self.trainingRegions], [])))
        self.samples = self.samplesInfo.get_samples(list(set(self.backgroundSampleNames + self.signalSampleNames + self.dataSampleNames)))

        self.trainingRegionsDict = {}
        for trainingRegion in self.trainingRegions:
            treeCutName = config.get(trainingRegion, 'treeCut') if config.has_option(trainingRegion, 'treeCut') else trainingRegion
            treeVarSet = config.get(trainingRegion, 'treeVarSet').strip()
            #systematics = [x for x in config.get('systematics', 'systematics').split(' ') if len(x.strip())>0]
            if config.has_option(trainingRegion, 'systematics'):
                systematicsString = config.get(trainingRegion, 'systematics').strip()
                if systematicsString.startswith('['):
                    systematics = eval(systematicsString)
                else:
                    systematics = systematicsString.split(' ')
            else:
                systematics = []
            mvaVars = config.get(treeVarSet, 'Nominal').split(' ')
            weightVars = []
            #for systematic in systematics:
            for syst in systematics: 
                systNameUp   = syst+'_UP'   if self.config.has_option('Weights',syst+'_UP')   else syst+'_Up'
                systNameDown = syst+'_DOWN' if self.config.has_option('Weights',syst+'_DOWN') else syst+'_Down'
                if self.config.has_option('Weights',systNameUp):
                    weightVars.append(self.config.get('Weights',systNameUp))
                if self.config.has_option('Weights',systNameDown):
                    weightVars.append(self.config.get('Weights',systNameDown))

            self.trainingRegionsDict[trainingRegion] = {
                    'cut': config.get('Cuts', treeCutName),
                    'vars': mvaVars,
                    'weightVars': weightVars,
                    }

        self.TrainCut = config.get('Cuts', 'TrainCut') 
        self.EvalCut = config.get('Cuts', 'EvalCut')

        self.splitFilesChunks = splitFilesChunks
        self.chunkNumber = chunkNumber
        self.splitFilesChunkSize = splitFilesChunkSize
        
        VHbbNameSpace=config.get('VHbbNameSpace','library')
        ROOT.gSystem.Load(VHbbNameSpace)
    
    def printInfo(self):
        print ("REGION:".ljust(24),"CUT:")
        for trainingRegion,trainingRegionInfo in self.trainingRegionsDict.iteritems():
            print (" > ",trainingRegion.ljust(20), trainingRegionInfo['cut'])

    def run(self):
        # ----------------------------------------------------------------------------------------------------------------------
        # cache samples
        # ----------------------------------------------------------------------------------------------------------------------
        for sampleToCache in [self.sampleIdentifier]:
            print ('*'*80)
            print (' ',sampleToCache)
            print ('*'*80)
            # prepare caches for training and evaluation samples
            treeCaches = []
            self.sampleTree = None

            # use all (sub)samples which come from the same files (sampleIdentifier)
            subsamples = [x for x in self.samples if x.identifier == sampleToCache]

            # list of branches to keep for use as MVA input variables
            branchListOfMVAVars = BranchList()
            for sample in subsamples:
                for trainingRegion,trainingRegionInfo in self.trainingRegionsDict.iteritems():
                    for additionalCut in [self.TrainCut, self.EvalCut]:
                        branchListOfMVAVars.addCut(trainingRegionInfo['vars'])
                    for weightVar in trainingRegionInfo['weightVars']:
                        branchListOfMVAVars.addCut(weightVar)
            branchListOfMVAVars.addCut(self.config.get('Weights', 'weightF'))
            mvaBranches = branchListOfMVAVars.getListOfBranches()

            # loop over all samples
            for sample in subsamples:

                # add cuts for all training regions
                for trainingRegion,trainingRegionInfo in self.trainingRegionsDict.iteritems():

                    # add cuts for training and evaluation
                    additionalCuts = [None] if sample.isData() else [self.TrainCut, self.EvalCut]
                    for additionalCut in additionalCuts:

                        # cuts
                        sampleCuts = [sample.subcut]
                        if additionalCut:
                            sampleCuts.append(additionalCut)
                        if trainingRegionInfo['cut']:
                            sampleCuts.append(trainingRegionInfo['cut'])

                        # add cache object
                        tc = TreeCache.TreeCache(
                            name='{region}_{sample}_{tr}'.format(region=trainingRegion, sample=sample.name, tr='TRAIN' if additionalCut==self.TrainCut else 'EVAL'),
                            sample=sample.name,
                            cutList=sampleCuts,
                            inputFolder=self.samplesPath,
                            splitFilesChunks=self.splitFilesChunks,
                            chunkNumber=self.chunkNumber,
                            splitFilesChunkSize=self.splitFilesChunkSize,
                            branches=mvaBranches,
                            config=self.config,
                            debug=True
                        )

                        # check if this part of the sample is already cached
                        isCached = tc.partIsCached()
                        if not isCached or self.force:
                            if isCached:
                                tc.deleteCachedFiles(chunkNumber=self.chunkNumber)
                            # for the first sample which comes from this files, load the tree
                            if not self.sampleTree:
                                self.sampleTree = SampleTree({'name': sample.identifier, 'folder': self.samplesPath}, splitFilesChunkSize=self.splitFilesChunkSize, chunkNumber=self.chunkNumber, config=self.config, saveMemory=True)
                            treeCaches.append(tc.setSampleTree(self.sampleTree).cache())

            if len(treeCaches) > 0:
                # run on the tree
                self.sampleTree.process()
            else:
                print ("nothing to do!")

# read arguments
argv = sys.argv
parser = OptionParser()
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                          help="Verbose mode.")
parser.add_option("-C", "--config", dest="config", default=[], action="append",
                      help="configuration file")
parser.add_option("-t","--trainingRegions", dest="trainingRegions", default='',
                      help="cut region identifier")
parser.add_option("-s","--sampleIdentifier", dest="sampleIdentifier", default='',
                      help="sample identifier (no subsample!)")
parser.add_option("-n","--splitFilesChunks", dest="splitFilesChunks", default='',
                      help="number of parts")
parser.add_option("-i","--chunkNumber", dest="chunkNumber", default='',
                      help="number of part to cache")
parser.add_option("-p","--splitFilesChunkSize", dest="splitFilesChunkSize", default='',
                      help="number of files per part")
parser.add_option("-f","--force", action="store_true", dest="force", default=False,
                      help="force overwriting of already cached files")
parser.add_option("-l","--fileList", dest="fileList", default="",
                      help="file list")
(opts, args) = parser.parse_args(argv)
if opts.config =="":
        opts.config = "config"

# Import after configure to get help message
from myutils import BetterConfigParser, ParseInfo

# load config
config = BetterConfigParser()
config.read(opts.config)

# initialize
trainingRegions = opts.trainingRegions.split(',')
splitFilesChunks = int(opts.splitFilesChunks) if len(opts.splitFilesChunks) > 0 else 1
chunkNumber = int(opts.chunkNumber) if len(opts.chunkNumber) > 0 else 1
splitFilesChunkSize = int(opts.splitFilesChunkSize) if len(opts.splitFilesChunkSize) > 0 else -1
ct = CacheTraining(config=config, sampleIdentifier=opts.sampleIdentifier, trainingRegions=trainingRegions, chunkNumber=chunkNumber, splitFilesChunks=splitFilesChunks, splitFilesChunkSize=splitFilesChunkSize, force=opts.force)
ct.printInfo()

# run training
ct.run()
