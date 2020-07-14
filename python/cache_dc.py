#!/usr/bin/env python
from __future__ import print_function

import sys, ROOT, warnings
ROOT.gROOT.SetBatch(True)
ROOT.v5.TFormula.SetMaxima(10000)
#suppres the EvalInstace conversion warning bug
warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='creating converter.*' )
from optparse import OptionParser
import json
from myutils import NewTreeCache as TreeCache
from myutils.sampleTree import SampleTree as SampleTree
from myutils.Datacard import Datacard
from myutils.BranchList import BranchList
from myutils.FileList import FileList
import resource
import os

# ------------------------------------------------------------------------------
# script to produce cached tree for datacards
# ------------------------------------------------------------------------------
class CacheDatacards(object):

    def __init__(self, config, regions, sampleToCache, splitFilesChunks=1, chunkNumber=1, splitFilesChunkSize=-1, forceRedo=False, fileList=None, verbose=False):
        self.verbose = verbose or ('XBBDEBUG' in os.environ)
        self.config = config
        self.regions = regions
        self.treeCaches = []
        self.sampleTree = None
        self.sampleToCache = sampleToCache
        self.forceRedo = forceRedo
        
        # settings which part of input files to process
        self.splitFilesChunkSize = splitFilesChunkSize
        self.splitFilesChunks = splitFilesChunks
        self.chunkNumber = chunkNumber
        self.fileList = FileList.decompress(fileList) if fileList else None

        # initialize Datacard objects
        self.dcMakers = [Datacard(config=self.config, region=region) for region in self.regions]
    
    # make a minimum list of samples which is needed to produce all the Datacard regions at the same time
    def getAllSamples(self):
        samples = []
        for dcMaker in self.dcMakers:
            for sample in dcMaker.getAllSamples():
                if len([x for x in samples if x.name == sample.name]) < 1:
                    samples.append(sample)
        return samples

    def prepare(self):
        if len(self.dcMakers) > 0:
            self.treeCaches = []
            self.sampleTree = None

            # cuts
            allSamples = self.getAllSamples() 
            subsamples = [x for x in allSamples if x.identifier == self.sampleToCache]

            # loop over all datacard regions
            for dcMaker in self.dcMakers:

                # loop over all subsamples (which come from the same root tree files)
                for sample in subsamples:

                    # combine subcut and systematics cut with logical AND
                    # systematics cuts are combined with logical OR, such that 1 cache file can be used for all the systematics

                    isData = (sample.type == 'DATA')
                    systematicsCuts = sorted(list(set([x['cachecut'] for x in dcMaker.getSystematicsList(isData=isData)])))
                    sampleCuts = {'AND': [sample.subcut, {'OR': systematicsCuts}]}
                    if self.verbose:
                        print (json.dumps(sampleCuts, sort_keys=True, indent=8, default=str))

                    # make list of branches to keep in root file
                    branchList = BranchList(sample.subcut)
                    branchList.addCut([x['cachecut'] for x in dcMaker.getSystematicsList()])
                    branchList.addCut([x['cut'] for x in dcMaker.getSystematicsList()])
                    branchList.addCut([x['var'] for x in dcMaker.getSystematicsList()])
                    branchList.addCut([x['weight'] for x in dcMaker.getSystematicsList()])
                    branchList.addCut(self.config.get('Weights', 'weightF'))
                    branchList.addCut(eval(self.config.get('Branches', 'keep_branches')))
                    branchesToKeep = branchList.getListOfBranches()

                    # arbitrary (optional) name for the output tree, used for print-out (the TreeCache object has no idea what it is doing, e.g. dc, plot etc.)
                    cacheName = 'dc:{region}_{sample}'.format(region=dcMaker.getRegion(), sample=sample.name) 
                    
                    # add cache object
                    tc = TreeCache.TreeCache(
                        name=cacheName,
                        sample=sample.name,
                        cutList=sampleCuts,
                        cutSequenceMode='TREE',
                        branches=branchesToKeep,
                        inputFolder=dcMaker.path,
                        splitFilesChunks=self.splitFilesChunks,
                        chunkNumber=self.chunkNumber,
                        splitFilesChunkSize=self.splitFilesChunkSize,
                        fileList=self.fileList,
                        config=self.config,
                        debug=self.verbose
                    )

                    # check if this part of the sample is already cached
                    isCached = tc.partIsCached() 
                    print ("check if sample \x1b[34m{sample}\x1b[0m part {part} is cached:".format(sample=sample.name, part=self.chunkNumber), isCached)
                    if not isCached or self.forceRedo:
                        if isCached:
                            tc.deleteCachedFiles(chunkNumber=self.chunkNumber)

                        # for the first sample which comes from this files, load the tree
                        if not self.sampleTree:
                            self.sampleTree = SampleTree({'name': sample.identifier, 'folder': dcMaker.path}, splitFilesChunkSize=self.splitFilesChunkSize, chunkNumber=self.chunkNumber, config=self.config, saveMemory=True)
                            if not self.sampleTree or not self.sampleTree.tree:
                                print ("\x1b[31mERROR: creation of sample tree failed!!\x1b[0m")
                                raise Exception("CreationOfSampleTreeFailed")

                            # consistency check on the file list at submission time and now
                            fileListNow = self.sampleTree.getSampleFileNameChunk(self.chunkNumber)
                            if self.fileList and (sorted(self.fileList) != sorted(fileListNow)):
                                print ("\x1b[31mERROR: sample files have changed between submission and run of the job!\x1b[0m")
                                raise Exception("SampleFilesHaveChanged")

                        # connect the TreeCache object to the input sampleTree and add it to the list of cached trees 
                        self.treeCaches.append(tc.setSampleTree(self.sampleTree).cache())
        else:
            print("WARNING: no datacard regions added, nothing to do.")
        return self

    def run(self):
        if len(self.treeCaches) > 0:
            # run on the tree
            self.sampleTree.process()
        else:
            print ("nothing to do!")


if __name__ == "__main__":
    # read arguments
    argv = sys.argv
    parser = OptionParser()
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                              help="Verbose mode.")
    parser.add_option("-C", "--config", dest="config", default=[], action="append",
                          help="configuration file")
    parser.add_option("-t", "--regions", dest="regions", default='',
                          help="cut regions identifier")
    parser.add_option("-s", "--sampleIdentifier", dest="sampleIdentifier", default='',
                          help="sample identifier (no subsample!)")
    parser.add_option("-n", "--splitFilesChunks", dest="splitFilesChunks", default='',
                          help="number of chunks the cached file is split into")
    parser.add_option("-i", "--chunkNumber", dest="chunkNumber", default='',
                          help="number of part to cache")
    parser.add_option("-p", "--splitFilesChunkSize", dest="splitFilesChunkSize", default='',
                          help="number of files per part")
    parser.add_option("-f", "--force", action="store_true", dest="force", default=False,
                          help="force overwriting of already cached files")
    parser.add_option("-l", "--fileList", dest="fileList", default="",
                          help="file list")
    (opts, args) = parser.parse_args(argv)
    if opts.config == "":
            opts.config = "config"

    # Import after configure to get help message
    from myutils import BetterConfigParser, mvainfo, ParseInfo

    # load config
    config = BetterConfigParser()
    config.read(opts.config)

    # if no region is given in argument, run it for all of them
    regionsListString = opts.regions if len(opts.regions.strip())>0 else config.get('LimitGeneral', 'List')
    regions = [x.strip() for x in regionsListString.split(',') if len(x.strip()) > 0]
    
    # init and run
    cacheDC = CacheDatacards(
            config=config, 
            regions=regions, 
            sampleToCache=opts.sampleIdentifier, 
            chunkNumber=int(opts.chunkNumber) if len(opts.chunkNumber) > 0 else 1, 
            splitFilesChunks=int(opts.splitFilesChunks) if len(opts.splitFilesChunks) > 0 else 1,
            splitFilesChunkSize=int(opts.splitFilesChunkSize) if len(opts.splitFilesChunkSize) > 0 else -1, 
            forceRedo=opts.force, 
            fileList=opts.fileList, 
            verbose=opts.verbose) 

    if cacheDC.prepare():
        cacheDC.run()
    else:
        raise Exception("CacheDCinitFailed")
