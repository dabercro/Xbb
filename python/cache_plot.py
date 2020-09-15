#!/usr/bin/env python
from __future__ import print_function
from optparse import OptionParser
import ROOT
ROOT.gROOT.SetBatch(True)

from myutils import NewTreeCache as TreeCache
from myutils.sampleTree import SampleTree as SampleTree
from myutils import BetterConfigParser, ParseInfo
from myutils.FileList import FileList
from myutils.BranchList import BranchList
import os,sys

class CachePlot(object):

    def __init__(self, config, sampleIdentifier, regions, splitFilesChunks=1, chunkNumber=1, splitFilesChunkSize=-1, forceRedo=False, fileList=None):
        self.config = config
        self.sampleIdentifier = sampleIdentifier
        self.regions = list(set(regions))
        self.forceRedo = forceRedo

        self.sampleTree = None
        self.samplesPath = self.config.get('Directories', 'plottingSamples')
        self.samplesInfo = ParseInfo(samples_path=self.samplesPath, config=self.config) 
        self.sampleFilesFolder = self.config.get('Directories', 'samplefiles')

        self.sampleNames = list(eval(self.config.get('Plot_general', 'samples')))
        self.dataNames = list(eval(self.config.get('Plot_general', 'Data')))
        self.samples = self.samplesInfo.get_samples(self.sampleNames + self.dataNames)

        self.regionsDict = {}
        for region in self.regions:
            treeCut = config.get('Cuts', region)
            self.regionsDict[region] = {'cut': treeCut}
        self.splitFilesChunkSize = splitFilesChunkSize
        self.splitFilesChunks = splitFilesChunks
        self.chunkNumber = chunkNumber
        self.fileList = FileList.decompress(fileList) if fileList else None
    
        VHbbNameSpace=config.get('VHbbNameSpace','library')
        returnCode = ROOT.gSystem.Load(VHbbNameSpace)
        if returnCode != 0:
            print ("\x1b[31mERROR: loading VHbbNameSpace failed with code %d\x1b[0m"%returnCode)
        else:
            print ("INFO: loaded VHbbNameSpace: %s"%VHbbNameSpace)

    def printInfo(self):
        print ("REGION:".ljust(24),"CUT:")
        for region,regionInfo in self.regionsDict.iteritems():
            print (" > ",region.ljust(20), regionInfo['cut'])

    def run(self):

        # keep additional branches for plotting
        try:
            keepBranchesPlot = eval(self.config.get('Branches', 'keep_branches_plot'))
        except:
            keepBranchesPlot = []
        try:
            keepBranchesPlot += eval(self.config.get('Branches', 'keep_branches'))
        except:
            pass

        # also keep some branches which might be used later in variables definition and weights
        try:
            for section in self.config.sections():
                try:
                    if section.startswith('plotDef:') and self.config.has_option(section, 'relPath'):
                        keepBranchesPlot.append(self.config.get(section, 'relPath'))
                except Exception as e:
                    print("\x1b[31mWARNING: config error in:", section, "=>", e, "\x1b[0m")
        except Exception as e2:
            print("\x1b[31mERROR: config file contains an error! automatic selection of branches to keep will not work!\x1b[0m")
            print(e2)
        try:
            keepBranchesPlot.append(self.config.get('Weights', 'weightF'))
        except:
            pass
        # plotting region cut
        for region,regionInfo in self.regionsDict.iteritems():
            keepBranchesPlot.append(regionInfo['cut'])
        keepBranchesPlotFinal = BranchList(keepBranchesPlot).getListOfBranches()
        print("KEEP:", keepBranchesPlotFinal)


        # ----------------------------------------------------------------------------------------------------------------------
        # cache samples
        # ----------------------------------------------------------------------------------------------------------------------
        for sampleToCache in [self.sampleIdentifier]:
            print ('*'*80)
            print (' ',sampleToCache)
            print ('*'*80)
            # prepare caches for training and evaluation samples
            treeCaches = []
            sampleTree = None

            # for all (sub)samples which come from the same files (sampleIdentifier)
            subsamples = [x for x in self.samples if x.identifier == sampleToCache]
            for sample in subsamples:

                # add cuts for all training regions
                for region,regionInfo in self.regionsDict.iteritems():

                    configSection = 'Plot:%s'%region
                    
                    # cuts
                    sampleCuts = [sample.subcut]
                    if regionInfo['cut']:
                        sampleCuts.append(regionInfo['cut'])
                    if self.config.has_option(configSection, 'Datacut'):
                        sampleCuts.append(self.config.get(configSection, 'Datacut'))
                    if self.config.has_option('Plot_general','addBlindingCut'):
                        sampleCuts.append(self.config.has_option('Plot_general', 'addBlindingCut'))

                    # arbitrary (optional) name for the output tree, used for print-out (the TreeCache object has no idea what it is doing, e.g. dc, plot etc.)
                    cacheName = 'plot:{region}_{sample}'.format(region=region, sample=sample.name)

                    # add cache object
                    tc = TreeCache.TreeCache(
                        name=cacheName,
                        sample=sample.name,
                        cutList=sampleCuts,
                        inputFolder=self.samplesPath,
                        splitFilesChunks=self.splitFilesChunks,
                        chunkNumber=self.chunkNumber,
                        splitFilesChunkSize=self.splitFilesChunkSize,
                        fileList=self.fileList,
                        branches=keepBranchesPlotFinal,
                        config=self.config,
                        debug=True
                    )

                    # check if this part of the sample is already cached
                    isCached = tc.partIsCached()
                    if not isCached or self.forceRedo:
                        if isCached:
                            tc.deleteCachedFiles(chunkNumber=self.chunkNumber)

                        # for the first sample which comes from this files, load the tree
                        if not self.sampleTree:
                            self.sampleTree = SampleTree({'name': sample.identifier, 'folder': self.samplesPath}, splitFilesChunkSize=self.splitFilesChunkSize, chunkNumber=self.chunkNumber, config=self.config, saveMemory=True)
                            if not self.sampleTree or not self.sampleTree.tree:
                                print("\x1b[31mERROR: creation of sample tree failed!!\x1b[0m")
                                raise Exception("CreationOfSampleTreeFailed")
                            # consistency check on the file list at submission time and now
                            fileListNow = self.sampleTree.getSampleFileNameChunk(self.chunkNumber)
                            if self.fileList and (sorted(self.fileList) != sorted(fileListNow)):
                                print("\x1b[31mERROR: sample files have changed between submission and run of the job!\x1b[0m")
                                raise Exception("SampleFilesHaveChanged")

                        treeCaches.append(tc.setSampleTree(self.sampleTree).cache())
                    else:
                        print ("INFO: already cached!",tc, "(",tc.hash,")")

            if len(treeCaches) > 0:
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
    parser.add_option("-t","--regions", dest="regions", default='',
                          help="cut regions identifier")
    parser.add_option("-s","--sampleIdentifier", dest="sampleIdentifier", default='',
                          help="sample identifier (no subsample!)")
    parser.add_option("-n","--splitFilesChunks", dest="splitFilesChunks", default='',
                          help="number of chunks the cached file is split into")
    parser.add_option("-i","--chunkNumber", dest="chunkNumber", default='',
                          help="number of part to cache")
    parser.add_option("-p","--splitFilesChunkSize", dest="splitFilesChunkSize", default='',
                          help="number of files per part")
    parser.add_option("-f","--force", action="store_true", dest="force", default=False,
                          help="force overwriting of already cached files")
    parser.add_option("-l","--fileList", dest="fileList", default="",
                          help="file list")
    (opts, args) = parser.parse_args(argv)
    if opts.config == "":
            opts.config = "config"

    # Import after configure to get help message
    from myutils import BetterConfigParser, mvainfo, ParseInfo

    # load config
    config = BetterConfigParser()
    vhbbPlotDef = opts.config[0].split('/')[0]+'/vhbbPlotDef.ini'
    opts.config.append(vhbbPlotDef)
    config.read(opts.config)

    # initialize
    regions = opts.regions.split(',')
    splitFilesChunks = int(opts.splitFilesChunks) if len(opts.splitFilesChunks) > 0 else 1
    chunkNumber = int(opts.chunkNumber) if len(opts.chunkNumber) > 0 else 1
    splitFilesChunkSize = int(opts.splitFilesChunkSize) if len(opts.splitFilesChunkSize) > 0 else -1
    ct = CachePlot(config=config, sampleIdentifier=opts.sampleIdentifier, regions=regions, chunkNumber=chunkNumber, splitFilesChunks=splitFilesChunks, splitFilesChunkSize=splitFilesChunkSize, forceRedo=opts.force, fileList=opts.fileList)
    ct.printInfo()

    # run training
    ct.run()

