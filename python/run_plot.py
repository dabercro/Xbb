#!/usr/bin/env python
from __future__ import print_function
from optparse import OptionParser
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.v5.TFormula.SetMaxima(10000)

from myutils import NewTreeCache as TreeCache
from myutils.sampleTree import SampleTree as SampleTree
from myutils import BetterConfigParser, ParseInfo
from myutils.NewStackMaker import NewStackMaker as StackMaker
from myutils.samplesclass import Sample
from myutils.FileLocator import FileLocator
from myutils.XbbTools import XbbTools
import os,sys

class PlotHelper(object):

    def __init__(self, config, region, vars = None, title=None, sampleIdentifier=None):
        self.debug = 'XBBDEBUG' in os.environ
        self.config = config
        self.region = region
        self.vars = vars
        self.title = title if title and len(title)>0 else None
        self.sampleIdentifiers = sampleIdentifier.split(',') if sampleIdentifier and len(sampleIdentifier) > 0 else None

        # VHbb namespace
        VHbbNameSpace=config.get('VHbbNameSpace','library')
        returnCode = ROOT.gSystem.Load(VHbbNameSpace)
        if returnCode != 0:
            print ("\x1b[31mERROR: loading VHbbNameSpace failed with code %d\x1b[0m"%returnCode)
        else:
            print ("INFO: loaded VHbbNameSpace: %s"%VHbbNameSpace)

        # input/output paths
        self.samplesPath = config.get('Directories', 'plottingSamples')
        self.samplesInfo = ParseInfo(samples_path=self.samplesPath, config=self.config)
        self.sampleFilesFolder = config.get('Directories', 'samplefiles')
        self.plotPath = config.get('Directories', 'plotpath')

        # plot regions
        self.configSection='Plot:%s'%region
        self.dataOverBackground = self.config.has_option('Plot_general', 'plotDataOverBackground') and eval(self.config.get('Plot_general', 'plotDataOverBackground'))

        # variables
        if self.vars and type(self.vars) == list:
            self.vars = [x.strip() for x in self.vars if len(x.strip()) > 0] 

        # if variables not specified in command line, read from config
        if not self.vars or len(self.vars) < 1:
            varListFromConfig = self.config.get(self.configSection, 'vars').split(',')
            print ("VARS::", self.configSection, " => ", varListFromConfig)
            self.vars = [x.strip() for x in varListFromConfig if len(x.strip()) > 0]

        # resolve plot variables (find plot section name if ROOT expression is given)
        self.vars = [XbbTools.resolvePlotVariable(var, self.config) for var in vars]

        # additional cut to only plot a subset of the region
        self.subcut = None
        if self.config.has_option(self.configSection, 'subcut'):
            self.subcut = self.config.get(self.configSection, 'subcut')
            print("INFO: use cut:", self.subcut)

        # additional global blinding cut:
        self.addBlindingCut = None
        if self.config.has_option('Plot_general','addBlindingCut'): #contained in plots, cut on the event number
            self.addBlindingCut = self.config.get('Plot_general','addBlindingCut')
            print ('adding add. blinding cut:', self.addBlindingCut)

        # load samples
        if self.config.has_section(self.configSection):
            # read data from region definition
            if self.config.has_option(self.configSection, 'Datas'):
                self.data = eval(self.config.get(self.configSection, 'Datas')) # read the data corresponding to each CR (section)
            elif self.config.has_option(self.configSection, 'Data'):
                self.data = eval(self.config.get(self.configSection, 'Data')) # read the data corresponding to each CR (section)
            else:
                self.data = eval(self.config.get('Plot_general', 'Data'))
        else:
            # use default datasets
            self.data = eval(self.config.get('Plot_general', 'Data'))
        self.mc = eval(self.config.get('Plot_general', 'samples')) # read the list of mc samples
        self.total_lumi = eval(self.config.get('General', 'lumi'))
        self.signalRegion = False

        self.dataSamples = self.samplesInfo.get_samples(self.data)
        self.mcSamples = self.samplesInfo.get_samples(self.mc)

        # filter samples used in the plot
        if self.sampleIdentifiers:
            self.dataSamples = [x for x in self.dataSamples if x.identifier in self.sampleIdentifiers]
            self.mcSamples =   [x for x in self.mcSamples   if x.identifier in self.sampleIdentifiers]

        self.groupDict = eval(self.config.get('Plot_general', 'Group'))
        self.subcutPlotName = ''
        self.histogramStacks = {}

    def getSampleGroup(self, sample):
        if sample.name in self.groupDict:
            if self.debug and sample.group != self.groupDict[sample.name]:
                print("\x1b[41m\x1b[33mDEBUG: overwrite group from", sample.group, " to ", self.groupDict[sample.name], "\x1b[0m")
            return self.groupDict[sample.name]
        else:
            return sample.group

    def prepare(self):
        print ("INFO: starting plot for region \x1b[34m{region}\x1b[0m, variables:".format(region=region))
        for var in self.vars:
            print ("  > {var}".format(var=var))

        self.histogramStacks = {}
        for var in self.vars:
            self.histogramStacks[var] = StackMaker(self.config, var, self.region, self.signalRegion, None, '_'+self.subcutPlotName, title=self.title)

        fileLocator = FileLocator(config=self.config, useDirectoryListingCache=True)

        # add DATA + MC samples
        for sample in self.dataSamples + self.mcSamples:
            
            # cuts
            sampleCuts = [sample.subcut]
            if self.config.has_option('Cuts', self.region):
                sampleCuts.append(self.config.get('Cuts', self.region))
            if self.config.has_option(self.configSection, 'Datacut'):
                sampleCuts.append(self.config.get(self.configSection, 'Datacut'))
            if self.addBlindingCut:
                sampleCuts.append(self.addBlindingCut)
            
            # get sample tree from cache
            tc = TreeCache.TreeCache(
                    sample=sample,
                    cutList=sampleCuts,
                    inputFolder=self.samplesPath,
                    config=config,
                    fileLocator=fileLocator
                )
            sampleTree = tc.getTree()

            if sampleTree:
                groupName = self.getSampleGroup(sample) 
                print (" > found the tree, #entries = ", sampleTree.tree.GetEntries())
                print ("   > group =", groupName)
                print (" > now adding the tree for vars=", self.vars)
                
                # add the sample tree for all the variables
                for var in self.vars:
                    self.histogramStacks[var].addSampleTree(sample=sample, sampleTree=sampleTree, groupName=groupName, cut=self.subcut if self.subcut else '1')
            else:
                print ("\x1b[31mERROR: sampleTree not available for ", sample,", run caching again!!\x1b[0m")
                raise Exception("CachedTreeMissing")
        return self

    def run(self):
        # draw
        for var in self.vars:
            self.histogramStacks[var].Draw(outputFolder=self.plotPath, prefix='{region}__{var}_'.format(region=self.region, var=var), dataOverBackground=self.dataOverBackground)
            if self.config.has_option('Plot_general', 'drawNormalizedPlots') and eval(self.config.get('Plot_general', 'drawNormalizedPlots')):
                self.histogramStacks[var].Draw(outputFolder=self.plotPath, prefix='comp_{region}__{var}_'.format(region=self.region, var=var), normalize=True)
        return self

    def getHistogramStack(self, var):
        if var in self.vars and var in self.histogramStacks:
            return self.histogramStacks[var]
        else:
            return None

if __name__ == "__main__":
    # read arguments
    argv = sys.argv
    parser = OptionParser()
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                              help="Verbose mode.")
    parser.add_option("-C", "--config", dest="config", default=[], action="append",
                          help="configuration file")
    parser.add_option("-r","--regions", dest="regions", default='',
                          help="cut region identifiers, separated by comma")
    parser.add_option("-p","--vars", dest="vars", default='',
                          help="plot variables, separated by comma")
    parser.add_option("-t","--title", dest="title", default='',
                          help="plot title")
    parser.add_option("-s","--sampleIdentifier", dest="sampleIdentifier", default='',
                                      help="sample identifier (no subsample!)")
    (opts, args) = parser.parse_args(argv)
    if opts.config == "":
            opts.config = ["config"]

    # Import after configure to get help message
    from myutils import BetterConfigParser, mvainfo, ParseInfo

    # load config
    config = BetterConfigParser()
    config.read(opts.config)

    # run plotter
    regions = opts.regions.split(',')
    vars = opts.vars.split(',')
    for region in regions:
        plotter = PlotHelper(config=config, region=region, vars=vars, title=opts.title, sampleIdentifier=opts.sampleIdentifier)
        plotter.prepare()
        plotter.run()

