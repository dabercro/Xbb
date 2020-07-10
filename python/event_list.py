#! /usr/bin/env python
from __future__ import print_function
import ROOT
ROOT.gROOT.SetBatch(True)
from optparse import OptionParser
from myutils.XbbConfig import XbbConfigReader, XbbConfigTools
from myutils.sampleTree import SampleTree as SampleTree
from myutils.BranchList import BranchList
import sys


if __name__ == "__main__":
    # read arguments
    argv = sys.argv
    parser = OptionParser()
    parser.add_option("-T", "--tag", dest="tag", help="config tag")
    parser.add_option("-R", "--regions", dest="regions", default="", help="regions")
    (opts, args) = parser.parse_args(argv)

    config      = XbbConfigTools(XbbConfigReader.read(opts.tag))
    inputFolder = config.get('Directories', 'dcSamples')
    logFolder   = config.get('Directories', 'tagDir')
    
    config.loadNamespaces()

    regions = config.getDatacardRegions() if len(opts.regions) < 1 else config.parseCommaSeparatedList(opts.regions) 
    for region in regions:
        dataSamples = eval(config.get('dc:'+region, 'data'))
        for dataSample in dataSamples:
            sampleTree = SampleTree({'name': dataSample, 'folder': inputFolder}, config=config)
           
            outputFileName = logFolder + '/' + region + '_' + dataSample + '.txt'
            print("save event list to:", outputFileName)

            treePlayer = sampleTree.tree.GetPlayer()
            treePlayer.SetScanRedirect(True)
            treePlayer.SetScanFileName(outputFileName)

            branchList = BranchList(["run","event"])
            regionCut  = config.get('Cuts', config.get('dc:'+region, 'cut') if config.has_option('dc:'+region, 'cut') else region)
            branchList.addCut(regionCut)
            sampleTree.enableBranches(branchList.getListOfBranches())

            sampleTree.tree.Scan("run:event", regionCut, "colsize=16")


