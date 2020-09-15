#!/usr/bin/env python
import ROOT
import sys
ROOT.gROOT.SetBatch(True)
from myutils import BetterConfigParser, ParseInfo
import glob
from multiprocessing import Process
from myutils.sampleTree import SampleTree as SampleTree
import argparse
import itertools
from myutils.BranchList import BranchList
from myutils.XbbConfig import XbbConfigReader, XbbConfigTools

parser = argparse.ArgumentParser(description='getExtWeightsNew: compute stitching coefficients for samples with overlap')
parser.add_argument('--samples', dest='samples', help='list of sample identifiers')
parser.add_argument('--cuts', dest='cuts', help='cuts', default='')
parser.add_argument('--from', dest='fromFolder', help='folder name to use, e.g. PREPout', default='SYSout')
parser.add_argument('--fc', dest='fc', help='fc', default='')
parser.add_argument('--prune', dest='prune', help='remove contributions with fraction lower than this', default='0.0')
parser.add_argument('-T', dest='tag', help='config tag')
args = parser.parse_args()

# 1) ./getExtWeights.py configname
# 2) ./getExtWeight.py configname verify

def getEventCount(config, sampleIdentifier, cut="1", sampleTree=None, sample=None):
    if not sampleTree:
        sampleTree = SampleTree({
                'name': sampleIdentifier, 
                'folder': config.get('Directories',args.fromFolder).strip()
            }, config=config)
    h1 = ROOT.TH1D("h1", "h1", 1, 0, 2)
    scaleToXs = sampleTree.getScale(sample)
    #nEvents = sampleTree.tree.Draw("1>>h1", "(" + cut + ")*genWeight*%1.6f"%scaleToXs, "goff")
    nEvents = sampleTree.tree.Draw("1>>h1", cut, "goff")
    nEventsWeighted = h1.GetBinContent(1)
    #print("DEBUG:", sampleIdentifier, cut, " MC events:", nEvents, " (weighted:", nEventsWeighted, ")")
    h1.Delete()
    return nEvents


# load config
config = XbbConfigReader.read(args.tag)
sampleInfo = ParseInfo(samples_path=config.get('Directories', args.fromFolder), config=config)
mcSamples = sampleInfo.get_samples(XbbConfigTools(config).getMC())

pruneThreshold = float(args.prune)

sampleGroups = []
for x in args.samples.split(','):
    sampleGroups.append(x.split('+'))

sampleCuts = args.cuts.strip().split(',')
if args.fc != '':
    cutGroups = [x.strip(',').split(',') for x in args.fc.strip(';').split(';')]
    # cartesian product
    sampleCuts = list(set(['&&'.join(x) for x in itertools.product(*cutGroups)]))
    print sampleCuts
    print "=> len=",len(sampleCuts)

countDict = {}
for sampleGroup in sampleGroups:
    count = 0
    for sampleIdentifier in sampleGroup:
        print "\x1b[32m",sampleIdentifier,"\x1b[0m"
        countDict[sampleIdentifier] = {}

        samples_matching = [x for x in mcSamples if x.identifier == sampleIdentifier]
        if len(samples_matching) > 0:
            sample = samples_matching[0]

            sampleTree = SampleTree({'sample': sample, 'folder': config.get('Directories',args.fromFolder).strip()}, config=config)
            print "CUT=",sampleCuts,":"
            for sampleCut in sampleCuts:
                sampleCount = getEventCount(config, sampleIdentifier, sampleCut, sampleTree=sampleTree, sample=sample)
                print sampleIdentifier,sampleCut,"\x1b[34m=>",sampleCount,"\x1b[0m"
                if sampleCut in countDict[sampleIdentifier]:
                    print "duplicate!!", sampleIdentifier, sampleCut, countDict[sampleIdentifier][sampleCut]
                    raise Exception("duplicate")
                countDict[sampleIdentifier][sampleCut] = sampleCount
        else:
            print sampleIdentifier
            raise Exception("SampleMissing")

# pruning
for sampleIdentifier,counts in countDict.iteritems():
    specialweight = ''
    for cut,cutCount in counts.iteritems():
        totalCount = 0
        for sample,sampleCounts in countDict.iteritems():
            if cut in sampleCounts:
                totalCount += sampleCounts[cut]

        if cutCount>0 and totalCount>0:
            fraction = 1.0*cutCount/totalCount
            if fraction < pruneThreshold:
                print "\x1b[31mfor "+sampleIdentifier + " set count("+cut+") from %d"%counts[cut] + " to 0\x1b[0m"
                counts[cut] = 0

print "result:"
print countDict
print "-"*80

for sampleIdentifier,counts in countDict.iteritems():
    specialweight = ''
    for cut,cutCount in counts.iteritems():
        totalCount = 0
        for sample,sampleCounts in countDict.iteritems():
            if cut in sampleCounts:
                totalCount += sampleCounts[cut]

        if cutCount>0:
            if cutCount==totalCount:
                if len(specialweight) > 0:
                    specialweight += ' + '
                specialweight += '(' + cut + ')'
            else:
                if len(specialweight) > 0:
                    specialweight += ' + '
                specialweight += '(' + cut + ')*%1.5f'%(1.0*cutCount/totalCount)
    print sampleIdentifier, specialweight

print "-"*80

with open("./stitching.csv", "w") as of:
    for sampleIdentifier,counts in countDict.iteritems():
        for cut,cutCount in counts.iteritems():
            of.write("{s},{c},{v}\n".format(s=sampleIdentifier,c=cut,v=cutCount))

