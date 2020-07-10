#! /usr/bin/env python
from __future__ import print_function
from optparse import OptionParser
import sys
import os
from myutils.XbbConfig import XbbConfigReader, XbbConfigTools
from myutils import ParseInfo

argv = sys.argv
parser = OptionParser()
parser.add_option("-T","--tag", dest="tag", default='', help="config tag")
parser.add_option("-R","--regions", dest="regions", default='', help="region")
(opts, args) = parser.parse_args(argv)

config = XbbConfigTools(config=XbbConfigReader.read(opts.tag))

if len(opts.regions) > 0:
    dnnRegions = opts.regions.split(',')
else:
    dnnRegions = config.get('MVALists', 'List_for_submitscript', '').strip().split(',')

plotVarList = []
for dnnRegion in dnnRegions:
    print("-"*10, dnnRegion, "-"*80)
    dnnVarset = config.get(dnnRegion, 'treeVarSet') 
    dnnVars   = config.get(dnnVarset, 'Nominal').strip().split(' ') 
    for dnnVar in dnnVars:
        found = False
        for section in config.getPlotVariableSections():
            plotVar = config.get(section, 'relPath').strip().replace(' ', '')
            if plotVar == dnnVar or ('('+plotVar+')' == dnnVar):
                print("match:", dnnVar, " --> ", section.split(':')[1])
                plotVarList.append(section.split(':')[1].strip())
                found = True
            if found:
                break
        if not found:
            print("\x1b[31mmissing:", dnnVar, "\x1b[0m")
    print("\ncommand for plotting:")
    print("--vars", ','.join(sorted(list(set(plotVarList)))))


