#! /usr/bin/env python
from __future__ import print_function
import ROOT
ROOT.gROOT.SetBatch(True)
from myutils.XbbConfig import XbbConfigReader, XbbConfigTools
from myutils.sampleTree import SampleTree as SampleTree
from myutils.BranchList import BranchList

config     = XbbConfigTools(XbbConfigReader.read("Zvv2017"))
sampleTree = SampleTree({'name': 'MET', 'folder': config.get('Directories', 'dcSamples')}, config=config)
variables  = ["H_pt","MET_Pt","H_pt/MET_Pt"]

# enable only explicitly used branches
sampleTree.enableBranches(BranchList(variables).getListOfBranches())

# create TTReeFormula's
for variable in variables:
    sampleTree.addFormula(variable)

# loop over events
for event in sampleTree:
    print(sampleTree.tree.GetReadEntry(), ", ".join([x + "=%1.4f"%sampleTree.evaluate(x) for x in variables])) 
    if sampleTree.tree.GetReadEntry() > 98:
        break
