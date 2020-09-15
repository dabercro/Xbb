#!/usr/bin/env python
from __future__ import print_function
import ROOT
import numpy as np
import array
import os
import fnmatch
from BranchTools import Collection
from BranchTools import AddCollectionsModule


class DoubleBTagWeightsFromCSV(AddCollectionsModule):

    def __init__(self, year, branchName="bTagWeightDoubleB", fileName="data/btag/deepak8v2_bbvslight.csv"):
        super(DoubleBTagWeightsFromCSV, self).__init__()
        self.version    = 3
        self.year       = year
        self.branchName = branchName
        self.fileName   = fileName
        self.CSVformat  = [('year','<i4'), ('wp','S5'), ('ptmin', 'S8'), ('ptmax', 'S8'), ('sf', 'f8'), ('err_low', 'f8'), ('err_high', 'f8')]
        self.bins       = {'mp': lambda x: (x>=0.8 and x < 0.97), 'hp': lambda x: (x>=0.97)}
        self.SF         = None
        self.jetPtName  = "FatJet_Pt"
        self.taggerName = "FatJet_deepTagMD_bbvsLight"

        self.applyToSamples         = ['DYJets*','ZJets*','WJets*','DYBJets*','ZBJets*','WBJets*','ZH*','WH*','WplusH*','WminusH*','ZZ*','WZ*','WW*']
        self.applyToHeavyFlavorOnly = ['DYJets*','ZJets*','WJets*','DYBJets*','ZBJets*','WBJets*','ZZ*','WZ*','WW*'] 

        self.isHeavyFlavorProcess = "(nGenBpt25eta2p6>0)" 

    def read_pt(self, x):
        if x.strip().lower() == 'inf':
            return 999999.9
        else:
            return float(x)

    def get_jet_sf(self, pt, score, variation=0):
        if self.SF is None:
            return 1.0
        wpBinMatches = [k for k,v in self.bins.items() if v(score)]
        if len(wpBinMatches) > 1:
            raise Exception("DoubleBscoreBinsNotExclusive")
        elif len(wpBinMatches) == 1:
            ptBinMatches = [ptBin[3:] for ptBin in self.ptWpBins if ptBin[0] == wpBinMatches[0] and pt >= ptBin[1] and pt < ptBin[2] ]
            if len(ptBinMatches) > 1:
                raise Exception("DoubleBptBinsNotExclusive")
            elif len(ptBinMatches) == 1:
                if variation < 0:
                    return ptBinMatches[0][0] + variation*ptBinMatches[0][1]
                elif variation > 0:
                    return ptBinMatches[0][0] + variation*ptBinMatches[0][2]
                else:
                    return ptBinMatches[0][0]
            else:
                return 1.0
        else:
            return 1.0


    def customInit(self, initVars):
        self.sample = initVars['sample']
        self.sampleTree = initVars['sampleTree']
        self.isData = self.sample.type == 'DATA'
        self.applyForThisSample = False
        self.checkHeavyFlavor = False
        if not self.isData:
            self.SF = [x for x in np.genfromtxt(self.fileName,  dtype=self.CSVformat, delimiter=None, skip_header=1) if int(x[0]) == self.year]

            self.ptWpBins = [ [x[1], self.read_pt(x[2]), self.read_pt(x[3]), x[4], x[5], x[6], x[1] + "_pt" + x[2] + "to" + x[3]] for x in self.SF]

            self.addBranch(self.branchName)
            for ptWpBin in self.ptWpBins:
                for UD in ['Up','Down']:
                    self.addBranch(self.branchName + '_{ptWpBinName}_{UD}'.format(ptWpBinName=ptWpBin[6], UD=UD))

            self.sampleTree.addFormula(self.isHeavyFlavorProcess)
            #fnmatch.fnmatch
            if any([fnmatch.fnmatch(self.sample.identifier,x) for x in self.applyToSamples]): 
                self.applyForThisSample = True
            if any([fnmatch.fnmatch(self.sample.identifier,x) for x in self.applyToHeavyFlavorOnly]): 
                self.checkHeavyFlavor = True

        if self.applyForThisSample:
            print("INFO: Double btag SFs are applied for this sample!")
            if self.checkHeavyFlavor:
                print("INFO: -> only for heavy flavor events!")
        else:
            print("INFO: Double btag SFs are not applied for this sample!")


    def applies(self):
        return self.applyForThisSample and (not self.checkHeavyFlavor or self.sampleTree.evaluate(self.isHeavyFlavorProcess))

    def processEvent(self, tree):

        if not self.hasBeenProcessed(tree) and not self.isData:
            self.markProcessed(tree)

            jetPt = getattr(tree, self.jetPtName)
            jetScore = getattr(tree, self.taggerName)

            jetSelection = lambda x,j : (x.FatJet_lepFilter[j]>0 and abs(x.FatJet_eta[j]) < 2.5)

            # Nominal
            sf = 1.0
            if self.applies():
                for i in range(tree.nFatJet):
                    if jetSelection(tree, i): 
                        sf *= self.get_jet_sf(jetPt[i], jetScore[i])
            self._b(self.branchName)[0] = sf

            # variations
            for ptWpBin in self.ptWpBins:
                for UD in ['Up','Down']:
                    sf = 1.0
                    if self.applies():
                        for i in range(tree.nFatJet):
                            if jetSelection(tree, i): 
                                wpBinMatches = [k for k,v in self.bins.items() if v(jetScore[i])]
                                if jetPt[i] >= ptWpBin[1] and jetPt[i] < ptWpBin[2] and len(wpBinMatches) == 1 and wpBinMatches[0] == ptWpBin[0]:
                                    variation = 1.0 if UD == 'Up' else -1.0
                                else:
                                    variation = 0
                                sf *= self.get_jet_sf(jetPt[i], jetScore[i], variation=variation)

                    self._b(self.branchName + '_{ptWpBinName}_{UD}'.format(ptWpBinName=ptWpBin[6], UD=UD))[0] = sf



if __name__ == "__main__":
    doubleBtagSF = DoubleBTagWeightsFromCSV(year=2017)
    doubleBtagSF.customInit(initVars={'sample': type('obj', (object,), {'type' : 'MC'})})
    print(doubleBtagSF.get_jet_sf(270, 0.88))
    print(doubleBtagSF.get_jet_sf(270, 0.88, 1.0))
    print(doubleBtagSF.get_jet_sf(270, 0.88, -1.0))
    print(doubleBtagSF.get_jet_sf(650, 0.98))
    print(doubleBtagSF.get_jet_sf(650, 0.98, 1.0))
    print(doubleBtagSF.get_jet_sf(650, 0.98, -1.0))
