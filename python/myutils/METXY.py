#!/usr/bin/env python
from __future__ import print_function
import ROOT
from BranchTools import Collection
from BranchTools import AddCollectionsModule
import array
import os
import math
import numpy as np
from XbbConfig import XbbConfigTools

# MET X/Y correction 
class METXY(AddCollectionsModule):

    def __init__(self, year):
        super(METXY, self).__init__()
        self.debug = 'XBBDEBUG' in os.environ
        self.year  = int(year)
        self.quickloadWarningShown = False

    def customInit(self, initVars):
        self.sampleTree  = initVars['sampleTree']
        self.sample      = initVars['sample']
        self.config      = initVars['config']
        self.xbbConfig   = XbbConfigTools(self.config)
        self.systematics = self.xbbConfig.getJECuncertainties(step='METXY')
        self.METsystematics = [x for x in self.systematics if 'jerReg' not in x] + ['unclustEn']

        # load METXYCorr_Met_MetPhi from VHbb namespace
        VHbbNameSpace = self.config.get('VHbbNameSpace', 'library')
        ROOT.gSystem.Load(VHbbNameSpace)

        self.MET_Pt     = array.array('f', [0.0])
        self.MET_Phi    = array.array('f', [0.0])
        self.sampleTree.tree.SetBranchAddress("MET_Pt", self.MET_Pt)
        self.sampleTree.tree.SetBranchAddress("MET_Phi", self.MET_Phi)

        self.addBranch("MET_Pt_uncorrected")
        self.addBranch("MET_Phi_uncorrected")

        if self.sample.isMC():
            self.MET_Pt_syst  = {}
            self.MET_Phi_syst = {}
            for syst in self.METsystematics:
                self.MET_Pt_syst[syst] = {}
                self.MET_Phi_syst[syst] = {}
                for Q in self._variations(syst):
                    self.MET_Pt_syst[syst][Q]  = array.array('f', [0.0])
                    self.MET_Phi_syst[syst][Q] = array.array('f', [0.0])
                    self.sampleTree.tree.SetBranchAddress("MET_pt_"+syst+Q, self.MET_Pt_syst[syst][Q])
                    self.sampleTree.tree.SetBranchAddress("MET_phi_"+syst+Q, self.MET_Phi_syst[syst][Q])

                    self.addBranch("MET_pt_uncorrected_"+syst+Q)
                    self.addBranch("MET_phi_uncorrected_"+syst+Q)


    def processEvent(self, tree):
        if not self.hasBeenProcessed(tree):
            self.markProcessed(tree)

            # backup uncorrected branches
            self._b('MET_Pt_uncorrected')[0]  = tree.MET_Pt
            self._b('MET_Phi_uncorrected')[0] = tree.MET_Phi

            MET_Pt_corrected, MET_Phi_corrected = ROOT.VHbb.METXYCorr_Met_MetPhi(tree.MET_Pt, tree.MET_Phi, tree.run, self.year, self.sample.isMC(), tree.PV_npvs)

            # overwrite MET_Pt, MET_Phi branches
            self.MET_Pt[0]  = MET_Pt_corrected
            self.MET_Phi[0] = MET_Phi_corrected

            if self.sample.isMC():
                for syst in self.METsystematics:
                    for Q in self._variations(syst):
                        
                        # backup uncorrected branches
                        self._b("MET_pt_uncorrected_"+syst+Q)[0]  = self.MET_Pt_syst[syst][Q][0]
                        self._b("MET_phi_uncorrected_"+syst+Q)[0] = self.MET_Phi_syst[syst][Q][0]

                        MET_Pt_corrected, MET_Phi_corrected = ROOT.VHbb.METXYCorr_Met_MetPhi(self.MET_Pt_syst[syst][Q][0], self.MET_Phi_syst[syst][Q][0], tree.run, self.year, self.sample.isMC(), tree.PV_npvs)

                        # overwrite MET_Pt, MET_Phi branches
                        self.MET_Pt_syst[syst][Q][0]  = MET_Pt_corrected
                        self.MET_Phi_syst[syst][Q][0] = MET_Phi_corrected
            
            # formulas by default reload the branch content when evaluating the first instance of the object!
            # SetQuickLoad(1) turns off this behavior
            for formulaName, treeFormula in self.sampleTree.formulas.items():
                if 'MET' in formulaName:
                    if not self.quickloadWarningShown:
                        self.quickloadWarningShown = True
                        print("INFO: SetQuickLoad(1) called for formula:", formulaName)
                        print("INFO: -> EvalInstance(0) on formulas will not re-load branches but will take values from memory, which might have been modified by this module.") 
                    treeFormula.SetQuickLoad(1)


