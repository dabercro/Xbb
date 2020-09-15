#!/usr/bin/env python
import ROOT
import numpy as np
import array
import os
from BranchTools import Collection
from BranchTools import AddCollectionsModule
from XbbTools import XbbTools
from XbbConfig import XbbConfigTools

class isBoosted(AddCollectionsModule):

    def __init__(self, branchName='isBoosted', cutName='all_BOOST', useFlags=False, flags=None):
        super(isBoosted, self).__init__()
        self.branchName = branchName
        self.cutName = cutName
        self.version = 4
        self.useFlags = useFlags
        self.flags    = ['resolvedCR','resolvedSR','boostedCR','boostedSR'] if flags is None else flags
        self.variations = self._variations("isBoosted") 

    # returns cut string with variables replaced by their systematic variations
    def getSystVarCut(self, cut, syst, UD):
        replacementRulesList = XbbTools.getReplacementRulesList(self.config, syst)
        systVarCut           = XbbTools.getSystematicsVariationTemplate(cut, replacementRulesList)
        systVarCut           = systVarCut.replace('{syst}', syst).replace('{UD}', UD)
        return systVarCut

    def customInit(self, initVars):
        self.sample     = initVars['sample']
        self.sampleTree = initVars['sampleTree']
        self.config     = initVars['config']
        self.xbbConfig  = XbbConfigTools(self.config)

        self.boostedCut = XbbTools.sanitizeExpression(self.config.get('Cuts', self.cutName), config=self.config)
        self.systVarCuts = {}

        #self.systematics = sorted(list(set(sum([eval(self.config.get('LimitGeneral', x)) for x in ['sys_cr', 'sys_BDT', 'sys_Mjj']], []))))
        self.systematics = self.xbbConfig.getJECuncertainties(step='BoostedFlags') + ['unclustEn']

        # Nominal
        self.addIntegerBranch(self.branchName)
        self.sampleTree.addFormula(self.boostedCut)

        # systematic variations
        if self.sample.isMC():
            for syst in self.systematics:
                for UD in self.variations: 
                    systVarBranchName = self._v(self.branchName, syst, UD)
                    self.addIntegerBranch(systVarBranchName)
                    self.systVarCuts[systVarBranchName] = self.getSystVarCut(self.boostedCut, syst=syst, UD=UD)
                    self.sampleTree.addFormula(self.systVarCuts[systVarBranchName])

        if self.useFlags:
            # CR/SR flags
            self.flagCuts = {k: XbbTools.sanitizeExpression(self.config.get('Cuts', k), config=self.config) for k in self.flags} 
            self.systVarFlagCuts = {k:{} for k in self.flags}

            for k in self.flags:
                self.sampleTree.addFormula(self.flagCuts[k])
                self.addIntegerBranch(k)

            # systematic variations
            if self.sample.isMC():
                for k in self.flags:
                    for syst in self.systematics:
                        for UD in self.variations: 
                            systVarflagName = self._v(k, syst, UD)
                            self.addIntegerBranch(systVarflagName)
                            self.systVarFlagCuts[k][systVarflagName] = XbbTools.sanitizeExpression(self.getSystVarCut(self.flagCuts[k], syst=syst, UD=UD), config=self.config)
                            self.sampleTree.addFormula(self.systVarFlagCuts[k][systVarflagName])


    def processEvent(self, tree):
        # if current entry has not been processed yet
        if not self.hasBeenProcessed(tree):
            self.markProcessed(tree)

            # Nominal
            b = int(self.sampleTree.evaluate(self.boostedCut))
            self._b(self._v(self.branchName))[0] = 1 if b > 0 else 0 

            # systematic variations
            if self.sample.isMC():
                for syst in self.systematics:
                    for UD in self.variations: 
                        systVarBranchName = self._v(self.branchName, syst, UD)
                        b = int(self.sampleTree.evaluate(self.systVarCuts[systVarBranchName]))
                        self._b(systVarBranchName)[0] = 1 if b > 0 else 0 

            if self.useFlags:
                # CR/SR flags
                for k in self.flags:
                    b =  int(self.sampleTree.evaluate(self.flagCuts[k]))
                    self._b(self._v(k))[0] = 1 if b > 0 else 0
                if self.sample.isMC():
                    for k in self.flags:
                        for syst in self.systematics:
                            for UD in self.variations: 
                                systVarflagName = self._v(k, syst, UD)
                                b = int(self.sampleTree.evaluate(self.systVarFlagCuts[k][systVarflagName]))
                                self._b(systVarflagName)[0] = 1 if b > 0 else 0


