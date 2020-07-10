#!/usr/bin/env python
import ROOT
from BranchTools import Collection
from BranchTools import AddCollectionsModule
import array
import os
import numpy as np

# apply x/y/phi corrections to data
class METcorr(AddCollectionsModule):

    def __init__(self, debug=False, year=None, debugEvents=[], backupPreviousCorrection=True):
        self.backupPreviousCorrection = backupPreviousCorrection
        self.debug = debug or 'XBBDEBUG' in os.environ
        self.year = year if type(year)==str else str(year)

    def customInit(self, initVars):
        self.sampleTree = initVars['sampleTree']
        self.isData = initVars['sample'].isData()
        self.sample = initVars['sample']

        # new branches to write
        self.addBranch("METcorr_x")
        self.addBranch("METcorr_y")
        self.MET_Phi = array.array('f', [-99.0]*1)
        self.MET_Pt  = array.array('f', [-99.0]*1)
  
        self.sampleTree.tree.SetBranchAddress("MET_Phi", self.MET_Phi)
        self.sampleTree.tree.SetBranchAddress("MET_Pt", self.MET_Pt)

        if self.backupPreviousCorrection:
            self.addBranch("MET_PhiOld")
            self.addBranch("MET_PtOld")

        self.corr = {
        "2016": {
            "2016MC": {"METxcorr":"-(-0.195191*npv -0.170948)", "METycorr":"-(-0.0311891*npv +0.787627)"},
            "2016B": {"METxcorr":"-(-0.0478335*npv -0.108032)", "METycorr":"-(0.125148*npv +0.355672)"},
            "2016C": {"METxcorr":"-(-0.0916985*npv +0.393247)", "METycorr":"-(0.151445*npv +0.114491)"},
            "2016D": {"METxcorr":"-(-0.0581169*npv +0.567316)", "METycorr":"-(0.147549*npv +0.403088)"},
            "2016E": {"METxcorr":"-(-0.065622*npv +0.536856)", "METycorr":"-(0.188532*npv +0.495346)"},
            "2016F": {"METxcorr":"-(-0.0313322*npv +0.39866)", "METycorr":"-(0.16081*npv +0.960177)"},
            "2016G": {"METxcorr":"-(0.040803*npv -0.290384)", "METycorr": "-(0.0961935*npv +0.666096)"},
            "2016H": {"METxcorr":"-(0.0330868*npv -0.209534)", "METycorr":"-(0.141513*npv +0.816732)"} 
            },
         "2018": {
            "2018MC": {"METxcorr":"-(0.296713*npv -0.141506)", "METycorr":"-(0.115685*npv +0.0128193)"},
            "2018A": {"METxcorr":"-(0.362865*npv -1.94505)", "METycorr":"-(0.0709085*npv -0.307365)"},
            "2018B": {"METxcorr":"-(0.492083*npv -2.93552)", "METycorr":"-(0.17874*npv -0.786844)"},
            "2018C": {"METxcorr":"-(0.521349*npv -1.44544)", "METycorr":"-(0.118956*npv -1.96434)"},
            "2018D": {"METxcorr":"-(0.531151*npv -1.37568)", "METycorr":"-(0.0884639*npv -1.57089)"}
            }
        }

    def processEvent(self, tree):
        if not self.hasBeenProcessed(tree):
            self.markProcessed(tree)

            if self.backupPreviousCorrection:           
                self._b("MET_PhiOld")[0]  = self.MET_Phi[0]
                self._b("MET_PtOld")[0]   = self.MET_Pt[0]

            self._b("METcorr_x")[0] = -999.0
            self._b("METcorr_y")[0] = -999.0

            debugEvent = [tree.run, tree.event] in self.debugEvents
            if debugEvent:
                print "DEBUG-EVENT:", tree.run, tree.event

            runnb = tree.run 
            runera = self.get_era(runnb)
            npv = tree.PV_npvs

            if(npv>100): npv=100;

            MET_p4 = ROOT.TLorentzVector()
            MET_p4.SetPtEtaPhiM(tree.MET_Pt, 0.0, tree.MET_Phi, 0.0)
            MET = MET_p4.E() 
            MET_x = MET_p4.X()
            MET_y = MET_p4.Y()
            MET_phi = MET_p4.Phi()

            if (runera != -1):

                METxcorr = eval(self.corr[self.year][runera]["METxcorr"])
                METycorr = eval(self.corr[self.year][runera]["METycorr"])

                METcorr_x = MET * np.cos(MET_phi) + METxcorr 
                METcorr_y = MET * np.sin(MET_phi) + METycorr
                METcorr_pt = np.sqrt(METcorr_x * METcorr_x + METcorr_y * METcorr_y)

                if (METcorr_x == 0) and (METcorr_y > 0): METcorr_phi = np.pi 
                elif (METcorr_x == 0) and (METcorr_y < 0): METcorr_phi = -np.pi 
                elif (METcorr_x > 0): METcorr_phi = np.arctan(METcorr_y/METcorr_x) 
                elif (METcorr_x < 0) and (METcorr_y > 0): METcorr_phi = np.arctan(METcorr_y/METcorr_x)+ np.pi 
                elif (METcorr_x < 0) and (METcorr_y < 0): METcorr_phi = np.arctan(METcorr_y/METcorr_x)- np.pi 
                else: METcorr_phi = 0 

                self._b("METcorr_x")[0] = METcorr_x
                self._b("METcorr_y")[0] = METcorr_y
                self.MET_Phi[0] = METcorr_phi
                self.MET_Pt[0]  = METcorr_pt


    def get_era(self, runnb):

        runera =-1

        # Data        
        if (self.isData) and (self.year == "2016") and (runnb >= 272007) and (runnb <= 275376): runera = "2016B"
        if (self.isData) and (self.year == "2016") and (runnb >= 275657) and (runnb <= 276283): runera = "2016C"
        if (self.isData) and (self.year == "2016") and (runnb >= 276315) and (runnb <= 276811): runera = "2016D"
        if (self.isData) and (self.year == "2016") and (runnb >= 276831) and (runnb <= 277420): runera = "2016E"
        if (self.isData) and (self.year == "2016") and (runnb >= 277772) and (runnb <= 278808): runera = "2016F"
        if (self.isData) and (self.year == "2016") and (runnb >= 278820) and (runnb <= 280385): runera = "2016G"
        if (self.isData) and (self.year == "2016") and (runnb >= 280919) and (runnb <= 284044): runera = "2016H"

        # MC
        if (not self.isData) and (self.year == "2016"): runera = "2016MC"

        if (self.year == "2018"):
            try:
                if (self.isData):                        
                    if (runnb >=315252 and runnb<=316995): runera = "2018A"
                    if (runnb >=316998 and runnb<=319312): runera = "2018B"
                    if (runnb >=319313 and runnb<=320393): runera = "2018C"
                    if (runnb >=320394 and runnb<=325273): runera = "2018D"
                if not (self.isData):
                    runera = "2018MC"
            except:
                print "Run number {} not found for {}".format(runnb,self.year)
                            
        return runera

    def afterProcessing(self):
        print "MET corrected!"
