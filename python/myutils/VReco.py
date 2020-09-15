#!/usr/bin/env python
import ROOT
from BranchTools import Collection
from BranchTools import AddCollectionsModule
import array
import os
import numpy as np
from XbbConfig import XbbConfigTools

# do jet/lepton selection and skimming
class VReco(AddCollectionsModule):

    def __init__(self, debug=False, replaceNominal=False, puIdCut=6, jetIdCut=4):
        self.debug = debug or 'XBBDEBUG' in os.environ
        self.replaceNominal = replaceNominal
        self.puIdCut = puIdCut 
        self.jetIdCut = jetIdCut
        super(VReco, self).__init__()
        self.version = 2

    def customInit(self, initVars):
        self.sampleTree = initVars['sampleTree']
        self.isData     = initVars['sample'].isData()
        self.sample     = initVars['sample']
        self.config     = initVars['config']
        self.xbbConfig  = XbbConfigTools(self.config)

        self.systematics = self.xbbConfig.getJECuncertainties(step='VReco') + ['unclustEn']

        #if self.replaceNominal:
        self.allVariations = ['Nominal']
        #else:
        #    self.allVariations = []
       
        # jet and MET systematics for MC
        if not self.isData:
            self.allVariations += self.systematics 

        for syst in self.allVariations:
            self.addVsystematics(syst)

        self.addBranch("selLeptonWln_pt")
        self.addBranch("selLeptonWln_eta")
        self.addBranch("selLeptonWln_phi")
        self.addBranch("selLeptonWln_mass")

    def addVsystematics(self, syst):
        if self._isnominal(syst): 
            # replace values from post-processor / selection
            if self.replaceNominal:
                self.addBranch("V_pt")
                self.addBranch("V_eta")
                self.addBranch("V_phi")
                self.addBranch("V_mass")
                self.addBranch("V_mt")
            self.addBranch("MET_sig30")
            self.addBranch("MET_sig30puid")
        else:
            for UD in self._variations(syst):
                self.addBranch(self._v("V_mt", syst, UD))
                self.addBranch(self._v("V_pt", syst, UD))
                self.addBranch(self._v("V_eta", syst, UD))
                self.addBranch(self._v("V_phi", syst, UD))
                self.addBranch(self._v("V_mass", syst, UD))
                if syst not in ['jerReg']:
                    self.addBranch(self._v("MET_sig30", syst, UD))
                    self.addBranch(self._v("MET_sig30puid", syst, UD))

    def processEvent(self, tree):
        if not self.hasBeenProcessed(tree):
            self.markProcessed(tree)
            self._b("selLeptonWln_pt")[0] = -99.0
            self._b("selLeptonWln_eta")[0] = -99.0
            self._b("selLeptonWln_phi")[0] = -99.0
            self._b("selLeptonWln_mass")[0] = -99.0

            for syst in self.allVariations:
                directions = [''] if syst.lower() == 'nominal' else ['Up','Down']
                for UD in self._variations(syst):
                    if not self._isnominal(syst) or self.replaceNominal:
                        self._b(self._v("V_mt", syst, UD))[0] = -1.0

                    if tree.Vtype == 0 or tree.Vtype == 1:
                        lep1 = ROOT.TLorentzVector()
                        lep2 = ROOT.TLorentzVector()
                        i1 = tree.vLidx[0]
                        i2 = tree.vLidx[1]
                        if tree.Vtype==1:
                            lep1.SetPtEtaPhiM(tree.Electron_pt[i1], tree.Electron_eta[i1], tree.Electron_phi[i1], tree.Electron_mass[i1])
                            lep2.SetPtEtaPhiM(tree.Electron_pt[i2], tree.Electron_eta[i2], tree.Electron_phi[i2], tree.Electron_mass[i2])
                        else:
                            lep1.SetPtEtaPhiM(tree.Muon_pt[i1], tree.Muon_eta[i1], tree.Muon_phi[i1], tree.Muon_mass[i1])
                            lep2.SetPtEtaPhiM(tree.Muon_pt[i2], tree.Muon_eta[i2], tree.Muon_phi[i2], tree.Muon_mass[i2])
                        V = lep1 + lep2
                    elif  tree.Vtype == 2 or tree.Vtype == 3:
                        i1 = tree.vLidx[0]
                        if tree.Vtype==3:
                            sel_lepton_pt   = tree.Electron_pt[i1]
                            sel_lepton_eta  = tree.Electron_eta[i1]
                            sel_lepton_phi  = tree.Electron_phi[i1]
                            sel_lepton_mass = tree.Electron_mass[i1]
                        else:
                            sel_lepton_pt   = tree.Muon_pt[i1]
                            sel_lepton_eta  = tree.Muon_eta[i1]
                            sel_lepton_phi  = tree.Muon_phi[i1]
                            sel_lepton_mass = tree.Muon_mass[i1]

                        MET = ROOT.TLorentzVector()
                        Lep = ROOT.TLorentzVector()
                        if syst.lower()=='nominal' or syst=='jerReg':
                            MET.SetPtEtaPhiM(tree.MET_Pt, 0.0, tree.MET_Phi, 0.0)
                        else:
                            MET.SetPtEtaPhiM(getattr(tree, "MET_pt_{syst}{UD}".format(syst=syst, UD=UD)), 0.0, getattr(tree, "MET_phi_{syst}{UD}".format(syst=syst, UD=UD)), 0.0)
                        Lep.SetPtEtaPhiM(sel_lepton_pt, sel_lepton_eta, sel_lepton_phi, sel_lepton_mass)
                        cosPhi12 = (Lep.Px()*MET.Px() + Lep.Py()*MET.Py()) / (Lep.Pt() * MET.Pt())
                        if not self._isnominal(syst) or self.replaceNominal:
                            self._b(self._v("V_mt", syst,UD))[0] = ROOT.TMath.Sqrt(2*Lep.Pt()*MET.Pt() * (1 - cosPhi12))

                        V = MET + Lep

                        if self._isnominal(syst):
                            self._b("selLeptonWln_pt")[0] = Lep.Pt()
                            self._b("selLeptonWln_eta")[0] = Lep.Eta()
                            self._b("selLeptonWln_phi")[0] = Lep.Phi()
                            self._b("selLeptonWln_mass")[0] = Lep.M()
                    elif tree.Vtype == 4:
                        MET = ROOT.TLorentzVector()
                        if syst.lower()=='nominal':
                            MET.SetPtEtaPhiM(tree.MET_Pt, 0.0, tree.MET_Phi, 0.0)
                        else:
                            MET.SetPtEtaPhiM(getattr(tree, "MET_pt_{syst}{UD}".format(syst=syst, UD=UD)), 0.0, getattr(tree, "MET_phi_{syst}{UD}".format(syst=syst, UD=UD)), 0.0)
                        V = MET
                    else:
                        V = None
                        if (not self._isnominal(syst)) or self.replaceNominal:
                            self._b(self._v("V_pt", syst, UD))[0] = -1.0
                            self._b(self._v("V_eta", syst, UD))[0] = -1.0
                            self._b(self._v("V_phi", syst, UD))[0] = -1.0
                            self._b(self._v("V_mass", syst, UD))[0] = -1.0

                    if V is not None:
                        if (not self._isnominal(syst)) or self.replaceNominal:
                            self._b(self._v("V_pt", syst, UD))[0] = V.Pt()
                            self._b(self._v("V_eta", syst, UD))[0] = V.Eta()
                            self._b(self._v("V_phi", syst, UD))[0] = V.Phi()
                            self._b(self._v("V_mass", syst, UD))[0] = V.M()

                    # MET significance (approx.)
                    if syst != 'jerReg':
                        HTsum30 = 0
                        HTsum30puid = 0
                        if syst.lower() == 'nominal' or syst == 'unclustEn':
                            jetPt = tree.Jet_Pt
                        else:
                            jetPt = getattr(tree, "Jet_pt_{syst}{UD}".format(syst=syst, UD=UD))
                        for i in range(tree.nJet):
                            if jetPt[i]>30 and tree.Jet_lepFilter[i]>0 and tree.Jet_jetId[i] > self.jetIdCut:
                                HTsum30 += jetPt[i]
                                if tree.Jet_puId[i]>self.puIdCut or jetPt[i]>50.0:
                                    HTsum30puid += jetPt[i]
                        
                        if syst.lower() == 'nominal':
                            metPt = tree.MET_Pt
                        else:
                            metPt = getattr(tree, "MET_pt_{syst}{UD}".format(syst=syst, UD=UD))

                        self._b(self._v("MET_sig30", syst, UD))[0] = (metPt / np.sqrt(HTsum30)) if HTsum30 > 0 else -1.0 
                        self._b(self._v("MET_sig30puid", syst, UD))[0] = (metPt / np.sqrt(HTsum30puid)) if HTsum30puid > 0 else -1.0

