#!/usr/bin/env python
import ROOT
from BranchTools import Collection
from BranchTools import AddCollectionsModule
import array
import os
import itertools

from vLeptons import vLeptonSelector

# do jet/lepton selection and skimming
class VHbbSelection(AddCollectionsModule):

    # original 2017: puIdCut=0, jetIdCut=-1
    # now:           puIdCut=6, jetIdCut=4
    def __init__(self, debug=False, year="none", channels=["Wln","Zll","Znn"], puIdCut=6, jetIdCut=4, debugEvents=[], isoWen=0.06, isoWmn=0.06, isoZmm=0.25, isoZee=0.15, recomputeVtype=False, btagMaxCut=None, idWmn="default", idWmnCut=1, idWen="default", idWenCut=1, idZmm="default", idZmmCut=1, idZee="default", idZeeCut=1):
        super(VHbbSelection, self).__init__()
        self.debug = debug or 'XBBDEBUG' in os.environ
        self.version = 6
        self.stats = {}
        self.year = year
        self.channels = channels
        self.debugEvents = debugEvents
        self.puIdCut = puIdCut
        self.jetIdCut = jetIdCut
        self.isoWen = isoWen
        self.isoWmn = isoWmn
        self.isoZee = isoZee
        self.isoZmm = isoZmm
        self.idWmn    = idWmn
        self.idWmnCut = idWmnCut
        self.idWen    = idWen
        self.idWenCut = idWenCut
        self.idZmm    = idZmm
        self.idZmmCut = idZmmCut
        self.idZee    = idZee
        self.idZeeCut = idZeeCut
        self.recomputeVtype = recomputeVtype
        self.btagMaxCut = btagMaxCut
        # only use puId below this pT:
        self.puIdMaxPt = 50.0
        # run QCD estimation analysis with inverted isolation cuts (CAREFUL!)
        if self.recomputeVtype:
            print("\x1b[45m\x1b[37m" + "-"*160 + "\x1b[0m")
            print("\x1b[45m\x1b[37m" + "-"*160 + "\x1b[0m")
            print("\x1b[45m\x1b[37m RECOMPUTE Vtype, custom definitions might be used!\x1b[0m") 
            print("\x1b[45m\x1b[37m" + "-"*160 + "\x1b[0m")
            print("\x1b[45m\x1b[37m" + "-"*160 + "\x1b[0m")

    def customInit(self, initVars):
        self.sampleTree = initVars['sampleTree']
        self.isData = initVars['sample'].isData()
        self.sample = initVars['sample']
        self.config = initVars['config']

        # settings
        # originally Electron_mvaFall17V1Iso_WP80 was used for 2017, which was called Electron_mvaFall17Iso_WP80
        # consistently available for all 3 years: Electron_mvaFall17V2Iso
        self.electronID = {
                    1: {
                                "2018": "Electron_mvaFall17V2Iso_WP80",
                                "2017": "Electron_mvaFall17V2Iso_WP80",
                                #"2016": "Electron_mvaSpring16GP_WP80",
                                "2016": "Electron_mvaFall17V2Iso_WP80",
                            },
                    2: {
                                "2018": "Electron_mvaFall17V2Iso_WP90",
                                "2017": "Electron_mvaFall17V2Iso_WP90",
                                #"2016": "Electron_mvaSpring16GP_WP90",
                                "2016": "Electron_mvaFall17V2Iso_WP90",
                            }
                    }
       
        # default lepton IDs for backward compatibility
        if self.idWen == "default":
            self.idWen = self.electronID[1][self.year]
        if self.idZee == "default":
            self.idZee = self.electronID[2][self.year]
        if self.idWmn == "default":
            # cut-based muon ID
            self.idWmn = "Muon_tightId"
        if self.idZmm == "default":
            self.idZmm = None 


        # updated to july 2018 Jet/MET recommendations
        # reference: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
        if self.year == "2016":
            self.metFilters = ["Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter"]
        elif self.year in ["2017","2018"]:
            self.metFilters = ["Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter", "Flag_ecalBadCalibFilter"]
        if self.isData:
            self.metFilters.append("Flag_eeBadScFilter")

        # main tagger (-> hJidx)
        self.taggerName = "Jet_btagDeepB"

        ## alternative jet selections to check (->hJidx_*)
        self.jetDefinitions = []

        # additional jet definitions for DeepJet
        if self.year in ["2016","2017", "2018"]:
            self.jetDefinitions = [
                    {'taggerName': 'Jet_btagDeepB'},
                    #{'taggerName': 'Jet_btagDeepFlavB', 'indexName':'DeepFlavB'},
                    ]
        else:
            self.jetDefinitions = [
                    {'taggerName': 'Jet_btagDeepB'},
                    ]

        self.btagWPs = {
                       "2018": {
                      'Jet_btagDeepB': {
                            'loose':  0.1241,
                            'medium': 0.4184,
                            'tight':  0.7527,
                            'none':   -1.0,
                            },
                      'Jet_btagDeepFlavB': {
                            'loose':  0.0494, 
                            'medium': 0.2770,
                            'tight':  0.7264,
                            'none':   -1.0,
                            },
                        },
                        "2017": {
                        'Jet_btagDeepB': {
                            'loose':  0.1522,
                            'medium': 0.4941,
                            'tight':  0.8001,
                            'none':   -1.0,
                            },
                        'Jet_btagDeepFlavB': {
                            'loose':  0.0521,
                            'medium': 0.3033,
                            'tight':  0.7489,
                            'none':   -1.0,
                            },
                        },
                        "2016": {
                        'Jet_btagDeepB': {
                            'loose':  0.2217,
                            'medium': 0.6321,
                            'tight':  0.8953,
                            'none':   -1.0,
                            },
                        'Jet_btagDeepFlavB': {
                            'loose':  0.0614,
                            'medium': 0.3093,
                            'tight':  0.7221,
                            'none':   -1.0,
                            },
                        },
                   }
        # for main tagger, for others use self.btagWPs below
        self.btagWP = self.btagWPs[self.year][self.taggerName]

        if self.year == "2018":
            self.HltPaths = {
                        'Znn': ['HLT_PFMET120_PFMHT120_IDTight'],
                        'Wln': ['HLT_Ele32_WPTight_Gsf','HLT_IsoMu24'],
                        'Zll': ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL'],
                        }
        elif self.year == "2017":
            self.HltPaths = {
                        'Znn': ['HLT_PFMET120_PFMHT120_IDTight','HLT_PFMET120_PFMHT120_IDTight_PFHT60'],
                        'Wln': ['HLT_Ele32_WPTight_Gsf_L1DoubleEG','HLT_IsoMu27'],
                        'Zll': ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL'],
                        }
        elif self.year == "2016":
            self.HltPaths = {
                        'Znn': ['HLT_PFMET110_PFMHT110_IDTight','HLT_PFMET120_PFMHT120_IDTight','HLT_PFMET170_NoiseCleaned','HLT_PFMET170_HBHECleaned','HLT_PFMET170_HBHE_BeamHaloCleaned'],
                        'Wln': ['HLT_Ele27_WPTight_Gsf','HLT_IsoMu24','HLT_IsoTkMu24'],
                        'Zll': ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL','HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL','HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ','HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ','HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ'],
                        }
        else:
            raise Execption("unknown year")

        # Vtype cut on datasets
        self.leptonFlav = {
                'DoubleMuon': 0,
                'DoubleEG': 1,
                'SingleMuon': 2,
                'SingleElectron': 3,
                'MET': 4,
                }

        if self.year=='2018':
            if "Zll" in self.channels:
                self.leptonFlav['EGamma']=1
            elif "Wlv" in self.channels:
                self.leptonFlav['EGamma']=3

        self.systematics = ['jer','jerReg','jesAbsoluteStat','jesAbsoluteScale','jesAbsoluteFlavMap','jesAbsoluteMPFBias','jesFragmentation','jesSinglePionECAL','jesSinglePionHCAL','jesFlavorQCD','jesRelativeJEREC1','jesRelativeJEREC2','jesRelativeJERHF','jesRelativePtBB','jesRelativePtEC1','jesRelativePtEC2','jesRelativePtHF','jesRelativeBal','jesRelativeFSR','jesRelativeStatFSR','jesRelativeStatEC','jesRelativeStatHF','jesPileUpDataMC','jesPileUpPtRef','jesPileUpPtBB','jesPileUpPtEC1','jesPileUpPtEC2','jesPileUpPtHF','jesPileUpMuZero','jesPileUpEnvelope','jesTotal']
        self.METsystematics = [x for x in self.systematics if 'jerReg' not in x] + ['unclustEn']
        self.systematicsBoosted = [x for x in self.systematics if 'jerReg' not in x] + ['jms', 'jmr']

        self.cutFlow = [0] * 16

        # new branches to write
        self.addIntegerBranch("isZee")
        self.addIntegerBranch("isZmm")
        self.addIntegerBranch("isWmunu")
        self.addIntegerBranch("isWenu")
        self.addIntegerBranch("isZnn")

        # selected lepton indices
        self.addIntegerVectorBranch("vLidx", default=-1, length=2)
        
        # selected jet indices
        defaultJetDefinitions = [x for x in self.jetDefinitions if 'indexName' not in x or x['indexName']=='']
        assert (len(defaultJetDefinitions) == 1)

        for jetDefinition in self.jetDefinitions:
            indexName = jetDefinition['indexName'] if 'indexName' in jetDefinition else '' 
            self.addIntegerVectorBranch("hJidx{suffix}".format(suffix=indexName), default=-1, length=2)
            if self.sample.isMC():
                for syst in self.systematics:
                    for UD in ['Up','Down']:
                        self.addIntegerVectorBranch("hJidx{suffix}_{syst}{UD}".format(suffix=indexName, syst=syst, UD=UD), default=-1, length=2)

        # selected fatjet index
        self.addIntegerBranch("Hbb_fjidx")
        if self.sample.isMC():
            for syst in self.systematicsBoosted:
                for UD in ['Up','Down']:
                    self.addIntegerBranch("Hbb_fjidx_{syst}{UD}".format(syst=syst, UD=UD))

        self.addBranch("lepMetDPhi")
        self.addBranch("V_mt")
        self.addBranch("V_pt")
        self.addBranch("V_eta")
        self.addBranch("V_phi")
        self.addBranch("V_mass")

        print "DEBUG: sample identifier:", self.sample.identifier, " lep flav", self.leptonFlav, " -> ", self.leptonFlav[self.sample.identifier] if self.sample.identifier in self.leptonFlav else "UNDEFINED LEPTON FLAVOR!!"

    def HighestPtGoodElectronsOppCharge(self, tree, min_pt, max_rel_iso, idcut, etacut, isOneLepton, electronID=None):
        indices = []
        for i in range(tree.nElectron):
            passIso        = tree.Electron_pfRelIso03_all[i] < max_rel_iso 
            passID         = getattr(tree, electronID)[i] >= idcut if electronID is not None else True
            passAcceptance = abs(tree.Electron_eta[i]) < etacut and tree.Electron_pt[i] > min_pt
            if passAcceptance and passIso and passID: 
                if len(indices) < 1:
                    indices.append(i)
                    if isOneLepton:
                        break
                else:
                    if tree.Electron_charge[i] * tree.Electron_charge[indices[0]] < 0:
                        indices.append(i)
                        break
        return indices

    def HighestPtGoodMuonsOppCharge(self, tree, min_pt, max_rel_iso, idcut, etacut, isOneLepton, muonID=None):
        indices = []
        for i in range(tree.nMuon):
            passIso        = tree.Muon_pfRelIso04_all[i] < max_rel_iso 
            passID         = getattr(tree, muonID)[i] >= idcut if muonID is not None else True

            passAcceptance = abs(tree.Muon_eta[i]) < etacut and tree.Muon_pt[i] > min_pt

            if passAcceptance and passIso and passID: 
                if len(indices) < 1:
                    indices.append(i)
                    if isOneLepton:
                        break
                else:
                    if tree.Muon_charge[i] * tree.Muon_charge[indices[0]] < 0:
                        indices.append(i)
                        break
        return indices

    def HighestTaggerValueFatJet(self, tree, ptCut, msdCut, etaCut, taggerName, systList=None, Vphi=None, Vphi_syst=None):
        fatJetMaxTagger = []
        Hbb_fjidx       = []
        Pt              = []
        Msd             = []
        systematics     = systList if systList is not None else [None]
        nSystematics    = len(systematics)
            
        Pt_nom  = tree.FatJet_Pt
        Msd_nom = tree.FatJet_Msoftdrop
        eta     = tree.FatJet_eta
        phi     = tree.FatJet_phi
        lepFilter = tree.FatJet_lepFilter
        tagger  = getattr(tree, taggerName) 

        for syst in systematics: 
            Hbb_fjidx.append(-1)
            fatJetMaxTagger.append(-999.9)

            if syst is None:
                Pt.append(Pt_nom)
                Msd.append(Msd_nom)
            elif syst in ['jmrUp','jmrDown','jmsUp','jmsDown']:
                # for jms,jmr the branches contain the mass itself
                Pt.append(Pt_nom)
                Msd.append(getattr(tree, "FatJet_msoftdrop_"+syst))
            else:
                # for all other variations, the branches contain a multiplicative factor... duh
                Pt_scale  = getattr(tree, "FatJet_pt_{syst}".format(syst=syst)) 
                Pt.append([Pt_nom[i] * Pt_scale[i] for i in range(len(Pt_nom))])
                Msd_scale = getattr(tree, "FatJet_msoftdrop_{syst}".format(syst=syst)) 
                Msd.append([Msd_nom[i] * Msd_scale[i] for i in range(len(Msd_nom))])

        for i in range(tree.nFatJet):
            for j in range(nSystematics):
                selectedVphi = Vphi_syst[systematics[j]] if (Vphi_syst is not None and systematics[j] in Vphi_syst) else Vphi
                #if Pt[j][i] > ptCut and abs(eta[i]) < etaCut and Msd[j][i] > msdCut and lepFilter[i]:
                if Pt[j][i] > ptCut and abs(eta[i]) < etaCut and Msd[j][i] > msdCut:
                    if selectedVphi is None or abs(ROOT.TVector2.Phi_mpi_pi(selectedVphi-phi[i]))>1.57:
                        if tagger[i] > fatJetMaxTagger[j]:
                            fatJetMaxTagger[j] = tagger[i]
                            Hbb_fjidx[j] = i

        # return list for systematic variations given and scalar otherwise
        if systList is not None:
            return Hbb_fjidx
        else:
            return Hbb_fjidx[0]

    # returns indices of the two second highest b-tagged jets passing the selection, if at least two exist
    # indices are sorted by b-tagger, descending
    def HighestTaggerValueBJets(self, tree, j1ptCut, j2ptCut, taggerName, puIdCut=0, jetIdCut=-1, maxJetPtCut=0, systList=None):
        jets = []
        indices = []
        nJet = tree.nJet
        systematics = systList if systList is not None else [None]

        jetTagger = getattr(tree, taggerName)
        Eta       = tree.Jet_eta
        jetId     = tree.Jet_jetId
        lepFilter = tree.Jet_lepFilter
        puId      = tree.Jet_puId
        PtReg_nom = tree.Jet_PtReg
        Pt_nom    = tree.Jet_Pt

        # momenta with sys variations applied
        PtReg = []
        Pt    = []
        nSystematics = len(systematics)
        for syst in systematics:
            indices.append([])
            if syst is None:
                PtReg.append(PtReg_nom)
                Pt.append(Pt_nom)
            elif syst == 'jerRegUp':
                PtReg.append(tree.Jet_PtRegUp)
                Pt.append(Pt_nom)
            elif syst == 'jerRegDown':
                PtReg.append(tree.Jet_PtRegDown)
                Pt.append(Pt_nom)
            else:
                Pt_syst    = getattr(tree, "Jet_pt_"+syst) 
                PtReg_syst = [PtReg_nom[i] * Pt_syst[i] / Pt_nom[i] for i in range(nJet)]
                Pt.append(Pt_syst)
                PtReg.append(PtReg_syst)

        for i in range(nJet):
            for j in range(nSystematics):
                pass_acceptance            = PtReg[j][i] > j1ptCut and abs(Eta[i]) < 2.5
                pass_jetIDandPUIDlepFilter = (puId[i] > puIdCut or Pt[j][i] > self.puIdMaxPt) and jetId[i] > jetIdCut and lepFilter[i]
                if pass_acceptance and pass_jetIDandPUIDlepFilter: 
                    if len(indices[j]) < 1:
                        indices[j].append(i)
                    else:
                        if jetTagger[i] > jetTagger[indices[j][0]]:
                            indices[j][0] = i
        for i in range(nJet):
            for j in range(nSystematics):
                if len(indices[j]) < 1 or i == indices[j][0]:
                    continue
                pass_acceptance            = PtReg[j][i] > j2ptCut and abs(Eta[i]) < 2.5
                pass_jetIDandPUIDlepFilter = (puId[i] > puIdCut or Pt[j][i] > self.puIdMaxPt) and jetId[i] > jetIdCut and lepFilter[i]

                if pass_acceptance and pass_jetIDandPUIDlepFilter:
                    if len(indices[j]) < 2:
                        indices[j].append(i)
                    else:
                        if jetTagger[i] > jetTagger[indices[j][1]]:
                            indices[j][1] = i

        for j in range(nSystematics):
            if len(indices[j]) > 1:
                if jetTagger[indices[j][1]] > jetTagger[indices[j][0]]:
                    indices = [indices[j][1], indices[j][0]]

            if len(indices[j]) == 2 and max(PtReg[j][indices[j][0]],PtReg[j][indices[j][1]]) < maxJetPtCut:
                indices[j] = []

        if systList is not None:
            return indices
        else:
            return indices[0]

    def processEvent(self, tree):

        if not self.hasBeenProcessed(tree):
            self.markProcessed(tree)

            for n in ["isZmm", "isZee", "isWenu", "isWmunu", "isZnn"]:
                self._b(n)[0] = 0
            self._b("lepMetDPhi")[0] = -99.0
            self._b("V_mt")[0] = -99.0
            self._b("V_pt")[0] = -99.0
            self._b("V_eta")[0] = -99.0
            self._b("V_phi")[0] = -99.0
            self._b("V_mass")[0] = -99.0
            self._b("hJidx")[0] = -1
            self._b("hJidx")[1] = -1
            self._b("vLidx")[0] = -1
            self._b("vLidx")[1] = -1

            self.cutFlow[0] += 1

            debugEvent = [tree.run, tree.event] in self.debugEvents
            if debugEvent:
                print "DEBUG-EVENT:", tree.run, tree.event

            # TRIGGER
            triggerPassed = {k: any([getattr(tree, x) for x in v if hasattr(tree,x)])  for k,v in self.HltPaths.items()}
            #print("------------------------------------")
            #for k,v in self.HltPaths.items():
            #    print(k,v)
            #    for x in v:
            #        print(x, hasattr(tree,x))
            #        if hasattr(tree,x):
            #            print(getattr(tree,x))

            #print(triggerPassed, "triggerPassed")
            #print(self.channels, "self.channels")
            #if not ((triggerPassed['0-lep'] and "Znn" in self.channels) or (triggerPassed['1-lep'] and "Wln" in self.channels) or (triggerPassed['2-lep'] and  "Zll" in self.channels)):

            if debugEvent:
                print "triggers:", triggerPassed

            Vtype = tree.Vtype
            if self.recomputeVtype:
                leptonSelector = vLeptonSelector(tree, config=self.config, isoZmm=self.isoZmm, isoZee=self.isoZee, isoWmn=self.isoWmn, isoWen=self.isoWen)
                customVtype = leptonSelector.getVtype()
                if customVtype != Vtype:
                    self.count("recomputeVtype_mismatch")
                else:
                    self.count("recomputeVtype_match")
                self.count("recomputeVtype_Vtype%d"%customVtype)
                Vtype = customVtype

            isZnn = Vtype == 4 and triggerPassed['Znn']
            isWln = Vtype in [2,3] and triggerPassed['Wln']
            isZll = Vtype in [0,1] and triggerPassed['Zll']

            if not any([isZll, isWln, isZnn]):
                return False
            self.cutFlow[1] += 1

            # LEPTONS
            if self.sample.identifier not in self.leptonFlav or self.leptonFlav[self.sample.identifier] == Vtype:
                if Vtype == 0:
                    good_muons_2lep = self.HighestPtGoodMuonsOppCharge(tree, min_pt=20.0, max_rel_iso=self.isoZmm, idcut=self.idZmmCut, etacut=2.4, isOneLepton=False, muonID=self.idZmm)
                    if len(good_muons_2lep) > 1:
                        self._b("isZmm")[0] = 1
                        self._b("vLidx")[0] = good_muons_2lep[0]
                        self._b("vLidx")[1] = good_muons_2lep[1]
                    elif debugEvent:
                        print "DEBUG-EVENT: 2 mu event, but less than 2 good muons -> discard"
                elif Vtype == 1:
                    good_elecs_2lep = self.HighestPtGoodElectronsOppCharge(tree, min_pt=20.0, max_rel_iso=self.isoZee, idcut=self.idZeeCut, etacut=2.5, isOneLepton=False, electronID=self.idZee)
                    if len(good_elecs_2lep) > 1:
                        self._b("isZee")[0] = 1
                        self._b("vLidx")[0] = good_elecs_2lep[0]
                        self._b("vLidx")[1] = good_elecs_2lep[1]
                    elif debugEvent:
                        print "DEBUG-EVENT: 2 e event, but less than 2 good electrons -> discard"
                elif Vtype == 2:
                    good_muons_1lep = self.HighestPtGoodMuonsOppCharge(tree, min_pt=25.0, max_rel_iso=self.isoWmn, idcut=self.idWmnCut, etacut=2.4, isOneLepton=True, muonID=self.idWmn)
                    if len(good_muons_1lep) > 0:
                        self._b("isWmunu")[0] = 1
                        self._b("vLidx")[0] = good_muons_1lep[0]
                        self._b("vLidx")[1] = -1
                    elif debugEvent:
                        print "DEBUG-EVENT: 1 mu event, but no good muon found"
                elif Vtype == 3:
                    good_elecs_1lep = self.HighestPtGoodElectronsOppCharge(tree, min_pt=30.0, max_rel_iso=self.isoWen, idcut=self.idWenCut, etacut=2.5, isOneLepton=True, electronID=self.idWen)
                    if len(good_elecs_1lep) > 0:
                        self._b("isWenu")[0] = 1
                        self._b("vLidx")[0] = good_elecs_1lep[0]
                        self._b("vLidx")[1] = -1
                    elif debugEvent:
                        print "DEBUG-EVENT: 1 e event, but no good electron found"
                elif Vtype == 4:
                        passMetFilters = all([getattr(tree, x) for x in self.metFilters]) 
                        if tree.MET_Pt > 170.0 and passMetFilters:
                            self._b("isZnn")[0] = 1
                        else:
                            if debugEvent:
                                print "DEBUG-EVENT: 0 lep event, but MET criteria not passed!"
                            return False
                else:
                    return False
            else:
                return False
            self.cutFlow[2] += 1

            # CHANNEL
            if (self._b("isZmm")[0] or self._b("isZee")[0]) and not ("Zll" in self.channels):
                return False
            #elif (self._b("isWmunu")[0] or self._b("isWenu")[0]) and not ("Wln" in self.channels or "Znn" in self.channels):
            # update: now 1-lepton events are only needed for 1-lepton channel (before they were also needed for 0-lepton channel for tt CR)
            elif (self._b("isWmunu")[0] or self._b("isWenu")[0]) and not ("Wln" in self.channels):
                return False
            elif (self._b("isZnn")[0]) and not ("Znn" in self.channels):
                return False
            self.cutFlow[3] += 1
            
            # VECTOR BOSON
            Vphi_syst = {}
            any_v_sysvar_passes = False
            if self._b("isZee")[0] or self._b("isZmm")[0]:
                lep1 = ROOT.TLorentzVector()
                lep2 = ROOT.TLorentzVector()
                i1 = self._b("vLidx")[0]
                i2 = self._b("vLidx")[1]
                if self._b("isZee")[0]:
                    lep1.SetPtEtaPhiM(tree.Electron_pt[i1], tree.Electron_eta[i1], tree.Electron_phi[i1], tree.Electron_mass[i1])
                    lep2.SetPtEtaPhiM(tree.Electron_pt[i2], tree.Electron_eta[i2], tree.Electron_phi[i2], tree.Electron_mass[i2])
                elif self._b("isZmm")[0]:
                    lep1.SetPtEtaPhiM(tree.Muon_pt[i1], tree.Muon_eta[i1], tree.Muon_phi[i1], tree.Muon_mass[i1])
                    lep2.SetPtEtaPhiM(tree.Muon_pt[i2], tree.Muon_eta[i2], tree.Muon_phi[i2], tree.Muon_mass[i2])
                V = lep1 + lep2
            elif self._b("isWenu")[0] or self._b("isWmunu")[0]:
                i1 = self._b("vLidx")[0]
                if self._b("isWenu")[0]:
                    sel_lepton_pt   = tree.Electron_pt[i1]
                    sel_lepton_eta  = tree.Electron_eta[i1]
                    sel_lepton_phi  = tree.Electron_phi[i1]
                    sel_lepton_mass = tree.Electron_mass[i1]
                elif self._b("isWmunu")[0]:
                    sel_lepton_pt   = tree.Muon_pt[i1]
                    sel_lepton_eta  = tree.Muon_eta[i1]
                    sel_lepton_phi  = tree.Muon_phi[i1]
                    sel_lepton_mass = tree.Muon_mass[i1]
                self._b("lepMetDPhi")[0] = abs(ROOT.TVector2.Phi_mpi_pi(sel_lepton_phi - tree.MET_Phi))
                #if self._b("lepMetDPhi")[0] > 2.0:
                #    return False

                MET = ROOT.TLorentzVector()
                Lep = ROOT.TLorentzVector()
                MET.SetPtEtaPhiM(tree.MET_Pt, 0.0, tree.MET_Phi, 0.0)
                Lep.SetPtEtaPhiM(sel_lepton_pt, sel_lepton_eta, sel_lepton_phi, sel_lepton_mass)
                cosPhi12 = (Lep.Px()*MET.Px() + Lep.Py()*MET.Py()) / (Lep.Pt() * MET.Pt())
                self._b("V_mt")[0] = ROOT.TMath.Sqrt(2*Lep.Pt()*MET.Pt() * (1 - cosPhi12))

                V = MET + Lep

                if self.sample.isMC():
                    MET_syst = ROOT.TLorentzVector()
                    for syst in self.METsystematics:
                        for UD in self._variations(syst):
                            MET_syst.SetPtEtaPhiM(getattr(tree, 'MET_pt_'+syst+UD), 0.0, getattr(tree, 'MET_phi_'+syst+UD), 0.0)
                            V_syst = MET_syst + Lep
                            Vphi_syst[syst+UD] = V_syst.Phi()
                            if V_syst.Pt() > 150.0:
                                any_v_sysvar_passes = True

            elif self._b("isZnn")[0]:
                MET = ROOT.TLorentzVector()
                MET.SetPtEtaPhiM(tree.MET_Pt, 0.0, tree.MET_Phi, 0.0)
                if self.sample.isMC():
                    for syst in self.METsystematics:
                        for UD in self._variations(syst):
                            Vphi_syst[syst+UD] = getattr(tree, 'MET_phi_'+syst+UD)
                V = MET
            else:
                self.count("fail_lepton_selection")
                return False


            self._b("V_pt")[0] = V.Pt()
            self._b("V_eta")[0] = V.Eta()
            self._b("V_phi")[0] = V.Phi()
            self._b("V_mass")[0] = V.M()
            
            if (self._b("isZee")[0] or self._b("isZmm")[0]) and self._b("V_pt")[0] < 50.0:
                return False
            elif self._b("isZnn")[0] and self._b("V_pt")[0] < 150.0:
                return False
            elif (self._b("isWenu")[0] or self._b("isWmunu")[0]) and (self._b("V_pt")[0] < 150.0 and not any_v_sysvar_passes):
                return False
            
            self.cutFlow[4] += 1


            # JETS
            if self._b("isZnn")[0]:
                j1ptCut = 35.0
                j2ptCut = 35.0
                maxJetPtCut = 60.0
                j1BtagName = 'loose'
            elif self._b("isWmunu")[0] or self._b("isWenu")[0]:
                j1ptCut = 25.0
                j2ptCut = 25.0
                maxJetPtCut = 0.0
                j1BtagName = 'loose'
            elif self._b("isZmm")[0] or self._b("isZee")[0]:
                j1ptCut = 20.0
                j2ptCut = 20.0
                maxJetPtCut = 0.0
                j1BtagName = 'none'
            else:
                return False
            j2BtagName = 'none'

            # btagMaxCut can overwrite the default b-tag max cut
            if self.btagMaxCut is not None and len(self.btagMaxCut) > 0:
                j1BtagName = self.btagMaxCut

            any_taggerPassed_syst = False
            for jetDefinition in self.jetDefinitions:

                indexName  = jetDefinition['indexName'] if 'indexName' in jetDefinition else ''
                taggerName = jetDefinition['taggerName']
                j1Btag     = self.btagWPs[self.year][taggerName][j1BtagName]
                j2Btag     = self.btagWPs[self.year][taggerName][j2BtagName]

                jetTagger  = getattr(tree, taggerName)
                taggerSelection = lambda x: (
                        len(x) == 2 and 
                        jetTagger[x[0]] >= j1Btag and
                        jetTagger[x[1]] >= j2Btag 
                    )

                # -> nominal
                selectedJets = self.HighestTaggerValueBJets(tree, j1ptCut, j2ptCut, self.taggerName, puIdCut=self.puIdCut, jetIdCut=self.jetIdCut, maxJetPtCut=maxJetPtCut, systList=None)
                taggerPassed = taggerSelection(selectedJets) 
                if taggerPassed:
                    self._b("hJidx"+indexName)[0] = selectedJets[0]
                    self._b("hJidx"+indexName)[1] = selectedJets[1]
                    any_taggerPassed_syst         = True
                    self.count("pass_tagger_"+taggerName)
                else:
                    if debugEvent:
                        print("DEBUG-EVENT: did not pass jet-selection.")
                        print("DEBUG-EVENT: selected jets:", selectedJets)
                        print("DEBUG-EVENT: passed btag cuts:", taggerPassed)
                    self._b("hJidx"+indexName)[0] = -1
                    self._b("hJidx"+indexName)[1] = -1

                # -> systematic variations
                if self.sample.isMC():
                    systList     = ["".join(x) for x in itertools.product(self.systematics, ["Up", "Down"])]
                    nSystematics = len(systList)

                    # this returns a list of lists of jet indices for all systematic variations
                    selectedJets = self.HighestTaggerValueBJets(tree, j1ptCut, j2ptCut, self.taggerName, puIdCut=self.puIdCut, jetIdCut=self.jetIdCut, maxJetPtCut=maxJetPtCut,  systList=systList)

                    # apply pre-selection on tagger
                    for j in range(nSystematics):
                        hJidxName_syst = "hJidx{suffix}_{syst}".format(suffix=indexName, syst=systList[j])
                        if taggerSelection(selectedJets[j]):
                            self._b(hJidxName_syst)[0] = selectedJets[j][0]
                            self._b(hJidxName_syst)[1] = selectedJets[j][1]
                            any_taggerPassed_syst      = True
                        else:
                            self._b(hJidxName_syst)[0] = -1
                            self._b(hJidxName_syst)[1] = -1 

            self.cutFlow[5] += 1

            #print tree.Hbb_fjidx,tree.FatJet_pt[tree.Hbb_fjidx] if(tree.Hbb_fjidx>-1) else "outside",V.Pt()

            # BOOSTED JET selection
            any_boostedJet_syst = False

            # nominal
            Hbb_fjidx = self.HighestTaggerValueFatJet(tree, ptCut=250.0, msdCut=50.0, etaCut=2.5, taggerName='FatJet_deepTagMD_bbvsLight', systList=None, Vphi=V.Phi())
            self._b("Hbb_fjidx")[0] = Hbb_fjidx
            if Hbb_fjidx > -1:
                any_boostedJet_syst = True

            # systematics
            if self.sample.isMC():
                systList     = ["".join(x) for x in itertools.product(self.systematicsBoosted, ["Up", "Down"])]
                nSystematics = len(systList)

                # return list of jet indices for systematic variations
                Hbb_fjidx_syst = self.HighestTaggerValueFatJet(tree, ptCut=250.0, msdCut=50.0, etaCut=2.5, taggerName='FatJet_deepTagMD_bbvsLight', systList=systList, Vphi=V.Phi(), Vphi_syst=Vphi_syst)
                for j in range(nSystematics):
                    fJidxName_syst = "Hbb_fjidx_{syst}".format(syst=systList[j])
                    self._b(fJidxName_syst)[0] = Hbb_fjidx_syst[j]
                    if Hbb_fjidx_syst[j] > -1:
                        any_boostedJet_syst = True

            #boostedPass = (tree.Hbb_fjidx>-1 and tree.FatJet_Pt[tree.Hbb_fjidx]>250 and V.Pt()>250)  #AC: selection for boosted analysis
            boostedPass = (any_boostedJet_syst and V.Pt()>250)

            if boostedPass:
                self.count("pass_boosted")

            # discard event if none of the jet selections finds two Higgs candidate jets
            #if not (boostedPass or defaultTaggerPassed or any([jets['isSelected'] for jets in self.jetDefinitions])):
            if not (boostedPass or any_taggerPassed_syst): 
                return False

            self.cutFlow[6] += 1


            # yield in the end
            self.cutFlow[7] += 1
            if debugEvent:
                print "DEBUG-EVENT: event passed!!!"

        return True

    def afterProcessing(self):
        super(VHbbSelection, self).afterProcessing()

        print "cut-flow:"
        print "  beginning          ", self.cutFlow[0]
        print "  HLT                ", self.cutFlow[1]
        print "  Leptons            ", self.cutFlow[2]
        print "  Channel            ", self.cutFlow[3]
        print "  Vector boson       ", self.cutFlow[4]
        print "  Jets               ", self.cutFlow[5]
        print "  Vpt                ", self.cutFlow[6]
        print "  end                ", self.cutFlow[7]

        print "efficiency:", (1.0*self.cutFlow[7]/self.cutFlow[0]) if self.cutFlow[0] > 0 else "-"

        #print "jet selections:", self.jetDefinitions


