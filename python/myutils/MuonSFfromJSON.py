from __future__ import division
import json
import os
import array
from JsonTable import JsonTable
from vLeptons import vLeptonSelector
from BranchTools import Collection
from BranchTools import AddCollectionsModule
import math

class MuonSFfromJSON(AddCollectionsModule):

    def __init__(self, jsonFiles=None, branchName='muonSF', channel='None', year=None):
        super(MuonSFfromJSON, self).__init__()
        self.jsonFiles = jsonFiles
        self.debug = 'XBBDEBUG' in os.environ
        self.branchName = branchName

        self.year = year
        self.tablenames = {
        2016:  {'Zll':{ 
                    'idSf':'SF_MuIDTightBCDEF',
                    'isoSf':'NUM_LooseRelIso_DEN_LooseID',
                    'eff_IsoMu8_mc':'Mu8LegMC',
                    'eff_IsoMu17_mc':'Mu17LegMC',
                    'eff_IsoMu8_data':'Mu8Leg',
                    'eff_IsoMu17_data':'Mu17Leg'
                    },
                'Wlv':{
                    'idSf':'SF_MuIDTight', # important to add BCDEF / GH for Run dependent lumi
                    'isoSf':'SF_MuIsoTight',
                    'triggerSf':'SF_MuTrigger'},
                'lumi':{
                    'BCDEF': 20.1/36.4,
                    'GH': 16.3/36.4}
                },
        2017:  {'Zll':{
                    'idSf':'NUM_LooseID_DEN_genTracks',
                    'isoSf':'NUM_LooseRelIso_DEN_LooseID',
                    'eff_IsoMu8_mc':'Mu8LegMC',
                    'eff_IsoMu17_mc':'Mu17LegMC',
                    'eff_IsoMu8_data':'Mu8Leg',
                    'eff_IsoMu17_data':'Mu17Leg'
                    },
                'Wlv':{
                    'idSf':'NUM_TightID_DEN_genTracks',
                    'isoSf':'NUM_UltraTightIso4_DEN_TightIDandIPCut',
                    'triggerSf':'NUM_IsoMu27_DEN_empty'}
                },
        2018:   {'Zll':{
                    'idSf':'NUM_LooseID_DEN_TrackerMuons2018',
                    'isoSf':'NUM_LooseRelIso_DEN_LooseID2018',
                    'eff_IsoMu8_mc':'Mu8LegMC',
                    'eff_IsoMu17_mc':'Mu17LegMC',
                    'eff_IsoMu8_data':'Mu8Leg',
                    'eff_IsoMu17_data':'Mu17Leg'},
                'Wlv':{
                    'idSf':'TightID',
                    'isoSf':'TightISO',
                    'triggerSf':'TRGIsoMu24'}
                }
        }

        # load JOSN files
        self.jsonTable = JsonTable(jsonFiles)
        self.channel = channel 
        if self.channel== 'Zll':
            self.idSf = self.jsonTable.getEtaPtTable(self.tablenames[self.year][self.channel]['idSf'], 'abseta_pt')
            self.isoSf = self.jsonTable.getEtaPtTable(self.tablenames[self.year][self.channel]['isoSf'], 'abseta_pt')
            self.eff_IsoMu8_mc = self.jsonTable.getEtaPtTable(self.tablenames[self.year][self.channel]['eff_IsoMu8_mc'], 'abseta_pt') 
            self.eff_IsoMu17_mc = self.jsonTable.getEtaPtTable(self.tablenames[self.year][self.channel]['eff_IsoMu17_mc'], 'abseta_pt') 
            self.eff_IsoMu8_data = self.jsonTable.getEtaPtTable(self.tablenames[self.year][self.channel]['eff_IsoMu8_data'], 'abseta_pt') 
            self.eff_IsoMu17_data = self.jsonTable.getEtaPtTable(self.tablenames[self.year][self.channel]['eff_IsoMu17_data'], 'abseta_pt') 
            #self.triggerSf = ; not implemented yet
        elif self.channel == 'Wlv':
            if self.year == 2016:
                for corr in self.tablenames[self.year][self.channel].keys():
                    for lumi_sec in self.tablenames[self.year]['lumi'].keys():
                        setattr(self, corr+lumi_sec, self.jsonTable.getEtaPtTable(self.tablenames[self.year][self.channel][corr]+lumi_sec, 'abseta_pt')) 
            else:
                for corr in self.tablenames[self.year][self.channel].keys():
                    setattr(self, corr, self.jsonTable.getEtaPtTable(self.tablenames[self.year][self.channel][corr], 'abseta_pt'))

        else: 
            print "Channel not defined!"
            raise Exception("ChannelNotDefined")

        self.systVariations = [None, 'Down', 'Up']

    def customInit(self, initVars):
        sample = initVars['sample']
        self.isData = sample.type == 'DATA'
        self.config = initVars['config']

        # prepare buffers for new branches to add
        self.branches = []
        self.branchBuffers = {}
        self.lastEntry = -1
        if not self.isData:
            for n in ['', '_Id', '_Iso', '_trigger']:
                self.addVectorBranch(self.branchName + n, default=1.0, length=3)

    def processEvent(self, tree):
        # if current entry has not been processed yet
        if not self.hasBeenProcessed(tree) and not self.isData:
            self.markProcessed(tree)

            #TODO check the assigned weights for events with more than 2 letpns
            Vtype = tree.Vtype
            vLidx = []
            lep_eta = []
            lep_pt = []
         
            if self.channel == "Wlv" and Vtype == 2:
               vLidx = [tree.vLidx[0]]
               lep_pt = [tree.Muon_pt[vLidx[0]]]
               lep_eta = [tree.Muon_eta[vLidx[0]]]
            elif self.channel == "Zll" and Vtype == 0:
               vLidx = [tree.vLidx[0],tree.vLidx[1]]
               lep_pt = [tree.Muon_pt[vLidx[0]], tree.Muon_pt[vLidx[1]]]
               lep_eta = [tree.Muon_eta[vLidx[0]], tree.Muon_eta[vLidx[1]]]

            self.computeSF(
                    weight_trigg=self._b(self.branchName + '_trigger'),
                    weight_Id=self._b(self.branchName + '_Id'),
                    weight_Iso=self._b(self.branchName + '_Iso'),
                    weight_SF=self._b(self.branchName),
                    lep_eta=lep_eta,
                    lep_pt=lep_pt,
                    lep_n=len(vLidx)
                    )

        return True

    def computeSF(self, weight_trigg, weight_Id, weight_Iso, weight_SF, lep_eta, lep_pt, lep_n):
        '''Computes the trigger, IdIso (including separated variations in eta) and final event SF'''
        # require two electrons

        if lep_n == 1 or lep_n == 2:
        #Calculating the trigger and IdIso weights and Down/Up variations
            for i, syst in enumerate(self.systVariations):
                if self.year == 2016:
                    weight_trigg[i] = 0.0
                    weight_Id[i] = 0.0
                    weight_Iso[i] = 0.0

                    for corr in self.tablenames[self.year][self.channel].keys():
                        for lumi_sec in self.tablenames[self.year]['lumi'].keys():
                            w = self.getSf(corr+lumi_sec, lep_eta, lep_pt, lep_n, syst=syst)
                            if corr.startswith("triggerSf"):
                                weight_trigg[i] += self.tablenames[self.year]['lumi'][lumi_sec] * w
                            if corr.startswith("idSf"):
                                weight_Id[i] += self.tablenames[self.year]['lumi'][lumi_sec] * w
                            if corr.startswith("isoSf"):
                                weight_Iso[i] += self.tablenames[self.year]['lumi'][lumi_sec] * w

                else:
                    weight_trigg[i] = self.getTriggerSf(lep_eta, lep_pt, lep_n, syst=syst)
                    weight_Id[i] = self.getSf("idSf", lep_eta, lep_pt, lep_n)
                    weight_Iso[i] = self.getSf("isoSf", lep_eta, lep_pt, lep_n)
                weight_SF[i] = weight_trigg[i] * weight_Id[i] * weight_Iso[i]

            #print "trigg", weight_trigg
            #print "Id", weight_Id
            #print "Iso", weight_Iso
            #print "SF", weight_SF

        #This is when an event has 0 or more than 2 lepton
        else:
            for i, syst in enumerate(self.systVariations):
                weight_trigg[i] = 1.
                weight_Id[i] = 1.
                weight_Iso[i] = 1.
                weight_SF[i] = 1.

    def getSf(self, corr, eta, pt, len_n, syst=None):
        SF = 1.
        for i in range(len_n):
            SF = SF * self.jsonTable.find(getattr(self, corr), eta[i], pt[i], syst=syst)
        return SF


#    def getTriggerSf(self, eta1, pt1, eta2, pt2, syst=None):
#        leg1 = 1.0 #not implemented yet
#        leg2 = 1.0 
#        #define efficiency for MC and 
#        effData_leg8 = []
#        effData_leg17= []
#        effMC_leg8= []
#        effMC_leg17 = []
#        #if self.debug:
#        #    print "leg1: eta:", eta1, " pt:", pt1, "->", leg1
#        #    print "leg2: eta:", eta2, " pt:", pt2, "->", leg2
#        #    print "-->", leg1*leg2
#        return leg1*leg2

    def computeEventSF_fromleg(self, eff_mu8_l1, eff_mu17_l1, eff_mu8_l2, eff_mu17_l2):
        #returns event efficiency and relative uncertainty
        eff_event = [1.,0.]
        eff_event[0] = ((eff_mu8_l1[0]**2 * eff_mu17_l2[0] + eff_mu8_l2[0]**2 * eff_mu17_l1[0])/(eff_mu8_l1[0] + eff_mu8_l2[0]))
        #relative uncertainty down
        #print('nominal', eff_event[0])
        uncert_down = (abs(((eff_mu8_l1[1]**2 * eff_mu17_l2[1] + eff_mu8_l2[1]**2 * eff_mu17_l1[1])/(eff_mu8_l1[1] + eff_mu8_l2[1])) - eff_event[0])/eff_event[0])
        #print('down',(abs(((eff_mu8_l1[0]**2 * eff_mu17_l2[0] + eff_mu8_l2[0]**2 * eff_mu17_l1[0])/(eff_mu8_l1[0] + eff_mu8_l2[0])) - eff_event[0])))
        #relative uncertainty up
        uncert_up = (abs(((eff_mu8_l1[2]**2 * eff_mu17_l2[2] + eff_mu8_l2[2]**2 * eff_mu17_l1[2])/(eff_mu8_l1[2] + eff_mu8_l2[2])) - eff_event[0])/eff_event[0])
        eff_event[1]  = (uncert_down+uncert_up)/2.
        return eff_event


    def getTriggerSf(self, eta, pt, len_n, syst=None):
        triggSF = 1.
        if self.channel == 'Zll':
            if len_n==2:
                eff_mu8_l1_mc  = self.jsonTable.findvalerr(self.eff_IsoMu8_mc, eta[0], pt[0], self.systVariations)
                eff_mu17_l1_mc  = self.jsonTable.findvalerr(self.eff_IsoMu17_mc, eta[0], pt[0], self.systVariations)
                eff_mu8_l1_data = self.jsonTable.findvalerr(self.eff_IsoMu8_data, eta[0], pt[0], self.systVariations)
                eff_mu17_l1_data = self.jsonTable.findvalerr(self.eff_IsoMu17_data, eta[0], pt[0], self.systVariations)

                eff_mu8_l2_mc  = self.jsonTable.findvalerr(self.eff_IsoMu8_mc, eta[1], pt[1], self.systVariations)
                eff_mu17_l2_mc  = self.jsonTable.findvalerr(self.eff_IsoMu17_mc, eta[1], pt[1], self.systVariations)
                eff_mu8_l2_data = self.jsonTable.findvalerr(self.eff_IsoMu8_data, eta[1], pt[1], self.systVariations)
                eff_mu17_l2_data = self.jsonTable.findvalerr(self.eff_IsoMu17_data, eta[1], pt[1], self.systVariations)
                #print("new event triggerSf")
                #print(eta[0], pt[0])
                #print(eff_mu8_l1_mc)
                #print(eff_mu8_l1_data)

                EffData = self.computeEventSF_fromleg(eff_mu8_l1_data,eff_mu17_l1_data,eff_mu8_l2_data,eff_mu17_l2_data)
                EffMC = self.computeEventSF_fromleg(eff_mu8_l1_mc,eff_mu17_l1_mc,eff_mu8_l2_mc,eff_mu17_l2_mc)

                #print(EffData)
                #print(EffMC)
                
                if syst==None:
                    triggSF = (EffData[0]/EffMC[0]) 
                elif syst=='Down':
                    triggsf = (EffData[0]/EffMC[0])
                    triggSF = (1-math.sqrt(EffData[1]**2+ EffMC[1]**2))*triggsf
                elif syst=='Up':
                    triggsf = (EffData[0]/EffMC[0])
                    triggSF = (1+math.sqrt(EffData[1]**2+ EffMC[1]**2))*triggsf

        elif self.channel == 'Wlv': 
           for i in range(len_n):
                triggSF = triggSF * self.jsonTable.find(self.triggerSf[i], eta[i], pt[i], syst=syst)
            #print "leg", i, " : eta:", eta[i], " pt:", pt[i], "->", self.jsonTable.find(self.triggerSf[i], eta[i], pt[i], syst=syst)
        return triggSF

if __name__ == "__main__":
    sfObject = MuonSFfromJSON([
            '/work/krgedia/CMSSW_10_1_0/src/Xbb/python/data/Zll/Muons/RunBCDEF_SF_ID.json',
            '/work/krgedia/CMSSW_10_1_0/src/Xbb/python/data/Zll/Muons/RunBCDEF_SF_ISO.json',
        ],channel="Wlv",year=2017)
    #print sfObject.getIdSf([0.5],[121] ,1)
    #print sfObject.getIdSf(-0.5, 42)
    #print sfObject.getIsoSf(1.5, 21)
    #print sfObject.getIsoSf(-1.5, 21)

    sfObject = MuonSFfromJSON([
            '/work/krgedia/CMSSW_10_1_0/src/Xbb/python/data/Zll/Muons/Eff_IsoMu17Cut_NUM_IsoMu17leg_DEN_LooseRelIso_PAR_newpt_etaMC2018.json',
            '/work/krgedia/CMSSW_10_1_0/src/Xbb/python/data/Zll/Muons/Eff_IsoMu17Cut_NUM_IsoMu17leg_DEN_LooseRelIso_PAR_newpt_eta_DATA2018.json',
            '/work/krgedia/CMSSW_10_1_0/src/Xbb/python/data/Zll/Muons/Eff_IsoMu8Cut_NUM_IsoMu8leg_DEN_LooseRelIso_PAR_newpt_etaMC2018.json',
            '/work/krgedia/CMSSW_10_1_0/src/Xbb/python/data/Zll/Muons/Eff_IsoMu8Cut_NUM_IsoMu8leg_DEN_LooseRelIso_PAR_newpt_eta_DATA2018.json',
            '/work/krgedia/CMSSW_10_1_0/src/Xbb/python/data/Zll/Muons/RunABCD2018_SF_ID.json'
        ],channel="Zll",year=2018)

    #print sfObject.getIdSf([-1.86865234375,-1.86865234375],[121.0,65.0] , 2, None)
    #print sfObject.getTriggerSf([-1.86865234375,0.5],[121.0,121] , 2, None)
    #print sfObject.getTriggerSf([-1.86865234375,0.5],[121.0,121] , 2, 'Up')
    #print sfObject.getTriggerSf([-1.86865234375,0.5],[121.0,121] , 2, 'Down')
    #print sfObject.getTriggerSf([0.5,0.5],[65,85] , 2, 'Up')
    #print sfObject.getTriggerSf([0.5,0.5],[65,85] , 2, 'Down')



