#!/usr/bin/env python
import ROOT
from BranchTools import Collection
from BranchTools import AddCollectionsModule
import array
import numpy as np

# create branches for gen level particles with HEPPY naming (only for those used)
class HeppyStyleGen(AddCollectionsModule):

    def __init__(self, debug=False):
        self.debug = debug
        super(HeppyStyleGen, self).__init__()
        self.version = 4

    def customInit(self, initVars):
        self.sampleTree = initVars['sampleTree']
        self.isData = initVars['sample'].isData()
        if not self.isData:
            self.addCollection(Collection('GenVbosons', ['pt','pdgId','GenPartIdx'], maxSize=40))
            self.addCollection(Collection('GenTop', ['pt','GenPartIdx'], maxSize=4))
            self.addCollection(Collection('GenHiggsBoson', ['pt','GenPartIdx'], maxSize=4))
        
            self.branchBuffers['VtypeSim'] = array.array('i', [0])
            self.branches.append({'name': 'VtypeSim', 'formula': self.getBranch, 'arguments': 'VtypeSim', 'type': 'i'})

            self.addCollection(Collection('GenBs',['pt','eta','phi','genPartIdx'], maxSize=20))
            self.addVectorBranch("GenJetAK8_nBhadrons", default=0, branchType='i', length=100, leaflist="GenJetAK8_nBhadrons[nGenJetAK8]/i")
            self.addVectorBranch("GenJetAK8_nBhadrons2p4", default=0, branchType='i', length=100, leaflist="GenJetAK8_nBhadrons2p4[nGenJetAK8]/i")
            self.addVectorBranch("GenJet_nBhadrons", default=0, branchType='i', length=100, leaflist="GenJet_nBhadrons[nGenJet]/i")
            self.addVectorBranch("GenJet_nBhadrons2p4", default=0, branchType='i', length=100, leaflist="GenJet_nBhadrons2p4[nGenJet]/i")

            #Sum$(GenBs_pt>25&&abs(GenBs_eta)<2.6)
            self.addIntegerBranch("nGenBpt25eta2p6")
            self.addIntegerBranch("nGenBpt20eta2p6")
    
    def processEvent(self, tree):
        if not self.hasBeenProcessed(tree):
            self.markProcessed(tree)

            if not self.isData:
                higgsParticles = []
                vParticles = []
                tParticles = []

                genParticles = self.sampleTree.getCollection('GenPart')
                for idx in range(genParticles.size()):
                    if genParticles.pdgId[idx] == 25 and genParticles.genPartIdxMother[idx] == 0:
                        higgsParticles.append({'pt': genParticles.pt[idx], 'GenPartIdx': idx})
                    if genParticles.pdgId[idx] in [23, 24, -24]:
                        vParticles.append({'pt': genParticles.pt[idx], 'pdgId': genParticles.pdgId[idx], 'daughters': [], 'GenPartIdx': idx})
                    if genParticles.pdgId[idx] in [-6, 6] and genParticles.genPartIdxMother[idx] == 0:
                        tParticles.append({'pt': genParticles.pt[idx], 'GenPartIdx': idx})

                    for vParticle in vParticles:
                        if genParticles.genPartIdxMother[idx] == vParticle['GenPartIdx']:
                            vParticle['daughters'].append({'pdgId': genParticles.pdgId[idx], 'GenPartIdx': idx})


                # only keep last V which decays into leptons
                oldLen = len(vParticles)
                vParticles = [x for x in vParticles if len(x['daughters'])>0 and x['pdgId']!=x['daughters'][0]['pdgId']]
               
                # write collections
                self.collections['GenHiggsBoson'].fromList(higgsParticles)
                self.collections['GenVbosons'].fromList(vParticles)
                self.collections['GenTop'].fromList(tParticles)

                # get VtypeSim
                VtypeSim = -1
                if len(vParticles) == 1:
                    vParticle = vParticles[0]
                    if vParticle['pdgId'] == 23 and len(vParticle['daughters'])>1 and abs(vParticle['daughters'][0]['pdgId'])==abs(vParticle['daughters'][1]['pdgId']):
                        if abs(vParticle['daughters'][0]['pdgId']) == 11:
                            VtypeSim = 1
                        elif abs(vParticle['daughters'][0]['pdgId']) == 13:
                            VtypeSim = 0
                        elif abs(vParticle['daughters'][0]['pdgId']) == 15:
                            VtypeSim = 5
                        elif abs(vParticle['daughters'][0]['pdgId']) in [12,14,16]:
                            VtypeSim = 4
                    elif abs(vParticle['pdgId']) == 24 and len(vParticle['daughters'])==2:
                        if (abs(vParticle['daughters'][0]['pdgId']) == 11 and abs(vParticle['daughters'][1]['pdgId']) == 12) or (abs(vParticle['daughters'][0]['pdgId']) == 12 and abs(vParticle['daughters'][1]['pdgId']) == 11):
                            VtypeSim = 3
                        elif (abs(vParticle['daughters'][0]['pdgId']) == 13 and abs(vParticle['daughters'][1]['pdgId']) == 14) or (abs(vParticle['daughters'][0]['pdgId']) == 14 and abs(vParticle['daughters'][1]['pdgId']) == 13):
                            VtypeSim = 2
                elif len(vParticles) == 2:
                    VtypeSim = -2
                self.branchBuffers['VtypeSim'][0] = VtypeSim

                # decay chains
                mothers  = [genParticles.genPartIdxMother[idx] for idx in range(genParticles.size()) ]
                products = [idx for idx in range(genParticles.size()) if idx not in mothers] 

                chains = []
                for product in products:
                    # find last B
                    idx = product
                    foundB = False
                    while not foundB:
                        if ((((abs(genParticles.pdgId[idx]) // 100) % 10) == 5) or (((abs(genParticles.pdgId[idx]) // 1000) % 10) == 5)) and genParticles.status[idx] == 2:
                            foundB = True
                        else:
                            if genParticles.genPartIdxMother[idx] < 0:
                                break
                            idx = genParticles.genPartIdxMother[idx]
                    if foundB:
                        chain = [idx]
                        while True:
                            idx = genParticles.genPartIdxMother[idx]
                            if idx < 0:
                                break
                            chain.append(idx)
                        chains.append(chain)

                # filter duplicates
                unique_chains = [list(x) for x in set(tuple(x) for x in chains)] 
                delIndices = []
                for k in range(len(unique_chains)):
                    for j in range(len(unique_chains)):
                        if j != k:
                            if len(unique_chains[k]) > len(unique_chains[j]) and unique_chains[k][-len(unique_chains[j]):] == unique_chains[j]:
                                delIndices.append(j)
                delIndices = list(set(delIndices))
                delIndices.sort(reverse=True)
                for k in delIndices:
                    del unique_chains[k]

                # fill collections
                bHadrons = []
                for unique_chain in unique_chains:
                    bHadrons.append({'pt': genParticles.pt[unique_chain[0]], 'eta': genParticles.eta[unique_chain[0]], 'phi': genParticles.phi[unique_chain[0]], 'genPartIdx': unique_chain[0]})
                self.collections['GenBs'].fromList(bHadrons)

                #GenJetAK8_nBhadrons counting
                bHadronsMatched = []
                for i in range(tree.nGenJetAK8):
                    nB = 0
                    nB2p4 = 0
                    for j in range(len(bHadrons)):
                        dPhi = ROOT.TVector2.Phi_mpi_pi(tree.GenJetAK8_phi[i]-bHadrons[j]['phi'])
                        dEta = tree.GenJetAK8_eta[i] - bHadrons[j]['eta']
                        dR = np.sqrt(dPhi*dPhi + dEta*dEta)
                        if dR < 0.8 and j not in bHadronsMatched:
                            nB += 1
                            bHadronsMatched.append(j)
                            if abs(bHadrons[j]['eta']) < 2.4:
                                nB2p4 += 1
                    self._b("GenJetAK8_nBhadrons")[i] = nB
                    self._b("GenJetAK8_nBhadrons2p4")[i] = nB2p4
                
                #GenJet_nBhadrons counting
                bHadronsMatched = []
                for i in range(tree.nGenJet):
                    nB = 0
                    nB2p4 = 0
                    for j in range(len(bHadrons)):
                        dPhi = ROOT.TVector2.Phi_mpi_pi(tree.GenJet_phi[i]-bHadrons[j]['phi'])
                        dEta = tree.GenJet_eta[i] - bHadrons[j]['eta']
                        dR = np.sqrt(dPhi*dPhi + dEta*dEta)
                        if dR < 0.8 and j not in bHadronsMatched:
                            nB += 1
                            bHadronsMatched.append(j)
                            if abs(bHadrons[j]['eta']) < 2.4:
                                nB2p4 += 1
                    self._b("GenJet_nBhadrons")[i] = nB
                    self._b("GenJet_nBhadrons2p4")[i] = nB2p4

                self._b("nGenBpt25eta2p6")[0] = len([x for x in bHadrons if x['pt']>25 and abs(x['eta'])<2.6])
                self._b("nGenBpt20eta2p6")[0] = len([x for x in bHadrons if x['pt']>20 and abs(x['eta'])<2.6])

        return True

