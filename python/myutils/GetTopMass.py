# Copied from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/TopQuarkAnalysis/SingleTop/src/TopProducer.cc?revision=1.9&view=markup
# and converted from C++ to Python

import numpy as np
import array
from pdgId import pdgId
import sys
import os
from ROOT import *
from math import pow
from math import sqrt
from math import acos
from math import cos
from math import sin
import sys
import numpy as np
from BranchTools import Collection
from BranchTools import AddCollectionsModule

class GetTopMass(AddCollectionsModule):

    def __init__(self, sample=None, nano=False, propagateJES = False, METmethod = 2, useHJets = 0, minbTag = 0.5, branchName='top_mass', addBoostSystematics=False, addMinMax=False, puIdCut=6, jetIdCut=4):
        super(GetTopMass, self).__init__()
        self.version = 2
        self.nano = nano
        self.lastEntry = -1
        self.branchName = branchName
        self.METmethod = METmethod #Check the getJetPtMass function for the different methods
        self.useHJets = useHJets # Using only H jets if useHJets == 1
        self.minbTag = minbTag
        self.propagateJES = propagateJES
        self.addBoostSystematics = addBoostSystematics
        self.addMinMax = addMinMax
        self.puIdCut = puIdCut
        self.jetIdCut = jetIdCut

    def dependencies(self):
        return {'VHbbSelection': 5}

    def customInit(self, initVars):
        self.sample = initVars['sample']
        self.config = initVars['config']

        # Branch names from config
        self.brJetPtReg = "Jet_PtReg"
        self.brJetPtNom = "Jet_Pt"

        self.Jet_btag = self.config.get('General', 'Jet_btag')
        self.tagidx = self.config.get('General', 'hJidx')

        if self.config.has_option('General', 'eIdx') and  self.config.has_option('General', 'muIdx'):
            self.eIdx = self.config.get('General', 'eIdx')
            self.muIdx = self.config.get('General', 'muIdx')
        else:
            print "WARNING: eIdx and muIdx not defined in [General]! Using default lepton index: \'vLidx\' "
            self.eIdx = "vLidx"
            self.muIdx = "vLidx"

        self.dataset = self.config.get('General','dataset')

        self.METpt = 'MET_Pt'
        self.METphi = 'MET_Phi'

        if self.sample.isMC():
            self.brJetMassNom = "Jet_mass_nom"
        else:
            self.brJetMassNom = "Jet_mass"

        # nominal
        self.systematics = [None]

        ##########
        # add JES/JER systematics
        ##########
        # only defined for nano at the moment

        if self.propagateJES and self.nano:
            self.jetSystematics = ['jer','jerReg','jesAbsoluteStat','jesAbsoluteScale','jesAbsoluteFlavMap','jesAbsoluteMPFBias','jesFragmentation','jesSinglePionECAL','jesSinglePionHCAL','jesFlavorQCD','jesRelativeJEREC1','jesRelativeJEREC2','jesRelativeJERHF','jesRelativePtBB','jesRelativePtEC1','jesRelativePtEC2','jesRelativePtHF','jesRelativeBal','jesRelativeFSR','jesRelativeStatFSR','jesRelativeStatEC','jesRelativeStatHF','jesPileUpDataMC','jesPileUpPtRef','jesPileUpPtBB','jesPileUpPtEC1','jesPileUpPtEC2','jesPileUpPtHF','jesPileUpMuZero','jesPileUpEnvelope','jesTotal','unclustEn'] 

            if self.addBoostSystematics:
                self.jetSystematics+= ['jms']
                self.jetSystematics+= ['jmr']
            if self.sample.isMC():
                self.systematics += self.jetSystematics

        for syst in self.systematics:
            for Q in self._variations(syst):
                self.addBranch(self._v(self.branchName, syst, Q))
                self.addBranch(self._v(self.branchName + '_noReg', syst, Q))

    def processEvent(self, tree):

        if not self.hasBeenProcessed(tree):
            self.markProcessed(tree)

            # leptons/MET
            treeMETpt, treeMETphi = self.getMET(tree)
            lep = TLorentzVector()
            met = TLorentzVector()
            if not self.nano:
                lep.SetPtEtaPhiM(tree.vLeptons_new_pt[0], tree.vLeptons_new_eta[0], tree.vLeptons_new_phi[0], tree.vLeptons_new_mass[0])
                met.SetPtEtaPhiM(tree.met_pt, tree.met_eta, tree.met_phi, tree.met_mass)
            else:
                if tree.Vtype == 2:
                    mu1Idx = getattr(tree,self.muIdx)[0]
                    lep.SetPtEtaPhiM(tree.Muon_pt[mu1Idx], tree.Muon_eta[mu1Idx], tree.Muon_phi[mu1Idx], tree.Muon_mass[mu1Idx])
                if tree.Vtype == 3:
                    e1Idx = getattr(tree,self.eIdx)[0]
                    lep.SetPtEtaPhiM(tree.Electron_pt[e1Idx], tree.Electron_eta[e1Idx], tree.Electron_phi[e1Idx], tree.Electron_mass[e1Idx])
                met.SetPtEtaPhiM(treeMETpt, 0, treeMETphi, 0)

            # compute top mass for all variations
            for syst in self.systematics:
                for UD in self._variations(syst):

                    # MET
                    treeMETpt, treeMETphi = self.getMET(tree, syst, UD)
                    met.SetPtEtaPhiM(treeMETpt, 0, treeMETphi, 0)

                    # jets
                    Jet_PtReg, Jet_MassReg, Jet_Pt = self.getJetPtMasses(tree, syst=syst, UD=UD)

                    # closest jets
                    cJidx, cHJidx = self.closestJetIdx(tree, lep, syst=syst, UD=UD, Jet_PtReg=Jet_PtReg, Jet_MassReg=Jet_MassReg, Jet_Pt=Jet_Pt)
                    closestIdx    = cJidx if self.useHJets == 0 else cHJidx

                    if closestIdx != - 1:
                       cJet = TLorentzVector()
                       cJet.SetPtEtaPhiM(Jet_PtReg[closestIdx], tree.Jet_eta[closestIdx], tree.Jet_phi[closestIdx], Jet_MassReg[closestIdx])

                       top_mass = self.computeTopMass_new(lep, met, cJet, self.METmethod)

                       cJet.SetPtEtaPhiM(tree.Jet_pt[closestIdx], tree.Jet_eta[closestIdx], tree.Jet_phi[closestIdx], tree.Jet_mass[closestIdx])

                       top_mass_noReg = self.computeTopMass_new(lep, met, cJet, self.METmethod)

                    elif closestIdx == -1:
                       top_mass = -99
                       top_mass_noReg = -99

                    self._b(self._v(self.branchName, syst, UD))[0] = top_mass
                    self._b(self._v(self.branchName + '_noReg', syst, UD))[0] = top_mass_noReg

            return True

    def getResolvedJetIndicesFromTree(self, tree, syst=None, UD=None):
        indexNameSyst = (self.tagidx + '_' + syst + UD) if not self._isnominal(syst) else self.tagidx
        if hasattr(tree, indexNameSyst):
            self.count('_debug_resolved_idx_syst_exists')
            return getattr(tree, indexNameSyst)
        else:
            self.count('_debug_resolved_idx_fallback_nom')
            self.count('_debug_resolved_idx_fallback_nom_for_'+str(syst))
            return getattr(tree, self.tagidx)

    ##
    # \Function EquationSolver:
    #
    # Solves 3rd degree equations
    #
    # \Author A. Orso M. Iorio
    #
    #
    # \version  $Id: EquationSolver.h,v 1.1 2013/02/27 12:18:42 degrutto Exp $
    #
    def EquationSolve(self, a, b, c, d):


      result = []



      if (a != 0):

        q = (3*a*c-b*b)/(9*a*a)
        r = (9*a*b*c - 27*a*a*d - 2*b*b*b)/(54*a*a*a)
        Delta = q*q*q + r*r

        rho=0.
        theta=0.

        if( Delta<=0):
          rho = sqrt(-(q*q*q))

          theta = acos(r/rho)

          s = complex(sqrt(-q)*cos(theta/3.0),sqrt(-q)*sin(theta/3.0))
          t = complex(sqrt(-q)*cos(-theta/3.0),sqrt(-q)*sin(-theta/3.0))

        if(Delta>0):
          #print r, sqrt(Delta)
          if (r+sqrt(Delta) > 0):
              s = complex(pow((r+sqrt(Delta)),(1./3)),0)
          else:
              s = complex(-pow(abs((r+sqrt(Delta))),(1./3)),0)
          if (r-sqrt(Delta) > 0):
              t = complex(pow((r-sqrt(Delta)),(1./3)),0)
          else:
              t = complex(-pow(abs((r-sqrt(Delta))),(1./3)),0)

        i = complex(0.,1.0)

        x1 = s+t+complex(-b/(3.0*a),0)
        x2 = (s+t)*complex(-0.5,0)-complex(b/(3.0*a),0)+(s-t)*i*complex(sqrt(3)/2.0,0)
        x3 = (s+t)*complex(-0.5,0)-complex(b/(3.0*a),0)-(s-t)*i*complex(sqrt(3)/2.0,0)

        if(abs(x1.imag)<0.0001): result.append(x1.real)
        if(abs(x2.imag)<0.0001): result.append(x2.real)
        if(abs(x3.imag)<0.0001): result.append(x3.real)

        #print x1,x2,x3
        return result
      else:
          return result


      return result

    def getNu4Momentum(self, TLepton, TMET, debug=False):

      branch = -1

      Lepton = TLorentzVector()
      Lepton.SetPxPyPzE(TLepton.Px(), TLepton.Py(), TLepton.Pz(), TLepton.E());
      MET = TLorentzVector()
      MET.SetPxPyPzE(TMET.Px(), TMET.Py(), 0., TMET.E());

      mW = 80.38;

      result = []

      MisET2 = (MET.Px()*MET.Px() + MET.Py()*MET.Py());
      mu = (mW*mW)/2 + MET.Px()*Lepton.Px() + MET.Py()*Lepton.Py();
      a  = (mu*Lepton.Pz())/(Lepton.Energy()*Lepton.Energy() - Lepton.Pz()*Lepton.Pz());
      a2 = pow(a,2);
      b  = (pow(Lepton.Energy(),2.)*(MisET2) - pow(mu,2.))/(pow(Lepton.Energy(),2) - pow(Lepton.Pz(),2));
      pz1 = 0.
      pz2 = 0.
      pznu = 0.
      nNuSol = 0

      p4nu_rec = TLorentzVector()
      p4W_rec = TLorentzVector()
      p4b_rec = TLorentzVector()
      p4Top_rec = TLorentzVector()
      p4lep_rec = TLorentzVector()

      p4lep_rec.SetPxPyPzE(Lepton.Px(),Lepton.Py(),Lepton.Pz(),Lepton.Energy());

      #print a2,b
      if(a2-b > 0 ):
        root = sqrt(a2-b);
        pz1 = a + root;
        pz2 = a - root;
        nNuSol = 2;

        pznu = pz1;
        branch = 1
        if abs(pz2)<abs(pz1):
            pznu = pz2
            branch = 2

        Enu = sqrt(MisET2 + pznu*pznu);

        p4nu_rec.SetPxPyPzE(MET.Px(), MET.Py(), pznu, Enu);

        result.append(p4nu_rec);

      else:

        ptlep = Lepton.Pt()
        pxlep=Lepton.Px()
        pylep=Lepton.Py()
        metpx=MET.Px()
        metpy=MET.Py()

        EquationA = 1.
        EquationB = -3.*pylep*mW/(ptlep)
        EquationC = mW*mW*(2*pylep*pylep)/(ptlep*ptlep)+mW*mW-4*pxlep*pxlep*pxlep*metpx/(ptlep*ptlep)-4*pxlep*pxlep*pylep*metpy/(ptlep*ptlep)
        EquationD = 4.*pxlep*pxlep*mW*metpy/(ptlep)-pylep*mW*mW*mW/ptlep

        solutions = self.EquationSolve(EquationA,EquationB,EquationC,EquationD)

        solutions2 = self.EquationSolve(EquationA,-EquationB,EquationC,-EquationD);

        deltaMin = 14000*14000
        zeroValue = -mW*mW/(4*pxlep)
        minPx=0
        minPy=0

        for i in range(len(solutions)):
            if(solutions[i]<0 ): continue
            p_x = (solutions[i]*solutions[i]-mW*mW)/(4*pxlep)
            p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x -mW*ptlep*solutions[i])/(2*pxlep*pxlep)
            Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy)

            if(Delta2< deltaMin and Delta2 > 0):
                deltaMin = Delta2
                minPx=p_x
                minPy=p_y
                branch = 3

        for i in range(len(solutions2)):
            if(solutions2[i]<0 ): continue
            p_x = (solutions2[i]*solutions2[i]-mW*mW)/(4*pxlep)
            p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x +mW*ptlep*solutions2[i])/(2*pxlep*pxlep)
            Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy)
            if(Delta2< deltaMin and Delta2 > 0):
                deltaMin = Delta2
                minPx=p_x
                minPy=p_y
                branch = 4

        pyZeroValue= ( mW*mW*pxlep + 2*pxlep*pylep*zeroValue)
        delta2ZeroValue= (zeroValue-metpx)*(zeroValue-metpx) + (pyZeroValue-metpy)*(pyZeroValue-metpy)

        if(deltaMin==14000*14000):
            if debug:
                branch = 6
                return TLorentzVector(0,0,0,0), branch
            else:
                return TLorentzVector(0,0,0,0)

        if(delta2ZeroValue < deltaMin):
            deltaMin = delta2ZeroValue
            minPx=zeroValue
            minPy=pyZeroValue
            branch = 5

        mu_Minimum = (mW*mW)/2 + minPx*pxlep + minPy*pylep
        a_Minimum  = (mu_Minimum*Lepton.Pz())/(Lepton.Energy()*Lepton.Energy() - Lepton.Pz()*Lepton.Pz())
        pznu = a_Minimum

        Enu = sqrt(minPx*minPx+minPy*minPy + pznu*pznu)
        p4nu_rec.SetPxPyPzE(minPx, minPy, pznu , Enu)
        result.append(p4nu_rec)

      if debug:
          return result[0], branch
      else:
          return result[0]

#    def computeTopMass(self, lep, met, jets, jetsIdx):
    def computeTopMass(self, lep, met, jets):

        #neutrino = self.getNu4Momentum(lep, met)
        bjet = TLorentzVector()
        minDR = 99
        closeIdx = -1
        for i, jet in enumerate(jets):
            if (jet.Pt() > 30): # and jet.bTagCSV > CSVL and jet.puID > 0 and jet.Id > 0? how?
                dR = jet.DeltaR(lep)
                if (dR < minDR):
                    minDR = dR
                    bjet = jet
#                    closeIdx = jetsIdx[i]
        if (bjet.Pt() <= 0):
            return -99
#        print "Old closest idx: {}".format(closeIdx)
        top = lep + met + bjet
        return top.M()

    def getJetPtMasses(self, tree, syst=None, UD=None):

        JetPtReg = np.array(getattr(tree, self.brJetPtReg))
        JetPtNom = np.array(getattr(tree,self.brJetPtNom))
        JetMassNom = np.array(getattr(tree,self.brJetMassNom))

        if self._isnominal(syst) or syst.startswith('unclustEn'):
           JetPt      = JetPtReg
           JetMass    = JetMassNom * JetPtReg /JetPtNom
        elif 'jerReg' in syst:
           JetPtSys   = np.array(getattr(tree,"Jet_PtReg" + UD))
           JetPt      = JetPtSys
           JetMass    = JetMassNom * JetPtSys / JetPtNom
        else:
           JetPtSys   = np.array(getattr(tree,"Jet_pt" + "_" + syst + UD))
           JetMassSys = np.array(getattr(tree, "Jet_mass" + "_" + syst + UD))
           JetPt      = JetPtSys * JetPtReg / JetPtNom
           JetMass    = JetMassSys * JetPtReg / JetPtNom
        return JetPt, JetMass, JetPtNom


#    def getJetPtMass(self, tree, i, variation=None):
#
#        JetPtReg = getattr(tree, self.brJetPtReg)[i]
#        JetPtNom = getattr(tree,self.brJetPtNom)[i] 
#        JetMassNom = getattr(tree,self.brJetMassNom)[i] 
#
#        if variation is None or variation.startswith('unclustEn'):
#           JetPt =  JetPtReg
#           JetMass = JetMassNom * JetPtReg /JetPtNom
#
#        elif 'jerReg' in variation:
#           Q = "Up" if "Up" in variation else "Down"
#           JetPtSys = getattr(tree,"Jet_PtReg" + Q)[i]
#
#           JetPt =  JetPtSys
#           JetMass = JetMassNom * JetPtSys / JetPtNom
#
#        else:
#           JetPtSys = getattr(tree,"Jet_pt" + "_" + variation)[i]
#           JetMassSys = getattr(tree, "Jet_mass" + "_" + variation)[i] 
#
#           JetPt =  JetPtSys * JetPtReg / JetPtNom
#           JetMass = JetMassSys * JetPtReg / JetPtNom
#           print "{}: {}, {}: {}".format("Jet_pt" + "_" + variation,JetPtSys,"Jet_mass" + "_" + variation,JetMassSys)
#           print "JetPt {}, JetMass: {}".format(JetPt,JetMass)
#        return JetPt, JetMass

    def getMET(self, tree, syst=None, UD=None):
        if not self._isnominal(syst) and not syst.startswith('jerReg'):
            treeMETpt = getattr(tree, "MET_pt_" + syst + UD)
            treeMETphi = getattr(tree, "MET_phi_" + syst + UD)
        else:
            treeMETpt = getattr(tree, self.METpt)
            treeMETphi = getattr(tree, self.METphi)
        return treeMETpt, treeMETphi

    def closestJetIdx(self, tree, lep, syst=None, UD=None, Jet_PtReg=None, Jet_MassReg=None, Jet_Pt=None):
        #neutrino = self.getNu4Momentum(lep, met)

        minDR = 999
        minHDR = 999

        ClosestJIdx = -1
        ClosestHJIdx = -1

        # reuse pt/mass already obtained from tree if available
        # if one is missing, recompute them again
        if Jet_PtReg is None or Jet_MassReg is None or Jet_Pt is None:
            Jet_PtReg_r, Jet_MassReg_r, Jet_Pt_r= self.getJetPtMasses(tree, syst, UD)
            Jet_PtReg   = Jet_PtReg_r if Jet_PtReg is None else Jet_PtReg
            Jet_MassReg = Jet_MassReg_r if Jet_MassReg is None else Jet_MassReg
            Jet_Pt      = Jet_Pt_r if Jet_Pt is None else Jet_Pt

        # b-tagger
        jetBtags = getattr(tree, self.Jet_btag)

        Jet_lepFilter = tree.Jet_lepFilter
        Jet_jetId     = tree.Jet_jetId
        Jet_puId      = tree.Jet_puId
        Jet_eta       = tree.Jet_eta
        Jet_phi       = tree.Jet_phi
        nJet          = tree.nJet

        # selected jet indices
        hJidx0, hJidx1 = self.getResolvedJetIndicesFromTree(tree, syst, UD)

        for i in range(nJet):
           #if temPt < 30 or tembTag < self.minbTag or not (tree.Jet_lepFilter): continue

           if Jet_PtReg[i] >= 30 and jetBtags[i] >= self.minbTag and tree.Jet_lepFilter[i] and tree.Jet_jetId[i] > self.jetIdCut and (tree.Jet_puId[i]>self.puIdCut or Jet_Pt[i] >50.0):

               tempJet = TLorentzVector()
               tempJet.SetPtEtaPhiM(Jet_PtReg[i],tree.Jet_eta[i],tree.Jet_phi[i], Jet_MassReg[i])
               dR = tempJet.DeltaR(lep)
    #           print "Jetidx: {}, dR = {}".format(i,dR)

               if dR < minDR:
                  minDR = dR
                  ClosestJIdx = i

               if (i == hJidx0 or i == hJidx1) and dR < minHDR:
                  minHDR = dR
                  ClosestHJIdx = i

#        print "hJidx0 = {}, hJidx1 = {}".format(self.hJidx0,self.hJidx1)
        return ClosestJIdx, ClosestHJIdx

            
    def computeTopMass_new(self, lep, met, jet, METmethod):
        #neutrino = self.getNu4Momentum(lep, met)

        new_top = TLorentzVector()
        jet_transverse = TLorentzVector()
        lep_transverse = TLorentzVector()

        if METmethod == 1:
           jet_transverse.SetPxPyPzE(jet.Px(),jet.Py(),0,sqrt((jet.M())**2 + (jet.Pt())**2))   
           lep_transverse.SetPxPyPzE(lep.Px(),lep.Py(),0,sqrt((lep.M())**2 + (lep.Pt())**2))   

           new_top = jet_transverse + lep_transverse + met
        elif METmethod == 2:
           neutrino = TLorentzVector()
           neutrino = self.getNu4Momentum(lep, met)
           new_top = lep + jet + neutrino
        elif METmethod == 3:
           new_top = lep + jet + met

        return new_top.M()


