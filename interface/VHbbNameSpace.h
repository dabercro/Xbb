#include "TLorentzVector.h"
#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"
#include <iostream>


namespace VHbb {
    
#include "XYMETCorrection.h"

    double LUT5(int n, float n0, float n1, float n2, float n3, float n4, float n5) {
        if (n < 1)
            return n0;
        else if (n == 1)
            return n1;
        else if (n == 2)
            return n2;
        else if (n == 3)
            return n3;
        else if (n == 4)
            return n4;
        else if (n > 4)
            return n5;
    }

    double LUT8(int n, float n0, float n1, float n2, float n3, float n4, float n5, float n6, float n7) {
        if (n < 1)
            return n0;
        else if (n == 1)
            return n1;
        else if (n == 2)
            return n2;
        else if (n == 3)
            return n3;
        else if (n == 4)
            return n4;
        else if (n == 5)
            return n5;
        else if (n == 6)
            return n6;
        else if (n > 6)
            return n7;
    }

double DYBW_F100(int nGenJet) {
    //{3, 0.99}, {4, 1.08}, {5, 1.21}, {6, 1.24}, {7, 1.17}, {8, 1.09}
    if (nGenJet > 7) nGenJet = 8;
    if (nGenJet < 3) nGenJet = 3;
    if (nGenJet == 3)
        return 0.99;
    else if (nGenJet == 4)
        return 1.08;
    else if (nGenJet == 5)
        return 1.21;
    else if (nGenJet == 6)
        return 1.24;
    else if (nGenJet == 7)
        return 1.17;
    else
        return 1.08;
}

double DYBW_F200(int nGenJet) {
    //{3, 1.27}, {4, 1.31}, {5, 1.28}, {6, 1.34}, {7, 1.24}, {8, 1.15}
    if (nGenJet > 7) nGenJet = 8;
    if (nGenJet < 3) nGenJet = 3;
    if (nGenJet == 3)
        return 1.27;
    else if (nGenJet == 4)
        return 1.31;
    else if (nGenJet == 5)
        return 1.28;
    else if (nGenJet == 6)
        return 1.34;
    else if (nGenJet == 7)
        return 1.24;
    else
        return 1.15;
}

double DYBW_B100(int nGenJet) {
    //{3, 1.23}, {4, 1.33}, {5, 1.41}, {6, 1.37}, {7, 1.24}, {8, 1.10}
    if (nGenJet > 7) nGenJet = 8;
    if (nGenJet < 3) nGenJet = 3;
    if (nGenJet == 3)
        return 1.23;
    else if (nGenJet == 4)
        return 1.33;
    else if (nGenJet == 5)
        return 1.41;
    else if (nGenJet == 6)
        return 1.37;
    else if (nGenJet == 7)
        return 1.24;
    else
        return 1.10;
}

double DYBW_B200(int nGenJet) {
    //{3, 1.46}, {4, 1.51}, {5, 1.61}, {6, 1.47}, {7, 1.41}, {8, 1.25}
    if (nGenJet > 7) nGenJet = 8;
    if (nGenJet < 3) nGenJet = 3;
    if (nGenJet == 3)
        return 1.46;
    else if (nGenJet == 4)
        return 1.51;
    else if (nGenJet == 5)
        return 1.61;
    else if (nGenJet == 6)
        return 1.47;
    else if (nGenJet == 7)
        return 1.41;
    else
        return 1.25;
}

double VptNLOweightDY(int nGenJets, int nB, double Vpt) {
    double SF=1.0;
    if (Vpt > 500) Vpt=500.0;
    if (Vpt < 50.0) Vpt=50.0;
    if (nGenJets==0&&nB==0) SF=1.41667257-0.001230725*Vpt;
    if (nGenJets==1&&nB==0) SF=1.393552626-0.001233994*Vpt;
    if (nGenJets==2&&nB==0) SF=1.347529697-0.001160407*Vpt;
    if (nGenJets==3&&nB==0) SF=1.071661786-0.000555197*Vpt;
    if (nGenJets==4&&nB==0) SF=0.729464542-2.05154e-05*Vpt;
    if (nGenJets==5&&nB==0) SF=0.679199746-0.000232446*Vpt;
    if (nGenJets==0&&nB==1) SF=1.393949937-0.001378941*Vpt;
    if (nGenJets==1&&nB==1) SF=1.424499504-0.001661213*Vpt;
    if (nGenJets==2&&nB==1) SF=1.288121414-0.001039063*Vpt;
    if (nGenJets==3&&nB==1) SF=1.090142309-0.000854487*Vpt;
    if (nGenJets==4&&nB==1) SF=0.861623061-0.000975446*Vpt;
    if (nGenJets==5&&nB==1) SF=0.376231418+0.000346819*Vpt;
    if (nGenJets==0&&nB==2) SF=1.681719504-0.004419806*Vpt;
    if (nGenJets==1&&nB==2) SF=1.546154404-0.002363307*Vpt;
    if (nGenJets==2&&nB==2) SF=1.453219102-0.002123416*Vpt;
    if (nGenJets==3&&nB==2) SF=1.078354952-0.000859408*Vpt;
    if (nGenJets==4&&nB==2) SF=0.949781804-0.001050688*Vpt;
    if (nGenJets==5&&nB==2) SF=0.709471838-0.000694619*Vpt;    
    if (SF < 0.0) SF=0;
    return SF;
}


double SoverSBWeight(double BDT, int channel) {

    //Check to which bin the BDT belongs to
    int BDTbin = -1;
    if (BDT <= -0.8)                    BDTbin = 0;
    if (BDT > -0.8 && BDT <= -0.68)     BDTbin = 0;
    if (BDT > -0.68 && BDT <= -0.56)    BDTbin = 1;
    if (BDT > -0.56 && BDT <= -0.44)    BDTbin = 2;
    if (BDT > -0.44 && BDT <= -0.32)    BDTbin = 3;
    if (BDT > -0.32 && BDT <= -0.02)    BDTbin = 4;
    if (BDT > -0.2 && BDT <= -0.08)     BDTbin = 5;
    if (BDT > -0.08 && BDT <= 0.04)     BDTbin = 6;
    if (BDT > 0.04 && BDT <= 0.16)      BDTbin = 7;
    if (BDT > 0.16 && BDT <= 0.28)      BDTbin = 8;
    if (BDT > 0.28 && BDT <= 0.4)       BDTbin = 9;
    if (BDT > 0.4 && BDT <= 0.52)       BDTbin = 10;
    if (BDT > 0.52 && BDT <= 0.64)      BDTbin = 11;
    if (BDT > 0.64 && BDT <= 0.76)      BDTbin = 12;
    if (BDT > 0.76 && BDT <= 0.88)      BDTbin = 13;
    if (BDT > 0.88)                     BDTbin = 14;

    double SB_Zuu_low[15]={2.76737387847e-05,
                           6.64750589345e-05,
                           0.000149447397438,
                           0.000364720697321,
                           0.000704082928445,
                           0.00114156054762,
                           0.00178676530053,
                           0.00289993913949,
                           0.00473330069259,
                           0.00714558106403,
                           0.0109339081032,
                           0.0184934984159,
                           0.0326417501844,
                           0.0664331803746,
                           0.155709437265
    };

    double SB_Zee_low[15]={4.38820281278e-05,
                           6.59406637058e-05,
                           0.000151812163876,
                           0.00037181791516,
                           0.000651510194471,
                           0.00109292029164,
                           0.00180371097475,
                           0.00289255480983,
                           0.00440456389377 ,
                           0.00717194663279,
                           0.011191506056,
                           0.0192027774896,
                           0.0313549416529,
                           0.0627899701193,
                           0.188193463981
    };

    double SB_Zuu_high[15]={0.000150118338319,
                            0.000177184418437,
                            0.000381785086736,
                            0.00081621546275,
                            0.00132284030661,
                            0.00315465582298,
                            0.00513913238438,
                            0.00871858305194,
                            0.014172742012,
                            0.0265281704166,
                            0.0408929588652,
                            0.0762144976363,
                            0.138083647382,
                            0.204843885382,
                            0.321478566847
    };

    double SB_Zee_high[15]={7.252018604e-08,
                            0.000204092490389,
                            0.00045852950862,
                            0.000812207870653,
                            0.00187548258053,
                            0.00255126940987,
                            0.00515682235312,
                            0.00931094133369,
                            0.0140552030406,
                            0.0230148836113,
                            0.0413999559234,
                            0.0668556717614,
                            0.131102613516,
                            0.242663806491,
                            0.331429498934
    };

    if (channel == 1) return SB_Zee_low[BDTbin];
    if (channel == 2) return SB_Zuu_low[BDTbin];
    if (channel == 3) return SB_Zee_high[BDTbin];
    if (channel == 4) return SB_Zuu_high[BDTbin];

    return 0;

    }



  double deltaPhi(double phi1, double phi2) {
    double result = phi1 - phi2;
    if (result > TMath::Pi()) {
      result -= 2*TMath::Pi();
    } else if (result <= -TMath::Pi()) {
      result += 2*TMath::Pi();
    }
    return result;
  }

  double deltaR(double eta1, double phi1, double eta2, double phi2) {
    double deta = eta1 - eta2;
    double dphi = deltaPhi(phi1, phi2);
    return TMath::Sqrt(deta*deta + dphi*dphi);
  }

  double HJetVarByCSV(double CSVj1, double Varj1, double CSVj2, double Varj2, bool pickhighestBtag) {
    if ((CSVj1 > CSVj2) && pickhighestBtag) {
      return Varj1;
    } else {
      return Varj2;
    }
  }

  double Hmass(double V_eta, double V_phi, double V_pt,
               double hJet1_eta, double hJet1_phi, double hJet1_pt,
               double hJet2_eta, double hJet2_phi, double hJet2_pt) {

    TVector3 V(1, 1, 1);
    V.SetPtEtaPhi(V_pt, V_eta, V_phi);

    TVector3 H1(1, 1, 1);
    H1.SetPtEtaPhi(hJet1_pt, hJet1_eta, hJet1_phi);
    H1.SetMag(1 / sin(H1.Theta()));

    TVector3 H2(1, 1, 1);
    H2.SetPtEtaPhi(hJet2_pt, hJet2_eta, hJet2_phi);
    H2.SetMag(1 / sin(H2.Theta()));

    TVector3 n1(H1), n2(H2);

    double det = n1.Px()*n2.Py() - n2.Px()*n1.Py();

    H1.SetMag((-n2.Py()*V.Px() + n2.Px()*V.Py()) / (sin(n1.Theta())*det));
    H2.SetMag((n1.Py()*V.Px() - n1.Px()*V.Py()) / (sin(n2.Theta())*det));

    double mass = TMath::Sqrt(TMath::Power(H1.Mag() + H2.Mag(), 2) - TMath::Power((H1+H2).Mag(), 2));

    return mass;
  }

  double Hmass_comb(double hJet1_eta, double hJet1_phi, double hJet1_pt, double hJet1_mass,
		    double hJet2_eta, double hJet2_phi, double hJet2_pt, double hJet2_mass) {

    TLorentzVector H1, H2;

    H1.SetPtEtaPhiM(hJet1_pt, hJet1_eta, hJet1_phi, hJet1_mass);
    H2.SetPtEtaPhiM(hJet2_pt, hJet2_eta, hJet2_phi, hJet2_mass);

    return (H1+H2).M();
  }

  double Hmass_3j(double h_eta, double h_phi, double h_pt, double h_mass,
		  double aJet_eta, double aJet_phi, double aJet_pt, double aJet_mass) {

    TLorentzVector H, H3;

    H.SetPtEtaPhiM(h_pt, h_eta, h_phi, h_mass);
    H3.SetPtEtaPhiM(aJet_pt, aJet_eta, aJet_phi, aJet_mass);

    return (H+H3).M();
  }

  double GetEnergy(double pt, double eta, double mass) {
    TLorentzVector m1;
    m1.SetPtEtaPhiM(pt, eta, 0, mass);
    return m1.Energy();
  }

  double HVMass(double pt, double eta, double phi, double mass, double pt2, double eta2, double phi2, double mass2) {
    TLorentzVector m1, m2, msum;
    m1.SetPtEtaPhiM(pt, eta, phi, mass);
    m2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);
    msum = m1 + m2;
    return msum.M();
  }

  double ANGLELZ(double pt, double eta, double phi, double mass, double pt2, double eta2, double phi2, double mass2) {
    TLorentzVector m1, m2, msum;
    m1.SetPtEtaPhiM(pt, eta, phi, mass);
    m2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);
    msum = m1 + m2;

    TVector3 bZ =  msum.BoostVector();
    m1.Boost(-bZ);
    m2.Boost(-bZ);

    TVector3 b1;
    if (int(pt) % 2) {
      b1 =  m2.BoostVector();
    } else {
      b1 =  m1.BoostVector();
    }

    double cosTheta = b1.Dot(msum.BoostVector()) / (b1.Mag() * msum.BoostVector().Mag());

    return cosTheta;
  }

  double ANGLEHB(double pt, double eta, double phi, double e, double pt2, double eta2, double phi2, double e2) {
    TLorentzVector m1, m2, msum;
    m1.SetPtEtaPhiE(pt, eta, phi, e);
    m2.SetPtEtaPhiE(pt2, eta2, phi2, e2);
    msum = m1 + m2;

    TVector3 bZ =  msum.BoostVector();
    m1.Boost(-bZ);
    m2.Boost(-bZ);

    TVector3 b1;
    if (int(pt) % 2) {
      b1 =  m2.BoostVector();
    } else {
      b1 =  m1.BoostVector();
    }

    double cosTheta = b1.Dot(msum.BoostVector()) / (b1.Mag() * msum.BoostVector().Mag());

    return cosTheta;
  }

  double ANGLEHBWithM(double pt, double eta, double phi, double mass, double pt2, double eta2, double phi2, double mass2) {
    TLorentzVector m1, m2, msum;
    m1.SetPtEtaPhiM(pt, eta, phi, mass);
    m2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);
    msum = m1 + m2;

    TVector3 bZ =  msum.BoostVector();
    m1.Boost(-bZ);
    m2.Boost(-bZ);

    TVector3 b1;
    if (int(pt) % 2) {
      b1 = m2.BoostVector();
    } else {
      b1 = m1.BoostVector();
    }

    double cosTheta = b1.Dot(msum.BoostVector()) / (b1.Mag() * msum.BoostVector().Mag());

    return cosTheta;
  }

  double metCorSysShift(double met, double metphi, int Nvtx, int EVENT_run) {
    double metx = met * cos(metphi);
    double mety = met * sin(metphi);
    double px = 0.0, py = 0.0;

    if (EVENT_run != 1) {
      //pfMEtSysShiftCorrParameters_2012runAplusBvsNvtx_data
      px = +1.68804e-01 + 3.37139e-01*Nvtx;
      py = -1.72555e-01 - 1.79594e-01*Nvtx;
    } else {
      //pfMEtSysShiftCorrParameters_2012runAplusBvsNvtx_mc
      px = +2.22335e-02 - 6.59183e-02*Nvtx;
      py = +1.52720e-01 - 1.28052e-01*Nvtx;
    }

    metx -= px;
    mety -= py;

    return std::sqrt(metx*metx + mety*mety);
  }

  double metphiCorSysShift(double met, double metphi, int Nvtx, int EVENT_run) {
    double metx = met * cos(metphi);
    double mety = met * sin(metphi);
    double px = 0.0, py = 0.0;

    if (EVENT_run != 1) {
      //pfMEtSysShiftCorrParameters_2012runAplusBvsNvtx_data
      px = +1.68804e-01 + 3.37139e-01*Nvtx;
      py = -1.72555e-01 - 1.79594e-01*Nvtx;
    } else {
      //pfMEtSysShiftCorrParameters_2012runAplusBvsNvtx_mc
      px = +2.22335e-02 - 6.59183e-02*Nvtx;
      py = +1.52720e-01 - 1.28052e-01*Nvtx;
    }

    metx -= px;
    mety -= py;

    if (metx == 0.0 && mety == 0.0) return 0.0;

    double phi1 = std::atan2(mety, metx);
    double phi2 = std::atan2(mety, metx) - 2.0*M_PI;

    if (std::abs(phi1 - metphi) < std::abs(phi2 - metphi) + 0.5*M_PI) {
      return phi1;
    } else {
      return phi2;
    }
  }

  TVector2 metType1Reg(double met, double metphi, double corr1, double corr2,
                       double pt1, double eta1, double phi1, double e1,
                       double pt2, double eta2, double phi2, double e2) {

    double metx = met * cos(metphi);
    double mety = met * sin(metphi);

    TLorentzVector j1, j2;
    j1.SetPtEtaPhiE(pt1, eta1, phi1, e1);
    j2.SetPtEtaPhiE(pt2, eta2, phi2, e2);

    metx += (j1.Px() * (1-corr1)) + (j2.Px() * (1-corr2));
    mety += (j1.Py() * (1-corr1)) + (j2.Py() * (1-corr2));

    TVector2 corrMET(metx, mety);

    return corrMET;
  }

  double metType1Phi(double met, double metphi, double corr1, double corr2,
                     double pt1, double eta1, double phi1, double e1,
                     double pt2, double eta2, double phi2, double e2) {
    return metType1Reg(met, metphi, corr1, corr2, pt1, eta1, phi1, e1, pt2, eta2, phi2, e2).Phi();
  }

  double metType1Et(double met, double metphi, double corr1, double corr2,
                    double pt1, double eta1, double phi1, double e1,
                    double pt2, double eta2, double phi2, double e2) {
    return metType1Reg(met, metphi, corr1, corr2, pt1, eta1, phi1, e1, pt2, eta2, phi2, e2).Mod();
  }

  double met_MPF(double met, double metphi, double pt, double phi) {
    return 1. + (met * pt * std::cos(deltaPhi(metphi, phi)) / (pt*pt));
  }

  double resolutionBias(double eta) {
    if (eta < 0.5) return 0.052;
    if (eta < 1.1) return 0.057;
    if (eta < 1.7) return 0.096;
    if (eta < 2.3) return 0.134;
    if (eta < 5) return 0.28;
    // Nominal
    return 0;
  }

  double evalJERBias(double ptreco, double ptgen, double eta1) {
    double eta = fabs(eta1);
    double cor = 1;

    // Limit the effect to the core.
    if ((fabs(ptreco-ptgen) / ptreco) < 0.5) {
      cor = (ptreco + resolutionBias(eta) * (ptreco-ptgen)) / ptreco;
    }

    if (ptgen > 0.) {
      return ptreco * cor;
    } else {
      return ptreco;
    }
  }

  double evalEt(double pt, double eta, double phi, double e) {
    TLorentzVector j;
    j.SetPtEtaPhiE(pt, eta, phi, e);
    return j.Et();
  }

  double evalEtFromPtEtaPhiM(double pt, double eta, double phi, double m) {
    TLorentzVector j;
    j.SetPtEtaPhiM(pt, eta, phi, m);
    return j.Et();
  }

  double evalMt(double pt, double eta, double phi, double e) {
    TLorentzVector j;
    j.SetPtEtaPhiE(pt, eta, phi, e);
    return j.Mt();
  }

  double evalMtFromPtEtaPhiM(double pt, double eta, double phi, double m) {
    TLorentzVector j;
    j.SetPtEtaPhiM(pt, eta, phi, m);
    return j.Mt();
  }

  double weight_QCD(int nGenHiggsBoson, int nGenTop, int nGenVbosons, double lheHT, int GenVbosons_pdgId) {
    // QCD reweighting applied to W+Jets samples.
    if (nGenHiggsBoson==0 && nGenTop==0 && nGenVbosons==1 && lheHT>100 && abs(GenVbosons_pdgId)==24) {
      return (lheHT>100 && lheHT<200)*(1.459/1.21) + (lheHT>200 && lheHT<400)*(1.434/1.21) + (lheHT>400 && lheHT<600)*(1.532/1.21) + (lheHT>600)*(1.004/1.21);
    }
    return 1.;
  }

  double weight_EWK(int nGenHiggsBoson, int nGenTop, int nGenVbosons, double GenVbosons_pt, int VtypeSim, int GenVbosons_pdgId) {
    // EWK reweighting applied to W+Jets and Z+Jets samples.
    if (nGenHiggsBoson==0 && nGenTop==0 && nGenVbosons==1 && GenVbosons_pt>100 && GenVbosons_pt<3000) {
      if ((VtypeSim==0 || VtypeSim==1 || VtypeSim==4 || VtypeSim==5) && GenVbosons_pdgId==23) {
        return -0.1808051 + 6.04146*TMath::Power(GenVbosons_pt + 759.098, -0.242556);
      } else if ((VtypeSim==2 || VtypeSim==3) && GenVbosons_pdgId==24) {
        return -0.830041 + 7.93714*TMath::Power(GenVbosons_pt + 877.978, -0.213831);
      }
    }
    return 1.;
  }

  double weight_LOtoNLO(int nGenHiggsBoson, int nGenTop, int nGenVbosons, int GenVbosons_pdgId, double deta_jj, int nGenBJets) {
    // LO to NLO reweighting applied to Drell-Yan samples.
    if (nGenHiggsBoson==0 && nGenTop==0 && nGenVbosons==1 && GenVbosons_pdgId==23 && deta_jj<5) {
      if (nGenBJets < 1) {
        return 0.935422 + 0.0403162*deta_jj - 0.0089026*deta_jj*deta_jj + 0.0064324*deta_jj*deta_jj*deta_jj - 0.000212443*deta_jj*deta_jj*deta_jj*deta_jj;
      } else if (nGenBJets == 1) {
        return 0.962415 + 0.0329463*deta_jj - 0.0414479*deta_jj*deta_jj + 0.0240993*deta_jj*deta_jj*deta_jj - 0.00278271*deta_jj*deta_jj*deta_jj*deta_jj;
      } else if (nGenBJets >= 2) {
        return (0.721265 - 0.105643*deta_jj - 0.0206835*deta_jj*deta_jj + 0.00558626*deta_jj*deta_jj*deta_jj)*TMath::Exp(0.450244*deta_jj);
      }
    }
    return 1.;
  }

  double weight_TTbar_genPt(int nGenTop, double GenTop_pt1, double GenTop_pt2) {
    if (nGenTop == 2) {
      double sf_top1 = exp(0.0615 - 0.0005*GenTop_pt1);
      double sf_top2 = exp(0.0615 - 0.0005*GenTop_pt2);
      return sqrt(sf_top1 * sf_top2);
    }
    return 1.;
  }

  double weight_TTbar_nJetCentral(int nGenTop, int nJetsCentral) {
    if (nGenTop == 2) {
      double sf_njets = 1.2217 * (0.924773 - 0.0212496*nJetsCentral);
      return sf_njets;
    }
    return 1.;
  }

   double weight_TTbar_nJetCentral_Up(int nGenTop, int nJetsCentral) {
    if (nGenTop == 2) {
      double sf_njets = 1.2222 * (0.9058043 - 0.01751612*nJetsCentral);
      return sf_njets;
    }
    return 1.;
  }

  double weight_TTbar_nJetCentral_Down(int nGenTop, int nJetsCentral) {
    if (nGenTop == 2) {
      double sf_njets = 1.2213 * (0.9437417 - 0.02498308*nJetsCentral);
      return sf_njets;
    }
    return 1.;
  }

  double weight_EWK_VH(int is_ZH, double GenVbosons_pt) {
    // Apply the EWK signal correction to ZH signal sample.
    if (is_ZH) {
      if (GenVbosons_pt>0 && GenVbosons_pt<20) {
        return 0.963285466565;
      } else if (GenVbosons_pt>20 && GenVbosons_pt<40) {
        return 0.960945354673;
      } else if (GenVbosons_pt>40 && GenVbosons_pt<60) {
        return 0.958125845181;
      } else if (GenVbosons_pt>60 && GenVbosons_pt<80) {
        return 0.956040178567;
      } else if (GenVbosons_pt>80 && GenVbosons_pt<100) {
        return 0.954829397636;
      } else if (GenVbosons_pt>100 && GenVbosons_pt<120) {
        return 0.955710296883;
      } else if (GenVbosons_pt>120 && GenVbosons_pt<140) {
        return 0.958933532209;
      } else if (GenVbosons_pt>140 && GenVbosons_pt<160) {
        return 0.960145045414;
      } else if (GenVbosons_pt>160 && GenVbosons_pt<180) {
        return 0.958798002291;
      } else if (GenVbosons_pt>180 && GenVbosons_pt<200) {
        return 0.954417774192;
      } else if (GenVbosons_pt>200 && GenVbosons_pt<220) {
        return 0.950510442261;
      } else if (GenVbosons_pt>220 && GenVbosons_pt<240) {
        return 0.943227377306;
      } else if (GenVbosons_pt>240 && GenVbosons_pt<260) {
        return 0.93876408527;
      } else if (GenVbosons_pt>260 && GenVbosons_pt<280) {
        return 0.929695202634;
      } else if (GenVbosons_pt>280 && GenVbosons_pt<300) {
        return 0.924760963363;
      } else if (GenVbosons_pt>300 && GenVbosons_pt<320) {
        return 0.915490487966;
      } else if (GenVbosons_pt>320 && GenVbosons_pt<340) {
        return 0.906099906197;
      } else if (GenVbosons_pt>340 && GenVbosons_pt<360) {
        return 0.901215165553;
      } else if (GenVbosons_pt>360 && GenVbosons_pt<380) {
        return 0.891758657639;
      } else if (GenVbosons_pt>380 && GenVbosons_pt<400) {
        return 0.882966044972;
      } else if (GenVbosons_pt>400 && GenVbosons_pt<420) {
        return 0.87385231926;
      } else if (GenVbosons_pt>420 && GenVbosons_pt<440) {
        return 0.868965659921;
      } else if (GenVbosons_pt>440 && GenVbosons_pt<460) {
        return 0.866066608534;
      } else if (GenVbosons_pt>460 && GenVbosons_pt<480) {
        return 0.854245464035;
      } else if (GenVbosons_pt>480 && GenVbosons_pt<500) {
        return 0.846544787867;
      } else if (GenVbosons_pt>500) {
        return 0.846544787867;
      }
    }
    return 1.;
  }

  double ptWeightDY(double lheV_pt, double sign = 1.) {
    double SF = 1.;
    if (50. < lheV_pt && lheV_pt < 100.) {
      return SF = 0.873885 + 0.00175853*sign*lheV_pt;
    } else if (lheV_pt > 100) {
      return SF = 1.10651 - 0.000705265*sign*lheV_pt;
    } else {
      return SF;
    }
  }

  // weights correction for EWK NLO correction
  double ptWeightZllH(double genHPt) {
    double SF = 1.;
    if (genHPt > 1.) SF = 0.999757 - 0.000174527*genHPt;
    return (SF > 0) ? SF : 0;
  }

  // weights correction for EWK NLO correction from Atlas
  double ewkAtlas8TeVZllH(double genHPt, double genZPt) {
    if (genHPt < 1.) return 1.;

    double hll8_contents[95] = {0.000664024, -0.00357095, -0.00767076, -0.00967366, -0.0134844, -0.0157148, -0.0181885, -0.0209647, -0.0232788, -0.0252373, -0.0265634, -0.0275069, -0.0285776, -0.0281683, -0.0294206, -0.0299975, -0.0308047, -0.0311716, -0.030913, -0.0324821, -0.0323192, -0.0324639, -0.0319356, -0.0322621, -0.0331146, -0.0338905, -0.0345189, -0.0358591, -0.0358407, -0.040018, -0.0396389, -0.0407177, -0.0445103, -0.0441406, -0.0471215, -0.0463301, -0.0513777, -0.0536773, -0.0546446, -0.0568508, -0.0590333, -0.0612157, -0.0633981, -0.0655805, -0.067763, -0.0699454, -0.0721278, -0.0743103, -0.0764927, -0.0786751, -0.0808575, -0.08304, -0.0852224, -0.0874048, -0.0895872, -0.0917697, -0.0939521, -0.0961345, -0.098317, -0.100499, -0.102682, -0.104864, -0.107047, -0.109229, -0.111412, -0.113594, -0.115776, -0.117959, -0.120141, -0.122324, -0.124506, -0.126689, -0.128871, -0.131053, -0.133236, -0.135418, -0.137601, -0.139783, -0.141965, -0.144148, -0.14633, -0.148513, -0.150695, -0.152878, -0.15506, -0.157242, -0.159425, -0.161607, -0.16379, -0.165972, -0.168155, -0.170337, -0.172519, -0.174702, -0.176884};

    int corrBin = int((genZPt-25.) / 5.);
    if (corrBin < 0) corrBin = 0;
    if (corrBin > 94) corrBin = 94;

    double SF = 1. + hll8_contents[corrBin];
    return SF;
  }

  double minCSVold(int EVENT_run, double hJet_csv0, double hJet_csv1, double hJet_csvOld0, double hJet_csvOld1) {
    if (EVENT_run < 2) {
      return std::min(hJet_csvOld0, hJet_csvOld1);
    } else {
      return std::min(hJet_csv0, hJet_csv1);
    }
  }

  double maxCSVold(int EVENT_run, double hJet_csv0, double hJet_csv1, double hJet_csvOld0, double hJet_csvOld1) {
    if (EVENT_run < 2) {
      return std::max(hJet_csvOld0, hJet_csvOld1);
    } else {
      return std::max(hJet_csv0, hJet_csv1);
    }
  }

  double eleEff(double pt, double eta) {
    if (pt >= 20.0 && pt < 25.0 && eta < -1.57000005245 && eta > -2.5) return 0.955429255962;
    if (pt >= 20.0 && pt < 25.0 && eta < -1.44000005722 && eta > -1.57000005245) return 0.927380621433;
    if (pt >= 20.0 && pt < 25.0 && eta < -0.800000011921 && eta > -1.44000005722) return 0.752811849117;
    if (pt >= 20.0 && pt < 25.0 && eta < 0.0 && eta > -0.800000011921) return 0.758898377419;
    if (pt >= 20.0 && pt < 25.0 && eta < 0.800000011921 && eta > 0.0) return 0.755395233631;
    if (pt >= 20.0 && pt < 25.0 && eta < 1.44000005722 && eta > 0.800000011921) return 0.765102207661;
    if (pt >= 20.0 && pt < 25.0 && eta < 1.57000005245 && eta > 1.44000005722) return 0.889689266682;
    if (pt >= 20.0 && pt < 25.0 && eta < 2.5 && eta > 1.57000005245) return 0.960859179497;
    if (pt >= 25.0 && pt < 30.0 && eta < -1.57000005245 && eta > -2.5) return 0.973749041557;
    if (pt >= 25.0 && pt < 30.0 && eta < -1.44000005722 && eta > -1.57000005245) return 0.95562505722;
    if (pt >= 25.0 && pt < 30.0 && eta < -0.800000011921 && eta > -1.44000005722) return 0.854979991913;
    if (pt >= 25.0 && pt < 30.0 && eta < 0.0 && eta > -0.800000011921) return 0.859907507896;
    if (pt >= 25.0 && pt < 30.0 && eta < 0.800000011921 && eta > 0.0) return 0.856850981712;
    if (pt >= 25.0 && pt < 30.0 && eta < 1.44000005722 && eta > 0.800000011921) return 0.854210674763;
    if (pt >= 25.0 && pt < 30.0 && eta < 1.57000005245 && eta > 1.44000005722) return 0.961790204048;
    if (pt >= 25.0 && pt < 30.0 && eta < 2.5 && eta > 1.57000005245) return 0.984032571316;
    if (pt >= 30.0 && pt < 35.0 && eta < -1.57000005245 && eta > -2.5) return 0.975343048573;
    if (pt >= 30.0 && pt < 35.0 && eta < -1.44000005722 && eta > -1.57000005245) return 0.971944391727;
    if (pt >= 30.0 && pt < 35.0 && eta < -0.800000011921 && eta > -1.44000005722) return 0.916224777699;
    if (pt >= 30.0 && pt < 35.0 && eta < 0.0 && eta > -0.800000011921) return 0.918354153633;
    if (pt >= 30.0 && pt < 35.0 && eta < 0.800000011921 && eta > 0.0) return 0.914323210716;
    if (pt >= 30.0 && pt < 35.0 && eta < 1.44000005722 && eta > 0.800000011921) return 0.918937265873;
    if (pt >= 30.0 && pt < 35.0 && eta < 1.57000005245 && eta > 1.44000005722) return 0.949005782604;
    if (pt >= 30.0 && pt < 35.0 && eta < 2.5 && eta > 1.57000005245) return 0.978435099125;
    if (pt >= 35.0 && pt < 40.0 && eta < -1.57000005245 && eta > -2.5) return 0.977416038513;
    if (pt >= 35.0 && pt < 40.0 && eta < -1.44000005722 && eta > -1.57000005245) return 0.96676659584;
    if (pt >= 35.0 && pt < 40.0 && eta < -0.800000011921 && eta > -1.44000005722) return 0.949674129486;
    if (pt >= 35.0 && pt < 40.0 && eta < 0.0 && eta > -0.800000011921) return 0.955099225044;
    if (pt >= 35.0 && pt < 40.0 && eta < 0.800000011921 && eta > 0.0) return 0.95435655117;
    if (pt >= 35.0 && pt < 40.0 && eta < 1.44000005722 && eta > 0.800000011921) return 0.949129760265;
    if (pt >= 35.0 && pt < 40.0 && eta < 1.57000005245 && eta > 1.44000005722) return 0.966977357864;
    if (pt >= 35.0 && pt < 40.0 && eta < 2.5 && eta > 1.57000005245) return 0.974912643433;
    if (pt >= 40.0 && pt < 45.0 && eta < -1.57000005245 && eta > -2.5) return 0.976379692554;
    if (pt >= 40.0 && pt < 45.0 && eta < -1.44000005722 && eta > -1.57000005245) return 0.96400809288;
    if (pt >= 40.0 && pt < 45.0 && eta < -0.800000011921 && eta > -1.44000005722) return 0.967367231846;
    if (pt >= 40.0 && pt < 45.0 && eta < 0.0 && eta > -0.800000011921) return 0.976113498211;
    if (pt >= 40.0 && pt < 45.0 && eta < 0.800000011921 && eta > 0.0) return 0.970276534557;
    if (pt >= 40.0 && pt < 45.0 && eta < 1.44000005722 && eta > 0.800000011921) return 0.965076506138;
    if (pt >= 40.0 && pt < 45.0 && eta < 1.57000005245 && eta > 1.44000005722) return 0.962332487106;
    if (pt >= 40.0 && pt < 45.0 && eta < 2.5 && eta > 1.57000005245) return 0.975438177586;
    if (pt >= 45.0 && pt < 50.0 && eta < -1.57000005245 && eta > -2.5) return 0.974642693996;
    if (pt >= 45.0 && pt < 50.0 && eta < -1.44000005722 && eta > -1.57000005245) return 0.976564764977;
    if (pt >= 45.0 && pt < 50.0 && eta < -0.800000011921 && eta > -1.44000005722) return 0.975869596004;
    if (pt >= 45.0 && pt < 50.0 && eta < 0.0 && eta > -0.800000011921) return 0.984652400017;
    if (pt >= 45.0 && pt < 50.0 && eta < 0.800000011921 && eta > 0.0) return 0.976347208023;
    if (pt >= 45.0 && pt < 50.0 && eta < 1.44000005722 && eta > 0.800000011921) return 0.974825322628;
    if (pt >= 45.0 && pt < 50.0 && eta < 1.57000005245 && eta > 1.44000005722) return 0.968582391739;
    if (pt >= 45.0 && pt < 50.0 && eta < 2.5 && eta > 1.57000005245) return 0.973198652267;
    if (pt >= 50.0 && pt < 200.0 && eta < -1.57000005245 && eta > -2.5) return 0.976521909237;
    if (pt >= 50.0 && pt < 200.0 && eta < -1.44000005722 && eta > -1.57000005245) return 0.968289792538;
    if (pt >= 50.0 && pt < 200.0 && eta < -0.800000011921 && eta > -1.44000005722) return 0.980251729488;
    if (pt >= 50.0 && pt < 200.0 && eta < 0.0 && eta > -0.800000011921) return 0.988307952881;
    if (pt >= 50.0 && pt < 200.0 && eta < 0.800000011921 && eta > 0.0) return 0.983105957508;
    if (pt >= 50.0 && pt < 200.0 && eta < 1.44000005722 && eta > 0.800000011921) return 0.976873219013;
    if (pt >= 50.0 && pt < 200.0 && eta < 1.57000005245 && eta > 1.44000005722) return 0.970499277115;
    if (pt >= 50.0 && pt < 200.0 && eta < 2.5 && eta > 1.57000005245) return 0.972752988338;
    return 1.;
  }

  double mueEff(int Vtype, double eta0, double eta1, double pt0, double pt1) {
    if (Vtype == 0) {
      return 1.;
    } else if (Vtype == 1) {
      return eleEff(pt0, eta0) * eleEff(pt1, eta1);
    } else {
      return 1.;
    }
  }

  float weight2(int i) {
    if (i < 0 || i > 51) return 1;

    double mc2[52] = {
      4.8551E-07,
      1.74806E-06,
      3.30868E-06,
      1.62972E-05,
      4.95667E-05,
      0.000606966,
      0.003307249,
      0.010340741,
      0.022852296,
      0.041948781,
      0.058609363,
      0.067475755,
      0.072817826,
      0.075931405,
      0.076782504,
      0.076202319,
      0.074502547,
      0.072355135,
      0.069642102,
      0.064920999,
      0.05725576,
      0.047289348,
      0.036528446,
      0.026376131,
      0.017806872,
      0.011249422,
      0.006643385,
      0.003662904,
      0.001899681,
      0.00095614,
      0.00050028,
      0.000297353,
      0.000208717,
      0.000165856,
      0.000139974,
      0.000120481,
      0.000103826,
      8.88868E-05,
      7.53323E-05,
      6.30863E-05,
      5.21356E-05,
      4.24754E-05,
      3.40876E-05,
      2.69282E-05,
      2.09267E-05,
      1.5989E-05,
      4.8551E-06,
      2.42755E-06,
      4.8551E-07,
      2.42755E-07,
      1.21378E-07,
      4.8551E-08
    };

    double data2[52] = {
      4.7027e-05,
      0.000281565,
      0.00028437,
      0.00038727,
      0.000569421,
      0.000952123,
      0.00319069,
      0.0203182,
      0.0699736,
      0.130068,
      0.180077,
      0.198876,
      0.174006,
      0.118772,
      0.06317,
      0.026531,
      0.00902068,
      0.00258006,
      0.000659883,
      0.000164919,
      4.46606e-05,
      1.44451e-05,
      5.83791e-06,
      2.78026e-06,
      1.40517e-06,
      7.0225e-07,
      3.36679e-07,
      1.53294e-07,
      6.60997e-08,
      2.69735e-08,
      1.04154e-08,
      3.80539e-09,
      1.31553e-09,
      4.30311e-10,
      1.3318e-10,
      3.90006e-11,
      1.08063e-11,
      2.83309e-12,
      7.02782e-13,
      1.64952e-13,
      3.66335e-14,
      7.69806e-15,
      1.53064e-15,
      2.87972e-16,
      5.12673e-17,
      8.63513e-18,
      1.37688e-18,
      2.04286e-19,
      3.72485e-20,
      0,
      0,
      0
    };

    return data2[i] / mc2[i];
  }

  float weight2_up(int i) {
    if (i < 0 || i > 51) return 1;

    double mc2[52] = {
      4.8551E-07,
      1.74806E-06,
      3.30868E-06,
      1.62972E-05,
      4.95667E-05,
      0.000606966,
      0.003307249,
      0.010340741,
      0.022852296,
      0.041948781,
      0.058609363,
      0.067475755,
      0.072817826,
      0.075931405,
      0.076782504,
      0.076202319,
      0.074502547,
      0.072355135,
      0.069642102,
      0.064920999,
      0.05725576,
      0.047289348,
      0.036528446,
      0.026376131,
      0.017806872,
      0.011249422,
      0.006643385,
      0.003662904,
      0.001899681,
      0.00095614,
      0.00050028,
      0.000297353,
      0.000208717,
      0.000165856,
      0.000139974,
      0.000120481,
      0.000103826,
      8.88868E-05,
      7.53323E-05,
      6.30863E-05,
      5.21356E-05,
      4.24754E-05,
      3.40876E-05,
      2.69282E-05,
      2.09267E-05,
      1.5989E-05,
      4.8551E-06,
      2.42755E-06,
      4.8551E-07,
      2.42755E-07,
      1.21378E-07,
      4.8551E-08
    };

    double data2_P[52] = {
      4.02952e-05,
      0.000254497,
      0.000280989,
      0.000335717,
      0.00050839,
      0.000746426,
      0.00186153,
      0.0100804,
      0.0445114,
      0.0982319,
      0.150143,
      0.184585,
      0.184905,
      0.148911,
      0.0954775,
      0.0489179,
      0.020324,
      0.00701388,
      0.00208482,
      0.000564834,
      0.000150948,
      4.35228e-05,
      1.47302e-05,
      6.1078e-06,
      2.96971e-06,
      1.54517e-06,
      8.03286e-07,
      4.03918e-07,
      1.94127e-07,
      8.88486e-08,
      3.86863e-08,
      1.60215e-08,
      6.31045e-09,
      2.36387e-09,
      8.42156e-10,
      2.8534e-10,
      9.19462e-11,
      2.81777e-11,
      8.21254e-12,
      2.27642e-12,
      6.00106e-13,
      1.50456e-13,
      3.58753e-14,
      8.13565e-15,
      1.75469e-15,
      3.59936e-16,
      7.02215e-17,
      1.30278e-17,
      2.29864e-18,
      3.87514e-19,
      7.22359e-20,
      7.82805e-22
    };

    return data2_P[i] / mc2[i];
  }

  float weight2_down(int i) {
    if (i < 0 || i > 51) return 1;

    double mc2[52] = {
      4.8551E-07,
      1.74806E-06,
      3.30868E-06,
      1.62972E-05,
      4.95667E-05,
      0.000606966,
      0.003307249,
      0.010340741,
      0.022852296,
      0.041948781,
      0.058609363,
      0.067475755,
      0.072817826,
      0.075931405,
      0.076782504,
      0.076202319,
      0.074502547,
      0.072355135,
      0.069642102,
      0.064920999,
      0.05725576,
      0.047289348,
      0.036528446,
      0.026376131,
      0.017806872,
      0.011249422,
      0.006643385,
      0.003662904,
      0.001899681,
      0.00095614,
      0.00050028,
      0.000297353,
      0.000208717,
      0.000165856,
      0.000139974,
      0.000120481,
      0.000103826,
      8.88868E-05,
      7.53323E-05,
      6.30863E-05,
      5.21356E-05,
      4.24754E-05,
      3.40876E-05,
      2.69282E-05,
      2.09267E-05,
      1.5989E-05,
      4.8551E-06,
      2.42755E-06,
      4.8551E-07,
      2.42755E-07,
      1.21378E-07,
      4.8551E-08
    };

    double data2_M[52] = {
      5.48234e-05,
      0.00031101,
      0.000294382,
      0.000447642,
      0.000648535,
      0.00130872,
      0.00642627,
      0.0386726,
      0.102155,
      0.165772,
      0.205682,
      0.198571,
      0.146521,
      0.0819434,
      0.0350992,
      0.0118095,
      0.00325042,
      0.000781694,
      0.000181477,
      4.58499e-05,
      1.41184e-05,
      5.55293e-06,
      2.58233e-06,
      1.26088e-06,
      6.01759e-07,
      2.73079e-07,
      1.16867e-07,
      4.70691e-08,
      1.78332e-08,
      6.35519e-09,
      2.13024e-09,
      6.71619e-10,
      1.99164e-10,
      5.55507e-11,
      1.45733e-11,
      3.59601e-12,
      8.34591e-13,
      1.82189e-13,
      3.74081e-14,
      7.22453e-15,
      1.31237e-15,
      2.24241e-16,
      3.60383e-17,
      5.44906e-18,
      7.76605e-19,
      1.04021e-19,
      1.24685e-20,
      0,
      0,
      0,
      0,
      0
    };

    return data2_M[i] / mc2[i];
  }

  double ptWeightQCD(int nGenVbosons, double lheHT, int GenVbosons_pdgId) {
    double SF = 1.;

    if (lheHT > 100 && nGenVbosons == 1) {
      if (GenVbosons_pdgId == 23) {
        // Z Boson
        SF = (lheHT > 100 && lheHT < 200) * 1.588 * (280.35/409.860000) +
             (lheHT > 200 && lheHT < 400) * 1.438 * (77.67/110.880000) +
             (lheHT > 400 && lheHT < 600) * 1.494 * (10.73/13.189) +
             (lheHT > 600) * 1.139 * (4.116/4.524300);
      } else if (abs(GenVbosons_pdgId) == 24) {
        // W Bosons
        SF = (lheHT > 100 && lheHT < 200) * 1.588 * (1345 / (1.23*1.29e3)) +
             (lheHT > 200 && lheHT < 400) * 1.438 * (359.7 / (1.23*3.86e2)) +
             (lheHT > 400 && lheHT < 600) * 1.494 * (48.91 / (1.23*47.9)) +
             (lheHT > 600) * 1.139 * (18.77 / (1.23*19.9));
      }
    }

    return (SF > 0) ? SF : 0;
  }

// weights correction for EWK NLO correction (for ZllHbb only !!!)
double ptWeightEWK_Zllv2(int isDY ,double GenVbosons_pt){
    double SF = 1.;
    //for Z options
    if (isDY > 0){
        if (GenVbosons_pt > 100. && GenVbosons_pt < 3000) SF = -0.1808051+6.04146*(TMath::Power((GenVbosons_pt+759.098),-0.242556));
    }
    return SF>0?SF:0;
}

// weights correction for EWK NLO correction (for ZllHbb only !!!)
double ptWeightEWK_Zll(int nGenVbosons,double GenVbosons_pt,int VtypeSim, int nGenTop, int nGenHiggsBoson){
    double SF = 1.;
    if (nGenVbosons ==1 && nGenTop == 0 && nGenHiggsBoson == 0)
    {
        if (VtypeSim == 0 || VtypeSim == 1 || VtypeSim == 4 || VtypeSim == 5)
        {
            {
                //for Z options
                if (GenVbosons_pt > 100. && GenVbosons_pt < 3000) SF = -0.1808051+6.04146*(TMath::Power((GenVbosons_pt+759.098),-0.242556));
            }
        }
    }
    return SF>0?SF:0;
}
  // weights correction for EWK NLO correction
  double ptWeightEWK(int nGenVbosons, double GenVbosons_pt, int VtypeSim, int GenVbosons_pdgId) {
    double SF = 1.;

    if (nGenVbosons == 1) {
      if (VtypeSim == 0 || VtypeSim == 1 || VtypeSim == 4 || VtypeSim == 5) {
        if (GenVbosons_pdgId == 23) {
          // Z Boson
          if (GenVbosons_pt > 100. && GenVbosons_pt < 3000) SF = -0.1808051 + 6.04146*TMath::Power(GenVbosons_pt + 759.098, -0.242556);
        }
      } else if (VtypeSim == 2 || VtypeSim == 3) {
        // W Bosons
        if (GenVbosons_pdgId == 24 || GenVbosons_pdgId == -24) {
          if (GenVbosons_pt > 100. && GenVbosons_pt < 3000) SF = -0.830041 + 7.93714*TMath::Power(GenVbosons_pt + 877.978,-0.213831);
        }
      }
    }

    return (SF > 0) ? SF : 0;
  }

double LOtoNLOWeightEtabb(double etabb){

    double SF = 1.;
    if(etabb < 5){

       //From 8/11/16
       SF =   0.940679 + 0.0306119*etabb -0.0134403*etabb*etabb + 0.0132179*etabb*etabb*etabb -0.00143832*etabb*etabb*etabb*etabb;

    }

    return SF;

}

double LOtoNLOWeightBjetSplitEtabb(double etabb, int njets){

    double SF = 1.;
    if(etabb < 5){
        if(njets < 1){

        //0b jets
        SF =   0.935422 + 0.0403162*etabb -0.0089026*etabb*etabb +0.0064324*etabb*etabb*etabb -0.000212443*etabb*etabb*etabb*etabb;
        }else if(njets == 1){
        //1b jets
        SF =   0.962415 +0.0329463*etabb -0.0414479*etabb*etabb +0.0240993*etabb*etabb*etabb -0.00278271*etabb*etabb*etabb*etabb;
        }else if(njets >=2){
        //2b jets
        SF =   (0.721265 -0.105643*etabb -0.0206835*etabb*etabb +0.00558626*etabb*etabb*etabb)*TMath::Exp(0.450244*etabb);
        }
    }

    return SF;

}

double LOtoNLOWeightBjetSplitEtabb2017(double etabb, int njets, double Vpt){

    double SF = 1.;
    if(etabb < 5){
        if (Vpt < 150){
                //low pt
            if(njets < 1){
                //0b jets
                SF =   1.02053 - 0.0600883*etabb + 0.0217384*etabb*etabb;
                }else if(njets == 1){
                //1b jets
                SF =   1.05412 - 0.122462*etabb + 0.0395305*etabb*etabb;
                }else if(njets >=2){
                //2b jets
                SF =   1.01633 - 0.0862970*etabb + 0.0362290*etabb*etabb;
                }

                //high pt
            }else if (Vpt > 150){
                //0b jets
                SF =   1.00853 - 0.0268244*etabb + 0.0101112*etabb*etabb;
                }else if(njets == 1){
                //1b jets
                SF =   (34.6575 + 2.00480*etabb + 1.54075*etabb*etabb)*0.0315321*TMath::Exp(-0.198796*etabb);
                }else if(njets >=2){
                //2b jets
                SF =   1.05759 - 0.0946755*etabb + 0.0239830*etabb*etabb;
            }
        }
    return SF;
}

double LOtoNLOWeightBjetSplitEtabb2017v2(double etabb, int njets){
    double SF = 1.;
    if (etabb > 5.0) etabb = 5.0;
    if(njets < 1){
        //0b jets
        SF = (0.958 + 0.0286 * etabb + 0.0014156 * etabb * etabb);
    } else if(njets == 1){
        //1b jets
        SF = ((0.972 - 0.264 * etabb + 0.026919 * etabb * etabb) * TMath::Exp(0.29901 * etabb));
    } else if(njets >=2){
        //2b jets
        SF = (0.81 + 0.1493 * etabb - 0.000965976 * etabb * etabb); 
    }
    return SF;
}

// residual deltaEta correction after Vpt reweighting is applied
// etabb = abs(GenJet_eta[0]-GenJet_eta[1])
//2b: 8.967e-01 + 7.218e-02 * etabb + 4.231e-02 * etabb*etabb - 1.092e-02 * etabb*etabb*etabb
//1b: 9.596e-01 + 1.806e-02 * etabb + 4.633e-03 * etabb*etabb - 8.210e-05 * etabb*etabb*etabb
//0b: 9.096e-01 + 8.609e-02 * etabb - 1.938e-02 * etabb*etabb + 2.080e-03 * etabb*etabb*etabb

// VHbb::LOtoNLOWeightBjetSplitEtabb2017ZJetsResidual(abs(GenJet_eta[0]-GenJet_eta[1]), Sum$(GenJet_pt>25 && abs(GenJet_eta)<2.4 && GenJet_hadronFlavour==5))

double LOtoNLOWeightBjetSplitEtabb2017ZJetsResidual(double etabb, int njets){
    double SF = 1.;
    if (etabb > 5.0) etabb = 5.0;
    if(njets < 1){
        //0b jets
        SF = 9.096e-01 + 8.609e-02 * etabb - 1.938e-02 * etabb*etabb + 2.080e-03 * etabb*etabb*etabb;
    } else if(njets == 1){
        //1b jets
        SF = 9.596e-01 + 1.806e-02 * etabb + 4.633e-03 * etabb*etabb - 8.210e-05 * etabb*etabb*etabb; 
    } else if(njets >=2){
        //2b jets
        SF = 8.967e-01 + 7.218e-02 * etabb + 4.231e-02 * etabb*etabb - 1.092e-02 * etabb*etabb*etabb; 
    }
    return SF;
}

double LOtoNLOWeightBjetSplitVpt2017(double lhevpt, int nbjets) {
    if (lhevpt<100.0) lhevpt=100.0;
    if (lhevpt>500.0) lhevpt=500.0;
    double SF = 1.0;
    if (nbjets < 1) {
        SF = 1.056-9.480e-04*lhevpt;
    } else if (nbjets == 1) {
        SF = 1.05-1.235e-03*lhevpt;
    } else if (nbjets > 1) {
        SF = 1.15-2.068e-03*lhevpt;
    }
    return SF;
}


int checkgen(int Jet_mcIdx0, int Jet_mcIdx1, int length){

   if (Jet_mcIdx0 == -1 || Jet_mcIdx1 == -1 || Jet_mcIdx0 >= length || Jet_mcIdx1 >= length){
       return 0;
   }else{
       return Jet_mcIdx0;
   }

}

double LOtoNLOWeightgenEtabb(double etabb){

    double SF = 1.;
        if(etabb < 5 && etabb > 0.000001){
            //GeneralCuts wQCD
            //SF = 0.889908 + 0.0185636*etabb +0.0353842*etabb*etabb + 0.00142895*etabb*etabb*etabb -0.000413957*etabb*etabb*etabb*etabb;
            //GeneralCuts noQCD
            //SF = 0.889875 + 0.0185131*etabb +0.0355196*etabb*etabb + 0.00138352*etabb*etabb*etabb -0.00041181*etabb*etabb*etabb*etabb;
            //Z+light
            //w QCD
            SF =  0.870811+ 0.0386754*etabb + 0.0479672*etabb*etabb - 0.011603*etabb*etabb*etabb + 0.00252448*etabb*etabb*etabb*etabb;
        }
    return SF;

}

double LOtoNLOWeightgenEtabbPtJ(double etabb){

    double SF = 1.;
        if(etabb < 5 && etabb > 0.000001){
            //SF = 0.860459 + 0.0191792*etabb + 0.0546046*etabb*etabb - 0.00623299*etabb*etabb*etabb + 0.000123637*etabb*etabb*etabb*etabb;
            //Z+light
            //w QCD
            SF = 0.841774 + 0.0725955*etabb + 0.0202381*etabb*etabb + 0.0016097*etabb*etabb*etabb + 0.0000391441*etabb*etabb*etabb*etabb;
        }

    return SF;

}

//int FirstGenJet(Float_t GenJet_eta[], Float_t GenJet_phi[], Float_t Jet_eta, Float_t Jet_phi, int gensize = 15){
//
//    double dR = 999;
//    int ind = -1;
//
//    for (int k = 0; k < gensize; ++k){
//        if (deltaR(GenJet_eta[k], GenJet_phi[k], Jet_eta, Jet_phi) < dR){
//           dR = deltaR(GenJet_eta[k], GenJet_phi[k], Jet_eta[l], Jet_phi[l]);
//           ind = k;
//        }
//    }
//    if dR <
//
//}

float puWeight_ichep(float I){

//int i = std::nearbyint(I);
int i = (int) I;

double puw[38]={0.00026954692859,
                0.00834201570744,
                0.0146551644559,
                0.0267593718187,
                0.045546866213,
                0.0351739785649,
                0.0389242914756,
                0.131305445658,
                0.337172065156,
                0.645651266938,
                0.908493559723,
                1.09739459449,
                1.26533979742,
                1.41586705341,
                1.53524109717,
                1.49032150282,
                1.39123249386,
                1.45737145535,
                1.38588035333,
                1.4464676134,
                1.23282560995,
                1.08387439504,
                1.02663484504,
                0.97464906611,
                0.829692105572,
                0.758597838706,
                0.554540190093,
                0.442016886186,
                0.291492210196,
                0.203297812168,
                0.131804681706,
                0.0834846669777,
                0.0680713812995,
                0.0497105409204,
                0.0496692836227,
                0.0581182475108,
                0.089209326104,
                0.0941178579122,
               };

if (i < 0) return 1.;
if (i > 37) return puw[37];

return puw[i];

}

float puWeight_ichep_up(float I){

int i = (int) I;
//int i = std::nearbyint(I);

double puw[38]={0.000168728884,
                0.006518139648,
                0.012767671099,
                0.021849922433,
                0.041129948359,
                0.031318966078,
                0.031506395408,
                0.069139953987,
                0.201682923142,
                0.431573069513,
                0.676688315382,
                0.885061147726,
                1.02772732515,
                1.1496216224,
                1.26358384716,
                1.24734104964,
                1.19950685946,
                1.30515668683,
                1.28849179873,
                1.39604961327,
                1.24030942651,
                1.13881318563,
                1.12938083469,
                1.13608683795,
                1.04512825001,
                1.04965380671,
                0.848278572582,
                0.748408550686,
                0.548182746816,
                0.426699751655,
                0.308798449343,
                0.217654143858,
                0.197782762841,
                0.16222993513 ,
                0.183774325453,
                0.245159931575,
                0.426559360962,
                0.491832630157,
               };

if (i < 0) return 1.;
if (i > 37) return puw[37];

return puw[i];

}

float puWeight_ichep_down(float I){

int i = (int) I;
//int i = std::nearbyint(I);

double puw[38]={0.000394948025924,
                0.010589291412,
                0.0168294422346,
                0.0324871095383,
                0.0512240404177,
                0.0392771251619,
                0.0585204632322,
                0.252037472998,
                0.543974217454,
                0.93479398165,
                1.17260170322,
                1.36170067318,
                1.5793616475,
                1.74509976454,
                1.86131835195,
                1.75449370357,
                1.57204882742,
                1.57954789714,
                1.44148253688,
                1.43725104611,
                1.16760012792,
                0.975000220575,
                0.8640144709,
                0.750072693636,
                0.574219048746,
                0.469714734739,
                0.306667353503,
                0.217186088114,
                0.126764763282,
                0.0784330110012,
                0.0451854756945,
                0.0252751549695,
                0.0180049218223,
                0.0113901195301,
                0.00987729209672,
                0.0103547199185,
                0.0158597867448,
                0.0210296659761,
               };

if (i < 0) return 1.;
if (i > 37) return puw[37];

return puw[i];

}

int argmax8(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8) {
        std::vector<double> p = {p1,p2,p3,p4,p5,p6,p7,p8};
            return std::distance(p.begin(), std::max_element(p.begin(), p.end()));
}
double max8(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8) {
        std::vector<double> p = {p1,p2,p3,p4,p5,p6,p7,p8};
        return *std::max_element(p.begin(), p.end());
}

float MindRDoubleSingleBtag(double FatJet_eta, double FatJet_phi, double Jet1_eta, double Jet1_phi, int Jet1_ind,double Jet2_eta, double Jet2_phi, int Jet2_ind) {
    double dR1 = 999;
    double dR2 = 999;
    double dRMin = 999;

    if (Jet1_ind != -1){
        dR1 = deltaR(FatJet_eta,FatJet_phi,Jet1_eta,Jet1_phi);
    };
    if (Jet2_ind != -1){
        dR2 = deltaR(FatJet_eta,FatJet_phi,Jet2_eta,Jet2_phi);
    };

   dRMin = min(dR1,dR2); 
   if (dRMin == 999){
       return 0;} 
   else{return dRMin;}
}

float MaxdRDoubleSingleBtag(double FatJet_eta, double FatJet_phi, double Jet1_eta, double Jet1_phi, int Jet1_ind,double Jet2_eta, double Jet2_phi, int Jet2_ind) {
    double dR1 = 999;
    double dR2 = 999;
    double dRMax = 999;

    if (Jet1_ind != -1){
        dR1 = deltaR(FatJet_eta,FatJet_phi,Jet1_eta,Jet1_phi);
    };
    if (Jet2_ind != -1){
        dR2 = deltaR(FatJet_eta,FatJet_phi,Jet2_eta,Jet2_phi);
    };

   dRMax = max(dR1,dR2); 
   if (dRMax == 999){
       return 0;} 
   else{return dRMax;}

}

double m_inv2(double pt1, double eta1, double phi1, double m1, double pt2, double eta2, double phi2, double m2) {
    TLorentzVector v1;
    TLorentzVector v2;
    v1.SetPtEtaPhiM(pt1,eta1,phi1,m1);
    v2.SetPtEtaPhiM(pt2,eta2,phi2,m2);
    return (v1+v2).M();
}

double m_inv3(double pt1, double eta1, double phi1, double m1, double pt2, double eta2, double phi2, double m2, double pt3, double eta3, double phi3, double m3) {
    TLorentzVector v1;
    TLorentzVector v2;
    TLorentzVector v3;
    v1.SetPtEtaPhiM(pt1,eta1,phi1,m1);
    v2.SetPtEtaPhiM(pt2,eta2,phi2,m2);
    v3.SetPtEtaPhiM(pt3,eta3,phi3,m3);
    TLorentzVector r;
    r = v1 + v2 + v3;
    return r.M();
}

 double pt (double pt1, double phi1, double pt2, double phi2) {

   double x = pt1 * TMath::Cos(phi1) + pt2 * TMath::Cos(phi2);
   double y = pt1 * TMath::Sin(phi1) + pt2 * TMath::Sin(phi2);

   return TMath::Sqrt(x*x + y*y);

 }

}

float triggerMet(float x) {
  return ((0.5+(0.5*erf((x-77.7235)/47.6344)))*0.994699)-0.00492778;
}

float triggerMetUp(float x) {
  return ((0.5+(0.5*erf((x-71.6922)/47.6344)))*0.998982)-0.00492778;
}

float triggerMetDown(float x) {
  return ((0.5+(0.5*erf((x-83.7549)/47.6344)))*0.990415)-0.00492778;
}



