#!/usr/bin/env python
import ROOT
import numpy as np
import array
import os
from BranchTools import Collection
from BranchTools import AddCollectionsModule

class LOtoNLOweight(AddCollectionsModule):

    def __init__(self, branchName='weightLOtoNLO', year=2016, addV3=False, addV4=False):
        super(LOtoNLOweight, self).__init__()
        self.branchName = branchName
        self.year = int(year)
        self.version = 5

        self.addV3 = addV3
        self.addV4 = addV4

        self.ZjetsToNuNuNLOweightNjet2017 = [
                [0,0,1.228032,0.006932,-0.001342,0.000038],
                [0,1,0.918723,0.009699,-0.000621,0.000052],
                [0,2,0.610806,0.014035,-0.000130,0.000075],
                [0,3,0.450058,0.025563,0.000270,0.000135],
                [0,4,0.508454,0.047272,-0.000388,0.000246],
                [0,5,0.451616,0.094075,-0.000509,0.000487],
                [1,0,1.126002,0.014546,-0.001122,0.000062],
                [1,1,0.884043,0.017420,-0.000675,0.000069],
                [1,2,0.543071,0.022312,-0.000016,0.000087],
                [1,3,0.455269,0.037464,0.000053,0.000142],
                [1,4,0.399175,0.066277,0.000092,0.000248],
                [1,5,0.414872,0.118091,-0.000196,0.000411],
                [2,0,1.216332,0.025716,-0.001479,0.000110],
                [2,1,0.910763,0.023878,-0.000651,0.000097],
                [2,2,0.588476,0.027266,-0.000040,0.000107],
                [2,3,0.518375,0.040573,-0.000144,0.000154],
                [2,4,0.380038,0.067387,0.000255,0.000253],
                [2,5,0.315496,0.119609,0.000194,0.000443]
            ]
        self.WjetsToLNuNLOweightNjet2017 = [
                [0,0,1.205752,0.003731,-0.001030,0.000023],
                [0,1,1.271092,0.004954,-0.001083,0.000027],
                [0,2,1.051401,0.007353,-0.000513,0.000038],
                [0,3,0.716065,0.011309,0.000036,0.000056],
                [0,4,0.604683,0.021541,0.000221,0.000103],
                [0,5,0.597507,0.042707,-0.000073,0.000194],
                [1,0,1.357121,0.032714,-0.001750,0.000181],
                [1,1,1.300590,0.028820,-0.001623,0.000151],
                [1,2,1.012177,0.034249,-0.000642,0.000171],
                [1,3,0.708395,0.046308,-0.000186,0.000224],
                [1,4,0.741707,0.077736,-0.000826,0.000364],
                [1,5,0.750060,0.133655,-0.000920,0.000577],
                [2,0,1.330875,0.096135,-0.002060,0.000531],
                [2,1,1.205912,0.054905,-0.001236,0.000295],
                [2,2,1.034273,0.057342,-0.000573,0.000284],
                [2,3,0.855016,0.067325,-0.000550,0.000323],
                [2,4,0.572950,0.098955,0.000412,0.000465],
                [2,5,0.441597,0.171852,0.000354,0.000793]
            ]
        self.DYjetsToLLNLOweightNjet2017 = [
                [0,0,1.41667257,0.009776957,-0.001230725,6.5397E-05],
                [0,1,1.393552626,0.007503729,-0.001233994,4.53404E-05],
                [0,2,1.347529697,0.010084719,-0.001160407,5.50119E-05],
                [0,3,1.071661786,0.015302059,-0.000555197,7.72272E-05],
                [0,4,0.729464542,0.024818469,-2.05154E-05,0.000118966],
                [0,5,0.679199746,0.047661377,-0.000232446,0.000218842],
                [1,0,1.393949937,0.044851716,-0.001378941,0.000308352],
                [1,1,1.424499504,0.031378919,-0.001661213,0.000194317],
                [1,2,1.288121414,0.038578486,-0.001039063,0.000212253],
                [1,3,1.090142309,0.053589808,-0.000854487,0.00026954],
                [1,4,0.861623061,0.078617749,-0.000975446,0.000367765],
                [1,5,0.376231418,0.139875467,0.000346819,0.000616674],
                [2,0,1.681719504,0.183653773,-0.004419806,0.00120732],
                [2,1,1.546154404,0.073798177,-0.002363307,0.000443821],
                [2,2,1.453219102,0.059564281,-0.002123416,0.000333749],
                [2,3,1.078354952,0.068476135,-0.000859408,0.000353418],
                [2,4,0.949781804,0.095093814,-0.001050688,0.000471687],
                [2,5,0.709471838,0.160303941,-0.000694619,0.000742524]
            ]
        self.ZjetsToNuNuNLOweight2017V4 = [
            [
                [1.22677725339,0.00311575545047,-0.00112963223927,1.29282259793e-05],
                [1.29049031827,0.00461358046641,-0.00124756450272,1.94460346591e-05],
                [1.3636793139,0.00380872777152,-0.00147972942279,1.57968703605e-05],
                [1.38027222023,0.00423764008367,-0.0015447688251,1.76107315882e-05],
                [1.40308681015,0.00484583833964,-0.00161204424483,2.02253690332e-05],
                [1.43389363026,0.00571708611515,-0.00169512827031,2.39816294723e-05],
                [1.46158117198,0.00694411641241,-0.00175257156669,2.9361789442e-05],
                [1.49538841797,0.00864418260506,-0.00183752191762,3.68273172387e-05],
                [1.50673458945,0.0111005152962,-0.0018229532664,4.78534425631e-05],
                [1.5194558797,0.020088611979,-0.0018006013307,9.06177722205e-05],
                [1.53032921477,0.0147428635733,-0.0014809810881,6.78293429486e-05]
            ],
            [
                [1.16736127488,0.0121284431641,-0.00107476823455,4.98385067402e-05],
                [1.24110889013,0.0125519094754,-0.00128726343266,5.15656976743e-05],
                [1.24818352644,0.0139319992556,-0.0014163561105,5.72458255563e-05],
                [1.28446754174,0.0157218011294,-0.00152302151278,6.50078558206e-05],
                [1.3115449472,0.0176953960807,-0.0016501961541,7.29410162652e-05],
                [1.32636417558,0.0217764785483,-0.00163611932605,9.12341112311e-05],
                [1.34862919668,0.0384974157419,-0.00163294420101,0.000165453542905],
                [1.45334247929,0.0494593511794,-0.00204285127658,0.000213486537985],
                [1.493770589,0.0477271519291,-0.00203499950982,0.000207551408471],
                [1.48657664816,0.0656064096975,-0.00187600001609,0.000290582178555],
                [1.60806190872,0.0992119167923,-0.00255268993508,0.000456569756199]
            ],
            [
                [0.97687993581,0.0165063765063,-0.000790248462335,6.69412274479e-05],
                [1.0602744594,0.0178802047799,-0.00101398560128,7.23711369095e-05],
                [1.17901054066,0.0216841821638,-0.00129990153026,8.74063731256e-05],
                [1.2714700415,0.0257095173323,-0.0015415569768,0.000103859490534],
                [1.29002891386,0.0311320734108,-0.00154839673828,0.000126614718642],
                [1.42815664871,0.0399962299922,-0.00193721071673,0.000162485130853],
                [1.40042261369,0.0530796123521,-0.00185747235945,0.000220476663638],
                [1.44260749302,0.0730886547618,-0.00206568571207,0.000306464240409],
                [1.17947501595,0.158249999217,-0.00162620925245,0.000693126503846],
                [1.28595132842,0.120696407702,-0.00250134055258,0.000498755812061],
                [0.7557556876,0.245657709045,-0.00054160616066,0.00107978539421]
            ]
        ]

    def customInit(self, initVars):
        self.sample = initVars['sample']
        # this is just to separate V+jets from other samples, LOtoNLO weight does not use distinction between 
        # different V+jets components from sampleIndex, but explicit formula used in derivation of weight
        self.sampleIndexBase = self.sample.index[0] if type(self.sample.index) == list else self.sample.index
        self.sampleIndexBase = int(self.sampleIndexBase - (self.sampleIndexBase % 10))

        self.config = initVars['config']
        if not self.sample.isData():
            self.addBranch(self.branchName)
            self.addBranch(self.branchName + '_LHEVptV2b')
            self.addBranch(self.branchName + '_LHEVptV2b_p0_Up')
            self.addBranch(self.branchName + '_LHEVptV2b_p0_Down')
            self.addBranch(self.branchName + '_LHEVptV2b_p1_Up')
            self.addBranch(self.branchName + '_LHEVptV2b_p1_Down')

            self.addBranch(self.branchName + '_LHEVptV5')
            self.addBranch(self.branchName + '_LHEVptV5_p0_Up')
            self.addBranch(self.branchName + '_LHEVptV5_p0_Down')
            self.addBranch(self.branchName + '_LHEVptV5_p1_Up')
            self.addBranch(self.branchName + '_LHEVptV5_p1_Down')

            # those are TESTS only
            if self.addV3:
                self.addBranch(self.branchName + '_LHEVptV3')
                for njet in [0,1,2,3,4,5]:
                    self.addBranch(self.branchName + '_LHEVptV3_njet%d_p0_Up'%njet)
                    self.addBranch(self.branchName + '_LHEVptV3_njet%d_p0_Down'%njet)
                    self.addBranch(self.branchName + '_LHEVptV3_njet%d_p1_Up'%njet)
                    self.addBranch(self.branchName + '_LHEVptV3_njet%d_p1_Down'%njet)
            if self.addV4:
                self.addBranch(self.branchName + '_LHEVptV4')
                for deltaEtaBin in [0,1,2,3,4,5,6,7,8,9,10]:
                    self.addBranch(self.branchName + '_LHEVptV4_deta%d_p0_Up'%deltaEtaBin)
                    self.addBranch(self.branchName + '_LHEVptV4_deta%d_p0_Down'%deltaEtaBin)
                    self.addBranch(self.branchName + '_LHEVptV4_deta%d_p1_Up'%deltaEtaBin)
                    self.addBranch(self.branchName + '_LHEVptV4_deta%d_p1_Down'%deltaEtaBin)

            self.addBranch(self.branchName + '_2016')

            self.sampleTree = initVars['sampleTree']
            self.nAddJetsFormula       = "Sum$(GenJet_pt>30&&abs(GenJet_eta)<2.4)-2"
            self.deltaEtaFormula       = "min(int(abs(GenJet_eta[0]-GenJet_eta[1])/0.5),10)" 
            self.absGenDeltaEtaFormula = "abs(GenJet_eta[0]-GenJet_eta[1])"
            self.nGenBJetsFormula      = "Sum$(GenJet_pt>25 && abs(GenJet_eta)<2.4 && GenJet_hadronFlavour==5)"
            self.sampleTree.addFormula(self.nAddJetsFormula)
            self.sampleTree.addFormula(self.deltaEtaFormula)
            self.sampleTree.addFormula(self.absGenDeltaEtaFormula)
            self.sampleTree.addFormula(self.nGenBJetsFormula)

    def processEvent(self, tree):
        # if current entry has not been processed yet
        if not self.sample.isData() and not self.hasBeenProcessed(tree):
            self.markProcessed(tree)
            self._b(self.branchName)[0]                        = 1.0
            self._b(self.branchName + '_LHEVptV2b')[0]         = 1.0
            self._b(self.branchName + '_LHEVptV2b_p0_Up')[0]   = 1.0
            self._b(self.branchName + '_LHEVptV2b_p0_Down')[0] = 1.0
            self._b(self.branchName + '_LHEVptV2b_p1_Up')[0]   = 1.0
            self._b(self.branchName + '_LHEVptV2b_p1_Down')[0] = 1.0
            self._b(self.branchName + '_LHEVptV5')[0]          = 1.0
            self._b(self.branchName + '_LHEVptV5_p0_Up')[0]    = 1.0
            self._b(self.branchName + '_LHEVptV5_p0_Down')[0]  = 1.0
            self._b(self.branchName + '_LHEVptV5_p1_Up')[0]    = 1.0
            self._b(self.branchName + '_LHEVptV5_p1_Down')[0]  = 1.0
            self._b(self.branchName + '_2016')[0]              = 1.0

            # SFs below are tests only
            if self.addV3:
                self._b(self.branchName + '_LHEVptV3')[0]         = 1.0
                for njet in [0,1,2,3,4,5]:
                    for v in ['p0_Up', 'p0_Down', 'p1_Up', 'p1_Down']:
                        self._b(self.branchName + '_LHEVptV3_njet%d_%s'%(njet, v))[0] = 1.0
            if self.addV4: 
                self._b(self.branchName + '_LHEVptV4')[0]         = 1.0
                for n in [0,1,2,3,4,5,6,7,8,9,10]:
                    for v in ['p0_Up', 'p0_Down', 'p1_Up', 'p1_Down']:
                        self._b(self.branchName + '_LHEVptV4_deta%d_%s'%(n, v))[0] = 1.0

            if self.applies(tree):

                # V5
                # only derived for 2017 so far, can probably be applied to 2018 (but NOT to 2016)
                #LOtoNLOweight2017V5(self, vpt, nb, deltaEta, sampleIdentifier, var0=0.0, var1=0.0):
                
                vpt         = tree.LHE_Vpt
                nb          = self.sampleTree.evaluate(self.nGenBJetsFormula)
                gendeltaeta = self.sampleTree.evaluate(self.absGenDeltaEtaFormula)

                self._b(self.branchName + '_LHEVptV5')[0]         = self.LOtoNLOweight2017V5(vpt, nb, gendeltaeta, self.sample.identifier)
                self._b(self.branchName + '_LHEVptV5_p0_Up')[0]   = self.LOtoNLOweight2017V5(vpt, nb, gendeltaeta, self.sample.identifier, var0=1.0) 
                self._b(self.branchName + '_LHEVptV5_p0_Down')[0] = self.LOtoNLOweight2017V5(vpt, nb, gendeltaeta, self.sample.identifier, var0=-1.0)
                self._b(self.branchName + '_LHEVptV5_p1_Up')[0]   = self.LOtoNLOweight2017V5(vpt, nb, gendeltaeta, self.sample.identifier, var1=1.0)
                self._b(self.branchName + '_LHEVptV5_p1_Down')[0] = self.LOtoNLOweight2017V5(vpt, nb, gendeltaeta, self.sample.identifier, var1=-1.0)

                # <V5
                if tree.hJidx[0] > -1 and tree.hJidx[1] > -1:
                    etabb = abs(tree.Jet_eta[tree.hJidx[0]]-tree.Jet_eta[tree.hJidx[1]])
                else:
                    etabb = 0.0
                #njets = self.sampleIndex % 10

                if self.year == 2017:
                    # apply only one of them!

                    # eta bb derived from 2017 DY
                    self._b(self.branchName)[0] = 1.153 * self.LOtoNLOWeightBjetSplitEtabb2017(etabb, nb) 

                    # eta bb derived from 2016 DY
                    self._b(self.branchName + '_2016')[0] = 1.153 * self.LOtoNLOWeightBjetSplitEtabb(etabb, nb) 


    def applies(self, tree):
        isVJets = False
        sampleCat = self.sampleIndexBase

        # sync with AT: DYJetsToLL_M-4to50 not reweighted

        # Z+jets normal, W+jets normal, W+jets b-enriched
        if sampleCat in [4000,4100,4200,4300,4400,4500,4600,4700,5000,5100,5300,5400,11000,11100,11200,11300,11400,11500,11600,11700,15000,15100,15200,15300,15400,15500,15600]:
            isVJets = True
        
        # Z+jets b-enriched
        if sampleCat in [14000,12000,12100,12200,14100,14200,16000,16100,16200,16300]:
            isVJets = True

        return isVJets

    def LOtoNLOWeightBjetSplitEtabb(self, etabb, njets):
        SF = 1.
        if etabb < 5:
            if njets < 1:
                SF = 0.935422 + 0.0403162*etabb -0.0089026*etabb*etabb +0.0064324*etabb*etabb*etabb -0.000212443*etabb*etabb*etabb*etabb
            elif njets == 1:
                SF = 0.962415 +0.0329463*etabb -0.0414479*etabb*etabb +0.0240993*etabb*etabb*etabb -0.00278271*etabb*etabb*etabb*etabb
            elif njets >= 2:
                SF = (0.721265 -0.105643*etabb -0.0206835*etabb*etabb +0.00558626*etabb*etabb*etabb)*np.exp(0.450244*etabb)
        return SF

    def LOtoNLOWeightBjetSplitEtabb2017(self, etabb, njets):
        SF = 1.0
        etabb = min(etabb, 5.0)
        if njets < 1:
            SF = (0.958 + 0.0286 * etabb + 0.0014156 * etabb * etabb)
        elif njets == 1:
            SF = ((0.972 - 0.264 * etabb + 0.026919 * etabb * etabb) * np.exp(0.29901 * etabb))
        else:
            SF = (0.81 + 0.1493 * etabb - 0.000965976 * etabb * etabb)
        return SF
    
    def LOtoNLOWeightBjetSplitVpt2017(self, vpt, njets):
        SF = 1.0
        vpt = max(min(vpt,500.0),50.0)
        if njets < 1:
            SF = 1.1596-7.1563e-04*(vpt-5.000e+01)-1.1169e-06*(vpt-5.000e+01)**2 
        elif njets == 1:
            SF = 1.1153-7.3720e-04*(vpt-5.000e+01)-1.7232e-06*(vpt-5.000e+01)**2
        else:
            SF = 1.1667-8.9528e-04*(vpt-5.000e+01)-2.1150e-06*(vpt-5.000e+01)**2
        return SF
    
    def LOtoNLOWeightBjetSplitVpt2017preserveNormalization(self, vpt, njets):
        SF = 1.0
        vpt = max(min(vpt,500.0),50.0)
        if njets < 1:
            SF = 0.883*(1.1596-7.1563e-04*(vpt-5.000e+01)-1.1169e-06*(vpt-5.000e+01)**2) 
        elif njets == 1:
            SF = 0.926*(1.1153-7.3720e-04*(vpt-5.000e+01)-1.7232e-06*(vpt-5.000e+01)**2)
        else:
            SF = 0.887*(1.1667-8.9528e-04*(vpt-5.000e+01)-2.1150e-06*(vpt-5.000e+01)**2)
        return SF

    ## each of the 0b/1b/2b components relative contribution is taken from LO
    #def LOtoNLOWeightBjetSplitVpt2017V2_HFcompositionFromLO(self, vpt, njets, sampleIdentifier, var=0.0):
    #    SF = 1.0
    #    vpt = max(min(vpt,500.0),50.0)
    #    # exclude NLO samples itself
    #    if 'amc' not in sampleIdentifier:
    #        if any([x in sampleIdentifier for x in ['WJets', 'WBJets']]):
    #            if njets < 1:
    #                SF = 1.210 - (1.013e-3 + var*0.010e-3)*vpt 
    #            elif njets == 1:
    #                SF = 1.098 - (1.230e-3 + var*0.051e-3)*vpt 
    #            else:
    #                SF = 1.087 - (0.871e-3 + var*0.102e-3)*vpt
    #        elif any([x in sampleIdentifier for x in ['ZJets', 'ZBJets']]):
    #            pass
    #        elif any([x in sampleIdentifier for x in ['DYJets', 'DYBJets']]):
    #            if njets < 1:
    #                SF = 1.290 - (1.167e-3 + var*0.014e-3)*vpt 
    #            elif njets == 1:
    #                SF = 1.122 - (1.065e-3 + var*0.054e-3)*vpt 
    #            else:
    #                SF = 1.110 - (1.345e-3 + var*0.095e-3)*vpt
    #    return SF

    # relative 0b/1b/2b contributions from NLO, scaled to NNLO (LO x NNLO-k factor)
    def LOtoNLOWeightBjetSplitVpt2017V2(self, vpt, njets, sampleIdentifier, var0=0.0, var1=0.0):
        SF = 1.0
        vpt = max(min(vpt,500.0),50.0)
        # exclude NLO samples itself
        if 'amc' not in sampleIdentifier:
            if any([x in sampleIdentifier for x in ['WJets', 'WBJets']]):
                if njets < 1:
                    SF = (1.209 + var0*0.002) - (1.013e-3 + var1*0.010e-3)*vpt 
                elif njets == 1:
                    SF = (1.241 + var0*0.009) - (1.391e-3 + var1*0.057e-3)*vpt 
                else:
                    SF = (1.087 + var0*0.019) - (0.871e-3 + var1*0.103e-3)*vpt
            elif any([x in sampleIdentifier for x in ['ZJets', 'ZBJets']]):
                if njets < 1:
                    SF = (1.185 + var0*0.0036) - (1.248e-3 + var1*0.015e-3)*vpt 
                elif njets == 1:
                    SF = (1.109 + var0*0.0145) - (1.259e-3 + var1*0.062e-3)*vpt 
                else:
                    SF = (1.017 + var0*0.0253) - (1.173e-3 + var1*0.104e-3)*vpt
            elif any([x in sampleIdentifier for x in ['DYJets', 'DYBJets']]):
                if njets < 1:
                    SF = (1.288 + var0*0.002) - (1.166e-3 + var1*0.014e-3)*vpt 
                elif njets == 1:
                    SF = (1.231 + var0*0.006) - (1.168e-3 + var1*0.059e-3)*vpt 
                else:
                    SF = (1.219 + var0*0.013) - (1.478e-3 + var1*0.105e-3)*vpt
            if SF < 0:
                SF = 0
        return SF

    # normalizations from NLO  (not rescaled to LOxk-factor)
    # derived from NanoAOD without applying analysis selections
    def LOtoNLOWeightBjetSplitVpt2017V2b(self, vpt, njets, sampleIdentifier, var0=0.0, var1=0.0):
        SF = 1.0
        vpt = max(min(vpt,500.0),50.0)
        # exclude NLO samples itself
        if 'amc' not in sampleIdentifier:
            if any([x in sampleIdentifier for x in ['WJets', 'WBJets']]):
                if njets < 1:
                    SF = (1.374 + var0*1.457e-3) - (1.221e-3 + var1*4.986e-6)*vpt
                elif njets == 1:
                    SF = (1.202 + var0*7.968e-3) - (1.179e-3 + var1*3.067e-5)*vpt
                else:
                    SF = (1.000 + var0*1.406e-2) - (6.495e-4 + var1*5.384e-5)*vpt 
            elif any([x in sampleIdentifier for x in ['ZJets', 'ZBJets']]):
                if njets < 1:
                    SF = (1.365 + var0*1.127e-3) - (1.425e-3 + var1*4.218e-6)*vpt
                elif njets == 1:
                    SF = (1.273 + var0*6.169e-3) - (1.401e-3 + var1*2.528e-5)*vpt 
                else:
                    SF = (1.150 + var0*9.340e-3) - (1.222e-3 + var1*3.775e-5)*vpt 
            elif any([x in sampleIdentifier for x in ['DYJets', 'DYBJets']]):
                if njets < 1:
                    SF = (1.360 + var0*1.943e-3) - (1.400e-3 + var1*8.141e-6)*vpt
                elif njets == 1:
                    SF = (1.257 + var0*7.418e-3) - (1.371e-3 + var1*3.061e-5)*vpt 
                else:
                    SF = (1.177 + var0*1.195e-2) - (1.357e-3 + var1*4.697e-5)*vpt 
            if SF < 0:
                SF = 0
        return SF
    
    
    def LOtoNLOWeight2017V3(self, vpt, nb, njets, sampleIdentifier, var0=0.0, var1=0.0):
        SF = 1.0
        vpt = max(min(vpt,500.0),50.0)
        if 'amc' not in sampleIdentifier:
            sfTable = []
            if any([x in sampleIdentifier for x in ['ZJets', 'ZBJets']]):
                sfTable = self.ZjetsToNuNuNLOweightNjet2017
            elif any([x in sampleIdentifier for x in ['WJets', 'WBJets']]):
                sfTable = self.WjetsToLNuNLOweightNjet2017
            elif any([x in sampleIdentifier for x in ['DYJets', 'DYBJets']]):
                sfTable = self.DYjetsToLLNLOweightNjet2017

            nb = max(min(nb,2),0)
            njets = max(min(njets,5),0)
            for l in sfTable: 
                if l[0]==nb and l[1]==njets:
                    SF = (l[2] + var0*l[3]) + (l[4] + var1*l[5])*vpt
                    break
            if SF < 0:
                SF = 0
        return SF

    def LOtoNLOWeight2017V4(self, vpt, nb, deltaEtaBin, sampleIdentifier, var0=0.0, var1=0.0):
        SF = 1.0
        vpt = max(min(vpt,500.0),50.0)
        if 'amc' not in sampleIdentifier:
            nb = max(min(nb,2),0)
            deltaEtaBin = max(min(deltaEtaBin,10),0)
            coeffs = [1.0, 0.0, 0.0, 0.0]
            norm = 1.0
            if any([x in sampleIdentifier for x in ['ZJets', 'ZBJets']]):
                coeffs = self.ZjetsToNuNuNLOweight2017V4[nb][deltaEtaBin]
                # this factor is to keep total normalization fixed at LO * NNLO-k-factor 
                norm = 0.86365
            SF = norm * ((coeffs[0] + var0*coeffs[1]) + (coeffs[2] + var1*coeffs[3])*vpt)
            if SF < 0.01:
                SF = 0.01
        return SF


    # LOtoNLOweight2017V5
    #
    # 1) LHE_Vpt correction, normalization changing from LO(pdfwgt=True) to NLO
    #    important: LO samples with pdfwgt=False need an additional correction 
    # 2) residual deltaEta correction, normalization preserving. Only for ZJetsToNuNu, flat for other samples
    #
    # INPUTS:
    # vpt:              LHE_Vpt
    # nb:               Sum$(GenJet_pt>25 && abs(GenJet_eta)<2.4 && GenJet_hadronFlavour==5) 
    # deltaEta:         abs(GenJet_eta[0]-GenJet_eta[1]) 
    # sampleIdentifier: sampleIdentifier
    # var0/1:           0.0 for nominal, 1.0 for up, -1.0 for down
    #
    # OUTPUT:
    # MC SF to be applied to LO(pdfwgt=True) V+jets samples
    def LOtoNLOweight2017V5(self, vpt, nb, deltaEta, sampleIdentifier, var0=0.0, var1=0.0):
        SF = 1.0

        if any([x in sampleIdentifier for x in ['ZJets', 'ZBJets']]):
            vpt = max(min(vpt,400.0),0.0)

            if nb < 1:
                SF = (1.688e+00 + var0*1.727e-03)+(-1.785e-03 + var1*7.239e-06)*vpt
            elif nb == 1:
                SF = (1.575e+00 + var0*6.579e-03)+(-1.754e-03 + var1*2.721e-05)*vpt
            else:
                SF = (1.424e+00 + var0*1.035e-02)+(-1.539e-03 + var1*4.161e-05)*vpt

            x = max(min(abs(deltaEta),4.0),0.0) 
            if nb < 1:
                if vpt>=150.0 and vpt<170.0:
                    SF *= (8.902e-01)+(1.171e-01)*x+(-4.028e-02)*x**2+(5.413e-03)*x**3
                elif vpt>=170.0 and vpt<200.0:
                    SF *= (8.960e-01)+(1.125e-01)*x+(-4.521e-02)*x**2+(6.914e-03)*x**3
                elif vpt>=200.0 and vpt<250.0:
                    SF *= (9.055e-01)+(1.016e-01)*x+(-4.262e-02)*x**2+(6.732e-03)*x**3
                elif vpt>=250.0 and vpt<300.0:
                    SF *= (9.258e-01)+(6.067e-02)*x+(-2.653e-02)*x**2+(4.067e-03)*x**3
                elif vpt>=300.0 and vpt<400.0:
                    SF *= (9.446e-01)+(8.139e-02)*x+(-4.740e-02)*x**2+(8.163e-03)*x**3
            elif nb == 1:
                if vpt>=150.0 and vpt<170.0:
                    SF *= (9.509e-01)+(4.284e-02)*x+(-2.355e-02)*x**2+(6.062e-03)*x**3 
                elif vpt>=170.0 and vpt<200.0:
                    SF *= (9.476e-01)+(6.023e-02)*x+(-4.144e-02)*x**2+(9.121e-03)*x**3
                elif vpt>=200.0 and vpt<250.0:
                    SF *= (9.424e-01)+(8.237e-02)*x+(-4.904e-02)*x**2+(9.142e-03)*x**3
                elif vpt>=250.0 and vpt<300.0:
                    SF *= (9.838e-01)+(-5.520e-03)*x+(-1.425e-02)*x**2+(4.843e-03)*x**3
                elif vpt>=300.0 and vpt<400.0:
                    SF *= (1.012e+00)+(5.707e-03)*x+(-3.566e-02)*x**2+(1.008e-02)*x**3
            else:
                if vpt>=150.0 and vpt<170.0:
                    SF *= (8.635e-01)+(4.481e-02)*x+(8.009e-02)*x**2+(-1.883e-02)*x**3
                elif vpt>=170.0 and vpt<200.0:
                    SF *= (8.942e-01)+(6.551e-02)*x+(2.721e-02)*x**2+(-7.706e-03)*x**3
                elif vpt>=200.0 and vpt<250.0:
                    SF *= (8.709e-01)+(8.808e-02)*x+(1.561e-02)*x**2+(-6.648e-03)*x**3
                elif vpt>=250.0 and vpt<300.0:
                    SF *= (9.174e-01)+(2.229e-02)*x+(4.194e-02)*x**2+(-1.010e-02)*x**3
                elif vpt>=300.0 and vpt<400.0:
                    SF *= (9.712e-01)+(-1.882e-02)*x+(5.351e-02)*x**2+(-1.159e-02)*x**3

        elif any([x in sampleIdentifier for x in ['WJets', 'WBJets']]):
            vpt = max(min(vpt, 500.0), 100.0)

            # LHE_Vpt correction, normalization changing from LO(pdfwgt=True) to NLO
            # LO samples with pdfwgt=False need an additional correction before this one can be applied
            if nb < 1:
                SF = (1.628e+00 + var0*4.637e-03)+(-1.339e-03 + var1*2.052e-05)*vpt
            elif nb == 1:
                SF = (1.586e+00 + var0*2.657e-02)+(-1.531e-03 + var1*1.118e-04)*vpt
            else:
                SF = (1.440e+00 + var0*4.865e-02)+(-9.250e-04 + var1*2.034e-04)*vpt

        elif any([x in sampleIdentifier for x in ['DYJets', 'DYBJets']]):
            vpt = max(min(vpt, 500.0), 50.0)

            # LHE_Vpt correction, normalization changing from LO(pdfwgt=True) to NLO
            # LO samples with pdfwgt=False need an additional correction before this one can be applied
            if nb < 1:
                SF = (1.650e+00 + var0*2.112e-03)+(-1.707e-03 + var1*1.960e-05)*vpt  
            elif nb == 1:
                SF = (1.534e+00 + var0*9.535e-03)+(-1.485e-03 + var1*8.048e-05)*vpt
            else:
                SF = (1.519e+00 + var0*1.868e-02)+(-1.916e-03 + var1*1.396e-04)*vpt
        
        return SF


