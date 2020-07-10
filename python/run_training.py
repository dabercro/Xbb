#!/usr/bin/env python
from __future__ import print_function
from optparse import OptionParser
import ROOT
ROOT.gROOT.SetBatch(True)
from myutils import NewTreeCache as TreeCache
from myutils.sampleTree import SampleTree as SampleTree
from myutils import BetterConfigParser, ParseInfo
import resource
import os
import sys
import pickle
import glob
import shutil
import math

class MvaTrainingHelper(object):

    def __init__(self, config, mvaName):
        self.config = config
        self.factoryname = config.get('factory', 'factoryname')
        self.factorysettings = config.get('factory', 'factorysettings')
        if config.has_option('Directories', 'trainingSamples'):
            self.samplesPath = config.get('Directories', 'trainingSamples')
        else:
            self.samplesPath = config.get('Directories', 'MVAin')
        self.samplesInfo = ParseInfo(samples_path=self.samplesPath, config=self.config) 

        self.sampleFilesFolder = config.get('Directories', 'samplefiles')

        self.treeVarSet = config.get(mvaName, 'treeVarSet')
        self.MVAtype = config.get(mvaName, 'MVAtype')
        self.MVAsettings = config.get(mvaName,'MVAsettings')
        self.mvaName = mvaName

        VHbbNameSpace = config.get('VHbbNameSpace', 'library')
        ROOT.gSystem.Load(VHbbNameSpace)

        # variables
        self.MVA_Vars = {}
        self.MVA_Vars['Nominal'] = config.get(self.treeVarSet, 'Nominal').strip().split(' ')

        # samples
        backgroundSampleNames = eval(config.get(mvaName, 'backgrounds'))
        signalSampleNames = eval(config.get(mvaName, 'signals'))
        self.samples = {
            'BKG': self.samplesInfo.get_samples(backgroundSampleNames),
            'SIG': self.samplesInfo.get_samples(signalSampleNames),
        }

        self.treeCutName = config.get(mvaName, 'treeCut') if config.has_option(mvaName, 'treeCut') else mvaName
        self.treeCut = config.get('Cuts', self.treeCutName)

        self.TrainCut = config.get('Cuts', 'TrainCut') 
        self.EvalCut = config.get('Cuts', 'EvalCut')
        print("TRAINING CUT:", self.TrainCut)
        print("EVAL CUT:", self.EvalCut)

        self.globalRescale = 2.0
        
        self.trainingOutputFileName = 'mvatraining_{factoryname}_{region}.root'.format(factoryname=self.factoryname, region=mvaName)
        print("INFO: MvaTrainingHelper class created.")


    def prepare(self):

        self.trainingOutputFile = ROOT.TFile.Open(self.trainingOutputFileName, "RECREATE")
        # ----------------------------------------------------------------------------------------------------------------------
        # create TMVA factory
        # ----------------------------------------------------------------------------------------------------------------------
        self.factory = ROOT.TMVA.Factory(self.factoryname, self.trainingOutputFile, self.factorysettings)
        if self.trainingOutputFile and self.factory:
            print ("INFO: initialized MvaTrainingHelper.", self.factory) 
        else:
            print ("\x1b[31mERROR: initialization of MvaTrainingHelper failed!\x1b[0m") 

        # ----------------------------------------------------------------------------------------------------------------------
        # add sig/bkg x training/eval trees
        # ----------------------------------------------------------------------------------------------------------------------
        try:
            addBackgroundTreeMethod = self.factory.AddBackgroundTree
            addSignalTreeMethod = self.factory.AddSignalTree
            self.dataLoader = None
        except:
            print("oh no..")
            # the DataLoader wants to be called '.'
            self.dataLoader = ROOT.TMVA.DataLoader(".")
            addBackgroundTreeMethod = self.dataLoader.AddBackgroundTree
            addSignalTreeMethod = self.dataLoader.AddSignalTree

        if self.config.has_option('Weights','useSpecialWeight') and eval(self.config.get('Weights','useSpecialWeight')):
            print("\x1b[31mERROR: specialweight cannot be used with TMVA training, set it to false and add the DY_specialWeight to weightF!!\x1b[0m")
            raise Exception("SpecialWeightNotSupported")
        
        # DEBUG: restrict memory
        # resource.setrlimit(resource.RLIMIT_AS, (4.0*1024*1024*1024, 5.0*1024*1024*1024))

        self.sampleTrees = []
        for addTreeFcn, samples in [
                    [addBackgroundTreeMethod, self.samples['BKG']],
                    [addSignalTreeMethod, self.samples['SIG']]
                ]:
            for sample in samples:
                print ('*'*80,'\n%s\n'%sample,'*'*80)
                for additionalCut in [self.TrainCut, self.EvalCut]:
                    # cuts
                    sampleCuts = [sample.subcut]
                    if additionalCut:
                        sampleCuts.append(additionalCut)
                    # cut from the mva region
                    if self.treeCut:
                        sampleCuts.append(self.treeCut)

                    tc = TreeCache.TreeCache(
                            sample=sample,
                            cutList=sampleCuts,
                            inputFolder=self.samplesPath,
                            config=self.config,
                            debug=True
                        )
                    sampleTree = tc.getTree()
                    sampleTree.tree.SetCacheSize(32*1024)

                    # prevent garbage collection
                    self.sampleTrees.append(sampleTree)
                    if sampleTree:
                        treeScale = sampleTree.getScale(sample) * self.globalRescale

                        # only non-empty trees can be added
                        if sampleTree.tree.GetEntries() > 0:
                            addTreeFcn(sampleTree.tree, treeScale, ROOT.TMVA.Types.kTraining if additionalCut == self.TrainCut else ROOT.TMVA.Types.kTesting)
                            print('max mem used = %d'%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
                    else:
                        print ("\x1b[31mERROR: TREE NOT FOUND:", sample.name, " -> not cached??\x1b[0m")
                        raise Exception("CachedTreeMissing")

        if self.dataLoader:
            for var in self.MVA_Vars['Nominal']:
                self.dataLoader.AddVariable(var, 'D')
        else:
            for var in self.MVA_Vars['Nominal']:
                self.factory.AddVariable(var, 'D')

        return self

    # ----------------------------------------------------------------------------------------------------------------------
    # backup old .xml and .info files 
    # ----------------------------------------------------------------------------------------------------------------------
    def backupOldFiles(self):
        success = False
        MVAdir = self.config.get('Directories','vhbbpath')+'/python/weights/'
        backupDir = MVAdir + 'backup/'
        try:
            os.makedirs(backupDir)
        except:
            pass
        freeNumber = 1
        try:
            lastUsedBackupDirectories = sorted(glob.glob(backupDir + '/v*/'), key=lambda x: int(x.strip('/').split('/')[-1][1:]), reverse=True)
            freeNumber = 1 + int(lastUsedBackupDirectories[0].strip('/').split('/')[-1][1:]) if len(lastUsedBackupDirectories) > 0 else 1
        except Exception as e:
            print("\x1b[31mERROR: creating backup of MVA files failed!", e, "\x1b[0m")
            freeNumber = -1
        if freeNumber > -1:
            try:
                fileNamesToBackup = glob.glob(MVAdir + self.factoryname+'_'+self.mvaName + '.*')
                fileNamesToBackup += glob.glob(MVAdir + '/../mvatraining_MVA_ZllBDT_*.root')
                os.makedirs(backupDir + 'v%d/'%freeNumber)
                for fileNameToBackup in fileNamesToBackup:
                    shutil.copy(fileNameToBackup, backupDir + 'v%d/'%freeNumber)
                success = True
            except Exception as e:
                print("\x1b[31mERROR: creating backup of MVA files failed!", e, "\x1b[0m")
        return success


    def run(self):
        backupFiles = False
        try:
            backupFiles = eval(self.config.get('MVAGeneral', 'backupWeights'))
        except:
            pass
        if backupFiles:
            print('backing up old BDT files')
            self.backupOldFiles()
        # ----------------------------------------------------------------------------------------------------------------------
        # Execute TMVA
        # ----------------------------------------------------------------------------------------------------------------------
        self.factory.Verbose()
        print('max mem used = %d'%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        print('Execute TMVA: factory.BookMethod("%s", "%s", "%s")'%(self.MVAtype, self.mvaName, self.MVAsettings))
        print('max mem used = %d'%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        weightF = self.config.get('Weights','weightF')
        try:
            self.factory.BookMethod(self.MVAtype, self.mvaName, self.MVAsettings)
            print("ROOT 5 style TMVA found")
            self.factory.SetSignalWeightExpression(weightF)
            self.factory.SetBackgroundWeightExpression(weightF)
        except:
            print("ROOT 6 style TMVA found, using data loader object!!! >_<")
            print(" weights dir:", ROOT.TMVA.gConfig().GetIONames().fWeightFileDir)
            print(" data loader:", self.dataLoader)
            print(" type:       ", self.MVAtype)
            print(" name:       ", self.mvaName)
            print(" settings:   ", self.MVAsettings)
            ROOT.TMVA.gConfig().GetIONames().fWeightFileDir = 'weights'
            self.dataLoader.SetSignalWeightExpression(weightF)
            self.dataLoader.SetBackgroundWeightExpression(weightF)
            self.factory.BookMethod(self.dataLoader, self.MVAtype, self.mvaName, self.MVAsettings)
        sys.stdout.flush()
        print('Execute TMVA: TrainAllMethods')
        print('max mem used = %d'%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        self.factory.TrainAllMethods()
        sys.stdout.flush()
        print('Execute TMVA: TestAllMethods')
        print('max mem used = %d'%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        self.factory.TestAllMethods()
        sys.stdout.flush()
        print('Execute TMVA: EvaluateAllMethods')
        print('max mem used = %d'%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        self.factory.EvaluateAllMethods()
        sys.stdout.flush()
        print('Execute TMVA: output.Write')
        print('max mem used = %d'%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        self.trainingOutputFile.Close()
        return self

    def printInfo(self):
        #WRITE INFOFILE
        MVAdir = self.config.get('Directories','vhbbpath')+'/python/weights/'
        infofile = open(MVAdir+self.factoryname+'_'+self.mvaName+'.info','w')
        print ('@DEBUG: output infofile name')
        print (infofile)

        info=mvainfo(self.mvaName)
        info.factoryname=self.factoryname
        info.factorysettings=self.factorysettings
        info.MVAtype=self.MVAtype
        info.MVAsettings=self.MVAsettings
        info.weightfilepath=MVAdir
        info.path=self.samplesPath
        info.varset=self.treeVarSet
        info.vars=self.MVA_Vars['Nominal']
        pickle.dump(info,infofile)
        infofile.close()

    def getExpectedSignificance(self, tree, nBins, xMin, xMax, power=1.0, rescaleSig=1.0, rescaleBkg=1.0):
        hSIG = ROOT.TH1D("hSig","hSig",nBins,xMin,xMax)
        hBKG = ROOT.TH1D("hBkg","hBkg",nBins,xMin,xMax)
        print("INFO: GetEntries() = ", tree.GetEntries())
        if power != 1.0:
            print("INFO: rescale BDT score with power ", power)
        for event in tree:
            if power != 1.0:
                x = (getattr(event, self.mvaName)-xMin)/(xMax-xMin)
                if x<0:
                    x=0
                if x>0.999999:
                    x=0.999999
                value = math.pow(x, power)*(xMax-xMin)+xMin
            else:
                value = max(min(getattr(event, self.mvaName),xMax-0.00001),xMin)

            weight = event.weight
            if event.classID == 1:
                hSIG.Fill(value, weight * rescaleSig)
            else:
                hBKG.Fill(value, weight * rescaleBkg)
        ssbSum = 0.0
        sSum = 0
        bSum = 0
        sbTableFormat = "{bin: <16}{signal: <16}{background: <16}{ssb: <16}"
        print("---- nBins =", nBins, " from ", xMin, "..", xMax, "-----")
        print(sbTableFormat.format(bin="bin", signal="signal", background="background", ssb="S/sqrt(S+B)"))
        for i in range(nBins):
            ssbSum += hSIG.GetBinContent(1+i)*hSIG.GetBinContent(1+i)/(hSIG.GetBinContent(1+i) + hBKG.GetBinContent(1+i)) if (hSIG.GetBinContent(1+i) + hBKG.GetBinContent(1+i)) > 0 else 0
            sSum += hSIG.GetBinContent(1+i)
            bSum += hBKG.GetBinContent(1+i)
            ssb = hSIG.GetBinContent(1+i)/math.sqrt(hSIG.GetBinContent(1+i) + hBKG.GetBinContent(1+i)) if (hSIG.GetBinContent(1+i) + hBKG.GetBinContent(1+i)) > 0 else 0
            print(sbTableFormat.format(bin=i, signal=round(hSIG.GetBinContent(1+i),2), background=round(hBKG.GetBinContent(1+i),2), ssb=round(ssb,3)))
        expectedSignificance = math.sqrt(ssbSum)
        print(sbTableFormat.format(bin="SUM", signal=round(sSum,1), background=round(bSum,1), ssb="\x1b[34mZ=%1.3f\x1b[0m"%expectedSignificance))
        print("-"*40)
        hSIG.Delete()
        hBKG.Delete()
        return expectedSignificance, sSum, bSum

    def estimateExpectedSignificance(self):
        print("INFO: open ", self.trainingOutputFileName)
        rootFile = ROOT.TFile.Open(self.trainingOutputFileName, "READ")
        print("INFO: ->", rootFile)
        testTree = rootFile.Get('./TestTree')

        # run a few tests with different binnings and rescaling of BDT score
        self.getExpectedSignificance(testTree, 15, -0.8, 1.0)
        self.getExpectedSignificance(testTree, 15, -0.8, 0.9)
        self.getExpectedSignificance(testTree, 15, -0.8, 0.8)
        self.getExpectedSignificance(testTree, 15, -0.8, 0.75)
        self.getExpectedSignificance(testTree, 15, -0.8, 0.7)
        self.getExpectedSignificance(testTree, 15, -0.8, 0.8, power=0.5)
        self.getExpectedSignificance(testTree, 15, -0.8, 0.8, power=0.33)
        self.getExpectedSignificance(testTree, 15, -0.8, 0.8, power=1.5)
        self.getExpectedSignificance(testTree, 15, -0.8, 0.8, power=2.0)

        # close to nominal binning
        print("---- ~nominal TEST -----")
        esTest, sTest, bTest = self.getExpectedSignificance(testTree, 15, -0.8, 0.8)
        print("---- ~nominal TRAINING (without correct normalization) -----")
        trainTree = rootFile.Get('./TrainTree')
        esTrain, sTrain, bTrain = self.getExpectedSignificance(trainTree, 15, -0.8, 0.8)

        # the tree ./TrainTree contains the input events for training AFTER re-balancing the classes
        # therefore for SIG/BKG separately the normalization is fixed to the one of the TEST events
        rescaleSig = 1.0*sTest/sTrain
        rescaleBkg = 1.0*bTest/bTrain
        print("---- ~nominal TRAINING -----")
        trainTree = rootFile.Get('./TrainTree')
        esTrain, sTrain, bTrain = self.getExpectedSignificance(trainTree, 15, -0.8, 0.8, rescaleSig=rescaleSig, rescaleBkg=rescaleBkg)

    def getbdtHistogram(self, tree):
        hSIG = ROOT.TH1D("hSig","TMVA overtraining check for classifier: %s"%self.mvaName,40,-1,1)
        hBKG = ROOT.TH1D("hBkg","TMVA overtraining check for classifier: %s"%self.mvaName,40,-1,1)
        print("INFO: GetEntries() = ", tree.GetEntries())
        for event in tree:
            value = getattr(event, self.mvaName)

            if event.classID == 1:
                hSIG.Fill(value)
            else:
                hBKG.Fill(value)
        return [hSIG, hBKG]

    def setTMVASyle(self):
        # style
        self.hSIGtest.SetLineColor(ROOT.TColor.GetColor("#0000ee"))
        self.hSIGtest.SetLineWidth(1)
        self.hSIGtest.SetFillStyle(1001)
        self.hSIGtest.SetFillColor(ROOT.TColor.GetColor("#7d99d1"))
        self.hSIGtest.SetTitle("TMVA overtraining check for classifier: %s"%self.mvaName )

        self.hBKGtest.SetLineColor(ROOT.TColor.GetColor("#ff0000"))
        self.hBKGtest.SetLineWidth(1)
        self.hBKGtest.SetFillStyle(3554)
        self.hBKGtest.SetFillColor(ROOT.TColor.GetColor("#ff0000"))
        self.hBKGtest.SetTitle(self.hSIGtest.GetTitle())


        self.hSIGtrain.SetMarkerColor(self.hSIGtest.GetLineColor())
        self.hSIGtrain.SetMarkerSize(0.7)
        self.hSIGtrain.SetMarkerStyle(20)
        self.hSIGtrain.SetLineWidth(1)
        self.hSIGtrain.SetLineColor(self.hSIGtest.GetLineColor())
        self.hSIGtrain.SetTitle(self.hSIGtest.GetTitle())

        self.hBKGtrain.SetMarkerColor(self.hBKGtest.GetLineColor())
        self.hBKGtrain.SetMarkerSize(0.7)
        self.hBKGtrain.SetMarkerStyle(20)
        self.hBKGtrain.SetLineWidth(1)
        self.hBKGtrain.SetLineColor(self.hBKGtest.GetLineColor())
        self.hBKGtrain.SetTitle(self.hSIGtest.GetTitle())

        TMVAStyle = ROOT.TStyle(ROOT.gROOT.GetStyle("Plain"))# // our style is based on Plain
        TMVAStyle.SetName("TMVA")
        TMVAStyle.SetTitle("TMVA style based on \"Plain\" with modifications defined in tmvaglob.C")
        ROOT.gROOT.GetListOfStyles().Add(TMVAStyle)
        ROOT.gROOT.SetStyle("TMVA")
         	
        TMVAStyle.SetLineStyleString( 5, "[52 12]" )
        TMVAStyle.SetLineStyleString( 6, "[22 12]" )
        TMVAStyle.SetLineStyleString( 7, "[22 10 7 10]" )
        
        UsePaperStyle = False
    
        #// the pretty color palette of old
        TMVAStyle.SetPalette((18 if UsePaperStyle else 1))
    
        #// use plain black on white colors
        TMVAStyle.SetFrameBorderMode(0)
        TMVAStyle.SetCanvasBorderMode(0)
        TMVAStyle.SetPadBorderMode(0)
        TMVAStyle.SetPadColor(0)
        TMVAStyle.SetFillStyle(0)
    
        TMVAStyle.SetLegendBorderSize(0)
    
        c_TitleBox = ROOT.TColor.GetColor( "#5D6B7D" )
        c_TitleText = ROOT.TColor.GetColor( "#FFFFFF" )
        c_TitleBorder = ROOT.TColor.GetColor( "#7D8B9D" )
        c_FrameFill = ROOT.TColor.GetColor( "#fffffd" )
        c_Canvas = ROOT.TColor.GetColor( "#f0f0f0" )


        TMVAStyle.SetTitleFillColor( c_TitleBox )
        TMVAStyle.SetTitleTextColor( c_TitleText )
        TMVAStyle.SetTitleBorderSize( 1 )
        TMVAStyle.SetLineColor( c_TitleBorder )
        if not UsePaperStyle:
            TMVAStyle.SetFrameFillColor( c_FrameFill )
            TMVAStyle.SetCanvasColor( c_Canvas )
    
        #// set the paper & margin sizes
        TMVAStyle.SetPaperSize(20,26)
        TMVAStyle.SetPadTopMargin(0.10)
        TMVAStyle.SetPadRightMargin(0.05)
        TMVAStyle.SetPadBottomMargin(0.11)
        TMVAStyle.SetPadLeftMargin(0.12)
    
        #// use bold lines and markers
        TMVAStyle.SetMarkerStyle(21)
        TMVAStyle.SetMarkerSize(0.3)
        TMVAStyle.SetHistLineWidth(2)
        TMVAStyle.SetLineStyleString(2,"[12 12]") #// postscript dashes
    
        #// do not display any of the standard histogram decorations
        TMVAStyle.SetOptTitle(1)
        TMVAStyle.SetTitleH(0.052)
    
        TMVAStyle.SetOptStat(0)
        TMVAStyle.SetOptFit(0)
    
        #// put tick marks on top and RHS of plots
        TMVAStyle.SetPadTickX(1)
        TMVAStyle.SetPadTickY(1)

    def nomrmaliseHist(self,hSIG, hBKG):
        if (hSIG.GetSumw2N() == 0):
            hSIG.Sumw2()
        if (hBKG and hBKG.GetSumw2N() == 0): 
            hBKG.Sumw2()
     
        if(hSIG.GetSumOfWeights()!=0): 
            dx = (hSIG.GetXaxis().GetXmax() - hSIG.GetXaxis().GetXmin())/hSIG.GetNbinsX()
            hSIG.Scale(1.0/hSIG.GetSumOfWeights()/dx)
        if (hBKG != 0 and hBKG.GetSumOfWeights()!=0):
            dx = (hBKG.GetXaxis().GetXmax() - hBKG.GetXaxis().GetXmin())/hBKG.GetNbinsX()
            hBKG.Scale( 1.0/hBKG.GetSumOfWeights()/dx )

    def drawOvertraining(self):

        #normalise histograms
        self.nomrmaliseHist(self.hSIGtest, self.hBKGtest)
        self.nomrmaliseHist(self.hSIGtrain, self.hBKGtrain)

        c = ROOT.TCanvas("canvas1", "TMVA comparison %s"%self.mvaName, 0, 200, 600, 468) 

        # frame limits (choose judicuous x range)
        nrms = 10
        xmin = ROOT.TMath.Max(ROOT.TMath.Min(self.hSIGtest.GetMean() - nrms*self.hSIGtest.GetRMS(), self.hBKGtest.GetMean() - nrms*self.hBKGtest.GetRMS() ),self.hSIGtest.GetXaxis().GetXmin() )
        xmax = ROOT.TMath.Min(ROOT.TMath.Max(self.hSIGtest.GetMean() + nrms*self.hSIGtest.GetRMS(), self.hBKGtest.GetMean() + nrms*self.hBKGtest.GetRMS()), self.hSIGtest.GetXaxis().GetXmax())
        ymin = 0
        maxMult = 1.3
        #maxMult = (htype == CompareType) ? 1.3 : 1.2
        ymax = ROOT.TMath.Max(self.hSIGtest.GetMaximum(), self.hBKGtest.GetMaximum())*maxMult
        ymax = ROOT.TMath.Max(ymax,ROOT.TMath.Max(self.hSIGtrain.GetMaximum(), self.hBKGtrain.GetMaximum())*maxMult)
        #print ('ymax is', ymax)
        #print (self.hSIGtest.GetMaximum())
        #print (self.hBKGtest.GetMaximum())
        #print (self.hSIGtrain.GetMaximum())
        #print (self.hBKGtrain.GetMaximum())
        #sys.exit()
   
        # build a frame
        nb = 500
        hFrameName = "frame" + self.mvaName
        #o = ROOT.gROOT.FindObject(hFrameName)
        frame = ROOT.TH2F(hFrameName, self.hSIGtest.GetTitle(), nb, xmin, xmax, nb, ymin, ymax )
        frame.GetXaxis().SetTitle(self.mvaName + " response")
        frame.GetYaxis().SetTitle("(1/N) dN^{ }/^{ }dx")

        #TMVAGlob.SetFrameStyle( frame )
        frame.SetLabelOffset( 0.012, "X" )
        frame.SetLabelOffset( 0.012, "Y" )
        frame.GetXaxis().SetTitleOffset( 1.25 )
        frame.GetYaxis().SetTitleOffset( 1.22 )
        frame.GetXaxis().SetTitleSize( 0.045)
        frame.GetYaxis().SetTitleSize( 0.045)
        frame.GetXaxis().SetLabelSize( 0.04)
        frame.GetYaxis().SetLabelSize( 0.04)

        #// global style settings
        ROOT.gPad.SetTicks()
        ROOT.gPad.SetLeftMargin  ( 0.108)
        ROOT.gPad.SetRightMargin ( 0.050)
        ROOT.gPad.SetBottomMargin( 0.120)
   
        # eventually: draw the frame
        frame.Draw()  
    
        c.GetPad(0).SetLeftMargin(0.105 )
        frame.GetYaxis().SetTitleOffset( 1.2 )

        # Draw legend               
        legend = ROOT.TLegend(c.GetLeftMargin(), 1 - c.GetTopMargin() - 0.12, c.GetLeftMargin() + 0.40, 1 - c.GetTopMargin() )
        legend.SetFillStyle(1)
        legend.AddEntry(self.hSIGtest,"Signal"     + " (test sample)", "F")
        legend.AddEntry(self.hBKGtest,"Background" + " (test sample)", "F")
        legend.SetBorderSize(1)
        legend.SetMargin(0.2)
        legend.Draw("same")

        legend2= ROOT.TLegend( 1 - c.GetRightMargin() - 0.42, 1 - c.GetTopMargin() - 0.12, 1 - c.GetRightMargin(), 1 - c.GetTopMargin() )
        legend2.SetFillStyle(1)
        legend2.SetBorderSize(1)
        legend2.AddEntry(self.hSIGtrain,"Signal (training sample)","P")
        legend2.AddEntry(self.hBKGtrain,"Background (training sample)","P")
        legend2.SetMargin( 0.1 )
        legend2.Draw("same")

        self.setTMVASyle()

        self.hSIGtest.Draw('samehist')
        self.hBKGtest.Draw('samehist')
        self.hSIGtrain.Draw('e1same')
        self.hBKGtrain.Draw('e1same')

        #perform K-S test
        print("--- Perform Kolmogorov-Smirnov tests")
        #//Double_t kolS = sig->KolmogorovTest( self.hSIGtrain, "X" );
        #//Double_t kolB = bgd->KolmogorovTest( bgdOv, "X" );
        kolS = self.hSIGtest.KolmogorovTest( self.hSIGtrain);
        kolB = self.hBKGtest.KolmogorovTest( self.hBKGtrain);
        print ("--- Goodness of signal (background) consistency: " + str(kolS) + " (" + str(kolB) + ")")

        probatext = "Kolmogorov-Smirnov test: signal (background) probability = % 5.3g (%5.3g)"% (kolS, kolB)
        tt = ROOT.TText(0.12, 0.74, probatext)
        tt.SetNDC()
        tt.SetTextSize(0.032)
        tt.AppendPad()

        # redraw axes
        frame.Draw("sameaxis")

        #/text for overflows
        nbin = self.hSIGtest.GetNbinsX()
        dxu  = self.hSIGtest.GetBinWidth(0)
        dxo  = self.hSIGtest.GetBinWidth(nbin+1)
        uoflow =  "U/O-flow (S,B): (%.1f, %.1f)%% / (%.1f, %.1f)%%"% (self.hSIGtest.GetBinContent(0)*dxu*100, self.hBKGtest.GetBinContent(0)*dxu*100, self.hSIGtest.GetBinContent(nbin+1)*dxo*100, self.hBKGtest.GetBinContent(nbin+1)*dxo*100)
        t = ROOT.TText( 0.975, 0.115, uoflow )
        t.SetNDC()
        t.SetTextSize( 0.030 )
        t.SetTextAngle( 90 )
        t.AppendPad()    
   
        # update canvas
        c.Update()


        MVAdir = self.config.get('Directories','vhbbpath')+'/python/weights/'
        c.SaveAs(MVAdir+'overtraining%s.pdf'%self.mvaName)
        print ('I saved the canvase in', MVAdir+'overtraining%s.pdf'%self.mvaName)

    def saveOvertrainingPlots(self):
        print("INFO: open ", self.trainingOutputFileName)
        rootFile = ROOT.TFile.Open(self.trainingOutputFileName, "READ")
        print("INFO: ->", rootFile)

        self.hSIGtest = rootFile.Get('./Method_%s/%s/MVA_%s_S'%(self.mvaName,self.mvaName,self.mvaName))
        self.hBKGtest = rootFile.Get('./Method_%s/%s/MVA_%s_B'%(self.mvaName,self.mvaName,self.mvaName))
        self.hSIGtrain = rootFile.Get('./Method_%s/%s/MVA_%s_Train_S'%(self.mvaName,self.mvaName,self.mvaName))
        self.hBKGtrain = rootFile.Get('./Method_%s/%s/MVA_%s_Train_B'%(self.mvaName,self.mvaName,self.mvaName))
        print("./Method_%s/%s/MVA_%s_Train_B"%(self.mvaName,self.mvaName,self.mvaName))
        self.drawOvertraining()



# read arguments
argv = sys.argv
parser = OptionParser()
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                          help="Verbose mode.")
parser.add_option("-C", "--config", dest="config", default=[], action="append",
                      help="configuration file")
parser.add_option("-t","--trainingRegions", dest="trainingRegions", default='',
                      help="cut region identifier")
parser.add_option("-s", "--expectedSignificance" ,action="store_true", dest="expectedSignificance", default=False,
                          help="Compute estimate for expected significance (without systematics)")
(opts, args) = parser.parse_args(argv)
if opts.config =="":
        opts.config = ["config"]

# Import after configure to get help message
from myutils import BetterConfigParser, mvainfo, ParseInfo

# load config
config = BetterConfigParser()
config.read(opts.config)

# initialize
trainingRegions = opts.trainingRegions.split(',')
if len(trainingRegions) > 1:
    print ("ERROR: not implemented!")
    exit(1)
for trainingRegion in trainingRegions:
    th = MvaTrainingHelper(config=config, mvaName=trainingRegion)
    if opts.expectedSignificance:
        th.estimateExpectedSignificance()
    else:
        th.prepare().run()
        th.printInfo()
        try:
            th.estimateExpectedSignificance()
        except:
            pass
        try:
            th.saveOvertrainingPlots()
        except: 
            pass

