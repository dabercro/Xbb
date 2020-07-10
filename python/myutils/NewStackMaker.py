from __future__ import print_function
import ROOT 
ROOT.gROOT.SetBatch(True)
import TdrStyles
import os
import sys
import array
import time
import subprocess
import math
import numpy as np
from copy import deepcopy

from Ratio import getRatio
from NewHistoMaker import NewHistoMaker as HistoMaker
from sampleTree import SampleTree as SampleTree
from XbbTools import XbbTools

# ------------------------------------------------------------------------------
# produces histograms from trees with HistoMaker, groups them, and draws a
# stacked histogram
# ------------------------------------------------------------------------------
class NewStackMaker:

    def readConfig(self, section, option, default=None):
        if self.config.has_section(section) and self.config.has_option(section, option):
            return eval(self.config.get(section, option))
        else:
            return default
    
    def readConfigStr(self, section, option, default=None):
        if self.config.has_section(section) and self.config.has_option(section, option):
            return self.config.get(section, option)
        else:
            return default

    def __init__(self, config, var, region, SignalRegion, setup=None, subcut='', title='CMS', configSectionPrefix='Plot'):
        self.debug = 'XBBDEBUG' in os.environ
        self.garbage = []
        self.config = config
        self.saveShapes = True
        self.var = var
        self.region = region
        self.outputFolder = '.'
        self.configSectionPrefix = configSectionPrefix
        self.configSection = '{prefix}:{region}'.format(prefix=self.configSectionPrefix, region=region)
        if not self.config.has_section(self.configSection):
            print("\x1b[31mWARNING: no section '%s' in config, using defaults.\x1b[0m"%('Plot:%s'%region))
            self.configSection = None
            #raise Exception("ConfigError")

        self.plotVarSection = 'plotDef:%s'%var
        self.hasPlotVarSection = self.config.has_section(self.plotVarSection)

        self.dataGroupName = 'DATA'
        self.anaTag        = self.config.get("Analysis", "tag")
        self.subcut        = subcut
        self.forceLog      = None
        self.normalize     = self.readConfig(self.configSection, 'Normalize', False)
        self.log           = self.readConfig(self.configSection, 'log', False)
        self.blind         = self.readConfig(self.configSection, 'blind', False)

        self.mcUncertaintyLegend = self.readConfig('Plot_general', 'mcUncertaintyLegend', 'MC uncert. (stat.)')
        self.drawMCErrorForNormalizedPlots = self.readConfig('Plot_general', 'drawMCErrorForNormalizedPlots', False)
        self.asimovData    = self.readConfig('Plot_general', 'AsimovData', False)

        self.xAxis = self.readConfigStr(self.plotVarSection, 'xAxis', '')
        self.yAxis = self.readConfigStr(self.plotVarSection, 'yAxis', '')

        if self.hasPlotVarSection: 
            if self.config.has_option(self.plotVarSection, 'log'):
                self.log = eval(self.config.get(self.plotVarSection, 'log'))

        self.is2D = True if self.yAxis else False
        self.typLegendDict = self.readConfig('Plot_general','typLegendDict', {})
        self.legendEntries = []
        self.plotLabels = {}
        if setup is None:
            self.setup = [x.strip() for x in self.config.get('Plot_general', 'setupLog' if self.log else 'setup').split(',') if len(x.strip()) > 0]
        else:
            self.setup = setup

        self.plotTitle = self.readConfig('Plot_general', 'title', title) 

        self.rebin = 1
        self.histogramOptions = {
                'rebin': 1,
                'var': self.var,
                }
        self.additionalTextLines = [""]
        if self.config.has_option('Plot_general', 'additionalText'):
            aText =  eval(self.config.get('Plot_general', 'additionalText'))
            if type(aText) == list:
                self.additionalTextLines += aText
            else:
                self.additionalTextLines.append(aText)
        if self.config.has_section(self.configSection) and self.config.has_option(self.configSection, 'additionalText'):
            aText = eval(self.config.get(self.configSection, 'additionalText'))
            if type(aText) == list:
                self.additionalTextLines += aText
            else:
                self.additionalTextLines.append(aText)
        if self.hasPlotVarSection:
            if self.config.has_option(self.plotVarSection, 'additionalText'):
                aText = eval(self.config.get(self.plotVarSection, 'additionalText'))
                if type(aText) == list:
                    self.additionalTextLines += aText
                else:
                    self.additionalTextLines.append(aText)


        # general event by event weight which is applied to all samples
        #  for special plots of weights itself, weightF can be defined in the plot definition
        if self.config.has_option('plotDef:%s'%self.var,'weightF'):
            self.histogramOptions['weight'] = self.config.get('plotDef:%s'%self.var,'weightF')
        elif self.config.has_option('plotDef:%s'%self.var,'weight'):
            self.histogramOptions['weight'] = self.config.get('plotDef:%s'%self.var,'weight')
        elif self.config.has_option('Weights','weightF'):
            self.histogramOptions['weight'] = self.config.get('Weights','weightF')
        else:
            self.histogramOptions['weight'] = None

        optionNames = {
                    'treeVar': 'relPath',
                    'rebin': 'rebin',
                    'xAxis': 'xAxis',
                    'yAxis': 'yAxis',
                    'drawOption': 'drawOption',
                    'drawOptionData': 'drawOptionData',
                    'draw': 'draw',
                    'binList': 'binList',
                    'plotEqualSize': 'plotEqualSize',
                    'rebinFlat': 'rebinFlat',
                    'visualizeBlindCutThreshold': 'visualizeBlindCutThreshold',
                    'min': ['min', 'minX'],
                    'max': ['max', 'maxX'],
                    'blindCut': 'blindCut',
                    'minX': ['minX', 'min'],
                    'minY': ['minY', 'min'],
                    'minZ': ['minZ'],
                    'maxX': ['maxX', 'max'],
                    'maxY': ['maxY', 'max'],
                    'maxZ': ['maxZ'],
                    'nBins': ['nBins', 'nBinsX'],
                    'nBinsX': ['nBinsX', 'nBins'],
                    'nBinsY': ['nBinsY', 'nBins'],
                    'fractions': 'fractions',
                }
        numericOptions = ['rebin', 'min', 'minX', 'minY', 'maxX', 'maxY', 'nBins', 'nBinsX', 'nBinsY', 'minZ', 'maxZ']
        evalOptions = ['binList', 'plotEqualSize','fractions','rebinFlat']
        for optionName, configKeys in optionNames.iteritems():
            # use the first available option from the config, first look in region definition, afterwards in plot definition
            configKeysList = configKeys if type(configKeys) == list else [configKeys]
            for configKey in configKeysList:
                if self.config.has_section(self.configSection) and self.config.has_option(self.configSection, configKey):
                    self.histogramOptions[optionName] = self.config.get(self.configSection, configKey)
                    #print(self.configSection, configKey,self.histogramOptions[optionName])
                    break
                elif self.config.has_option('plotDef:%s'%var, configKey):
                    self.histogramOptions[optionName] = self.config.get('plotDef:%s'%var, configKey)
                    #print('plotDef:%s'%var, configKey,self.histogramOptions[optionName])
                    break
            # convert numeric options to float/int
            if optionName in numericOptions and optionName in self.histogramOptions and type(self.histogramOptions[optionName]) == str:
                self.histogramOptions[optionName] = float(self.histogramOptions[optionName]) if ('.' in self.histogramOptions[optionName] or 'e' in self.histogramOptions[optionName]) else int(self.histogramOptions[optionName])
       
        if self.config.has_section('plotDef:%s'%var) and self.config.has_option('plotDef:%s'%var, 'rebinMethod'):
            if 'binList' in self.histogramOptions:
                del self.histogramOptions['binList']
                print("DEBUG: rebinMethod present, ignoring any bin lists given.")

        # evaluate options given as python code
        for evalOption in evalOptions: 
            if evalOption in self.histogramOptions and type(evalOption) == str:
                self.histogramOptions[evalOption] = eval(self.histogramOptions[evalOption])

        # region/variable specific blinding cut
        if self.config.has_option(self.configSection, 'blindCuts'):
            blindCuts = eval(self.config.get(self.configSection, 'blindCuts'))
            if self.var in blindCuts:
                self.histogramOptions['blindCut'] = blindCuts[self.var]
                if '{var}' in self.histogramOptions['blindCut']:
                    self.histogramOptions['blindCut'] = self.histogramOptions['blindCut'].format(var=self.histogramOptions['treeVar'])
                print("\x1b[31mINFO: for region {region} var {var} using the blinding cut: {cut}\x1b[0m".format(region=self.region, var=self.var, cut=self.histogramOptions['blindCut']))

        self.groups = {}
        self.histograms = []
        self.legends = {}
        self.plotTexts = {}
        self.collectedObjects = []
        self.dataTitle = 'Data'
        self.maxRatioUncert = 0.5
        self.lumi = self.config.get('General','lumi')
        self.ratioError = None
        self.ratioPlot = None
        if SignalRegion:
            self.maxRatioUncert = 1000.
        self.outputTeX = False
        self.outputFileTemplate = "{outputFolder}/{prefix}.{ext}"
        try:
            self.outputFileFormats = [x.strip() for x in config.get('Plot_general','outputFormats').split(',') if len(x.strip())>0] 
        except:
            self.outputFileFormats = ["png"]
        
        if 'tex' in self.outputFileFormats:
            self.outputFileFormats = [x for x in self.outputFileFormats if x != 'tex']
            self.outputTeX = True

        self.plotTextMarginLeft = 0.16

        if 'treeVar' not in self.histogramOptions: 
            print ("ERROR: treeVar ", self.var)
            raise Exception ("config error")
        if self.debug:
            print ("INFO: StackMaker initialized!", self.histogramOptions['treeVar'], " min=", self.histogramOptions['minX'], " max=", self.histogramOptions['maxX'], "nBins=", self.histogramOptions['nBins'] if 'nBins' in self.histogramOptions else '-')

    def setPlotText(self, text):
        if type(text) == list:
            self.additionalTextLines = text
        elif type(text) == str:
            self.additionalTextLines = [text]
        else:
            print("ERROR: can't set plto text, unknown type:", type(text))

    # ------------------------------------------------------------------------------
    # draw text
    # ------------------------------------------------------------------------------
    @staticmethod
    def myText(txt="CMS Preliminary", ndcX=0.0, ndcY=0.0, size=0.8, color=None):
        ROOT.gPad.Update()
        text = ROOT.TLatex()
        text.SetNDC()
        if color:
            text.SetTextColor(color)
        text.SetTextSize(text.GetTextSize()*size)
        text.DrawLatex(ndcX,ndcY,txt)
        return text

    # ------------------------------------------------------------------------------
    # create histogram out of a tree
    # ------------------------------------------------------------------------------
    def addSampleTree(self, sample, sampleTree, groupName, cut='1'):
        if 'treeVar' not in self.histogramOptions:
            print("ERROR:", self.var, " not found.")
        print ("INFO: var=", self.var, "-> treeVar=\x1b[34m", self.histogramOptions['treeVar'] , "\x1b[0m add sample \x1b[34m", sample,"\x1b[0m from sampleTree \x1b[34m", sampleTree, "\x1b[0m to group \x1b[34m", groupName, "\x1b[0m")
        histogramOptions = self.histogramOptions.copy()
        histogramOptions['group'] = groupName

        self.histoMaker = HistoMaker(self.config, sample=sample, sampleTree=sampleTree, histogramOptions=histogramOptions) 
        sampleHistogram = self.histoMaker.getHistogram(cut)
        self.histograms.append({
            'name': sample.name,
            'histogram': sampleHistogram,
            'group': groupName,
            'signal': sample.type=='SIG'
            })

        if self.config.has_option('Plot_general', 'drawWeightSystematicError'):
            histogramOptions_Up   = histogramOptions.copy()
            histogramOptions_Down = histogramOptions.copy()
            try:
                histogramOptions_Up['weight']   = self.config.get('Weights', self.config.get('Plot_general', 'drawWeightSystematicError') + '_UP')
                histogramOptions_Down['weight'] = self.config.get('Weights', self.config.get('Plot_general', 'drawWeightSystematicError') + '_DOWN')
            except:
                histogramOptions_Up['weight']   = self.config.get('Weights', self.config.get('Plot_general', 'drawWeightSystematicError') + '_Up')
                histogramOptions_Down['weight'] = self.config.get('Weights', self.config.get('Plot_general', 'drawWeightSystematicError') + '_Down')
            histoMaker_Up   = HistoMaker(self.config, sample=sample, sampleTree=sampleTree, histogramOptions=histogramOptions_Up)
            histoMaker_Down = HistoMaker(self.config, sample=sample, sampleTree=sampleTree, histogramOptions=histogramOptions_Down)
            self.histograms[-1]['histogram_Up']   = histoMaker_Up.getHistogram(cut)
            self.histograms[-1]['histogram_Down'] = histoMaker_Down.getHistogram(cut)

            # set up/down normalizations to nominal
            if self.config.has_option('Plot_general', 'drawWeightSystematicErrorNormalized') and eval(self.config.get('Plot_general', 'drawWeightSystematicErrorNormalized')):
                if self.histograms[-1]['histogram_Up'].Integral() > 0:
                    self.histograms[-1]['histogram_Up'].Scale(self.histograms[-1]['histogram'].Integral()/self.histograms[-1]['histogram_Up'].Integral())
                if self.histograms[-1]['histogram_Down'].Integral() > 0:
                    self.histograms[-1]['histogram_Down'].Scale(self.histograms[-1]['histogram'].Integral()/self.histograms[-1]['histogram_Down'].Integral())

    # add object to collection
    def addObject(self, object):
        self.collectedObjects.append(object)

    # ------------------------------------------------------------------------------
    # create canvas and load default style 
    # ------------------------------------------------------------------------------
    def initializeSplitCanvas(self):
        TdrStyles.tdrStyle()

        # initialize canvas
        self.canvas = ROOT.TCanvas(self.var, '', 600, 600)
        self.canvas.SetFillStyle(4000)
        self.canvas.SetFrameFillStyle(1000)
        self.canvas.SetFrameFillColor(0)
        self.canvas.SetTopMargin(0.035)
        self.canvas.SetBottomMargin(0.12)

        self.pads = {}
        if 'fractions' in self.histogramOptions and self.histogramOptions['fractions']:
            self.pads['oben'] = ROOT.TPad('oben', 'oben', 0, 0.4, 1.0, 1.0)
            self.pads['oben'].SetBottomMargin(0)
            self.pads['oben'].SetFillStyle(4000)
            self.pads['oben'].SetFrameFillStyle(1000)
            self.pads['oben'].SetFrameFillColor(0)
            self.pads['fractions'] = ROOT.TPad('fractions', 'fractions', 0, 0.0, 1.0, 0.2)
            self.pads['fractions'].SetTopMargin(0.)
            self.pads['fractions'].SetBottomMargin(0.35)
            self.pads['fractions'].SetFillStyle(4000)
            self.pads['fractions'].SetFrameFillStyle(1000)
            self.pads['fractions'].SetFrameFillColor(0)
            self.pads['unten'] = ROOT.TPad('unten', 'unten', 0, 0.2, 1.0, 0.4)
            self.pads['unten'].SetTopMargin(0.)
            self.pads['unten'].SetBottomMargin(0.0)
            self.pads['unten'].SetFillStyle(4000)
            self.pads['unten'].SetFrameFillStyle(1000)
            self.pads['unten'].SetFrameFillColor(0)

            self.pads['oben'].Draw()
            self.pads['fractions'].Draw()
            self.pads['unten'].Draw()

        else:
            self.pads['oben'] = ROOT.TPad('oben', 'oben', 0, 0.3, 1.0, 1.0)
            self.pads['oben'].SetBottomMargin(0)
            self.pads['oben'].SetFillStyle(4000)
            self.pads['oben'].SetFrameFillStyle(1000)
            self.pads['oben'].SetFrameFillColor(0)
            self.pads['unten'] = ROOT.TPad('unten', 'unten', 0, 0.0, 1.0, 0.3)
            self.pads['unten'].SetTopMargin(0.)
            self.pads['unten'].SetBottomMargin(0.35)
            self.pads['unten'].SetFillStyle(4000)
            self.pads['unten'].SetFrameFillStyle(1000)
            self.pads['unten'].SetFrameFillColor(0)

            self.pads['oben'].Draw()
            self.pads['unten'].Draw()

        self.pads['oben'].cd()
        return self.canvas

    # ------------------------------------------------------------------------------
    # create canvas and load default style 
    # ------------------------------------------------------------------------------
    def initializeCanvas(self):
        TdrStyles.tdrStyle()

        self.canvas = ROOT.TCanvas(self.var+'Comp','',600,600)
        self.canvas.SetFillStyle(4000)
        self.canvas.SetFrameFillStyle(1000)
        self.canvas.SetFrameFillColor(0)
        self.pads = {}
        if 'fractions' in self.histogramOptions and self.histogramOptions['fractions']:
            self.pads['oben'] = ROOT.TPad('oben', 'oben', 0, 0.3, 1.0, 1.0)
            self.pads['oben'].SetBottomMargin(0)
            self.pads['oben'].SetFillStyle(4000)
            self.pads['oben'].SetFrameFillStyle(1000)
            self.pads['oben'].SetFrameFillColor(0)
            self.pads['oben'].SetBottomMargin(0.0)
            self.pads['fractions'] = ROOT.TPad('fractions', 'fractions', 0, 0.0, 1.0, 0.3)
            self.pads['fractions'].SetTopMargin(0.)
            self.pads['fractions'].SetBottomMargin(0.35)
            self.pads['fractions'].SetFillStyle(4000)
            self.pads['fractions'].SetFrameFillStyle(1000)
            self.pads['fractions'].SetFrameFillColor(0)

            self.pads['oben'].Draw()
            self.pads['fractions'].Draw()

        if self.is2D:
            ROOT.gPad.SetTopMargin(0.05)
            ROOT.gPad.SetBottomMargin(0.13)
            ROOT.gPad.SetLeftMargin(0.17)
            ROOT.gPad.SetRightMargin(0.16)
        return self.canvas

    # ------------------------------------------------------------------------------
    # data/MC ratio
    # ------------------------------------------------------------------------------
    def drawRatioPlot(self, dataHistogram, mcHistogram, yAxisTitle="Data/MC", same=False, ratioRange=[0.5,1.75], mcHistogram_Up=None, mcHistogram_Down=None):

        if not self.useSplitCanvas or 'unten' not in self.pads:
            return

        self.pads['unten'].cd()
        ROOT.gPad.SetTicks(1,1)

        self.legends['ratio'] = ROOT.TLegend(0.39, 0.85, 0.93, 0.97)
        self.legends['ratio'].SetLineWidth(2)
        self.legends['ratio'].SetBorderSize(0)
        self.legends['ratio'].SetFillColor(0)
        self.legends['ratio'].SetFillStyle(4000)
        self.legends['ratio'].SetTextSize(0.075)
        self.legends['ratio'].SetNColumns(2)

        # convert TGraphAsymmErrors to TH1D
        if type(dataHistogram) == ROOT.TGraphAsymmErrors:
            print("INFO: converting TGraphAsymmErrors to TH1D for ratio plot...")
            convertedDataHistogram = ROOT.TH1D("data_th1d","data_th1d",self.histogramOptions['nBins'],self.histogramOptions['minX'], self.histogramOptions['maxX'])
            pointX = array.array('d', [0.0, 0.0])
            pointY = array.array('d', [0.0, 0.0])
            for i in range(dataHistogram.GetN()):
                dataHistogram.GetPoint(i, pointX, pointY)
                convertedDataHistogram.SetBinContent(1+i, pointY[0])
                convertedDataHistogram.SetBinError(1+i, 0.5*(dataHistogram.GetErrorYhigh(i)+dataHistogram.GetErrorYlow(i)))
                convertedDataHistogram.SetDrawOption("EP")
            dataHistogram = convertedDataHistogram

        # draw ratio plot
        try:
            self.ratioPlot, error = getRatio(dataHistogram, reference=mcHistogram, min=self.histogramOptions['minX'], max=self.histogramOptions['maxX'], yTitle=yAxisTitle, maxUncertainty=self.maxRatioUncert, restrict=True, yRange=ratioRange)
        except:
             self.ratioPlot = None
             error = None

        if mcHistogram is not None and self.ratioPlot is not None:
            ksScore = dataHistogram.KolmogorovTest(mcHistogram)
            chiScore = dataHistogram.Chi2Test(mcHistogram, "UWCHI2/NDF")
        else:
            ksScore = -1
            chiScore = -1
        print ("INFO: data/MC ratio, KS test:", ksScore, " chi2:", chiScore)
        try:
            self.ratioPlot.SetStats(0)
            self.ratioPlot.GetXaxis().SetTitle(self.xAxis)
            self.ratioError = ROOT.TGraphErrors(error)
            self.ratioError.SetFillColor(ROOT.kGray+3)
            self.ratioError.SetFillStyle(3013)
            self.ratioPlot.Draw("E1 SAME" if same else "E1")
            self.ratioError.Draw('SAME2')
            self.ratioPlot.SetDirectory(0)

            # print table and export ratio histogram with stats+MCstats error combined into .root file
            if self.debug:
                print("INFO: DATA/MC ratios:")
                print("      bin     ratio    error(stat+MCstat)")
                shapesName = "ratioWithStatsPlusMCstatsError_" + self.config.get('Plot_general', 'suffix') + '.root' if self.config.has_option('Plot_general', 'suffix') else "ratioWithStatsPlusMCstatsError.root"
                outputFileName = self.outputFileTemplate.format(outputFolder=self.outputFolder, prefix=self.prefix, prefixSeparator='_' if len(self.prefix)>0 else '', var=self.var,  ext=shapesName)
                print("INFO: \x1b[33mratio with stats+MCstats error\x1b[0m saved to:", outputFileName)
                ratioShapesFile = ROOT.TFile.Open(outputFileName, "RECREATE")
                #thRatioWithTotalError = ROOT.TH1D("rte", "", self.ratioPlot.GetXaxis().GetNbins(), self.ratioPlot.GetXaxis().GetXmin(), self.ratioPlot.GetXaxis().GetXmax())
                thRatioWithTotalError = self.ratioPlot.Clone()
                thRatioWithTotalError.SetDirectory(ratioShapesFile)
                for i in range(self.ratioPlot.GetXaxis().GetNbins()):
                    # mind the different bin number conventions for TH1D and TGraphErrors!
                    statError   = self.ratioPlot.GetBinError(1+i)
                    mcStatError = self.ratioError.GetErrorY(i)
                    totalError = math.sqrt(statError*statError + mcStatError*mcStatError)
                    print("     ", self.ratioPlot.GetXaxis().GetBinCenter(1+i), self.ratioPlot.GetBinContent(1+i), totalError)
                    thRatioWithTotalError.SetBinContent(1+i, self.ratioPlot.GetBinContent(1+i))
                    thRatioWithTotalError.SetBinError(1+i, totalError)
                thRatioWithTotalError.Write()
                ratioShapesFile.Close()

            if self.config.has_option('Plot_general', 'drawWeightSystematicError') and mcHistogram_Up is not None and mcHistogram_Down is not None:
                self.ratioPlot_Up, error_Up     = getRatio(dataHistogram, reference=mcHistogram_Up, min=self.histogramOptions['minX'], max=self.histogramOptions['maxX'], yTitle=yAxisTitle, maxUncertainty=self.maxRatioUncert, restrict=True, yRange=ratioRange)
                self.ratioPlot_Down, error_Down = getRatio(dataHistogram, reference=mcHistogram_Down, min=self.histogramOptions['minX'], max=self.histogramOptions['maxX'], yTitle=yAxisTitle, maxUncertainty=self.maxRatioUncert, restrict=True, yRange=ratioRange)

                px = []
                py = []
                ex = []
                eyh = []
                eyl = []

                x = array.array('d', [0.0])
                y = array.array('d', [0.0])

                for i in range(self.ratioError.GetN()):
                    self.ratioError.GetPoint(i, x, y)
                    px.append(x[0])
                    py.append(self.ratioPlot.GetBinContent(1+i))
                    ex.append(self.ratioError.GetErrorX(i))

                    var_up   = self.ratioPlot_Up.GetBinContent(1+i)
                    var_down = self.ratioPlot_Down.GetBinContent(1+i)

                    var_max  = max(var_up, var_down)
                    var_min  = min(var_up, var_down)

                    eyh.append(abs(var_max-self.ratioPlot.GetBinContent(1+i)))
                    eyl.append(abs(var_min-self.ratioPlot.GetBinContent(1+i)))

                a_px  = array.array('d', px)
                a_py  = array.array('d', py)
                a_ex  = array.array('d', ex)
                a_eyh = array.array('d', eyh)
                a_eyl = array.array('d', eyl)

                if self.config.has_option('Plot_general', 'drawWeightSystematicErrorHatchesSpacing'):
                    ROOT.gStyle.SetHatchesSpacing(float(self.config.get('Plot_general', 'drawWeightSystematicErrorHatchesSpacing')))
                self.ratioWeightSystematicError = ROOT.TGraphAsymmErrors(len(px), a_px, a_py, a_ex, a_ex, a_eyl, a_eyh)
                self.ratioWeightSystematicError.SetFillColor(eval(self.config.get('Plot_general', 'drawWeightSystematicErrorColor')) if self.config.has_option('Plot_general', 'drawWeightSystematicErrorColor') else ROOT.kBlue)
                self.ratioWeightSystematicError.SetFillStyle(3144)
                self.ratioWeightSystematicError.Draw('SAME2')
                self.ratioPlot.Draw("E1 SAME")
                self.legends['ratio'].AddEntry(self.ratioWeightSystematicError, self.config.get('Plot_general', 'drawWeightSystematicError') ,"f")

        except Exception as e:
            print ("\x1b[31mERROR: with ratio histogram!", e, "\x1b[0m")

        # draw blinded region in the ratio plot with hashed area
        if 'blindCut' in self.histogramOptions:

            print("DEBUG: blinding cut will be applied:", self.histogramOptions['blindCut'], self.histogramOptions['treeVar'], self.histogramOptions['blindCut'].startswith(self.histogramOptions['treeVar']))
            blindCutExpr = self.histogramOptions['blindCut'].replace(self.histogramOptions['treeVar'], '{var}')
            print("DEBUG: final expression is:", blindCutExpr)
            if self.ratioPlot is not None:

                nBinsOriginal      = dataHistogram.GetXaxis().GetNbins()
                binListBlindRegion = array.array('d', [dataHistogram.GetXaxis().GetBinLowEdge(1+i) for i in range(nBinsOriginal)] + [dataHistogram.GetXaxis().GetXmax()])
                blindedRegion      = ROOT.TH1D("blind", "blind", nBinsOriginal, binListBlindRegion) 

                for i in range(nBinsOriginal):

                    if len(self.histoMaker.originalBins) > 0:
                        # if plot is drawn with 'equal size' option, then have to look at original bin boundaries 
                        binLowEdgeValue = self.histoMaker.originalBins[i][0] + 1e-6
                        binUpEdgeValue = self.histoMaker.originalBins[i][1] - 1e-6
                    else:
                        # otherwise, take bins directly from histogram
                        binLowEdgeValue = dataHistogram.GetXaxis().GetBinLowEdge(1+i) + 1e-6
                        binUpEdgeValue = dataHistogram.GetXaxis().GetBinUpEdge(1+i) - 1e-6
                    try:
                        exprLow = ROOT.gInterpreter.ProcessLine(blindCutExpr.replace('{var}',str(binLowEdgeValue)))
                        exprUp  = ROOT.gInterpreter.ProcessLine(blindCutExpr.replace('{var}',str(binUpEdgeValue))) 
                    except Exception as e:
                        print("EXCEPTION: could not draw blind region in ratio plot: ", e)
                        exprLow = True
                        exprUp = True
                    value = 1.0
                    if not exprLow or not exprUp: 
                        # draw hashed area as an error bar
                        error = 2.0
                    else:
                        error = 0.0
                    blindedRegion.SetBinContent(1+i, value)
                    blindedRegion.SetBinError(1+i, error)
                    if self.debug:
                        print("DEBUG:", i, binLowEdgeValue, binUpEdgeValue, error) 
                blindedRegion.SetFillColor(ROOT.kRed)
                blindedRegion.SetFillStyle(3018)
                blindedRegion.SetMarkerSize(0)
                blindedRegion.Draw("SAME E2")
                self.addObject(blindedRegion)

        if self.ratioPlot is not None:
            self.m_one_line = ROOT.TLine(self.histogramOptions['minX'], 1, self.histogramOptions['maxX'], 1)
            self.m_one_line.SetLineStyle(ROOT.kSolid)
            self.m_one_line.Draw("Same")

            self.legends['ratio'].AddEntry(self.ratioError, self.mcUncertaintyLegend,"f")
            self.legends['ratio'].Draw() 
            if not self.blind:
                self.addObject(self.myText("#chi^{2}_{ }#lower[0.1]{/^{}#it{dof} = %.2f}"%(chiScore), self.plotTextMarginLeft, 0.895, 1.55))
                t0 = ROOT.TText()
                t0.SetTextSize(ROOT.gStyle.GetLabelSize()*2.4)
                t0.SetTextFont(ROOT.gStyle.GetLabelFont())
                if not self.log:
                    t0.DrawTextNDC(0.1059, 0.96, "0")

    def drawSampleLegend(self, groupedHistograms, theErrorGraph, normalize=False):
        if 'oben' in self.pads and self.pads['oben']:
            self.pads['oben'].cd()

        legendX1 = float(self.config.get('Plot_general','legendX1')) if self.config.has_option('Plot_general','legendX1') else 0.45
        legendY1 = float(self.config.get('Plot_general','legendY1')) if self.config.has_option('Plot_general','legendY1') else 0.6
        legendX2 = float(self.config.get('Plot_general','legendX2')) if self.config.has_option('Plot_general','legendX2') else 0.92
        legendY2 = float(self.config.get('Plot_general','legendY2')) if self.config.has_option('Plot_general','legendY2') else 0.92

        self.legends['left'] = ROOT.TLegend(legendX1,legendY1,legendX1+0.638*(legendX2-legendX1),legendY2)
        self.legends['left'].SetLineWidth(2)
        self.legends['left'].SetBorderSize(0)
        self.legends['left'].SetFillColor(0)
        self.legends['left'].SetFillStyle(4000)
        self.legends['left'].SetTextFont(62)
        self.legends['left'].SetTextSize(0.035)
        self.legends['right'] = ROOT.TLegend(legendX1+0.49*(legendX2-legendX1), legendY1,legendX2,legendY2)
        self.legends['right'].SetLineWidth(2)
        self.legends['right'].SetBorderSize(0)
        self.legends['right'].SetFillColor(0)
        self.legends['right'].SetFillStyle(4000)
        self.legends['right'].SetTextFont(62)
        self.legends['right'].SetTextSize(0.035)


        nLeft = 0
        nRight = 0
        
        legendEntries = []
        if self.dataGroupName in groupedHistograms:
            self.legends['left'].AddEntry(groupedHistograms[self.dataGroupName], self.dataTitle, 'P')
            nLeft += 1
        groupNames = list(set([groupName for groupName, groupHistogram in groupedHistograms.iteritems()]))
        groupNamesOrdered = self.setup + sorted([x for x in groupNames if x not in self.setup])

        numLegendEntries = len(groupNames) + 2
        if self.config.has_option('Plot_general', '__modNumLegentries'):
            numLegendEntries += int(self.config.get('Plot_general', '__modNumLegentries'))

        for itemPosition, groupName in enumerate(groupNamesOrdered): 
            if groupName != self.dataGroupName and groupName in groupedHistograms:
                legendEntryName = self.typLegendDict[groupName] if groupName in self.typLegendDict else groupName
                if itemPosition < numLegendEntries/2.-2:
                    self.legends['left'].AddEntry(groupedHistograms[groupName], legendEntryName, 'F')
                    nLeft += 1
                else:
                    self.legends['right'].AddEntry(groupedHistograms[groupName], legendEntryName, 'F')
                    nRight += 1
            elif groupName not in groupedHistograms:
                print("WARNING: histogram group not found:", groupName)
        if theErrorGraph and not normalize:
            self.legends['right'].AddEntry(theErrorGraph, self.mcUncertaintyLegend, "fl")
            nRight += 1
        self.canvas.Update()

        for legendEntry in self.legendEntries:
            if nLeft < nRight:
                pos = 'left'
                nLeft += 1
            else:
                pos = 'right'
                nRight += 1
            if len(legendEntry) > 2:
                self.legends[pos].AddEntry(legendEntry[0],legendEntry[1],legendEntry[2])
            else:
                self.legends[pos].AddEntry(legendEntry[0],legendEntry[1])


        ROOT.gPad.SetTicks(1,1)
        self.legends['left'].SetFillColor(0)
        self.legends['left'].SetBorderSize(0)
        self.legends['right'].SetFillColor(0)
        self.legends['right'].SetBorderSize(0)
        self.legends['left'].Draw()
        self.legends['right'].Draw()

    def drawPlotTexts(self):
        if 'oben' in self.pads and self.pads['oben']:
            self.pads['oben'].cd()
        posY = 0.88
        if type(self.plotTitle) == list:
            size = 1.04
            for plotTitleLine in self.plotTitle:
                self.addObject(self.myText(plotTitleLine, self.plotTextMarginLeft + (0.03 if self.is2D else 0),posY,size))
                posY -= 0.05
                size *= 0.77
        else:
            self.addObject(self.myText(self.plotTitle,self.plotTextMarginLeft+(0.03 if self.is2D else 0),posY,1.04))
            posY -= 0.05
        print ('self.lumi is', self.lumi)
        try:
            self.addObject(self.myText("#sqrt{s} = %s, L = %.2f fb^{-1}"%(self.anaTag, (float(self.lumi)/1000.0)), self.plotTextMarginLeft+(0.03 if self.is2D else 0), posY, 0.75))
        except Exception as e:
            print ("WARNING: exception while adding text: ", e)
            pass
        dataNames = list(set([histogram['name'] for histogram in self.histograms if histogram['group'] == self.dataGroupName]))
        addFlag = ''
        isZee = False
        isZmm = False
        for data_ in dataNames:
            if 'DoubleEG' in data_:
                isZee = True
            if 'DoubleMuon' in data_:
                isZmm = True
        if ('Zee' in dataNames and 'Zmm' in dataNames) or (isZee and isZmm):
            addFlag = 'Z(l^{-}l^{+})H(b#bar{b})'
        elif 'Zee' in dataNames or isZee:
            addFlag = 'Z(e^{-}e^{+})H(b#bar{b})'
        elif 'Zmm' in dataNames or isZmm:
            addFlag = 'Z(#mu^{-}#mu^{+})H(b#bar{b})'
        elif 'Znn' in dataNames:
            addFlag = 'Z(#nu#nu)H(b#bar{b})'
        elif 'Wmn' in dataNames:
            addFlag = 'W(#mu#nu)H(b#bar{b})'
        elif 'Wen' in dataNames:
            addFlag = 'W(e#nu)H(b#bar{b})'
        #self.addObject(self.myText(addFlag, self.plotTextMarginLeft+(0.03 if self.is2D else 0), 0.78))

        try:
            for labelName, label in self.plotLabels.iteritems():
                self.addObject(self.myText(label['text'], label['x'], label['y'], label['size']))
        except:
            pass

        try:
            for j, additionalTextLine in enumerate(self.additionalTextLines):
                self.addObject(self.myText(additionalTextLine, self.plotTextMarginLeft, 0.73-0.03*j, 0.6))
        except Exception as e:
            print(e)

        #print 'Add Flag %s' %self.addFlag2
        #if self.addFlag2:
        #    tAddFlag2 = self.myText(self.addFlag2,0.17,0.73)

    # adds TH1 objects together
    @staticmethod
    def sumHistograms(histograms, outputName='summed_histogram'):
        summedHistogram = None
        for histogram in histograms:
            if summedHistogram:
                summedHistogram.Add(histogram)
            else:
                summedHistogram = histogram.Clone(outputName)
        return summedHistogram

    # ------------------------------------------------------------------------------
    # draw the stacked histograms 
    # dataOverBackground=False => draw DATA/MC
    # ------------------------------------------------------------------------------
    def Draw(self, outputFolder='./', prefix='', normalize=False, dataOverBackground=False, ratioRange=None):

        if ratioRange is None:
            if self.config.has_option('Plot_general', 'ratioRange'):
                ratioRange = eval(self.config.get('Plot_general', 'ratioRange'))
                print("INFO: using custom range for ratio plot:", ratioRange)
            else:
                ratioRange = [0.5,1.75]

        self.is2D = any([isinstance(h['histogram'], ROOT.TH2) for h in self.histograms])
        self.outputFolder = outputFolder
        self.prefix = prefix

        # MC histograms, defined in setup
        dataGroupName = self.dataGroupName
        mcHistogramGroups          = list(set([histogram['group'] for histogram in self.histograms if histogram['group']!=dataGroupName]))
        mcHistogramGroupsToPlot    = sorted(mcHistogramGroups, key=lambda x: self.setup.index(x) if x in self.setup else 9999)
        mcHistogramGroupsUndefined = [x for x in mcHistogramGroups if x not in self.setup]
        if len(mcHistogramGroupsUndefined) > 0:
            print("\x1b[97m\x1b[41mWARNING: some MC samples are not defined in 'setup' definition for plots: \x1b[0m")
            for hiddenGroup in mcHistogramGroupsUndefined:
                print(" > ", hiddenGroup, " is not defined in setup")
        if dataGroupName in mcHistogramGroupsToPlot:
            raise Exception("DATA contained in MC groups!")


        # rebinning
        rebinMethod = self.config.get(self.plotVarSection, 'rebinMethod').lower() if self.config.has_option(self.plotVarSection, 'rebinMethod') else 'none'
        print("DEBUG: rebin method:", rebinMethod)
        if rebinMethod == 'flatsignal':
            signalHistograms = [histogram['histogram'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot and histogram['signal']]
            signalSum = NewStackMaker.sumHistograms(histograms=signalHistograms, outputName="signalHistogramsSum")
            signalIntegral = signalSum.Integral()

            cumulatedSignal = 0.0
            nBinsTarget = 15
            newBinBoundaries = [float(self.histogramOptions['min'])]
            for i in range(self.histogramOptions['nBins']):
                cumulatedSignal += signalSum.GetBinContent(1+i)
                if cumulatedSignal >= (len(newBinBoundaries))*(signalIntegral/nBinsTarget):
                    newBinBoundaries.append(signalSum.GetXaxis().GetBinUpEdge(1+i))
            if len(newBinBoundaries) < nBinsTarget + 1:
                newBinBoundaries.append(float(self.histogramOptions['max']))
            else:
                newBinBoundaries[nBinsTarget] = float(self.histogramOptions['max'])
            print("INFO: new bin boundaries have been computed, to be put to [plotDef:...] binList=")
            print("INFO: ", "[" + ", ".join(["%1.5f"%x for x in newBinBoundaries]) + "]") 
            print("INFO: plotting is skipped for this variable, please change binList in the config and rerun")
            print("INFO: to ensure correct MC errors in the histograms.")
            return
        elif rebinMethod != 'none':
            if rebinMethod == 'arctansignal':
                # Round[Table[(Pi/2 + ArcTan[x])/Pi, {x, -3, 11}], 0.0001]
                fractions = np.array([0.1024, 0.1476, 0.25, 0.5, 0.75, 0.8524, 0.8976, 0.922, 0.9372, 0.9474, 0.9548, 0.9604, 0.9648, 0.9683, 0.9711])
            elif rebinMethod == 'gausssignal':
                # Table[Exp[-(x - 0.7)^2/0.3^2], {x, 0, 0.99, 1/15}]
                fractions = np.array([0.0043, 0.0116, 0.0282, 0.0622, 0.1241, 0.2245, 0.3679, 0.5461, 0.7344, 0.8948, 0.9877, 0.9877, 0.8948, 0.7344, 0.5461])
            elif rebinMethod.strip().startswith('['):
                fractions = np.array(eval(rebinMethod))
            else:
                print("ERROR: not valid:", rebinMethod)
                raise Exception("RebinMethodError")
            cumulativeFractions = np.cumsum(fractions/np.sum(fractions))
            print("DEBUG: c=", cumulativeFractions)

            signalHistograms = [histogram['histogram'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot and histogram['signal']]
            signalSum = NewStackMaker.sumHistograms(histograms=signalHistograms, outputName="signalHistogramsSum")
            signalIntegral = signalSum.Integral()

            cumulatedSignal = 0.0
            nBinsTarget = 15
            newBinBoundaries = [float(self.histogramOptions['min'])]
            for i in range(self.histogramOptions['nBins']):
                cumulatedSignal += signalSum.GetBinContent(1+i)/signalIntegral
                if len(newBinBoundaries) < nBinsTarget + 1 and cumulatedSignal >= cumulativeFractions[len(newBinBoundaries)-1]: 
                    newBinBoundaries.append(signalSum.GetXaxis().GetBinUpEdge(1+i))
            if len(newBinBoundaries) < nBinsTarget + 1:
                newBinBoundaries.append(float(self.histogramOptions['max']))
            else:
                newBinBoundaries[nBinsTarget] = float(self.histogramOptions['max'])
            print("INFO: new bin boundaries have been computed, to be put to [plotDef:...] binList=")
            print("INFO: ", "[" + ", ".join(["%1.5f"%x for x in newBinBoundaries]) + "]") 
            print("INFO: plotting is skipped for this variable, please change binList in the config and rerun")
            print("INFO: to ensure correct MC errors in the histograms.")
            return


        # group ("sum") MC+DATA histograms 
        groupedHistograms = {}
        histogramGroups = list(set([histogram['group'] for histogram in self.histograms]))
        for histogramGroup in histogramGroups:
            histogramsInGroup = [histogram['histogram'] for histogram in self.histograms if histogram['group'] == histogramGroup]
            groupedHistograms[histogramGroup] = NewStackMaker.sumHistograms(histograms=histogramsInGroup, outputName="group_" + histogramGroup)
            try:
                groupedHistograms[histogramGroup].SetStats(0)
            except:
                pass

        # [histogram['histogram'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot]

        mcHistogramList = [histogram['histogram'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot]

        self.useSplitCanvas = not (normalize or self.is2D or dataGroupName not in histogramGroups or len(mcHistogramList) < 1)
        c = self.initializeSplitCanvas() if self.useSplitCanvas else self.initializeCanvas()

        # add summed MC histograms to stack
        allStack = ROOT.THStack(self.var, '')
        colorDict = eval(self.config.get('Plot_general', 'colorDict'))
        maximumNormalized = 0
        for groupName in mcHistogramGroupsToPlot[::-1]:
            if groupName in groupedHistograms and groupName != dataGroupName:
                if groupName in colorDict:
                    if normalize:
                        groupedHistograms[groupName].SetFillColor(0)
                        groupedHistograms[groupName].SetLineColor(colorDict[groupName])
                        groupedHistograms[groupName].SetLineWidth(3)
                    else:
                        groupedHistograms[groupName].SetFillColor(colorDict[groupName])
                if groupedHistograms[groupName]:
                    groupedHistograms[groupName].SetStats(0)
                    if normalize and groupedHistograms[groupName].Integral()>0:
                        groupedHistograms[groupName].Scale(1./groupedHistograms[groupName].Integral())
                        if groupedHistograms[groupName].GetMaximum() > maximumNormalized:
                            maximumNormalized = groupedHistograms[groupName].GetMaximum()
                    if not normalize:
                        groupedHistograms[groupName].SetLineColor(ROOT.kBlack)
                    allStack.Add(groupedHistograms[groupName])

        # draw options
        drawOption = "hist" if not normalize else "histnostack"
        if self.is2D:
            drawOption = self.histogramOptions['drawOption'] if 'drawOption' in self.histogramOptions else 'colz'

        # draw stack/sum
        if self.is2D:
            self.plotLabels['topright1'] = {'text': '2Dplot', 'x': 0.69, 'y': 0.9, 'size': 0.7}
            if 'draw' in self.histogramOptions and self.histogramOptions['draw'].strip() == self.dataGroupName:
                allStack = NewStackMaker.sumHistograms(histograms=[histogram['histogram'] for histogram in self.histograms if histogram['group'] == self.dataGroupName], outputName='summedDataHistograms')
                self.plotLabels['topright1']['text'] = 'DATA'
            elif 'draw' in self.histogramOptions and self.histogramOptions['draw'].strip().upper() == 'RATIO':
                allStackData = NewStackMaker.sumHistograms(histograms=[histogram['histogram'] for histogram in self.histograms if histogram['group'] == self.dataGroupName], outputName='summedDataHistograms')
                allStackMc = NewStackMaker.sumHistograms(histograms=[histogram['histogram'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot], outputName='summedMcHistograms')
                allStack = allStackData.Clone('summedRatioHistogram')
                allStack.Sumw2()
                allStack.Divide(allStackMc)
                self.plotLabels['topright1']['text'] = 'DATA/MC'
                try:
                    unityPos = 0.33
                    try:
                        if 'minZ' in self.histogramOptions and 'maxZ' in self.histogramOptions:
                            unityPos = (1.0-self.histogramOptions['minZ'])/(self.histogramOptions['maxZ']-self.histogramOptions['minZ'])
                    except:
                        pass
                    stops = array.array('d', [0.0, unityPos, 1.0])
                    reds = array.array('d', [0.34, 0.63, 1.0])
                    greens = array.array('d', [0.34, 0.63, 0])
                    blues = array.array('d', [0.63, 0.63, 0])
                    nc = 20
                    ROOT.TColor.CreateGradientColorTable(len(stops), stops, reds, greens, blues, nc)
                    ROOT.gStyle.SetNumberContours(nc)
                except Exception as e:
                    print("\x1b[31m",e,"\x1b[0m")
            else:
                allStack = NewStackMaker.sumHistograms(histograms=[histogram['histogram'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot], outputName='summedMcHistograms')
                self.plotLabels['topright1']['text'] = 'MC'
            allStack.SetStats(0)
            allStack.SetTitle('')
        ROOT.gPad.SetLogy(0)

        if 'fractions' in self.histogramOptions and self.histogramOptions['fractions']:
            self.pads['fractions'].cd()
            allStackCopy = allStack.Clone() 
            histList = allStackCopy.GetHists()    
            for binNumber in range(self.histogramOptions['nBins']):
                binTotal = sum([h.GetBinContent(1+binNumber) for h in histList])
                if binTotal > 0:
                    scale = 100.0/binTotal
                    print("DEBUG: scale by", scale, " in bin ", binNumber)
                    for h in histList:
                        h.SetBinContent(1+binNumber,h.GetBinContent(1+binNumber) * scale)
                        h.SetBinError(1+binNumber,h.GetBinError(1+binNumber) * scale)
            ROOT.gStyle.SetHistTopMargin(0.)
            allStackCopy.SetMaximum(100.0)
            allStackCopy.Draw(drawOption)
            allStackCopy.GetYaxis().SetTitle("fraction [%]")
            allStackCopy.GetXaxis().SetTitle(self.xAxis)
            allStackCopy.GetXaxis().SetTitleSize(0.16)
            allStackCopy.GetXaxis().SetLabelSize(0.16)
            allStackCopy.GetYaxis().SetTitleSize(0.16)
            allStackCopy.GetYaxis().SetTitleOffset(0.3)
            allStackCopy.GetYaxis().SetLabelSize(0.16)
            allStackCopy.SetMaximum(100.001)
            self.garbage.append(allStackCopy)
            self.garbage.append(histList)

            self.pads['oben'].cd()

        allStack.Draw(drawOption)
        self.garbage.append(allStack)

        if self.config.has_option('Plot_general', 'drawSignal') and eval(self.config.get('Plot_general', 'drawSignal')):
            factor = 1.0*eval(self.config.get('Plot_general', 'drawSignal'))
            print("DEBUG: draw signal line, multiplied by:", factor)
            signalHistograms = [histogram['histogram'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot and histogram['signal']]
            signalSum = NewStackMaker.sumHistograms(histograms=signalHistograms, outputName="signalHistogramsSum")
            signalSum.SetLineColor(eval(self.config.get('Plot_general', 'drawSignalColor')) if self.config.has_option('Plot_general', 'drawSignalColor') else ROOT.kRed)
            signalSum.SetFillStyle(0)
            signalSum.Scale(factor)
            signalSum.Draw("hist same")

            if factor != 1.0:
                self.legendEntries.append([signalSum, "signal (%1.1f x)"%factor])
            else:
                self.legendEntries.append([signalSum, "signal"])

        # set axis titles
        if self.is2D:
            yTitle = self.yAxis
        else:
            try:
                if dataGroupName in groupedHistograms and not groupedHistograms[dataGroupName].GetSumOfWeights() % 1 == 0.0:
                    yTitle = 'S/(S+B) weighted entries'
                else:
                    yTitle = 'Entries'
                if '/' not in yTitle:
                    if allStack and allStack.GetXaxis():
                        if 'GeV' in self.xAxis:
                            yAppend = '%.0f' %(allStack.GetXaxis().GetBinWidth(1))
                        else:
                            yAppend = '%.2f' %(allStack.GetXaxis().GetBinWidth(1))
                        yTitle = '%s / %s' %(yTitle, yAppend)
                        if 'GeV' in self.xAxis:
                            yTitle += ' GeV'
            except:
                yTitle = "Entries"
        self.yTitle = yTitle
        if allStack and allStack.GetXaxis():
            allStack.GetYaxis().SetTitle(self.yTitle if self.yTitle else '-')
            allStack.GetXaxis().SetRangeUser(self.histogramOptions['minX'], self.histogramOptions['maxX'])
            if not self.is2D:
                allStack.GetYaxis().SetRangeUser(0,90000)
        if not self.is2D:
            if normalize:
                Ymax = maximumNormalized * 1.5
            else:
                Ymax = max(allStack.GetMaximum(), groupedHistograms[dataGroupName].GetMaximum() if dataGroupName in groupedHistograms else 0) * 1.7
            if self.log and not normalize:
                allStack.SetMinimum(0.1)
                Ymax = Ymax*ROOT.TMath.Power(10,1.2*(ROOT.TMath.Log(1.2*(Ymax/0.2))/ROOT.TMath.Log(10)))*(0.2*0.1)
                ROOT.gPad.SetLogy()
            else:
                ROOT.gPad.SetLogy(0)
            allStack.SetMaximum(Ymax)
        if (not normalize and dataGroupName in groupedHistograms) and not self.is2D:
            if allStack and allStack.GetXaxis():
                allStack.GetXaxis().SetLabelOffset(999)
                allStack.GetXaxis().SetLabelSize(0)
        else:
            allStack.GetXaxis().SetTitle(self.xAxis)


        if self.is2D:
            if 'minZ' in self.histogramOptions and 'maxZ' in self.histogramOptions:
                allStack.GetZaxis().SetRangeUser(self.histogramOptions['minZ'], self.histogramOptions['maxZ'])

        # draw DATA
        if dataGroupName in groupedHistograms and not self.is2D:
            drawOption = self.histogramOptions['drawOptionData'] if 'drawOptionData' in self.histogramOptions else 'PE'
            if allStack and allStack.GetXaxis():
                drawOption += ',SAME'

            if self.asimovData:
                dataHistogram = NewStackMaker.sumHistograms(histograms=[histogram['histogram'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot], outputName='asimovDataHistogram')
                for i in range(dataHistogram.GetXaxis().GetNbins()):
                    dataHistogram.SetBinError(1+i, math.sqrt(dataHistogram.GetBinContent(1+i)))
            else:
                dataHistogram = groupedHistograms[dataGroupName]
            if normalize and dataHistogram.Integral() > 0:
                dataHistogram.Scale(1./dataHistogram.Integral())
            dataHistogram.Draw(drawOption)

            # adjust axis range to data if MC not present and put axis labels to data histogram
            if len(mcHistogramList) < 1:
                dataHistogram.GetYaxis().SetRangeUser(0, dataHistogram.GetBinContent(dataHistogram.GetMaximumBin()) * 1.5)
                dataHistogram.GetXaxis().SetTitle(self.xAxis)
                dataHistogram.GetYaxis().SetTitle(self.yTitle)

        # draw total entry number
        if not self.is2D and 'drawOption' in self.histogramOptions and 'TEXT' in self.histogramOptions['drawOption'].upper():
            mcHistogram0 = NewStackMaker.sumHistograms(histograms=[histogram['histogram'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot], outputName='summedMcHistograms0')
            try:
                mcHistogram0.SetLineColor(ROOT.kBlack)
                mcHistogram0.SetMarkerColor(ROOT.kBlack)
                mcHistogram0.SetMarkerStyle(1)
                #mcHistogram0.SetMarkerSize(0.001)
                mcHistogram0.SetStats(0)
                mcHistogram0.Draw("SAME TEXT0")
                print("drawn total entry histogram!!!")
                binContents = ["%1.4f"%mcHistogram0.GetBinContent(1+j) for j in range(mcHistogram0.GetXaxis().GetNbins())]
                print("bin contents:", ", ".join(binContents))
            except:
                pass

        # draw ratio plot
        theErrorGraph = None
        mcHistogram   = None
        if not normalize and not self.is2D:
            if dataGroupName in groupedHistograms:

                if self.asimovData:
                    dataHistogram = NewStackMaker.sumHistograms(histograms=[histogram['histogram'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot], outputName='asimovDataHistogram')
                    for i in range(dataHistogram.GetXaxis().GetNbins()):
                        dataHistogram.SetBinError(1+i, math.sqrt(dataHistogram.GetBinContent(1+i)))
                else:
                    dataHistogram = groupedHistograms[dataGroupName]

                if dataOverBackground:
                    bkgHistogramList = [histogram['histogram'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot and not histogram['signal']]
                    sigHistogramList = [histogram['histogram'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot and histogram['signal']]
                    backgroundHistogram = NewStackMaker.sumHistograms(histograms=bkgHistogramList, outputName='summedBackgroundHistograms')
                    signalHistogram     = NewStackMaker.sumHistograms(histograms=sigHistogramList, outputName='summedSignalHistograms')
                    mcHistogram         = backgroundHistogram.Clone()
                    mcHistogram.Add(signalHistogram)

                    self.drawRatioPlot(dataHistogram, backgroundHistogram, yAxisTitle='Data/BKG', ratioRange=ratioRange)

                    # draw blue line for (S+B)/B expectation
                    expectedRatio = mcHistogram.Clone()
                    expectedRatio.Divide(backgroundHistogram)
                    expectedRatio.SetLineColor(ROOT.kBlue)
                    expectedRatio.SetLineWidth(2)
                    expectedRatio.SetFillStyle(0)
                    self.pads['unten'].cd()
                    expectedRatio.Draw("HIST;SAME")

                    # set up/down normalizations of groups to nominal
                    if self.config.has_option('Plot_general', 'drawWeightSystematicErrorGroupsNormalized') and eval(self.config.get('Plot_general', 'drawWeightSystematicErrorGroupsNormalized')):
                        print("INFO: normalize systematic variations per group")

                        groupNormalizations = {'Nom': {}, 'Up':{}, 'Down': {}}
                        for histogram in self.histograms:
                            if histogram['group'] in mcHistogramGroupsToPlot:
                                if histogram['group'] not in groupNormalizations['Nom']:
                                    groupNormalizations['Nom'][histogram['group']] = 0.0
                                    groupNormalizations['Up'][histogram['group']] = 0.0
                                    groupNormalizations['Down'][histogram['group']] = 0.0
                                groupNormalizations['Nom'][histogram['group']] += histogram['histogram'].Integral()
                                groupNormalizations['Up'][histogram['group']] += histogram['histogram_Up'].Integral()
                                groupNormalizations['Down'][histogram['group']] += histogram['histogram_Down'].Integral()

                        for histogram in  self.histograms:
                            if histogram['group'] in mcHistogramGroupsToPlot:
                                if groupNormalizations['Up'][histogram['group']] > 0:
                                    histogram['histogram_Up'].Scale(groupNormalizations['Nom'][histogram['group']]/groupNormalizations['Up'][histogram['group']])
                                if groupNormalizations['Down'][histogram['group']] > 0:
                                    histogram['histogram_Down'].Scale(groupNormalizations['Nom'][histogram['group']]/groupNormalizations['Down'][histogram['group']])

                    if self.config.has_option('Plot_general', 'drawWeightSystematicError'):
                        backgroundHistogram_Up   = NewStackMaker.sumHistograms(histograms=[histogram['histogram_Up'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot and not histogram['signal']], outputName='summedBackgroundHistograms_Up')
                        backgroundHistogram_Down = NewStackMaker.sumHistograms(histograms=[histogram['histogram_Down'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot and not histogram['signal']], outputName='summedBackgroundHistograms_Down')
                        self.drawRatioPlot(dataHistogram, backgroundHistogram, mcHistogram_Up=backgroundHistogram_Up, mcHistogram_Down=backgroundHistogram_Down, yAxisTitle='Data/BKG', same=True, ratioRange=ratioRange)
                    else:
                        self.drawRatioPlot(dataHistogram, backgroundHistogram, yAxisTitle='Data/BKG', same=True, ratioRange=ratioRange)

                else:
                    mcHistogram = NewStackMaker.sumHistograms(histograms=[histogram['histogram'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot], outputName='summedMcHistograms')
                    if self.config.has_option('Plot_general', 'drawWeightSystematicError'):
                        mcHistogram_Up   = NewStackMaker.sumHistograms(histograms=[histogram['histogram_Up'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot], outputName='summedMcHistograms_Up')
                        mcHistogram_Down = NewStackMaker.sumHistograms(histograms=[histogram['histogram_Down'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot], outputName='summedMcHistograms_Down')
                        self.drawRatioPlot(dataHistogram, mcHistogram, mcHistogram_Up=mcHistogram_Up, mcHistogram_Down=mcHistogram_Down, ratioRange=ratioRange)
                    else:
                        self.drawRatioPlot(dataHistogram, mcHistogram, ratioRange=ratioRange)

                if 'oben' in self.pads:
                    self.pads['oben'].cd()
            else:
                print("INFO: no DATA available")

            # draw MC error
            if len(mcHistogramList) > 0:
                mcHistogram = NewStackMaker.sumHistograms(histograms=mcHistogramList, outputName='summedMcHistograms')
                theErrorGraph = ROOT.TGraphErrors(mcHistogram)
                theErrorGraph.SetFillColor(ROOT.kGray+3)
                theErrorGraph.SetFillStyle(3013)
                theErrorGraph.Draw('SAME2')
                #    self.garbage.append(mcHistogramClone)
                #    self.garbage.append(theErrorGraph)

        elif normalize and not self.is2D and self.drawMCErrorForNormalizedPlots:
            # draw error bands for all individual groups
            self.errorGraphs = []
            colorDict = eval(self.config.get('Plot_general', 'colorDict'))
            histogramGroups = set(list([histogram['group'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot]))
            print(groupedHistograms)
            for histogramGroup in histogramGroups:
                theErrorGraph = ROOT.TGraphErrors(NewStackMaker.sumHistograms(histograms=[groupedHistograms[histogramGroup]], outputName='summedMcHistograms_'+histogramGroup))
                theErrorGraph.SetFillColor(colorDict[histogramGroup])
                theErrorGraph.SetFillStyle(3013)
                theErrorGraph.Draw('SAME2')
                self.errorGraphs.append(theErrorGraph)

        # draw legend
        if not self.is2D:
            self.drawSampleLegend(groupedHistograms, theErrorGraph, normalize=normalize)

        # draw various labels
        if not normalize:
            self.drawPlotTexts()

        # save to file
        for ext in self.outputFileFormats:
            print("INFO: save as ", ext)
            sys.stdout.flush()
            outputFileName = self.outputFileTemplate.format(outputFolder=outputFolder, prefix=prefix, prefixSeparator='_' if len(prefix)>0 else '', var=self.var, ext=ext)
            c.SaveAs(outputFileName)
            if os.path.isfile(outputFileName):
                print("INFO: saved as \x1b[34m", outputFileName, "\x1b[0m")
            else:
                print("\x1b[31mERROR: could not save canvas to the file:", outputFileName, "\x1b[0m")
        if self.outputTeX:
            outputFileName = self.outputFileTemplate.format(outputFolder=outputFolder, prefix=prefix, prefixSeparator='_' if len(prefix)>0 else '', var=self.var, ext='tex')
            #ROOT.gPad.Print(outputFileName)
            c.SaveAs(outputFileName)
            print("INFO: saved as \x1b[32mTeX \x1b[34m", outputFileName, "\x1b[0m")

        self.histoCounts = {'unweighted':{}, 'weighted': {}}

        # save shapes
        if not normalize and self.saveShapes:
            mcHistogram = NewStackMaker.sumHistograms(histograms=[histogram['histogram'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot], outputName='summedMcHistograms')
            print("MC:", mcHistogram)
            try:
                shapesName = "shapes.root"
                if self.config.has_option('Plot_general', 'suffix'):
                    shapesName = "shapes_" + self.config.get('Plot_general', 'suffix') + '.root'
                outputFileName = self.outputFileTemplate.format(outputFolder=outputFolder, prefix=prefix, prefixSeparator='_' if len(prefix)>0 else '', var=self.var,  ext=shapesName)
                print("INFO: \x1b[33mshapes\x1b[0m saved to:", outputFileName)
                shapesFile = ROOT.TFile.Open(outputFileName, "RECREATE")
                if dataGroupName and dataGroupName in groupedHistograms:
                    groupedHistograms[dataGroupName].SetDirectory(shapesFile)
                    ratioHistogram = groupedHistograms[dataGroupName].Clone()
                    ratioHistogram.SetName("ratio")
                    ratioHistogram.SetTitle("ratio")
                    ratioHistogram.Sumw2()
                    ratioHistogram.Divide(mcHistogram)
                    ratioHistogram.SetDirectory(shapesFile)
                mcHistogram.SetDirectory(shapesFile)

                for histogram in self.histograms:
                    histogram['histogram'].SetDirectory(shapesFile)

                try:
                    signals = self.readConfig(self.configSection, 'Signal', self.config.get('Plot_general','allSIG'))
                    print("INFO: signals=", signals) 
                    backgroundHistogram = NewStackMaker.sumHistograms(histograms=[histogram['histogram'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot and histogram['name'] not in signals], outputName='summedBackgroundHistograms')
                    signalHistogram = NewStackMaker.sumHistograms(histograms=[histogram['histogram'] for histogram in self.histograms if histogram['group'] in mcHistogramGroupsToPlot and histogram['name'] in signals], outputName='summedSignalHistograms')
                    backgroundHistogram.SetDirectory(shapesFile)
                    signalHistogram.SetDirectory(shapesFile)

                    sspb_sum = 0.0
                    q_a_sum = 0.0
                    s_sum = 0.0
                    b_sum = 0.0
                    for i in range(backgroundHistogram.GetXaxis().GetNbins()):
                        s = signalHistogram.GetBinContent(1+i)
                        b = backgroundHistogram.GetBinContent(1+i)
                        sspb = s/math.sqrt(s+b) if s+b > 0 else 0
                        sspb_sum += sspb*sspb
                        q_a = math.sqrt(2*((s+b)*math.log(1+s/b)-s)) if b>0 else 0
                        q_a_sum += q_a*q_a
                        s_sum += s
                        b_sum += b
                        print("{s: <12.3f}{b: <12.3f}{a: <12.3f}{q: <12.3f}".format(s=s, b=b, a=sspb, q=q_a))
                    print("-"*48)
                    q_a_sum = math.sqrt(q_a_sum)
                    sspb_sum = math.sqrt(sspb_sum)
                    print("{s: <12.3f}{b: <12.3f}{a: <12.3f}{q: <12.3f}".format(s=s_sum, b=b_sum, a=sspb_sum, q=q_a_sum))

                except Exception as e:
                    print("EXCEPTION:", e)
                shapesFile.Write()
            except Exception as e:
                print("ERROR: could not save shapes:", e)

        # print yield tables
        try:
            for histogram in self.histograms:
                self.histoCounts['weighted'][histogram['name']] = histogram['histogram'].Integral()
                self.histoCounts['unweighted'][histogram['name']] = histogram['histogram'].GetEntries()

            # print n events and weights for each sample
            keys = list(set(sorted([histogram['name'] for histogram in self.histograms])))
            keys.sort()
            for key in keys:
                print(key.ljust(50),("%d"%self.histoCounts['unweighted'][key]).ljust(10), ("%f"%self.histoCounts['weighted'][key]).ljust(10))

            ## print n events and weights for group
            self.groupCounts =  {'unweighted':{}, 'weighted': {}}
            for histogram in self.histograms:
                group_ = histogram['group']
                if not group_ in self.groupCounts['weighted']:
                    self.groupCounts['weighted'][group_] = histogram['histogram'].Integral()
                    self.groupCounts['unweighted'][group_] = histogram['histogram'].GetEntries()
                else:
                    self.groupCounts['weighted'][group_] += histogram['histogram'].Integral()
                    self.groupCounts['unweighted'][group_] += histogram['histogram'].GetEntries()

            print('--------------------')
            print('Printing the #events and yield by group')
            print('--------------------\n')

            #print('self.groupCounts[weighted] is',self.groupCounts['weighted'])
            keys = list(set(sorted([histogram['group'] for histogram in self.histograms])))
            keys.sort()
            for key in keys:
                print(key.ljust(50),("%d"%self.groupCounts['unweighted'][key]).ljust(10), ("%f"%self.groupCounts['weighted'][key]).ljust(10))

            wSumMC = sum([self.groupCounts['weighted'][key] for key in keys if key!=dataGroupName], 0)
            wSumData = sum([self.groupCounts['weighted'][key] for key in keys if key==dataGroupName], 0)

            print('--------------------')
            # makes no sense but keep this line for compatibility reason 
            wSum = wSumMC + wSumData
            print('total weighted events are:', wSum)

            print('total weighted MC events are:', wSumMC)
            print('total data yield is:', wSumData)
            print('data/MC ratio:', '%1.4f'%(wSumData/wSumMC))

        except Exception as e:
            print("ERROR: could not obtain counts (not implemented for TGraphAsymErrors yet)", e)

        # save data histogram
        if self.config.has_option('Plot_general','saveDataHistograms') and eval(self.config.get('Plot_general','saveDataHistograms')):
            if dataGroupName in groupedHistograms:
                outputFileName = self.outputFileTemplate.format(outputFolder=outputFolder, prefix='DATA_'+prefix, prefixSeparator='_' if len(prefix)>0 else '', var=self.var, ext="root")
                dataOutputFile = ROOT.TFile.Open(outputFileName, "recreate")
                groupedHistograms[dataGroupName].SetDirectory(dataOutputFile)
                dataOutputFile.Write()
                dataOutputFile.Close()
        print("INFO: stack created.")

