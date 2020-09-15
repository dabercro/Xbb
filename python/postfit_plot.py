#!/usr/bin/env python
from __future__ import print_function
from optparse import OptionParser
from collections import defaultdict
import ROOT
ROOT.gROOT.SetBatch(True)

from myutils.NewStackMaker import NewStackMaker as StackMaker
from myutils.samplesclass import Sample
from myutils.sample_parser import ParseInfo
from myutils.XbbTools import XbbTools

import os,sys
import array
import math

class PostfitPlotter(object):

    def __init__(self, config, region, vars=None, directory="shapes_fit_s", title=None, blind=True):
        self.config = config
        self.region = region if type(region) == list else [region]
        self.fitRegions = eval(self.config.get('Fit', 'regions'))
        self.dcRegion = [self.fitRegions[x] for x in self.region] 
        self.var = None 
        self.directory = directory
        self.shapesFile = None
        self.plotPath = config.get('Directories', 'plotpath')
        self.title = None # "CMS #scale[0.8]{work in progress}"
        self.blindBins = []
        self.blind = blind
        self.plotText = [""]
        if self.config.has_option('Fit', 'plotText'):
            self.plotText = eval(self.config.get('Fit', 'plotText'))

        # plot Data/B and (S+B)/B instead of Data/(S+B)
        self.dataOverBackground = self.config.has_option('Fit', 'plotDataOverBackground') and eval(self.config.get('Fit', 'plotDataOverBackground'))
        self.ratioRange = None 

        if self.config.has_section('Fit:'+self.region[0]):
            if self.config.has_option('Fit:'+self.region[0], 'var'):
                self.var = self.config.get('Fit:'+self.region[0], 'var')
            if self.blind and self.config.has_option('Fit:'+self.region[0], 'blindBins'):
                self.blindBins = eval(self.config.get('Fit:'+self.region[0], 'blindBins'))
            if self.config.has_option('Fit:'+self.region[0], 'plotText'):
                self.plotText += eval(self.config.get('Fit:'+self.region[0],'plotText'))
            if self.config.has_option('Fit:'+self.region[0], 'plotDataOverBackground'):
                self.dataOverBackground = eval(self.config.get('Fit:'+self.region[0], 'plotDataOverBackground'))
            if self.config.has_option('Fit:'+self.region[0], 'ratioRange'):
                self.ratioRange = eval(self.config.get('Fit:'+self.region[0], 'ratioRange'))

        self.samplesInfo = ParseInfo(samples_path=config.get('Directories', 'dcSamples'), config=self.config)
        self.sampleGroupDict = eval(self.config.get('LimitGeneral', 'Group')) if self.config.has_option('LimitGeneral', 'Group') else {}

        # make dictionary to get sample type ('SIG'/'BKG') for each sample group
        self.sampleGroupTypeDict = {}
        for x in self.samplesInfo:
            sampleGroup = self.getSampleGroup(x)
            if sampleGroup:
                if sampleGroup not in self.sampleGroupTypeDict:
                    self.sampleGroupTypeDict[sampleGroup] = []
                self.sampleGroupTypeDict[sampleGroup].append(x.type)

        for k,v in self.sampleGroupTypeDict.items():
            v_unique = list(set(v))
            if len(v_unique) > 1:
                print(k, v_unique)
                raise Exception("ConfigError")
            self.sampleGroupTypeDict[k] = v_unique[0] if len(v_unique) > 0 else None

    
    def getSampleGroup(self, sample):
        # if Group dictionary is used, prioritize it over group from config
        if sample.name in self.sampleGroupDict: 
            sampleGroup = self.sampleGroupDict[sample.name]
        else:
            sampleGroup = sample.group
        return sampleGroup

    def prepare(self):
        self.dcDict = eval(self.config.get('LimitGeneral','Dict'))
        self.reverseDcDict = {v:k for k,v in self.dcDict.items()}
        self.regionDict = eval(self.config.get('Fit', 'regions'))
        self.reverseRegionDict = {v:k for k,v in self.regionDict.items()}
        self.setup = eval(self.config.get('LimitGeneral','setup'))
        self.overwriteSignalDefinition = eval(self.config.get('LimitGeneral','overwriteSignalDefinition')) if self.config.has_option('LimitGeneral','overwriteSignalDefinition') else False
        if self.overwriteSignalDefinition and self.config.has_option('LimitGeneral','setupSignals'):
            self.setupSignals = eval(self.config.get('LimitGeneral','setupSignals'))
        else:
            self.setupSignals = [x for x in self.setup if x in self.sampleGroupTypeDict and self.sampleGroupTypeDict[x] == 'SIG']
        print("DEBUG: signals:", self.setupSignals)

        shapesFileName = self.config.get('Fit', 'FitDiagnosticsDump')
        self.shapesFile = ROOT.TFile.Open(shapesFileName, "READ")
        if self.shapesFile is None or self.shapesFile.IsZombie() or self.shapesFile.TestBit(ROOT.TFile.kRecovered):
            print("FILENAME:", shapesFileName)
            raise Exception("FileError")
        print("DEBUG: reverse dict:", self.reverseDcDict)

    # process is datacard name convention, convert it back to Xbb process convention if necessary
    def getNameFromDcname(self, process):
        return self.reverseDcDict[process] if process in self.reverseDcDict else process

    def getDcnameFromName(self, processName):
        return self.dcDict[processName]

    def getProcesses(self):
        return self.setup

    def getHistogramNames(self, process):
        processName = self.getNameFromDcname(process) 
        
        # use pre-fit signal strength for plotting if blind
        if self.blind and processName in self.setupSignals:
            print("INFO: \x1b[31mBLIND: prefit shapes are plotted for signal!\x1b[0m")
            histogramName = ['shapes_prefit/{dcRegion}/' + process, '{dcRegion}_prefit/' + process]
        else:
            fitType = 'postfit' if 'fit_s' in self.directory else 'prefit'
            histogramName = [self.directory + '/{dcRegion}/' + process, '{dcRegion}_' + fitType + '/' + process]

        print("DEBUG: get shape ", histogramName)
        return histogramName

    def getBins(self, region):
        nBins = -1
        bins = None
        if self.config.has_section('Fit:'+region) and self.config.has_option('Fit:'+region, 'nBins'):
            nBins   = eval(self.config.get('Fit:'+region, 'nBins'))
        try:
            varName = self.config.get('Fit:'+region, 'var')
            bins    = XbbTools.getPlotVariableBins(varName, config=self.config)
        except Exception as e:
            print("ERROR: Can't get plot variable information:", e)

        if bins is None:
            bins = XbbTools.uniformBins(nBins,0.0,1.0)

        return bins

    def getShape(self, process):

        histogramNameTemplates = self.getHistogramNames(process)
        histograms = []
        histogramBins = []
        for x in self.dcRegion:
            found = False
            for histogramNameTemplate in histogramNameTemplates:
                print("REGION:", x, histogramNameTemplate.format(dcRegion=x), " FROM ", self.shapesFile)
                histogram = self.shapesFile.Get(histogramNameTemplate.format(dcRegion=x))
                print("H:", histogramNameTemplate, histogram)
                if isinstance(histogram, ROOT.TGraphAsymmErrors) or isinstance(histogram, ROOT.TH1) or isinstance(histogram, ROOT.TH1D) or isinstance(histogram, ROOT.TH1F):
                    print("->OK")
                    found = True
                    break
            if not found:
                return None

            histograms.append(histogram)
            region = self.reverseRegionDict[x] 
            bins   = self.getBins(region)
            print("BINS:", bins, type(bins)) 
            histogramBins.append(bins)

        if len(histograms) == 1:
            print("fix bin boundaries:", process)
            histogram = ROOT.TH1D("h1", "", len(histogramBins[0])-1, histogramBins[0])
            histogram.SetDirectory(0)
            histogram.Sumw2()

            if isinstance(histograms[0], ROOT.TGraphAsymmErrors):
                histogram = histograms[0]
                for i in range(len(histogramBins[0])-1):
                    p_x = array.array('d', [0.0])
                    p_y = array.array('d', [0.0])
                    histogram.GetPoint(i, p_x, p_y)
                    histogram.SetPoint(i, 0.5*(histogramBins[0][i]+histogramBins[0][i+1]), p_y[0])
                    histogram.SetPointEXhigh(i, 0.5*(histogramBins[0][i+1]-histogramBins[0][i]))
                    histogram.SetPointEXlow(i, 0.5*(histogramBins[0][i+1]-histogramBins[0][i]))
            else:
                try:
                    for i in range(len(histogramBins[0])-1):
                        histogram.SetBinContent(i+1, histograms[0].GetBinContent(i+1))
                        histogram.SetBinError(i+1, histograms[0].GetBinError(i+1))
                except Exception as e:
                    print("EXCEPTION:",e)

        else:
            nBinsTotal = sum([len(x)-1 for x in histogramBins])
            histogram = ROOT.TH1D("h1", "", nBinsTotal, 0, nBinsTotal)
            histogram.SetDirectory(0)

            iBin = 1
            if isinstance(histograms[0], ROOT.TGraphAsymmErrors):
                histogram = ROOT.TGraphAsymmErrors()
                tlist = ROOT.TList()
                for i,h in enumerate(histograms):
                    h.Set(histogramBins[i])
                    tlist.Add(h)
                histogram.Merge(tlist)
                offset = 0
                for i in range(histogram.GetN()):
                    p_x = array.array('d', [0.0])
                    p_y = array.array('d', [0.0])
                    histogram.GetPoint(i, p_x, p_y)
                    if i == 0:
                        offset = p_x[0]
                    print("P:",i,p_x[0],p_y[0],"-->", i, p_y[0])
                    histogram.SetPoint(i, i + offset, p_y[0])
            else:
                for i,h in enumerate(histograms):
                    if h:
                        for n in range(histogramBins[i]):
                            histogram.SetBinContent(iBin, h.GetBinContent(1+n))
                            histogram.SetBinError(iBin, h.GetBinError(1+n))
                            iBin += 1

        #histogram = self.shapesFile.Get(histogramName)
        #print("DEBUG: -->", histogram, histogram.Integral())
        return histogram

    def setBinRange(self, histogram, nBins, rangeMin=0.0, rangeMax=1.0):
        newHistogram = ROOT.TH1D("nh", "", nBins, rangeMin, rangeMax)
        newHistogram.SetDirectory(0)
        newHistogram.Sumw2()
        for i in range(nBins):
            newHistogram.SetBinContent(i+1, histogram.GetBinContent(i+1))
            newHistogram.SetBinError(i+1, histogram.GetBinError(i+1))
        return newHistogram

    def getDcProcessName(self, process):
        return self.dcDict[process] if process in self.dcDict else process

    def run(self):

        # if variable definition not given explicitly
        if self.var is None:
            print("make auto variable...")
            self.var = "__auto"
            self.varSection = "plotDef:" + self.var
            if not self.config.has_section(self.varSection):
                self.config.add_section(self.varSection)
            nBins = next( (self.getShape(self.getDcProcessName(process)).GetXaxis().GetNbins() for process in self.setup), 15)
            self.config.set(self.varSection, 'nBins', '%d'%nBins)
            self.config.set(self.varSection, 'min', '0')
            self.config.set(self.varSection, 'max', '%d'%nBins) 
        print("INFO: create StackMaker...")

        # overwrite treeVar since we plot from histograms and not trees 
        self.varSection = "plotDef:" + self.var
        if not self.config.has_option(self.varSection, "relPath"):
            self.config.set(self.varSection, "relPath", "1")

        self.stack = StackMaker(self.config, self.var, self.region[0], True, self.setup, '_', title=self.title, configSectionPrefix="Fit")
        self.stack.setPlotText(self.plotText)

        # add MC
        print("INFO: setup = \x1b[31m", self.setup, "\x1b[0m")
        for process in self.setup:
            print("\x1b[32mINFO: getShape:", process, "\x1b[0m")
            histogram = self.getShape(self.getDcProcessName(process))
            print("\x1b[32mINFO: OK getShape:", process, "\x1b[0m")

            if histogram:
                self.stack.histograms.append({
                        'name':      process, 
                        'histogram': histogram, 
                        'group':     process,
                        'signal':    process in self.setupSignals
                    })
            else:
                print("\x1b[31mINFO: empty: ",process,"\x1b[0m")

        print("\x1b[32mINFO: all shapes!\x1b[0m")
        #for histogram in self.stack.histograms:
        #    print("DEBUG:", histogram['name'])
        #    for i in range(histogram['histogram'].GetXaxis().GetNbins()):
        #        print(" > ", histogram['histogram'].GetBinContent(1+i))

        # dump S/B
        #shapeS = self.getShape("total_signal")
        #shapeB = self.getShape("total_background")
        shapeS = StackMaker.sumHistograms(histograms=[histogram['histogram'] for histogram in self.stack.histograms if histogram['group'] in self.setupSignals], outputName='totalSignal') 
        shapeB = StackMaker.sumHistograms(histograms=[histogram['histogram'] for histogram in self.stack.histograms if histogram['group'] not in self.setupSignals], outputName='totalBackground') 
        print("name:", shapeB.GetName(), " / ", self.directory)
        print("file:", self.shapesFile)
        print("shape:", shapeS, shapeB)
        print(self.region,"-"*40)
        print( ("bin").ljust(5), ("S").ljust(8), ("B").ljust(8))
        sum_s_over_sqrtb = 0.0
        sum_s = 0.0
        sum_b = 0.0
        sum_med_z_a = 0.0
        if shapeS is not None:
            for i in range(shapeS.GetXaxis().GetNbins()):
                s = shapeS.GetBinContent(1+i)
                b = shapeB.GetBinContent(1+i)
                print( ("%d"%i).ljust(5), ("%1.2f"%s).ljust(8), ("%1.2f"%b).ljust(8))

                sum_s_over_sqrtb += s*s/b if b>0 else 0
                sum_s += s
                sum_b += b
                sum_med_z_a += 2.0* ( (s+b) * math.log(1.0 + s/b) - s) if s > 0 else 0.0
        sum_s_over_sqrtb = math.sqrt(sum_s_over_sqrtb)
        sum_med_z_a = math.sqrt(sum_med_z_a)
        print("TOTAL: s/sqrt(b):", "%1.3f"%sum_s_over_sqrtb, " s:",sum_s, " b:",sum_b, " med[Z]=sqrt(q_A)=", sum_med_z_a)
        
        # add DATA
        dataHistogram = self.getShape("data") 
        print("DEBUG: \x1b[32m data histogram=", dataHistogram, dataHistogram == None, "\x1b[0m")
        if dataHistogram == None:
            print("-> data histogram not found, trying alternative names")
            dataHistogram = self.getShape("data_obs") 

        # blind last few bins
        dataIntegral = 0
        try:
            if isinstance(dataHistogram, ROOT.TGraphAsymmErrors) or isinstance(dataHistogram, ROOT.TGraphErrors):
                for i in range(dataHistogram.GetN()):
                    dataHistogram.GetPoint(i, pointX, pointY)
                    dataIntegral += pointY[0]
                    if int(pointX[0]+1) in self.blindBins:
                        dataHistogram.SetPoint(i, -100, -100)
            else:
                dataIntegral = dataHistogram.Integral()
                for i in self.blindBins:
                    dataHistogram.SetBinContent(1+i,0)
                    dataHistogram.SetBinError(1+i,0)
        except Exception as e:
            print("ERROR: BlindDataHistogram",e)

        print("DATA:", dataIntegral, "MC:", sum_s+sum_b)

        # style data
        dataHistogram.SetMarkerColor(ROOT.kBlack)
        dataHistogram.SetMarkerStyle(20)

        self.stack.histograms.append({
                'name': 'DATA', 
                'histogram': dataHistogram,
                'group': 'DATA'
            })

        try:
            totalError = dataHistogram = self.getShape("TotalProcs")
            if totalError != None:
                self.stack.setTotalError(totalError)
        except:
            pass

        # draw
        self.stack.Draw(outputFolder=self.plotPath, prefix='{region}__{var}_'.format(region=self.region[0], var=self.directory), dataOverBackground=self.dataOverBackground, ratioRange=self.ratioRange)


if __name__ == "__main__":
    # read arguments
    argv = sys.argv
    parser = OptionParser()
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                              help="Verbose mode.")
    parser.add_option("-C", "--config", dest="config", default=[], action="append",
                          help="configuration file")
    parser.add_option("-r","--regions", dest="regions", default='',
                          help="cut region identifiers, separated by comma")
    parser.add_option("-t","--type", dest="fitType", default='',
                          help="shapes_prefit, shapes_fit_b, shapes_fit_s")
    parser.add_option("-u", "--unblind", action="store_true", dest="unblind", default=False,
                              help="Unblind")
    parser.add_option("-b","--bins", dest="bins", default='',
            help="for eff. scale factor computation: restrict to this list of comma separated bins (Starting at bin 1)")

    (opts, args) = parser.parse_args(argv)
    if opts.config == "":
            opts.config = ["config"]

    # Import after configure to get help message
    from myutils import BetterConfigParser, mvainfo, ParseInfo
    config = BetterConfigParser()
    config.read(opts.config)

    if len(opts.fitType) < 1:
        if config.has_option('Fit','FitType'):
            opts.fitType = config.get('Fit','FitType')
        else:
            opts.fitType = "shapes_prefit"

    # run plotter
    if len(opts.regions) < 1 or opts.regions.strip() == 'None':
        regions = eval(config.get('Fit', 'regions')).keys()
    else:
        regions = opts.regions.strip().split(',')
        for i in range(len(regions)):
            if '+' in regions[i]:
                regions[i] = regions[i].split('+')
        print("REGIONS:", regions)
    
    scaleFactorTable = []
    plotCommands = []
    for region in regions:
        print("INFO: ----------")
        print("INFO: region:", region)
        try:
            plotter = PostfitPlotter(config=config, region=region, directory=opts.fitType, blind=not opts.unblind)
            plotter.prepare()
            print("\x1b[32mINFO: run!\x1b[0m")
            plotter.run()
        except Exception as e:
            print("ERROR: plot:", e)

        print("INFO: compute SFs:")
        # COMPUTE effective scale factors
        if opts.fitType == 'shapes_fit_s':
            plotter_prefit = PostfitPlotter(config=config, region=region, directory="shapes_prefit", blind=not opts.unblind)
            plotter_prefit.prepare()
        
            regionSF = defaultdict(lambda: 1.0)

            for process in plotter.getProcesses():
                try:
                    histogram_postfit = plotter.getShape(plotter.getDcnameFromName(process))
                    histogram_prefit  = plotter_prefit.getShape(plotter.getDcnameFromName(process))

                    nBins_prefit = histogram_prefit.GetXaxis().GetNbins()
                    nBins_postfit = histogram_postfit.GetXaxis().GetNbins()
                    assert nBins_prefit == nBins_postfit

                    # restrict to bin list
                    if len(opts.bins.strip()) > 0:
                        binList = [int(x) for x in opts.bins.strip().split(',')]
                        for i in range(nBins_prefit):
                            if i+1 not in binList:
                                histogram_prefit.SetBinContent(i+1, 0)
                                histogram_postfit.SetBinContent(i+1, 0)

                    histogram_postfit.Rebin(nBins_postfit)
                    histogram_prefit.Rebin(nBins_prefit)
                    
                    print("\x1b[35m", region, process, histogram_prefit.GetBinContent(1), histogram_postfit.GetBinContent(1), "\x1b[0m")
                    histogram_postfit.Divide(histogram_prefit)

                    scaleFactorTable.append("{region} {process} prefit/postfit = {scale:.2f} +/- {error:.2f}".format(region=region, process=process, scale=histogram_postfit.GetBinContent(1), error=histogram_postfit.GetBinError(1))) 
                    regionSF[process] = histogram_postfit.GetBinContent(1)
                except Exception as e:
                    print("WARNING:", e)

            # VHbb specific plot commands for simplifity
            try:
                regionFormatted = ",".join(region) if type(region) == list else region
                plotCommand = "./submit.py -J runplot --parallel=8 --regions '{region}' --set='General.SF_TT={SF[TT]};General.SF_ZJets=[{SF[ZJets_0b_udsg]},{SF[ZJets_0b_c]},{SF[ZJets_1b]},{SF[ZJets_2b]}];General.SF_WJets=[{SF[WJets_0b_udsg]},{SF[WJets_0b_c]},{SF[WJets_1b]},{SF[WJets_2b]}]'".format(region=regionFormatted, SF=regionSF)
                plotCommands.append(plotCommand)
            except Exception as e:
                print("WARNING:", e)

    # PRINT effective scale factors
    if opts.fitType == 'shapes_fit_s':
        print("INFO: --------------------------------------------------------------------------------")
        print("INFO: effective scale factors (including combined effects of all uncertainties)")
        print("INFO: --------------------------------------------------------------------------------")

        for line in scaleFactorTable:
            print("INFO:", line)
        print("INFO: --------------------------------------------------------------------------------")
        for line in plotCommands:
            print("", line)
        print("INFO: --------------------------------------------------------------------------------")

