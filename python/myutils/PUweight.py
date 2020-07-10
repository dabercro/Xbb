#!/usr/bin/env python
import ROOT
from BranchTools import Collection
from BranchTools import AddCollectionsModule

class PUweight(AddCollectionsModule):

    def __init__(self, fileNameData, fileNameMC, fileNameDataUp=None, fileNameDataDown=None, nano=True, puWeightName='puWeight'):
        super(PUweight, self).__init__()
        self.version = 2
        self.nano = nano
        self.puWeightName = puWeightName
        self.fileNameData = fileNameData
        self.fileNameDataUp = fileNameDataUp
        self.fileNameDataDown = fileNameDataDown
        self.fileNameMC = fileNameMC
        self.ignoreHistogramMismatch = True

    def customInit(self, initVars):
        self.sample     = initVars['sample']
        self.sampleTree = initVars['sampleTree']
        self.config     = initVars['config']
        if self.sample.isMC():
            self.systematics = {
                    'Nominal': self.fileNameData, 
                    'Up': self.fileNameDataUp,
                    'Down': self.fileNameDataDown
                    }

            # load histograms from file
            rootFileData = {k: ROOT.TFile.Open(v, "read") if v else None for k,v in self.systematics.iteritems()}
            if not rootFileData['Nominal']:
                raise Exception("RootFileMissing")
            self.histogramData = {k: v.Get('pileup') if v else None for k,v in rootFileData.iteritems()}

            # if no file for MC is given, try to load MC pileup from 'autoPU' histogram (2017)
            if self.fileNameMC:
                rootFileMC = ROOT.TFile.Open(self.fileNameMC, "read")
                histogramMC = rootFileMC.Get('pileup')
            else:
                if 'autoPU' in self.sampleTree.histograms:
                    histogramMC = self.sampleTree.histograms['autoPU']
                    print "INFO: using autoPU histogram!"
                else:
                    print "\x1b[31mERROR: no MC PU file given and no autoPU histogram for MC found!\x1b[0m"
                    raise Exception("NoPUforMC")

            # normalize histograms
            for k,v in self.histogramData.iteritems():
                if v:
                    v.Scale(1.0/v.Integral())
            histogramMC.Scale(1.0/histogramMC.Integral())

            # check histogram compatibility
            for k,v in self.histogramData.iteritems():
                if v:
                    if not self.ignoreHistogramMismatch and (v.GetNbinsX()!=histogramMC.GetNbinsX() or v.GetXaxis().GetXmin()!=histogramMC.GetXaxis().GetXmin() or v.GetXaxis().GetXmax()!=histogramMC.GetXaxis().GetXmax()):
                        print "\x1b[31mERROR: ", k, "histograms not compatible! (bins, xmin, xmax)\x1b[0m"
                        print v.GetNbinsX(),histogramMC.GetNbinsX()
                        print v.GetXaxis().GetXmin(),histogramMC.GetXaxis().GetXmin()
                        print v.GetXaxis().GetXmax(),histogramMC.GetXaxis().GetXmax()
                        raise Exception("IncompatiblePUhistograms")
            
            if not self.histogramData['Nominal']:
                raise Exception("HistogramMissing")

            # determine PU range
            minPU = int(self.histogramData['Nominal'].GetXaxis().GetBinLowEdge(1))
            maxPU = int(self.histogramData['Nominal'].GetXaxis().GetBinLowEdge(self.histogramData['Nominal'].GetXaxis().GetNbins()))

            # compute table of PUweight vs. nTrueInt
            self.puWeightLUT = {}
            for k,v in self.histogramData.iteritems():
                self.puWeightLUT[k] = [1.0]*(maxPU+1)
            for pu in range(minPU, maxPU+1):
                nMC = histogramMC.GetBinContent(histogramMC.FindBin(pu))
                for k,v in self.histogramData.iteritems():
                    if v:
                        self.puWeightLUT[k][pu] = 1.0 * v.GetBinContent(v.FindBin(pu))/nMC  if nMC > 0 else 1.0
                print "PU = ",pu," weight = ", self.puWeightLUT['Nominal'][pu]
            self.puWeightLUTmin = minPU 
            self.puWeightLUTmax = maxPU 

            # add PU weight for all available systematic variations
            for k,v in self.histogramData.iteritems():
                if v:
                    self.addBranch(self.puWeightName + (('_'+k) if k != 'Nominal' else ''))
            print('dictionary_pu',self.puWeightLUT)

    def getPUweight(self, tree, syst=None):
        PUweight = 1
        pu = tree.Pileup_nTrueInt
        if pu < self.puWeightLUTmin and pu > 0:
            PUweight = self.puWeightLUT[syst][self.puWeightLUTmin]
        elif pu > self.puWeightLUTmax:
            PUweight = self.puWeightLUT[syst][self.puWeightLUTmax]
        else:
            PUweight = self.puWeightLUT[syst][int(pu)]
        return PUweight


    def processEvent(self, tree):
        if not self.hasBeenProcessed(tree) and not self.sample.isData():
            self.markProcessed(tree)
            for k,v in self.histogramData.iteritems():
                if v:
                    self._b(self.puWeightName + (('_'+k) if k != 'Nominal' else ''))[0] = self.getPUweight(tree, k)

