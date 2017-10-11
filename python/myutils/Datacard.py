#!/usr/bin/env python
from __future__ import print_function
import os, sys, ROOT, warnings, pickle
ROOT.gROOT.SetBatch(True)
from array import array
from math import sqrt
from copy import copy, deepcopy
#suppres the EvalInstace conversion warning bug
warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='creating converter.*' )
from optparse import OptionParser
from samplesclass import Sample
from sample_parser import ParseInfo 
import re
import json
import NewTreeCache as TreeCache
from sampleTree import SampleTree as SampleTree

class Datacard(object):
    def __init__(self, config, region, forceRedo=False):
        self.verbose = True
        self.forceRedo = forceRedo
        self.config = config
        VHbbNameSpace = config.get('VHbbNameSpace','library')
        returnCode = ROOT.gSystem.Load(VHbbNameSpace)
        if returnCode != 0:
            print ("\x1b[31mERROR: loading VHbbNameSpace failed with code %d\x1b[0m"%returnCode)
        else:
            print ("INFO: loaded VHbbNameSpace: %s"%VHbbNameSpace)
        
        self.region = region
        self.anaTag = config.get("Analysis", "tag")
        if self.anaTag not in ['7TeV', '8TeV', '13TeV']:
            raise Exception("anaTag %s unknown. Specify 8TeV or 7TeV or 13TeV in the general config"%self.anaTag)
        # Directories:
        self.Wdir=config.get('Directories','Wdir')
        self.vhbbpath=config.get('Directories','vhbbpath')
        self.samplesinfo=config.get('Directories','samplesinfo')
        self.path = config.get('Directories','dcSamples')
        self.cachedPath = self.config.get('Directories', 'tmpSamples')
        self.tmpPath = self.config.get('Directories', 'scratch')
        self.outpath=config.get('Directories','limits')
        self.optimisation= "" #TODO: opts.optimisation
        self.optimisation_training = False
        self.UseTrainSample = eval(config.get('Analysis','UseTrainSample'))
        if self.UseTrainSample:
            print ('Training events will be used')
        if not self.optimisation == '':
            print ('Preparing DC for BDT optimisaiton')
            self.optimisation_training = True
        print ('optimisation is', self.optimisation)
        try:
            os.stat(self.outpath)
        except:
            os.mkdir(self.outpath)
        
        # parse histogram config:
        self.treevar = config.get('dc:%s'%self.region,'var')
        print ('treevar is', self.treevar)
        self.name = config.get('dc:%s'%self.region,'wsVarName')
        if self.optimisation_training:
            self.treevar = self.optimisation+'.Nominal'
            self.name += '_'+ self.optimisation
            if self.UseTrainSample:
                self.name += '_Train'
        print ('again, treevar is', self.treevar)
        
        # set binning
        self.binning = {
                'nBins': int(config.get('dc:%s'%self.region,'range').split(',')[0]),
                'xMin': float(config.get('dc:%s'%self.region,'range').split(',')[1]),
                'xMax': float(config.get('dc:%s'%self.region,'range').split(',')[2]),
                }
        if self.verbose:
            print ("DEBUG: binning is ", self.binning)
        self.ROOToutname = config.get('dc:%s'%self.region,'dcName')
        if self.optimisation_training:
           self.ROOToutname += self.optimisation
           if self.UseTrainSample:
               self.ROOToutname += '_Train'
        self.RCut = config.get('dc:%s'%self.region,'cut')
        self.signals = eval('['+config.get('dc:%s'%self.region,'signal')+']') #TODO
        self.Datacardbin=config.get('dc:%s'%self.region,'dcBin')
        self.anType = config.get('dc:%s'%self.region,'type').upper()
        self.setup = eval(config.get('LimitGeneral','setup'))

        #new
        try:
            self.BDTmin = eval(config.get('LimitGeneral', 'BDTmin'))
        except:
            self.BDTmin = None

        self.keep_branches = eval(config.get('Branches', 'keep_branches'))
        
        #Systematics:
        if config.has_option('LimitGeneral', 'addSample_sys'):
            self.addSample_sys = eval(config.get('LimitGeneral','addSample_sys'))
            self.additionals = [self.addSample_sys[key] for key in self.addSample_sys]
        else:
            self.addSample_sys = None
            self.additionals = []

        #TODO: move to config!
        analysisSystematics = {
                'BDT': 'sys_BDT',
                'MJJ': 'sys_Mjj',
                'CR': 'sys_cr',
                }
        if self.anType in analysisSystematics:
            self.systematics = eval(config.get('LimitGeneral', analysisSystematics[self.anType]))
        else:
            print ("\x1b[31mEXIT: please specify if your datacards are BDT, Mjj or cr.\x1b[0m")
            raise Exception("InvalidDatacardSystematicsType")

        self.sysOptions = {}
        # define the options read directly from the config
        sysOptionNames = ['sys_cut_suffix', 'sys_weight_corr', 'decorrelate_sys_weight', 'sys_cut_include', 'sys_factor', 'sys_affecting', 'sys_lhe_affecting', 'rescaleSqrtN', 'toy', 'blind', 'addBlindingCut', 'change_shapes', 'Group', 'Dict', 'binstat',
                'rebin_active', 'ignore_stats', 'signal_inject', 'add_signal_as_bkg', 'systematicsnaming', 'weightF_sys', 'sample_sys_info', 'addSample_sys' 
                ]
        for sysOptionName in sysOptionNames:
            self.sysOptions[sysOptionName] = eval(config.get('LimitGeneral', sysOptionName)) if config.has_option('LimitGeneral', sysOptionName) else None
            print (" > \x1b[34m{name}\x1b[0m:{value}".format(name=sysOptionName.ljust(40), value=self.sysOptions[sysOptionName]))

        # read weights
        self.weightF = config.get('Weights','weightF')
        self.SBweight = None
        print ('before adding SBweight, weightF is',  self.weightF)
        if self.anType == 'MJJ':
            print ('Passed mJJ')
            if config.has_option('dc:%s'%self.region,'SBweight'):
                print ('passed config')
                self.SBweight = config.get('dc:%s'%self.region,'SBweight')
                self.weightF ='('+self.weightF+')*('+self.SBweight+')'
                print ('after adding SBweight, weightF is',  self.weightF)
            else:
                print ('NOT Passed config')

        self.TrainFlag = eval(config.get('Analysis','TrainFlag'))
        self.treecut = config.get('Cuts', self.RCut)
        
        # checks on read options
        #on control region cr never blind. Overwrite whatever is in the config
        if self.anType == 'CR':
            if self.sysOptions['blind']:
                print ('@WARNING: Changing blind to false since you are running for control region.')
            self.sysOptions['blind'] = False
            self.sysOptions['binstat'] = False

        if self.sysOptions['blind']:
            printc('red','', 'I AM BLINDED!')    
            
        if self.anType != 'BDT':
            if self.sysOptions['rebin_active']:
                print ('@WARNING: Changing rebin_active to false since you are running for control region.')
            self.sysOptions['rebin_active'] = False
        if self.sysOptions['add_signal_as_bkg']:
            self.setup.append(self.sysOptions['add_signal_as_bkg'])

        #--Setup--------------------------------------------------------------------
        #Assign Pt region for sys factors
        print ('Assign Pt region for sys factors')
        print ('================================\n')
        ptRegionsDict = {
                'HighPtLooseBTag': ['HighPtLooseBTag'],
                'HighPt': ['HighPt', 'ATLAS', 'MJJ'],
                'MedPt': ['MedPt'],
                'LowPt': ['lowPt'],
                'bin0': ['Vptbin0'],
                'bin1': ['Vptbin1'],
                'bin2': ['Vptbin2'],
                'bin150To200': ['Vptbin150To200'],
                'bin200ToInf': ['Vptbin200ToInf'],
                }
        self.ptRegion = [ptRegion for ptRegion,outputNames in ptRegionsDict.iteritems() if len([x for x in outputNames if x.upper() in self.ROOToutname.upper()])>0]
        print ("\x1b[34mptRegion:\x1b[0m", self.ptRegion)
       
        # TODO: move to config
        removeSystematicsDict = {
                'Zee': ['CMS_vhbb_eff_m_13TeV', 'CMS_vhbb_eff_m_trigger_Zll_13TeV', 'CMS_vhbb_eff_m_MVAID_Zll_13TeV', 'CMS_vhbb_eff_m_tracker_Zll_13TeV', 'CMS_vhbb_eff_m_ISO_Zll_13TeV'],
                'Zuu': ['CMS_vhbb_eff_e_13TeV', 'CMS_vhbb_eff_e_trigger_Zll_13TeV', 'CMS_vhbb_eff_e_MVAID_Zll_13TeV', 'CMS_vhbb_eff_e_MVAID_Zll_eta0_13TeV', 'CMS_vhbb_eff_e_MVAID_Zll_eta1_13TeV', 'CMS_vhbb_eff_e_tracker_Zll_13TeV'],
                }
        for outputName, removeSystematics in removeSystematicsDict.iteritems():
            if outputName in self.ROOToutname:
                self.sysOptions['weightF_sys'] = [x for x in self.sysOptions['weightF_sys'] if x not in removeSystematics]
        
         
        if self.TrainFlag:
            self.MC_rescale_factor = 2.
            print ('I RESCALE BY 2.0')
        else: 
            self.MC_rescale_factor = 1.

        #systematics up/down
        self.UD = ['Up','Down']

        print ('Parse the sample information')
        print ('============================\n')
        #Parse samples configuration
        self.samplesInfo = ParseInfo(self.samplesinfo, self.path)
        # get all the treeCut sets
        # create different sample Lists

        print ('Get the sample list')
        print ('===================\n')
        self.backgrounds = eval(config.get('dc:%s'%self.region,'background'))
        self.signals = eval(config.get('dc:%s'%self.region,'signal'))
        self.data_sample_names = eval(config.get('dc:%s'%self.region,'data'))

        #Systematics:
        if self.sysOptions['addSample_sys']: 
            self.additionals = [self.sysOptions['addSample_sys'][key] for key in self.sysOptions['addSample_sys']]
        else:
            self.additionals = []

        if self.sysOptions['sample_sys_info']:
            self.sample_sys_list = []#List of all the samples used for the sys. Those samples need to be skiped, except for corresponding sys
            #Extract list of sys samples
            for key, item in self.sysOptions['sample_sys_info'].iteritems():
                for item2 in item:
                    for sample_type in item2:
                        NOMsamplesys = sample_type[0]
                        noNom = False
                        for nomsample in NOMsamplesys: #This is for mergesyscachingdcsplit. Doesn't add the sys if nom is not present
                            print ('nomsample is', nomsample)
                            print ('signals+backgrounds are', self.signals+self.backgrounds)
                            if nomsample not in self.signals+self.backgrounds: noNom = True
                        if noNom: continue
                        DOWNsamplesys = sample_type[1]
                        UPsamplesys = sample_type[2]
                        self.sample_sys_list += DOWNsamplesys
                        self.sample_sys_list += UPsamplesys
            #self.additionals += self.sample_sys_list
        else:
            self.sample_sys_list = None
        #Create dictonary to "turn of" all the sample systematic (for nominal)
        self.sample_sys_dic = {}
        for sample_sys in self.sample_sys_list:
            self.sample_sys_dic[sample_sys] = False
        print ("\x1b[34msample_sys_list\x1b[0m =", self.sample_sys_list)
        
        self.MC_samples = self.signals + self.backgrounds + self.additionals
        
        self.samples = {
                'SIG': self.samplesInfo.get_samples(self.signals),
                'BKG': self.samplesInfo.get_samples(self.backgrounds),
                'DATA': self.samplesInfo.get_samples(self.data_sample_names),
                'ADD': self.samplesInfo.get_samples(self.additionals),
                }

        if self.verbose:
            print ('sample list')
            print ('===================\n')
            print (json.dumps(self.samples, sort_keys=True, indent=8, default=str))

        if self.verbose:
            print ('systematics')
            print ('===================\n')
            print (json.dumps(self.sysOptions['sys_affecting'], sort_keys=True, indent=8, default=str))

        #To prepare Histomaker
        all_samples_HM = self.samplesInfo.get_samples(self.signals+self.backgrounds+self.additionals)
        all_samples = self.samplesInfo.get_samples(self.MC_samples)

        self.systematicsDictionaryNominal = {
                'cut': self.treecut,
                'var': self.treevar,
                'name': self.name,
                'binning': self.binning,
                'weight': self.weightF,
                'countHisto': 'CountWeighted',
                'countbin': 0,
                'blind': self.sysOptions['blind'],
                'sysType': 'nominal',
                'SBweight': self.SBweight,
                'sample_sys_dic': self.sample_sys_dic,
                }
        
        # contains all systematicDictionaries
        self.systematicsList = [self.systematicsDictionaryNominal]

        print ('Assign the systematics')
        print ('======================\n')
        #shape systematics
        for syst in self.systematics:
            for Q in self.UD:
                systematicsDictionary = deepcopy(self.systematicsDictionaryNominal)
                systematicsDictionary.update({
                        'cut': self.getSystematicsCut(syst, Q),
                        'var': self.getSystematicsVar(syst, Q),
                        'weight': self.getSystematicsWeight(syst, Q),
                        'sysType': 'shape',
                    })
                self.systematicsList.append(systematicsDictionary)
        
        if self.verbose:
            print ('systematics dict')
            print ('===================\n')
            print (json.dumps(self.systematicsList, sort_keys=True, indent=8, default=str))
                    
    def getSystematicsList(self, isData=False):
        if isData:
            return [x for x in self.systematicsList if x['sysType'] == 'nominal']
        else:
            return self.systematicsList

    def getRegion(self):
        return self.region

    def cacheSample(self, sampleToCache):
        treeCaches = []
        self.sampleTree = None
        
        # cuts
        allSamples = sum([y for x,y in self.samples.iteritems()], [])
        subsamples = [x for x in allSamples if x.identifier == sampleToCache]
        for sample in subsamples:
            # combine subcut and systematics cut
            sampleCuts = ['(%s)&&(%s)'%(sample.subcut,x['cut']) for x in self.systematicsList]

            # add cache object
            tc = TreeCache.TreeCache(
                sample=sample.name,
                cutList=sampleCuts,
                cutSequenceMode='OR',
                inputFolder=self.path,
                tmpFolder=self.tmpPath,
                outputFolder=self.cachedPath,
                debug=True
            )

            # check if this part of the sample is already cached
            isCached = tc.isCached()
            print ("IS CACHED?", isCached)
            if not isCached or self.forceRedo:
                if isCached:
                    tc.deleteCachedFiles()

                # for the first sample which comes from this files, load the tree
                if not self.sampleTree:
                    self.sampleTree = SampleTree({'name': sample.identifier, 'folder': self.path})
                    if not self.sampleTree or not self.sampleTree.tree:
                        print ("\x1b[31mERROR: creation of sample tree failed!!\x1b[0m")
                        raise Exception("CreationOfSampleTreeFailed")
                treeCaches.append(tc.setSampleTree(self.sampleTree).cache())

        if len(treeCaches) > 0:
            # run on the tree
            self.sampleTree.process()
        else:
            print ("nothing to do!")

    def getSystematicsVar(self, syst, Q):
        treevar = self.treevar
        #replace tree variable
        if self.anType == 'BDT':
            if not 'UD' in syst:
                print ('treevar was', _treevar)
                treevar = treevar.replace('.Nominal','.%s_%s'%(syst,Q))
                print ('.nominal by','.%s_%s'%(syst,Q))
            else:
                treevar = treevar.replace('.nominal','.%s'%(syst.replace('UD',Q)))
                print ('.nominal by','.%s'%(syst.replace('UD',Q)))
        elif self.anType == 'MJJ':
            if not 'UD' in syst:
                print ('treevar was', _treevar)
                treevar = treevar.replace('_reg_mass','_reg_mass_corr%s%s'%(syst,Q))
                print ('_reg_mass','_reg_mass_corr%s%s'%(syst,Q))
            else:
                print ('\x1b[31m@ERROR: Why is there UD in sys ? Abort\x1b[0m')
                raise Exception("DcSystMjjContainsUD")
        return treevar 

    def getSystematicsWeight(self, syst, Q):
        weight = self.weightF
        if syst in self.sysOptions['sys_weight_corr']:
            weightName = '{sysWeightCorr}_{ud}'.format(sysWeightCorr=self.sysOptions['sys_weight_corr'][syst], ud=Q.upper())
            weight = self.config.get('Weights', weightName)
            print ('weightName is \x1b[34m{weightName}\x1b[0m'.format(weightName=weightName))
            print ('weight is \x1b[34m{weight}\x1b[0m'.format(weight=weight))
        return weight

    def getSystematicsCut(self, syst, Q):
        cut = self.treecut
        new_cut_list = self.sysOptions['sys_cut_suffix'][syst] if isinstance(self.sysOptions['sys_cut_suffix'][syst], list) else [self.sysOptions['sys_cut_suffix'][syst]] 
        for new_cut in new_cut_list:
            if not new_cut == 'nominal':
                old_str, new_str = new_cut.split('>')
                cut = cut.replace(old_str, new_str.replace('SYS', syst).replace('UD', Q).replace('?', Q))
        return cut
