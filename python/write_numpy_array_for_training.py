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
import numpy as np
import math
from copy import deepcopy
import gzip
import h5py
import json
import datetime

class SampleTreesToNumpyConverter(object):

    def __init__(self, config, mvaName, useSyst=True, useWeightSyst=True, testRun=False, includeData=False):
        self.mvaName = mvaName
        self.includeData = includeData
        VHbbNameSpace = config.get('VHbbNameSpace', 'library')
        ROOT.gSystem.Load(VHbbNameSpace)
        self.dataFormatVersion = 4
        self.sampleTrees = []
        self.config = config
        self.testRun = testRun
        if config.has_option('Directories', 'trainingSamples'):
            self.samplesPath = config.get('Directories', 'trainingSamples')
        else:
            self.samplesPath = config.get('Directories', 'MVAin')
        self.samplesInfo = ParseInfo(samples_path=self.samplesPath, config=self.config) 

        # region
        self.treeCutName = config.get(mvaName, 'treeCut') if config.has_option(mvaName, 'treeCut') else mvaName
        self.treeCut = config.get('Cuts', self.treeCutName)

        # split in train/eval sets
        self.trainCut = config.get('Cuts', 'TrainCut') 
        self.evalCut = config.get('Cuts', 'EvalCut')

        # rescale MC by 2 because of train/eval split
        self.globalRescale = 2.0

        # variables and systematics
        self.treeVarSet = config.get(mvaName, 'treeVarSet')
        self.MVA_Vars = {'Nominal': [x for x in config.get(self.treeVarSet, 'Nominal').strip().split(' ') if len(x.strip()) > 0]}

        self.weightSYS = []
        self.weightSYSweights = {}

        self.systematics = []
        if useSyst:
            print('INFO: use systematics in training!')
            self.systList = eval(self.config.get(mvaName, 'systematics')) if self.config.has_option(mvaName, 'systematics') else []
            if config.has_option(mvaName, 'systematics'):
                systematicsString = config.get(mvaName, 'systematics').strip()
                if systematicsString.startswith('['):
                    self.systList = eval(systematicsString)
                else:
                    self.systList = systematicsString.split(' ')
            else:
                self.systList = []

            for syst in self.systList:
                systNameUp   = syst+'_UP'   if self.config.has_option('Weights',syst+'_UP')   else syst+'_Up'
                systNameDown = syst+'_DOWN' if self.config.has_option('Weights',syst+'_DOWN') else syst+'_Down'

                self.systematics.append({
                    'name': syst,
                    'U': self.config.get('Weights', systNameUp),
                    'D': self.config.get('Weights', systNameDown),
                    })

        # default: signal vs. background
        self.sampleNames = {
                    'SIG_ALL': eval(self.config.get(mvaName, 'signals')) if self.config.has_option(mvaName, 'signals') else eval(self.config.get('Plot_general', 'allSIG')),
                    'BKG_ALL': eval(self.config.get(mvaName, 'backgrounds')) if self.config.has_option(mvaName, 'backgrounds') else  eval(self.config.get('Plot_general', 'allBKG')),
                }
        # for multi-output classifiers load dictionary from config
        self.categories = None
        if self.config.has_option(mvaName, 'classDict'):
            self.sampleNames = eval(self.config.get(mvaName, 'classDict'))
            #print('checking.........',eval(self.config.get(mvaName, 'classDict')))
            #print('checking.........',self.config.get(mvaName, 'classDict'))
            self.samples = {category: self.samplesInfo.get_samples(samples) for category,samples in self.sampleNames.iteritems()} 
            self.categories = self.samples.keys()
            print("classes dict:", self.sampleNames)
        elif self.config.has_option(mvaName, 'classes'):
            self.sampleNames = dict(eval(self.config.get(mvaName, 'classes')))
            self.categories = [x[0] for x in eval(self.config.get(mvaName, 'classes'))]
        self.samples = {category: self.samplesInfo.get_samples(samples) for category,samples in self.sampleNames.iteritems()}
        if not self.categories:
            self.categories = self.samples.keys()

        # DATA
        if self.config.has_option(mvaName, 'includeData'):
            self.includeData = eval(self.config.get(mvaName, 'includeData'))

        self.dataSamples = []
        if self.includeData:
            if not self.config.has_option(mvaName, 'data'):
                print("\x1b[31mERROR: in training.ini, the option 'data' has to be specified for the MVA\x1b[0m")
                raise Exception("ConfigError")
            self.dataSampleNames =  eval(self.config.get(mvaName, 'data')) if self.config.has_option(mvaName, 'data') else eval(self.config.get('Plot_general', 'Data')) 
            print("INFO: sample names for DATA are:", self.dataSampleNames)
            self.dataSamples = self.samplesInfo.get_samples(self.dataSampleNames)
            print("\x1b[32mINFO: added DATA:", [x.identifier for x in self.dataSamples],"\x1b[0m")
        else:
            print("INFO: DATA not added, use --include-data to add it to the h5 file")

        if self.testRun:
            print("\x1b[31mDEBUG: TEST-RUN, using only small subset of samples!\x1b[0m")


    def run(self):
        # ----------------------------------------------------------------------------------------------------------------------
        # add sig/bkg x training/testing trees
        # ----------------------------------------------------------------------------------------------------------------------
        categories = self.categories 
        if categories:
            print("categories:")
            for i,category in enumerate(categories):
                print(" ",i,":", category)
        datasetParts = {'train': self.trainCut, 'test': self.evalCut}

        systematics = self.systematics
        arrayLists = {datasetName:[] for datasetName in datasetParts.iterkeys()}
        #arrayLists_sys = {x: {datasetName:[] for datasetName in datasetParts.iterkeys()} for x in systematics}
        weightLists = {datasetName:[] for datasetName in datasetParts.iterkeys()}
        targetLists = {datasetName:[] for datasetName in datasetParts.iterkeys()}
        sampleIndexList = {datasetName:[] for datasetName in datasetParts.iterkeys()}

        weightListsSYS = {x: {datasetName:[] for datasetName in datasetParts.iterkeys()} for x in self.weightSYS} 
        
        # standard weight expression
        weightF = self.config.get('Weights','weightF')

        weightListSYStotal = {datasetName:[] for datasetName in datasetParts.iterkeys()}
        eventNumberStat = []

        for iCat,category in enumerate(categories):
            if self.testRun:
                self.samples[category] = self.samples[category][0:1]
            for j,sample in enumerate(self.samples[category]):
                print ('*'*80,'\n%s (category %d/%d sample %d/%d)\n'%(sample, iCat+1, len(categories), j+1, len(self.samples[category])),'*'*80)
                for datasetName, additionalCut in datasetParts.iteritems():
                    # cuts
                    sampleCuts = [sample.subcut]
                    if additionalCut:
                        sampleCuts.append(additionalCut)
                    # cut from the mva region
                    if self.treeCut:
                        sampleCuts.append(self.treeCut)

                    # get ROOT tree for selected sample & region cut
                    tc = TreeCache.TreeCache(
                            sample=sample,
                            cutList=sampleCuts,
                            inputFolder=self.samplesPath,
                            config=self.config,
                            debug=True
                        )
                    sampleTree = tc.getTree()
                    if sampleTree:
                        treeScale = sampleTree.getScale(sample) * self.globalRescale
                        print ('scale:', treeScale)
                        
                        # initialize numpy array
                        nSamples = sampleTree.GetEntries()
                        features = self.MVA_Vars['Nominal']
                        #features_sys = {x: self.MVA_Vars[x] for x in systematics} 
                        nFeatures = len(features) 
                        print('nFeatures:', nFeatures)
                        inputData = np.zeros((nSamples, nFeatures), dtype=np.float32)
                        #inputData_sys = {x: np.zeros((nSamples, nFeatures), dtype=np.float32) for x in systematics}

                        sampleIndex = sample.index

                        # initialize formulas for ROOT tree
                        for feature in features:
                            sampleTree.addFormula(feature)
                        #for k, features_s in features_sys.iteritems():
                        #    for feature in features_s:
                        #        sampleTree.addFormula(feature)
                        sampleTree.addFormula(weightF)
                        #for syst in self.weightSYS:
                        #    sampleTree.addFormula(self.weightSYSweights[syst])
                        for syst in self.systematics:
                            sampleTree.addFormula(syst['U'])
                            sampleTree.addFormula(syst['D'])

                        useSpecialWeight = self.config.has_option('Weights', 'useSpecialWeight') and eval(self.config.get('Weights', 'useSpecialWeight')) 
                        if useSpecialWeight:
                            sampleTree.addFormula(sample.specialweight)

                        sumOfWeights = 0.0
                        nMCevents = 0
                        # fill numpy array from ROOT tree
                        for i, event in enumerate(sampleTree):
                            for j, feature in enumerate(features):
                                inputData[i, j] = sampleTree.evaluate(feature)
                            # total weight comes from weightF (btag, lepton sf, ...) and treeScale to scale MC to x-section
                            eventWeight = sampleTree.evaluate(weightF)
                            specialWeight =  sampleTree.evaluate(sample.specialweight) if useSpecialWeight else 1.0 
                            totalWeight = treeScale * eventWeight * specialWeight 
                            weightLists[datasetName].append(totalWeight)
                            targetLists[datasetName].append(categories.index(category))
                            sampleIndexList[datasetName].append(sampleIndex)
                            sumOfWeights += totalWeight
                            nMCevents    += 1
                            
                            # add weights varied by (btag) systematics
                            #for syst in self.weightSYS:
                            #    weightListsSYS[syst][datasetName].append(treeScale * sampleTree.evaluate(self.weightSYSweights[syst]))
                            deltas = []
                            for syst in self.systematics:
                                delta_up   = sampleTree.evaluate(syst['U']) - eventWeight
                                delta_down = sampleTree.evaluate(syst['D']) - eventWeight
                                delta = 0.5 * (np.abs(delta_up) + np.abs(delta_down))
                                deltas.append(delta*delta)
                            totalDelta = np.sqrt(sum(deltas))

                            # convert to absolute error on total event weight
                            weightListSYStotal[datasetName].append(treeScale * totalDelta * specialWeight)

                            # fill systematics 
                            #for k, feature_s in features_sys.iteritems():
                            #    for j, feature in enumerate(feature_s):
                            #        inputData_sys[k][i,j] = sampleTree.evaluate(feature)
                        print("\x1b[43mINFO:", sample, ":", nMCevents, "MC events ->", sumOfWeights, "\x1b[0m")
                        eventNumberStat.append([sample, nMCevents, sumOfWeights])
                        arrayLists[datasetName].append(inputData)
                        #for sys in systematics:
                        #    arrayLists_sys[sys][datasetName].append(inputData_sys[sys])

                    else:
                        print ("\x1b[31mERROR: TREE NOT FOUND:", sample.name, " -> not cached??\x1b[0m")
                        raise Exception("CachedTreeMissing")

        if self.includeData:
            arrayListsData = []
            # all DATA 
            print("INFO: DATA=", self.dataSamples)
            for sample in self.dataSamples:
                print("INFO: add DATA sample:", sample.identifier)

                # cuts
                sampleCuts = [sample.subcut]
                # cut from the mva region
                if self.treeCut:
                    sampleCuts.append(self.treeCut)

                # get ROOT tree for selected sample & region cut
                tc = TreeCache.TreeCache(
                        sample=sample,
                        cutList=sampleCuts,
                        inputFolder=self.samplesPath,
                        config=self.config,
                        debug=True
                    )
                sampleTree = tc.getTree()
                if sampleTree:
                    # initialize numpy array
                    nSamples = sampleTree.GetEntries()
                    features = self.MVA_Vars['Nominal']
                    nFeatures = len(features) 
                    inputData = np.zeros((nSamples, nFeatures), dtype=np.float32)

                    # initialize formulas for ROOT tree
                    for feature in features:
                        sampleTree.addFormula(feature)

                    # fill numpy array from ROOT tree
                    for i, event in enumerate(sampleTree):
                        for j, feature in enumerate(features):
                            inputData[i, j] = sampleTree.evaluate(feature)

                    arrayListsData.append(inputData)

        ##systematics for training
        #puresystematics = deepcopy(systematics)
        #if 'Nominal' in puresystematics:
        #    puresystematics.remove('Nominal')
        puresystematics = [x['name'] for x in self.systematics]

        xSecs = {}
        SFs = {}
        for i,category in enumerate(categories):
            for j,sample in enumerate(self.samples[category]):
                if sample.name not in xSecs:
                    xSecs[sample.name] = sample.xsec
                if sample.name not in SFs:
                    SFs[sample.name] = sample.sf

        # concatenate all data from different samples
        self.data = {
                'train': {
                    'X': np.concatenate(arrayLists['train'], axis=0),
                    'y': np.array(targetLists['train'], dtype=np.float32),
                    'sample_weight': np.array(weightLists['train'], dtype=np.float32),
                    'sample_weight_error': np.array(weightListSYStotal['train'], dtype=np.float32),
                    'sample_index': np.array(sampleIndexList['train'], dtype=np.int32)
                    },
                'test': {
                    'X': np.concatenate(arrayLists['test'], axis=0), 
                    'y': np.array(targetLists['test'], dtype=np.float32), 
                    'sample_weight': np.array(weightLists['test'], dtype=np.float32),
                    'sample_weight_error': np.array(weightListSYStotal['test'], dtype=np.float32),
                    'sample_index': np.array(sampleIndexList['test'], dtype=np.int32)
                    },
                'category_labels': {idx: label for idx, label in enumerate(categories)},
                'meta': {
                    'version': self.dataFormatVersion,
                    'region': self.mvaName,
                    'cutName': self.treeCutName,
                    'cut': self.treeCut,
                    'trainCut': self.trainCut,
                    'testCut': self.evalCut,
                    'samples': self.sampleNames,
                    'xSecs': xSecs,
                    'scaleFactors': SFs,
                    'weightF': weightF,
                    'weightSYS': self.weightSYS,
                    'variables': ' '.join(self.MVA_Vars['Nominal']),
                    'systematics': puresystematics,
                    }
                }

        if self.includeData:
            self.data['data'] = {'X': np.concatenate(arrayListsData, axis=0)}

        ## add systematics variations
        #for sys in systematics:
        #    self.data['train']['X_'+sys] = np.concatenate(arrayLists_sys[sys]['train'], axis=0)
        #for syst in self.weightSYS:
        #    self.data['train']['sample_weight_'+syst] = np.array(weightListsSYS[syst]['train'], dtype=np.float32)

        if not os.path.exists("./dumps"):
            os.makedirs("dumps")
        baseName = './dumps/' + self.config.get('Configuration','channel') + self.config.get('General','dataset') + '_' + self.mvaName + '_' + datetime.datetime.now().strftime("%y%m%d")

        if config.has_option('MVAGeneral','ntupleVersion'):
            baseName += '_' + config.get('MVAGeneral','ntupleVersion')

        numpyOutputFileName = baseName + '.dmpz'
        hdf5OutputFileName = baseName + '.h5'
        print("INFO: saving output...")
        
        success = False
        try:
            if self.config.has_option(self.mvaName, 'writeNumpy') and eval(self.config.get(self.mvaName, 'writeNumpy')):
                self.saveAsPickledNumpy(numpyOutputFileName)
                success = True
        except Exception as e:
            print("ERROR: writing numpy array failed.", e)

        try:
            self.saveAsHDF5(hdf5OutputFileName)
            success = True
        except Exception as e:
            print("ERROR: writing HDF5 file failed.", e)

        eventNumberStat.sort(key = lambda x: x[0].identifier)
        print("MC stats:")
        for es in eventNumberStat:
            print(es[0].identifier.ljust(50), es[0].name.ljust(20), ("%d"%es[1]).ljust(10), ("%1.4f"%es[2]).ljust(10))
        print("*"*80)

        print("INFO: done." if success else "ERROR: no output file written")
        return success

    def saveAsPickledNumpy(self, outputFileName):
        with gzip.open(outputFileName, 'wb') as outputFile:
            pickle.dump(self.data, outputFile)
        print("written to:\x1b[34m", outputFileName, " \x1b[0m")

    def saveAsHDF5(self, outputFileName):
        # check arrays
        for k in ['train', 'test', 'data']:
            if k in self.data:
                lens = len(list(set([len(self.data[k][x]) for x in self.data[k].keys()])))
                if lens != 1:
                    print("ERROR: array length mismatch:", k, [[x,len(self.data[k][x])] for x in self.data[k].keys()])
                    raise Exception("LengthMismatch")

        # save to file
        f = h5py.File(outputFileName, 'w')
        for k in ['meta', 'category_labels']:
            f.attrs[k] = json.dumps(self.data[k].items())
        for k in ['train', 'test', 'data']:
            if k in self.data:
                for k2 in self.data[k].keys():
                    f.create_dataset(k + '/' + k2, data=self.data[k][k2], compression="gzip", compression_opts=9)
        f.close()
        print("written to:\x1b[34m", outputFileName, " \x1b[0m")


# read arguments
argv = sys.argv
parser = OptionParser()
parser.add_option("-C", "--config", dest="config", default=[], action="append",
                      help="configuration file")
parser.add_option("-T", "--tag", dest="tag", default='',
                      help="configuration tag")
parser.add_option("-t","--trainingRegions", dest="trainingRegions", default='',
                      help="cut region identifier")
parser.add_option("-S","--systematics", dest="systematics", default=0,
                      help="include systematics (0 for none, 1 for bdtVars, 2 for all (with btagWeights)")
parser.add_option("-x", "--test", dest="test", action="store_true", help="for debugging only!!!", default=False)
parser.add_option("-d", "--include-data", dest="include_data", action="store_true", help="include data in output file as additional dataset", default=False)
(opts, args) = parser.parse_args(argv)
if opts.config =="":
        opts.config = ["config"]

# Import after configure to get help message
from myutils import BetterConfigParser, mvainfo, ParseInfo

if len(opts.tag.strip()) > 1:
    config = BetterConfigParser()
    config.read("{tag}config/paths.ini".format(tag=opts.tag))
    configFiles = config.get("Configuration", "List").split(' ')
    opts.config = ["{tag}config/{file}".format(tag=opts.tag, file=x.strip()) for x in configFiles]
    print("reading config files:", opts.config)

sys = False
btagSys = False
if int(opts.systematics) > 0:
    sys = True
    if int(opts.systematics) > 1:
        btagSys = True
# load config
config = BetterConfigParser()
config.read(opts.config)
converter = SampleTreesToNumpyConverter(config, opts.trainingRegions, useSyst=sys, useWeightSyst=btagSys, testRun=opts.test, includeData=opts.include_data)
success = converter.run()
if not success:
    raise Exception("WriteTrainingDataFailed")
