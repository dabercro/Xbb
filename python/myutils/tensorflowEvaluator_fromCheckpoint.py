#!/usr/bin/env python
from __future__ import print_function
import ROOT
import numpy as np
from Jet import Jet
from BranchTools import Collection
from BranchTools import AddCollectionsModule
import sys
import os
import json
from myutils.XbbTools import XbbTools
from myutils.XbbTools import XbbMvaInputsList

# tfZllDNN repository has to be cloned inside the python folder
sys.path.append("..")
sys.path.append("../tfVHbbDNN/")
from tfVHbbDNN.tfDNNevaluator import TensorflowDNNEvaluator

# reads tensorflow checkpoints and applies DNN output to ntuples,
# including variations from systematic uncertainties
# needs: tensorflow >=1.4
class tensorflowEvaluator(AddCollectionsModule):

    def __init__(self, mvaName, nano=False, condition=None):
        self.mvaName = mvaName
        self.nano = nano
        self.debug = False
        self.condition = condition
        self.fixInputs = []
        super(tensorflowEvaluator, self).__init__()

    def customInit(self, initVars):
        self.config = initVars['config']
        self.sampleTree = initVars['sampleTree']
        self.sample = initVars['sample']

        if self.condition:
            self.sampleTree.addFormula(self.condition)

        self.hJidx = self.config.get('General', 'hJidx') if self.config.has_option('General', 'hJidx') else 'hJidx'

        self.scalerDump = self.config.get(self.mvaName, 'scalerDump') if self.config.has_option(self.mvaName, 'scalerDump') else None
        self.checkpoint = self.config.get(self.mvaName, 'checkpoint') if self.config.has_option(self.mvaName, 'checkpoint') else None
        if self.config.has_option(self.mvaName, 'branchName'):
            self.branchName = self.config.get(self.mvaName, 'branchName')
        elif self.checkpoint is not None:
            self.branchName = self.checkpoint.strip().replace('/model.ckpt','').replace('/','_')
            if self.branchName[0] in ['0','1','2','3','4','5','6','7','8','9']:
                self.branchName = 'DNN_' + self.branchName

        self.addDebugVariables = eval(self.config.get('Multi', 'evalAddDebugVariables')) if self.config.has_section('Multi') and self.config.has_option('Multi', 'evalAddDebugVariables') else False

        if self.checkpoint is None:
            print("\x1b[31mERROR: 'checkpoint' option missing for MVA config section [%s]! .../model.ckpt has to be specified to be able to restore classifier.\x1b[0m"%self.mvaName)
            raise Exception("CheckpointError")

        if self.scalerDump is not None and not os.path.isfile(self.scalerDump):
            self.scalerDump = None

        if self.config.has_option(self.mvaName, 'fixInputs'):
            self.fixInputs = eval(self.config.get(self.mvaName, 'fixInputs'))

        if os.path.isdir(self.checkpoint):
            self.checkpoint += '/model.ckpt'

        if not os.path.isfile(self.checkpoint + '.meta'):
            print("\x1b[31mERROR: can't restore from graph! .meta file not found in checkpoint:", self.checkpoint, "\x1b[0m")
            raise Exception("CheckpointError")

        # INFO file (with training parameters)
        if os.path.isfile(self.checkpoint + '.info'):
            with open(self.checkpoint + '.info', 'r') as infoFile:
                self.info = json.load(infoFile)
        else:
            print("WARNING: (optional) .info file not found in checkpoint!")
            self.info = {}

        # CLASSES
        self.classes = eval(self.config.get(self.mvaName, 'classes')) if self.config.has_option(self.mvaName, 'classes') else None
        if self.classes:
            self.nClasses = len(self.classes)
        else:
            try:
                self.nClasses = eval(self.config.get(self.mvaName, 'nClasses'))
            except Exception as e:
                self.nClasses = len(self.info["labels"].keys())
        if self.nClasses < 3:
            self.nClasses = 1
        if self.nClasses > 1:
            print("INFO: multi-class checkpoint found! number of classes =", self.nClasses)
        else:
            print("INFO: binary-classifier found!")

        # FEATURES
        self.featuresConfig = None
        self.featuresCheckpoint = None
        try:
            self.featuresConfig = self.config.get(self.config.get(self.mvaName, "treeVarSet"), "Nominal").strip().split(" ")
        except Exception as e:
            print("WARNING: could not get treeVarSet from config:", e)
        if 'variables' in self.info:
            self.featuresCheckpoint = self.info['variables']

        if self.featuresConfig is None and self.featuresCheckpoint is not None:
            self.features = self.featuresCheckpoint
        elif self.featuresConfig is not None and self.featuresCheckpoint is None:
            self.features = self.featuresConfig
        elif self.featuresConfig is None and self.featuresCheckpoint is None:
            raise Exception("NoInputFeaturesDefined")
        else:
            self.features = self.featuresCheckpoint
            if len(self.featuresConfig) != len(self.featuresCheckpoint):
                print("\x1b[31mWARNING: number of input features does not match!")
                print(" > classifier expects:", len(self.featuresCheckpoint))
                print(" > configuration has:", len(self.featuresConfig), "\x1b[0m")
                print("INFO: => feature list from checkpoint will be used.")
            else:
                if self.config.has_option(self.mvaName, 'forceInputFeaturesFromConfig') and eval(self.config.get(self.mvaName, 'forceInputFeaturesFromConfig')):
                    print("INFO: forceInputFeaturesFromConfig is enabled, features from configuration will be used")
                    self.features = self.featuresConfig

            print("INFO: list of input features:")
            print("INFO:", "config".ljust(40), "---->", "checkpoint")
            match = True
            for i in range(min(len(self.featuresConfig),len(self.featuresCheckpoint))):
                if self.featuresConfig[i] != self.featuresCheckpoint[i]:
                    print("\x1b[41m\x1b[37mINFO:", self.featuresConfig[i].ljust(40), "---->", self.featuresCheckpoint[i].ljust(40), " => MISMATCH!\x1b[0m")
                    match = False
                else:
                    print("INFO:\x1b[32m", self.featuresConfig[i].ljust(40), "---->", self.featuresCheckpoint[i], "(match)\x1b[0m")
            if match:
                print("INFO: => all input variables match!")
            else:
                print("INFO: some variables are not identically defined as for the training, please check!")
                if self.config.has_option(self.mvaName, 'forceInputFeaturesFromConfig') and eval(self.config.get(self.mvaName, 'forceInputFeaturesFromConfig')):
                    print("\x1b[31mWARNING: forceInputFeaturesFromConfig is enabled, features from configuration will be used although they could be incompatible with the features used during training.\x1b[0m")

            # fix input features at constant values: check if features given exist in checkpoint
            if self.fixInputs:
                for feature,value in self.fixInputs.items():
                    if feature not in self.featuresCheckpoint:
                        print("ERROR: can't fix input feature '", feature, "' to value ", value," => feature not found in checkpoint.")
                        raise Exception("ConfigError")

        self.featureList = XbbMvaInputsList(self.features, config=self.config)
        self.nFeatures   = self.featureList.length()

        # SIGNAL definition
        self.signalIndex = 0
        if self.config.has_option(self.mvaName, 'signalIndex'):
            self.signalIndex = eval(self.config.get(self.mvaName, 'signalIndex'))
        self.signalClassIds = [self.signalIndex]
        if self.classes:
            self.signalClassIds = [x for x,y in enumerate(self.classes) if y[0].startswith('SIG')]
            print("INFO: signals:", self.signalClassIds)

        print("INFO: number of classes:", self.nClasses if self.nClasses > 1 else 2)

        # systematics
        self.systematics = self.config.get(self.mvaName, 'systematics').split(' ') if self.config.has_option(self.mvaName, 'systematics') else self.config.get('systematics', 'systematics').split(' ')

        # create output branches
        self.dnnCollections = []
        for i in range(self.nClasses):
            collectionName = self.branchName if self.nClasses==1 else self.branchName + "_%d"%i
            if self.nClasses==1 or self.addDebugVariables:
                self.dnnCollection = Collection(collectionName, self.systematics, leaves=True)
                self.addCollection(self.dnnCollection)
                self.dnnCollections.append(self.dnnCollection)

        # create formulas for input variables
        self.inputVariables = {}
        for syst in self.systematics:
            systBase, UD              = XbbTools.splitSystVariation(syst, sample=self.sample)
            self.inputVariables[syst] = [XbbTools.sanitizeExpression(self.featureList.get(i, syst=systBase, UD=UD), self.config, debug=self.debug) for i in range(self.nFeatures)]
            for var in self.inputVariables[syst]:
                self.sampleTree.addFormula(var)

        # additional pre-computed values for multi-classifiers
        if self.nClasses > 1:
            self.dnnCollectionsMulti = {
                        'default': Collection(self.branchName,   self.systematics, leaves=True),
                    }
            if self.addDebugVariables:
                self.dnnCollectionsMulti.update({
                            'argmax' : Collection(self.branchName + "_argmax", self.systematics, leaves=True),
                            'max'    : Collection(self.branchName + "_max",    self.systematics, leaves=True),
                            'max2'   : Collection(self.branchName + "_max2",   self.systematics, leaves=True),
                            'signal' : Collection(self.branchName + "_signal", self.systematics, leaves=True),
                        })

            for k,v in self.dnnCollectionsMulti.items():
                self.addCollection(v)

        # create tensorflow graph
        self.ev = TensorflowDNNEvaluator(checkpoint=self.checkpoint, scaler=self.scalerDump)

    # called from main loop for every event
    def processEvent(self, tree):

        if not self.hasBeenProcessed(tree):
            self.markProcessed(tree)

            if self.condition is None or self.sampleTree.evaluate(self.condition):
                # fill input variables
                inputs = np.full((len(self.systematics), self.nFeatures), 0.0, dtype=np.float32)
                for j, syst in enumerate(self.systematics):
                    for i, var in enumerate(self.inputVariables[syst]):
                        if var in self.fixInputs:
                            inputs[j,i] = self.fixInputs[var]
                        else:
                            inputs[j,i] = self.sampleTree.evaluate(var)

                # use TensorflowDNNEvaluator
                probabilities = self.ev.eval(inputs)
            else:
                probabilities = np.full((len(self.systematics), self.nFeatures), -1.0, dtype=np.float32)

            # fill output branches
            if len(self.dnnCollections) == 1:
                # for signal/background store just one of the two, selected by self.signalIndex
                for i,dnnCollection in enumerate(self.dnnCollections):
                    for j, syst in enumerate(self.systematics):
                        dnnCollection[dnnCollection.name][j] = probabilities[j, self.signalIndex]
            else:
                if self.addDebugVariables:
                    # for multi-class, store all
                    for i,dnnCollection in enumerate(self.dnnCollections):
                        for j, syst in enumerate(self.systematics):
                            dnnCollection[dnnCollection.name][j] = probabilities[j, i]

                for j, syst in enumerate(self.systematics):
                    sortedValues = np.sort(probabilities[j])
                    argmaxP = np.argmax(probabilities[j])
                    maxP = min(sortedValues[-1], 0.9999)
                    if self.addDebugVariables:
                        self.dnnCollectionsMulti['argmax']._b()[j] = argmaxP
                        self.dnnCollectionsMulti['max']._b()[j]    = maxP
                        self.dnnCollectionsMulti['max2']._b()[j]   = min(sortedValues[-2], 0.9999)

                        # sum the total signal probability, useful in case of multiple signals
                        self.dnnCollectionsMulti['signal']._b()[j] = min(sum([probabilities[j, i] for i in self.signalClassIds]), 0.9999)

                    # default linearization method
                    self.dnnCollectionsMulti['default']._b()[j] = argmaxP + maxP

        return True

    def cleanUp(self):
        pass

