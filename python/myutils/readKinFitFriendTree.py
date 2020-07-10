#!/usr/bin/env python
from __future__ import print_function
import ROOT
from BranchTools import Collection
from BranchTools import AddCollectionsModule
from kinFitterXbb import kinFitterXbb

# this module is to read the friend trees created by the kinematic fit when
# the systematics are split to different jobs (due to memory leaks in kinfit code)
class readKinFitFriendTree(AddCollectionsModule):
    def __init__(self, directory, name='Events_f', systematics='Nominal'):
        self.directory = directory
        self.systematics = [(x if not self._isnominal(x) else None) for x in systematics.strip().split(' ')]
        self.systematicsData = [None]
        self.name = name
        super(readKinFitFriendTree, self).__init__()
        self.first = True
        self.enabled = True
        self.hash = '%d'%hash(','.join([str(y) for y in self.systematics])) 

    def customInit(self, initVars):
        self.tree = initVars['tree']
        self.sampleTree = initVars['sampleTree']
        self.sample = initVars['sample']

        if self.sample.isData():
            self.systematics = [x for x in self.systematics if x in self.systematicsData or self._isnominal(x)]

        if len(self.systematics) < 1:
            self.enabled = False
            return

        # each tree needs a unique alias
        # copy files to scratch first to keep number of concurrent connections to storage element as low as possible
        self.sampleTree.addFriend(self.directory, alias="Events_" + self.hash, copy=True)

        # unfortunately friend trees are not copied by CloneTree and getattr() also doesn't work in pyROOT :-(
        # pretty ugly hack
        self.kinFitter = kinFitterXbb(year=2016, jetIdCut=2, puIdCut=6)
        self.kinFitter.systematics = self.systematics
        self.kinFitter.customInit(initVars)
        self.branches = self.kinFitter.branches
        self.branchBuffers = self.kinFitter.branchBuffers
        # this will disable the fitting itself but will keep the buffers for new branches, which will be filled by this module instead!
        self.kinFitter.disable()
        # /pretty ugly hack


    def processEvent(self, tree):
        if not self.hasBeenProcessed(tree) and self.enabled:
            self.markProcessed(tree)

            if self.first:
                self.first = False
                # another ugly hack to make friend trees usable in python
                # create a formula to make the leaves of the friend trees accessible
                # after all friend trees have been added!
                for k in self.branchBuffers.keys():
                    self.sampleTree.addFormula(k)

            for k in self.branchBuffers.keys():
                self._b(k)[0] = int(self.sampleTree.evaluate(k)) if type(self._b(k)[0]) == int else self.sampleTree.evaluate(k)

    def finish(self):
        for t in self.sampleTree.outputTrees:
            fl = t['tree'].GetListOfFriends()
            if fl:
                fl.Delete()

