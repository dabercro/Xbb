#!/usr/bin/env python
from __future__ import print_function
import BetterConfigParser
import os
import inspect
import ROOT
from itertools import ifilter
from copytreePSI import filelist
from FileLocator import FileLocator
from sample_parser import ParseInfo
from XbbTools import XbbTools

# helper class to read full config object
# use:
#      config = XbbConfigReader.read('Zll2017')

class XbbConfigReader(object):

    def __init__(self):
        self.config = None

    @staticmethod
    def read(configTag):

        configDirectory = configTag + 'config/'
        debug = 'XBBDEBUG' in os.environ
        if debug:
            print("DEBUG: read configuration \x1b[35m", configTag, "\x1b[0m")

        # paths.ini contains list of config files to use, relative to configDirectory
        if configTag.endswith('.ini'):
            config = BetterConfigParser.BetterConfigParser()
            config.read(configTag)
            config.set('Configuration','__temporary_config', "{tag}.volatile.ini".format(tag=configTag))
        else:
            pathconfig = BetterConfigParser.BetterConfigParser()
            pathconfig.read(configDirectory + '/paths.ini')
            configFiles = [x.strip() for x in pathconfig.get('Configuration', 'List').split(' ') if x.strip() != 'volatile.ini']

            # read actual config
            config = BetterConfigParser.BetterConfigParser()
            for configFile in configFiles:
                if debug:
                    print("DEBUG: --> read configFile:", configFile)
                config.read(configDirectory + configFile)
            config.set('Configuration','__temporary_config', "{tag}config/volatile.ini".format(tag=configTag))
        if debug:
            print('DEBUG: \x1b[35m read', len(config.sections()), 'sections\x1b[0m')
        if not config.has_option('configuration','__self'):
            config.set('Configuration', '__self', configTag)

        return config

class XbbConfigTools(object):

    def __init__(self, config):
        self.config = config
        self.fileLocator = None
        self.samplesInfo = None

    def initFS(self, force=False): 
        if self.fileLocator is None or force:
            self.fileLocator = FileLocator(config=self.config)

    def fs(self):
        if self.fileLocator is None:
            self.initFS()
        return self.fileLocator

    def loadNamespaces(self):
        #default
        try:
            defaultNamespace = self.get('VHbbNameSpace','library')
            ROOT.gSystem.Load(defaultNamespace)
        except Exception as e:
            print(e)

    # list of DATA sample names
    def getData(self):
        return eval(self.config.get('Plot_general', 'Data'))

    # list of MC sample names
    def getMC(self):
        return eval(self.config.get('Plot_general', 'samples'))

    def getSamplesInfo(self):
        if self.samplesInfo is None:
            self.samplesInfo = ParseInfo(config=self.config)
        return self.samplesInfo

    # processed sample identifiers (may not be actually USED)
    def getSampleIdentifiers(self, filterList=None):
        s = self.getSamplesInfo().getSampleIdentifiers()
        if filterList is not None:
            s = XbbTools.filterSampleList(s, filterList)
        s.sort()
        return s

    # list of all sample names (data + mc)
    def getUsedSamples(self):
        return self.getMC() + self.getData()

    # get list of original file names: /store/...
    def getOriginalFileNames(self, sampleIdentifier):
        return filelist(self.config.get('Directories', 'samplefiles'), sampleIdentifier)

    # get list of file names (e.g. in SYSout folder)
    def getFileNames(self, sampleIdentifier, folder='SYSout'):
        self.initFS()
        try:
            originalFileNames = self.getOriginalFileNames(sampleIdentifier)
        except:
            originalFileNames = []
        samplePath = self.config.get('Directories', folder)
        fileNames = ["{path}/{subfolder}/{filename}".format(path=samplePath, subfolder=sampleIdentifier, filename=self.fileLocator.getFilenameAfterPrep(x)) for x in originalFileNames]
        return fileNames

    def parseCommaSeparatedList(self, listAsString):
        return [x.strip() for x in listAsString.split(',') if len(x.strip()) > 0]

    def parseSpaceSeparatedList(self, listAsString):
        return [x.strip() for x in listAsString.split(' ') if len(x.strip()) > 0]

    def getPlotRegions(self):
        return self.parseCommaSeparatedList(self.get('Plot_general', 'List'))
    
    def getDatacardRegions(self):
        return self.parseCommaSeparatedList(self.get('LimitGeneral', 'List'))

    def getTrainingRegions(self):
        return self.parseCommaSeparatedList(self.get('MVALists', 'List_for_submitscript'))

    def getPlotRegionCutName(self, plotRegion):
        configSection = 'Plot:' + plotRegion
        if self.has_option(configSection, 'Cut'):
            return self.get(configSection, 'Cut')
        else:
            return plotRegion 

    def getDatacardCutName(self, datacardRegion):
        configSection = 'dc:' + datacardRegion
        if self.has_option(configSection, 'Cut'):
            return self.get(configSection, 'Cut')
        else:
            return datacardRegion

    def getTrainingRegionCutName(self, trainingRegion):
        configSection = trainingRegion
        if self.has_option(configSection, 'Cut'):
            return self.get(configSection, 'Cut')
        elif self.has_option(configSection, 'treeCut'):
            return self.get(configSection, 'treeCut')
        else:
            return trainingRegion

    def getTrainingRegionVarSet(self, trainingRegion):
        return self.get(trainingRegion, 'treeVarSet')

    def getTrainingRegionVariables(self, trainingRegion):
        treeVarSet = self.getTrainingRegionVarSet(trainingRegion)
        return self.parseSpaceSeparatedList(self.get(treeVarSet, 'Nominal'))

    def getDatacardRegionType(self, datacardRegion):
        configSection = 'dc:' + datacardRegion
        if self.has_option(configSection, 'type'):
            datacardType = self.get(configSection, 'type')
            if datacardType.lower() not in ['bdt', 'cr', 'dnn', 'mjj']:
                print("ERROR: unknown datacard type:", datacardType)
                raise Exception("DatacardTypeUnknown")
            return datacardType
        else:
            raise Exception("DatacardTypeUndefined")
        return None
    
    def getDatacardRegionSignals(self, datacardRegion):
        configSection = 'dc:' + datacardRegion
        if self.has_option(configSection, 'signal'):
            return eval(self.get(configSection, 'signal'))
        else:
            raise Exception("DatacardTypeUndefined")
        return None
    
    def getDatacardRegionBackgrounds(self, datacardRegion):
        configSection = 'dc:' + datacardRegion
        if self.has_option(configSection, 'background'):
            return eval(self.get(configSection, 'background'))
        else:
            raise Exception("DatacardTypeUndefined")
        return None

    def getPlotVariables(self):
        return self.parseCommaSeparatedList(self.get('Plot_general', 'var'))

    def getPlotVariableDefinition(self, plotVariableName):
        configSection = 'plotDef:' + plotVariableName
        if self.has_section(configSection) and self.has_option(configSection, 'relPath'):
            return self.get(configSection, 'relPath')
        else:
            return None

    def getCutString(self, cutName):
        return self.get('Cuts', cutName)

    def sections(self):
        return self.config.sections()

    def has_section(self, section):
        return self.config.has_section(section)

    def has_option(self, section, option):
        return self.config.has_section(section) and self.config.has_option(section, option)

    def get(self, section, option, default=None):
        if default is not None:
            if self.has_option(section, option):
                return self.config.get(section, option)
            else:
                return default
        else:
            return self.config.get(section, option)

    def getPlotVariableSections(self):
        return [x for x in self.config.sections() if x.startswith('plotDef:')]

    def getJECuncertainties(self, step=None):
        if step is None:
            configOption = 'JEC'
        else:
            configOption = 'JEC_' + step
        if self.config.has_option('systematics', configOption):
            systematics = eval(self.config.get('systematics', configOption))
        elif  self.config.has_option('systematics', 'JEC'):
            systematics = eval(self.config.get('systematics', 'JEC'))
        else:
            # default
            #systematics = ['jer','jerReg','jesAbsoluteStat','jesAbsoluteScale','jesAbsoluteFlavMap','jesAbsoluteMPFBias','jesFragmentation','jesSinglePionECAL','jesSinglePionHCAL','jesFlavorQCD','jesRelativeJEREC1','jesRelativeJEREC2','jesRelativeJERHF','jesRelativePtBB','jesRelativePtEC1','jesRelativePtEC2','jesRelativePtHF','jesRelativeBal','jesRelativeFSR','jesRelativeStatFSR','jesRelativeStatEC','jesRelativeStatHF','jesPileUpDataMC','jesPileUpPtRef','jesPileUpPtBB','jesPileUpPtEC1','jesPileUpPtEC2','jesPileUpPtHF','jesPileUpMuZero','jesPileUpEnvelope','jesTotal']
            raise Exception("ConfigError: Specify the JEC list in [systematics]") 
        systematics = list(set(systematics))
        systematics.sort()
        return systematics

    def setList(self, setOptions):
        # escaping of semicolon
        semicolonEscapeSequence = '##SEMICOLON##'
        setOptions = setOptions.replace('\;', semicolonEscapeSequence) 
        prevSection = None
        for optValue in setOptions.split(';'):
            optValue = optValue.replace(semicolonEscapeSequence, ';').strip()
            syntaxOk = True
            try:
                if ':=' in optValue:
                    opt = optValue.split(':=')[0]
                    value = optValue.split(':=')[1]
                elif '=' in optValue:
                    splitParts = optValue.split('=')
                    if len(splitParts) > 2:
                        print("\x1b[31mWARNING: more than one equal sign found in expression, split at the first one! use ':=' to force split at another position!\x1b[0m")
                    opt = optValue.split('=')[0]
                    value = '='.join(optValue.split('=')[1:])
                elif optValue:
                    opt = optValue.split(':')[0]
                    value = optValue.split(':')[1]
            except Exception as e:
                print("ERROR:",e)
                print("ERROR: syntax error in:", optValue)
                print("ERROR: use ; to separate options and use \; to escape semicolons in case they are inside the value. Use := for assignment.")
                syntaxOk = False
                raise

            if syntaxOk:

                configSection = opt.split('.')[0]
                configOption  = opt.split('.')[1]

                if len(configSection.strip()) < 1:
                    if prevSection is None:
                        raise Exception("ConfigSetError")
                    else:
                        configSection = prevSection
                
                prevSection = configSection
                if not self.config.has_section(configSection):
                    self.config.add_section(configSection)
                if self.config.has_section(configSection) and self.config.has_option(configSection, configOption):
                    print("\x1b[31mCONFIG: SET", "{s}.{o}".format(s=configSection, o=configOption), "=", value, "\x1b[0m")
                else:
                    print("\x1b[31mCONFIG: ADD", "{s}.{o}".format(s=configSection, o=configOption), "=", value, "\x1b[0m")
                self.config.set(configSection, configOption, value)

    def formatSampleName(self, sampleIdentifier, maxlen=80, padding=False):
        if len(sampleIdentifier) > maxlen:
            s = sampleIdentifier[:maxlen-7] + '...' + sampleIdentifier[-4:]
        else:
            s = sampleIdentifier
        if padding:
            s = s.ljust(maxlen)
        return s

class XbbConfigChecker(object):
    def __init__(self, config):
        self.config = config
        self.errors = []

    def addError(self, category, message):
        self.errors.append([category, message])

    def printErrors(self):
        for e in self.errors:
            print(("\x1b[31m[%s]"%e[0]).ljust(16), e[1], "\x1b[0m")
        print("-"*80)
        if len(self.errors) > 0:
            print("\x1b[31m%d errors in total.\x1b[0m"%len(self.errors))
        else:
            print("%d errors in total."%len(self.errors))

    def getStatus(self):
        status = 0
        if len(self.errors) > 0:
            status = 1
        return status

    def checkPlotRegions(self):
        plotRegions = self.config.getPlotRegions()
        for plotRegion in plotRegions:
            try:
                cutName = self.config.getPlotRegionCutName(plotRegion)
                cutString = self.config.getCutString(cutName)

                print("plot region:", plotRegion)
                print("  ->", cutName)
                print("    ->", cutString)

                if cutString.count('(') != cutString.count(')'):
                    raise Exception("CutStringUnbalancedRoundBrackets")
                if cutString.count('[') != cutString.count(']'):
                    raise Exception("CutStringUnbalancedSquareBrackets")

            except Exception as e:
                self.addError('Plot region', plotRegion + ' ' + repr(e))

    def checkPlotVariables(self):
        plotVariables = self.config.getPlotVariables()
        for plotVariable in plotVariables:
            plotVariableDefinition = self.config.getPlotVariableDefinition(plotVariable)
            if plotVariableDefinition is None:
                self.addError('Plot variables', 'Variable not found: %s'%plotVariable)

    def checkDatacardRegions(self):
        datacardRegions      = self.config.getDatacardRegions()
        samples              = ParseInfo(config=self.config)
        availableSampleNames = [x.name for x in samples]
        for datacardRegion in datacardRegions:
            try:
                cutName = self.config.getDatacardCutName(datacardRegion)
                cutString = self.config.getCutString(cutName)
                print("datacard region:", datacardRegion)
                print("  ->", cutName)
                print("    ->", cutString)

                regionType = self.config.getDatacardRegionType(datacardRegion)
                print("  -> TYPE:", regionType)

                signals = self.config.getDatacardRegionSignals(datacardRegion)
                backgrounds = self.config.getDatacardRegionBackgrounds(datacardRegion)
                print("  -> SIG:\x1b[34m", signals, "\x1b[0m")
                print("  -> BKG:\x1b[35m", backgrounds, "\x1b[0m")

                for x in list(signals) + list(backgrounds):
                    if x not in availableSampleNames:
                         print("ERROR: not found: sample for datacard:", x)
                         raise Exception("SampleNotFound")

            except Exception as e:
                self.addError('datacard region', datacardRegion + ' ' + repr(e))

    def checkTrainingRegions(self):
        trainingRegions = self.config.getTrainingRegions()
        for trainingRegion in trainingRegions:
            print("training region:", trainingRegion)
            try:
                cutName = self.config.getTrainingRegionCutName(trainingRegion)
                cutString = self.config.getCutString(cutName)
                print("  ->", cutName)
                print("    ->", cutString)
                treeVarSet = self.config.getTrainingRegionVarSet(trainingRegion)
                print("  -> VARSET:\x1b[35m", treeVarSet, "\x1b[0m")
                variables = self.config.getTrainingRegionVariables(trainingRegion)
                variablesList = " ".join([(("\x1b[34m"+x+"\x1b[0m") if (i % 2) == 0 else ("\x1b[32m"+x+"\x1b[0m")) for i,x in enumerate(variables)])
                print("  -> VARS:", variablesList)

            except Exception as e:
                self.addError('Training region', trainingRegion + ' ' + repr(e))

    def checkSamples(self):
        samples         = ParseInfo(config=self.config)
        availableSampleNames = [x.name for x in samples]
        usedSampleNames = self.config.getUsedSamples()
        for sampleName in usedSampleNames:
            print(sampleName)
            if sampleName not in availableSampleNames:
                print("ERROR: not found sample:", sampleName)
                raise Exception("SampleNotFound")

    def checkConfigFiles(self):
        configFiles = [x.strip() for x in self.config.get('Configuration', 'List').split(' ') if x.strip() != 'volatile.ini' and len(x.strip()) > 0]
        for configFile in configFiles:
            filePath = self.config.get('Configuration','__self') + 'config/' + configFile
            print("  -> config file:", filePath)
            if not os.path.isfile(filePath):
                raise Exception("FileNotFound:"+filePath)

    # runs all methods starting with 'check'
    def checkAll(self):
        for f in ifilter(inspect.ismethod, [getattr(self, name) for name in dir(self) if name.startswith('check') and name != 'checkAll']):
            f()

