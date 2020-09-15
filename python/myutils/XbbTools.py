#!/usr/bin/env python
from __future__ import print_function
import json
import re
import fnmatch
import importlib
import ROOT
import array
import numpy as np

class XbbTools(object):

    @staticmethod
    def readJson(fileName):
        with open(fileName, 'r') as fp:
            res = json.load(fp)
        return res

    @staticmethod
    def splitParts(s):
        p_raw = [x for x in re.split('([+\-=/*,\(\)><^\|\&\!])', s) if len(x) > 0]
        p     = []
        ac    = []
        for x in p_raw:
            ac.append(x)
            t = ''.join(ac)
            if t.count('[')==t.count(']'):
                p.append(t)
                ac=[]
        if len(ac) > 0:
            raise Exception("UnbalancedBrackets")
        return p

    @staticmethod
    def joinParts(p):
        return ''.join(p)

    @staticmethod
    def parseReplacementRules(replacementRulesList):
        return [ x.split('>') for x in replacementRulesList ]

    @staticmethod
    def applyReplacementRules(nominalVariable, replacementRules):
        result = nominalVariable
        replacements = []
        # replace all matching parts by placeholders first to avoid multiple replacements of a single expression
        for x in replacementRules:
            placeholder = "###%d###"%len(replacements)
            newResult = result.replace(x[0], placeholder)
            if newResult != result:
                replacements.append([placeholder, x[1]])
                result = newResult
        # fill the placeholders
        for replacement in replacements:
            result = result.replace(replacement[0], replacement[1])
        return result

    @staticmethod
    def getSystematicsVariationTemplate(nominalVariable, replacementsList):
        stringParts         = XbbTools.splitParts(nominalVariable)
        replacementRules    = XbbTools.parseReplacementRules(replacementsList)
        replacedStringParts = [XbbTools.applyReplacementRules(x, replacementRules) for x in stringParts]
        replacedVariable    = XbbTools.joinParts(replacedStringParts)
        return replacedVariable

    @staticmethod
    def getReplacementRulesList(config, syst):
        replacementRulesDict = eval(config.get('LimitGeneral', 'sys_cut_suffix'))
        if syst in replacementRulesDict:
            if isinstance(replacementRulesDict[syst], list):
                return replacementRulesDict[syst]
            else:
                return [replacementRulesDict[syst]]
        else:
            if isinstance(replacementRulesDict['default'], list):
                return replacementRulesDict['default']
            else:
                return [replacementRulesDict['default']]

    @staticmethod
    def filterSampleList(sampleIdentifiers, samplesList):
        if samplesList and len([x for x in samplesList if x]) > 0:
            filteredList = []
            for expr in samplesList:
                if expr in sampleIdentifiers:
                    filteredList.append(expr)
                elif '*' in expr and not expr.startswith('!'):
                    for sampleIdentifier in sampleIdentifiers:
                        if fnmatch.fnmatch(sampleIdentifier, expr):
                            filteredList.append(sampleIdentifier)
                elif '*' in expr and expr.startswith('!'):
                    expr = expr[1:]
                    newList = []
                    for x in filteredList:
                        if not fnmatch.fnmatch(x, expr):
                            newList.append(x)
                    filteredList = newList
                elif expr.startswith('!'):
                    expr = expr[1:]
                    filteredList = [x for x in filteredList if x!=expr]

            filteredList = list(set(filteredList))
            return filteredList
        else:
            return sampleIdentifiers

    # this takes all elements from inputFilters and if they contain a *, replaces them with all matching elements from inputList, removing duplicates
    @staticmethod
    def filterList(inputList, inputFilters):
        if not type(inputFilters) == list:
            inputFilters = [inputFilters]
        matchingWildscard  = [x for x in inputList if any([fnmatch.fnmatch(x, y) for y in inputFilters])]
        nonWildcardFilters = [x for x in inputFilters if '*' not in x]
        return list(set(matchingWildscard + nonWildcardFilters))

    @staticmethod
    def filterListStrict(inputList, inputFilters):
        if not type(inputFilters) == list:
            inputFilters = [inputFilters]
        return [x for x in inputList if any([fnmatch.fnmatch(x, y) for y in inputFilters])]


    @staticmethod
    def parseSamplesList(samplesListString):
        return XbbTools.parseList(samplesListString, separator=',')

    @staticmethod
    def parseList(listString, separator=' '):
        return [x.strip() for x in listString.strip().split(separator) if len(x.strip()) > 0]

    @staticmethod
    def getSampleTreeFileNames(pathIN, folderName):
        filenames = open(pathIN+'/'+folderName+'.txt').readlines()

        ## search the folder containing the input files
        inputFiles = []

        for filename_ in filenames:
            if '.root' in filename_ :
                inputFiles.append(filename_.rstrip('\n'))

        return inputFiles

    @staticmethod
    def splitSystVariation(syst, sample=None):
        if syst.lower() == 'nominal' or (sample is not None and sample.isData()):
            systBase = None
            UD = None
        else:
            systBase = '_'.join(syst.split('_')[:-1])
            UD = syst.split('_')[-1]
        return systBase, UD

    @staticmethod
    def sanitizeExpression(expression, config, debug=False):
        if config.has_option('General', 'sanitizeExpression'):
            rules = eval(config.get('General', 'sanitizeExpression'))
            newExpression = expression
            for rule in rules:
                newExpression = newExpression.replace(rule[0], rule[1])
            if debug and newExpression != expression:
                print("DEBUG: sanitized expression:")
                print("DEBUG: from:", expression)
                print("DEBUG:   to:", newExpression)
            return newExpression
        else:
            return expression

    @staticmethod
    def getMvaSystematics(mvaName, config):
        return XbbTools.parseList(config.get(mvaName, 'systematics') if (config.has_section(mvaName) and config.has_option(mvaName, 'systematics')) else config.get('systematics', 'systematics'), separator=' ')

    @staticmethod
    def getModuleInfo(module, config):
        modulesInfo = []
        try:
            if '.' in module:
                section = module.split('.')[0]
                key = module.split('.')[1]
                if config.has_section(section) and config.has_option(section, key):
                    pyCode = config.get(section, key)
                    if pyCode.strip().startswith('['):
                        modulesList = eval(pyCode)
                        for submodule in modulesList:
                            modulesInfo += XbbTools.getModuleInfo(submodule, config=config)
                    else:
                        # import module from myutils
                        moduleName = pyCode.split('(')[0].split('.')[0].strip()
                        globals()[moduleName] = importlib.import_module(".{module}".format(module=moduleName), package="myutils")

                        # get object
                        wObject = eval(pyCode)
                        version = wObject.getVersion() if hasattr(wObject, "getVersion") else 0

                        modulesInfo.append([wObject, moduleName, version])
            else:
                modulesInfo.append([None, "unknown", -1])
        except Exception as e:
            print("\x1b[31mEXCEPTION: An exception occured while getting the module meta information (version etc) this might be due to an error in the module or config!\x1b[0m")
            print("\x1b[31mEXCEPTION:", e, "\x1b[0m")
            modulesInfo.append([None, "unknown", -2])
        return modulesInfo

    @staticmethod
    def readDictFromRootFile(rootFileName, treeName, key, value):
        dictResult = {}
        infoFile = ROOT.TFile.Open(rootFileName, "READ")
        infoTree = infoFile.Get(treeName)
        for ev in infoTree:
            dictResult["%s"%getattr(ev,key)] = getattr(ev,value)
        infoFile.Close()
        return dictResult

    @staticmethod
    def writeDictToRootFile(rootFileName, treeName, treeDescription, key, value, dictToWrite, fileLocator):
        pathOUT = '/'.join(rootFileName.split('/')[:-1])
        if not fileLocator.exists(pathOUT):
            fileLocator.makedirs(pathOUT)
        if fileLocator.exists(rootFileName):
            fileLocator.rm(rootFileName)
        infoFile = ROOT.TFile.Open(rootFileName, "NEW")
        infoTree = ROOT.TTree(treeName, treeDescription)
        infoTree.SetDirectory(infoFile)
        key_s = ROOT.std.string()
        value_i = array.array( 'i', [ 0 ])
        infoTree.Branch(key, key_s)
        infoTree.Branch(value, value_i, value + "/I")
        for k,v in dictToWrite.items():
            key_s.replace(0, ROOT.std.string.npos, k)
            value_i[0] = v
            infoTree.Fill()
        infoFile.Write()
        infoFile.Close()

    @staticmethod
    def resolvePlotVariable(varName, config):
        if config.has_section('plotDef:' + varName):
            return varName
        else:
            for section in config.sections():
                if section.startswith('plotDef:'):
                    if config.has_option(section, 'relPath') and config.get(section, 'relPath').strip() == varName:
                        return section.split('plotDef:')[1]
            print("\x1b[31mERROR: plot variable", varName, " could not be resolved!\x1b[0m")
            return varName

    @staticmethod
    def getPlotVariableBins(varName, config):
        plotDefSection = 'plotDef:' + XbbTools.resolvePlotVariable(varName, config)
        if config.has_section(plotDefSection):
            if config.has_option(plotDefSection, "binList"):
                return np.array(eval(config.get(plotDefSection, "binList")))
            else:
                binDef = {x:eval(config.get(plotDefSection, x)) for x in ['nBins','min','max']}
                return np.linspace(binDef['min'], binDef['max'], binDef['nBins']+1)

    @staticmethod
    def uniformBins(nBins, xMin, xMax):
        return np.linspace(xMin, xMax, nBins+1)

class XbbMvaInputsList(object):

    def __init__(self, expressions, config=None):
        self.config = config
        self.expressions = expressions

    def get(self, i, syst=None, UD=None):
        nominal = self.expressions[i]
        if syst is None:
            return nominal
        else:
            replacementRulesList = XbbTools.getReplacementRulesList(self.config, syst)
            sysVar = XbbTools.getSystematicsVariationTemplate(nominal, replacementRulesList)
            sysVar = sysVar.replace('{syst}', syst).replace('{UD}', UD)
            return sysVar

    def length(self):
        return len(self.expressions)
