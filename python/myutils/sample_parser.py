import os, sys, warnings
from copy import copy
from optparse import OptionParser
from BetterConfigParser import BetterConfigParser
from samplesclass import Sample
import fnmatch
import traceback

class ParseInfo:
    '''Class containing a list of Sample. Is filled during the prep stage.'''
    def __init__(self,samples_config=None,samples_path=None,config=None):
        '''
        Methode filling a list of Sample "self._samplelist = []" contained in the class. 
        "sample_path" contains the path where the samples are stored (PREPin). 
        "samples_config" is the "samples_nosplit.cfg" file. Depending of the variable "run_on_files" defined in "samples_nosplit.cfg", 
        the sample list are generated from the input folder (PREPin) or the list in "samples_nosplit.cfg" '''

        self.debug = 'XBBDEBUG' in os.environ

        # legacy behavior not supported anymore (INSTEAD: pass config object)
        if not config:
            print("\x1b[41m\x1b[37mERROR: ParseInfo() is now required to have config= argument passed\x1b[0m")
            traceback.print_exc(file=sys.stdout)
            raise Exception("Error")
        elif self.debug:
            print "DEBUG: config object passed to sample parser."

        lumi = float(config.get('General','lumi'))

        self._samplelist = []
        self.__fileslist = []

        configSamples = [x for x in config.sections() if config.has_option(x, 'sampleName')]
        if self.debug:
            print "DEBUG:", len(configSamples), " samples found."
            
        for sample in configSamples:
            sampleName = config.get(sample, 'sampleName')
            sampleType = config.get(sample,'sampleType')
            cut = config.get(sample, 'cut') if config.has_option(sample, 'cut') else '1'

            specialweight = config.get(sample, 'specialweight') if config.has_option(sample, 'specialweight') else "1"
            fullname = config.get(sample, 'sampleName')
            mergeCachingSize = int(config.get(sample, 'mergeCachingSize')) if config.has_option(sample, 'mergeCachingSize') else -1

            #fill the sample
            newsample = Sample(sampleName, sampleType)
            newsample.addtreecut = cut
            newsample.identifier = sample
            newsample.weightexpression = '' 
            newsample.specialweight = specialweight
            newsample.lumi = lumi
            newsample.prefix = '' 
            newsample.FullName = sampleName 

            if mergeCachingSize > 0:
                newsample.mergeCachingSize = mergeCachingSize
            if config.has_option(sample, 'skipParts'):
                newsample.skipParts = eval(config.get(sample, 'skipParts'))

            newsample.index = eval(config.get(sample, 'sampleIndex')) if config.has_option(sample, 'sampleIndex') else -999999

            # add and fill all the subsamples
            if config.has_option(sample,'subsamples') and eval(config.get(sample,'subsamples')):
                subgroups = eval((config.get(sample,'sampleGroup')))
                try:
                    subnames = eval((config.get(sample, 'subnames')))
                except:
                    # create subnames automatically based on subgroup name to avoid duplication
                    try:
                        shortname = config.get(sample, 'shortname').strip()
                    except:
                        # use full name if no short name given
                        shortname = sampleName
                    # if suffixes explicitly given use them
                    if config.has_option(sample, 'suffixes'):
                        suffixes = eval(config.get(sample, 'suffixes'))
                        subnames = [shortname + '_' + x for x in suffixes]
                    # otherwise use sample group as suffix
                    else:
                        subnames = [shortname + '_' + x for x in subgroups]
                subcuts = eval((config.get(sample, 'subcuts')))

                if sampleType != 'DATA':
                    subxsecs = eval((config.get(sample, 'xSec')))
                    if type(subxsecs) != list:
                        subxsecs = [subxsecs] * len(subcuts)

                    if len(list(set(subxsecs))) == 1:
                        newsample.xsec = [subxsecs[0]]
                    else:
                        print "\x1b[31mWARNING: different cross sections for the sub-samples of", sampleName, " are you sure you want to do this?\x1b[0m"
                    subsfs = eval((config.get(sample, 'SF'))) if config.has_option(sample, 'SF') else [1.0]*len(subxsecs)

                try:
                    subspecialweights = eval((config.get(sample, 'specialweight')))
                    if len(subspecialweights) < 2:
                        subspecialweights = []
                        print "\x1b[31mWARNING: specialweight not defined for subsamples but for full sample only!\x1b[0m"
                except:
                    subspecialweights = []

                subindices = None
                newsamples = []
                for i,cut in enumerate(subcuts):
                    newsubsample = copy(newsample)
                    newsubsample.subsample = True
                    newsubsample.name = subnames[i]
                    newsubsample.subcut = subcuts[i]
                    newsubsample.group = subgroups[i]
                    if sampleType != 'DATA':
                        newsubsample.sf = float(subsfs[i])
                        newsubsample.xsec = float(subxsecs[i])
                    if len(subspecialweights) == len(subcuts):
                        newsubsample.specialweight = subspecialweights[i] 

                    if type(newsample.index) == list:
                        newsubsample.index = newsample.index[i]
                    if config.has_option(sample, 'offsets'):
                        sampleIndexOffsets = eval(config.get(sample, 'offsets'))
                        if type(sampleIndexOffsets) == list and len(sampleIndexOffsets) == len(subcuts):
                            newsubsample.index += sampleIndexOffsets[i]
                        else:
                            print "sample:", newsample.name
                            print "=>", type(sampleIndexOffsets), len(sampleIndexOffsets), len(subcuts)
                            print "\x1b[31mERROR: sampleIndex offset does not match subcuts\x1b[0m"
                            raise Exception("ConfigError")

                    newsamples.append(newsubsample)

                self._samplelist.extend(newsamples)
                self._samplelist.append(newsample)
            else:
                if sampleType != 'DATA':
                    newsample.xsec = eval((config.get(sample,'xSec')))    
                    newsample.sf = eval((config.get(sample, 'SF'))) if config.has_option(sample, 'SF') else 1.0
                newsample.group = config.get(sample,'sampleGroup')
                self._samplelist.append(newsample)

    def __iter__(self):
        for sample in self._samplelist:
            if sample.active:
                yield sample

    def get_sample(self, samplename):
        '''return the sample whose name matches the sample.name'''
        for sample in self._samplelist:
            if sample.name == samplename:
                return sample
        return None
    
    def get_samples(self, samplenames=''):
        '''Samplenames is list of the samples names. Returns a list of samples corresponding to the names'''
        samples = []
        thenames = []
        #for splitted samples use the identifier. There is always only one. if list, they are all true
        if (len(samplenames)>0 and self.checkSplittedSampleName(samplenames[0])):
          print "The samples is splitted"
          for sample in self._samplelist:
                  if (sample.subsample): continue #avoid multiple submissions from subsamples
                  print '@DEBUG: samplenames ' + samplenames[0]
                  print '@DEBUG: sample identifier ' + sample.identifier
                  if sample.identifier == samplenames[0]:
                          samples.append(sample)
                          thenames.append(sample.name)
        #else check the name
        else:
            if not isinstance(samplenames, (list, tuple)): 
                samplenames = [x.strip() for x in samplenames.split(',')]

            # filter empty sample names
            samplenames = [x for x in samplenames if len(x.strip()) > 0]

            for samplename in samplenames:
                found = False
                for sample in self._samplelist:
                    if sample.name == samplename:
                        found = True
                        samples.append(sample)
                        thenames.append(sample.name)
                if not found:
                    print "\x1b[31mERROR: sample not found:", samplename, "\x1b[0m"
                    print "requested:", samplenames
                    print "existing:", [x.name for x in self._samplelist]
                    raise Exception("SampleMissing")
                
        return samples

    # return list of all identifiers
    def getSampleIdentifiers(self):
        return list(set([x.identifier for x in self]))

    #
    def getFullSample(self, sampleIdentifier):
        fullSamples = [x for x in self if x.identifier == sampleIdentifier and x.subsample == False]
        if len(fullSamples) != 1:
            print("ERROR: sample not found/ not unique! :", fullSamples)
            raise Exception("SampleMissing")
        return fullSamples[0]
    
    def getSubsamples(self, sampleIdentifier):
        return [x for x in self if x.identifier == sampleIdentifier and x.subsample == True]

    # return list of all samples with matching identifier
    def find(self, identifier, subsamples=False):
        return [x for x in self if (
                    ('*' in identifier and fnmatch.fnmatch(x.identifier, identifier)) or 
                    ('*' not in identifier and x.identifier==identifier)
                ) and (not x.subsample or subsamples)
            ]

    @staticmethod
    def filterIdentifiers(sampleIdentifiers, samples):
        usedIdentifiers = [x.identifier for x in samples]
        return [x for x in sampleIdentifiers if x in usedIdentifiers]

    # DEPRECATED
    #it checks whether filename is a splitted sample or is a pure samples and returns the file name without the _#
    def checkSplittedSample(self, filename):
            try:
                    isinstance( eval(filename[filename.rfind('_')+1:] ) , int )
                    print '@DEBUG: fileName in CHECKSPLITTEDSAMPLE : ' + filename
                    print '@DEBUG: return in CHECKSPLITTEDSAMPLE : ' + filename[:filename.rfind('_')]
                    return filename[:filename.rfind('_')]
            except:
                    return filename
    # DEPRECATED
    #bool
    def checkSplittedSampleName(self,filename):
            #print '### CHECKSPLITTEDSAMPLENAME ###',filename
            # if there is an underscore in the filename
            if ( filename.rfind('_') > 0. ) :
                    try:
                            return isinstance( eval(filename[filename.rfind('_')+1:] ) , int )
                    except:
                            return False
            else:
                    return False

    
if __name__ == '__main__':

    # read config

    #configSamples = ['WminusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8', 'WplusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8', 'ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8', 'ggZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8', 'ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8', 'ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8', 'ZZ_TuneCP5_13TeV-pythia8', 'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8', 'WZ_TuneCP5_13TeV-pythia8', 'WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8', 'WW_TuneCP5_13TeV-pythia8', 'WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8',.....]
    ## this is what the XbbConfigReader module is doing:
    pathconfig = BetterConfigParser()
    pathconfig.read('/mnt/t3nfs01/data01/shome/krgedia/CMSSW_10_1_0/src/Xbb/python/Wlv2018config/paths.ini') #parent class 'ConfigParser' method
    configFiles = pathconfig.get('Configuration', 'List').split(' ') 
    config = BetterConfigParser()
    print ('configFiles parsed', configFiles)
    for configFile in configFiles:
        config.read('Wlv2017config/' + configFile)

    #print(config.get('Weights','weightF'))
    #config = XbbConfigReader.read('Zvv2017')

    #inputFile = 'root://t3dcachedb03.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/berger_p2/VHbb/VHbbPostNano2017/V5/Zvv/rerun/v4j/eval/ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8/tree_aa5e971734ef4e885512748d534e6937ff03dc61feed21b6772ba943_000000_000000_0000_9_a6c5a52b56e5e0c7ad5aec31429c8926bf32cf39adbe087f05cfb323.root'
    path = 'root://t3dcachedb03.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/krgedia/VHbb/Wlv/VHbbPostNano2018/mva/17jan20v4' 
    #samplefiles = '../samples/VHbbPostNano2017_V5/merged_Zvv2017/' 
    #samplesinfo = 'Zvv2017config/samples_nosplit.ini' 
    info = ParseInfo(samples_path=path, config=config)
    
    #print('_samplelist ',info._samplelist)
    
    #sample = [x for x in info if x.identifier == 'ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8'][0]
    #print('sample ',sample)
    
    # read sample
    #sampleTree = SampleTree([inputFile], config=config)
    #print 'sampleTree', sampleTree

