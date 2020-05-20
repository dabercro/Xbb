from __future__ import print_function
import xml.etree.ElementTree
import subprocess
import fnmatch
import hashlib
import json
from BatchSystem import BatchSystem

class BatchSystemHTCondor(BatchSystem):
    
    def __init__(self, config=None, interactive=False, local=False, configFile=None):
        super(BatchSystemHTCondor, self).__init__(interactive=interactive, local=local, configFile=configFile)
        self.name = 'HTCondor'
        self.config = config
        self.noBatch = False
        self.templateFileName = 'batch/condor/template.sub'
        self.template = None
        self.runLocally = local
        self.condorBatchGroups = {}

    def loadTemplate(self):
        with open(self.templateFileName, 'r') as templateFile:
            self.template = templateFile.read()

    def getJobNames(self):
        p = subprocess.Popen(["condor_q", "-nobatch"], stdout=subprocess.PIPE, shell=True)
        out, err = p.communicate()
        return out.split("\n")

    def submit(self, job, repDict):
        self.nJobsProcessed += 1
        self.submitPreprocess(job, repDict)

        runScript = self.getRunScriptCommand(repDict)
        logPaths = self.getLogPaths(repDict)

        if not self.template:
            self.loadTemplate()

        # -----------------------------------------------------------------------------
        # CONDOR
        # -----------------------------------------------------------------------------
        firstFileOfBatch = False
        isBatched = 'batch' in repDict and not self.noBatch
        if isBatched:
            if repDict['batch'] not in self.condorBatchGroups:
                # first file of batch -> make new submit file
                firstFileOfBatch = True
                self.condorBatchGroups[repDict['batch']] = '%(task)s_%(timestamp)s_%(batch)s'%(repDict)
            # use existing submit file and append
            dictHash = self.condorBatchGroups[repDict['batch']]
        else:
            # create a new submit file
            dictHash = '%(task)s_%(timestamp)s'%(repDict) + '_%x'%hash('%r'%repDict)

        condorDict = {
            'runscript': runScript.split(' ')[0],
            'arguments': ' '.join(runScript.split(' ')[1:]),
            'output': logPaths['out'],
            'log': logPaths['log'],
            'error': logPaths['error'],
            'queue': 'workday',
        }
        submitFileName = 'condor_{hash}.sub'.format(hash=dictHash)

        # append to existing bath
        if isBatched:
            with open(submitFileName, 'a') as submitFile:
                submitFile.write("\n")
                submitFile.write(self.template.format(**condorDict))
            command = None
        else:
            with open(submitFileName, 'w') as submitFile:
                submitFile.write(self.template.format(**condorDict))
            command = 'condor_submit {submitFileName}'.format(submitFileName=submitFileName)
        print("COMMAND:\x1b[35m", runScript, "\x1b[0m")
        return self.run(command, runScript, repDict)

    def submitQueue(self):
        for batchName, submitFileIdentifier in self.condorBatchGroups.iteritems():
            submitFileName = 'condor_{identifier}.sub'.format(identifier=submitFileIdentifier)
            command = 'condor_submit {submitFileName}  -batch-name {batchName}'.format(submitFileName=submitFileName, batchName=batchName)
            if self.interactive:
                print("SUBMIT:\x1b[34m", command, "\x1b[0m\n(press ENTER to run it and continue)")
                answer = raw_input().strip()
                if answer.lower() in ['no', 'n', 'skip']:
                    continue
            else:
                print("the command is ", command)
            subprocess.call([command], shell=True)

