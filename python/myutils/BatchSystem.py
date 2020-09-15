from __future__ import print_function
import xml.etree.ElementTree
import subprocess
import fnmatch
import hashlib
import json
import importlib

class BatchJob(object):
    def __init__(self, jobName='', submitCommand='', log=''):
        self.data = {
                'jobName': jobName,
                'submitCommand': submitCommand,
                'log': log,
                'err': '',
                'outputFiles': []
                }

    def setProperty(self, prop, value):
        self.data[prop] = value

    def toDict(self):
        return self.data

# batch system base class / factory
class BatchSystem(object):
    # if command returns 0, batch system is installed
    autoDetectionCommnds = {
                'SLURM': 'command -v squeue', 
                'SGE': 'command -v qstat',
                'HTCondor': 'command -v condor_q',
                }

    def __init__(self, interactive=False, local=False, configFile=None):
        self.name = 'undefined'
        self.nJobsProcessed = 0
        self.nJobsSubmitted = 0
        self.nJobsSkipped = 0
        self.interactive = interactive
        self.submittedJobs = []
        self.runLocally = local
        #self.configFile = None if interactive else configFile
        self.configFile = configFile

    # factory
    @staticmethod
    def create(config, interactive=False, local=False, configFile=None):
        whereToLaunch = config.get('Configuration', 'whereToLaunch').strip()

        # automatic detection
        if whereToLaunch.lower() == 'auto':
            batchSystemsFound = []
            for batchSystem, detectionCommand in BatchSystem.autoDetectionCommnds.items():
                try:
                    if subprocess.call(detectionCommand, shell=True) == 0:
                        batchSystemsFound.append(batchSystem)
                except:
                    pass
            if len(batchSystemsFound) < 1:
                print("\x1b[31mERROR: no batch system could be detected! Configuration.whereToLaunch = auto will not work.\x1b[0m")
                raise Exception("NoBatchSystemFound")
            elif len(batchSystemsFound) > 1:
                print("\x1b[31mERROR: multiple batch systems found:", ",".join(batchSystemsFound), " use Configuration.whereToLaunch to specify which one to use, not auto.\x1b[0m")
                raise Exception("MultipleBatchSystemsFound")
            else:
                whereToLaunch = batchSystemsFound[0]
                print("INFO: auto-detected batch system:", whereToLaunch)

        # for backward compatibility
        if whereToLaunch == 'PSI':
            moduleName = 'BatchSystemSGE'
        elif whereToLaunch == 'lxplus-condor':
            moduleName = 'BatchSystemHTCondor'
        else:
            moduleName = 'BatchSystem' + whereToLaunch

        print("INFO: using batch system \x1b[34m", moduleName, "\x1b[0m")
        batchSystemClass = getattr(importlib.import_module(".{module}".format(module=moduleName), package="myutils"), moduleName)
        return batchSystemClass(config, interactive=interactive, local=local, configFile=configFile)

    # helper functions
    def getName(self):
        return self.name
    
    def getJobs(self):
        raise Exception("not implemented")

    def get_single_job_name(self, jobPrefix, sampleIdentifier, fileNumber, jobSuffix):
        return jobPrefix + '_' + sampleIdentifier + '_part%d'%(fileNumber) + '_' + jobSuffix

    def job_for_file_exists(self, jobList, jobPrefix, sampleIdentifier, fileNumber, jobSuffix):
        # 1 file per job
        if self.get_single_job_name(jobPrefix, sampleIdentifier, fileNumber, jobSuffix) in jobList:
            return True

        # multiple files per job
        jobName = jobPrefix + '_' + sampleIdentifier + '_part*_files*'
        for runningJob in jobList:
            if fnmatch.fnmatch(runningJob, jobName):
                filesSpec = runningJob.split('_files')[1].split('_')[0]
                filesFrom = int(filesSpec.split('to')[0])
                filesTo = int(filesSpec.split('to')[1])
                if filesFrom <= fileNumber and fileNumber <= filesTo:
                    return True
        return False

    def getNJobsSubmitted(self):
        return self.nJobsSubmitted

    # -----------------------------------------------------------------------------
    # prepare RUN SCRIPT
    # (this is independent of batch system or local)
    # -----------------------------------------------------------------------------
    def getRunScriptCommand(self, repDict):
        runScript = 'runAll.sh %(job)s %(config)s '%(repDict)
        runScript += repDict['task'] + ' ' + repDict['nprocesses'] + ' ' + repDict['job_id'] + ' ' + repDict['additional']

        # add named arguments to run script
        if 'arguments' in repDict:
            for argument, value in repDict['arguments'].iteritems():
                runScript += (' --{argument}={value} '.format(argument=argument, value=value)) if len('{value}'.format(value=value)) > 0 else ' --{argument} '.format(argument=argument)

        # add path to config file (Can be dumped combined config instead of single files!)
        if self.configFile:
            runScript += " --configFile=" + self.configFile
        
        if self.interactive or self.runLocally:
            runScript += " --noretry"

        return runScript

    def getLogPaths(self, repDict):
        # log=batch system log, out=stdout of process
        logPaths = {
                'log': '%(logpath)s/%(task)s_%(timestamp)s_%(job)s_%(en)s_%(additional)s.log' %(repDict),
                'error': '%(logpath)s/%(task)s_%(timestamp)s_%(job)s_%(en)s_%(additional)s.err' %(repDict),
                'out': '%(logpath)s/%(task)s_%(timestamp)s_%(job)s_%(en)s_%(additional)s.out' %(repDict),
                'config': "%(logpath)s/%(timestamp)s_%(task)s.config" %(repDict),
            }
        return logPaths

    def submitPreprocess(self, job, repDict):
        repDict['job'] = job
        repDict['name'] = '%(job)s_%(en)s%(task)s' %repDict

    def runShell(self, commands):
        try:
            return subprocess.check_output(commands, shell=True)
        except:
            return None

    def run(self, command, runScript='', repDict={}, getJobIdFn=None):
        # -----------------------------------------------------------------------------
        # RUN command
        # -----------------------------------------------------------------------------
        if command:
            batchJob = BatchJob(repDict['name'], command, self.getLogPaths(repDict)['out'])
            batchJob.setProperty('repDict', repDict)

            if self.interactive or self.runLocally:
                if self.runLocally:
                    answer = 'l'
                else:
                    print("SUBMIT:\x1b[34m", command, "\x1b[0m\n(press ENTER to run it and continue, \x1b[34ml\x1b[0m to run it locally, \x1b[34md\x1b[0m for debug mode, \x1b[34ma\x1b[0m to run all jobs locally and \x1b[34ms\x1b[0m to submit the remaining jobs)")
                    answer = raw_input().strip()
                if answer.lower() in ['no', 'n']:
                    self.nJobsSkipped += 1
                    return
                elif answer.lower() == 's':
                    self.interactive = False
                    self.nJobsSubmitted += 1
                elif answer.lower() in ['l', 'local', 'a','d']:
                    if answer.lower() == 'a':
                        self.runLocally = True
                    print("run locally")
                    command = 'sh {runscript}'.format(runscript=runScript)
                    if answer.lower() == 'd':
                        command = "XBBDEBUG=1 " + command
                else:
                    self.nJobsSubmitted += 1

            else:
                print("the command is ", command)
                self.nJobsSubmitted += 1

            if getJobIdFn and not (self.interactive or self.runLocally):
                try:
                    stdOutput = subprocess.check_output([command], shell=True)
                except Exception as e:
                    print("\x1b[31mERROR: SUBMISSION FAILED!", e,"\x1b[0m")
                    stdOutput = ""

                try:
                    jobId = getJobIdFn(stdOutput)
                    batchJob.setProperty('id', jobId)
                except Exception as e:
                    print(e)
            else:
                subprocess.call([command], shell=True)

            self.submittedJobs.append(batchJob)
            return batchJob
    
    def submitQueue(self):
        pass

    def dumpSubmittedJobs(self, fileName):
        with open(fileName, 'w') as outfile:
            json.dump([x.toDict() for x in self.submittedJobs], outfile)

    def cancelJob(self, job):
        print("WARNING: cancel job not implemented for this batch system!")


