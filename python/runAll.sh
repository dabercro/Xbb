#!/bin/bash

#====================================================================
#
#        FILE: runAll.sh
#
#       USAGE: runAll.sh sample energy task
#
# DESCRIPTION: Script to be launched in the batch system.
#              Can also be used, with some care, to run locally.
#
#      AUTHOR: VHbb team
#              ETH Zurich
#
#====================================================================

hostname

export SCRAM_ARCH=slc6_amd64_gcc630
export CMSSW_VER=CMSSW_10_1_0

source /cvmfs/cms.cern.ch/cmsset_default.sh

if [ -f xbb_condor.tgz ]
then

    scramv1 project CMSSW $CMSSW_VER

    cd $CMSSW_VER

    tar -xf ../xbb_condor.tgz

    cd src
    eval `scram runtime -sh`

fi

cd ${CMSSW_BASE}/src/Xbb/python/
echo "cd ${CMSSW_BASE}/src/Xbb/python/"
eval `scram runtime -sh`
STARTTIME=$(date +%s.%N)
ulimit -c 0
ulimit -S -c 0

# Fix Python escape sequence bug.
export TERM=""

#-------------------------------------------------
# Parse Input Arguments

# kept for backwards compatibility, please use named arguments for new code!!!

sample=$1     # The sample to run on. It must match a sampleName in samples_nosplit.ini.
tag=$2        # The analysis configuration tag, e.g. 13TeV.
task=$3       # The task to perform.
nprocesses=$4 # Dummy variable used to shift the other parameters by +1. FIXME: Remove this argument?
job_id=$5     # Needed for split step and train optimisation. FIXME: It does not have a unique meaning.
bdt_params=$6 # The set of hyperparameters for BDT optimisation.
filelist=$7   # Needed to run the prep and sys step with a limited number of files per job.
region=$8     # plotting region

# Verify the number of input arguments.
if [ $# -lt 3 ]; then
    echo "RuntimeError: At least 3 arguments are required, e.g."
    echo "./runAll.sh sample tag task"
    exit
elif [ $task = "mva_opt" -a $# -lt 5 ]; then
    echo "RuntimeError: At least 5 arguments are required for BDT hyperparameter optimisation, e.g."
    echo "./runAll.sh sample tag mva_opt job_id bdt_params"
    exit
fi

echo "Task: $task"
echo "Pwd: $PWD"
echo "Host: $HOSTNAME"
echo "User: $USER"
echo "StartTime: "`date +"%Y-%m-%d %T"`
echo


#-------------------------------------------------
# parse named input arguments
# todo: pass everything as named argument
force="0"
unblind="0"
friend="0"
join="0"
noretry="0"
while [ $# -gt 0 ]; do
  case "$1" in
    --inputDir=*)
      inputDir="${1#*=}"
      ;;
    --outputDir=*)
      outputDir="${1#*=}"
      ;;
    --trainingRegions=*)
      trainingRegions="${1#*=}"
      ;;
    --regions=*)
      regions="${1#*=}"
      ;;
    --process=*)
      process="${1#*=}"
      ;;
    --force)
      force="1"
      ;;
    --unblind)
      unblind="1"
      ;;
    --friend)
      friend="1"
      ;;
    --join)
      join="1"
      ;;
    --noretry)
      noretry="1"
      ;;
    --expectedSignificance)
      expectedSignificance="1"
      ;;
    --verbose)
      verbose="1"
      ;;
    --vars=*)
      vars="${1#*=}"
      ;;
    --sampleIdentifier=*)
      sampleIdentifier="${1#*=}"
      ;;
    --sampleName=*)
      sampleName="${1#*=}"
      ;;
    --splitFilesChunkSize=*)
      splitFilesChunkSize="${1#*=}"
      ;;
    --chunkNumber=*)
      chunkNumber="${1#*=}"
      ;;
    --splitFilesChunks=*)
      splitFilesChunks="${1#*=}"
      ;;
    --fileList=*)
      fileList="${1#*=}"
      ;;
    --limit=*)
      limit="${1#*=}"
      ;;
    --addCollections=*)
      addCollections="${1#*=}"
      ;;
    --configFile=*)
      configFile="${1#*=}"
      ;;
    *)
      ;;
  esac
  shift
done


#-------------------------------------------------
# Setup Environment

# Change to the job submission directory when using lxbatch.
if [ -n "${LS_SUBCWD-}" ]; then
    cd $LS_SUBCWD
fi

echo "Parsing files in ${tag}config..."
echo


#-------------------------------------------------
# Run Task

if [ $task = "prep" ]; then
    # deprecated, use 'run' for this
    runCommand="python ./prepare_environment_with_config.py"
    if [ "$limit" ]; then runCommand="${runCommand} --limit ${limit}"; fi
    
elif [ $task = "run" ] || [ $task = "sysnew" ]; then
    runCommand="python ./sys_new.py"
    if [ "$limit" ]; then runCommand="${runCommand} --limit ${limit}"; fi
    if [ "$addCollections" ]; then runCommand="${runCommand} --addCollections ${addCollections}"; fi

elif [ $task = "cachetraining" ]; then
    runCommand="python ./cache_training.py --trainingRegions ${trainingRegions} --splitFilesChunkSize ${splitFilesChunkSize} --splitFilesChunks ${splitFilesChunks} --chunkNumber ${chunkNumber}"

elif [ $task = "runtraining" ]; then
    runCommand="python ./run_training.py --trainingRegions ${trainingRegions}"
    if [ "$expectedSignificance" = "1" ]; then runCommand="${runCommand} --expectedSignificance"; fi

elif [ $task = "dnn" ]; then
    runCommand="python tfZllDNN/train.py -c tfZllDNN/config.cfg -i ${trainingRegions} -l"
    config_filenames=()
    unset configFile

elif [ $task = "hadd" ]; then
    runCommand="python ./hadd.py --chunkNumber ${chunkNumber}"

elif [ $task = "cacheplot" ]; then
    runCommand="python ./cache_plot.py --regions ${regions} --splitFilesChunkSize ${splitFilesChunkSize} --splitFilesChunks ${splitFilesChunks} --chunkNumber ${chunkNumber}"

elif [ $task = "runplot" ]; then
    if [ -z "$vars" ]; then 
        runCommand="python ./run_plot.py --regions ${regions}";
    else
        runCommand="python ./run_plot.py --regions ${regions} --vars ${vars}";
    fi
elif [ $task = "postfitplot" ]; then
    if [ -z "$regions" ]; then
        runCommand="python ./postfit_plot.py"
    else
        runCommand="python ./postfit_plot.py --regions ${regions}"
    fi

elif [ $task = "cachedc" ]; then
    runCommand="python ./cache_dc.py --regions ${regions} --splitFilesChunkSize ${splitFilesChunkSize} --splitFilesChunks ${splitFilesChunks} --chunkNumber ${chunkNumber}"

elif [ $task = "rundc" ]; then
    if [ -z "$chunkNumber" ]; then
        runCommand="python ./run_dc.py --regions ${regions}"
    else
        runCommand="python ./run_dc.py --regions ${regions}  --chunkNumber ${chunkNumber}";
    fi

elif [ $task = "mergedc" ]; then
    runCommand="python ./merge_dc.py --regions ${regions}";

elif [ $task = "export_h5" ] || [ $task = "export_hdf5" ]; then
    runCommand="python ./write_numpy_array_for_training.py -t ${trainingRegions}"

elif [ $task = "make_skims" ]; then
    runCommand="python ./make_skims.py --regions ${regions}";

fi

# add standard arguments, print command and run 
if [ "$runCommand" ]; then
    if [ "$fileList" ]; then
        runCommand="${runCommand} --fileList ${fileList}"
    fi
    if [ "$force" = "1" ]; then
        runCommand="${runCommand} --force"
    fi
    if [ "$friend" = "1" ]; then
        runCommand="${runCommand} --friend"
    fi
    if [ "$join" = "1" ]; then
        runCommand="${runCommand} --join"
    fi
    if [ "$sampleIdentifier" ]; then
        runCommand="${runCommand} --sampleIdentifier ${sampleIdentifier}"
    fi
    if [ "$verbose" = "1" ]; then
        runCommand="${runCommand} --verbose"
    fi
    
    if [ "$configFile" ]; then 
        runCommand="${runCommand} --config ${configFile}"
    else
        runCommand="${runCommand} ${config_filenames[@]}"
    fi
    if [ "$inputDir" ]; then
        runCommand="${runCommand} --inputDir ${inputDir}"
    fi
    if [ "$outputDir" ]; then
        runCommand="${runCommand} --outputDir ${outputDir}"
    fi
    if [ "$unblind" = "1" ]; then
        runCommand="${runCommand} --unblind"
    fi

    echo "$runCommand"
    eval "$runCommand"
fi


EXITCODE=$?
echo "--------------------------------------------------------------------------------"
echo "exit code: $EXITCODE"
echo "EndTime: "`date +"%Y-%m-%d %T"`
ENDTIME=$(date +%s.%N)
DIFFTIME=$(echo "($ENDTIME - $STARTTIME)/60" | bc)
echo "duration (real time): $DIFFTIME minutes"
if [ "$EXITCODE" -ne "0" ]; then
    if [ "$noretry" = "1" ]; then
        echo "--- STOP ---"
    else
        echo "--- RETRY ---"
        if [ "$runCommand" ]; then
            echo "$runCommand"
            eval "$runCommand"
            EXITCODE=$?
            echo "exit code: $EXITCODE"
            ENDTIME=$(date +%s.%N)
            DIFFTIME=$(echo "($ENDTIME - $STARTTIME)/60" | bc)
            echo "duration (real time): $DIFFTIME minutes, including retry"
        fi
    fi
else
    echo "--- OK ---"
fi
echo
echo "Exiting runAll.sh"
echo
exit $EXITCODE
