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

# If on Condor, run this setup

if [ -f xbb_condor.tgz ]
then

    export SCRAM_ARCH=slc6_amd64_gcc630
    export CMSSW_VER=CMSSW_10_1_0

    source /cvmfs/cms.cern.ch/cmsset_default.sh

    scramv1 project CMSSW $CMSSW_VER

    cd $CMSSW_VER

    tar -xf ../xbb_condor.tgz

    cd src
    eval `scram runtime -sh`

fi


cd ${CMSSW_BASE}/src/Xbb/python/
echo "cd ${CMSSW_BASE}/src/Xbb/python/"
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

# The job submission environment.
whereToLaunch=$(
python - << END
import myutils
parser = myutils.BetterConfigParser()
parser.read('${tag}config/paths.ini')
print parser.get('Configuration', 'whereToLaunch')
END
)

# The configuration file names, formatted as arguments to the task scripts.
config_filenames=( $(
python - << END
import myutils
parser = myutils.BetterConfigParser()
parser.read('${tag}config/paths.ini')
print parser.get('Configuration', 'List')
END
) )

echo "Configuration Files: ${config_filenames[@]}"
echo

for (( i=0; i<${#config_filenames[@]}; i++ )); do
    config_filenames[$i]="--config ${tag}config/${config_filenames[$i]}"
done

# The log file path.
logpath=$(
python - << END
import myutils
parser = myutils.BetterConfigParser()
parser.read('${tag}config/paths.ini')
try:
  parser.read('${tag}config/volatile.ini')
except:
  pass
print parser.get('Directories', 'logpath')
END
)

if [ ! -d $logpath ]; then
    mkdir -p $logpath
fi

# The MVA list, only compute for the eval steps!
if [[ $task = *"eval"* ]]; then
MVAList=$(
python << END
import myutils
parser = myutils.BetterConfigParser()
parser.read('${tag}config/paths.ini')
configList = ['${tag}config/' + x for x in parser.get('Configuration', 'List').strip().split(' ')]
config= myutils.BetterConfigParser()
config.read(configList)
print config.get('MVALists', 'List_for_submitscript')
END
)
fi

#-------------------------------------------------
# Run Task

if [ $task = "prep" ]; then
    runCommand="python ./prepare_environment_with_config.py"
    if [ "$limit" ]; then runCommand="${runCommand} --limit ${limit}"; fi
    
elif [ $task = "sysnew" ] || [ $task = "run" ]; then
    runCommand="python ./sys_new.py"
    if [ "$limit" ]; then runCommand="${runCommand} --limit ${limit}"; fi
    if [ "$addCollections" ]; then runCommand="${runCommand} --addCollections ${addCollections}"; fi

elif [ $task = "singleprep" ]; then
    echo "python ./prepare_environment_with_config.py --samples $sample ${config_filenames[@]} --filelist ...(${#filelist} char)"
    python ./prepare_environment_with_config.py --samples $sample ${config_filenames[@]} --filelist $filelist

elif [ $task = "mergesingleprep" ]; then
    echo "python ./myutils/mergetreePSI.py --samples $sample ${config_filenames[@]}"
    python ./myutils/mergetreePSI.py --samples $sample ${config_filenames[@]}

elif [ $task = "trainReg" ]; then
    echo "python ./trainRegression.py --config ${tag}config/regression.ini ${config_filenames[@]}"
    python ./trainRegression.py --config ${tag}config/regression.ini ${config_filenames[@]}

elif [ $task = "sys" ]; then
    runCommand="python ./write_regression_systematics.py --samples ${sampleIdentifier}"
    if [ "$limit" ]; then runCommand="${runCommand} --limit ${limit}"; fi

elif [ $task = "vars" ]; then
    echo "python ./write_newVariables.py --samples $sample ${config_filenames[@]}"
    python ./write_newVariables.py --samples $sample ${config_filenames[@]}

elif [ $task = "singlesys" ]; then
    echo "python ./write_regression_systematics.py --samples $sample ${config_filenames[@]} --filelist ...(${#filelist} char)"
    python ./write_regression_systematics.py --samples $sample ${config_filenames[@]} --filelist $filelist

elif [ $task = "mergesinglesys" ]; then
    echo "python ./myutils/mergetreePSI.py --samples $sample ${config_filenames[@]}  --mergesys True"
    python ./myutils/mergetreePSI.py --samples $sample ${config_filenames[@]} --mergesys True

elif [ $task = "reg" ]; then
    echo "python ./only_regression.py --samples $sample ${config_filenames[@]}"
    python ./only_regression.py --samples $sample ${config_filenames[@]}

elif [ $task = "eval" ]; then
    runCommand="python  ./evaluateMVA.py --discr $MVAList"

elif [ $task = "eval_scikit" ]; then
    runCommand="python  ./evaluateMVA_scikit.py --discr $MVAList"

elif [ $task = "singleeval" ]; then
    echo "python ./evaluateMVA.py --discr $MVAList --samples $sample ${config_filenames[@]} --filelist ...(${#filelist} char)"
    python ./evaluateMVA.py --discr $MVAList --samples $sample ${config_filenames[@]} --filelist $filelist

elif [ $task = "mergesingleeval" ]; then
    echo "python ./myutils/mergetreePSI.py --samples $sample ${config_filenames[@]}  --mergeeval True"
    python ./myutils/mergetreePSI.py --samples $sample ${config_filenames[@]} --mergeeval True

elif [ $task = "syseval" ]; then
    echo "python ./write_regression_systematics.py --samples $sample ${config_filenames[@]}"
    python ./write_regression_systematics.py --samples $sample ${config_filenames[@]}
    echo "python ./evaluateMVA.py --discr $MVAList --samples $sample ${config_filenames[@]}"
    python ./evaluateMVA.py --discr $MVAList --samples $sample ${config_filenames[@]}

elif [ $task = "train" ] || [ $task = "splitsubcaching" ]; then
    echo "python ./train.py --training $sample ${config_filenames[@]} --setting $bdt_params --local True"
    python ./train.py --training $sample ${config_filenames[@]} --setting $bdt_params --local True

elif [ $task = "cachetraining" ]; then
    runCommand="python ./cache_training.py --trainingRegions ${trainingRegions} --splitFilesChunkSize ${splitFilesChunkSize} --splitFilesChunks ${splitFilesChunks} --chunkNumber ${chunkNumber}"

elif [ $task = "runtraining" ]; then
    runCommand="python ./run_training.py --trainingRegions ${trainingRegions}"
    if [ "$expectedSignificance" = "1" ]; then runCommand="${runCommand} --expectedSignificance"; fi

elif [ $task = "runtraining_scikit" ]; then
    runCommand="python ./run_training_scikit.py --trainingRegions ${trainingRegions}"

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

elif [ $task = "mva_opt" ] || [ $task = "splitsubcaching" ]; then
    echo "python ./train.py --training $sample ${config_filenames[@]} --setting $bdt_params --local True"
    python ./train.py --training $sample ${config_filenames[@]} --setting $bdt_params --local True

elif [ $task = "singleplot" ]; then
    echo "./tree_stack.py --region $region ${config_filenames[@]} --filelist ...(${#filelist} char)"
    ./tree_stack.py --region $region ${config_filenames[@]} --filelist $filelist

elif [ $task = "mergesingleplot" ]; then
    echo "./tree_stack.py --region $region ${config_filenames[@]} --mergeplot True"
    ./tree_stack.py --region $region ${config_filenames[@]} --filelist $sample --mergeplot True

elif [ $task = "plot" ]; then
    echo "python ./tree_stack.py --region $sample ${config_filenames[@]} --settings $bdt_params"
    python ./tree_stack.py --region $sample ${config_filenames[@]} --settings $bdt_params

elif [ $task = "mergecachingplotvar" ]; then
    echo "python ./tree_stack.py --region $sample ${config_filenames[@]} --settings $bdt_params --mergecachingplot"
    python ./tree_stack.py --region $sample ${config_filenames[@]} --settings $bdt_params --mergecachingplot

elif [ $task = "splitvarplot" ]; then
    echo "python ./tree_stack.py --region $sample ${config_filenames[@]} --settings $bdt_params"
    python ./tree_stack.py --region $sample ${config_filenames[@]} --settings $bdt_params

elif [ $task = "splitcaching" ]; then
    echo "./tree_stack.py --region $sample ${config_filenames[@]} --settings $bdt_params --mergeplot False"
    ./tree_stack.py --region $sample ${config_filenames[@]} --settings $bdt_params --mergeplot False

elif [ $task = "mergecaching" ]; then
    echo "./tree_stack.py --region $sample ${config_filenames[@]} --settings $bdt_params --mergeplot False  --filelist $filelist"
    ./tree_stack.py --region $sample ${config_filenames[@]} --settings $bdt_params --mergeplot False --filelist $filelist

elif [ $task = "mergecachingplot" ]; then
    echo "python ./tree_stack.py --region $sample ${config_filenames[@]} --settings $bdt_params --mergecachingplot"
    python ./tree_stack.py --region $sample ${config_filenames[@]} --settings $bdt_params --mergecachingplot

elif [ $task = "mergesubcaching" ]; then
    echo "./train.py --training $sample ${config_filenames[@]} --settings $bdt_params --mergeplot False  --filelist $filelist"
    ./train.py --training $sample ${config_filenames[@]} --setting $bdt_params --mergeplot False  --filelist $filelist
elif [ $task = "mergesubcachingtrain" ]; then
    echo "./train.py --training $sample ${config_filenames[@]} --settings $bdt_params --mergeplot False"
    ./train.py --training $sample ${config_filenames[@]} --setting $bdt_params --mergeplot False --mergecachingplot
elif [ $task = "mergesyscaching" ]; then
    #echo "./tree_stack.py --region $sample ${config_filenames[@]} --settings $bdt_params --mergeplot False  --filelist $filelist"
    #./tree_stack.py --region $sample ${config_filenames[@]} --settings $bdt_params --mergeplot False --filelist $filelist
    echo "python ./workspace_datacard.py --variable $sample ${config_filenames[@]} --settings $bdt_params --mergeplot False  --filelist $filelist"
    python ./workspace_datacard.py --variable $sample ${config_filenames[@]} --settings $bdt_params --mergeplot False --filelist $filelist

elif [ $task = "mergesyscachingdc" ]; then
    echo "python ./workspace_datacard.py --variable $sample ${config_filenames[@]} --mergecachingplot True"
    python ./workspace_datacard.py --variable $sample ${config_filenames[@]} --mergecachingplot True

elif [ $task = "mergesyscachingdcsplit" ]; then
    echo "python ./workspace_datacard.py --variable $sample ${config_filenames[@]} --mergecachingplot True --settings $bdt_params"
    python ./workspace_datacard.py --variable $sample ${config_filenames[@]} --mergecachingplot True --settings $bdt_params
elif [ $task = "mergesyscachingdcmerge" ]; then
    echo "python ./workspace_datacard.py --variable $sample ${config_filenames[@]} --mergecachingplot True --settings $bdt_params"
    python ./workspace_datacard.py --variable $sample ${config_filenames[@]} --mergecachingplot True --settings $bdt_params

elif [ $task = "checksingleprep" ] || [ $task = "checksinglesys" ] || [ $task = "checksingleeval" ] || [ $task = "checksingleplot" ]; then
    if [[ $region ]]; then
        echo "./myutils/check_singlestep.py --region $region ${config_filenames[@]} --task $task  --sample $sample --filelist ...(${#filelist} char)"
        ./myutils/check_singlestep.py --region $region ${config_filenames[@]} --task $task --filelist $filelist --sample $sample
    else
        echo "./myutils/check_singlestep.py ${config_filenames[@]} --task $task  --sample $sample --filelist ...(${#filelist} char)"
        ./myutils/check_singlestep.py ${config_filenames[@]} --task $task --sample $sample --filelist $filelist
    fi
    exitcode=$?
    # echo $exitcode
    exit $exitcode
elif [ $task = "dc" ]; then
    echo "python ./workspace_datacard.py --variable $sample ${config_filenames[@]}"
    python ./workspace_datacard.py --variable $sample ${config_filenames[@]}
elif [ $task = "splitcachingdc" ]; then
    echo "python ./workspace_datacard.py --variable $sample ${config_filenames[@]} --settings $bdt_params"
    python ./workspace_datacard.py --variable $sample ${config_filenames[@]} --settings $bdt_params

elif [ $task = "split" ]; then
    echo "python ./split_tree.py --samples $sample ${config_filenames[@]} --max-events $job_id"
    python ./split_tree.py --samples $sample ${config_filenames[@]} --max-events $job_id

elif [ $task = "stack" ]; then
    echo "python ./manualStack.py --config ${config_filenames[@]}"
    python ./manualStack.py ${config_filenames[@]}

elif [ $task = "plot_sys" ]; then
    echo "python ./plot_systematics.py ${config_filenames[@]}"
    python ./plot_systematics.py ${config_filenames[@]}
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
