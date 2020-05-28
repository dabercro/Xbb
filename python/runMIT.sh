#! /bin/bash

infile=$1
config=$2

export SCRAM_ARCH=slc6_amd64_gcc630
export CMSSW_VER=CMSSW_10_1_0

source /cvmfs/cms.cern.ch/cmsset_default.sh

scramv1 project CMSSW $CMSSW_VER

cd $CMSSW_VER

tar -xf ../xbb_condor.tgz

cd src
eval `scram runtime -sh`

cp ../../*.ini Xbb/python/.

cd ../..

cd $CMSSW_VER/src/Xbb/python

filelist=$(python -c "from myutils.FileList import FileList; print(FileList.compress('$infile'))")
outfile=$(python -c "from myutils.FileLocator import FileLocator; from myutils.BetterConfigParser import BetterConfigParser; config = BetterConfigParser(); config.read('$config'); print(FileLocator(config).getFilenameAfterPrep('$infile'))")

./runAll.sh run_ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8_part0 Zvv2018 run 1 noid --inputDir=PREPin --sampleIdentifier=ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8 --addCollections=Prep.VHbb --fileList=$filelist --outputDir=PREPout --noretry

./runAll.sh sysnew_ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8_part0 Zvv2018 sysnew 1 noid  --inputDir=SYSin  --sampleIdentifier=ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8  --addCollections=Sys.all --fileList=$filelist --outputDir=SYSout --noretry
