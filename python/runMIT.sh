#! /bin/bash

tag=$1
sample=$2
filelist=$3
part=$4

export SCRAM_ARCH=slc6_amd64_gcc630
export CMSSW_VER=CMSSW_10_1_0

source /cvmfs/cms.cern.ch/cmsset_default.sh

scramv1 project CMSSW $CMSSW_VER

cd $CMSSW_VER

tar -xf ../xbb_condor.tgz

cd src
eval `scram runtime -sh`

cd Xbb/python

./runAll.sh run_${sample}_part$part $tag run 1 noid --inputDir=PREPin --sampleIdentifier=$sample --addCollections=Prep.VHbb --fileList=$filelist --outputDir=PREPout --noretry

./runAll.sh sysnew_${sample}_part$part $tag sysnew 1 noid  --inputDir=SYSin  --sampleIdentifier=$sample  --addCollections=Sys.all --fileList=$filelist --outputDir=SYSout --noretry

# Will need weights
# ./runAll.sh eval_${sample}_part$part $tag eval 1 noid  --inputDir=MVAin  --sampleIdentifier=$sample  --outputDir=MVAout  --fileList=$filelist  --noretry
