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

if [ "$(expr substr $tag 1 3)" = "Zll" ]
then

    ./runAll.sh run_${sample}_part$part $tag run 1 noid --inputDir=KINFITin --sampleIdentifier=$sample  --addCollections=KinematicFit.fitter  --fileList=$filelist  --outputDir=KINFITout --noretrxy

fi

if [ "$(./submitMIT.py -T $tag -c 'Directories:CONDORin')" = "$(./submitMIT.py -T $tag -c 'Directories:MVAout')" ]
then

    ./runAll.sh run_${sample}_part$part $tag run 1 noid --inputDir=MVAin --sampleIdentifier=$sample --addCollections=Eval.VH --fileList=$filelist --outputDir=MVAout --noretry

fi

pwd
find $(./submitMIT.py -T $tag -c 'Directories:CONDORin') -type f
