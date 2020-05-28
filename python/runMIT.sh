#! /bin/bash

filelist=

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

./runAll.sh run_ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8_part0 Zvv2018 run 1 noid --inputDir=PREPin --sampleIdentifier=ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8 --addCollections=Prep.VHbb --fileList=base64:eNoVzL0OgjAQAOAn0usdYnCEOMAgGEMYWC6QNC2J9pr+aHh7w/ZNH8QkQYMJkj14u0e2mzER7LqCS9m/dYSpXdenxNQvToAUVjAhgTFzy+0oTcPzKH3uMz+QSsZi1BN7+Vlt2O/JbksFr+y6rs4pfxxWx1MP9+/1dIjKAS8EpBQRMqEq8AZKKQUpaM10DiLpDxdMNyw= --outputDir=PREPout --noretry

./runAll.sh sysnew_ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8_part0 Zvv2018 sysnew 1 noid  --inputDir=SYSin  --sampleIdentifier=ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8  --addCollections='VHbbCommon.EWKweights,VHbbCommon.BTagSFDeepCSV_18,VHbbCommon.HiggsReco,VHbbCommon.isGGZH,VHbbCommon.vLeptons,Sys.newBranches,Sys.fill1,Sys.fill2,Sys.fill3,VHbbCommon.LOtoNLOweight,VHbbCommon.TTweights,VHbbCommon.DYspecialWeight,VHbbCommon.JetPUIDSF,Sys.metTriggerSF,Sys.PUweight,VHbbCommon.DoubleBtagSF,VHbbCommon.HeppyStyleGen,VHbbCommon.METcorr' --fileList=base64:eNoVzL0OgjAQAOAn0usdYnCEOMAgGEMYWC6QNC2J9pr+aHh7w/ZNH8QkQYMJkj14u0e2mzER7LqCS9m/dYSpXdenxNQvToAUVjAhgTFzy+0oTcPzKH3uMz+QSsZi1BN7+Vlt2O/JbksFr+y6rs4pfxxWx1MP9+/1dIjKAS8EpBQRMqEq8AZKKQUpaM10DiLpDxdMNyw= --outputDir=SYSout --noretry
