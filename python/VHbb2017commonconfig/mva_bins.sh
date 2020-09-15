#!/bin/bash

# example: VHbb2017commonconfig/mva_bins.sh "arctans(15)"
# arctan shape binning with S(number of signal events) bins, max 15

REBINMETHOD=$1

echo "bin boundaries for method ${REBINMETHOD}"
echo ""
echo "dc:SR_med_Znn_0j.rebin_list="`./submit.py -J runplot -T Zvv2017 -F tmp --regions SR_med_Znn_0j --vars DNN --set="plotDef:DNN.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`
echo "dc:SR_med_Znn_ge1j.rebin_list="`./submit.py -J runplot -T Zvv2017 -F tmp --regions SR_med_Znn_ge1j --vars DNN --set="plotDef:DNN.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`
echo "dc:SR_high_Znn.rebin_list="`./submit.py -J runplot -T Zvv2017 -F tmp --regions SR_high_Znn --vars DNN --set="plotDef:DNN.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`
echo "dc:SR_high_Znn_BOOST.rebin_list="`./submit.py -J runplot -T Zvv2017 -F tmp --regions SR_high_Znn_BOOST --vars BDT_Zvv_BOOSTFinal_wdB --set="plotDef:BDT_Zvv_BOOSTFinal_wdB.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`

echo "dc:SR_med_Wmn.rebin_list="`./submit.py -J runplot -T Wlv2017 -F tmp --regions SR_med_Wmn --vars DNN --set="plotDef:DNN.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`
echo "dc:SR_med_Wen.rebin_list="`./submit.py -J runplot -T Wlv2017 -F tmp --regions SR_med_Wen --vars DNN --set="plotDef:DNN.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`
echo "dc:SR_high_Wmn.rebin_list="`./submit.py -J runplot -T Wlv2017 -F tmp --regions SR_high_Wmn --vars DNN --set="plotDef:DNN.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`
echo "dc:SR_high_Wen.rebin_list="`./submit.py -J runplot -T Wlv2017 -F tmp --regions SR_high_Wen --vars DNN --set="plotDef:DNN.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`
echo "dc:SR_high_Wmn_BOOST.rebin_list="`./submit.py -J runplot -T Wlv2017 -F tmp --regions SR_high_Wmn_BOOST --vars BDT_Wlv_BOOSTFinal_wdB --set="plotDef:BDT_Wlv_BOOSTFinal_wdB.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`
echo "dc:SR_high_Wen_BOOST.rebin_list="`./submit.py -J runplot -T Wlv2017 -F tmp --regions SR_high_Wen_BOOST --vars BDT_Wlv_BOOSTFinal_wdB --set="plotDef:BDT_Wlv_BOOSTFinal_wdB.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`

echo "dc:SR_low_Zmm.rebin_list="`./submit.py -J runplot -T Zll2017 -F tmp --regions SR_low_Zmm --vars DNNlow --set="plotDef:DNNlow.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`
echo "dc:SR_low_Zee.rebin_list="`./submit.py -J runplot -T Zll2017 -F tmp --regions SR_low_Zee --vars DNNlow --set="plotDef:DNNlow.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`
echo "dc:SR_med_Zmm_0j.rebin_list="`./submit.py -J runplot -T Zll2017 -F tmp --regions SR_med_Zmm_0j --vars DNN --set="plotDef:DNN.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`
echo "dc:SR_med_Zee_0j.rebin_list="`./submit.py -J runplot -T Zll2017 -F tmp --regions SR_med_Zee_0j --vars DNN --set="plotDef:DNN.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`
echo "dc:SR_med_Zmm_ge1j.rebin_list="`./submit.py -J runplot -T Zll2017 -F tmp --regions SR_med_Zmm_ge1j --vars DNN --set="plotDef:DNN.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`
echo "dc:SR_med_Zee_ge1j.rebin_list="`./submit.py -J runplot -T Zll2017 -F tmp --regions SR_med_Zee_ge1j --vars DNN --set="plotDef:DNN.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`
echo "dc:SR_high_Zmm.rebin_list="`./submit.py -J runplot -T Zll2017 -F tmp --regions SR_high_Zmm --vars DNN --set="plotDef:DNN.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`
echo "dc:SR_high_Zee.rebin_list="`./submit.py -J runplot -T Zll2017 -F tmp --regions SR_high_Zee --vars DNN --set="plotDef:DNN.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`
echo "dc:SR_high_Zmm_BOOST.rebin_list="`./submit.py -J runplot -T Zll2017 -F tmp --regions SR_high_Zmm_BOOST --vars BDT_Zll_BOOSTFinal_wdB --set="plotDef:BDT_Zll_BOOSTFinal_wdB.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`
echo "dc:SR_high_Zee_BOOST.rebin_list="`./submit.py -J runplot -T Zll2017 -F tmp --regions SR_high_Zee_BOOST --vars BDT_Zll_BOOSTFinal_wdB --set="plotDef:BDT_Zll_BOOSTFinal_wdB.rebinMethod=${REBINMETHOD}" --local 2>1 | grep -A1 "new bin boundaries" | tail -n1 | sed 's/INFO:  //g'`
