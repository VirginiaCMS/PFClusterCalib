#!/bin/bash
#
# Run all PFCluster calibration steps in one shot.
# This script also serves as a source of documentation.
#
# Input root files must be produced with ncuAnalysis/PFClusterCalib CMSSW
# module and are expected to be found inside input/.
#
# NOTE: files input/ntuple_photongun_*.root were obtained with CMSSW_7_2_3_patch1.
# NOTE: in this script, there is a place with hardcoded "ntuple_photongun_nopu".
#
# NOTE: in this script, the execution of rootlogon.C is enforced.
#

# stop on first error
set -e

# draw/save distributions of input variables + deltaR + pfE/mcE + MC truth variables
python draw_inputs.py &

# draw/save fractions of various PFCluster sizes vs pT
python auxiliary/draw_pfsize.py &

# draw/save evolution of fitting function parameters with pT
python auxiliary/draw_fit_params.py &

wait

ntuples=`ls input/*.root`

mkdir -p output

echo "Training semi-parametric MVAs with GBRLikelihood:"

for infile in $ntuples; do
    # extract file name without extension
    fname="${infile##*/}"
    fname="${fname%.root}"

    # remove previous result, if any
    rm -f output/training_results_${fname}.root

    # do not use nVtx as input for the no-pileup MC
    useNumVtx=true
    [ "${fname%_rereco}" == "ntuple_photongun_nopu" ] && useNumVtx=false

    for pfSize in 1 2 3; do
        suff=" "
        [ "$pfSize" == 3 ] && suff="+"

        # ECAL barrel
        echo "   EB, pfSize=${pfSize}${suff}, useNumVtx=${useNumVtx}: ${infile} ..."
        echo "
            .x rootlogon.C
            .x train_one.cc+(\"${infile}\", \"output/training_results_${fname}.root\", \"ws_mva_EB_pfSize${pfSize}\", false, ${pfSize}, ${useNumVtx})
            .q" | root -b -l &>output/train_${fname}_EB_pfSize${pfSize}.log

        # ECAL endcaps
        echo "   EE, pfSize=${pfSize}${suff}, useNumVtx=${useNumVtx}: ${infile} ..."
        echo "
            .x rootlogon.C
            .x train_one.cc+(\"${infile}\", \"output/training_results_${fname}.root\", \"ws_mva_EE_pfSize${pfSize}\", true, ${pfSize}, ${useNumVtx})
            .q" | root -b -l &>output/train_${fname}_EE_pfSize${pfSize}.log
    done
done

echo "Evaluating outputs from semi-parametric MVAs:"

# compile eval_one.cc
echo "
   .x rootlogon.C
   .L eval_one.cc+
   .q" | root -b -l

# prepare common CINT commands for the next "for" block
cmd="vector<string> fnames"$'\n'
for name in $ntuples; do
    name="${name##*/}"
    cmd="${cmd}fnames.push_back(\"${name%.root}\")"$'\n'
done

for infile in $ntuples; do
    # extract file name without extension
    fname="${infile##*/}"
    fname="${fname%.root}"

    echo "   ${infile} ..."

    echo "
        .x rootlogon.C
        $cmd
        .x eval_one.cc+(\"${infile}\", \"output/friend_${fname}.root\", fnames)
        .q" | root -b -l &

done

wait

# draw/save distributions with achieved energy resolutions
python draw_results.py &

# draw/save some slices with achieved energy resolutions
python draw_slices.py &

# draw/save some train vs test slices
python draw_overtraining.py &

wait

echo "Finished."
