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

ntuples=`ls /afs/cern.ch/work/k/konush/public/PFClusterCalib_Apr2015/input/*.root`

mkdir -p output

echo "Training semi-parametric MVAs with GBRLikelihood:" 1>&2

for infile in $ntuples; do
    # extract file name without extension
    fname="${infile##*/}"
    fname="${fname%.root}"

    # do not use nVtx as input for the no-pileup MC
    useNumVtx=true
    [ "${fname%_rereco}" == "ntuple_photongun_nopu" ] && useNumVtx=false

    # remove previous result, if any
    rm -f output/training_results_${fname}.root

    echo "
        .x rootlogon.C
        .x train.cc+(\"${infile}\", \"output/training_results_${fname}.root\", ${useNumVtx})
        .q" | root -b -l >output/train_${fname}.log
done

echo "Evaluating outputs from semi-parametric MVAs:" 1>&2

# compile eval_one.cc
echo "
   .x rootlogon.C
   .L eval.cc+
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

    echo "   ${infile} ..." 1>&2

    echo "
        .x rootlogon.C
        $cmd
        .x eval.cc+(\"${infile}\", \"output/friend_${fname}.root\", fnames)
        .q" | root -b -l &

done

wait

# draw/save distributions with achieved energy resolutions
python draw_results.py &

# draw/save some slices with achieved energy resolutions
python draw_slices.py &

# draw/save some train vs test slices
python draw_overtraining.py &

# draw/save real vs estimated energy resolutions.
python draw_mva_pars.py &

wait

echo "Finished."
