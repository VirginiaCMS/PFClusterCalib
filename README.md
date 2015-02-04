# Calibration of PFClusters in ECAL (CMS experiment at the LHC)

Input root files must be produced with ncuAnalysis/PFClusterCalib CMSSW module
and are expected to be found inside input/.

GBRLikelihood module of CMSSW must be installed for this analysis. Instructions:

    cd $CMSSW_BASE/src__
    git clone https://github.com/bendavid/GBRLikelihood.git HiggsAnalysis/GBRLikelihood__
    scram b

To do everything in one shot, execute ./runall.sh.
