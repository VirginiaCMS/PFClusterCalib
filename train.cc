/* Trainer of semi-parametric MVAs.
 */

#include <vector>

#include <TCut.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooConstVar.h>
#include <RooWorkspace.h>

// GBRLikelihood
#include <RooRevCBExp.h>
#include <RooGausDoubleExp.h>
#include <HybridGBRForest.h>
#include <RooHybridBDTAutoPdf.h>

// prints a message and exits gracefully
#define FATAL(msg) do { fprintf(stderr, "FATAL: %s\n", msg); gSystem->Exit(1); } while (0)

using namespace RooFit;

void train_one(const char* infile, const char* outfile, bool isEE, int pfSize, bool useNumVtx,
               double ptMin = -1, double ptMax = -1)
{
   /* Trains one MVA.
    *
    * pfSize = 1: train only on 1x1 PFClusters;
    * pfSize = 2: train only on 1x2 PFClusters;
    * pfSize = any other value: train on all PFClusters, excluding 1x1 and 1x2;
    *
    * useNumVtx = if true, nVtx branch will be used as MVA input;
    *
    * [ptMin, ptMax) = take only events from this particular pfPt region;
    * negative ptMin/ptMax = no lower/upper limit.
    */

   fprintf(stderr, "   %s, pfSize=%i%s, useNumVtx=%i, ptMin=%.1f, ptMax=%.1f: %s ...\n",
          isEE ? "EE" : "EB", pfSize, pfSize > 2 ? "+" : " ", (int)useNumVtx, ptMin, ptMax, infile);

   // input variables + target variable
   RooArgList allvars;

   allvars.addOwned(*new RooRealVar("var1", "pfE",                0));
   allvars.addOwned(*new RooRealVar("var2", "pfIEtaIX",           0));
   allvars.addOwned(*new RooRealVar("var3", "pfIPhiIY",           0));

//    if (pfSize != 1)
//       allvars.addOwned(*new RooRealVar("var4", "pfE1x3/pfE",    0));
//
//    if (pfSize != 1 && pfSize != 2) {
//       allvars.addOwned(*new RooRealVar("var5", "pfE2x2/pfE",    0));
//       allvars.addOwned(*new RooRealVar("var6", "pfE2x5Max/pfE", 0));
//       allvars.addOwned(*new RooRealVar("var7", "pfE3x3/pfE",    0));
//       allvars.addOwned(*new RooRealVar("var8", "pfE5x5/pfE",    0));
//    }

   if (useNumVtx)
      allvars.addOwned(*new RooRealVar("nVtx", "nVtx", 0));

   if (isEE) {
      allvars.addOwned(*new RooRealVar("varEE1", "ps1E/pfE", 0));
      allvars.addOwned(*new RooRealVar("varEE2", "ps2E/pfE", 0));
   }

   // input variables only
   RooArgList invars = allvars;

   // target variable
   // NOTE: preshower energy is not subtracted
   // NOTE: limits were evaluated with draw_inputs.py
   RooRealVar* target = new RooRealVar("target", "log(mcE/pfE)", 0., -0.336, 0.916);
   allvars.addOwned(*target);

   // variables corresponding to regressed parameters
   RooRealVar mean("mean", "", 0.);
   RooRealVar sigma("sigma", "", 0.1);
   RooRealVar alphaL("alphaL", "", 1.2);
   RooRealVar alphaR("alphaR", "", 2.0);
   RooRealVar powerR("powerR", "", 5);

   mean.setConstant(false);
   sigma.setConstant(false);
   alphaL.setConstant(false);
   alphaR.setConstant(false);
   powerR.setConstant(false);

   // non-parametric functions for each regressed parameter
   RooGBRFunctionFlex funcMean("funcMean", "");
   RooGBRFunctionFlex funcSigma("funcSigma", "");
   RooGBRFunctionFlex funcAlphaL("funcAlphaL", "");
   RooGBRFunctionFlex funcAlphaR("funcAlphaR", "");
   RooGBRFunctionFlex funcPowerR("funcPowerR", "");

   // mapping of input variables to non-parametric functions
   RooGBRTargetFlex tgtMean("tgtMean", "", funcMean, mean, invars);
   RooGBRTargetFlex tgtSigma("tgtSigma", "", funcSigma, sigma, invars);
   RooGBRTargetFlex tgtAlphaL("tgtAlphaL", "", funcAlphaL, alphaL, invars);
   RooGBRTargetFlex tgtAlphaR("tgtAlphaR", "", funcAlphaR, alphaR, invars);
   RooGBRTargetFlex tgtPowerR("tgtPowerR", "", funcPowerR, powerR, invars);

   // parameters' bounds
   RooRealConstraint limMean("limMean", "", tgtMean, -0.336, 0.916);
   RooRealConstraint limSigma("limSigma", "", tgtSigma, 0.001, 0.4);
   RooRealConstraint limAlphaL("limAlphaL", "", tgtAlphaL, 0.2, 7.);
   RooRealConstraint limAlphaR("limAlphaR", "", tgtAlphaR, 0.2, 7.);
   RooRealConstraint limPowerR("limPowerR", "", tgtPowerR, 1.01, 100.);

   // Gaussian + left exponential tail + right power-law or exponential tail
   RooAbsPdf* pdf;
   if (pfSize == 1 || pfSize == 2)
      pdf = new RooGausDoubleExp("pdfGausDoubleExp", "", *target, limMean, limSigma, limAlphaL, limAlphaR);
   else
      pdf = new RooRevCBExp("pdfRevCBExp", "", *target, limMean, limSigma, limAlphaL, limAlphaR, limPowerR);

   // list of mapped functions to regress
   RooArgList tgts;
   tgts.add(tgtMean);
   tgts.add(tgtSigma);
   tgts.add(tgtAlphaL);
   tgts.add(tgtAlphaR);

   if (pfSize != 1 && pfSize != 2)
      tgts.add(tgtPowerR);

   // list of pdfs
   std::vector<RooAbsReal*> pdfs;
   pdfs.push_back(pdf);

   // open file and get tree with the inputs and the target
   TFile* fi = TFile::Open(infile);
   if (!fi || fi->IsZombie())
      FATAL("TFile::Open() failed");

   TTree* tree = dynamic_cast<TTree*>(fi->Get("ntuplizer/PFClusterTree"));
   if (!tree) FATAL("TFile::Get() failed");

   // create a memory-resident friend TTree with linear event numbers
   if (!gROOT->cd()) FATAL("TROOT::cd() failed");
   TTree evtree("ntuplizer/PFClusterTree", "Trivial event numbers");
   evtree.SetAutoFlush(0);
   evtree.SetAutoSave(0);
   Long64_t event;
   evtree.Branch("event", &event);
   for (event = 0; event < tree->GetEntriesFast(); event++)
      evtree.Fill();
   tree->AddFriend(&evtree);

   // pre-filtering cuts
   TCut cuts = (isEE ? "abs(pfEta) > 1.479" : "abs(pfEta) < 1.479");
   cuts += "pfE/mcE > 0.4";      // NOTE: evaluated with draw_inputs.py
   cuts += "pfPhoDeltaR < 0.03"; // NOTE: evaluated with draw_inputs.py
   cuts += "event % 2 == 0";     // NOTE: take only even tree entries

   if (pfSize == 1)
      cuts += "pfSize5x5_ZS == 1";
   else if (pfSize == 2)
      cuts += "pfSize5x5_ZS == 2";
   else
      cuts += "pfSize5x5_ZS >= 3";

   if (ptMin > 0)
      cuts += TString::Format("pfPt >= %f", ptMin);
   if (ptMax > 0)
      cuts += TString::Format("pfPt < %f", ptMax);

   // per-event weight
   // NOTE: title is used for per-event weights and selection cuts
   RooRealVar weightvar("weightvar", "", 1.);
   weightvar.SetTitle(cuts);

   // list of training datasets
   RooDataSet* dataset = RooTreeConvert::CreateDataSet("data", tree, allvars, weightvar);
   std::vector<RooAbsData*> datasets;
   datasets.push_back(dataset);

   // minimum event weight per tree
   std::vector<double> minweights;
   minweights.push_back(200);

   // dummies
   RooConstVar etermconst("etermconst", "", 0.);
   RooRealVar r("r", "", 1.);
   r.setConstant(true);

   // training
   RooHybridBDTAutoPdf bdtpdfdiff("bdtpdfdiff", "", tgts, etermconst, r, datasets, pdfs);
   if (pfSize == 1 || pfSize == 2)
      bdtpdfdiff.SetMinCutSignificance(1.);
   else
      bdtpdfdiff.SetMinCutSignificance(5.);
   //bdtpdfdiff.SetPrescaleInit(100);
   bdtpdfdiff.SetShrinkage(0.1);
   bdtpdfdiff.SetMinWeights(minweights);
   bdtpdfdiff.SetMaxNodes(750);
   bdtpdfdiff.TrainForest(1e+6); // NOTE: valid training will stop at ~100-500 trees

   // unique name of output workspace
   TString wsname = TString::Format("ws_mva_%s_pfSize%i", isEE ? "EE" : "EB", pfSize);
   if (ptMin > -0.5)
      wsname += TString::Format("_ptMin%.1f", ptMin);
   if (ptMax > -0.5)
      wsname += TString::Format("_ptMax%.1f", ptMax);

   // save output to file
   RooWorkspace* ws = new RooWorkspace(wsname);
   ws->import(*pdf);
   ws->writeToFile(outfile, false); // false = update output file, not recreate

   // NOTE: no memory cleanup for simplicity
}

void train(const char* infile, const char* outfile, bool useNumVtx)
{
   // Steering function.

   // EB vs EE
   for (int i = 0; i < 2; i++) {
      bool isEE = (i == 0 ? false : true);

      train_one(infile, outfile, isEE, 1, useNumVtx);
      train_one(infile, outfile, isEE, 2, useNumVtx);

      train_one(infile, outfile, isEE, 3, useNumVtx, 0, 5);
      train_one(infile, outfile, isEE, 3, useNumVtx, 4, 20);
      train_one(infile, outfile, isEE, 3, useNumVtx, 16, -1);
   }
}
