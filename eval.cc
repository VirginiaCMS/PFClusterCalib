/* Maker of TTree's friend with outputs from semi-parametric MVAs.
 */

#include <vector>
#include <string>

#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

// prints a message and exits gracefully
#define FATAL(msg) do { fprintf(stderr, "FATAL: %s\n", msg); gSystem->Exit(1); } while (0)

void eval(const char* infile, const char* outfile, std::vector<std::string> fnames)
{
   /* Main function.
    *
    * fnames = array with names of input ntuples.
    */

   // open file and get TTree with the inputs
   TFile* fi = TFile::Open(infile);
   if (!fi || fi->IsZombie())
      FATAL("TFile::Open() failed");

   TTree* intree = dynamic_cast<TTree*>(fi->Get("ntuplizer/PFClusterTree"));
   if (!intree) FATAL("TFile::Get() failed");

   // variables to be associated with the input tree branches
   Int_t nVtx;
   Int_t pfSize5x5_ZS;
   Int_t pfIEtaIX, pfIPhiIY;
   float pfE, pfPt, pfEta; //, pfPhi;
//    float pfE1x3, pfE2x2, pfE2x5Max, pfE3x3, pfE5x5;
   float ps1E, ps2E;

   // associate tree branches with variables
   intree->SetBranchAddress("pfE",   &pfE);
   intree->SetBranchAddress("pfPt",  &pfPt);
   intree->SetBranchAddress("pfEta", &pfEta);
//    intree->SetBranchAddress("pfPhi", &pfPhi);
   intree->SetBranchAddress("pfIEtaIX", &pfIEtaIX);
   intree->SetBranchAddress("pfIPhiIY", &pfIPhiIY);

   intree->SetBranchAddress("pfSize5x5_ZS", &pfSize5x5_ZS);

//    intree->SetBranchAddress("pfE1x3",    &pfE1x3);
//    intree->SetBranchAddress("pfE2x2",    &pfE2x2);
//    intree->SetBranchAddress("pfE2x5Max", &pfE2x5Max);
//    intree->SetBranchAddress("pfE3x3",    &pfE3x3);
//    intree->SetBranchAddress("pfE5x5",    &pfE5x5);

   intree->SetBranchAddress("nVtx", &nVtx);

   intree->SetBranchAddress("ps1E", &ps1E);
   intree->SetBranchAddress("ps2E", &ps2E);

   // number of ntuples
   size_t nent = fnames.size();
   if (nent < 1 || nent > 99) FATAL("fnames.size() not in range 1-99");

   // prepare output tree
   TFile* fo = TFile::Open(outfile, "RECREATE");
   if (!fo || fo->IsZombie())
      FATAL("TFile::Open() failed");

   TDirectory* dir = fo->mkdir("ntuplizer");
   if (!dir) FATAL("TFile::mkdir() failed");
   if (!dir->cd()) FATAL("TDirectory::cd() failed");

   TTree* outtree = new TTree("PFClusterTree", "Outputs from semi-parametric MVAs");

   // array of variables to be associated with the output tree branches
   float mean[99], sigma[99], alphaL[99], alphaR[99], powerR[99];

   // associate variables with the output tree branches
   for (size_t i = 0; i < nent; i++) {
      outtree->Branch(Form("mva_mean_%s", fnames[i].c_str()), &mean[i]);
      outtree->Branch(Form("mva_sigma_%s", fnames[i].c_str()), &sigma[i]);
      outtree->Branch(Form("mva_alphaL_%s", fnames[i].c_str()), &alphaL[i]);
      outtree->Branch(Form("mva_alphaR_%s", fnames[i].c_str()), &alphaR[i]);
      outtree->Branch(Form("mva_powerR_%s", fnames[i].c_str()), &powerR[i]);
   }

   // semi-parametric MVAs' inputs and outputs
   RooRealVar* invar[99][2][9][33];  // [mva number][EB or EE][pfSize][varnum]
   RooAbsReal* mvaMean[99][2][9];
   RooAbsReal* mvaSigma[99][2][9];
   RooAbsReal* mvaAlphaL[99][2][9];
   RooAbsReal* mvaAlphaR[99][2][9];
   RooAbsReal* mvaPowerR[99][2][9];

   // get trainings
   for (size_t i = 0; i < nent; i++) {
      TFile f(Form("output/training_results_%s.root", fnames[i].c_str()));
      if (f.IsZombie()) FATAL("TFile::Open() failed");

      for (int iBE = 0; iBE < 2; iBE++) // barrel vs endcaps
         for (int iS = 0; iS < 5; iS++) {  // pfSize = 1 vs 2 vs 3 and bigger (pfPt-sliced)
            const char* det = (iBE == 0 ? "EB" : "EE");

            int pfSize = (iS < 2 ? iS + 1 : 3);
            TString wsname = TString::Format("ws_mva_%s_pfSize%i", det, pfSize);

            // pfPt slices
            double ptMin = -1, ptMax = -1;
            if (iS == 2) {
               ptMin = 0;
               ptMax = 5;
            } else if (iS == 3) {
               ptMin = 4;
               ptMax = 20;
            } else if (iS == 4) {
               ptMin = 16;
               ptMax = -1;
            }

            if (ptMin > -0.5)
               wsname += TString::Format("_ptMin%.1f", ptMin);
            if (ptMax > -0.5)
               wsname += TString::Format("_ptMax%.1f", ptMax);

            RooWorkspace* ws = dynamic_cast<RooWorkspace*>(f.Get(wsname));
            if (!ws) FATAL("TFile::Get() failed");

            invar[i][iBE][iS][0] = ws->var("var1");       // pfE
            invar[i][iBE][iS][1] = ws->var("var2");       // pfIEtaIX
            invar[i][iBE][iS][2] = ws->var("var3");       // pfIPhiIY

//             if (iS > 0)
//                invar[i][iBE][iS][3] = ws->var("var4");    // pfE1x3/pfE
//
//             if (iS > 1) {
//                invar[i][iBE][iS][4] = ws->var("var5");    // pfE2x2/pfE
//                invar[i][iBE][iS][5] = ws->var("var6");    // pfE2x5Max/pfE
//                invar[i][iBE][iS][6] = ws->var("var7");    // pfE3x3/pfE
//                invar[i][iBE][iS][7] = ws->var("var8");    // pfE5x5/pfE
//             }

            invar[i][iBE][iS][8] = ws->var("nVtx");       // nVtx

            if (iBE == 1) {
               invar[i][iBE][iS][9] = ws->var("varEE1");  // ps1E/pfE
               invar[i][iBE][iS][10] = ws->var("varEE2"); // ps2E/pfE
            }

            mvaMean[i][iBE][iS] = ws->function("limMean");
            mvaSigma[i][iBE][iS] = ws->function("limSigma");
            mvaAlphaL[i][iBE][iS] = ws->function("limAlphaL");
            mvaAlphaR[i][iBE][iS] = ws->function("limAlphaR");

            if (iS > 1)
               mvaPowerR[i][iBE][iS] = ws->function("limPowerR");
         } // iS loop
   } // training number

   // loop over events
   for (Long64_t ev = 0; ev < intree->GetEntriesFast(); ev++) {
      if (intree->GetEntry(ev) <= 0)
         FATAL("TTree::GetEntry() failed");

      // 0=ECAL barrel vs 1=ECAL endcaps
      int iBE = (fabs(pfEta) < 1.479) ? 0 : 1;

      if (pfSize5x5_ZS <= 0) FATAL("pfSize5x5_ZS <= 0");

      // pfSize category
      int iS = (pfSize5x5_ZS > 2 ? 2 : pfSize5x5_ZS - 1);

      // pfPt slice category
      if (iS == 2) {
         if (pfPt >= 4.5 && pfPt < 18) iS = 3;
         else if (pfPt >= 18) iS = 4;
      }

      for (size_t i = 0; i < nent; i++) {
         *invar[i][iBE][iS][0] = pfE;
         *invar[i][iBE][iS][1] = pfIEtaIX;
         *invar[i][iBE][iS][2] = pfIPhiIY;

//          if (iS > 0)
//             *invar[i][iBE][iS][3] = pfE1x3/pfE;
//
//          if (iS > 1) {
//             *invar[i][iBE][iS][4] = pfE2x2/pfE;
//             *invar[i][iBE][iS][5] = pfE2x5Max/pfE;
//             *invar[i][iBE][iS][6] = pfE3x3/pfE;
//             *invar[i][iBE][iS][7] = pfE5x5/pfE;
//          }

         // NULL if not available
         if (invar[i][iBE][iS][8])
            *invar[i][iBE][iS][8] = nVtx;

         if (iBE == 1) {
            *invar[i][iBE][iS][9] = ps1E/pfE;
            *invar[i][iBE][iS][10] = ps2E/pfE;
         }

         mean[i] = TMath::Exp(mvaMean[i][iBE][iS]->getVal());
         sigma[i] = mvaSigma[i][iBE][iS]->getVal();
         alphaL[i] = mvaAlphaL[i][iBE][iS]->getVal();
         alphaR[i] = mvaAlphaR[i][iBE][iS]->getVal();
         powerR[i] = 0;

         if (iS > 1)
            powerR[i] = mvaPowerR[i][iBE][iS]->getVal();
      }

      outtree->Fill();

   } // event loop

   // flush caches
   if (!dir->cd()) FATAL("TDirectory::cd() failed");
   outtree->Write("", TObject::kOverwrite);

   // cleanup
   delete intree;
   delete outtree;
   delete fi;
   delete fo;
}
