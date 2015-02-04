/* Maker of TTree's friend with outputs from semi-parametric MVAs.
 */

#include <vector>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

// prints a message and exits gracefully
#define FATAL(msg) do { fprintf(stderr, "FATAL: %s\n", msg); gSystem->Exit(1); } while (0)

void eval_one(const char* infile, const char* outfile, std::vector<std::string> fnames)
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
   float pfE, pfEta, pfPhi;
   float pfE1x3, pfE2x2, pfE2x5Max, pfE3x3, pfE5x5;
   float ps1E, ps2E;

   // associate tree branches with variables
   intree->SetBranchAddress("pfE",   &pfE);
   intree->SetBranchAddress("pfEta", &pfEta);
   intree->SetBranchAddress("pfPhi", &pfPhi);

   intree->SetBranchAddress("pfSize5x5_ZS", &pfSize5x5_ZS);

   intree->SetBranchAddress("pfE1x3",    &pfE1x3);
   intree->SetBranchAddress("pfE2x2",    &pfE2x2);
   intree->SetBranchAddress("pfE2x5Max", &pfE2x5Max);
   intree->SetBranchAddress("pfE3x3",    &pfE3x3);
   intree->SetBranchAddress("pfE5x5",    &pfE5x5);

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
   RooRealVar* invar[99][2][3][11];  // [mva number][EB or EE][pfSize][varnum]
   RooAbsReal* mvaMean[99][2][3];
   RooAbsReal* mvaSigma[99][2][3];
   RooAbsReal* mvaAlphaL[99][2][3];
   RooAbsReal* mvaAlphaR[99][2][3];
   RooAbsReal* mvaPowerR[99][2][3];

   // get trainings
   for (size_t i = 0; i < nent; i++) {
      TFile f(Form("output/training_results_%s.root", fnames[i].c_str()));
      if (f.IsZombie()) FATAL("TFile::Open() failed");

      for (int iBE = 0; iBE < 2; iBE++) // barrel vs endcaps
         for (int iS = 0; iS < 3; iS++) {  // pfSize = 1 vs 2 vs 3 and bigger
            const char* det = (iBE == 0 ? "EB" : "EE");

            RooWorkspace* ws = dynamic_cast<RooWorkspace*>(f.Get(Form("ws_mva_%s_pfSize%i", det, iS + 1)));
            if (!ws) FATAL("TFile::Get() failed");

            invar[i][iBE][iS][0] = ws->var("var1");       // pfE
            invar[i][iBE][iS][1] = ws->var("var2");       // pfEta
            invar[i][iBE][iS][2] = ws->var("var3");       // pfPhi

            if (iS > 0)
               invar[i][iBE][iS][3] = ws->var("var4");    // pfE1x3/pfE

            if (iS > 1) {
               invar[i][iBE][iS][4] = ws->var("var5");    // pfE2x2/pfE
               invar[i][iBE][iS][5] = ws->var("var6");    // pfE2x5Max/pfE
               invar[i][iBE][iS][6] = ws->var("var7");    // pfE3x3/pfE
               invar[i][iBE][iS][7] = ws->var("var8");    // pfE5x5/pfE
            }

            invar[i][iBE][iS][8] = ws->var("nVtx");       // nVtx

            if (iBE == 1) {
               invar[i][iBE][iS][9] = ws->var("varEE1");  // ps1E/pfE
               invar[i][iBE][iS][10] = ws->var("varEE2"); // ps2E/pfE
            }

            mvaMean[i][iBE][iS] = ws->function("limMean");
            mvaSigma[i][iBE][iS] = ws->function("limSigma");

            if (iS < 2) { // 1x1 and 1x2
               mvaAlphaL[i][iBE][iS] = ws->function("limAlphaL");
               mvaAlphaR[i][iBE][iS] = ws->function("limAlphaR");
            }
            else {
               mvaAlphaR[i][iBE][iS] = ws->function("limAlphaR");
               mvaPowerR[i][iBE][iS] = ws->function("limPowerR");
            }
         }
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

      for (size_t i = 0; i < nent; i++) {
         *invar[i][iBE][iS][0] = pfE;
         *invar[i][iBE][iS][1] = pfEta;
         *invar[i][iBE][iS][2] = pfPhi;

         if (iS > 0)
            *invar[i][iBE][iS][3] = pfE1x3/pfE;

         if (iS > 1) {
            *invar[i][iBE][iS][4] = pfE2x2/pfE;
            *invar[i][iBE][iS][5] = pfE2x5Max/pfE;
            *invar[i][iBE][iS][6] = pfE3x3/pfE;
            *invar[i][iBE][iS][7] = pfE5x5/pfE;
         }

         // NULL if not available
         if (invar[i][iBE][iS][8])
            *invar[i][iBE][iS][8] = nVtx;

         if (iBE == 1) {
            *invar[i][iBE][iS][9] = ps1E/pfE;
            *invar[i][iBE][iS][10] = ps2E/pfE;
         }

         mean[i] = mvaMean[i][iBE][iS]->getVal();
         sigma[i] = mvaSigma[i][iBE][iS]->getVal();
         alphaL[i] = 0;
         alphaR[i] = 0;
         powerR[i] = 0;

         if (iS < 2) { // 1x1 and 1x2
            alphaL[i] = mvaAlphaL[i][iBE][iS]->getVal();
            alphaR[i] = mvaAlphaR[i][iBE][iS]->getVal();
         }
         else {
            alphaR[i] = mvaAlphaR[i][iBE][iS]->getVal();
            powerR[i] = mvaPowerR[i][iBE][iS]->getVal();
         }
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
