/* Helps draw_mva_pars.py to perform CPU-intensive tasks: fills and fits
 * histograms.
 */

#include <vector>

#include <TF1.h>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TString.h>
#include <TGraphErrors.h>

// prints a message and exits gracefully
#define FATAL(msg) do { fprintf(stderr, "FATAL: %s\n", msg); gSystem->Exit(1); } while (0)

using namespace std;

// global variables
vector<float> gDataMcPt;     // mcPt
vector<float> gDataPfEta;    // pfEta
vector<float> gDataResol;    // pfE/mcE * [MVA's mean]
vector<float> gDataExpWidth; // [MVA's sigma]/[MVA's mean]

TGraphErrors* grMeanVsMean = NULL;
TGraphErrors* grSigmaVsSigma = NULL;

//______________________________________________________________________________
void SetBranchAddress(TTree* tree, const char* bname, void* ptr)
{
   /* Associates tree branch with address of a variable.
    */

   // verify branch existence
   if (!tree->GetBranch(bname))
      FATAL(Form("tree branch \"%s\" does not exist", bname));

   // activate this branch
   tree->SetBranchStatus(bname, 1);

   // associate tree branch with address and check return value
   Int_t ret = tree->SetBranchAddress(bname, ptr);
   if (ret != 0 && ret != 4)
      FATAL("TTree::SetBranchAddress() returned bad code");
}

//______________________________________________________________________________
void fill_arrays(const char* infile, const char* friendname, const char* mva_name)
{
   // Fills global variables-arrays.

   // open root file
   TFile* fi = TFile::Open(infile);
   if (!fi || fi->IsZombie())
      FATAL("TFile::Open() failed");

   // get TTree with PFClusters
   TTree* tree = dynamic_cast<TTree*>(fi->Get("ntuplizer/PFClusterTree"));
   if (!tree) FATAL("TFile::Get() failed");

   // add branches with outputs from MVAs
   if (!tree->AddFriend("ntuplizer/PFClusterTree", friendname))
      FATAL("TTree::AddFriend() failed");

   // disable all branches by default
   tree->SetBranchStatus("*", 0);

   // variables to be associated with the input tree branches
   float mcE, mcPt, pfE, pfEta;
   float mva_mean, mva_sigma;

   // associate tree branches with variables
   SetBranchAddress(tree, "mcE",   &mcE);
   SetBranchAddress(tree, "mcPt",  &mcPt);
   SetBranchAddress(tree, "pfE",   &pfE);
   SetBranchAddress(tree, "pfEta", &pfEta);
   SetBranchAddress(tree, Form("mva_mean_%s", mva_name), &mva_mean);
   SetBranchAddress(tree, Form("mva_sigma_%s", mva_name), &mva_sigma);

   // cleanup from previous execution
   gDataMcPt.clear();
   gDataPfEta.clear();
   gDataResol.clear();
   gDataExpWidth.clear();

   // loop over events and collect data
   for (Long64_t ev = 1; ev < tree->GetEntriesFast(); ev += 2) {// NOTE: take only test events
      if (tree->GetEntry(ev) <= 0)
         FATAL("TTree::GetEntry() failed");

      gDataMcPt.push_back(mcPt);
      gDataPfEta.push_back(pfEta);
      gDataResol.push_back(pfE/mcE * mva_mean);

      // width of pfE/mcE distribution
      gDataExpWidth.push_back(mva_sigma/mva_mean);
   }
}

//______________________________________________________________________________
void MeanSigma(vector<float> &numbers, double &mean, double &sigma)
{
   /* Evaluates average (mean) and dispersion (sigma) of numbers "numbers".
    *
    * Mean and sigma are recalculated iteratively several times. During each
    * calculation, a region [mean - 3*sigma, mean + 3*sigma] is used, where
    * "mean" and "sigma" are taken from a previous iteration.
    */

   // zero-order iteration
   mean = 0;
   sigma = 0;
   for (size_t i = 0; i < numbers.size(); i++) {
      mean += numbers[i];
      sigma += numbers[i] * numbers[i];
   }
   mean /= numbers.size();
   sigma = sqrt(sigma/numbers.size() - mean*mean);

   // iterations
   for (int c = 0; c < 1000; c++) {
      double mean_prev = mean;
      double sigma_prev = sigma;

      double xmin = mean - 3*sigma;
      double xmax = mean + 3*sigma;

      // evaluate peak position and width
      mean = 0;
      sigma = 0;
      int nent = 0;
      for (size_t i = 0; i < numbers.size(); i++) {
         if (numbers[i] < xmin || numbers[i] > xmax) continue;

         mean += numbers[i];
         sigma += numbers[i] * numbers[i];
         nent++;
      }
      mean /= nent;
      sigma = sqrt(sigma/nent - mean*mean);

      // break when converged
      if (fabs(mean - mean_prev) <= 1e-6 * fabs(mean) &&
          fabs(sigma - sigma_prev) <= 1e-6 * fabs(sigma))
         return;
   }

   const char* fmt = "nent_total=%lli, mean=%f, sigma=%f\n";
   fprintf(stderr, fmt, numbers.size(), mean, sigma);
   FATAL("mean and/or sigma did not converged");
}

//______________________________________________________________________________
void fit_slices_real(vector<float>& x, vector<float>& y, int blockSize,
                     const char* title, const char* xtitle)
{
    /* Fits distributions of sorted blocks of data points, result is given in
     * grMeanVsMean and grSigmaVsSigma.
     *
     * NOTE: sigma = width/position.
     */

   // cleanup from previous execution
   if (grMeanVsMean) delete grMeanVsMean;
   if (grSigmaVsSigma) delete grSigmaVsSigma;

   grMeanVsMean = new TGraphErrors();
   grSigmaVsSigma = new TGraphErrors();

   size_t siz = x.size();
   if (siz < 1) FATAL("x.size() < 1");

   // sort by X axis
   size_t* ind = new size_t[siz]; // allocate on the heap to avoid segfaults
   TMath::Sort(siz, &x.front(), ind, false);

   // NOTE: last block is excluded if it has less than 0.5 * blockSize entries
   int nblocks = TMath::Nint(round(((float)siz)/blockSize));

   TCanvas* c = NULL;
   vector<TObject*> todel(64);

   // loop over blocks of ordered data
   for (int b = 0; b < nblocks; b++) {
      vector<float> bx;
      vector<float> by;

      // fill separate arrays with current block data
      for (size_t i = b * blockSize; i < (size_t) (b + 1) * blockSize; i++) {
         if (i >= siz) break;
         bx.push_back(x[ind[i]]);
         by.push_back(y[ind[i]]);
      }

      double meanX, sigmaX;
      double meanY, sigmaY;

      // mean and sigma in the block and Etrue/Erec
      MeanSigma(bx, meanX, sigmaX);
      MeanSigma(by, meanY, sigmaY);

      // fill histogram
      TH1* h = new TH1D("h", "", 200, 0.65, 1.2);
      for (size_t i = 0; i < by.size(); i++)
         h->Fill(by[i]);

      // create new canvas, if necessary
      if (b % 9 == 0) {
         if (c) {
            c->SaveAs(Form("output/plots_mva_pars/%s.png", c->GetTitle()));

            // memory cleanup
            delete c;
            for (size_t k = 0; k < todel.size(); k++)
               delete todel[k];
            todel.clear();
         }

         TString cname = TString::Format("fits_%s_blk%03dto%03d", title, b + 1, b + 9);
         c = new TCanvas(cname, cname, 1000, 700);

         c->SetLeftMargin(0);
         c->SetRightMargin(0);
         c->SetTopMargin(0);
         c->SetBottomMargin(0);
         c->Divide(3, 3);
      }

      c->cd(b % 9 + 1);
      gPad->SetLeftMargin(0.12);
      gPad->SetRightMargin(0.02);
      gPad->SetTopMargin(0.08);
      gPad->SetBottomMargin(0.08);

      h->SetTitle(Form("%s = (%.4f #pm %.2g)%%", xtitle, meanX * 100, sigmaX * 100));
      h->SetXTitle("E^{rec}/E^{gen}");
      h->SetYTitle("Entries");
      h->SetTitleOffset(1.6, "Y");

      h->Sumw2(true);
      h->SetLineColor(kBlack);
      h->Draw();

      // Gaussian + left power-law tail + right exponential tail
      const char* expr = "[0] * ( (x-[1])/[2] > -[3] ? "
         "( (x-[1])/[2] < [5] ? exp(-(x-[1])^2/(2*[2]*[2])) : exp(0.5*[5]*[5] - [5]*(x-[1])/[2]) ) : "
         "([4]/[3])^[4] * exp(-0.5*[3]^2) * (-(x-[1])/[2]-[3]+[4]/[3])^(-[4]) )";

      TF1* fit = new TF1("fit", expr, 0.65, 1.2);
      fit->SetLineWidth(1);
      fit->SetNpx(2000);

      fit->SetParameters(h->GetMaximum(), meanY, sigmaY, 1.5, 5, 1.5);
      fit->SetParLimits(0, 0.33 * h->GetMaximum(), 3 * h->GetMaximum());
      fit->SetParLimits(1, 0.65, 1.2);
      fit->SetParLimits(2, 0.33 * sigmaY, 1.5 * sigmaY);
      fit->SetParLimits(3, 0, 10);
      fit->SetParLimits(4, 1.01, 100);
      fit->SetParLimits(5, 0, 10);

      h->Fit(fit, "QEM", "same", 0.65, 1.2);  // improves convergence
      h->Fit(fit, "QEML", "same", 0.65, 1.2);

      grMeanVsMean->SetPoint(b, meanX, fit->GetParameter(1));
      grMeanVsMean->SetPointError(b, sigmaX, fit->GetParError(1));

      grSigmaVsSigma->SetPoint(b, meanX, fit->GetParameter(2)/fit->GetParameter(1));
      grSigmaVsSigma->SetPointError(b, sigmaX, fit->GetParError(2)/fit->GetParameter(1));

      todel.push_back(h);
      todel.push_back(fit);
   } // block loop

   // save the very last canvas
   if (c) {
      c->SaveAs(Form("output/plots_mva_pars/%s.png", c->GetTitle()));

      // memory cleanup
      delete c;
      for (size_t k = 0; k < todel.size(); k++)
         delete todel[k];
   }

   delete ind;
}

//______________________________________________________________________________
void fit_slices(int type, int blockSize, const char* title, const char* xtitle,
                double pt1=0, double pt2=0)
{
   // Steers work of fit_slices_real().

   // mcPt region, EB
   if (type == 0) {
      vector<float> x, y;

      for (size_t i = 0; i < gDataMcPt.size(); i++) {
         if (fabs(gDataPfEta[i]) > 1.479) continue;
         if (gDataMcPt[i] < pt1 || gDataMcPt[i] >= pt2) continue;

         x.push_back(gDataExpWidth[i]);
         y.push_back(gDataResol[i]);
      }

      fit_slices_real(x, y, blockSize, title, xtitle);
   }

   // mcPt region, EE
   else if (type == 1) {
      vector<float> x, y;

      for (size_t i = 0; i < gDataMcPt.size(); i++) {
         if (fabs(gDataPfEta[i]) < 1.479) continue;
         if (gDataMcPt[i] < pt1 || gDataMcPt[i] >= pt2) continue;

         x.push_back(gDataExpWidth[i]);
         y.push_back(gDataResol[i]);
      }

      fit_slices_real(x, y, blockSize, title, xtitle);
   }

   else
      FATAL("invalid type");
}
