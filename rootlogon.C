
void rootlogon()
{
   /* This function must be called by root on start for train.cc and eval.cc be compilable.
    */

   // keep current working directory clean from .d and .so files
   gSystem->SetBuildDir("output", true);

   if (gSystem->Getenv("CMSSW_BASE")) {
      gSystem->AddIncludePath("-I$ROOFITSYS/include");
      gSystem->AddIncludePath("-I$CMSSW_BASE/src/HiggsAnalysis/GBRLikelihood/interface");

      gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisGBRLikelihood.so");
   }
}
