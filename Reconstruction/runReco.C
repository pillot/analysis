void runReco(const char *filename="raw.root")
{
  /////////////////////////////////////////////////////////////////////////////////////////
  //
  // Reconstruction script for RAW data
  //
  /////////////////////////////////////////////////////////////////////////////////////////

  gROOT->Macro(TString::Format("rec.C(\"%s\")", filename));

  Info("runReco", ">>>>>>> Merging TreeR");
  gROOT->Macro("MergeTreeR.C");

  Info("runReco", ">>>>>>> Merging TreeD");
  gROOT->Macro("MergeTreeD.C");
}
