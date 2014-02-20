
{
gSystem->Load("libVMC.so");
gSystem->Load("libTree.so");
gSystem->Load("libPhysics.so");
gSystem->Load("libMatrix.so");
gSystem->Load("libMinuit.so");
gSystem->Load("libXMLParser.so");
gSystem->Load("libGui.so");
gSystem->Load("libSTEERBase.so");
gSystem->Load("libESD.so");
gSystem->Load("libAOD.so");
gSystem->Load("libANALYSIS.so");
gSystem->Load("libANALYSISalice.so");
gSystem->Load("libPWG3base.so");

  gSystem->SetIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT/PWG3/base");
  
  gStyle->SetPalette(1);
  gStyle->SetFillColor(0);
}
