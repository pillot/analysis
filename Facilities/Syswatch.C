void Syswatch(const char* file1="syswatch.log", 
	      const char* what="pI.fMemResident:id2",
	      const char* condition="id0==7 && id2>0")
{
  // for muon reco, use what="pI.fMemResident:id2" and condition="id0==7 && id2>0"
  // for filtering train, use what="pI.fMemResident:id0" and condition="id1==1001"
  // for analysis what="pI.fMemResident:id0" and condition="id1==0" (id1==0 for 1st analysis in the train)
  //
  // cpu times : what=(pI.fCpuUser-pIOld.fCpuUser)/(stampSec-stampOldSec):id2
  
  TTree t1;
  Double_t xmin(30);
  Double_t xmax(2000);
  Double_t ymax(5000);

  new TCanvas;

  t1.ReadFile(gSystem->ExpandPathName(file1));

  t1.SetLineColor(1);
  /*
  TH2* h = new TH2F("hframe","hframe",100,0,xmax,100,0,ymax);
  
  h->Draw("");
  
  t1.Draw(what,condition,"LPSAME");
  */
  t1.Draw(what,condition,"LP");
  
  new TCanvas;
  
  t1.Draw("pI.fMemVirtual:id2",condition,"LP");

//  TGraph *g1 = (TGraph*)gPad->GetPrimitive("Graph");
//
//  g1->SetName("g1");
//
//  g1->Fit("pol1","","",xmin,1000);

  new TCanvas;

//  t1.Draw("stampSec-stampOldSec:id2",condition,"PROF");
  t1.Draw("(pI.fCpuUser-pIOld.fCpuUser)/(stampSec-stampOldSec):id2",condition,"PROF");

//  new TCanvas;
//  
//  t1.Draw("pI.fMemResident-pIOld.fMemResident:id2",condition,"PROF");
  
}
