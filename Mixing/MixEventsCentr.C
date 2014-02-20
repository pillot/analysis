#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <Riostream.h>

#include <TCanvas.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TLine.h>
#include <TRandom.h>
#include <TList.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TGraphErrors.h>
#include <TParticle.h>
#include <TClonesArray.h>
#endif

Double_t DimuMass(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t); 
Double_t DimuRap(Double_t, Double_t);

void MixEventsCentr(TString dir, Int_t Run, TString cent, Bool_t display = kTRUE){
  
  //---------------------------------------------------------
  // Open File
  //---------------------------------------------------------
  
  TString fname = Form("%s/PoolCentr_%d.root", dir.Data(), Run);
  TFile *f = new TFile(fname.Data());
  TClonesArray *muplus = (TClonesArray*) f->Get(Form("muplus_%s",cent.Data()));
  TClonesArray *muminus = (TClonesArray*) f->Get(Form("muminus_%s",cent.Data()));
  
  TH1D *hDimuPM = (TH1D*) f->Get(Form("hDimuPM_%s",cent.Data()));
  TH1D *hDimuPP = (TH1D*) f->Get(Form("hDimuPP_%s",cent.Data()));
  TH1D *hDimuMM = (TH1D*) f->Get(Form("hDimuMM_%s",cent.Data()));
  
  TH1D *hDimuPtPM = (TH1D*) f->Get(Form("hDimuPtPM_%s",cent.Data()));
  TH1D *hDimuPtPP = (TH1D*) f->Get(Form("hDimuPtPP_%s",cent.Data()));
  TH1D *hDimuPtMM = (TH1D*) f->Get(Form("hDimuPtMM_%s",cent.Data()));
  
  TH1D *hDimuYPM = (TH1D*) f->Get(Form("hDimuYPM_%s",cent.Data()));
  TH1D *hDimuYPP = (TH1D*) f->Get(Form("hDimuYPP_%s",cent.Data()));
  TH1D *hDimuYMM = (TH1D*) f->Get(Form("hDimuYMM_%s",cent.Data()));
  
  //---------------------------------------------------------
  // Histos
  //---------------------------------------------------------
  
  TString hName = Form("hDimuPM_%s_Mix",cent.Data());
  TH1D *hDimuPM_Mix = new TH1D(hName.Data(),hName.Data(),200,0.,10.);
  hDimuPM_Mix->SetDirectory(0);
  hDimuPM_Mix->Sumw2();
  hName = Form("hDimuPP_%s_Mix",cent.Data());
  TH1D *hDimuPP_Mix = new TH1D(hName.Data(),hName.Data(),200,0.,10.);
  hDimuPP_Mix->SetDirectory(0);
  hDimuPP_Mix->Sumw2();
  hName = Form("hDimuMM_%s_Mix",cent.Data());
  TH1D *hDimuMM_Mix = new TH1D(hName.Data(),hName.Data(),200,0.,10.);
  hDimuMM_Mix->SetDirectory(0);
  hDimuMM_Mix->Sumw2();
  
  hName = Form("hDimuPtPM_%s_Mix",cent.Data());
  TH1D *hDimuPtPM_Mix = new TH1D(hName.Data(),hName.Data(),100,0.,20.);
  hDimuPtPM_Mix->SetDirectory(0);
  hDimuPtPM_Mix->Sumw2();
  hName = Form("hDimuPtPP_%s_Mix",cent.Data());
  TH1D *hDimuPtPP_Mix = new TH1D(hName.Data(),hName.Data(),100,0.,20.);
  hDimuPtPP_Mix->SetDirectory(0);
  hDimuPtPP_Mix->Sumw2();
  hName = Form("hDimuPtMM_%s_Mix",cent.Data());
  TH1D *hDimuPtMM_Mix = new TH1D(hName.Data(),hName.Data(),100,0.,20.);
  hDimuPtMM_Mix->SetDirectory(0);
  hDimuPtMM_Mix->Sumw2();
  
  hName = Form("hDimuYPM_%s_Mix",cent.Data());
  TH1D *hDimuYPM_Mix = new TH1D(hName.Data(),hName.Data(),100,-5.,0.);
  hDimuYPM_Mix->SetDirectory(0);
  hDimuYPM_Mix->Sumw2();
  hName = Form("hDimuYPP_%s_Mix",cent.Data());
  TH1D *hDimuYPP_Mix = new TH1D(hName.Data(),hName.Data(),100,-5.,0.);
  hDimuYPP_Mix->SetDirectory(0);
  hDimuYPP_Mix->Sumw2();
  hName = Form("hDimuYMM_%s_Mix",cent.Data());
  TH1D *hDimuYMM_Mix = new TH1D(hName.Data(),hName.Data(),100,-5.,0.);
  hDimuYMM_Mix->SetDirectory(0);
  hDimuYMM_Mix->Sumw2();
  
  //---------------------------------------------------------
  // Output file
  //---------------------------------------------------------
  
  TString fnameout = Form("%s/MixEvCentr_%s_%d.root", dir.Data(), cent.Data(), Run);
  TFile *fout = new TFile(fnameout.Data(),"RECREATE");
  
  //---------------------------------------------------------
  // Output Tree
  //---------------------------------------------------------
  
//  TTree *treemix = new TTree("TreePb","TreePb");
  Double_t DimuonMassPM,DimuonMassPP, DimuonMassMM;
  Double_t DimuonYPM, DimuonYPP, DimuonYMM;
  Double_t DimuonPtPM, DimuonPtPP, DimuonPtMM;
/*  
  treemix->Branch("DimuonMassPM",&DimuonMassPM,"DimuonMassPM/D");
  treemix->Branch("DimuonRapPM",&DimuonYPM,"DimuonRapPM/D");
  treemix->Branch("DimuonPtPM",&DimuonPtPM,"DimuonPtPM/D");
  treemix->Branch("DimuonMassPP",&DimuonMassPP,"DimuonMassPP/D");
  treemix->Branch("DimuonRapPP",&DimuonYPP,"DimuonRapPP/D");
  treemix->Branch("DimuonPtPP",&DimuonPtPP,"DimuonPtPP/D");
  treemix->Branch("DimuonMassMM",&DimuonMassMM,"DimuonMassMM/D");
  treemix->Branch("DimuonRapMM",&DimuonYMM,"DimuonRapMM/D");
  treemix->Branch("DimuonPtMM",&DimuonPtMM,"DimuonPtMM/D");
*/  
  //---------------------------------------------------------
  // Ev mixing
  //---------------------------------------------------------
  
  Int_t n_muplus = muplus->GetEntries();
  Int_t n_muminus = muminus->GetEntries();
  
  printf("Max number of PM pairs in %s = %ld (nmuplus = %d nmuminus = %d)\n",cent.Data(),((Long_t)n_muplus)*((Long_t)n_muminus),n_muplus,n_muminus);
  printf("N.+- data = %d\n", (Int_t)hDimuPM->GetEntries());
  Int_t nmaxpairs = 100 * hDimuPM->GetEntries();
  printf("Ncreated pairs (PM) = %d\n", nmaxpairs);
  
  for (Int_t i = 0; i < nmaxpairs; i++) {
    
    if ((i+1)%10000 == 0) cout << Form("\rCreating pair... %d / %d (%3.0f%%)", i+1, nmaxpairs, 100.*(i+1.)/nmaxpairs) << flush;
    
    Int_t nmu1 = gRandom->Integer(n_muplus);
    Int_t nmu2 = gRandom->Integer(n_muminus);
    TParticle *mu1 = (TParticle*) muplus->UncheckedAt(nmu1);
    TParticle *mu2 = (TParticle*) muminus->UncheckedAt(nmu2);
    DimuonMassPM = DimuMass(mu1->Energy(),mu1->Px(),mu1->Py(),mu1->Pz(),mu2->Energy(),mu2->Px(),mu2->Py(),mu2->Pz());
    DimuonPtPM = TMath::Sqrt((mu1->Px()+mu2->Px())*(mu1->Px()+mu2->Px())+(mu1->Py()+mu2->Py())*(mu1->Py()+mu2->Py())); 
    DimuonYPM = DimuRap((mu1->Energy()+mu2->Energy()),(mu1->Pz()+mu2->Pz()));
    hDimuPM_Mix->Fill(DimuonMassPM);
    hDimuPtPM_Mix->Fill(DimuonPtPM);
    hDimuYPM_Mix->Fill(DimuonYPM);
    
    nmu1 = gRandom->Integer(n_muplus);
    nmu2 = gRandom->Integer(n_muplus);
    mu1 = (TParticle*) muplus->UncheckedAt(nmu1);
    mu2 = (TParticle*) muplus->UncheckedAt(nmu2);
    DimuonMassPP = DimuMass(mu1->Energy(),mu1->Px(),mu1->Py(),mu1->Pz(),mu2->Energy(),mu2->Px(),mu2->Py(),mu2->Pz());
    DimuonPtPP = TMath::Sqrt((mu1->Px()+mu2->Px())*(mu1->Px()+mu2->Px())+(mu1->Py()+mu2->Py())*(mu1->Py()+mu2->Py())); 
    DimuonYPP = DimuRap((mu1->Energy()+mu2->Energy()),(mu1->Pz()+mu2->Pz()));
    hDimuPP_Mix->Fill(DimuonMassPP);
    hDimuPtPP_Mix->Fill(DimuonPtPP);
    hDimuYPP_Mix->Fill(DimuonYPP);
    
    nmu1 = gRandom->Integer(n_muminus);
    nmu2 = gRandom->Integer(n_muminus);
    mu1 = (TParticle*) muminus->UncheckedAt(nmu1);
    mu2 = (TParticle*) muminus->UncheckedAt(nmu2);
    DimuonMassMM = DimuMass(mu1->Energy(),mu1->Px(),mu1->Py(),mu1->Pz(),mu2->Energy(),mu2->Px(),mu2->Py(),mu2->Pz());
    DimuonPtMM = TMath::Sqrt((mu1->Px()+mu2->Px())*(mu1->Px()+mu2->Px())+(mu1->Py()+mu2->Py())*(mu1->Py()+mu2->Py())); 
    DimuonYMM = DimuRap((mu1->Energy()+mu2->Energy()),(mu1->Pz()+mu2->Pz()));
    hDimuMM_Mix->Fill(DimuonMassMM);
    hDimuPtMM_Mix->Fill(DimuonPtMM);
    hDimuYMM_Mix->Fill(DimuonYMM);
    
//    treemix->Fill();
    
  }
  cout << Form("\rCreating pair... %d / %d (100%%)", nmaxpairs, nmaxpairs) << endl;
  
  //---------------------------------------------------------
  // Scale histos and make ratios
  //---------------------------------------------------------
  
  Double_t normPM=(Double_t)(hDimuPM->Integral(1,200))/((Double_t)hDimuPM_Mix->Integral(1,200));
  printf("NormPM = %f\n",normPM);
  hDimuPM_Mix->Scale(normPM);
  hDimuPtPM_Mix->Scale(normPM);
  hDimuYPM_Mix->Scale(normPM);
  
  TH1D *ratioPM = new TH1D("ratioPM","ratioPM",200,0.,10.);
  ratioPM->SetDirectory(0);
  ratioPM->Sumw2();
  if(hDimuPM->GetEntries()!=0)ratioPM->Divide(hDimuPM_Mix,hDimuPM);
  TH1D *ratioPtPM = new TH1D("ratioPtPM","ratioPtPM",100,0.,20.);
  ratioPtPM->SetDirectory(0);
  ratioPtPM->Sumw2();
  if(hDimuPtPM->GetEntries()!=0)ratioPtPM->Divide(hDimuPtPM_Mix,hDimuPtPM);
  TH1D *ratioYPM = new TH1D("ratioYPM","ratioYPM",100,-5.,0.);
  ratioYPM->SetDirectory(0);
  ratioYPM->Sumw2();
  if(hDimuYPM->GetEntries()!=0)ratioYPM->Divide(hDimuYPM_Mix,hDimuYPM);
  
  Double_t normPP=(Double_t)(hDimuPP->Integral(1,200))/((Double_t)hDimuPP_Mix->Integral(1,200));
  printf("NormPP = %f\n",normPP);
  hDimuPP_Mix->Scale(normPP);
  hDimuPtPP_Mix->Scale(normPP);
  hDimuYPP_Mix->Scale(normPP);
  
  TH1D *ratioPP = new TH1D("ratioPP","ratioPP",200,0.,10.);
  ratioPP->SetDirectory(0);
  ratioPP->Sumw2();
  if(hDimuPP->GetEntries()!=0)ratioPP->Divide(hDimuPP_Mix,hDimuPP);
  TH1D *ratioPtPP = new TH1D("ratioPtPP","ratioPtPP",100,0.,20.);
  ratioPtPP->SetDirectory(0);
  ratioPtPP->Sumw2();
  if(hDimuYPP->GetEntries()!=0)ratioPtPP->Divide(hDimuPtPP_Mix,hDimuPtPP);
  TH1D *ratioYPP = new TH1D("ratioYPP","ratioYPP",100,-5.,0.);
  ratioYPP->SetDirectory(0);
  ratioYPP->Sumw2();
  if(hDimuYPP->GetEntries()!=0)ratioYPP->Divide(hDimuYPP_Mix,hDimuYPP);
  
  Double_t normMM=(Double_t)(hDimuMM->Integral(1,200))/((Double_t)hDimuMM_Mix->Integral(1,200));
  printf("NormMM = %f\n",normMM);
  hDimuMM_Mix->Scale(normMM);
  hDimuPtMM_Mix->Scale(normMM);
  hDimuYMM_Mix->Scale(normMM);
  
  TH1D *ratioMM = new TH1D("ratioMM","ratioMM",200,0.,10.);
  ratioMM->SetDirectory(0);
  ratioMM->Sumw2();
  if(hDimuMM->GetEntries()!=0)ratioMM->Divide(hDimuMM_Mix,hDimuMM);
  TH1D *ratioPtMM = new TH1D("ratioPtMM","ratioPtMM",100,0.,20.);
  ratioPtMM->SetDirectory(0);
  ratioPtMM->Sumw2();
  if(hDimuPtMM->GetEntries()!=0)ratioPtMM->Divide(hDimuPtMM_Mix,hDimuPtMM);
  TH1D *ratioYMM = new TH1D("ratioYMM","ratioYMM",100,-5.,0.);
  ratioYMM->SetDirectory(0);
  ratioYMM->Sumw2();
  if(hDimuYMM->GetEntries()!=0)ratioYMM->Divide(hDimuYMM_Mix,hDimuYMM);
  
  //---------------------------------------------------------
  // Output file
  //---------------------------------------------------------
  
//  treemix->Write();
  hDimuPM->Write();
  hDimuPP->Write();
  hDimuMM->Write();
  hDimuPtPM->Write();
  hDimuPtPP->Write();
  hDimuPtMM->Write();
  hDimuYPM->Write(); 
  hDimuYPP->Write(); 
  hDimuYMM->Write(); 
  hDimuPM_Mix->Write();
  hDimuPP_Mix->Write();
  hDimuMM_Mix->Write();
  hDimuPtPM_Mix->Write();
  hDimuPtPP_Mix->Write();
  hDimuPtMM_Mix->Write();
  hDimuYPM_Mix->Write(); 
  hDimuYPP_Mix->Write(); 
  hDimuYMM_Mix->Write(); 
  ratioPM->Write();
  ratioPtPM->Write();
  ratioYPM->Write();
  ratioPP->Write();
  ratioPtPP->Write();
  ratioYPP->Write();
  ratioMM->Write();
  ratioPtMM->Write();
  ratioYMM->Write();
  fout->Close();
  printf("Writing output file %s \n",fnameout.Data());
  
  //---------------------------------------------------------
  // Draw
  //---------------------------------------------------------
  
  if (display) {
    
    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    TString name = Form("PM_%s_%d.root", cent.Data(), Run);
    TCanvas *cPM = new TCanvas(name.Data(),name.Data(),20,20,1200,600);
    cPM->Divide(3,2);
    cPM->cd(1);
    hDimuPM->SetLineColor(4);
    hDimuPM->DrawClone();
    hDimuPM_Mix->SetLineColor(2);
    hDimuPM_Mix->DrawClone("same");
    cPM->cd(4);
    ratioPM->DrawClone();
    TLine l(0.,1,10.,1.);
    l.DrawClone();
    cPM->cd(2);
    hDimuPtPM->SetLineColor(4);
    hDimuPtPM->DrawClone();
    hDimuPtPM_Mix->SetLineColor(2);
    hDimuPtPM_Mix->DrawClone("same");
    cPM->cd(5);
    ratioPtPM->DrawClone();
    TLine lPt(0.,1,20.,1.);
    lPt.DrawClone();
    cPM->cd(3);
    hDimuYPM->SetLineColor(4);
    hDimuYPM->DrawClone();
    hDimuYPM_Mix->SetLineColor(2);
    hDimuYPM_Mix->DrawClone("same");
    cPM->cd(6);
    ratioYPM->DrawClone();
    TLine lY(-5.,1,0.,1.);
    lY.DrawClone();
    
    name = Form("PP_%s_%d.root", cent.Data(), Run);
    TCanvas *cPP = new TCanvas(name.Data(),name.Data(),40,40,1200,600);
    cPP->Divide(3,2);
    cPP->cd(1);
    hDimuPP->SetLineColor(4);
    hDimuPP->DrawClone();
    hDimuPP_Mix->SetLineColor(2);
    hDimuPP_Mix->DrawClone("same");
    cPP->cd(4);
    ratioPP->DrawClone();
    l.DrawClone();
    cPP->cd(2);
    hDimuPtPP->SetLineColor(4);
    hDimuPtPP->DrawClone();
    hDimuPtPP_Mix->SetLineColor(2);
    hDimuPtPP_Mix->DrawClone("same");
    cPP->cd(5);
    ratioPtPP->DrawClone();
    lPt.DrawClone();
    cPP->cd(3);
    hDimuYPP->SetLineColor(4);
    hDimuYPP->DrawClone();
    hDimuYPP_Mix->SetLineColor(2);
    hDimuYPP_Mix->DrawClone("same");
    cPP->cd(6);
    ratioYPP->DrawClone();
    lY.DrawClone();
    
    name = Form("MM_%s_%d.root", cent.Data(), Run);
    TCanvas *cMM = new TCanvas(name.Data(),name.Data(),60,60,1200,600);
    cMM->Divide(3,2);
    cMM->cd(1);
    hDimuMM->SetLineColor(4);
    hDimuMM->DrawClone();
    hDimuMM_Mix->SetLineColor(2);
    hDimuMM_Mix->DrawClone("same");
    cMM->cd(4);
    ratioMM->DrawClone();
    l.DrawClone();
    cMM->cd(2);
    hDimuPtMM->SetLineColor(4);
    hDimuPtMM->DrawClone();
    hDimuPtMM_Mix->SetLineColor(2);
    hDimuPtMM_Mix->DrawClone("same");
    cMM->cd(5);
    ratioPtMM->DrawClone();
    lPt.DrawClone();
    cMM->cd(3);
    hDimuYMM->SetLineColor(4);
    hDimuYMM->DrawClone();
    hDimuYMM_Mix->SetLineColor(2);
    hDimuYMM_Mix->DrawClone("same");
    cMM->cd(6);
    ratioYMM->DrawClone();
    lY.DrawClone();
    
  }
  
  //---------------------------------------------------
  // clean memory
  //---------------------------------------------------
  
  delete hDimuPM_Mix;
  delete hDimuPP_Mix;
  delete hDimuMM_Mix;
  delete hDimuPtPM_Mix;
  delete hDimuPtPP_Mix;
  delete hDimuPtMM_Mix;
  delete hDimuYPM_Mix; 
  delete hDimuYPP_Mix; 
  delete hDimuYMM_Mix; 
  delete ratioPM;
  delete ratioPtPM;
  delete ratioYPM;
  delete ratioPP;
  delete ratioPtPP;
  delete ratioYPP;
  delete ratioMM;
  delete ratioPtMM;
  delete ratioYMM;
//  delete treemix;
  f->Close();
  
}

//________________________________________________________________________
Double_t DimuMass(Double_t e1, Double_t px1, Double_t py1, Double_t pz1,
	          Double_t e2, Double_t px2, Double_t py2, Double_t pz2) 
{
  Double_t massrec = TMath::Sqrt((e1+e2)*(e1+e2)-((px1+px2)*(px1+px2)+
						  (py1+py2)*(py1+py2)+(pz1+pz2)*(pz1+pz2)));
  return massrec;
}

//________________________________________________________________________
Double_t DimuRap(Double_t e, Double_t pz) 
{
  Double_t rap = 999;
  if (e > TMath::Abs(pz)) rap = 0.5*TMath::Log((e+pz)/(e-pz));
  return rap;
}

