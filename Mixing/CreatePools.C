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
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TGraphErrors.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#endif

Double_t DimuMass(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t); 
Double_t DimuRap(Double_t, Double_t);

void CreatePools(TString indir, TString outdir, Int_t Run){
  
  Double_t muMass = TDatabasePDG::Instance()->GetParticle("mu-")->Mass(); // GeV
  
  //---------------------------------------------------------
  // Load files
  //---------------------------------------------------------
  
  TChain *c = new TChain("fTreeSingleMu");
  c->Add(Form("%s/CS_SingleMuon_000%d.root", indir.Data(), Run));
  
  //---------------------------------------------------------
  // Create dimuon histos
  //---------------------------------------------------------
  
  const Int_t nCent = 11;
  Double_t centBinRange[nCent][2] = {{0., 10.}, {10., 20.}, {20., 40.}, {40., 60.}, {60., 80.},
    {20., 30.}, {30., 40.}, {40., 50.}, {50., 60.}, {60., 70.}, {70., 80.}};
  TString centBinName[nCent] = {"010", "1020", "2040", "4060", "6080", "2030", "3040", "4050", "5060", "6070", "7080"};
  TH1D *hDimuPM_cent[nCent];
  TH1D *hDimuPP_cent[nCent];
  TH1D *hDimuMM_cent[nCent];
  TH1D *hDimuPtPM_cent[nCent];
  TH1D *hDimuPtPP_cent[nCent];
  TH1D *hDimuPtMM_cent[nCent];
  TH1D *hDimuYPM_cent[nCent];
  TH1D *hDimuYPP_cent[nCent];
  TH1D *hDimuYMM_cent[nCent];
  TString hName;
  for (Int_t i = 0; i < nCent; i++) {
    
    hName = Form("hDimuPM_%s", centBinName[i].Data());
    hDimuPM_cent[i] = new TH1D(hName.Data(),hName.Data(),200,0.,10.);
    hDimuPM_cent[i]->SetDirectory(0);
    hName = Form("hDimuPP_%s", centBinName[i].Data());
    hDimuPP_cent[i] = new TH1D(hName.Data(),hName.Data(),200,0.,10.);
    hDimuPP_cent[i]->SetDirectory(0);
    hName = Form("hDimuMM_%s", centBinName[i].Data());
    hDimuMM_cent[i] = new TH1D(hName.Data(),hName.Data(),200,0.,10.);
    hDimuMM_cent[i]->SetDirectory(0);
    
    hName = Form("hDimuPtPM_%s", centBinName[i].Data());
    hDimuPtPM_cent[i] = new TH1D(hName.Data(),hName.Data(),100,0.,20.);
    hDimuPtPM_cent[i]->SetDirectory(0);
    hName = Form("hDimuPtPP_%s", centBinName[i].Data());
    hDimuPtPP_cent[i] = new TH1D(hName.Data(),hName.Data(),100,0.,20.);
    hDimuPtPP_cent[i]->SetDirectory(0);
    hName = Form("hDimuPtMM_%s", centBinName[i].Data());
    hDimuPtMM_cent[i] = new TH1D(hName.Data(),hName.Data(),100,0.,20.);
    hDimuPtMM_cent[i]->SetDirectory(0);
    
    hName = Form("hDimuYPM_%s", centBinName[i].Data());
    hDimuYPM_cent[i] = new TH1D(hName.Data(),hName.Data(),100,-5.,0.);
    hDimuYPM_cent[i]->SetDirectory(0);
    hName = Form("hDimuYPP_%s", centBinName[i].Data());
    hDimuYPP_cent[i] = new TH1D(hName.Data(),hName.Data(),100,-5.,0.);
    hDimuYPP_cent[i]->SetDirectory(0);
    hName = Form("hDimuYMM_%s", centBinName[i].Data());
    hDimuYMM_cent[i] = new TH1D(hName.Data(),hName.Data(),100,-5.,0.);
    hDimuYMM_cent[i]->SetDirectory(0);
    
  }
  
  //---------------------------------------------------------
  // Load branches
  //---------------------------------------------------------
  
  UInt_t PassPhysicsSelection;
  c->SetBranchAddress("PassPhysicsSelection",&PassPhysicsSelection);
  
  UInt_t L0TrigInp;
  c->SetBranchAddress("L0TrigInp",&L0TrigInp);
  
  Int_t ZDCAccept;
  c->SetBranchAddress("ZDCAccept",&ZDCAccept);
  
  Float_t V0Cent;
  c->SetBranchAddress("V0Cent",&V0Cent);
  
  Float_t IPVx, IPVy, IPVz;
  c->SetBranchAddress("IPVx",&IPVx);
  c->SetBranchAddress("IPVy",&IPVy);
  c->SetBranchAddress("IPVz",&IPVz);
  
  Int_t NumOfMuonTracks;
  c->SetBranchAddress("NumOfMuonTracks",&NumOfMuonTracks);
  
  Int_t IsMuon;
  c->SetBranchAddress("IsMuon",&IsMuon);
  
  Int_t MatchTrig;
  c->SetBranchAddress("MatchTrig",&MatchTrig);
  
  Float_t Eta;
  c->SetBranchAddress("Eta",&Eta);
  
  Float_t RAbs;
  c->SetBranchAddress("RAbs",&RAbs);
  
  Float_t Charge;
  c->SetBranchAddress("Charge",&Charge);
  
  Float_t Px, Py, Pz;
  c->SetBranchAddress("Px",&Px);
  c->SetBranchAddress("Py",&Py);
  c->SetBranchAddress("Pz",&Pz);
  
  //---------------------------------------------------------
  // Loop on entries 
  //---------------------------------------------------------
  
  Int_t NEntries = c->GetEntries();
  printf("Entries = %d\n",(Int_t)c->GetEntries());
  
  Int_t imuplus[nCent];
  Int_t imuminus[nCent];
  TClonesArray *muplus_cent[nCent];
  TClonesArray *muminus_cent[nCent];
  for (Int_t i = 0; i < nCent; i++) {
    muplus_cent[i] = new TClonesArray("TParticle",0);
    muminus_cent[i] = new TClonesArray("TParticle",0);
    imuplus[i] = 0;
    imuminus[i] = 0;
  }
  
  Int_t nMuon = 0;
  Int_t nTotMuon = 0;
  TObjArray currentMu;
  for (Int_t i = 0; i < NEntries; i++) {
    
    if ((i+1)%10000 == 0) cout << Form("\rEvent processing... %d / %d (%3.0f%%)", i+1, NEntries, 100.*(i+1.)/NEntries) << flush;
    c->GetEntry(i);
    
    // select events
    if (NumOfMuonTracks == 0) continue;
    if (!(PassPhysicsSelection&1)) continue;
    if (!((L0TrigInp&4096) && (L0TrigInp&2) && (L0TrigInp&4))) continue;
    if (!ZDCAccept) continue;
    
    // starting a new event
    if (nTotMuon == 0) nTotMuon = NumOfMuonTracks;
    
    // increment the muon counter
    nMuon++;
    
    // muon energy
    Double_t E = TMath::Sqrt(Px*Px + Py*Py + Pz*Pz + muMass*muMass);
    
    // select muons in acceptance matching the trigger
    if (IsMuon > 0 && MatchTrig > 0 && Eta > -4. && Eta < -2.5 && RAbs > 17.622 && RAbs < 89.) {
      
      // fill pools
      TParticle* muon = 0x0;
      for (Int_t iCent = 0; iCent < nCent; iCent++) {
	
	if (V0Cent > centBinRange[iCent][0] && V0Cent <= centBinRange[iCent][1]) {
	  
	  if (Charge > 0.5) {
	    
	    muon = new((*(muplus_cent[iCent]))[imuplus[iCent]++]) TParticle(-13,0,-1,-1,-1,-1,Px,Py,Pz,E,IPVx,IPVy,IPVz,0.);
	    
	  } else if (Charge < -0.5) {
	    
	    muon = new((*(muminus_cent[iCent]))[imuminus[iCent]++]) TParticle(13,0,-1,-1,-1,-1,Px,Py,Pz,E,IPVx,IPVy,IPVz,0.);
	    
	  }
	  
	}
	
      }
      
      // stock muons of current event
      if (muon) currentMu.AddLast(muon);
      
    }
    
    // continue until we reach the last track of the current event
    if (nMuon < nTotMuon) continue;
    
    // build dimuons and fill histos
    for (Int_t imu1 = 0; imu1 < currentMu.GetEntries(); imu1++) {
      
      TParticle* mu1 = static_cast<TParticle*>(currentMu.UncheckedAt(imu1));
      Int_t pdg1 = mu1->GetPdgCode();
      
      for (Int_t imu2 = imu1+1; imu2 < currentMu.GetEntries(); imu2++) {
	
	TParticle* mu2 = static_cast<TParticle*>(currentMu.UncheckedAt(imu2));
	
	// select dimuons
	Double_t dimuY = DimuRap((mu1->Energy()+mu2->Energy()),(mu1->Pz()+mu2->Pz()));
	if(dimuY < -4. || dimuY > -2.5) continue;
	
	// fill histos
	Int_t pdg2 = mu2->GetPdgCode();
	Double_t dimuMass = DimuMass(mu1->Energy(),mu1->Px(),mu1->Py(),mu1->Pz(),mu2->Energy(),mu2->Px(),mu2->Py(),mu2->Pz());
	Double_t dimuPt = TMath::Sqrt((mu1->Px()+mu2->Px())*(mu1->Px()+mu2->Px())+(mu1->Py()+mu2->Py())*(mu1->Py()+mu2->Py()));
	for (Int_t iCent = 0; iCent < nCent; iCent++) {
	  
	  if (V0Cent > centBinRange[iCent][0] && V0Cent <= centBinRange[iCent][1]) {
	    
	    if (pdg1 < 0 && pdg2 < 0) {
	      
	      hDimuPP_cent[iCent]->Fill(dimuMass);
	      hDimuPtPP_cent[iCent]->Fill(dimuPt);
	      hDimuYPP_cent[iCent]->Fill(dimuY);
	      
	    } else if  (pdg1 > 0 && pdg2 > 0) {
	      
	      hDimuMM_cent[iCent]->Fill(dimuMass);
	      hDimuPtMM_cent[iCent]->Fill(dimuPt);
	      hDimuYMM_cent[iCent]->Fill(dimuY);
	      
	    } else {
	      
	      hDimuPM_cent[iCent]->Fill(dimuMass);
	      hDimuPtPM_cent[iCent]->Fill(dimuPt);
	      hDimuYPM_cent[iCent]->Fill(dimuY);
	      
	    }
	    
	  }
	  
	}
	
      }
      
    }
    
    // prepare for the next event
    currentMu.Clear();
    nTotMuon = 0;
    nMuon = 0;
    
  }
  cout << Form("\rEvent processing... %d / %d (100%%)", NEntries, NEntries) << endl;
  
  //---------------------------------------------------
  // Save histos
  //---------------------------------------------------
  
  TString fnameout = Form("%s/PoolCentr_%d.root", outdir.Data(), Run);
  TFile *fout = new TFile(fnameout.Data(), "RECREATE");
  TString name;
  for (Int_t i = 0; i < nCent; i++) {
    name = Form("muplus_%s", centBinName[i].Data());
    muplus_cent[i]->Write(name.Data(),1);
    name = Form("muminus_%s", centBinName[i].Data());
    muminus_cent[i]->Write(name.Data(),1);
    hDimuPM_cent[i]->Write();
    hDimuPP_cent[i]->Write();
    hDimuMM_cent[i]->Write();
    hDimuPtPM_cent[i]->Write();
    hDimuPtPP_cent[i]->Write();
    hDimuPtMM_cent[i]->Write();
    hDimuYPM_cent[i]->Write();
    hDimuYPP_cent[i]->Write();
    hDimuYMM_cent[i]->Write();
  }
  fout->Close();
  printf("Writing file %s\n",fnameout.Data());
  
  //---------------------------------------------------
  // clean memory
  //---------------------------------------------------
  
  for (Int_t i = 0; i < nCent; i++) {
    delete hDimuPM_cent[i];
    delete hDimuPP_cent[i];
    delete hDimuMM_cent[i];
    delete hDimuPtPM_cent[i];
    delete hDimuPtPP_cent[i];
    delete hDimuPtMM_cent[i];
    delete hDimuYPM_cent[i];
    delete hDimuYPP_cent[i];
    delete hDimuYMM_cent[i];
    delete muplus_cent[i];
    delete muminus_cent[i];
  }
  delete c;
  
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

