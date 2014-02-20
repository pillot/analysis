/*
 *  runGenTunerLoop.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 08/03/13.
 *  Copyright 2013 SUBATECH. All rights reserved.
 *
 */


//______________________________________________________________________________
void runGenTunerLoop(TString smode = "local", TString inputFileName = "AliAOD.root", Int_t nStep)
{
  /// run the generator tuner in a loop
  
  if (nStep <= 0) return;
  
  // prepare trending plots
  TGraphErrors *gOldPtParam[6];
  TGraphErrors *gOldPtParamMC[6];
  TGraphErrors *gNewPtParam[6];
  for (Int_t i = 0; i < 6; i++) {
    gOldPtParam[i] = new TGraphErrors(nStep);
    gOldPtParam[i]->SetNameTitle(Form("gOldPtParam%d",i), Form("p%d",i));
    gOldPtParamMC[i] = new TGraphErrors(nStep);
    gOldPtParamMC[i]->SetNameTitle(Form("gOldPtParamMC%d",i), Form("p%d",i));
    gNewPtParam[i] = new TGraphErrors(nStep);
    gNewPtParam[i]->SetNameTitle(Form("gNewPtParam%d",i), Form("p%d",i));
  }
  TGraphErrors *gOldYParam[8];
  TGraphErrors *gOldYParamMC[8];
  TGraphErrors *gNewYParam[8];
  for (Int_t i = 0; i < 8; i++) {
    gOldYParam[i] = new TGraphErrors(nStep);
    gOldYParam[i]->SetNameTitle(Form("gOldYParam%d",i), Form("p%d",i));
    gOldYParamMC[i] = new TGraphErrors(nStep);
    gOldYParamMC[i]->SetNameTitle(Form("gOldYParamMC%d",i), Form("p%d",i));
    gNewYParam[i] = new TGraphErrors(nStep);
    gNewYParam[i]->SetNameTitle(Form("gNewYParam%d",i), Form("p%d",i));
  }
  
  TString resume = "";
  for (Int_t i = 0; i < nStep; i++) {
    
    // resume or not
    TString inFileName = Form("Results_step%d.root",i);
    if (resume != "n") {
      if (!gSystem->AccessPathName(inFileName.Data())) {
	if (resume != "y") {
	  cout<<"Results already exist. Do you want to resume? [y/n] (if not previous results will be deleted) "<<flush;
	  do {resume.Gets(stdin,kTRUE);} while (resume != "y" && resume != "n");
	  if (resume == "n") gSystem->Exec("rm -f Results_step*.root");
	}
      } else resume = "n";
    }
    
    // run the generator tuner
    if (resume != "y") {
      if (i == 0)
	gSystem->Exec(Form("root -b -q $WORK/Macros/MuonEfficiency/GenTuner/runGenTuner.C\\(\\\"%s\\\",\\\"%s\\\",%d\\)",
			   smode.Data(), inputFileName.Data(), i));
      else
	gSystem->Exec(Form("root -b -q runGenTuner.C\\(\\\"%s\\\",\\\"%s\\\",%d,\\\'k\\\'\\)",
			   smode.Data(), inputFileName.Data(), i));
    }
    
    // get the new generator parameters and fill trending plots
    TFile *inFile = TFile::Open(inFileName.Data(),"READ");
    if (inFile && inFile->IsOpen()) {
      TF1 *fOldPtFunc = static_cast<TF1*>(inFile->FindObjectAny("fPtFunc"));
      TF1 *fOldPtFuncMC = static_cast<TF1*>(inFile->FindObjectAny("fPtFuncMC"));
      TF1 *fNewPtFunc = static_cast<TF1*>(inFile->FindObjectAny("fPtFuncNew"));
      for (Int_t j = 0; j < 6; j++) {
	if (fOldPtFunc) {
	  gOldPtParam[j]->SetPoint(i, i, fOldPtFunc->GetParameter(j));
	  //gOldPtParam[j]->SetPointError(i, 0., fOldPtFunc->GetParError(j));
	}
	if (fOldPtFuncMC) {
	  gOldPtParamMC[j]->SetPoint(i, i, fOldPtFuncMC->GetParameter(j));
	  //gOldPtParamMC[j]->SetPointError(i, 0., fOldPtFuncMC->GetParError(j));
	}
	if (fNewPtFunc) {
	  gNewPtParam[j]->SetPoint(i, i+1, fNewPtFunc->GetParameter(j));
	  //gNewPtParam[j]->SetPointError(i, 0., fNewPtFunc->GetParError(j));
	}
      }
      TF1 *fOldYFunc = static_cast<TF1*>(inFile->FindObjectAny("fYFunc"));
      TF1 *fOldYFuncMC = static_cast<TF1*>(inFile->FindObjectAny("fYFuncMC"));
      TF1 *fNewYFunc = static_cast<TF1*>(inFile->FindObjectAny("fYFuncNew"));
      for (Int_t j = 0; j < 8; j++) {
	if (fOldYFunc) {
	  gOldYParam[j]->SetPoint(i, i, fOldYFunc->GetParameter(j));
	  //gOldYParam[j]->SetPointError(i, 0., fOldYFunc->GetParError(j));
	}
	if (fOldYFuncMC) {
	  gOldYParamMC[j]->SetPoint(i, i, fOldYFuncMC->GetParameter(j));
	  //gOldYParamMC[j]->SetPointError(i, 0., fOldYFuncMC->GetParError(j));
	}
	if (fNewYFunc) {
	  gNewYParam[j]->SetPoint(i, i+1, fNewYFunc->GetParameter(j));
	  //gNewYParam[j]->SetPointError(i, 0., fNewYFunc->GetParError(j));
	}
      }
      inFile->Close();
    }
    
  }
  
// display trending plots
  TCanvas *cPtParams = new TCanvas("cPtParams", "cPtParams", 600, 400);
  cPtParams->Divide(3,2);
  for (Int_t i = 0; i < 6; i++) {
    cPtParams->cd(i+1);
    gOldPtParam[i]->SetMarkerStyle(kFullDotMedium);
    gOldPtParam[i]->SetMarkerColor(4);
    gOldPtParam[i]->SetLineColor(4);
    gOldPtParam[i]->Draw("ap");
    gOldPtParam[i]->GetXaxis()->SetLimits(-1., nStep+1.);
    gOldPtParamMC[i]->SetMarkerStyle(kFullDotMedium);
    gOldPtParamMC[i]->SetMarkerColor(3);
    gOldPtParamMC[i]->SetLineColor(3);
    gOldPtParamMC[i]->Draw("p");
    gNewPtParam[i]->SetMarkerStyle(kFullDotMedium);
    gNewPtParam[i]->SetMarkerColor(2);
    gNewPtParam[i]->SetLineColor(2);
    gNewPtParam[i]->Draw("p");
  }
  TCanvas *cYParams = new TCanvas("cYParams", "cYParams", 800, 400);
  cYParams->Divide(4,2);
  for (Int_t i = 0; i < 8; i++) {
    cYParams->cd(i+1);
    gOldYParam[i]->SetMarkerStyle(kFullDotMedium);
    gOldYParam[i]->SetMarkerColor(4);
    gOldYParam[i]->SetLineColor(4);
    gOldYParam[i]->Draw("ap");
    gOldYParam[i]->GetXaxis()->SetLimits(-1., nStep+1.);
    gOldYParamMC[i]->SetMarkerStyle(kFullDotMedium);
    gOldYParamMC[i]->SetMarkerColor(3);
    gOldYParamMC[i]->SetLineColor(3);
    gOldYParamMC[i]->Draw("p");
    gNewYParam[i]->SetMarkerStyle(kFullDotMedium);
    gNewYParam[i]->SetMarkerColor(2);
    gNewYParam[i]->SetLineColor(2);
    gNewYParam[i]->Draw("p");
  }
  
  // print and plot last results and save trending plots
  TString inFileName = Form("Results_step%d.root",nStep-1);
  TFile *inFile = TFile::Open(inFileName.Data(),"UPDATE");
  if (inFile && inFile->IsOpen()) {
    TF1 *fPtFuncMC = static_cast<TF1*>(inFile->FindObjectAny("fPtFuncMC"));
    TF1 *fYFuncMC = static_cast<TF1*>(inFile->FindObjectAny("fYFuncMC"));
    if (fPtFuncMC && fYFuncMC) {
      Double_t *param = fPtFuncMC->GetParameters();
      printf("\npT parameters for single muon generator:\n");
      printf("Double_t p[6] = {%g, %g, %g, %g, %g, %g};\n",
	     param[0], param[1], param[2], param[3], param[4], param[5]);
      param = fYFuncMC->GetParameters();
      printf("\ny parameters for single muon generator:\n");
      printf("Double_t p[8] = {%g, %g, %g, %g, %g, %g, %g, %g};\n\n",
	     param[0], param[1], param[2], param[3], param[4], param[5], param[6], param[7]);
    }
    TCanvas *cRes = static_cast<TCanvas*>(inFile->FindObjectAny("cRes"));
    if (cRes) cRes->DrawClone();
    TCanvas *cRat = static_cast<TCanvas*>(inFile->FindObjectAny("cRat"));
    if (cRat) cRat->DrawClone();
    cPtParams->Write(0x0, TObject::kOverwrite);
    cYParams->Write(0x0, TObject::kOverwrite);
  }
  inFile->Close();
  
}

