/* EvtMixNormalization.C
 *
 * Edited by Antoine Lardeux
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <Riostream.h>

#include <TCanvas.h>
#include <TTree.h>
#include <TChain.h>
#include <TMath.h>
#include <TFile.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TList.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <TDirectory.h>
#endif


void MergeCentrality(Int_t i, TString name, Int_t iBin, TH1F* hMergeDimu[12], TH1F* hDimu=0x0, TH1F* hDimuMix=0x0);
TH1F* Substract(TString name="", TH1F* hDimu=0x0, TH1F* hDimuMix=0x0);
void DrawRatio(TString name="", TCanvas* cRatio=0x0, TString sign="", TH1F* hDimu=0x0, TH1F* hDimuMix=0x0, Bool_t zoom=kTRUE);


//---------------------------------------------------------------------------
void EvtMixNormalization(TString dirOut = "SL",
						 Bool_t zoom = kTRUE,
						 TString dir = "../Task/Outputs")
{
    
  // Open File Mix
  TFile *fMix = new TFile(Form("%s/%s/MixOutput_MergeNoNorme.root",dir.Data(),dirOut.Data()), "read");
  if (!fMix || !fMix->IsOpen()) 
  {
	printf("cannot open file MixOutput_MergeNoNorme.root \n");
	return;
  }
  TList *listMix = static_cast<TList*>(fMix->FindObjectAny("cOut")); 
  if (!listMix) 
  {
	printf("cannot find listMix cOut\n");
	return;
  }
  
  // Init
  const Int_t nCent=9;
  const Int_t nCentMerge=6;
  const Int_t npt=24;
  const Int_t ny=10;
  TString centBinName[nCent] = {"010", "1020", "2030", "3040", "4050", "5060", "6070", "7080", "8090"};
  TString MergeBinName[nCentMerge] = {"090", "020", "2040", "4060", "6090", "4090"};
  TH1F *hMergeCent[12];
  Double_t normPM=0.;
  TH1F* hRacc;
  Int_t yshift = 0;
  
  TList *lMixNorm = new TList();
  TList *lSignal = new TList();
  TList *lRacc = new TList();
  TList *lRatio = new TList();
  
  
  // loop on cases pt & y
  TList* var = new TList();
  var->SetOwner(kTRUE);
  var->Add(new TObjString("pt"));
  var->Add(new TObjString("y"));
  TIter nextvar(var);
  TObjString* vartmp;
  while ( ( vartmp = static_cast<TObjString*>(nextvar()) ) )
  {
	printf("######  Observable : %s\n",vartmp->String().Data());
	if (vartmp->String()=="y") yshift = npt*(nCent+nCentMerge);
	
	// loop on pt and y bins
	for (Int_t iBin=0; iBin<npt; iBin++) {
	  if (vartmp->String()=="y" && iBin==ny) break;
	  
	  // loop on centrality
	  for (Int_t i=0; i<nCent; i++) {
		
		printf("--  Cent : %s  (i==%i)\n",centBinName[i].Data(),i);
		
		TString hname = Form("%s_%i_%s",vartmp->String().Data(), iBin, centBinName[i].Data()) ;

	    // Get histos
		TH1F *hDimuPM = static_cast<TH1F*> (listMix->FindObject(Form("hDimuPM_%s",hname.Data() ) ));
		TH1F *hDimuPM_Mix = static_cast<TH1F*> (listMix->FindObject(Form("hDimuPM_%s_Mix",hname.Data() ) ));
		TH1F *hDimuPP = static_cast<TH1F*> (listMix->FindObject( Form("hDimuPP_%s",hname.Data() ) ));
		TH1F *hDimuPP_Mix = static_cast<TH1F*> (listMix->FindObject(Form("hDimuPP_%s_Mix",hname.Data() ) ));
		TH1F *hDimuMM = static_cast<TH1F*> (listMix->FindObject( Form("hDimuMM_%s",hname.Data() ) ));
		TH1F *hDimuMM_Mix = static_cast<TH1F*> (listMix->FindObject(Form("hDimuMM_%s_Mix",hname.Data() ) ));
		hDimuPM->Sumw2();
		hDimuPM_Mix->Sumw2();		
		hDimuPP->Sumw2();
		hDimuPP_Mix->Sumw2();
		hDimuMM->Sumw2();
		hDimuMM_Mix->Sumw2();
		if (!hDimuPM || !hDimuPM_Mix || !hDimuPP || !hDimuPP_Mix || !hDimuMM || !hDimuMM_Mix) {
		  printf("cannot get histos\n");
		  return;
		}
		
		Int_t nbin = hDimuPM->GetNbinsX();
		Double_t xmin = hDimuPM->GetXaxis()->GetBinLowEdge(1);
		Double_t xmax = hDimuPM->GetXaxis()->GetBinUpEdge(nbin);
		Double_t binWidth = (xmax-xmin)/nbin;
		
		// empty bins -> error = 1
		for (Int_t b=1; b < nbin; b++) {
		  Double_t sigerrPM = hDimuPM->GetBinError(b);
		  Double_t sigerrPM_Mix = hDimuPM_Mix->GetBinError(b);
		  Double_t sigerrPP = hDimuPP->GetBinError(b);
		  Double_t sigerrPP_Mix = hDimuPP_Mix->GetBinError(b);
		  Double_t sigerrMM = hDimuMM->GetBinError(b);
		  Double_t sigerrMM_Mix = hDimuMM_Mix->GetBinError(b);
		  if (sigerrPM == 0.) hDimuPM->SetBinError(b,1.);
		  if (sigerrPM_Mix == 0.) hDimuPM_Mix->SetBinError(b,1.);
		  if (sigerrPP == 0.) hDimuPP->SetBinError(b,1.);
		  if (sigerrPP_Mix == 0.) hDimuPP_Mix->SetBinError(b,1.);
		  if (sigerrMM == 0.) hDimuMM->SetBinError(b,1.);
		  if (sigerrMM_Mix == 0.) hDimuMM_Mix->SetBinError(b,1.);
		}
		
		// Racceptance factor
		if (i==0) {
		  TH1F* hRaccTmp = new TH1F( Form("hRaccTmp_%s",hname.Data()), Form("hRaccTmp_%s",hname.Data()),nbin,xmin,xmax);
		  hRaccTmp->Sumw2();
		  
		  Int_t binlimit = hRaccTmp->FindBin(3.);
		  Double_t npmMixRacc, nppMix, nmmMix, errnpmMixRacc, errnppMix, errnmmMix;
		  for (Int_t j=1; j < binlimit; j++) {  //for (Int_t j=1; j < hDimuPM_Mix->GetNbinsX()+1; j++)
			npmMixRacc = hDimuPM_Mix->GetBinContent(j) / binWidth;
			errnpmMixRacc = hDimuPM_Mix->GetBinError(j) / binWidth;
			nppMix = hDimuPP_Mix->GetBinContent(j) / binWidth;
			errnppMix = hDimuPP_Mix->GetBinError(j) / binWidth;
			nmmMix = hDimuMM_Mix->GetBinContent(j) / binWidth;
			errnmmMix = hDimuMM_Mix->GetBinError(j) / binWidth;
			if (nppMix==0. || nmmMix==0.) continue; 
			Double_t nppmmMix = 2*TMath::Sqrt(nppMix*nmmMix);
			
			Double_t errnppmmMix = 2 * TMath::Sqrt(  ((nppMix/nmmMix) * errnmmMix*errnmmMix)  +  ((nmmMix/nppMix) * errnppMix*errnppMix)  ); 
			Double_t normRacc = npmMixRacc/nppmmMix;
			Double_t errnormRacc = TMath::Sqrt( ((errnppmmMix/nppmmMix)*(errnppmmMix/nppmmMix)) + ((errnpmMixRacc/npmMixRacc)*(errnpmMixRacc/npmMixRacc)) )*normRacc;
			hRaccTmp->SetBinContent(j,normRacc);
			hRaccTmp->SetBinError(j,0.); //errnormRacc
		  }
		  for (Int_t j=binlimit; j < nbin+1; j++) {
			hRaccTmp->SetBinContent(j,1.);
			hRaccTmp->SetBinError(j,0.);
		  }

		  // Draw Racceptance
		  TCanvas* c2 = new TCanvas(Form("cRacc_%s",hname.Data()),Form("Racc_%s",hname.Data()),800,800);
		  gPad->SetTopMargin(0.03);
		  gPad->SetRightMargin(0.03);
		  gPad->SetLeftMargin(0.11);
		  c2->SetTickx(1);
		  c2->SetTicky(1);
		  c2->cd();
		  hRaccTmp->Draw();
		  
		  hRacc = (TH1F*) hRaccTmp->Clone( Form("Racc_%s",hname.Data()));
		  
		  lRacc->Add(hRacc);
		  
		  delete hRaccTmp;
		  delete c2;
		}
		
		
		// histo n++n--
		TH1F* hnppmm = new TH1F( Form("hnppmm_%s",hname.Data()), Form("hnppmm_%s",hname.Data()),nbin,xmin,xmax);
		hnppmm->Sumw2();
		
		// Fill histo n++n--
		for (Int_t b=1; b<nbin+1; b++) {
		  Double_t npp = hDimuPP->GetBinContent(b);
		  Double_t errnpp = hDimuPP->GetBinError(b);
		  Double_t nmm = hDimuMM->GetBinContent(b);
		  Double_t errnmm = hDimuMM->GetBinError(b);
		  Double_t nppmm=0., errnppmm=1.;
		  if (npp==0. || nmm==0.) {hnppmm->SetBinContent(b,nppmm);hnppmm->SetBinError(b,errnppmm);}
		  else {
			nppmm = 2 * TMath::Sqrt(npp*nmm);
			errnppmm = 2 * TMath::Sqrt(  ((npp/nmm) * errnmm*errnmm)  +  ((nmm/npp) * errnpp*errnpp)  ); 
			hnppmm->SetBinContent(b,nppmm);
			hnppmm->SetBinError(b,errnppmm);
		  }
		}
		// Apply R factor on histo n++n--
		hnppmm->Multiply(hRacc);
				
		// Compute norm PM
		Int_t binD=hDimuPM_Mix->FindBin(2.);
		Int_t binU=hDimuPM_Mix->FindBin(8.);
		Double_t npmMix = hDimuPM_Mix->Integral(binD,binU);
		Double_t nppmmMain = hnppmm->Integral(binD,binU);
		normPM = nppmmMain/npmMix;
		delete hnppmm;
		
		
		// apply Norm
		hDimuPM_Mix->Scale(normPM);	
		hDimuPP_Mix->Scale(normPM);
		hDimuMM_Mix->Scale(normPM);
		
		
		// Write histos Raw and Mix normalized
		hDimuPM->SetName(Form("hSignal_%s",hname.Data()));
		lMixNorm->Add(hDimuPM);
		//lMixNorm->Add(hDimuPP);
		//lMixNorm->Add(hDimuMM);
		//lMixNorm->Add(hDimuPM_Mix);
		//lMixNorm->Add(hDimuPP_Mix);
		//lMixNorm->Add(hDimuMM_Mix);
		
		
		// merge centrality bins
		MergeCentrality(i, vartmp->String(), iBin, hMergeCent, hDimuPM, hDimuPM_Mix);
		
		
		// Substract
		TH1F* hSignal = (TH1F*) Substract(Form("hSignal_%s",hname.Data()),hDimuPM,hDimuPM_Mix);
		lSignal->Add(hSignal);
		if (i==nCent-1) {
		  for (Int_t j=0; j<nCentMerge*2; j+=2) {
			TH1F* hdimuPMmerge = (TH1F*) hMergeCent[j]->Clone(Form("hSignal_%s_%i_%s",vartmp->String().Data(), iBin, MergeBinName[j/2].Data()));
			lMixNorm->AddAt(hdimuPMmerge,(iBin*(nCent+nCentMerge))+j/2+yshift);
			//TH1F* hdimuPM_Mix_merge = (TH1F*) hMergeCent[j+1]->Clone(Form("hSignal_%s_%i_%s_Mix",vartmp->String().Data(), iBin, MergeBinName[j/2].Data()));
			//lMixNorm->AddAt(hdimuPM_Mix_merge,(iBin*14)+j+1);
			TH1F* hSignalMerge = Substract(Form("hSignal_%s_%i_%s",vartmp->String().Data(), iBin, MergeBinName[j/2].Data()),hMergeCent[j],hMergeCent[j+1]);
			lSignal->AddAt(hSignalMerge,(iBin*(nCent+nCentMerge))+j/2+yshift);
		  }
		}
	  
		
		// Ratio
		TCanvas *cRatioTmp = new TCanvas(hname, "ratioTmp",1200,800);
		TH1F *hDimuPM_tmp = (TH1F*) hDimuPM->Clone();
		TH1F *hDimuPM_Mix_tmp = (TH1F*) hDimuPM_Mix->Clone();
		DrawRatio(Form("Ratio_PM_%s",hname.Data()),cRatioTmp,"PM",hDimuPM_tmp,hDimuPM_Mix_tmp,zoom); 
		TCanvas* cRatio =  (TCanvas*) cRatioTmp->Clone("Ratio");
		delete cRatioTmp;
		lRatio->Add(cRatio);
		if (i==nCent-1) {
		  for (Int_t j=0; j<nCentMerge*2; j+=2) {
			TCanvas *cRatioMergeTmp = new TCanvas(hname, "ratioTmp",1200,800);
			DrawRatio(Form("Ratio_PM_%s_%i_%s",vartmp->String().Data(), iBin, MergeBinName[j/2].Data()),cRatioMergeTmp,"PM",hMergeCent[j],hMergeCent[j+1],zoom); 
			TCanvas* cRatioMerge =  (TCanvas*) cRatioMergeTmp->Clone("Ratio");
			delete cRatioMergeTmp;
		    lRatio->Add(cRatioMerge);
		  }
		}
		//TCanvas* cRatioPP = DrawRatio(Form("Ratio_PP_%s",hname.Data()),"PP",hDimuPP,hDimuPP_Mix,zoom); 
		//lRatio->(cRatioPP);
		//delete cRatioPP;
		//TCanvas* cRatioMM = DrawRatio(Form("Ratio_MM_%s",hname.Data()),"MM",hDimuMM,hDimuMM_Mix,zoom); 
		//lRatio->(cRatioMM);
		//delete cRatioMM;

	  }
	}
  }
    
  
  TFile *fNorm = new TFile(Form("%s/MixNorm.root",dirOut.Data()),"update");  
  lMixNorm->Write("MixNorm",1);
  fNorm->Close();
  delete fNorm;
  
  TFile *fRacc = new TFile(Form("%s/Racc.root",dirOut.Data()),"update");  
  lRacc->Write("Racc",1);
  fRacc->Close();
  delete fRacc;
  
  TFile *fSignal = new TFile(Form("%s/Signal.root",dirOut.Data()),"update");  
  lSignal->Write("Signal",1);
  fSignal->Close();
  delete fSignal;
   
  TFile *fRatio = new TFile(Form("%s/Ratio_zoom%i.root",dirOut.Data(),zoom),"update");  
  lRatio->Write("Ratio",1);
  fRatio->Close();
  delete fRatio;
  
  
  delete lMixNorm;
  delete lRacc;
  delete lRatio;
  delete lSignal;
   
  delete hRacc;
  delete var;
  
  fMix->Close();
  
}

//---------------------------------------------------------------------------
TH1F* Substract(TString name, TH1F* hDimu, TH1F* hDimuMix) 
{
  TH1F* hSignal = (TH1F*) hDimu->Clone(Form("%s",name.Data()));
  for (Int_t b=1; b < hDimu->GetNbinsX(); b++) {
	Double_t sigerr = hSignal->GetBinError(b);
	if (sigerr == 0.) hSignal->SetBinError(b,1.);
  }
  hSignal->Add(hDimuMix, -1);   
  
  return hSignal;
}

//---------------------------------------------------------------------------
void DrawRatio(TString name, TCanvas* cRatio, TString sign, TH1F* hDimu, TH1F* hDimuMix, Bool_t zoom) 
{
  //zoom
  if(zoom) {
	hDimu->SetAxisRange(2.,4.99);
	hDimuMix->SetAxisRange(2.,4.99);
  }
  
  // Compute ratio
  TH1F *ratio = (TH1F*) hDimu->Clone( Form("Ratio_%s_%s", sign.Data(), name.Data()));
  if(hDimuMix->GetEntries()!=0) ratio->Divide(hDimu,hDimuMix);
  
  // Draw
  Float_t fracOfHeight = 0.35;
  Float_t rightMargin = 0.03;
  //TCanvas *cRatio = new TCanvas(name, "ratio",1200,800);
  cRatio->Divide(1,2,0,0);
  
  cRatio->cd(1)->SetLogy();
  gPad->SetPad(0., fracOfHeight, 0.99, 0.99);
  gPad->SetTopMargin(0.03);
  gPad->SetRightMargin(rightMargin);
  
  hDimu->SetTitle("");
  hDimu->SetStats(0);
  //hDimu->SetMinimum(0.);
  //hDimu->SetMaximum(1.04);
  hDimu->GetXaxis()->SetTitle("Inv Mass, Pt, Y");
  hDimu->GetYaxis()->SetLabelSize(0.051);
  hDimu->GetYaxis()->SetTitleSize(0.051);
  hDimu->GetYaxis()->SetTitleOffset(0.8);
  hDimu->GetYaxis()->SetTitle("counts #");
  hDimu->SetMarkerStyle(20);
  hDimu->SetMarkerSize(0.4);
  hDimu->SetMarkerColor(4);
  hDimu->SetLineColor(4);
  hDimu->DrawClone("E");
  
  hDimuMix->SetMarkerStyle(20);
  hDimuMix->SetMarkerSize(0.4);
  hDimuMix->SetMarkerColor(2);
  hDimuMix->SetLineColor(2);
  hDimuMix->DrawClone("Esame");
  
  TLegend *legend = new TLegend (0.8, 0.8, 0.95, 0.95);
  legend->AddEntry(hDimuMix, " Mix", "lep");
  legend->AddEntry(hDimu, " Real", "lep");
  legend->Draw("same");
  
  cRatio->cd(2);
  gPad->SetPad(0., 0., 0.99, fracOfHeight);
  gPad->SetRightMargin(rightMargin);
  gPad->SetBottomMargin(0.08/fracOfHeight);
  gPad->SetGridy();
  
  ratio->SetName("mixOverReal");
  ratio->SetStats(0);
  ratio->SetTitle("");
  ratio->SetLineStyle(1);
  ratio->SetLineColor(1); 
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerSize(0.4);
  ratio->SetMarkerColor(1);
  ratio->GetXaxis()->SetLabelSize(0.11);
  ratio->GetXaxis()->SetTitleSize(0.12);
  ratio->GetXaxis()->SetTitle("GeV/c   ");
  ratio->GetXaxis()->SetTitleOffset(-0.6);
  ratio->GetXaxis()->SetLabelFont(22);
  ratio->GetXaxis()->SetTitleFont(22);
  
  ratio->GetYaxis()->SetLabelSize(0.09);
  ratio->GetYaxis()->SetTitleSize(0.11);
  ratio->GetYaxis()->SetTitle("Ratio");
  ratio->GetYaxis()->SetTitleOffset(0.37);
  ratio->GetYaxis()->SetLabelFont(22);
  ratio->GetYaxis()->SetLabelFont(22);
  ratio->SetMinimum(0.75);   //0.86
  ratio->SetMaximum(1.25);   //1.14
  ratio->DrawClone("E");
  
  TLegend *legend2 = new TLegend (0.60, 0.9, 0.9, 1.0);
  legend2->AddEntry(ratio, " mix / real ", "lep");
  //legend2->Draw("same");
  
  gStyle->SetFillColor(0);
  cRatio->Update();

}

//---------------------------------------------------------------------------
void MergeCentrality(Int_t i, TString name, Int_t iBin, TH1F* hMergeDimu[12], TH1F* hDimu, TH1F* hDimuMix)
{
  
  if (i==0) { // 0-90 & 0-20
	hMergeDimu[0] = (TH1F*) hDimu->Clone( Form("hSignal_%s_%i_090",name.Data(), iBin));
	hMergeDimu[1] = (TH1F*) hDimuMix->Clone( Form("hSignal_%s_%i_090_Mix",name.Data(), iBin));
	hMergeDimu[2] = (TH1F*) hDimu->Clone( Form("hSignal_%s_%i_020",name.Data(), iBin));
	hMergeDimu[3] = (TH1F*) hDimuMix->Clone( Form("hSignal_%s_%i_020_Mix",name.Data(), iBin));
  }
  else {
	if (i==1) {
	  hMergeDimu[2]->Add(hDimu);
	  hMergeDimu[3]->Add(hDimuMix);
	}
	hMergeDimu[0]->Add(hDimu);
	hMergeDimu[1]->Add(hDimuMix);
  }
  
  if (i==2) { // 20-40
	hMergeDimu[4] = (TH1F*) hDimu->Clone( Form("hSignal_%s_%i_2040",name.Data(), iBin));
	hMergeDimu[5] = (TH1F*) hDimuMix->Clone( Form("hSignal_%s_%i_2040_Mix",name.Data(), iBin));
  }
  else if (i==3) {
	hMergeDimu[4]->Add(hDimu);
	hMergeDimu[5]->Add(hDimuMix);
  }
  
  if (i==4) { // 40-60 & 40-90
	hMergeDimu[6] = (TH1F*) hDimu->Clone( Form("hSignal_%s_%i_4060",name.Data(), iBin));
	hMergeDimu[7] = (TH1F*) hDimuMix->Clone( Form("hSignal_%s_%i_4060_Mix",name.Data(), iBin));
	hMergeDimu[10] = (TH1F*) hDimu->Clone( Form("hSignal_%s_%i_4090",name.Data(), iBin));
	hMergeDimu[11] = (TH1F*) hDimuMix->Clone( Form("hSignal_%s_%i_4090_Mix",name.Data(), iBin));
  }
  else if (i>4) {
	if (i==5) {
	  hMergeDimu[6]->Add(hDimu);
	  hMergeDimu[7]->Add(hDimuMix);
	}
	hMergeDimu[10]->Add(hDimu);
	hMergeDimu[11]->Add(hDimuMix);
  }
  
  if (i==6) { // 60-90
	hMergeDimu[8] = (TH1F*) hDimu->Clone( Form("hSignal_%s_%i_6090",name.Data(), iBin));
	hMergeDimu[9] = (TH1F*) hDimuMix->Clone( Form("hSignal_%s_%i_6090_Mix",name.Data(), iBin));
  }
  else if (i>6) {
	hMergeDimu[8]->Add(hDimu);
	hMergeDimu[9]->Add(hDimuMix);
  }
  
}

