#include "TFile.h"
#include "THnSparse.h"
#include "TH1D.h"
#include "TCanvas.h"

Int_t firstMultRange[3] = {1,8}; // limits included
//Int_t firstMultRange[3] = {9,13}; // limits included

void DrawEff(TH1D* numerator, TH1D* denominator, TString name);

//---------------------------------------------------------------------------------------------
void ComputeEfficiencies(TString fileNameData = "AnalysisResults.root")
{
  // WARNING: Ntrk axis contains underflow. When projecting over Ntrk, if all THnSparse entries are projected,
  // underflow is counted in the number of entries. Otherwise, if some entries of THnSparse are not projected,
  // because of some restricted range on other axis, underflow is not counted in the number of entries.
  // In any case, underflow is not counted in the computation of mean, RMS, ...
  // When projecting over other axis, underflow in Ntrk is alway projected and counted in the number of entries,
  // unless the range on the Ntrk axis is restricted.

  TFile *file = new TFile(fileNameData.Data(), "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file %s \n",fileNameData.Data());
    return;
  }
  THnSparse *hNtrkCorrVsCuts = static_cast<THnSparse*>(file->FindObjectAny("hNtrkCorrVsCuts"));
  if (!hNtrkCorrVsCuts) return;

  // trigger+PS efficiency
  //hNtrkCorrVsCuts->GetAxis(2)->SetRangeUser(244540, 300000);
  //hNtrkCorrVsCuts->GetAxis(7)->SetRangeUser(1, 1); // |Vtx MC| < 10 cm
  TH1D* hNch = hNtrkCorrVsCuts->Projection(1,"e");
  hNch->SetName("hNch");
  hNtrkCorrVsCuts->GetAxis(3)->SetRangeUser(1, 1);
  TH1D* hNchPS = hNtrkCorrVsCuts->Projection(1,"e");
  hNchPS->SetName("hNchPS");
  hNtrkCorrVsCuts->GetAxis(3)->SetRange();
  hNtrkCorrVsCuts->GetAxis(9)->SetRangeUser(1, 1);
  TH1D* hNchINELpos = hNtrkCorrVsCuts->Projection(1,"e");
  hNchINELpos->SetName("hNchINELpos");
  hNtrkCorrVsCuts->GetAxis(3)->SetRangeUser(1, 1);
  TH1D* hNchINELposPS = hNtrkCorrVsCuts->Projection(1,"e");
  hNchINELposPS->SetName("hNchINELposPS");
  hNtrkCorrVsCuts->GetAxis(3)->SetRange();
  hNtrkCorrVsCuts->GetAxis(9)->SetRange();
  printf("\n------ trigger+PS efficiency ------\n");
  printf("N(INEL>=0 + PS) / N(INEL>=0) = %f\n",hNchPS->GetEntries()/hNch->GetEntries());
  printf("N(INEL>0 + PS) / N(INEL>0) = %f\n",hNchINELposPS->GetEntries()/hNchINELpos->GetEntries());
  DrawEff(hNchPS, hNch, "hTrigEffvsNch");

  // vtxQA efficiency
  //hNtrkCorrVsCuts->GetAxis(7)->SetRange(); // |Vtx MC| < 10 cm
  //hNtrkCorrVsCuts->GetAxis(6)->SetRangeUser(1, 1); // |Vtx| < 10 cm
  hNtrkCorrVsCuts->GetAxis(3)->SetRangeUser(1, 1);
  hNtrkCorrVsCuts->GetAxis(5)->SetRangeUser(1, 1);
  TH1D* hNchPSVtxQA = hNtrkCorrVsCuts->Projection(1,"e");
  hNchPSVtxQA->SetName("hNchPSVtxQA");
  hNtrkCorrVsCuts->GetAxis(9)->SetRangeUser(1, 1);
  TH1D* hNchINELposPSVtxQA = hNtrkCorrVsCuts->Projection(1,"e");
  hNchINELposPSVtxQA->SetName("hNchINELposPSVtxQA");
  hNtrkCorrVsCuts->GetAxis(3)->SetRange();
  hNtrkCorrVsCuts->GetAxis(5)->SetRange();
  hNtrkCorrVsCuts->GetAxis(9)->SetRange();
  printf("\n------ vtxQA efficiency ------\n");
  printf("N(INEL>=0 + PS + vtxQA) / N(INEL>=0 + PS) = %f\n",hNchPSVtxQA->GetEntries()/hNchPS->GetEntries());
  printf("N(INEL>0 + PS + vtxQA) / N(INEL>0 + PS) = %f\n",hNchINELposPSVtxQA->GetEntries()/hNchINELposPS->GetEntries());
  DrawEff(hNchPSVtxQA, hNchPS, "hVtxQAEffvsNch");

  // |Zvtx| < 10 cm efficiency
  hNtrkCorrVsCuts->GetAxis(3)->SetRangeUser(1, 1);
  hNtrkCorrVsCuts->GetAxis(5)->SetRangeUser(1, 1);
  hNtrkCorrVsCuts->GetAxis(6)->SetRangeUser(1, 1);
  TH1D* hNchPSVtxQAZvtx = hNtrkCorrVsCuts->Projection(1,"e");
  hNchPSVtxQAZvtx->SetName("hNchPSVtxQAZvtx");
  hNtrkCorrVsCuts->GetAxis(9)->SetRangeUser(1, 1);
  TH1D* hNchINELposPSVtxQAZvtx = hNtrkCorrVsCuts->Projection(1,"e");
  hNchINELposPSVtxQAZvtx->SetName("hNchINELposPSVtxQAZvtx");
  hNtrkCorrVsCuts->GetAxis(3)->SetRange();
  hNtrkCorrVsCuts->GetAxis(5)->SetRange();
  hNtrkCorrVsCuts->GetAxis(6)->SetRange();
  hNtrkCorrVsCuts->GetAxis(9)->SetRange();
  printf("\n------ |Zvtx| < 10 cm efficiency ------\n");
  printf("N(INEL>=0 + PS + vtxQA + |Zvtx|<10cm) / N(INEL>=0 + PS + vtxQA) = %f\n",hNchPSVtxQAZvtx->GetEntries()/hNchPSVtxQA->GetEntries());
  printf("N(INEL>0 + PS + vtxQA + |Zvtx|<10cm) / N(INEL>0 + PS + vtxQA) = %f\n",hNchINELposPSVtxQAZvtx->GetEntries()/hNchINELposPSVtxQA->GetEntries());
  DrawEff(hNchPSVtxQAZvtx, hNchPSVtxQA, "hZvtxEffvsNch");

  // INEL=0 contamination
  printf("\n------ INEL=0 contamination ------\n");
  printf("N(INEL=0 + PS) / N(INEL>=0 + PS) = %f\n",hNchPS->GetBinContent(1)/hNchPS->GetEntries());
  printf("N(INEL=0 + PS + vtxQA) / N(INEL>=0 + PS + vtxQA) = %f\n",hNchPSVtxQA->GetBinContent(1)/hNchPSVtxQA->GetEntries());

  // --------------------------- efficiencies for first multiplicity bin
  printf("\n====== efficiencies for first multiplicity bin ======\n");

  // trigger+PS efficiency
  hNtrkCorrVsCuts->GetAxis(7)->SetRangeUser(1, 1); // |Vtx MC| < 10 cm (NtrkCorr = -999 if |Vtx| > 10 cm)
  //hNtrkCorrVsCuts->GetAxis(5)->SetRangeUser(1, 1);
  TH1D* hNtrk = hNtrkCorrVsCuts->Projection(0,"e");
  hNtrk->SetName("hNtrk");
  hNtrkCorrVsCuts->GetAxis(3)->SetRangeUser(1, 1);
  TH1D* hNtrkPS = hNtrkCorrVsCuts->Projection(0,"e");
  hNtrkPS->SetName("hNtrkPS");
  hNtrkCorrVsCuts->GetAxis(3)->SetRange();
  hNtrkCorrVsCuts->GetAxis(9)->SetRangeUser(1, 1);
  TH1D* hNtrkINELpos = hNtrkCorrVsCuts->Projection(0,"e");
  hNtrkINELpos->SetName("hNtrkINELpos");
  hNtrkCorrVsCuts->GetAxis(3)->SetRangeUser(1, 1);
  TH1D* hNtrkINELposPS = hNtrkCorrVsCuts->Projection(0,"e");
  hNtrkINELposPS->SetName("hNtrkINELposPS");
  hNtrkCorrVsCuts->GetAxis(3)->SetRange();
  hNtrkCorrVsCuts->GetAxis(9)->SetRange();
  hNtrkCorrVsCuts->GetAxis(0)->SetRangeUser(firstMultRange[0],firstMultRange[1]);
  TH1D* hNch1stTrkBin = hNtrkCorrVsCuts->Projection(1,"e");
  hNch1stTrkBin->SetName("hNch1stTrkBin");
  hNtrkCorrVsCuts->GetAxis(3)->SetRangeUser(1, 1);
  TH1D* hNch1stTrkBinPS = hNtrkCorrVsCuts->Projection(1,"e");
  hNch1stTrkBinPS->SetName("hNch1stTrkBinPS");
  hNtrkCorrVsCuts->GetAxis(0)->SetRange();
  hNtrkCorrVsCuts->GetAxis(3)->SetRange();
  printf("\n------ trigger+PS efficiency ------\n");
  printf("N(Nch in [%d,%d] + INEL>=0 + PS) / N(Nch in [%d,%d] + INEL>=0) = %f\n",firstMultRange[0],firstMultRange[1],firstMultRange[0],firstMultRange[1],hNchPS->Integral(firstMultRange[0]+1,firstMultRange[1]+1)/hNch->Integral(firstMultRange[0]+1,firstMultRange[1]+1));
  printf("N(Ntrk in [%d,%d] + INEL>=0 + PS) / N(Ntrk in [%d,%d] + INEL>=0) = %f\n",firstMultRange[0],firstMultRange[1],firstMultRange[0],firstMultRange[1],hNtrkPS->Integral(firstMultRange[0]+1,firstMultRange[1]+1)/hNtrk->Integral(firstMultRange[0]+1,firstMultRange[1]+1));
  printf("N(Ntrk in [%d,%d] + INEL>0 + PS) / N(Ntrk in [%d,%d] + INEL>0) = %f\n",firstMultRange[0],firstMultRange[1],firstMultRange[0],firstMultRange[1],hNtrkINELposPS->Integral(firstMultRange[0]+1,firstMultRange[1]+1)/hNtrkINELpos->Integral(firstMultRange[0]+1,firstMultRange[1]+1));
  printf("------\n");
  hNch->GetXaxis()->SetRangeUser(firstMultRange[0],firstMultRange[1]);
  hNchPS->GetXaxis()->SetRangeUser(firstMultRange[0],firstMultRange[1]);
  printf("<Nch>(Nch in [%d,%d] + INEL>=0 + PS) / <Nch>(Nch in [%d,%d] + INEL>=0) = %f\n",firstMultRange[0],firstMultRange[1],firstMultRange[0],firstMultRange[1],hNchPS->GetMean()/hNch->GetMean());
  hNch->GetXaxis()->SetRange();
  hNchPS->GetXaxis()->SetRange();
  printf("<Nch>(Ntrk in [%d,%d] + INEL>=0 + PS) / <Nch>(Ntrk in [%d,%d] + INEL>=0) = %f\n",firstMultRange[0],firstMultRange[1],firstMultRange[0],firstMultRange[1],hNch1stTrkBinPS->GetMean()/hNch1stTrkBin->GetMean());
  hNch1stTrkBin->GetXaxis()->SetRangeUser(1,hNch1stTrkBin->GetXaxis()->GetXmax());
  hNch1stTrkBinPS->GetXaxis()->SetRangeUser(1,hNch1stTrkBinPS->GetXaxis()->GetXmax());
  printf("<Nch>(Ntrk in [%d,%d] + INEL>0 + PS) / <Nch>(Ntrk in [%d,%d] + INEL>0) = %f\n",firstMultRange[0],firstMultRange[1],firstMultRange[0],firstMultRange[1],hNch1stTrkBinPS->GetMean()/hNch1stTrkBin->GetMean());
  hNch1stTrkBin->GetXaxis()->SetRange();
  hNch1stTrkBinPS->GetXaxis()->SetRange();
  hNtrk->GetXaxis()->SetRangeUser(firstMultRange[0],firstMultRange[1]);
  hNtrkPS->GetXaxis()->SetRangeUser(firstMultRange[0],firstMultRange[1]);
  printf("<Ntrk>(Ntrk in [%d,%d] + INEL>=0 + PS) / <Ntrk>(Ntrk in [%d,%d] + INEL>=0) = %f\n",firstMultRange[0],firstMultRange[1],firstMultRange[0],firstMultRange[1],hNtrkPS->GetMean()/hNtrk->GetMean());
  hNtrk->GetXaxis()->SetRange();
  hNtrkPS->GetXaxis()->SetRange();
  hNtrkINELpos->GetXaxis()->SetRangeUser(firstMultRange[0],firstMultRange[1]);
  hNtrkINELposPS->GetXaxis()->SetRangeUser(firstMultRange[0],firstMultRange[1]);
  printf("<Ntrk>(Ntrk in [%d,%d] + INEL>0 + PS) / <Ntrk>(Ntrk in [%d,%d] + INEL>0) = %f\n",firstMultRange[0],firstMultRange[1],firstMultRange[0],firstMultRange[1],hNtrkINELposPS->GetMean()/hNtrkINELpos->GetMean());
  hNtrkINELpos->GetXaxis()->SetRange();
  hNtrkINELposPS->GetXaxis()->SetRange();
  DrawEff(hNtrkPS, hNtrk, "hTrigEffvsNtrk");
  DrawEff(hNtrkINELposPS, hNtrkINELpos, "hTrigEffINELposvsNtrk");

  // vtxQA efficiency
  hNtrkCorrVsCuts->GetAxis(3)->SetRangeUser(1, 1);
  hNtrkCorrVsCuts->GetAxis(5)->SetRangeUser(1, 1);
  TH1D* hNtrkPSVtxQA = hNtrkCorrVsCuts->Projection(0,"e");
  hNtrkPSVtxQA->SetName("hNtrkPSVtxQA");
  hNtrkCorrVsCuts->GetAxis(9)->SetRangeUser(1, 1);
  TH1D* hNtrkINELposPSVtxQA = hNtrkCorrVsCuts->Projection(0,"e");
  hNtrkINELposPSVtxQA->SetName("hNtrkINELposPSVtxQA");
  hNtrkCorrVsCuts->GetAxis(9)->SetRange();
  hNtrkCorrVsCuts->GetAxis(0)->SetRangeUser(firstMultRange[0],firstMultRange[1]);
  TH1D* hNch1stTrkBinPSVtxQA = hNtrkCorrVsCuts->Projection(1,"e");
  hNch1stTrkBinPSVtxQA->SetName("hNch1stTrkBinPSVtxQA");
  hNtrkCorrVsCuts->GetAxis(0)->SetRange();
  hNtrkCorrVsCuts->GetAxis(3)->SetRange();
  hNtrkCorrVsCuts->GetAxis(5)->SetRange();
  printf("\n------ vtxQA efficiency ------\n");
  printf("N(Nch in [%d,%d] + INEL>=0 + PS + vtxQA) / N(Nch in [%d,%d] + INEL>=0 + PS) = %f\n",firstMultRange[0],firstMultRange[1],firstMultRange[0],firstMultRange[1],hNchPSVtxQA->Integral(firstMultRange[0]+1,firstMultRange[1]+1)/hNchPS->Integral(firstMultRange[0]+1,firstMultRange[1]+1));
  printf("N(Ntrk in [%d,%d] + INEL>=0 + PS + vtxQA) / N(Ntrk in [%d,%d] + INEL>=0 + PS) = %f\n",firstMultRange[0],firstMultRange[1],firstMultRange[0],firstMultRange[1],hNtrkPSVtxQA->Integral(firstMultRange[0]+1,firstMultRange[1]+1)/hNtrkPS->Integral(firstMultRange[0]+1,firstMultRange[1]+1));
  printf("N(Ntrk in [%d,%d] + INEL>0 + PS + vtxQA) / N(Ntrk in [%d,%d] + INEL>0 + PS) = %f\n",firstMultRange[0],firstMultRange[1],firstMultRange[0],firstMultRange[1],hNtrkINELposPSVtxQA->Integral(firstMultRange[0]+1,firstMultRange[1]+1)/hNtrkINELposPS->Integral(firstMultRange[0]+1,firstMultRange[1]+1));
  printf("------\n");
  hNchPS->GetXaxis()->SetRangeUser(firstMultRange[0],firstMultRange[1]);
  hNchPSVtxQA->GetXaxis()->SetRangeUser(firstMultRange[0],firstMultRange[1]);
  printf("<Nch>(Nch in [%d,%d] + INEL>=0 + PS + vtxQA) / <Nch>(Nch in [%d,%d] + INEL>=0 + PS) = %f\n",firstMultRange[0],firstMultRange[1],firstMultRange[0],firstMultRange[1],hNchPSVtxQA->GetMean()/hNchPS->GetMean());
  hNchPS->GetXaxis()->SetRange();
  hNchPSVtxQA->GetXaxis()->SetRange();
  printf("<Nch>(Ntrk in [%d,%d] + INEL>=0 + PS + vtxQA) / <Nch>(Ntrk in [%d,%d] + INEL>=0 + PS) = %f\n",firstMultRange[0],firstMultRange[1],firstMultRange[0],firstMultRange[1],hNch1stTrkBinPSVtxQA->GetMean()/hNch1stTrkBinPS->GetMean());
  hNch1stTrkBinPS->GetXaxis()->SetRangeUser(1,hNch1stTrkBinPS->GetXaxis()->GetXmax());
  hNch1stTrkBinPSVtxQA->GetXaxis()->SetRangeUser(1,hNch1stTrkBinPSVtxQA->GetXaxis()->GetXmax());
  printf("<Nch>(Ntrk in [%d,%d] + INEL>0 + PS + vtxQA) / <Nch>(Ntrk in [%d,%d] + INEL>0 + PS) = %f\n",firstMultRange[0],firstMultRange[1],firstMultRange[0],firstMultRange[1],hNch1stTrkBinPSVtxQA->GetMean()/hNch1stTrkBinPS->GetMean());
  hNch1stTrkBinPS->GetXaxis()->SetRange();
  hNch1stTrkBinPSVtxQA->GetXaxis()->SetRange();
  hNtrkPS->GetXaxis()->SetRangeUser(firstMultRange[0],firstMultRange[1]);
  hNtrkPSVtxQA->GetXaxis()->SetRangeUser(firstMultRange[0],firstMultRange[1]);
  printf("<Ntrk>(Ntrk in [%d,%d] + INEL>=0 + PS + vtxQA) / <Ntrk>(Ntrk in [%d,%d] + INEL>=0 + PS) = %f\n",firstMultRange[0],firstMultRange[1],firstMultRange[0],firstMultRange[1],hNtrkPSVtxQA->GetMean()/hNtrkPS->GetMean());
  hNtrkPS->GetXaxis()->SetRange();
  hNtrkPSVtxQA->GetXaxis()->SetRange();
  hNtrkINELposPS->GetXaxis()->SetRangeUser(firstMultRange[0],firstMultRange[1]);
  hNtrkINELposPSVtxQA->GetXaxis()->SetRangeUser(firstMultRange[0],firstMultRange[1]);
  printf("<Ntrk>(Ntrk in [%d,%d] + INEL>0 + PS + vtxQA) / <Ntrk>(Ntrk in [%d,%d] + INEL>0 + PS) = %f\n",firstMultRange[0],firstMultRange[1],firstMultRange[0],firstMultRange[1],hNtrkINELposPSVtxQA->GetMean()/hNtrkINELposPS->GetMean());
  hNtrkINELposPS->GetXaxis()->SetRange();
  hNtrkINELposPSVtxQA->GetXaxis()->SetRange();
  DrawEff(hNtrkPSVtxQA, hNtrkPS, "hVtxQAEffvsNtrk");

  // INEL=0 contamination
  printf("\n------ INEL=0 contamination ------\n");
  printf("N(Ntrk in [%d,%d] + INEL=0 + PS) / N(Ntrk in [%d,%d] + INEL>=0 + PS) = %f\n",firstMultRange[0],firstMultRange[1],firstMultRange[0],firstMultRange[1],hNch1stTrkBinPS->GetBinContent(1)/hNch1stTrkBinPS->GetEntries());
  printf("N(Ntrk in [%d,%d] + INEL=0 + PS + vtxQA) / N(Ntrk in [%d,%d] + INEL>=0 + PS + vtxQA) = %f\n",firstMultRange[0],firstMultRange[1],firstMultRange[0],firstMultRange[1],hNch1stTrkBinPSVtxQA->GetBinContent(1)/hNch1stTrkBinPSVtxQA->GetEntries());
}

//---------------------------------------------------------------------------------------------
void DrawEff(TH1D* numerator, TH1D* denominator, TString name)
{
  TH1D* hEff = static_cast<TH1D*>(numerator->Clone(name.Data()));
  hEff->Divide(denominator);
  TCanvas* cEff = new TCanvas(name.Data(),name.Data());
  cEff->Divide(1,2);
  cEff->cd(1);
  denominator->SetLineColor(4);
  denominator->Draw();
  numerator->SetLineColor(2);
  numerator->Draw("sames");
  cEff->cd(2);
  hEff->Draw();
}
