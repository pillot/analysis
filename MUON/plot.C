#include <riostream>
#include <TMath.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TRandom.h>

/*
  FinalizeTrack(*track);

  AliMUONTrackParam* lastTrackParam = (AliMUONTrackParam*) track->GetTrackParamAtCluster()->Last();
  AliMUONVCluster* cluster2 = lastTrackParam->GetClusterPtr();
  Double_t x2 = cluster2->GetX();
  Double_t y2 = cluster2->GetY();
  Double_t z2 = cluster2->GetZ();
  AliMUONTrackParam* previousTrackParam = (AliMUONTrackParam*) track->GetTrackParamAtCluster()->Before(lastTrackParam);
  AliMUONVCluster* cluster1 = previousTrackParam->GetClusterPtr();
  Double_t x1 = cluster1->GetX();
  Double_t y1 = cluster1->GetY();
  Double_t z1 = cluster1->GetZ();
  cout<<"estimated non bending slope = "<<(x1 - x2) / (z1 - z2)<<endl;
  cout<<"estimated non bending impact param = "<<x2 - z2 * (x1 - x2) / (z1 - z2)<<endl;
  Double_t bendingImpact = y2 - z2 * (y1 - y2) / (z1 - z2);
  Double_t bendingMomentum = AliMUONTrackExtrap::GetBendingMomentumFromImpactParam(bendingImpact);
  cout<<"estimated bending momentum = "<<bendingMomentum<<endl;
  AliMUONTrackParam* firstTrackParam = (AliMUONTrackParam*) track->GetTrackParamAtCluster()->First();
  Double_t recoBendingMomentum = 1. / firstTrackParam->GetInverseBendingMomentum();
  cout<<"reconstructed bending momentum = "<<recoBendingMomentum<<endl;
  const TMatrixD& kParamCov = firstTrackParam->GetCovariances();
  Double_t errRecoBendingMomentum = TMath::Sqrt(kParamCov(4,4)) * recoBendingMomentum * recoBendingMomentum;
  cout<<"reconstructed bending momentum error = "<<errRecoBendingMomentum<<endl;
  Double_t errBendingMomentum = - TMath::Sqrt(z1 * z1 * cluster2->GetErrY2() + z2 * z2 * cluster1->GetErrY2()) / (z1 - z2) / bendingImpact * bendingMomentum;
  cout<<"estimated bending momentum error = "<<errBendingMomentum<<endl;
  Double_t errVertexBendingMomentum = - bendingMomentum / bendingImpact * AliMUONReconstructor::GetRecoParam()->GetBendingVertexDispersion();
  cout<<"estimated bending momentum vertex error = "<<errVertexBendingMomentum<<endl;
*/

void plot() {
  
  ifstream nonBendSlopeFile;
  nonBendSlopeFile.open("NonBendingSlope.txt", ifstream::in);
  ifstream nonBendImpactFile;
  nonBendImpactFile.open("NonBendingImpactParam.txt", ifstream::in);
  ifstream estimatedBendPFile;
  estimatedBendPFile.open("EstimatedBendingMomentum.txt", ifstream::in);
  ifstream errEstimatedBendPFile;
  errEstimatedBendPFile.open("ErrEstimatedBendingMomentum.txt", ifstream::in);
  ifstream errVertexEstimatedBendPFile;
  errVertexEstimatedBendPFile.open("ErrVertexEstimatedBendingMomentum.txt", ifstream::in);
  ifstream reconstructedBendPFile;
  reconstructedBendPFile.open("ReconstructedBendingMomentum.txt", ifstream::in);
  ifstream errReconstructedBendPFile;
  errReconstructedBendPFile.open("ErrReconstructedBendingMomentum.txt", ifstream::in);
  
  Double_t nonBendSlope, nonBendImpact, estimatedBendP, errEstimatedBendP, errVertexEstimatedBendP, reconstructedBendP, errReconstructedBendP;
  Double_t simuBendP1, simuBendP2, simuBendP3, simuBendPTot, simuBendPReco;
  
  TH1F *hNonBendSlope = new TH1F("hNonBendSlope", "non bending slope", 200, -1, 1);
  TH1F *hNonBendingImpact = new TH1F("hNonBendingImpact", "non bending impact param", 400, -200, 200);
  TH1F *hEstimatedBendP = new TH1F("hEstimatedBendP", "estimated bending momentum", 2000, -100, 100);
  TH1F *hReconstructedBendP = new TH1F("hReconstructedBendP", "reconstructed bending momentum", 2000, -100, 100);
  TH1F *hErrEstimatedBendP = new TH1F("hErrEstimatedBendP", "estimated error/reconstructed bending momentum", 200, -0.1, 0.1);
  TH1F *hDiffBendP = new TH1F("hDiffBendP", "(estimated - reconstructed)/reconstructed bending momentum", 200, -1, 1);
  TH1F *hSimuDiffBendP1 = new TH1F("hSimuDiffBendP1", "(simulated - reconstructed)/reconstructed bending momentum", 200, -1, 1);
  TH1F *hSimuDiffBendP2 = new TH1F("hSimuDiffBendP2", "(simulated - reconstructed)/reconstructed bending momentum", 200, -1, 1);
  TH1F *hSimuDiffBendP3 = new TH1F("hSimuDiffBendP3", "(simulated - reconstructed)/reconstructed bending momentum", 200, -1, 1);
  TH1F *hSimuDiffBendTot = new TH1F("hSimuDiffBendPTot", "(simulated - reconstructed)/reconstructed bending momentum", 200, -1, 1);
  TH1F *hSimuDiffBendPReco = new TH1F("hSimuDiffBendPReco", "(simulated - reconstructed)/reconstructed bending momentum", 200, -1, 1);
  
  while (!reconstructedBendPFile.eof()) {
    nonBendSlopeFile>>nonBendSlope;
    nonBendImpactFile>>nonBendImpact;
    estimatedBendPFile>>estimatedBendP;
    errEstimatedBendPFile>>errEstimatedBendP;
    errVertexEstimatedBendPFile>>errVertexEstimatedBendP;
    reconstructedBendPFile>>reconstructedBendP;
    errReconstructedBendPFile>>errReconstructedBendP;
    
    
    hNonBendSlope->Fill(nonBendSlope);
    hNonBendingImpact->Fill(nonBendImpact);
    hEstimatedBendP->Fill(estimatedBendP);
    hErrEstimatedBendP->Fill(errEstimatedBendP/reconstructedBendP);
    hReconstructedBendP->Fill(reconstructedBendP);
    hDiffBendP->Fill((estimatedBendP-reconstructedBendP)/reconstructedBendP);
    
    for (Int_t i=0; i<100; i++) {
      simuBendP1 = 0.98 * (reconstructedBendP + gRandom->Gaus(0.,errEstimatedBendP));
      simuBendP2 = 0.98 * (reconstructedBendP + gRandom->Gaus(0., errVertexEstimatedBendP));
      simuBendP3 = 0.98 * (reconstructedBendP + gRandom->Gaus(0., 0.07*reconstructedBendP));
      simuBendPTot = 0.98 * (reconstructedBendP + gRandom->Gaus(0.,errEstimatedBendP) +
//                          gRandom->Gaus(0., errVertexEstimatedBendP) +
                          gRandom->Gaus(0., 0.07*reconstructedBendP));
      simuBendPReco = 0.98 * (reconstructedBendP + gRandom->Gaus(0.,errReconstructedBendP));
      
      hSimuDiffBendP1->Fill((simuBendP1-reconstructedBendP)/reconstructedBendP);
      hSimuDiffBendP2->Fill((simuBendP2-reconstructedBendP)/reconstructedBendP);
      hSimuDiffBendP3->Fill((simuBendP3-reconstructedBendP)/reconstructedBendP);
      hSimuDiffBendPTot->Fill((simuBendPTot-reconstructedBendP)/reconstructedBendP);
      hSimuDiffBendPReco->Fill((simuBendPReco-reconstructedBendP)/reconstructedBendP);
    }
  }
  
  TCanvas *c = new TCanvas();
  c->Divide(3,2);
  c->cd(1);
  hNonBendSlope->Draw();
  c->cd(2);
  hNonBendingImpact->Draw();
  c->cd(3);
  hEstimatedBendP->Draw();
  c->cd(4);
  hReconstructedBendP->Draw();
  c->cd(5);
  hDiffBendP->Draw();
  hSimuDiffBendP1->SetLineColor(2);
  hSimuDiffBendP1->Scale(0.01);
  hSimuDiffBendP1->Draw("same");
//  hSimuDiffBendP2->SetLineColor(3);
//  hSimuDiffBendP2->Scale(0.01);
//  hSimuDiffBendP2->Draw("same");
  hSimuDiffBendP3->SetLineColor(6);
  hSimuDiffBendP3->Scale(0.01);
  hSimuDiffBendP3->Draw("same");
//  hSimuDiffBendPReco->SetLineColor(5);
//  hSimuDiffBendPReco->Scale(0.01);
//  hSimuDiffBendPReco->Draw("same");
  hSimuDiffBendPTot->SetLineColor(4);
  hSimuDiffBendPTot->Scale(0.01);
  hSimuDiffBendPTot->Draw("same");
  c->cd(6);
  hErrEstimatedBendP->Draw();
  
  
}
