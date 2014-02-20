#include <Riostream.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TStopwatch.h>
#include <TClonesArray.h>

#include <AliMUONClusterInfo.h>

//MUONClusterInfo.C before

void residuals(char *filein="clusterInfo.root") {
  TFile *fin= new TFile(filein);
  TTree *treeClInfo=(TTree*)fin->Get("clusterInfoTree");

  TClonesArray histosResX("TH1F",11);
  TClonesArray histosResY("TH1F",11);
  Int_t nbin=100;
  Double_t maxX=0.5, minX=-maxX; //cm
  Double_t maxY=0.07,minY=-maxY; //cm

  for(Int_t i=0;i<11;i++){
    char nameX[20], nameY[20];
    char titleX[50],titleY[50];
    if(i==10){
      sprintf(nameX,"hresXchAll");
      sprintf(titleX,"Cluster-track Residuals: X (all);x (cm)");
      sprintf(nameY,"hresYchAll");
      sprintf(titleY,"Cluster-track Residuals: Y (all);y (cm)");
    }else {
      sprintf(nameX,"hresXch%d",i);
      sprintf(titleX,"Cluster-track Residuals: X (ch%d);x (cm)",i);
      sprintf(nameY,"hresYch%d",i);
      sprintf(titleY,"Cluster-track Residuals: Y (ch%d);y (cm)",i);
    }

    new (histosResX[i]) TH1F(nameX,titleX,nbin,minX,maxX);
    new (histosResY[i]) TH1F(nameY,titleY,nbin,minY,maxY);

  }

  TStopwatch watch;
  watch.Start();
  Int_t entries=(Int_t)treeClInfo->GetEntries();
  cout<<"Entries = "<<entries<<endl;

  AliMUONClusterInfo *aliClInfo=0;

  TBranch *brClInfo=treeClInfo->GetBranch("clusterInfo");
 
  brClInfo->SetAddress(&aliClInfo);

  for(Int_t i=0;i<entries; i++){

    //cout<<"Dentro al for "<<i<<"\t";

    brClInfo->GetEntry(i);

    //cout<<aliClInfo->GetChamberId()<<"\t"<<aliClInfo->GetClusterX()<<"\t"<<aliClInfo->GetTrackX()<<"\t"<<aliClInfo->GetClusterY()<<"\t"<<aliClInfo->GetTrackY()<<endl;

    ((TH1F*)histosResX[10])->Fill(aliClInfo->GetClusterX()-aliClInfo->GetTrackX());
    ((TH1F*)histosResY[10])->Fill(aliClInfo->GetClusterY()-aliClInfo->GetTrackY());

    Int_t chId=aliClInfo->GetChamberId();

    switch (chId){
    case 0:
      
      ((TH1F*)histosResX[0])->Fill(aliClInfo->GetClusterX()-aliClInfo->GetTrackX());
      ((TH1F*)histosResY[0])->Fill(aliClInfo->GetClusterY()-aliClInfo->GetTrackY());
      break;
     

    case 1:
      
      ((TH1F*)histosResX[1])->Fill(aliClInfo->GetClusterX()-aliClInfo->GetTrackX());
      ((TH1F*)histosResY[1])->Fill(aliClInfo->GetClusterY()-aliClInfo->GetTrackY());
      break;

    case 2:
      
      ((TH1F*)histosResX[2])->Fill(aliClInfo->GetClusterX()-aliClInfo->GetTrackX());
      ((TH1F*)histosResY[2])->Fill(aliClInfo->GetClusterY()-aliClInfo->GetTrackY());
      break;
    case 3:
     
      ((TH1F*)histosResX[3])->Fill(aliClInfo->GetClusterX()-aliClInfo->GetTrackX());
      ((TH1F*)histosResY[3])->Fill(aliClInfo->GetClusterY()-aliClInfo->GetTrackY());
      break;

    case 4:
      
      ((TH1F*)histosResX[4])->Fill(aliClInfo->GetClusterX()-aliClInfo->GetTrackX());
      ((TH1F*)histosResY[4])->Fill(aliClInfo->GetClusterY()-aliClInfo->GetTrackY());
      break;
    case 5:
      
      ((TH1F*)histosResX[5])->Fill(aliClInfo->GetClusterX()-aliClInfo->GetTrackX());
      ((TH1F*)histosResY[5])->Fill(aliClInfo->GetClusterY()-aliClInfo->GetTrackY());
      break;
    case 6:
      
      ((TH1F*)histosResX[6])->Fill(aliClInfo->GetClusterX()-aliClInfo->GetTrackX());
      ((TH1F*)histosResY[6])->Fill(aliClInfo->GetClusterY()-aliClInfo->GetTrackY());
      break; 
    case 7:
      
      ((TH1F*)histosResX[7])->Fill(aliClInfo->GetClusterX()-aliClInfo->GetTrackX());
      ((TH1F*)histosResY[7])->Fill(aliClInfo->GetClusterY()-aliClInfo->GetTrackY());
      break; 
    case 8:
      
      ((TH1F*)histosResX[8])->Fill(aliClInfo->GetClusterX()-aliClInfo->GetTrackX());
      ((TH1F*)histosResY[8])->Fill(aliClInfo->GetClusterY()-aliClInfo->GetTrackY());
      break; 
    case 9:
      
      ((TH1F*)histosResX[9])->Fill(aliClInfo->GetClusterX()-aliClInfo->GetTrackX());
      ((TH1F*)histosResY[9])->Fill(aliClInfo->GetClusterY()-aliClInfo->GetTrackY());
      break; 

    default:
      cout<<"Default"<<endl;
      break;
      //cout<<"chId = "<<chId<<"\t\t";

    }
    
  }

  Double_t sigmaX[11],sigmaY[11],sigmaErrX[11],sigmaErrY[11];
  Double_t var[11];

  //gaussian fit
  for(Int_t i=0;i<11;i++){
    var[i]=i;
    cout<<((TH1F*)histosResX[i])->GetName()<<endl;
    Double_t rmsX=((TH1F*)histosResX[i])->GetRMS();
    TF1 *gaussianX;
    gaussianX = new TF1("gaussianX","gaus",-2*rmsX,2*rmsX);
    ((TH1F*)histosResX[i])->Fit("gaussianX","R");
    sigmaX[i]=gaussianX->GetParameter(2);
    sigmaErrX[i]=gaussianX->GetParError(2);
    cout<<"sigma fit = "<<sigmaX[i]<<" cm"<<endl<<endl;
    
    cout<<((TH1F*)histosResY[i])->GetName()<<endl;
    Double_t rmsY=((TH1F*)histosResY[i])->GetRMS();
    TF1 *gaussianY = new TF1("gaussianY","gaus",-2*rmsY,2*rmsY);
    ((TH1F*)histosResY[i])->Fit("gaussianY","R");
    sigmaY[i]=gaussianY->GetParameter(2);
    sigmaErrY[i]=gaussianY->GetParError(2);
    cout<<"sigma fit = "<<sigmaY[i]<<" cm"<<endl<<endl;
  }
  

  watch.Stop();
  watch.Print();

  TH1F *hSigX=new TH1F("hSigX","#sigma_{fit} (X);Chamber number;res (cm)",10,0.,10.);
  TH1F *hSigY=new TH1F("hSigY","#sigma_{fit} (Y);Chamber number;res (cm)",10,0.,10.);

  for(Int_t i=0;i<10;i++){
    hSigX->Fill(var[i],sigmaX[i]);
    hSigY->Fill(var[i],sigmaY[i]);

  }

  hSigX->SetError(sigmaErrX);
  hSigY->SetError(sigmaErrY);
  (hSigX->GetYaxis())->SetRangeUser(0.,0.09);
  (hSigY->GetYaxis())->SetRangeUser(0.,0.008);
  hSigX->SetMarkerStyle(3);
  hSigY->SetMarkerStyle(3);

  /*
  TCanvas *a=new TCanvas("a","allX");
  a->cd();
  hresXall->Draw();

  TCanvas *b=new TCanvas("b","allY");
  b->cd();
  hresYall->Draw();


  TCanvas *c1=new TCanvas("c1","ch0");
  c1->Divide(1,2);
  c1->cd(1);
  hresXch0->Draw();
  c1->cd(2);
  hresYch0->Draw();

  TCanvas *c2=new TCanvas("c2","ch1");
  c2->Divide(1,2);
  c2->cd(1);
  hresXch1->Draw();
  c2->cd(2);
  hresYch1->Draw();

  TCanvas *c3=new TCanvas("c3","ch2");
  c3->Divide(1,2);
  c3->cd(1);
  hresXch2->Draw();
  c3->cd(2);
  hresYch2->Draw();

  TCanvas *c4=new TCanvas("c4","ch3");
  c4->Divide(1,2);
  c4->cd(1);
  hresXch3->Draw();
  c4->cd(2);
  hresYch3->Draw();

  TCanvas *c5=new TCanvas("c5","ch4");
  c5->Divide(1,2);
  c5->cd(1);
  hresXch4->Draw();
  c5->cd(2);
  hresYch4->Draw();

  TCanvas *c6=new TCanvas("c6","ch5");
  c6->Divide(1,2);
  c6->cd(1);
  hresXch5->Draw();
  c6->cd(2);
  hresYch5->Draw();

  TCanvas *c7=new TCanvas("c7","ch6");
  c7->Divide(1,2);
  c7->cd(1);
  hresXch6->Draw();
  c7->cd(2);
  hresYch6->Draw();

  TCanvas *c8=new TCanvas("c8","ch7");
  c8->Divide(1,2);
  c8->cd(1);
  hresXch7->Draw();
  c8->cd(2);
  hresYch7->Draw();

  TCanvas *c9=new TCanvas("c9","ch8");
  c9->Divide(1,2);
  c9->cd(1);
  hresXch8->Draw();
  c9->cd(2);
  hresYch8->Draw();

  TCanvas *c10=new TCanvas("c10","ch9");
  c10->Divide(1,2);
  c10->cd(1);
  hresXch9->Draw();
  c10->cd(2);
  hresYch9->Draw();

  */

  TCanvas *c1=new TCanvas("c1","X");
  c1->cd();
  hSigX->Draw("P");

  TCanvas *c2=new TCanvas("c2","Y");
  c2->cd();
  hSigY->Draw("P");


  TFile *f=new TFile("Resolution.root","recreate");
  f->cd();

  for(Int_t i=0;i<11;i++){
    histosResX[i]->Write();
    histosResY[i]->Write();
  }
  hSigX->Write();
  hSigY->Write();

  f->Close();
  cout<<"Plots in Resolution.root"<<endl;

}
