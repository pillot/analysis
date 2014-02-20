///==========================================================================
///
///    macro to plot centrality bin values 
///==========================================================================
///
#include <cstdlib>

const Int_t nbins;

//double centPercent[]={5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,99.99}; 
//double centPercent[]={10.,20.,30.,40.,50.,60.,70.,80.,90.,99.99}; 
//double centPercent[]={20.,40.,60.,80.,99.99};
//double centPercent[]={10.,20.,40.,80.,99.99}; 
//double centPercent[]={50.,99.99};
//double centPercent[]={80.,99.99};
//double centPercent[]={10., 20., 50., 80., 99.99};
//double centPercent[]={20.,40.,80.,99.99};
//double centPercent[]={10., 20., 30., 50., 80., 99.99};
//double centPercent[]={20.,60.,80.,99.99};
//double centPercent[]={50.,99.99};
//double centPercent[]={50.,80.,99.99};
double centPercent[]={60.,80.,99.99};

nbins = sizeof(centPercent)/sizeof(double);

TArrayI* binUp = new TArrayI(nbins);
TArrayF* Multbin = new TArrayF(nbins);

void AnalyzeCentralityB(Bool_t weight = kTRUE, Int_t mysignn=64,Int_t mymind=4,Int_t myr=662,Int_t mya=546)    
// mysignn=64+-5,Int_t mymind=4+-4,Int_t myr=662+-6,Int_t mya=546+-10)    
{
 TGraphErrors *gNpart=new TGraphErrors(0);
 gNpart->SetName("gNpart"); 
 TGraphErrors *gNcoll=new TGraphErrors(0);
 gNcoll->SetName("gNcoll"); 
 TGraphErrors *gtAA=new TGraphErrors(0);
 gtAA->SetName("gtAA"); 


  TString suffix=Form("GlauberMC_PbPb_ntuple_sigma%d_mind%d_r%d_a%d.root",mysignn,mymind,myr,mya);
  TFile *f=new TFile(Form("/Users/pillot/Work/Alice/Work/Macros/MuonPhysics/centrality/%s",suffix.Data()));
  
  TNtuple *nt=(TNtuple*)f->Get("nt_Pb_Pb");
  
  Float_t B;
  Float_t Npart;
  Float_t Ncoll;
  Float_t tAA;
  
  nt->SetBranchAddress("B",&B);
  nt->SetBranchAddress("Npart",&Npart);
  nt->SetBranchAddress("Ncoll",&Ncoll);
  nt->SetBranchAddress("tAA",&tAA);


  Int_t totev=nt->GetEntries();
 
 
  TH1F*hB=new TH1F("hB","hB",100000,0.,5000.);
  TH1F*hNpart=new TH1F("hNpart","hNpart",1000000,-0.5,999999.5);
  TH1F*hNcoll=new TH1F("hNcoll","hNcoll",100000,-0.5,99999.5);
  TH1F*htAA=new TH1F("htAA","htAA",100000,-0.5,39.5);
 
  
  for (Int_t iEvent = 0; iEvent < totev; iEvent++) {
    nt->GetEvent(iEvent);    
      if (weight) hNpart->Fill(Ncoll*Npart);
      else hNpart->Fill(Npart);
      hNcoll->Fill(Ncoll);
      hB->Fill(B);
      htAA->Fill(tAA);
     } 

 //---------------------------------------------------
  getCentrality(hB);
 //---------------------------------------------------


for (int j=0;j<nbins;j++) cout<<"bin di centralità n. "<<j<<" relativo alla centralità "<<centPercent[j]<<"B "<<Multbin->At(j)<<endl; 
 
 TH1F* hnpartcutb[nbins];
 char histtitp[100];
 for (Int_t i=0; i<binUp->GetSize(); i++) {
   sprintf(histtitp,"npartcutb%d",i);
   hnpartcutb[i] = new TH1F(histtitp,histtitp,1000000,-0.5,999999.5);
   hnpartcutb[i]->SetLineWidth(1);
   hnpartcutb[i]->SetStats(1);

 }

 TH1F* hncollcutb[nbins];
 char histtitc[100];
 for (Int_t i=0; i<binUp->GetSize(); i++) {
   sprintf(histtitc,"ncollcutb%d",i);
   hncollcutb[i] = new TH1F(histtitc,histtitc,10000,-0.5,9999.5);
   hncollcutb[i]->SetLineWidth(1);
   hncollcutb[i]->SetStats(1);

 }

 TH1F* htaacutb[nbins];
 char histtitt[100];
 for (Int_t i=0; i<binUp->GetSize(); i++) {
   sprintf(histtitt,"taacutb%d",i);
   htaacutb[i] = new TH1F(histtitt,histtitt,100000,-0.5,39.5);
   htaacutb[i]->SetLineWidth(1);
   htaacutb[i]->SetStats(1);

 }



 for (Int_t iEvent = 0; iEvent < totev; iEvent++) {
   nt->GetEvent(iEvent);
   for (int ibin=0; ibin<nbins; ibin++) {
     if (B<Multbin->At(ibin)) {
       if (weight) hnpartcutb[ibin]->Fill(Ncoll*Npart);
       else hnpartcutb[ibin]->Fill(Npart);
       hncollcutb[ibin]->Fill(Ncoll); 
       htaacutb[ibin]->Fill(tAA); 
       break;
     } 
   } 
 } 
 

   


 // superimpose cut histograms
// new TCanvas();
// hnpartcutb[0]->Draw("");
 for (Int_t icentr=0; icentr<nbins;icentr++) {
   hnpartcutb[icentr]->SetLineColor(icentr);
   hnpartcutb[icentr]->Draw("same");
 } 

 for (Int_t icentr=0; icentr<nbins;icentr++) {
   new TCanvas();
   hnpartcutb[icentr]->SetLineColor(icentr+1);
   hnpartcutb[icentr]->SetStats(1); 
   hnpartcutb[icentr]->Draw("");
   cout<<" Centrality bin "<<icentr<<" Mean N_{part}--> "<<hnpartcutb[icentr]->GetMean()<<" rms--> "<<hnpartcutb[icentr]->GetRMS()<<endl;
   gNpart->SetPoint(icentr,Float_t(icentr),hnpartcutb[icentr]->GetMean());
   gNpart->SetPointError(icentr,0,hnpartcutb[icentr]->GetRMS());
 }

 for (Int_t icentr=0; icentr<nbins;icentr++) {
   new TCanvas();
   hncollcutb[icentr]->SetLineColor(icentr+1);
   hncollcutb[icentr]->SetStats(1); 
   hncollcutb[icentr]->Draw("");
   cout<<" Centrality bin "<<icentr<<" Mean N_{coll}--> "<<hncollcutb[icentr]->GetMean()<<" rms--> "<<hncollcutb[icentr]->GetRMS()<<endl;
   gNcoll->SetPoint(icentr,Float_t(icentr),hncollcutb[icentr]->GetMean());
   gNcoll->SetPointError(icentr,0,hncollcutb[icentr]->GetRMS());
 }

 for (Int_t icentr=0; icentr<nbins;icentr++) {
   new TCanvas();
   htaacutb[icentr]->SetLineColor(icentr+1);
   htaacutb[icentr]->SetStats(1); 
   htaacutb[icentr]->Draw("");
   cout<<" Centrality bin "<<icentr<<" Mean t_{AA}--> "<<htaacutb[icentr]->GetMean()<<" rms--> "<<htaacutb[icentr]->GetRMS()<<endl;
   gtAA->SetPoint(icentr,Float_t(icentr),htaacutb[icentr]->GetMean());
   gtAA->SetPointError(icentr,0,htaacutb[icentr]->GetRMS());
 }



 // TString suffixhisto=Form("/home/atoia/GlauberNtuple/GlauberMC_PbPb_histoB_sigma%d_mind%d_r%d_a%d_%d.root",mysignn,mymind,myr,mya,nbins);
 //  const Char_t* filehistoname=suffixhisto.Data();
 //TFile*filefinal=new TFile(filehistoname,"recreate");
 TFile*filefinal=new TFile(Form("GlauberInfo_sigma%d_mind%d_r%d_a%d.root",mysignn,mymind,myr,mya),"recreate");
 
 gNpart->Write();
 gNcoll->Write();
 gtAA->Write();
 
 hNpart->Write();
 hNcoll->Write();
 hB    ->Write();
 htAA  ->Write();
  
 for (int i=0; i<nbins; i++) {
   hnpartcutb[i]->Write();
   hncollcutb[i]->Write();
   htaacutb[i]  ->Write();
 }

 filefinal->Close();

}


void getCentrality(TH1 *histNch, Float_t ff=1.0)
// histNch - histo of multiplicity distribution (better with binsize=1)x
// ff fraction of accepted events. All losses are assumed to occur in most
// peripheral bin
{

 //double sum= histNch->GetEntries() - histNch->GetBinContent(1);
 double sum= histNch->Integral(); 
 int nbinsx=histNch->GetNbinsX();
 double frac=0.;
 int ic=0;
 for (int ib=1;ib<=nbinsx;ib++){
// for (int ib=nbinsx;ib>0;ib--){
   frac += histNch->GetBinContent(ib)/sum*100.*ff;
   if(frac > centPercent[ic]){
     binUp->SetAt(ib,ic);
     Multbin->SetAt(histNch->GetBinCenter(ib),ic);
//     cout<<" centrality="<<centPercent[ic]<<"   mult <="<< histNch->GetBinCenter(ib) <<endl;
     cout<<" centrality="<<centPercent[ic]<<" impact parameter <="<< histNch->GetBinCenter(ib) <<endl;
     ic++;
   }
   if(ic==nbins) break;
 }
 
 printf(" \n float binUp[%i] = {",nbins);  
 // cout <<" \n float multCent[nbins] = {";

 for (int ic=nbins-1; ic>-1; ic--){
   cout<< binUp->At(ic);
   if (ic!=0) cout<<", ";
 }
 cout<<"};\n"<<endl;


 printf(" \n float multCent[%i] = {",nbins);  
 // cout <<" \n float multCent[nbins] = {";

 for (int ic=nbins-1; ic>-1; ic--){
   cout<< histNch->GetBinCenter(binUp->At(ic));
   if (ic!=0) cout<<", ";
 }
 cout<<"};\n"<<endl;
}
