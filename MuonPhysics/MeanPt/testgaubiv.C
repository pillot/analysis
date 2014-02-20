void testgaubiv(Double_t xmean =4.32592, Double_t ymean=3.61705, Double_t sigmax =0.261907, Double_t sigmay=0.205432, Double_t corr=0.96904){
  
// Parameters from Philippe's macro
// p0=4.32592+/-0.261907 
// n=3.61705+/-0.205432
// corr. = 0.96904
// which gave <pT> = 2.7766 +/- 0.0442

TF2 *fgau = new TF2("fgau","1./(2*3.141*[2]*[3]*sqrt(1-[4]*[4]))*exp((-1./(2*(1-[4]*[4])))*((x-[0])*(x-[0])/([2]*[2])+(y-[1])*(y-[1])/([3]*[3])-2*[4]*(x-[0])/[2]*(y-[1])/[3]))",3.,6.,2.,5.);

fgau->SetParameters(xmean, ymean, sigmax, sigmay, corr);

TCanvas *c1 = new TCanvas("c1","c1");
fgau->Draw("surf2");


TH2D *hgau = new TH2D("hgau","hgau",300,3.,6.,300,2.,5.);
hgau->FillRandom("fgau",10000000);
TCanvas *c2 = new TCanvas("c2","c2");

hgau->Draw("surf");

TH1D *hptmedio = new TH1D("hptmedio","hptmedio",1000,2.,4.);
TF1* fpt = new TF1("fpt","x/((1+(x/[0])**2)**[1])",0.,15.);
  for(Int_t i=0;i<hgau->GetNbinsX();i++){
   for(Int_t j=0;j<hgau->GetNbinsY();j++){
     Double_t weight = hgau->GetBinContent(i+1,j+1);
     Double_t p0 = 3.+i*0.01+0.005;
     Double_t n = 2.+j*0.01+0.005;
//     printf("p0, n = %f %f\n",p0,n);
     fpt->SetParameters(p0,n);
     hptmedio->Fill(fpt->Mean(0.,15.),weight);
   }
  }
  
TCanvas *c3 = new TCanvas("c3","c3");
 
hptmedio->Fit("gaus");
hptmedio->Draw();
printf ("Average pT from bivariate normal density : %5.4f +/- %5.4f (histo)\n",hptmedio->GetMean(),hptmedio->GetRMS());
printf ("Average pT from bivariate normal density : %5.4f +/- %5.4f (fit)\n",hptmedio->GetFunction("gaus")->GetParameter(1),hptmedio->GetFunction("gaus")->GetParameter(2));

}