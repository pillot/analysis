#include <cmath>
#include <vector>

#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TString.h>
#include <TH2D.h>
#include <TH1D.h>

//define all used histograms 
std::vector<TH2D*> h2chi2_ndf;
std::vector<TH1D*> hprob;
std::vector<TH1D*> hk3x;
std::vector<TH1D*> hk3y;
std::vector<TH1D*> hAsymm;

//_________________________________________________________________________________________________
void LoadHist()
{
  //for naming histograms
  int station[3] = {1, 2, 345};

  // common histograms
  for (int i = 0; i < 3; ++i) {
    h2chi2_ndf.push_back(new TH2D(Form("h2chi2_ndf_%d", i), Form("#chi^{2} vs ndf (St.%d)",station[i]), 16, -0.5 , 15.5, 1000, 0, 50));
    h2chi2_ndf[i]->SetDirectory(0);
    h2chi2_ndf[i]->GetXaxis()->SetTitle("ndf");
    h2chi2_ndf[i]->GetYaxis()->SetTitle("#chi^{2}");
  }

  for (int i = 0; i < 3; ++i) {
    hprob.push_back(new TH1D(Form("hprob_%d", i), Form("P-Value (St.%d)",station[i]), 1000, 0 , 1));
    hprob[i]->SetDirectory(0);
    hprob[i]->GetXaxis()->SetTitle("p-value");
  }

  for (int i = 0; i < 3; ++i) {
    hk3x.push_back(new TH1D(Form("hk3x_%d", i), Form("K_{3X} (St.%d)",station[i]), 1000, 0 , 1));
    hk3x[i]->SetDirectory(0);
    hk3x[i]->GetXaxis()->SetTitle("K_{3X}");
  }

  for (int i = 0; i < 3; ++i) {
    hk3y.push_back(new TH1D(Form("hk3y_%d", i), Form("K_{3Y} (St.%d)",station[i]), 1000, 0 , 1));
    hk3y[i]->SetDirectory(0);
    hk3y[i]->GetXaxis()->SetTitle("K_{3Y}");
  }

  for (int i = 0; i < 6; ++i) {
    hAsymm.push_back(new TH1D(Form("hAsymm_%d", i), Form("Asymmetry bending - nbending (St.%d)",station[i / 2]), 700, -0.7 , 0.7));
    hAsymm[i]->SetDirectory(0);
    hAsymm[i]->GetXaxis()->SetTitle("(NB-B)/(NB+B)");
  }

}

//_________________________________________________________________________________________________

void DelHist() {

  for (auto& hist : h2chi2_ndf) delete hist;
  h2chi2_ndf.clear();

  for (auto& hist : hprob) delete hist;
  hprob.clear();

  for (auto& hist : hk3x) delete hist;
  hk3x.clear();

  for (auto& hist : hk3y) delete hist;
  hk3y.clear();

  for (auto& hist : hAsymm) delete hist;
  hAsymm.clear();
}

//_________________________________________________________________________________________________
void FillInfoGraphErrAsymm(TGraphAsymmErrors* graph, const TString& titleX, const TString& titleY, const TString& Station, const TString& Cathode, bool primeFile) 
{
    TString file = primeFile ? "(1)" : "(2)";
    TString title = Form("%s vs %s [%s] [%s] %s",
                        titleY.Data(), titleX.Data(),
                        Station.Data(), Cathode.Data(),
                        file.Data());

    graph->SetTitle(title);
    graph->SetName(title);
    graph->GetXaxis()->SetLimits(0, 5001);
    graph->GetXaxis()->SetTitle(titleX);
    graph->GetYaxis()->SetTitle(titleY);
    graph->SetMarkerColor((primeFile) ? kRed : kBlue);
    graph->SetMarkerStyle(47);
    graph->SetMinimum(0);
    graph->SetMaximum(100);
}
//_________________________________________________________________________________________________
void FillInfoGraph(TGraph* graph, const TString& titleX, const TString& titleY, const TString& Station, const TString& Cathode, bool primeFile) 
{
    TString file = primeFile ? "(1)" : "(2)";
    TString title = Form("%s vs %s [%s] [%s] %s",
                        titleY.Data(), titleX.Data(),
                        Station.Data(), Cathode.Data(),
                        file.Data());

    graph->SetTitle(title);
    graph->SetName(title);
    graph->GetXaxis()->SetLimits(0, 5001);
    graph->GetXaxis()->SetTitle(titleX);
    graph->GetYaxis()->SetTitle(titleY);
    graph->SetMarkerColor((primeFile) ? kRed : kBlue);
    graph->SetMarkerStyle(47);
    if(titleY == "#mu") {
      graph->SetMinimum(-20);
      graph->SetMaximum(20);
    } else {
      graph->SetMinimum(0);
      graph->SetMaximum(100);
    }
}
//_________________________________________________________________________________________________
void FillInfoHist(TH1D* h, const TString& titleX, const TString& titleY, const TString& Station, const TString& Cathode, bool primeFile, bool normalize = false) 
{
    TString oldName = h->GetName();
    TString file = primeFile ? "(1)" : "(2)";
    TString title = Form("%s vs %s %s [%s] %s %s",
                        titleY.Data(), titleX.Data(),
                        oldName.Data(), Station.Data(), Cathode.Data(), file.Data());

    if (normalize && (h->Integral() != 0)) {
        h->Scale(1.0 / h->Integral());
    }

    h->SetTitle(title);
    h->SetName(title);

    h->GetXaxis()->SetTitle(titleX);
    h->GetYaxis()->SetTitle(titleY);
    h->SetLineColor((primeFile) ? kRed : kBlue);
}
//_________________________________________________________________________________________________
TF1* plotNoise(std::string sigma, double alpha, double gamma = 0., bool Asymm = false)
{
    TF1 *Func = nullptr;

    if (sigma == "MC") {
      if(!Asymm){
        if (!Func) {
            Func = new TF1("func", 
                          [](double *x, double *par) { 
                              double charge = x[0];
                              double result = (std::pow(charge / par[1], 1. / par[2]) + par[0]);
                              return (0.5 * (sqrt(std::round(result)) + par[4])); 
                          }, 
                          0, 5000, 4); 
        }
        double signalParam[3] = {14., 13., 1.5};
        Func->SetParameters(signalParam[0], signalParam[1], signalParam[2], alpha); 
        Func->SetLineColor(4); 
        Func->SetLineWidth(3);   
        Func->SetLineStyle(1);
        Func->SetNpx(1000); 
      } else {
        if (!Func) {
          Func = new TF1("func", 
                        [](double *x, double *par) { 
                            double charge = x[0];
                            double result = (std::pow(charge / par[1], 1. / par[2]) + par[0]);
                            return sqrt(0.25 * (std::round(result) + par[3]) + 0.25 * charge * charge * (TMath::Exp(8 * par[4] * par[4]) - TMath::Exp(4 * par[4] * par[4]))); 
                        }, 
                        0, 5000, 5); 
        }
        double sigmaAsymm = gamma * 0.055;
        double signalParam[3] = {14., 13., 1.5};
        Func->SetParameters(signalParam[0], signalParam[1], signalParam[2], alpha, sigmaAsymm); 
        Func->SetLineColor(4); 
        Func->SetLineWidth(3);   
        Func->SetLineStyle(1);
        Func->SetNpx(1000); 
      }
    } else if (sigma == "sADC") {
      if(!Asymm) {
        if (!Func) {
            Func = new TF1("func", [](double *x, double *par) { return par[0] * sqrt(x[0]); }, 0, 1500, 1); 
        }

        alpha = (alpha < 0.001 ) ?  1. : alpha ;
        Func->SetParameters(alpha);
        Func->SetLineColor(4); 
        Func->SetLineWidth(3);             
        Func->SetLineStyle(1);
        Func->SetNpx(1000); 
      } else {
        if (!Func) {
          Func = new TF1("func", 
                        [](double *x, double *par) { 
                            double charge = x[0];
                            double intrinsic = par[0] * sqrt(charge);
                            return sqrt(par[0] * par[0] * sqrt(charge) * sqrt(charge) + 0.25 * charge * charge * (TMath::Exp(8 * par[1] * par[1]) - TMath::Exp(4 * par[1] * par[1]))); 
                        }, 
                        0, 5000, 2); 
        }
        double sigmaAsymm = gamma * 0.055;
        alpha = (alpha < 0.001 ) ?  1. : alpha ;
        Func->SetParameters(alpha, sigmaAsymm); 
        Func->SetLineColor(4); 
        Func->SetLineWidth(3);   
        Func->SetLineStyle(1);
        Func->SetNpx(1000); 
      }
    } else {
        std::cerr << "Unknown sigma type!" << std::endl;
        return nullptr; 
    }
    return Func;
  }

//_________________________________________________________________________________________________
void plotSAME(std::vector<TH1D*> hist, const char *name, const char *title, bool all = false) 
{
    TCanvas* c = new TCanvas(name, title, 800, 600);
    int N = hist.size();
    int COL = 1;
    if(all){COL = 2;}

    c->Divide(COL,N/2);
    for (int i = 1; i < (N/2 + 1); ++i) {
      c->cd(i);
      gPad->SetLogy(); 
      hist[2*(i-1)]->SetLineColor(2); 
      hist[2*(i-1)]->Draw("HIST");
      hist[2*(i-1) + 1]->Draw("HIST SAME");
    }
    c->Update();
    c->Draw();
}

//_________________________________________________________________________________________________
void plot1D(std::vector<TH1D*> hist, const char *name, const char *title, bool all = false) 
{
    TCanvas* c = new TCanvas(name, title, 800, 600);
    int N = hist.size();
    int COL = 1;
    if(all){COL = 2;}

    c->Divide(COL, N);  
    for (int i = 1; i < N+1; ++i) {
      c->cd(i);
      gPad->SetLogy(); 
      hist[i-1]->Draw("HIST");
    }
    c->Update();
    c->Draw();
}

//_________________________________________________________________________________________________ 
void plot2D(std::vector<TH2D*> hist, const char *name, const char *title, bool all = false) 
{
    TCanvas* c = new TCanvas(name, title, 800, 600);
    int N = hist.size();
    int COL = 1;
    if(all){ COL = 2;}  

    c->Divide(COL, N);
    for (int i = 1; i < N+1; ++i) {
      c->cd(i);
      gPad->SetLogz(); 
      gStyle->SetPalette(kRainBow);
      hist[i-1]->Draw("COLZ");
    }
    c->Update();
    c->Draw();
}

//_________________________________________________________________________________________________  
void tGraphErrAsymm(std::vector<TGraphAsymmErrors*>& graphs1, std::vector<TGraphAsymmErrors*>& graphs2, const std::string& name, int station)
{

  TCanvas* c = new TCanvas(name.c_str(), name.c_str(), 1200, 800);
  c->Divide(1, 2);

  auto func = plotNoise("sADC", 1.0); 

  for (int i = 0; i < 2; ++i) {
    c->cd(i+1);
    gPad->SetGrid();

    auto* g1 = graphs1[2*station + i];
    auto* g2 = graphs2[2*station + i];

    g1->Draw("AP");
    g2->Draw("P SAME");
    func->Draw("L SAME");

    TLegend* leg = new TLegend(0.65, 0.70, 0.88, 0.88);
    leg->AddEntry(g1, g1->GetTitle(), "p");
    leg->AddEntry(g2, g2->GetTitle(), "p");
    leg->AddEntry(func, "ToyMC th. noise", "l");
    leg->Draw();
  }

  c->Update();
  c->Draw();
}

//_________________________________________________________________________________________________  
void tGraph(std::vector<TGraph*>& graphs1, std::vector<TGraph*>& graphs2, const std::string& name, int station)
{

  TCanvas* c = new TCanvas(name.c_str(), name.c_str(), 1200, 800);
  c->Divide(2, 2);

  TGraph* plots[4] = {
      graphs1[2 * station], graphs2[2 * station],
      graphs1[2 * station + 1], graphs2[2 * station + 1]
  };

  for (int i = 0; i < 4; ++i) {
      c->cd(i + 1);
      gPad->SetGrid();

      plots[i]->Draw("AP");

      TLegend* leg = new TLegend(0.65, 0.70, 0.88, 0.88);
      leg->AddEntry(plots[i], plots[i]->GetTitle(), "p");
      leg->Draw();
  }

  c->Update();
  c->Draw();
}

//_________________________________________________________________________________________________  
void tHist(std::vector<TH1D*>& hists1, std::vector<TH1D*>& hists2, const std::string& name, int station, int cathode)
{
  std::vector<TH1D*> subh1(hists1.begin() + 8 * (2 * station + cathode), hists1.begin() + 8 * (2 * station + cathode + 1));
  std::vector<TH1D*> subh2(hists2.begin() + 8 * (2 * station + cathode), hists2.begin() + 8 * (2 * station + cathode + 1));

  TCanvas* c = new TCanvas(name.c_str(), name.c_str(), 1200, 800);
  c->Divide(2, 4);

  for (int i = 0; i < 8; ++i) {
        c->cd(i + 1);
        gPad->SetGrid();
        gPad->SetLogy();

        double max1 = subh1[i]->GetMaximum();
        double max2 = subh2[i]->GetMaximum();
        double ymax = std::max(max1, max2) * 1.5; 

        subh1[i]->GetXaxis()->SetRangeUser(-50, 50); 
        subh1[i]->SetMaximum(ymax);
        subh1[i]->Draw("HIST");
        subh2[i]->Draw("HIST SAME");

        TLegend* leg = new TLegend(0.65, 0.70, 0.88, 0.88);
        leg->AddEntry(subh1[i], subh1[i]->GetTitle(), "l");
        leg->AddEntry(subh2[i], subh2[i]->GetTitle(), "l");
        leg->Draw();
  }

  c->Update();
  c->Draw();
}

//_________________________________________________________________________________________________ 
void tRatio(std::vector<TGraphAsymmErrors*>& graphs1, std::vector<TGraphAsymmErrors*>& graphs2, const std::string& name, int station, const std::string& theory) 
{

  double alpha = std::stod(theory);
  auto lambda = [alpha](double charge) {
    return alpha * sqrt(charge);
  };

  TCanvas* c = new TCanvas(name.c_str(), name.c_str(), 1000, 800);
  c->Divide(1, 2);
  for (int i = 0; i < 2; ++i) {

    c->cd(i+1);
    gPad->SetGrid();

    TGraphAsymmErrors* gi = graphs1[2 * station + i];
    TGraphAsymmErrors* gj = graphs2[2 * station + i];
    TGraphAsymmErrors* ratio = new TGraphAsymmErrors();
    TGraphAsymmErrors* ratio_simple = new TGraphAsymmErrors(); //draw (#sigma_{DATA} / #sigma_{TMC})

    TString Cathode = (i == 0) ? "Bending" : "NonBending";
    TString title = Form("#sigma_{th}(#sigma_{DATA} / #sigma_{TMC}) - %s", Cathode.Data());
    TString title_simple = Form("(#sigma_{DATA} / #sigma_{TMC}) - %s", Cathode.Data());

    ratio->SetTitle(title);
    ratio_simple->SetTitle(title_simple);

    ratio->GetXaxis()->SetTitle("ADC");
    ratio->GetYaxis()->SetTitle("#sigma");
    ratio->SetMarkerStyle(2);
    ratio_simple->SetMarkerStyle(2);
    ratio_simple->SetMarkerColor(kGreen);

    auto func = plotNoise("sADC", 1.0); 

    int Ni = gi->GetN();
    for (int k = 0; k < Ni; ++k) {
        double xi, yi;
        gi->GetPoint(k, xi, yi);
        if (yi == 0) continue;
        double yi_err = 0.5 * (gi->GetErrorYlow(k) + gi->GetErrorYhigh(k));

        int Nj = gj->GetN();
        for (int l = 0; l < Nj; ++l) {
            double xj, yj;
            gj->GetPoint(l, xj, yj);
            if (std::abs(xi - xj) < 1e-6 && yj != 0) {
                double yj_err = 0.5 * (gj->GetErrorYlow(l) + gj->GetErrorYhigh(l));
                double r = yi / yj;
                double rel_err = std::sqrt(std::pow(yi_err / yi, 2) + std::pow(yj_err / yj, 2));
                double r_err = r * rel_err * lambda(xi);
                double r_err_simple = r * rel_err;

                ratio->AddPoint(xi, (r * lambda(xi)));
                ratio_simple->AddPoint(xi, r);
                ratio->SetPointError(ratio->GetN() - 1, 0, 0, r_err, r_err);
                ratio_simple->SetPointError(ratio_simple->GetN() - 1, 0, 0, r_err, r_err);
                break;
            }
        }
    }

    ratio->Draw("AP");
    ratio_simple->Draw("P SAME");
    func->Draw("L SAME");
    TLegend* leg = new TLegend(0.65, 0.70, 0.88, 0.88);
    leg->AddEntry(ratio, ratio->GetTitle(), "p");
    leg->AddEntry(ratio_simple, ratio_simple->GetTitle(), "p");
    leg->AddEntry(func, "ToyMC th. noise", "l");
    leg->Draw();

  }
  c->Update();
  c->Draw();
}