#include "AliAnalysisMuMu.h"

#include "AliAnalysisManager.h"
#include "AliCounterCollection.h"
#include "AliHistogramCollection.h"
#include "AliLog.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TDatabasePDG.h"
#include "TH1.h"
#include "TH2.h"
#include "THashList.h"
#include "TList.h"
#include "TMath.h"
#include "TObjString.h"
#include "TPaveText.h"
#include <algorithm>
#include "AliLog.h"

ClassImp(AliAnalysisMuMu)

namespace
{
  Int_t GetNbins(Double_t xmin, Double_t xmax, Double_t xstep)
  {
    if ( TMath::AreEqualRel(xstep,0.0,1E-9) ) return 1;
    
    return TMath::Nint(TMath::Abs((xmax-xmin)/xstep));
  }
  
}

//_____________________________________________________________________________
AliAnalysisMuMu::AliAnalysisMuMu() : AliAnalysisTaskSE("AliAnalysisMuMu"),
fHistogramCollection(0),
fOutput(0),
fTriggerClasses(0),
fEventCounters(0),
fIsFromESD(kFALSE),
fSingleTrackCutNames(0x0),
fPairTrackCutNames(0x0),
fCentralityLimits(),
fCentralityNames(0x0),
fGlobalEventSelectionNames(0x0),
fAA(kFALSE)
{
  /// default ctor
}

//_____________________________________________________________________________
AliAnalysisMuMu::AliAnalysisMuMu(TList* triggerClasses, Bool_t fromESD) 
: AliAnalysisTaskSE(Form("AliAnalysisMuMu-from%s",fromESD ? "ESD":"AOD")),
fHistogramCollection(0),
fOutput(0),
fTriggerClasses(new THashList),
fEventCounters(0),
fIsFromESD(fromESD),
fSingleTrackCutNames(0x0),
fPairTrackCutNames(0x0),
fCentralityLimits(),
fCentralityNames(new TObjArray),
fGlobalEventSelectionNames(0x0),
fAA(kFALSE)
{
  // Constructor
  // Define input and output slots here
  // Slot #0 works are defined in TaskSE
  // Output slot #1 writes into a TList container
  
  fTriggerClasses->SetOwner(kTRUE);
  
  DefineOutput(1, TList::Class());
  
  if ( !triggerClasses ) 
  {
    DefineDefaultTriggerClasses();
  }
  else
  {
    TIter next(triggerClasses);
    TObjString* tname;
    while ( ( tname = static_cast<TObjString*>(next()) ) )
    {
      fTriggerClasses->Add(new TObjString(*tname));
      if ( tname->String().BeginsWith("CMB") )
      {
        fAA = kTRUE;
      }
    }
  }  
  
  DefineCentralityClasses();  
}

//_____________________________________________________________________________
AliAnalysisMuMu::~AliAnalysisMuMu()
{
  delete fTriggerClasses;
  delete fSingleTrackCutNames;
  delete fPairTrackCutNames;
  if(fOutput && ! AliAnalysisManager::GetAnalysisManager()->IsProofMode()) 
  {
    delete fOutput;     
  }
  delete fCentralityNames;
  delete fGlobalEventSelectionNames;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::AddGlobalEventSelection(const char* name)
{
  if (!fGlobalEventSelectionNames)
  {
    fGlobalEventSelectionNames = new TObjArray;
    fGlobalEventSelectionNames->SetOwner(kTRUE);
  }
  fGlobalEventSelectionNames->Add(new TObjString(name));
}
                                  
//_____________________________________________________________________________
UInt_t AliAnalysisMuMu::CodePairCutMask(UInt_t maskForOneOrBothTrack, UInt_t maskForTrackPair) const
{
  return ( ( maskForOneOrBothTrack << 16 ) | maskForTrackPair );
}

//_____________________________________________________________________________
void AliAnalysisMuMu::DecodePairCutMask(UInt_t code, UInt_t& maskForOneOrBothTrack, UInt_t& maskForTrackPair) const
{
  maskForOneOrBothTrack = ( ( code & 0xFFFF0000 ) >> 16 );
  maskForTrackPair = ( code & 0xFFFF );
}

//_____________________________________________________________________________
void AliAnalysisMuMu::DefineCentralityClasses()
{

  if ( fAA ) 
  {
    fCentralityLimits.push_back(10.0);
    fCentralityLimits.push_back(20.0);
    fCentralityLimits.push_back(30.0);
    fCentralityLimits.push_back(40.0);
    fCentralityLimits.push_back(50.0);
    fCentralityLimits.push_back(60.0);
    fCentralityLimits.push_back(70.0);
    fCentralityLimits.push_back(80.0);
    fCentralityLimits.push_back(90.0);
  }
//  fCentralityLimits.push_back(100.0);
  
  for ( std::vector<double>::size_type i = 0; i < fCentralityLimits.size(); ++i )
  {
    Double_t limit = fCentralityLimits[i];
    fCentralityNames->Add(new TObjString(CentralityName(limit)));
  }
  
  fCentralityNames->Add(new TObjString(DefaultCentralityName()));
  fCentralityNames->SetOwner(kTRUE);
}

//_____________________________________________________________________________
void AliAnalysisMuMu::DefineDefaultTriggerClasses()
{
  fTriggerClasses->Add(new TObjString("ANY"));
  
  fTriggerClasses->Add(new TObjString("CBEAMB-ABCE-NOPF-ALL"));
  fTriggerClasses->Add(new TObjString("CINT1B-ABCE-NOPF-ALL"));
  fTriggerClasses->Add(new TObjString("CMUS1B-ABCE-NOPF-MUON"));
  fTriggerClasses->Add(new TObjString("CINT1C-ABCE-NOPF-ALL"));
  fTriggerClasses->Add(new TObjString("CINT1A-ABCE-NOPF-ALL"));
  fTriggerClasses->Add(new TObjString("CMUS1C-ABCE-NOPF-MUON"));
  fTriggerClasses->Add(new TObjString("CMUS1A-ABCE-NOPF-MUON"));
  
  fTriggerClasses->Add(new TObjString("CINT1-B-NOPF-ALLNOTRD"));
  fTriggerClasses->Add(new TObjString("CINT1-AC-NOPF-ALLNOTRD"));
  fTriggerClasses->Add(new TObjString("CMUS1-B-NOPF-ALLNOTRD"));    
  fTriggerClasses->Add(new TObjString("CMUS1-AC-NOPF-ALLNOTRD"));
  
  fTriggerClasses->Add(new TObjString("CINT5-B-NOPF-ALLNOTRD"));
  fTriggerClasses->Add(new TObjString("CINT5-AC-NOPF-ALLNOTRD"));
  fTriggerClasses->Add(new TObjString("CMUS5-B-NOPF-ALLNOTRD"));    
  fTriggerClasses->Add(new TObjString("CMUS5-AC-NOPF-ALLNOTRD"));
  
  fTriggerClasses->Add(new TObjString("CINT7-B-NOPF-ALLNOTRD"));
  fTriggerClasses->Add(new TObjString("CMUSH7-B-NOPF-MUON"));	
  fTriggerClasses->Add(new TObjString("CMUL7-B-NOPF-MUON"));
  fTriggerClasses->Add(new TObjString("CMUU7-B-NOPF-MUON"));
}
                                  
//_____________________________________________________________________________
void AliAnalysisMuMu::AddPairCut(const char* cutName, UInt_t maskForOneOrBothTrack, UInt_t maskForTrackPair)
{
  /// Add a cut for a pair.
  /// maskForOneOrBothTrack is the mask of cuts that at least one track must satisfy
  /// maskForTrackPair is the mask of cuts that *both* tracks must satisfy.
  /// if maskForTrackPair is 0, then no (extra) condition is applied to the pair
  
  if ( !fPairTrackCutNames ) 
  {
    fPairTrackCutNames = new TObjArray;
    fPairTrackCutNames->SetOwner(kTRUE);
  }
  TObjString* oname = new TObjString(Form("p%s",cutName));
  if ( !maskForTrackPair ) maskForTrackPair = maskForOneOrBothTrack;
  oname->SetUniqueID(CodePairCutMask(maskForOneOrBothTrack,maskForTrackPair));
  fPairTrackCutNames->Add(oname);
}

//_____________________________________________________________________________
void AliAnalysisMuMu::AddSingleCut(const char* name, UInt_t mask)
{
  if ( !fSingleTrackCutNames ) 
  {
    fSingleTrackCutNames = new TObjArray;
    fSingleTrackCutNames->SetOwner(kTRUE);
  }
  TObjString* oname = new TObjString(Form("s%s",name));
  oname->SetUniqueID(mask);
  fSingleTrackCutNames->Add(oname);
}

//_____________________________________________________________________________
void AliAnalysisMuMu::BeautifyHistos()
{
  /// Put the titles, marker sizes, color, etc...
  
  TIter next(fHistogramCollection->CreateIterator());
  TH1* h;
  
  while ( ( h = static_cast<TH1*>(next()) ) )
  {
    TString name(h->GetName());
    
    if ( name.Contains("Plus") ) 
    {
      h->SetMarkerStyle(kFullCircle);
      h->SetMarkerColor(kBlue);
      h->SetLineColor(kBlue);        
    }
    else
    {
      h->SetMarkerStyle(kOpenCircle);
      h->SetMarkerColor(kRed);
      h->SetLineColor(kRed);        
    }      
  }
}

//_____________________________________________________________________________
void 
AliAnalysisMuMu::CreateSingleHisto(const char* physics,
                                   const char* triggerClassName,
                                   const char* hname, const char* htitle, 
                                   Int_t nbinsx, Double_t xmin, Double_t xmax,
                                   Int_t nbinsy, Double_t ymin, Double_t ymax) const
{  
  CreateHisto(fSingleTrackCutNames,physics,triggerClassName,hname,htitle,
              nbinsx,xmin,xmax,nbinsy,ymin,ymax);
}

//_____________________________________________________________________________
void 
AliAnalysisMuMu::CreatePairHisto(const char* physics,
                                 const char* triggerClassName,
                                 const char* hname, const char* htitle, 
                                 Int_t nbinsx, Double_t xmin, Double_t xmax,
                                 Int_t nbinsy, Double_t ymin, Double_t ymax) const
{
  CreateHisto(fPairTrackCutNames,physics,triggerClassName,hname,htitle,
              nbinsx,xmin,xmax,nbinsy,ymin,ymax);
}

//_____________________________________________________________________________
void 
AliAnalysisMuMu::CreateEventHisto(const char* physics,
                                  const char* triggerClassName,
                                  const char* hname, const char* htitle, 
                                  Int_t nbinsx, Double_t xmin, Double_t xmax,
                                  Int_t nbinsy, Double_t ymin, Double_t ymax) const
{
  AliDebug(1,Form("physics=%s trigger=%s hname=%s",physics,triggerClassName,hname));
  
  TIter next(fCentralityNames);
  TObjString* cent;
  
  while ( ( cent = static_cast<TObjString*>(next()) ) )
  {
    TH1* h(0x0);
    
    if ( nbinsy > 0 )
    {  
      h = new TH2F(hname,htitle,nbinsx,xmin,xmax,nbinsy,ymin,ymax);
    }
    else
    {
      h = new TH1F(hname,htitle,nbinsx,xmin,xmax);
    }
    
    fHistogramCollection->Adopt(physics,triggerClassName,cent->String().Data(),h); 
  }
}

//_____________________________________________________________________________
void 
AliAnalysisMuMu::CreateHisto(TObjArray* array,
                             const char* physics,
                             const char* triggerClassName,
                             const char* hname, const char* htitle, 
                             Int_t nbinsx, Double_t xmin, Double_t xmax,
                             Int_t nbinsy, Double_t ymin, Double_t ymax) const
{
  AliDebug(1,Form("physics=%s trigger=%s hname=%s",physics,triggerClassName,hname));

  TIter next(array);
  TObjString* tcn;
  while ( ( tcn = static_cast<TObjString*>(next()) ) )
  {
    TIter nextCent(fCentralityNames);
    TObjString* cent;
    
    while ( ( cent = static_cast<TObjString*>(nextCent()) ) )
    {
      TH1* h(0x0);
    
      if ( nbinsy > 0 )
      {  
        h = new TH2F(hname,htitle,nbinsx,xmin,xmax,nbinsy,ymin,ymax);
      }
      else
      {
        h = new TH1F(hname,htitle,nbinsx,xmin,xmax);
      }
    
      fHistogramCollection->Adopt(physics,triggerClassName,cent->String().Data(),tcn->String().Data(),h);
    }
  }      
}

//_____________________________________________________________________________
const char* 
AliAnalysisMuMu::DefaultCentralityName() const
{
  if ( fAA ) return "CENTX";
  else return "PP";
}

//_____________________________________________________________________________
const char* 
AliAnalysisMuMu::CentralityName(Double_t centrality) const
{
  if ( centrality > 0 && centrality <= 100.0 ) 
  {
    return Form("CENT%02d",TMath::Nint(centrality));
  }
  else
  {
    return DefaultCentralityName();
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMu::AssertHistogramCollection(const char* physics, const char* triggerClassName)
{
  // insure that a given set of histogram is created
  TH1* test = fHistogramCollection->Histo(physics,triggerClassName,DefaultCentralityName(),"Zvertex");
  if (!test) FillHistogramCollection(physics,triggerClassName);
}

//_____________________________________________________________________________
void AliAnalysisMuMu::FillHistogramCollection(const char* physics, const char* triggerClassName)
{
  AliDebug(1,Form("(%s,%s)",physics,triggerClassName));
  
  Double_t ptMin = 0;
  Double_t ptMax = 50;
  Int_t nbinsPt = GetNbins(ptMin,ptMax,0.25);
  Double_t pMin = 0;
  Double_t pMax = 300;
  Int_t nbinsP = GetNbins(pMin,pMax,1.0);
  Double_t etaMin = -6;
  Double_t etaMax = 1;
  Int_t nbinsEta = GetNbins(etaMin,etaMax,0.1);
//  Double_t yMin = -10;
//  Double_t yMax = 10;
//  Int_t nbinsY = GetNbins(yMin,yMax,0.1);
  
  CreateSingleHisto(physics,triggerClassName,"PtEtaMuPlus", "#mu+ P_{T} distribution vs Eta", nbinsEta,etaMin,etaMax, nbinsPt,ptMin,ptMax);
  CreateSingleHisto(physics,triggerClassName,"PtEtaMuMinus", "#mu- P_{T} distribution vs Eta",nbinsEta,etaMin,etaMax, nbinsPt,ptMin,ptMax);
  
  CreateSingleHisto(physics,triggerClassName,"Pt", "P_{T} distribution",nbinsPt,ptMin,ptMax);
  
  CreateSingleHisto(physics,triggerClassName,"PEtaMuPlus", "#mu+ P distribution",nbinsEta,etaMin,etaMax,nbinsP,pMin,pMax);
  CreateSingleHisto(physics,triggerClassName,"PEtaMuMinus", "#mu- P distribution",nbinsEta,etaMin,etaMax,nbinsP,pMin,pMax);  
  
  Double_t chi2min = 0;
  Double_t chi2max = 20;
  Int_t nbinchi2 = GetNbins(chi2min,chi2max,0.1);
  
  CreateSingleHisto(physics, triggerClassName, "Chi2MuPlus", "chisquare per NDF #mu+", nbinchi2, chi2min, chi2max);
  CreateSingleHisto(physics, triggerClassName, "Chi2MuMinus", "chisquare per NDF #mu-", nbinchi2, chi2min, chi2max);
    
  CreateSingleHisto(physics, triggerClassName, "Chi2", "chisquare per NDF", nbinchi2, chi2min, chi2max);
  CreateSingleHisto(physics, triggerClassName, "Chi2Pt2GeV", "chisquare per NDF (p_{T} > 2 GeV/c)", nbinchi2, chi2min, chi2max);
  
  Double_t nChHitMin = 0;
  Double_t nChHitMax = 20;
  Int_t nbinNChHit = GetNbins(nChHitMin,nChHitMax,1.);
  
  CreateSingleHisto(physics, triggerClassName, "NChHit", "Number of chamber hit", nbinNChHit, nChHitMin, nChHitMax);
  CreateSingleHisto(physics, triggerClassName, "NChHitPt2GeV", "Number of chamber hit (p_{T} > 2 GeV/c)", nbinNChHit, nChHitMin, nChHitMax);
  
  Double_t MinvMin = 0;
  Double_t MinvMax = 16;
  Int_t nMinvBins = GetNbins(MinvMin,MinvMax,0.05);
  
  CreatePairHisto(physics,triggerClassName,"MinvUSPt", "#mu+#mu- inv. mass vs Pt",nbinsPt,ptMin,ptMax,nMinvBins,MinvMin,MinvMax);
  CreatePairHisto(physics,triggerClassName,"MinvPPPt", "#mu+#mu+ inv. mass vs Pt",nbinsPt,ptMin,ptMax,nMinvBins,MinvMin,MinvMax);
  CreatePairHisto(physics,triggerClassName,"MinvMMPt", "#mu-#mu- inv. mass vs Pt",nbinsPt,ptMin,ptMax,nMinvBins,MinvMin,MinvMax);

  Double_t xmin = -40;
  Double_t xmax = +40;
  Int_t nbins = GetNbins(xmin,xmax,0.5);
  
  CreateEventHisto(physics,triggerClassName,"Zvertex","z vertex",nbins,xmin,xmax);  
  CreateEventHisto(physics,triggerClassName,"Xvertex","x vertex",nbins,xmin,xmax);  
  CreateEventHisto(physics,triggerClassName,"Yvertex","y vertex",nbins,xmin,xmax);  

  CreateEventHisto(physics,triggerClassName,"Nevents","number of events",2,-0.5,1.5);  

  xmin = -200;
  xmax = +200;
  nbins = GetNbins(xmin,xmax,1.0);
  
  CreateSingleHisto(physics,triggerClassName,"XYdcaMuPlus","#mu+ DCA non bending vs bending;dcaY;dcaX",nbins,xmin,xmax,nbins,xmin,xmax);  
  CreateSingleHisto(physics,triggerClassName,"XYdcaMuMinus","#mu- DCA non bending vs bending;dcaY;dcaX",nbins,xmin,xmax,nbins,xmin,xmax);  
  
  CreateSingleHisto(physics,triggerClassName,"dcaP23MuPlus","#mu+ DCA vs P for 2-3 degrees;P (GeV);DCA (cm)",nbinsP,pMin,pMax,nbins,xmin,xmax);  
  CreateSingleHisto(physics,triggerClassName,"dcaP23MuMinus","#mu- DCA vs P for 2-3 degrees;P (GeV);DCA (cm)",nbinsP,pMin,pMax,nbins,xmin,xmax);  
  CreateSingleHisto(physics,triggerClassName,"dcaP310MuPlus","#mu+ DCA vs P for 3-10 degrees;P (GeV);DCA (cm)",nbinsP,pMin,pMax,nbins,xmin,xmax);  
  CreateSingleHisto(physics,triggerClassName,"dcaP310MuMinus","#mu- DCA vs P for 3-10 degrees;P (GeV);DCA (cm)",nbinsP,pMin,pMax,nbins,xmin,xmax);  

  CreateSingleHisto(physics,triggerClassName,"dcaPwPtCut23MuPlus","#mu+ DCA vs P for 2-3 degrees with Pt Cut;P (GeV);DCA (cm)",nbinsP,pMin,pMax,nbins,xmin,xmax);  
  CreateSingleHisto(physics,triggerClassName,"dcaPwPtCut23MuMinus","#mu- DCA vs P for 2-3 degrees with Pt Cut;P (GeV);DCA (cm)",nbinsP,pMin,pMax,nbins,xmin,xmax);  
  CreateSingleHisto(physics,triggerClassName,"dcaPwPtCut310MuPlus","#mu+ DCA vs P for 3-10 degrees with Pt Cut;P (GeV);DCA (cm)",nbinsP,pMin,pMax,nbins,xmin,xmax);  
  CreateSingleHisto(physics,triggerClassName,"dcaPwPtCut310MuMinus","#mu- DCA vs P for 3-10 degrees with Pt Cut;P (GeV);DCA (cm)",nbinsP,pMin,pMax,nbins,xmin,xmax);  
  
  xmin = 0;
  xmax = 10000;
  nbins = GetNbins(xmin,xmax,10.0);
  
  CreateSingleHisto(physics,triggerClassName,"pdca23","P#timesDCA for 2-3 degrees;PDCA (GeVcm)",nbins,xmin,xmax);  
  CreateSingleHisto(physics,triggerClassName,"pdca23Pt2GeV","P#timesDCA for 2-3 degrees (p_{T} > 2 GeV/c);PDCA (GeVcm)",nbins,xmin,xmax);  
  CreateSingleHisto(physics,triggerClassName,"pdca310","P#timesDCA for 3-10 degrees;PDCA (GeVcm)",nbins,xmin,xmax);  
  CreateSingleHisto(physics,triggerClassName,"pdca310Pt2GeV","P#timesDCA for 3-10 degrees (p_{T} > 2 GeV/c);PDCA (GeVcm)",nbins,xmin,xmax);  
  
  //  xmin = 0;
//  xmax = 2000; 
//  nbins = GetNbins(xmin,xmax,5);
//  
//  CreateSingleHisto(physics,triggerClassName,"PDCAP23MuMinus","#mu- P.DCA vs P for 2-3 degrees;P (GeV);P.DCA (GeV.cm)",nbinsP,pMin,pMax,nbins,xmin,xmax);
//  CreateSingleHisto(physics,triggerClassName,"PDCAP23MuPlus","#mu+ P.DCA vs P for 2-3 degrees;P (GeV);P.DCA (GeV.cm)",nbinsP,pMin,pMax,nbins,xmin,xmax);
//  
//  CreateSingleHisto(physics,triggerClassName,"PDCAP310MuMinus","#mu- P.DCA vs P for 3-10 degrees;P (GeV);P.DCA (GeV.cm)",nbinsP,pMin,pMax,nbins,xmin,xmax);
//  CreateSingleHisto(physics,triggerClassName,"PDCAP310MuPlus","#mu+ P.DCA vs P for 3-10 degrees;P (GeV);P.DCA (GeV.cm)",nbinsP,pMin,pMax,nbins,xmin,xmax);

  if ( fIsFromESD ) 
  {
    xmin = -80;
    xmax = 80;
    nbins = GetNbins(xmin,xmax,0.1);
    CreateEventHisto(physics,triggerClassName,"V0Time","Mean Time V0C versus V0A;Time V0A (ns);Time V0C (ns)",nbins,xmin,xmax,nbins,xmin,xmax);
    
    xmin = 0;
    xmax = 20000;
    nbins = GetNbins(xmin,xmax,10);
    
    CreateEventHisto(physics,triggerClassName,"V0Amplitude","V0Amult+V0Cmult",nbins,xmin,xmax);
  }
  
  Double_t* vbins = new Double_t[fCentralityLimits.size()+1];
  
  vbins[0] = 0.0;
  
  for ( std::vector<double>::size_type i = 0; i < fCentralityLimits.size(); ++i ) 
  {  
    vbins[i+1] = fCentralityLimits[i];
  }
  
  TH1* h = new TH1F("Centrality","Centrality",fCentralityLimits.size(),vbins);
  
  delete[] vbins;
  
  fHistogramCollection->Adopt(physics,triggerClassName,h);                                                              

  xmin = 0;
  xmax = 5000;
  nbins = GetNbins(xmin,xmax,10);
  
  CreateEventHisto(physics,triggerClassName,"Tracklets","Number of tracklets",nbins,xmin,xmax);
  
}


//_____________________________________________________________________________
Double_t AliAnalysisMuMu::MuonMass2() const
{
  static Double_t m2 = 1.11636129640000012e-02; // using a constant here as the line below is a problem for CINT...
//  static Double_t m2 = TDatabasePDG::Instance()->GetParticle("mu-")->Mass()*TDatabasePDG::Instance()->GetParticle("mu-")->Mass();
  return m2;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::FinishTaskOutput()
{
  /// prune empty histograms BEFORE mergin, in order to save some bytes...
  
  if ( fHistogramCollection )
  {
    UInt_t size = fHistogramCollection->EstimateSize(kFALSE);
    
    AliInfo(Form("size before prune histograms = %5.1f MB",size/1024.0/1024.0));
    
    fHistogramCollection->PruneEmptyHistograms();  
  }
}

//_____________________________________________________________________________
TH1* AliAnalysisMuMu::Histo(const char* physics, const char* triggerClassName, const char* histoname)
{
  return fHistogramCollection ? fHistogramCollection->Histo(physics,triggerClassName,histoname) : 0x0;
}

//_____________________________________________________________________________
TH1* AliAnalysisMuMu::Histo(const char* physics, const char* histoname)
{
  return fHistogramCollection ? fHistogramCollection->Histo(physics,histoname) : 0x0;
}

//_____________________________________________________________________________
TH1* AliAnalysisMuMu::Histo(const char* physics,
                            const char* triggerClassName, 
                            const char* what,
                            const char* histoname)
{
  return fHistogramCollection ? fHistogramCollection->Histo(physics,triggerClassName,what,histoname) : 0x0;
}

//_____________________________________________________________________________
TH1* AliAnalysisMuMu::Histo(const char* physics,
                            const char* triggerClassName, 
                            const char* cent,
                            const char* what,
                            const char* histoname)
{
  return fHistogramCollection ? fHistogramCollection->Histo(physics,triggerClassName,cent,what,histoname) : 0x0;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::NotifyRun()
{
  AliInfo(Form("Run %09d File %s",fCurrentRunNumber,CurrentFileName()));
}

//_____________________________________________________________________________
void AliAnalysisMuMu::Plot(const char* physics, const char* triggerClassName)
{
  if ( !fHistogramCollection->Histo(physics,triggerClassName,"Zvertex") ) return;
  
  TH1* hzvertex = static_cast<TH1*>(fHistogramCollection->Histo(physics,triggerClassName,"Zvertex")->Clone());
  
  TCanvas* c1 = new
  TCanvas(Form("%s-%s",ClassName(),triggerClassName),
          Form("%s-%s",ClassName(),triggerClassName),
          10,10,510,510);
  
  TH1* hptmuplus(0x0);
  TH1* hptmuminus(0x0);
  
  if ( fHistogramCollection->Histo(physics,triggerClassName,"sALL","PtEtaMuPlus") )
  {
    hptmuplus = (static_cast<TH2*>(fHistogramCollection->Histo(physics,triggerClassName,"sALL","PtEtaMuPlus")))->ProjectionY();
    hptmuplus->SetDirectory(0);
  }
  
  if ( fHistogramCollection->Histo(physics,triggerClassName,"sALL","PtEtaMuMinus") )
  {
    hptmuminus = (static_cast<TH2*>(fHistogramCollection->Histo(physics,triggerClassName,"sALL","PtEtaMuMinus")))->ProjectionY();
    hptmuminus->SetDirectory(0);
  }
  
  c1->Divide(1,2);
  c1->cd(1)->SetLogy();
  
  if (hptmuplus)
  {
    hptmuplus->SetLineColor(2);
    hptmuplus->Draw("HISTE");
  }
  if (hptmuminus)
  {
    hptmuminus->SetLineColor(4);
    if ( hptmuplus )
    {
      hptmuminus->Draw("HISTESAME");      
    }
    else
    {
      hptmuminus->Draw("HISTE");
    }
  }
  
  TPaveText* text = new TPaveText(0.20,0.90,0.6,0.99,"NDC");
  
  Int_t nplus = hptmuplus ? TMath::Nint(hptmuplus->GetEntries()) : 0;
  Int_t nminus = hptmuminus ? TMath::Nint(hptmuminus->GetEntries()) : 0;
  
  text->AddText(Form("N #mu %d N #mu+ %d N #mu- %d",
                     nplus+nminus,nplus,nminus));
  
  text->Draw();
  
  c1->cd(2);
  hzvertex->Draw("HISTE");
}

//_____________________________________________________________________________
void 
AliAnalysisMuMu::Print(Option_t* /*opt*/) const
{
  TIter next(fSingleTrackCutNames);
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    cout << Form("SINGLE CUT %20s MASK %x",str->String().Data(),str->GetUniqueID()) << endl;
  }
  
  TIter next2(fPairTrackCutNames);
  while ( ( str = static_cast<TObjString*>(next2()) ) )
  {
    UInt_t single(0);
    UInt_t pair(0);
    DecodePairCutMask(str->GetUniqueID(),single,pair);
    
    cout << Form("PAIR CUT %20s UID %x SINGLE MASK %x PAIR MASK %x",str->String().Data(),str->GetUniqueID(),single,pair) << endl;
  }
  
}

//_____________________________________________________________________________
void
AliAnalysisMuMu::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  
  fOutput = dynamic_cast<TList*>(GetOutputData(1));
  
  if (!fOutput) return;
  
  fHistogramCollection = dynamic_cast<AliHistogramCollection*>(fOutput->FindObject("mumu"));
  
  if (!fHistogramCollection)
  {
    AliError("Could not find back histogram collection in output...");
    return;
  }
  
  fHistogramCollection->PruneEmptyHistograms();

  UInt_t size2 = fHistogramCollection->EstimateSize();

  AliInfo(Form("size after prune histograms = %5.1f MB",size2/1024.0/1024.0));
  
  BeautifyHistos();

  fHistogramCollection->Print("*Minv*");
  
  fEventCounters = dynamic_cast<AliCounterCollection*>(fOutput->FindObject("eventCounters"));
  
  if (!fEventCounters)
  {
    AliError("Could not find back counters in output...");
    return;
  }
  
  fEventCounters->Print();
}


//_____________________________________________________________________________
void AliAnalysisMuMu::UserExec(Option_t* opt)
{
  MuUserExec(opt);
  
  // Post output data.
  PostData(1, fOutput);      
}

//_____________________________________________________________________________
void AliAnalysisMuMu::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  
  fTriggerClasses->Print();  
  
  fOutput = new TList;
  fOutput->SetOwner(kTRUE);
  
  fHistogramCollection = new AliHistogramCollection("mumu");
  
  fOutput->Add(fHistogramCollection);  
  
  TString triggerKeys;
  TIter next(fTriggerClasses);
  TObjString* tname;
  
  while ( ( tname = static_cast<TObjString*>(next()) ) )
  {
    triggerKeys += tname->GetName();
    triggerKeys += "/";
  }
  
  triggerKeys += "other/any";
  
  // initialize event counters
  fEventCounters = new AliCounterCollection("eventCounters");
  TIter nextGSN(fGlobalEventSelectionNames);
  TObjString* str;
  TString eventRubric;
  while ( ( str = static_cast<TObjString*>(nextGSN()) ) )
  {
    if ( eventRubric.Length() > 0 ) eventRubric += "/"; 
    eventRubric += str->String();
  }
  
  fEventCounters->AddRubric("event", eventRubric.Data());
  
  fEventCounters->AddRubric("trigger", triggerKeys.Data());
  fEventCounters->AddRubric("run", 1000000);
  fEventCounters->Init();
  
  fOutput->Add(fEventCounters);  
  
  // Post output data.
  PostData(1, fOutput);      
}

