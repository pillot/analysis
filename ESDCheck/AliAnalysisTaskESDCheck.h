#ifndef AliAnalysisTaskESDCheck_h
#define AliAnalysisTaskESDCheck_h

//
// scan ESD tracks and produce histograms
// Author: Philippe Pillot
//

#include <TMatrixD.h>

class TH1F;
class TObjArray;
class AliESDEvent;
class AliMUONTrack;
class AliCounterCollection;

class AliAnalysisTaskESDCheck : public AliAnalysisTaskSE {
 public:
  
  enum StationSelect {k1245, k245, kAll};
  
  AliAnalysisTaskESDCheck(const char *name = "ESDScan");
  virtual ~AliAnalysisTaskESDCheck();
  
//  virtual void   LocalInit();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   NotifyRun();
  virtual void   Terminate(Option_t *);
  
  void SelectStation(StationSelect stSelect) {fStSelect = stSelect;}
  void SetFullPrintout(Bool_t flag = kTRUE) {fFullPrintout = flag;}
  void SelectCharge(Short_t charge = 0) {fSelectCharge = charge;}

private:
  
  /// Not implemented
  AliAnalysisTaskESDCheck(const AliAnalysisTaskESDCheck& rhs);
  /// Not implemented
  AliAnalysisTaskESDCheck& operator = (const AliAnalysisTaskESDCheck& rhs);
  
  Double_t ChangeThetaRange(Double_t theta);
  
  Bool_t SelectTrack(AliMUONTrack &track);
  
  void Cov2CovP(const TMatrixD &param, TMatrixD &cov);
  
  void PrintCurrentStat();
  void PrintTotalStat();
  
private:
  
  enum EESD { 
    kESDnTracks                 = 0,  ///< number of tracks
    kESDMatchTrig               = 1,  ///< number of tracks matched with trigger
    kESDMomentum                = 2,  ///< P distribution
    kESDPt                      = 3,  ///< Pt distribution
    kESDRapidity                = 4,  ///< rapidity distribution
    kESDChi2                    = 5,  ///< normalized chi2 distribution
    kESDProbChi2                = 6,  ///< distribution of probability of chi2
    
    kESDClusterHitMap           = 7,  ///< cluster position distribution in chamber i
    kESDnClustersPerTrack       = 17, ///< number of clusters per track
    kESDnClustersPerCh          = 18, ///< number of clusters per chamber per track
    kESDnClustersPerDE          = 19, ///< number of clusters per DE per track
    kESDClusterChargeInCh       = 20, ///< cluster charge distribution in chamber i
    kESDClusterChargePerChMean  = 30, ///< cluster charge per Ch: mean
    kESDClusterChargePerChSigma = 31, ///< cluster charge per Ch: dispersion
    kESDClusterChargePerDE      = 32, ///< cluster charge per DE: mean
    kESDClusterSizeInCh         = 33, ///< cluster size distribution in chamber i
    kESDClusterSizePerChMean    = 43, ///< cluster size per Ch: mean
    kESDClusterSizePerChSigma   = 44, ///< cluster size per Ch: dispersion
    kESDClusterSizePerDE        = 45, ///< cluster size per DE: mean
    
    kESDResidualXInCh           = 46, ///< cluster-track residual-X distribution in chamber i
    kESDResidualYInCh           = 56, ///< cluster-track residual-Y distribution in chamber i
    kESDResidualXPerChMean      = 66, ///< cluster-track residual-X per Ch: mean
    kESDResidualYPerChMean      = 67, ///< cluster-track residual-Y per Ch: mean
    kESDResidualXPerChSigma     = 68, ///< cluster-track residual-X per Ch: dispersion
    kESDResidualYPerChSigma     = 69, ///< cluster-track residual-Y per Ch: dispersion
    kESDResidualXPerDEMean      = 70, ///< cluster-track residual-X per DE: mean
    kESDResidualYPerDEMean      = 71, ///< cluster-track residual-Y per DE: mean
    kESDResidualXPerDESigma     = 72, ///< cluster-track residual-X per DE: dispersion
    kESDResidualYPerDESigma     = 73, ///< cluster-track residual-Y per DE: dispersion
    kESDLocalChi2XInCh          = 74, ///< local chi2-X distribution in chamber i
    kESDLocalChi2YInCh          = 84, ///< local chi2-Y distribution in chamber i
    kESDLocalChi2XPerChMean     = 94, ///< local chi2-X per Ch: mean
    kESDLocalChi2YPerChMean     = 95, ///< local chi2-Y per Ch: mean
    kESDLocalChi2XPerDEMean     = 96, ///< local chi2-X per DE: mean
    kESDLocalChi2YPerDEMean     = 97, ///< local chi2-Y per DE: mean
    kESDLocalChi2InCh           = 98, ///< local chi2-X distribution in chamber i
    kESDLocalChi2PerChMean      = 108, ///< local chi2 per Ch: mean
    kESDLocalChi2PerDEMean      = 109, ///< local chi2 per DE: mean
    
    kESDThetaX                  = 110, ///< thetaX distribution
    kESDThetaY                  = 111, ///< thetaY distribution
    
    kESDMomentumUncorrected     = 112, ///< P distribution uncorrected
    kESDMomentumUncorrectedG    = 113, ///< P graph uncorrected
    kESDMomentumErrorG          = 114, ///< P errors
    kESDMomentumRelativeErrorG  = 115, ///< P relative errors
    kESDMomentumRelativeErrorMinG=116, ///< P relative errors - low edge
    kESDMomentumRelativeErrorMaxG=117, ///< P relative errors - high edge
    kESDMomentumRecoError       = 118, ///< P errors from track parameter covariance matrix versus p (at vertex)
    kESDUncorrMomentumRecoError = 119, ///< P errors from track parameter covariance matrix versus p (uncorrected)
    
    kESDSign                    = 120, ///< track sign
    kESDDCA                     = 121, ///< DCA distribution
    kESDnChamberHitPerTrack     = 122, ///< number of chamber hit per track
    kESDPtRecoError             = 123, ///< Pt errors from track parameter covariance matrix versus pt (at vertex)
    kESDUncorrPtRecoError       = 124, ///< Pt errors from track parameter covariance matrix versus pt (uncorrected)
    kESDSlopeXRecoError         = 125, ///< slopeX errors from track parameter covariance matrix versus p (at vertex)
    kESDUncorrSlopeXRecoError   = 126, ///< slopeX errors from track parameter covariance matrix versus p (uncorrected)
    kESDSlopeYRecoError         = 127, ///< slopeY errors from track parameter covariance matrix versus p (at vertex)
    kESDUncorrSlopeYRecoError   = 128, ///< slopeY errors from track parameter covariance matrix versus p (uncorrected)
    kESDPDCARecoError           = 129, ///< p*DCA errors from track parameter covariance matrix versus p (at DCA)
    
    kESDnTotClustersPerCh       = 1000, ///< total number of associated clusters per chamber
    kESDnTotClustersPerDE       = 1001, ///< total number of associated clusters per DE
    kESDnTotFullClustersPerDE   = 1002, ///< total number of associated clusters containing pad info per DE
    kESDSumClusterChargePerDE   = 1003, ///< sum of cluster charge per DE
    kESDSumClusterSizePerDE     = 1004, ///< sum of cluster size per DE
    kESDSumResidualXPerDE       = 1005, ///< sum of cluster-track residual-X per DE
    kESDSumResidualYPerDE       = 1006, ///< sum of cluster-track residual-Y per DE
    kESDSumResidualX2PerDE      = 1007, ///< sum of cluster-track residual-X**2 per DE
    kESDSumResidualY2PerDE      = 1008, ///< sum of cluster-track residual-Y**2 per DE
    kESDSumLocalChi2XPerDE      = 1009, ///< sum of local chi2-X per DE
    kESDSumLocalChi2YPerDE      = 1010, ///< sum of local chi2-Y per DE
    kESDSumLocalChi2PerDE       = 1011  ///< sum of local chi2 per DE
  };
  
  TObjArray*  fList;       //!< List of output object for everybody
  TObjArray*  fListExpert; //!< List of output object for experts

  AliCounterCollection* fTrackCounters; //!< track statistics
  AliCounterCollection* fEventCounters; //!< event statistics
  
  StationSelect fStSelect; ///< Select the requested station to validate a track
  Bool_t fFullPrintout;    ///< Switch ON/OFF the printout of statics at each event
  Short_t fSelectCharge;   ///< Fill histograms only with negative/position tracks (0=all)
  Bool_t fOCDBLoaded;      //!< flag telling if the OCDB has been properly loaded or not
  
  static const Int_t fgkNTriggerClass;
  static const char* fgkTriggerClass[10];
  static const char* fgkTriggerShortName[11];
  
  Int_t fcurrentNEvnt;        //!< number of events with something in MUON in the current run
  Int_t fTotalNEvnt;          //!< total number of events with something in MUON
  Int_t fcurrentNTotEvnt;     //!< number of events in the current run
  Int_t fTotalNTotEvnt;       //!< total number of events
  Int_t fcurrentNTracks[3];   //!< tracker = [0], trigger = [1], matched = [2]
  Int_t fTotalNTracks[3];     //!< tracker = [0], trigger = [1], matched = [2]
  Int_t fCurrentNTrig[11];    //!< number if events in [trigger class] with something in MUON in the current run
  Int_t fTotalNTrig[11];      //!< number if events in [trigger class] in the current run
  Int_t fCurrentNTotTrig[11]; //!< total number if events in [trigger class] with something in MUON
  Int_t fTotalNTotTrig[11];   //!< total number if events in [trigger class]
  Int_t fCurrentStat[11][3];  //!< statistics of the current run [trigger class][track type]
  Int_t fTotalStat[11][3];    //!< statistics of all runs [trigger class][track type]
  
  Int_t fPreviousRun; //!< previous run number
  
//  TH1F* fTest;
  
  ClassDef(AliAnalysisTaskESDCheck, 1); // Single muon analysis
};

#endif

