#include <cmath>
#include <stdexcept>
#include <vector>

#include <gsl/span>

#include <TCanvas.h>
#include <TH1F.h>

#include "DetectorsCommonDataFormats/DetID.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/TrackReference.h"
#include "Steer/MCKinematicsReader.h"

using namespace std;
using namespace o2;
using namespace o2::steer;

constexpr double pi() { return 3.14159265358979323846; }
void CreateResiduals(std::vector<TH1*>& histos);
void FillResiduals(const TrackReference& mcTrackRef, const MCTrack& mcTrack, std::vector<TH1*>& histos);
void DrawResiduals(std::vector<TH1*>& histos);

//_________________________________________________________________________________________________
void CheckMC()
{
  std::vector<TH1*> residuals{};
  CreateResiduals(residuals);

  MCKinematicsReader mcReader{};
  if (!mcReader.initFromDigitContext("collisioncontext.root")) {
    throw invalid_argument("initialization of MCKinematicsReader failed");
  }

  for (int iSource = 0; iSource < (int) mcReader.getNSources(); ++iSource) {
    for (int iEvent = 0; iEvent < (int)mcReader.getNEvents(iSource); ++iEvent) {
      const auto& mcTracks = mcReader.getTracks(iSource, iEvent);
      for (int iTrack = 0; iTrack < (int)mcTracks.size(); ++iTrack) {
        const auto& mcTrack = mcTracks[iTrack];
        if (mcTrack.GetPdgCode() != 13 || !mcTrack.isPrimary()) {
          continue;
        }
        const auto mcTrackRefs = mcReader.getTrackRefs(iSource, iEvent, iTrack);
        if (mcTrackRefs.size() == 0 || mcTrackRefs[0].getDetectorId() != o2::detectors::DetID::MCH) {
          continue;
        }
        FillResiduals(mcTrackRefs[0], mcTrack, residuals);
      }
    }
  }

  DrawResiduals(residuals);
}

//_________________________________________________________________________________________________
void CreateResiduals(std::vector<TH1*>& histos)
{
  /// create histograms holding residuals between simulated and reconstructed variables
  histos.emplace_back(new TH1F("dp", "dp;dp (GeV/c)", 501, -50.1, 50.1));
  histos.emplace_back(new TH1F("dpT", "dpT;dp_{T} (GeV/c)", 501, -5.01, 5.01));
  histos.emplace_back(new TH1F("deta", "deta;deta", 1001, -0.101, 0.101));
  histos.emplace_back(new TH1F("dphi", "dphi;dphi (rad)", 1001, -0.5005, 0.5005));
}

//_________________________________________________________________________________________________
void FillResiduals(const TrackReference& mcTrackRef, const MCTrack& mcTrack, std::vector<TH1*>& histos)
{
  /// fill residuals between simulated variables at vertex and at first cluster
  double p = mcTrackRef.P();
  double eta = 0.5 * log((p + mcTrackRef.Pz()) / (p - mcTrackRef.Pz()));
  double phi = pi() + atan2(-mcTrackRef.Py(), -mcTrackRef.Px());
  histos[0]->Fill(p - mcTrack.GetP());
  histos[1]->Fill(mcTrackRef.Pt() - mcTrack.GetPt());
  histos[2]->Fill(eta - mcTrack.GetEta());
  histos[3]->Fill(phi - mcTrack.GetPhi());
}

//_________________________________________________________________________________________________
void DrawResiduals(std::vector<TH1*>& histos)
{
  /// draw residuals between simulated and reconstructed variables
  TCanvas* c = new TCanvas("residuals", "residuals", 30, 30, 600, 600);
  c->Divide(2, 2);
  for (int i = 0; i < 4; ++i) {
    c->cd(i + 1);
    gPad->SetLogy();
    histos[i]->Draw();
  }
}
