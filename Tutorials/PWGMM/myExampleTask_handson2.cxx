#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
//#include "Common/DataModel/McParticles.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct myExampleTask {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};

  Filter trackDCA = nabs(aod::track::dcaXY)<0.2f; ///// soit on commente cette ligne mais on doit 1)commenter l19, 2)dé-commenter l40, 3)changer l35 avec myCompleteTracks au lieu de myFilteredTracks

  using myCompleteTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels>;
  using myFilteredTracks = soa::Filtered<myCompleteTracks>;

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisPt{nBinsPt, 0, 10, "p_{T}"};
    const AxisSpec axisDeltaPt{100, -1.0, +1.0, "#Delta(p_{T})"};
    // create histograms
    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt});
    histos.add("ptResolution", "ptResolution", kTH2F, {axisPt, axisDeltaPt}); 
  }
  void process(aod::Collision const& , myFilteredTracks const& tracks, aod::McParticles const& )
  {
    histos.fill(HIST("eventCounter"), 0.5);
    for (auto& track : tracks) {
      if( track.tpcNClsCrossedRows() < 70 ) continue; //badly tracked
      //if( fabs(track.dcaXY()) > 0.2) continue; //doesn’t point to primary vertex
      histos.fill(HIST("etaHistogram"), track.eta());
      histos.fill(HIST("ptHistogram"), track.pt());
      if (track.has_mcParticle()) {  // Vérifier si la piste a une particule MC associée
        auto mcParticle = track.mcParticle();  // De-référencer la particule MC
        histos.fill(HIST("ptResolution"), track.pt(), track.pt() - mcParticle.pt()); // Remplir l'histogramme de résolution
      }
    }
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myExampleTask>(cfgc)};
}

