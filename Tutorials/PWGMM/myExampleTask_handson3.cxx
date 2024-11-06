#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct myExampleTask {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Filter trackDCA = nabs(aod::track::dcaXY)<0.2f;

  using myCompleteTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels>;
  using myFilteredTracks = soa::Filtered<myCompleteTracks>;

  void init(InitContext const&)
  {
    // Define axes
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisPt{nBinsPt, 0, 10, "p_{T}"};
    const AxisSpec axisDeltaPt{100, -1.0, +1.0, "#Delta(p_{T})"};
    // Create histograms
    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt});
    // Adding histograms for pions, kaons, and protons
    histos.add("ptHistogramPion", "ptHistogramPion", kTH1F, {axisPt});
    histos.add("ptHistogramKaon", "ptHistogramKaon", kTH1F, {axisPt});
    histos.add("ptHistogramProton", "ptHistogramProton", kTH1F, {axisPt});
    // Adding generated pT histograms
    histos.add("ptGeneratedPion", "ptGeneratedPion", kTH1F, {axisPt}); 
    histos.add("ptGeneratedKaon", "ptGeneratedKaon", kTH1F, {axisPt}); 
    histos.add("ptGeneratedProton", "ptGeneratedProton", kTH1F, {axisPt});
  }

  void processReconstructed(aod::Collision const& , myFilteredTracks const& tracks, aod::McParticles const& ){
    histos.fill(HIST("eventCounter"), 0.5);
    for (auto& track : tracks) {
        if (track.has_mcParticle()) {  // Vérifier si la piste a une particule MC associée
        auto mcParticle = track.mcParticle();  // De-référencer la particule MC
        // Check if the MC particle is a physical primary and within |y| < 0.5 (mid-rapidity)
        if(mcParticle.isPhysicalPrimary() && fabs(mcParticle.y())<0.5){ // do this in the context of the track ! (context matters!!!) 
            if(abs(mcParticle.pdgCode())==211) histos.fill(HIST("ptHistogramPion"), mcParticle.pt()); 
            if(abs(mcParticle.pdgCode())==321) histos.fill(HIST("ptHistogramKaon"), mcParticle.pt()); 
            if(abs(mcParticle.pdgCode())==2212) histos.fill(HIST("ptHistogramProton"), mcParticle.pt());
        }
        }
    }
  }
  PROCESS_SWITCH(myExampleTask, processReconstructed, "process reconstructed information", true);
  void processSimulated(aod::McParticles const& mcParticles) {
    for (const auto& mcParticle : mcParticles) {
        if(mcParticle.isPhysicalPrimary() && fabs(mcParticle.y())<0.5){ // watch out for context!!!
            if(abs(mcParticle.pdgCode())==211) histos.fill(HIST("ptGeneratedPion"), mcParticle.pt()); 
            if(abs(mcParticle.pdgCode())==321) histos.fill(HIST("ptGeneratedKaon"), mcParticle.pt()); 
            if(abs(mcParticle.pdgCode())==2212) histos.fill(HIST("ptGeneratedProton"), mcParticle.pt());    
        } 
    }
  }
  PROCESS_SWITCH(myExampleTask, processSimulated, "process pure simulation information", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myExampleTask>(cfgc)};
}
