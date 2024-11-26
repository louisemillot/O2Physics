#include "ReconstructionDataFormats/Track.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DataModel/DerivedExampleTable.h"

// Header pour le filtre et partition
#include "Framework/ASoAHelpers.h"
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct DerivedBasicConsumer {
  
  // Filtre sur la position Z du vertex de collision
  Filter collZfilter = nabs(aod::collision::posZ) < 10.0f;
  // Ajout d'un SliceCache pour gérer les slices de partitions
  SliceCache cache;
  // Définition des partitions pour les pistes associées et les pistes trigger
  Partition<aod::DrTracks> associatedTracks = aod::exampleTrackSpace::pt < 6.0f && aod::exampleTrackSpace::pt > 4.0f;
  Partition<aod::DrTracks> triggerTracks = aod::exampleTrackSpace::pt > 6.0f;

  /// Fonction pour calculer delta-phi
  Double_t ComputeDeltaPhi(Double_t phi1, Double_t phi2)
  {
    Double_t deltaPhi = phi1 - phi2;
    if (deltaPhi < -TMath::Pi() / 2.) {
      deltaPhi += 2. * TMath::Pi();
    }
    if (deltaPhi > 3 * TMath::Pi() / 2.) {
      deltaPhi -= 2. * TMath::Pi();
    }
    return deltaPhi;
  }

  // Histogram registry: un objet pour contenir vos histogrammes
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext const&)
  {
    // Définition des axes et des histogrammes
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisVertexZ{100, -20, 20, "Vtx Z (cm)"};
    const AxisSpec axisPt{100, 0, 20, "pT (GeV/c)"};
    const AxisSpec axisDeltaPhi{100, -0.5*TMath::Pi(), +1.5*TMath::Pi(), "#Delta#phi"};
    const AxisSpec axisDeltaEta{100, -1.0, +1.0, "#Delta#eta"};
    histos.add("eventCounter", "Event Counter", kTH1F, {axisCounter});
    histos.add("eventVertexZ", "Event Vertex Z Position", kTH1F, {axisVertexZ});
    histos.add("ptAssoHistogram", "Associated Tracks pT", kTH1F, {axisPt});
    histos.add("ptTrigHistogram", "Trigger Tracks pT", kTH1F, {axisPt});
    histos.add("correlationFunction", "correlationFunction", kTH1F, {axisDeltaPhi});
    histos.add("correlationFunction2d", "correlationFunction2d", kTH2F, {axisDeltaPhi, axisDeltaEta});
  }

  void process(soa::Filtered<aod::DrCollisions>::iterator const& collision, aod::DrTracks const& tracks)
  {
    // Remplir l'histogramme de comptage d'événements
    histos.fill(HIST("eventCounter"), 0.5);
    // Remplir l'histogramme de la position Z du vertex
    histos.fill(HIST("eventVertexZ"), collision.posZ());
    // Grouper les partitions pour cette collision spécifique
    auto assoTracksThisCollision = associatedTracks->sliceByCached(aod::exampleTrackSpace::drCollisionId, collision.globalIndex(), cache);
    auto trigTracksThisCollision = triggerTracks->sliceByCached(aod::exampleTrackSpace::drCollisionId, collision.globalIndex(), cache);
    // Remplir les histogrammes de QA pour les pistes associées
    for (auto& track : assoTracksThisCollision) {
      histos.fill(HIST("ptAssoHistogram"), track.pt());
    }
    // Remplir les histogrammes de QA pour les pistes trigger
    for (auto& track : trigTracksThisCollision) {
      histos.fill(HIST("ptTrigHistogram"), track.pt());
    }
    for (auto& [trigger, associated] : combinations(o2::soa::CombinationsFullIndexPolicy(trigTracksThisCollision, assoTracksThisCollision))) {
      histos.fill(HIST("correlationFunction"), ComputeDeltaPhi(trigger.phi(),associated.phi()));
      histos.fill(HIST("correlationFunction2d"), ComputeDeltaPhi(trigger.phi(),associated.phi()), trigger.eta() - associated.eta());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<DerivedBasicConsumer>(cfgc, TaskName{"derived-basic-consumer"})};
  return workflow;
}


