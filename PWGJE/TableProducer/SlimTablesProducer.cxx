#include "PWGJE/DataModel/SlimTables.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisHelpers.h"

#include "PWGJE/DataModel/Jet.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SlimTablesProducer {

  // --- Tables produites ---
  Produces<o2::aod::SlimTracks> slimTracks;
  Produces<o2::aod::SlimParticles> slimParticles;

  // ------------------------------
  // SlimTracks
  // ------------------------------
  void processTracks(aod::SlimTracks const& tracks)
  {
    for (auto& tr : tracks) {
      slimTracks(tr.px(), tr.py(), tr.pz());
    }
  }
  PROCESS_SWITCH(SlimTablesProducer, processTracks,
                 "Produce slim track table", true);

  // ------------------------------
  // Slim MC Particles
  // ------------------------------
  void processParticles(aod::SlimParticles const& particles)
  {
    for (auto& p : particles) {
      slimParticles(p.px(), p.py(), p.pz(), p.e());
    }
  }
  PROCESS_SWITCH(SlimTablesProducer, processParticles,
                 "Produce slim MC particle table", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SlimTablesProducer>(cfgc, TaskName{"slim-tables-producer"})
  };
}
