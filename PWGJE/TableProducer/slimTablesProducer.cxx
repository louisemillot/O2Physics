#include "PWGJE/DataModel/SlimTables.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisHelpers.h"

#include "PWGJE/DataModel/Jet.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct slimTablesProducer {

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
  PROCESS_SWITCH(slimTablesProducer, processTracks,
                 "Produce slim track table", true);

  // ------------------------------
  // Slim MC Particles
  // ------------------------------
  void processParticles(aod::SlimParticles const& particles)
  {
    for (auto& p : particles) {
      slimParticles(p.px(), p.py(), p.pz(), p.energy());
    }
  }
  PROCESS_SWITCH(slimTablesProducer, processParticles,
                 "Produce slim MC particle table", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<slimTablesProducer>(cfgc, TaskName{"slim-tables-producer"})
  };
}