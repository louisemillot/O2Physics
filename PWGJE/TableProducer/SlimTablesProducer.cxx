#include "PWGJE/DataModel/SlimTables.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisHelpers.h"
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SlimTablesProducer {

  void init(InitContext&)
  {
  }

  // --- Tables produites ---
  Produces<o2::aod::SlimTracks> slimTracks;
  Produces<o2::aod::SlimParticles> slimParticles;


  // ------------------------------
  // SlimTracks
  // ------------------------------
  void processTracks(aod::SlimTracks const& tracks)
  {
    for (auto& tr : tracks) {
      slimTracks(tr.pt());
    }
  }
  PROCESS_SWITCH(SlimTablesProducer, processTracks,"Produce slim track table", true);

  void processParticles(aod::McParticles const& parts)
  {
    for (auto& p : parts) {
      slimParticles(p.px(), p.py(), p.pz());
    }
  }
  PROCESS_SWITCH(SlimTablesProducer, processParticles, "produce slim particles", true);



};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SlimTablesProducer>(cfgc, TaskName{"slim-tables-producer"})
  };
}
