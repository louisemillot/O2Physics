#include "PWGJE/DataModel/SlimTables.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisHelpers.h"

#include "PWGJE/DataModel/Jet.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SlimTablesProducer {

  void init(InitContext&)
  {
  }

  // --- Tables produites ---
  Produces<o2::aod::SlimTracks> slimTracks;

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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SlimTablesProducer>(cfgc, TaskName{"slim-tables-producer"})
  };
}
