#include "PWGJE/DataModel/SlimTables.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"


using namespace o2::aod;

// ---------------- SlimTracksProducer ----------------
struct SlimTracksProducer {
  Produces<SlimTracks> slimTracks;

  void process(Tracks const& tracks) {
    for (auto& tr : tracks) {
      slimTracks(tr.px(), tr.py(), tr.pz());
    }
  }
};

// ---------------- SlimParticlesProducer ----------------
struct SlimParticlesProducer {
  Produces<SlimParticles> slimParticles;

  void process(Particles const& particles) {
    for (auto& p : particles) {
      slimParticles(p.px(), p.py(), p.pz(), p.E());
    }
  }
};

// ---------------- SlimCollisionProducer ----------------
struct SlimCollisionProducer {
  Produces<SlimCollision> slimCollision;

  void process(Collision const& coll) {
    slimCollision(coll.posZ(), coll.Centrality(), coll.EventSel(), coll.EventWeight());
  }
};

// ---------------- SlimMcCollisionProducer ----------------
struct SlimMcCollisionProducer {
  Produces<SlimMcCollision> slimMcCollision;

  void process(McCollision const& mcColl) {
    slimMcCollision(mcColl.posZ(), mcColl.centrality(), mcColl.eventSel(), mcColl.eventWeight());
  }
};
