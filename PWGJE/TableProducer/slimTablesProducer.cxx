// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file slimTablesProducer.cxx
/// \brief Task to produce a reduced version of Tables for tracks, collisions, mcparticles and mccollisions.
/// \author Millot Louise <louise.millot@cern.ch>

#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGJE/Core/JetDQUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetV0Utilities.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/EMCALMatchedCollisions.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/SlimTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGUD/Core/SGCutParHolder.h"
#include "PWGUD/Core/SGSelector.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "ReconstructionDataFormats/Vertex.h"
#include <CommonConstants/MathConstants.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SlimTablesProducer {

  void init(InitContext&)
  {
  }

  Produces<o2::aod::SlimCollisions> slimCollisions;
  // Produces<o2::aod::SlimMcCollisions> slimMcCollisions;
  // Produces<o2::aod::SlimTracks> slimTracks;
  // Produces<o2::aod::SlimParticles> slimParticles;
  // Preslice<aod::JetTracks> tracksPerCollision = aod::jtrack::collisionId;

  // void processCollision(aod::JetCollisions const& collisions)
  // {
  //   for (const auto& coll : collisions) {
  //     // slimCollisions(coll.posZ(), coll.centFT0C(), coll.centFT0M(), coll.weight(), coll.eventSel(), coll.trackOccupancyInTimeRange());
  //   }
  // }
  // PROCESS_SWITCH(SlimTablesProducer, processCollision, "Produce slim collision table", false);

  void processCollision(aod::Collisions::iterator const& collision)
  {
    slimCollisions(collision.posZ());
  }
  PROCESS_SWITCH(SlimTablesProducer, processCollision, "Produce slim collision table", false);

  // void processMcCollision(aod::JetMcCollisions const& mccollisions)
  // {
  //   for (const auto& mccoll : mccollisions) {
  //     slimMcCollisions(mccoll.posZ(), mccoll.centFT0M(), mccoll.weight(), mccoll.accepted(), mccoll.ptHard());
  //   }
  // }
  // PROCESS_SWITCH(SlimTablesProducer, processMcCollision, "Produce slim mc collision table", false);

  // void processTracks(aod::Collisions::iterator const& collision,
  //                    aod::Tracks const& tracks)
  // {
  //   auto tracksInCollision = tracks.sliceBy(trackPerColl, collision.globalIndex());
  //   for (const auto& trk : tracksInCollision) {
  //     // slimTracks(trk.collision(), trk.pt(), trk.eta(), trk.phi(), trk.dcaXY());
  //     slimTracks(trk.collision(), trk.pt(), trk.eta(), trk.phi(), trk.px(), trk.py(), trk.pz());
  //   }
  // }
  // PROCESS_SWITCH(SlimTablesProducer, processTracks, "Produce slim track table", true);

  // void processTracks(aod::JetCollisions::iterator const& coll,
  //                    aod::JetTracks const& tracks)
  // {
  //   auto tracksInColl = tracks.sliceBy(tracksPerCollision, coll.globalIndex());
  //   for (const auto& trk : tracksInColl) {
  //     // slimTracks(trk.collision(), trk.pt(), trk.eta(), trk.phi(), trk.dcaXY());
  //     slimTracks(trk.collisionId(), trk.pt(), trk.eta(), trk.phi(), trk.px(), trk.py(), trk.pz());
  //   }
  // }
  // PROCESS_SWITCH(SlimTablesProducer, processTracks, "Produce slim track table", true);

  // void processParticles(aod::JetParticles const& parts)
  // {
  //   for (const auto& p : parts) {
  //     slimParticles(p.mcCollision(), p.pt(), p.eta(), p.phi());
  //   }
  // }
  // PROCESS_SWITCH(SlimTablesProducer, processParticles, "Produce slim particles", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SlimTablesProducer>(cfgc)};
}
