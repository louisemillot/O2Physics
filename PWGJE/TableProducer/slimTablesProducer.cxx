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

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/ASoA.h>
// #include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <fairlogger/Logger.h>

#include <Rtypes.h>

#include <string>
#include <unordered_map>
#include <vector>

namespace o2::aod
{
namespace slimcollision
{
DECLARE_SOA_COLUMN(Weight, weight, float);
}
DECLARE_SOA_TABLE(SlimCollisions, "AOD", "SlimCollisions",
                  o2::soa::Index<>,
                  o2::aod::collision::PosZ,
                  o2::aod::collision::CollisionTime,
                  slimcollision::Weight);
using SlimCollision = SlimCollisions::iterator;
namespace slmccollision
{
DECLARE_SOA_COLUMN(McWeight, mcWeight, float);
}
DECLARE_SOA_TABLE(SlMcCollisions, "AOD", "SlMcCollisions",
                  o2::soa::Index<>,
                  o2::aod::mccollision::PosZ,
                  slmccollision::McWeight);
using SlMcCollision = SlMcCollisions::iterator;
namespace slimtracks
{
DECLARE_SOA_INDEX_COLUMN(SlimCollision, slimCollision);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(E, e, float);
} // namespace slimtracks
DECLARE_SOA_TABLE(SlimTracks, "AOD", "SlimTracks",
                  o2::soa::Index<>,
                  slimtracks::SlimCollisionId,
                  slimtracks::Px,
                  slimtracks::Py,
                  slimtracks::Pz,
                  slimtracks::E);
using SlimTrack = SlimTracks::iterator;
namespace slimparticles
{
DECLARE_SOA_INDEX_COLUMN(SlMcCollision, slMcCollision);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(E, e, float);
} // namespace slimparticles
DECLARE_SOA_TABLE(SlimParticles, "AOD", "SlimParticles",
                  o2::soa::Index<>,
                  slimparticles::SlMcCollisionId,
                  slimparticles::Px,
                  slimparticles::Py,
                  slimparticles::Pz,
                  slimparticles::E);
using SlimParticle = SlimParticles::iterator;
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SlimTablesProducer {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<bool> checkCentFT0M{"checkCentFT0M", false, "0: centFT0C as default, 1: use centFT0M estimator"};
  Configurable<float> centralityMin{"centralityMin", -999, ""};
  Configurable<float> centralityMax{"centralityMax", 999, ""};
  Configurable<float> minPt{"minPt", 0.15, "min pT to save"};
  Configurable<float> maxPt{"maxPt", 200.0, "max pT to save"};
  Configurable<float> minEta{"minEta", -0.9, "min eta to save"};
  Configurable<float> maxEta{"maxEta", 0.9, "max eta to save"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> trackDcaZmax{"trackDcaZmax", 99, "additional cut on dcaZ to PV for tracks; uniformTracks in particular don't cut on this at all"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "Event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections; other option: uniformTracks"};
  Configurable<int> minTPCNClsCrossedRows{"minTPCNClsCrossedRows", 80, "min TPC crossed rows"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false, "flag to choose to reject min. bias gap events; jet-level rejection can also be applied at the jet finder level for jets only, here rejection is applied for collision and track process functions for the first time, and on jets in case it was set to false at the jet finder level"};
  Configurable<bool> applyRCTSelections{"applyRCTSelections", true, "decide to apply RCT selections"};

  std::vector<int> eventSelectionBits;
  std::unordered_map<int, int> recoGlobalToSlim;
  // Service<o2::framework::O2DatabasePDG> pdgDatabase;
  int trackSelection = -1;
  bool doSumw2 = false;

  void init(InitContext&)
  {
    doSumw2 = skipMBGapEvents; // true or false : storage of square erros when jet-jet

    AxisSpec centralityAxis = {1200, -10., 110., "Centrality"};

    histos.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
    histos.add("h_pt_particles", "Constituent pT; pT (GeV/c); Counts", kTH1F, {{200, 0.0, 200.0}});
    histos.add("h2_centrality_collisions", "event status vs. centrality;entries;centrality", {HistType::kTH2F, {centralityAxis, {4, 0.0, 4.0}}}, doSumw2);
    auto hColl = histos.get<TH1>(HIST("h_collisions"));
    hColl->GetXaxis()->SetBinLabel(1, "All");
    hColl->GetXaxis()->SetBinLabel(2, "eventSelection");

    histos.add("h_mcCollMCD_counts_weight", "MC event status;event status;weighted entries", {HistType::kTH1F, {{5, 0.0, 5.0}}});
    auto hMCD = histos.get<TH1>(HIST("h_mcCollMCD_counts_weight"));
    hMCD->GetXaxis()->SetBinLabel(1, "All");
    hMCD->GetXaxis()->SetBinLabel(2, "eventSelectionBits + skipMBGapEvents + applyRCTSelections ");

    histos.add("h_mcCollMCP_counts_weight", "MC event status;event status;weighted entries", {HistType::kTH1F, {{7, 0.0, 7.0}}});
    auto hMCP = histos.get<TH1>(HIST("h_mcCollMCP_counts_weight"));
    hMCP->GetXaxis()->SetBinLabel(1, "All");
    hMCP->GetXaxis()->SetBinLabel(2, "mcColl + skipMBGapEvents + applyRCTSelections");
    hMCP->GetXaxis()->SetBinLabel(3, "Zvertex");

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
  }

  Produces<o2::aod::SlimCollisions> slimCollisions;
  Produces<o2::aod::SlMcCollisions> slimMcCollisions;
  Produces<o2::aod::SlimTracks> slimTracks;
  Produces<o2::aod::SlimParticles> slimParticles;

  Filter trackFilter = (aod::jtrack::pt >= minPt && aod::jtrack::pt < maxPt && aod::jtrack::eta > minEta && aod::jtrack::eta < maxEta);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut &&
                      (checkCentFT0M ? aod::jcollision::centFT0M : aod::jcollision::centFT0C) >= centralityMin &&
                      (checkCentFT0M ? aod::jcollision::centFT0M : aod::jcollision::centFT0C) < centralityMax);
  Filter mcCollisionFilter = (nabs(aod::jmccollision::posZ) < vertexZCut && aod::jmccollision::centFT0M >= centralityMin && aod::jmccollision::centFT0M < centralityMax); // no centFT0C for mccollisions, using centFT0M for both
  Filter particleCuts = (aod::jmcparticle::pt >= minPt && aod::jmcparticle::pt < maxPt && aod::jmcparticle::eta > minEta && aod::jmcparticle::eta < maxEta);

  Preslice<aod::JetTracksMCD> perCollisionTracks = aod::jtrack::collisionId;
  Preslice<aod::JetParticles> perMcCollisionParticles = aod::jmcparticle::mcCollisionId;

  void processData(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                   soa::Filtered<aod::JetTracks> const& tracks)
  {
    histos.fill(HIST("h_collisions"), 0.5);
    float centrality = -1.0;
    checkCentFT0M ? centrality = collision.centFT0M() : centrality = collision.centFT0C();
    histos.fill(HIST("h2_centrality_collisions"), centrality, 0.5, 1.0);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, false, applyRCTSelections)) {
      return;
    }
    histos.fill(HIST("h_collisions"), 1.5);
    slimCollisions(collision.posZ(), collision.collisionTime(), 1.0);
    auto ts = collision.collisionTime();
    auto slimCollIndex = slimCollisions.lastIndex();
    LOG(INFO) << "Collision PbPB globalindex = " << collision.globalIndex();
    LOG(INFO) << "Collision PbPb slimIndex = " << slimCollisions.lastIndex();
    LOG(INFO) << "Collision time PbPb = " << ts;
    for (const auto& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      float mass = jetderiveddatautilities::mPion;
      float p = track.pt() * std::cosh(track.eta());
      float energy = std::sqrt(p * p + mass * mass);
      slimTracks(slimCollIndex, track.px(), track.py(), track.pz(), energy);
    }
  }
  PROCESS_SWITCH(SlimTablesProducer, processData, "process collisions and tracks for data", false);

  // void processMCD(soa::Filtered<aod::JetCollisionsMCD>::iterator const& collision,
  //                 aod::JetMcCollisions const&, // join the weight
  //                 soa::Filtered<aod::JetTracksMCD> const& tracks)
  // {
  //   float eventWeight = collision.mcCollision_as<aod::JetMcCollisions>().weight();
  //   histos.fill(HIST("h_mcCollMCD_counts_weight"), 0.5, eventWeight);
  //   float centrality = -1.0;
  //   checkCentFT0M ? centrality = collision.centFT0M() : centrality = collision.centFT0C();
  //   histos.fill(HIST("h2_centrality_MCD"), centrality, 0.5, eventWeight);
  //   if (!collision.has_mcCollision()) {
  //     LOG(warning) << "MC collision not found for collision with global ID " << collision.globalIndex();
  //     return;
  //   }
  //   if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents, applyRCTSelections)) {
  //     return;
  //   }
  //   histos.fill(HIST("h_mcCollMCD_counts_weight"), 1.5, eventWeight);
  //   slimCollisions(collision.posZ());
  //   auto slimCollIndex = slimCollisions.lastIndex();
  //   recoGlobalToSlim[collision.globalIndex()] = slimCollIndex;
  //   for (const auto& track : tracks) {
  //     if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
  //       continue;
  //     }
  //     float mass = jetderiveddatautilities::mPion;
  //     float p = track.pt() * std::cosh(track.eta());
  //     float energy = std::sqrt(p * p + mass * mass);
  //     slimTracks(slimCollIndex, track.px(), track.py(), track.pz(), energy);
  //   }
  // }
  // PROCESS_SWITCH(SlimTablesProducer, processMCD, "process collisions and tracks for MCD", false);

  // void processMCP(soa::Filtered<o2::aod::JetMcCollisions>::iterator const& mcCollision,
  //                 soa::SmallGroups<aod::JetCollisionsMCD> const& collisions, // SmallGroups contains and access the rec collisions associated to the mc collision
  //                 soa::Filtered<aod::JetParticles> const& particles)
  // {
  //   float eventWeight = mcCollision.weight();
  //   float centrality = mcCollision.centFT0M(); // checkCentFT0M ? centrality = mccollision.centFT0M() : centrality = mccollision.centFT0C();
  //   histos.fill(HIST("h_mcCollMCP_counts_weight"), 0.5, eventWeight);
  //   histos.fill(HIST("h2_centrality_MCP"), centrality, 0.5, eventWeight);

  //   if (std::abs(mcCollision.posZ()) > vertexZCut) {
  //     return;
  //   }
  //   histos.fill(HIST("h_mcCollMCP_counts_weight"), 1.5, eventWeight);
  //   if (collisions.size() < 1) {
  //     return;
  //   }
  //   histos.fill(HIST("h_mcCollMCP_counts_weight"), 2.5, eventWeight);
  //   bool hasSel8Coll = false;
  //   int matchedSlimCollId = -1;
  //   for (auto const& collision : collisions) {
  //     if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents, applyRCTSelections)) { // look if the rec collision associated to the mc collision passes the event selection
  //       hasSel8Coll = true;
  //       int globalId = collision.globalIndex();
  //       if (recoGlobalToSlim.find(globalId) != recoGlobalToSlim.end()) { // find the globalId of collision (original AO2D) till the end
  //         matchedSlimCollId = recoGlobalToSlim[globalId];                // if globalId found in map -> get the corresponding slimCollId
  //         LOG(info) << "==== MATCH FOUND ====";
  //         LOG(info) << "global MC coll ID          = " << mcCollision.globalIndex();
  //         LOG(info) << "global coll ID        = " << globalId;
  //         LOG(info) << "Slim coll ID        = " << matchedSlimCollId;
  //         break; // on prend le premier valide
  //       }
  //     }
  //   }
  //   if (!hasSel8Coll) {
  //     return;
  //   }
  //   histos.fill(HIST("h_mcCollMCP_counts_weight"), 3.5, eventWeight);
  //   if (!jetderiveddatautilities::selectMcCollision(mcCollision, skipMBGapEvents, applyRCTSelections)) {
  //     return;
  //   }
  //   histos.fill(HIST("h_mcCollMCP_counts_weight"), 4.5, eventWeight);
  //   if (matchedSlimCollId < 0) {
  //     return;
  //   }
  //   slimMcCollisions(mcCollision.posZ(), matchedSlimCollId);
  //   auto slimMcCollIndex = slimMcCollisions.lastIndex();
  //   LOG(info) << "Slim Mc coll ID         = " << slimMcCollIndex;
  //   LOG(info) << "======================";
  //   for (const auto& particle : particles) {
  //     // auto pdgParticle = pdgDatabase->GetParticle(particle.pdgCode());
  //     // float massParticle = pdgParticle ? pdgParticle->Mass() : jetderiveddatautilities::mPion;
  //     // float energyParticle = std::sqrt(particle.px() * particle.px() + particle.py() * particle.py() + particle.pz() * particle.pz() + massParticle * massParticle);
  //     slimParticles(slimMcCollIndex, particle.px(), particle.py(), particle.pz(), particle.energy());
  //   }
  // }
  // PROCESS_SWITCH(SlimTablesProducer, processMCP, "process mccollisions and mcparticles for MCD", false);

  void processMC(soa::Filtered<aod::JetCollisionsMCD> const& collisions,
                 aod::JetMcCollisions const&, // join the weight
                 soa::Filtered<aod::JetTracksMCD> const& tracks,
                 soa::Filtered<aod::JetParticles> const& particles)
  {
    for (auto& collision : collisions) {
      float eventWeight = collision.weight();
      LOG(info) << "Processing collision with global ID " << collision.globalIndex();
      if (!collision.has_mcCollision()) {
        return;
      }
      auto mcColl = collision.mcCollision(); // corresponding MC coll
      histos.fill(HIST("h_mcCollMCD_counts_weight"), 0.5, eventWeight);
      histos.fill(HIST("h_mcCollMCP_counts_weight"), 0.5, eventWeight);
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents, applyRCTSelections)) {
        return;
      }
      if (!jetderiveddatautilities::selectMcCollision(mcColl, skipMBGapEvents, applyRCTSelections)) {
        return;
      }
      histos.fill(HIST("h_mcCollMCD_counts_weight"), 1.5, eventWeight);
      histos.fill(HIST("h_mcCollMCP_counts_weight"), 1.5, eventWeight);
      if (std::abs(mcColl.posZ()) > vertexZCut)
        return;
      histos.fill(HIST("h_mcCollMCP_counts_weight"), 2.5, eventWeight);

      float eventMCWeight = mcColl.weight();
      LOG(INFO) << "eventWeight =" << eventWeight;
      LOG(INFO) << "eventMCWeight =" << eventMCWeight;
      LOG(INFO) << "MC collision global ID = " << mcColl.globalIndex()
                << " coll global ID = " << collision.globalIndex();
      slimCollisions(collision.posZ(), collision.collisionTime(), eventWeight);
      auto slimCollIndex = slimCollisions.lastIndex();
      auto ts = collision.collisionTime();
      LOG(INFO) << "Collision pp globalindex = " << collision.globalIndex();
      LOG(INFO) << "Collision pp slimIndex = " << slimCollisions.lastIndex();
      LOG(INFO) << "Collision time pp = " << ts;
      auto slicedTracks = tracks.sliceBy(perCollisionTracks, collision.globalIndex()); // tracks associated to the rec collision
      for (const auto& track : slicedTracks) {
        if (!jetderiveddatautilities::selectTrack(track, trackSelection))
          continue;
        float mass = jetderiveddatautilities::mPion;
        float p = track.pt() * std::cosh(track.eta());
        float energy = std::sqrt(p * p + mass * mass);
        slimTracks(slimCollIndex, track.px(), track.py(), track.pz(), energy);
      }
      slimMcCollisions(mcColl.posZ(), eventMCWeight);
      auto slimMcCollIndex = slimMcCollisions.lastIndex();
      LOG(INFO) << "slimMcCollIndex = " << slimMcCollisions.lastIndex();
      auto slicedParticles = particles.sliceBy(perMcCollisionParticles, mcColl.globalIndex()); // particles associated to the mc collision
      for (const auto& particle : slicedParticles) {
        if (!particle.isPhysicalPrimary())
          continue;
        histos.fill(HIST("h_pt_particles"), particle.pt(), eventMCWeight);
        slimParticles(slimMcCollIndex, particle.px(), particle.py(), particle.pz(), particle.energy());
      }
    }
  }
  PROCESS_SWITCH(SlimTablesProducer, processMC, "process collisions & tracks, MCcollisions & particles for MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SlimTablesProducer>(cfgc)};
}
