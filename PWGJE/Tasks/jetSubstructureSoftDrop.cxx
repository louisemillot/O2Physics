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

// jet analysis tasks (subscribing to jet finder task)
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
//

#include <vector>
#include <utility>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/Core/JetSubstructureUtilities.h"

#include "EventFiltering/filterTables.h"


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct JetSubstructureTask {

  Produces<aod::ChargedSPs> jetSplittingsDataTable;
  Produces<aod::ChargedMCDetectorLevelSPs> jetSplittingsMCDTable;
  Produces<aod::ChargedMCParticleLevelSPs> jetSplittingsMCPTable;
  Produces<aod::ChargedEventWiseSubtractedSPs> jetSplittingsDataSubTable;
  // Produces<aod::ChargedEventWiseSubtractedSPs> jetSplittingsMCDSubTable;
  // Produces<aod::ChargedEventWiseSubtractedSPs> jetSplittingsMCPSubTable;

  Configurable<float> zCut{"zCut", 0.1, "soft drop z cut"};
  Configurable<float> beta{"beta", 0.0, "soft drop beta"};
  Configurable<float> kappa{"kappa", 1.0, "angularity kappa"};
  Configurable<float> alpha{"alpha", 1.0, "angularity alpha"};
  Configurable<bool> doPairBkg{"doPairBkg", true, "save bkg pairs"};
  Configurable<float> pairConstituentPtMin{"pairConstituentPtMin", 1.0, "pt cut off for constituents going into pairs"};
  Configurable<float> centralityMin{"centralityMin", -999, ""};
  Configurable<float> centralityMax{"centralityMax", 999, ""};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> trackQAEtaMin{"trackQAEtaMin", -0.9, "minimum eta acceptance for tracks in the processTracks QA"};
  Configurable<float> trackQAEtaMax{"trackQAEtaMax", 0.9, "maximum eta acceptance for tracks in the processTracks QA"};
  Configurable<float> trackQAPtMin{"trackQAPtMin", 0.15, "minimum pT acceptance for tracks in the processTracks QA"};
  Configurable<float> trackQAPtMax{"trackQAPtMax", 100.0, "maximum pT acceptance for tracks in the processTracks QA"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<float> ptLeadingTrackCut{"ptLeadingTrackCut", 5.0f, "Leading track cut"};
  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum track occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<int> trackOccupancyInTimeRangeMin{"trackOccupancyInTimeRangeMin", -999999, "minimum track occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false, "flag to choose to reject min. bias gap events; jet-level rejection can also be applied at the jet finder level for jets only, here rejection is applied for collision and track process functions for the first time, and on jets in case it was set to false at the jet finder level"};


  Service<o2::framework::O2DatabasePDG> pdg;
  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  std::vector<float> energyMotherVec;
  std::vector<float> ptLeadingVec;
  std::vector<float> ptSubLeadingVec;
  std::vector<float> thetaVec;
  std::vector<float> nSub;
  std::vector<float> pairJetPtVec;
  std::vector<float> pairJetEnergyVec;
  std::vector<float> pairJetThetaVec;
  std::vector<float> pairJetPerpCone1PtVec;
  std::vector<float> pairJetPerpCone1EnergyVec;
  std::vector<float> pairJetPerpCone1ThetaVec;
  std::vector<float> pairPerpCone1PerpCone1PtVec;
  std::vector<float> pairPerpCone1PerpCone1EnergyVec;
  std::vector<float> pairPerpCone1PerpCone1ThetaVec;
  std::vector<float> pairPerpCone1PerpCone2PtVec;
  std::vector<float> pairPerpCone1PerpCone2EnergyVec;
  std::vector<float> pairPerpCone1PerpCone2ThetaVec;
  float angularity;
  float leadingConstituentPt;
  float perpConeRho;

  HistogramRegistry registry;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  void init(InitContext const&)
  {
    registry.add("h2_jet_pt_jet_zg", ";#it{p}_{T,jet} (GeV/#it{c});#it{z}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_rg", ";#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_thetag", ";#it{p}_{T,jet} (GeV/#it{c});#it{theta}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_nsd", ";#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});

    registry.add("h2_jet_pt_part_jet_zg_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{z}_{g}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_part_jet_rg_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{R}_{g}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_part_jet_thetag_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{theta}_{g}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_part_jet_nsd_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{n}_{SD}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});

    registry.add("h2_jet_pt_jet_zg_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{z}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_rg_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_thetag_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{theta}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_nsd_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});

    registry.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
    registry.add("h_jet_pT", "jet pT", {HistType::kTH1F, {{4, 0.0, 4.0}}});
    registry.add("h_jets", ";Number of jets;Count", {HistType::kTH1F, {{1, 0.5, 1.5}}});
    registry.add("h_tracks_per_collision", "Tracks per Collision;Collision Index;Counts", {HistType::kTH1F, {{1, 0.5, 1.5}}});
    registry.add("h_collisionidex", "Collision Index;Collision Index;Counts", {HistType::kTH1F, {{1, 0.5, 1.5}}});

      

    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
  }

  Preslice<aod::JetTracks> TracksPerCollision = aod::jtrack::collisionId;
  Preslice<aod::JetTracksSub> TracksPerCollisionDataSub = aod::bkgcharged::collisionId;
  Preslice<aod::JetParticles> ParticlesPerMcCollision = aod::jmcparticle::mcCollisionId;


  template <bool isMCP, bool isSubtracted, typename T, typename U>
  void jetReclustering(T const& jet, U& splittingTable, double weight)
  {
    // LOGF(info, " Entering jetReclustering " );
    energyMotherVec.clear(); //to be sure its empty before filling
    ptLeadingVec.clear();
    ptSubLeadingVec.clear();
    thetaVec.clear();
    jetReclustered.clear();
    fastjet::ClusterSequenceArea clusterSeq(jetReclusterer.findJets(jetConstituents, jetReclustered)); // jetConstituents is a vector of fastjet::PseudoJet containing (tracks/constituents)
    jetReclustered = sorted_by_pt(jetReclustered); //sort jets by pt decreasing order 
    fastjet::PseudoJet daughterSubJet = jetReclustered[0]; //order 0 of reclustering 
    fastjet::PseudoJet parentSubJet1; //parent1 of daughter
    fastjet::PseudoJet parentSubJet2; //parent2 of daughter
    bool softDropped = false;
    auto nsd = 0.0;
    auto zg = -1.0;
    auto rg = -1.0;

    while (daughterSubJet.has_parents(parentSubJet1, parentSubJet2)) {//while daughter has parents, until we reach the end of reclustering
      if (parentSubJet1.perp() < parentSubJet2.perp()) {
        std::swap(parentSubJet1, parentSubJet2);//for 1 to be greather in pt than 2
      }
      std::vector<int32_t> tracks;
      std::vector<int32_t> candidates;
      std::vector<int32_t> clusters;
      for (const auto& constituent : sorted_by_pt(parentSubJet2.constituents())) {//for boucle over all constituent in parent2
        if (constituent.template user_info<fastjetutilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::track)) {//if constituent is a track of the subjet,
          tracks.push_back(constituent.template user_info<fastjetutilities::fastjet_user_info>().getIndex());//then the index is stored in tracks vector
        }
      }
      splittingTable(jet.globalIndex(), tracks, clusters, candidates, parentSubJet2.perp(), parentSubJet2.eta(), parentSubJet2.phi(), 0);// fill a table with jet devision information : index, tracks, etc
      auto z = parentSubJet2.perp() / (parentSubJet1.perp() + parentSubJet2.perp());
      auto theta = parentSubJet1.delta_R(parentSubJet2);
      energyMotherVec.push_back(daughterSubJet.e());
      ptLeadingVec.push_back(parentSubJet1.pt());
      ptSubLeadingVec.push_back(parentSubJet2.pt());
      thetaVec.push_back(theta);

      if (z >= zCut * TMath::Power(theta / (jet.r() / 100.f), beta)) {
        float thetag = theta / (jet.r() / 100.f);
        // LOGF(info, "Jet radius = %.2f", jet.r() / 100.f);
        if (!softDropped) {//if the splitting hasent been already softdropped softdrop=false
          zg = z;
          rg = theta;
          
          if constexpr (!isSubtracted && !isMCP) {
            // LOGF(info, " Entering if statement for histograms :" );
            registry.fill(HIST("h2_jet_pt_jet_zg"), jet.pt(), zg, weight);
            registry.fill(HIST("h2_jet_pt_jet_rg"), jet.pt(), rg, weight);
            registry.fill(HIST("h2_jet_pt_jet_thetag"), jet.pt(), thetag, weight);
          }
          if constexpr (!isSubtracted && isMCP) {
            registry.fill(HIST("h2_jet_pt_part_jet_zg_part"), jet.pt(), zg, weight);
            registry.fill(HIST("h2_jet_pt_part_jet_rg_part"), jet.pt(), rg, weight);
            registry.fill(HIST("h2_jet_pt_part_jet_thetag_part"), jet.pt(), thetag, weight);
          }
          if constexpr (isSubtracted && !isMCP) {
            // LOGF(info, " Entering if statement for histograms :" );
            registry.fill(HIST("h2_jet_pt_jet_zg_eventwiseconstituentsubtracted"),jet.pt(), zg, weight);
            registry.fill(HIST("h2_jet_pt_jet_rg_eventwiseconstituentsubtracted"),jet.pt(), rg, weight);
            registry.fill(HIST("h2_jet_pt_jet_thetag_eventwiseconstituentsubtracted"),jet.pt(), thetag, weight);
          }
          softDropped = true;//mark the splitting as true to avoid filling it again
        }
        nsd++;//step up the number of iterations
      }
      daughterSubJet = parentSubJet1; //following with the hardest branch
    }
    if constexpr (!isSubtracted && !isMCP) {
      // LOGF(info, " Entering if statement for histograms: " );
      registry.fill(HIST("h2_jet_pt_jet_nsd"), jet.pt(), nsd, weight);
    }
    if constexpr (!isSubtracted && isMCP) {
      registry.fill(HIST("h2_jet_pt_part_jet_nsd_part"), jet.pt(), nsd, weight);
    }
    if constexpr (isSubtracted && !isMCP) {
      registry.fill(HIST("h2_jet_pt_jet_nsd_eventwiseconstituentsubtracted"), jet.pt(), nsd, weight);
    }
  }

  template <bool isSubtracted, typename T, typename U, typename N>
  void analyseCharged(T const& jet, U const& tracksOfCollisions, N& splittingTable, double weight = 1.0) 
  {
    // LOGF(info, " Entering analyseCharged " );
    jetConstituents.clear();
    for (auto& jetConstituent : jet.template tracks_as<U>()) {
      
      fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
    }
    // nSub = jetsubstructureutilities::getNSubjettiness(jet, jet.template tracks_as<U>(), jet.template tracks_as<U>(), jet.template tracks_as<U>(), 2, fastjet::contrib::CA_Axes(), true, zCut, beta);
    jetReclustering<false, isSubtracted>(jet, splittingTable, weight);
  }

  void processDummy(aod::JetTracks const&)
  {
  }
  PROCESS_SWITCH(JetSubstructureTask, processDummy, "Dummy process function turned on by default", true);

  void processCollisions(aod::JetCollisions::iterator const& collision)
  {
    registry.fill(HIST("h_collisions"), 0.5);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collisions"), 2.5);
  }
  PROCESS_SWITCH(JetSpectraCharged, processCollisions, "collisions Data and MCD", true);

  void processChargedJetsData(aod::JetCollisions::iterator const& collision,
                              aod::JetTracks const& tracksOfCollisions,
                              soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets)
  { 
    // LOGF(info, " Entering processChargedJetsData :" );

    ///////////// leading track cut try : (because filter doesnt work)
    bool hasHighPtConstituent = false;
    for (auto& jet : jets){
      for (auto& jetConstituent : jet.tracks_as<aod::JetTracks>()) {
        if (jetConstituent.pt() >= ptLeadingTrackCut) {
          hasHighPtConstituent = true;
          break; // Sortir de la boucle dès qu'un constituant valide est trouvé
        }
      }
      // Si un jet contient un constituant avec un pt > au critère, on l'analyse
      if (hasHighPtConstituent) {
        analyseCharged<false>(jet, tracksOfCollisions, jetSplittingsDataTable);//attention je donne TOUTES les traces pas juste les traces des jets comme dans jetSubstructure.cxx (jetConstituentS,jet.tracks_as<aod::JetTracks>()) en fait on s'en fou car c'est deja filtré dans analysisCharged() c'est emme mieux de tout donner si jamais il y manque
      }
    }
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsData, "charged jet substructure", false);

  void processChargedJetsEventWiseSubData(aod::JetCollisions::iterator const& collision,
                                          soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents> const& jets,
                                          aod::JetTracksSub const& tracksOfCollisions)
  {
    

    bool hasHighPtConstituent = false;
    for (auto& jet : jets){
    registry.fill(HIST("h_jet_pT"), 0.5);
    // auto & jetConstituent0 = jet.tracks_as<aod::JetTracksSub>().iteratorAt(0)
      for (auto& jetConstituent : jet.tracks_as<aod::JetTracksSub>()) {
        if (jetConstituent.pt() >= ptLeadingTrackCut) {
          LOGF(info, "Jet with leading constituent pt = %.2f found", jetConstituent.pt());
          // LOGF(info, "test1 ");
          hasHighPtConstituent = true;
          break; // Sortir de la boucle dès qu'un constituant valide est trouvé
        }
      }
      // Si un jet contient un constituant avec un pt > au critère, on l'analyse
      if (hasHighPtConstituent) {
        // LOGF(info, "test2 ");
        registry.fill(HIST("h_jet_pT"), 1.5);
        analyseCharged<true>(jet, tracksOfCollisions, jetSplittingsDataSubTable);
        registry.fill(HIST("h_jet_pT"), 2.5);
        
      }
    }
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsEventWiseSubData, "eventwise-constituent subtracted charged jet substructure", false);

    void processChargedJetsMCD(aod::JetCollisionsMCD::iterator const& collision,
                               aod::JetMcCollisions const&, //join the weight
                               soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                               aod::JetTracks const& tracks)
  { 
    ///////////// leading track cut /////////////
    bool hasHighPtConstituent = false;
    for (auto& jet : jets){
      for (auto& jetConstituent : jet.tracks_as<aod::JetTracks>()) {
        if (jetConstituent.pt() >= ptLeadingTrackCut) {
          hasHighPtConstituent = true;
          break; // Sortir de la boucle dès qu'un constituant valide est trouvé
        }
      }
      // Si un jet contient un constituant avec un pt > au critère, on l'analyse
      if (hasHighPtConstituent) {
        analyseCharged<false>(jet, tracks, jetSplittingsMCDTable, collision.mcCollision().weight());
        LOGF(info, "processChargedJetsMCD: weight = %.4f", collision.mcCollision().weight());
      }
    }
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsMCD, "charged jet substructure", false);

    void processChargedJetsEventWiseSubMCD(aod::JetCollisionsMCD::iterator const& collision,
                                           aod::JetMcCollisions const&, //join the weight
                                           soa::Join<aod::ChargedMCDetectorLevelEventWiseSubtractedJets, aod::ChargedMCDetectorLevelEventWiseSubtractedJetConstituents> const& jets,
                                           aod::JetTracksSub const& tracks)
  { 
    ///////////// leading track cut /////////////
    bool hasHighPtConstituent = false;
    for (auto& jet : jets){
      for (auto& jetConstituent : jet.tracks_as<aod::JetTracksSub>()) {
        if (jetConstituent.pt() >= ptLeadingTrackCut) {
          hasHighPtConstituent = true;
          break; // Sortir de la boucle dès qu'un constituant valide est trouvé
        }
      }
      // Si un jet contient un constituant avec un pt > au critère, on l'analyse
      if (hasHighPtConstituent) {
        analyseCharged<true>(jet, tracks, jetSplittingsDataSubTable, collision.mcCollision().weight());
        LOGF(info, "processChargedJetsEventWiseSubMCD: weight = %.4f", collision.mcCollision().weight());
      }
    }
  }
  
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsEventWiseSubMCD, "eventwise-constituent subtracted MCD charged jet substructure", false);

  // void processChargedJetsMCP(aod::JetMcCollisions::iterator const& mcCollision,
  //                            typename soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
  //                            aod::JetParticles const& particles)
  // {
  //   for (auto& jet : jets){
  //     jetConstituents.clear();
  //     for (auto& jetConstituent : jet.template tracks_as<aod::JetParticles>()) {
  //       fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex(), static_cast<int>(JetConstituentStatus::track), pdg->Mass(jetConstituent.pdgCode()));
  //     }
  //     jetReclustering<true, false>(jet, jetSplittingsMCPTable);
  //   }
  // }
  // PROCESS_SWITCH(JetSubstructureTask, processChargedJetsMCP, "charged jet substructure on MC particle level", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  return WorkflowSpec{adaptAnalysisTask<JetSubstructureTask>(
    cfgc, TaskName{"jet-substructure-softdrop"})};
}