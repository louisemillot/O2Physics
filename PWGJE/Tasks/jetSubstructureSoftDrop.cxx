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
#include "PWGJE/Core/JetFindingUtilities.h"



#include "EventFiltering/filterTables.h"


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct JetSubstructureTask {

  using ChargedMCDMatchedJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>;
  using ChargedMCPMatchedJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>;
  using ChargedMCDMatchedJetsWeighted = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetEventWeights>;
  using ChargedMCPMatchedJetsWeighted = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetEventWeights>;

  Produces<aod::ChargedSPs> jetSplittingsDataTable;
  Produces<aod::ChargedMCDetectorLevelSPs> jetSplittingsMCDTable;
  Produces<aod::ChargedMCParticleLevelSPs> jetSplittingsMCPTable;
  Produces<aod::ChargedEventWiseSubtractedSPs> jetSplittingsDataSubTable;
  Produces<aod::ChargedMCDetectorLevelEventWiseSubtractedSPs> jetSplittingsMCDSubTable;
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
  Configurable<float> ptLeadingTrackCut{"ptLeadingTrackCut", 5.0f, "Leading track cut : minimum pT selection on jet constituent"};
  Configurable<float> ptLeadingTrackCutMax{"ptleadingConstituentPtMax", 9999.0, "maximum pT selection on jet constituent"};
  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};
  Configurable<float> pTHatAbsoluteMin{"pTHatAbsoluteMin", -99.0, "minimum value of pTHat"};
  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum track occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<int> trackOccupancyInTimeRangeMin{"trackOccupancyInTimeRangeMin", -999999, "minimum track occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false, "flag to choose to reject min. bias gap events; jet-level rejection can also be applied at the jet finder level for jets only, here rejection is applied for collision and track process functions for the first time, and on jets in case it was set to false at the jet finder level"};
  Configurable<float> jetEtaMin{"jetEtaMin", -0.7, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 0.7, "maximum jet pseudorapidity"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<bool> checkLeadConstituentPtForMcpJets{"checkLeadConstituentPtForMcpJets", false, "flag to choose whether particle level jets should have their lead track pt above leadingConstituentPtMin to be accepted; off by default, as leadingConstituentPtMin cut is only applied on MCD jets for the Pb-Pb analysis using pp MC anchored to Pb-Pb for the response matrix"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<bool> checkGeoMatched{"checkGeoMatched", true, "0: turn off geometry matching, 1: do geometry matching "};
  Configurable<bool> checkPtMatched{"checkPtMatched", false, "0: turn off pT matching, 1: do pT matching"};
  Configurable<bool> checkGeoPtMatched{"checkGeoPtMatched", false, "0: turn off geometry and pT matching, 1: do geometry and pT matching"};





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
    AxisSpec jetPtAxis = {200, 0., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec centralityAxis = {1200, -10., 110., "Centrality"};


    //Reclustering
    registry.add("h_leading_prong_pt", "Leading prong p_{T};p_{T} (GeV/c);Entries", {HistType::kTH1F, {jetPtAxis}});
    registry.add("h_subleading_prong_pt", "Subleading prong p_{T};p_{T} (GeV/c);Entries", {HistType::kTH1F, {jetPtAxis}});

    if (doprocessChargedJetsData || doprocessChargedJetsMCD || doprocessChargedJetsMCDWeighted) {
      registry.add("h2_jet_pt_jet_zg", ";#it{p}_{T,jet} (GeV/#it{c});#it{z}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
      registry.add("h2_jet_pt_jet_rg", ";#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
      registry.add("h2_jet_pt_jet_thetag", ";#it{p}_{T,jet} (GeV/#it{c});#it{theta}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
      registry.add("h2_jet_pt_jet_nsd", ";#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});
      registry.add("h_jet_thetag", ";#it{theta}_{g};Entries", {HistType::kTH1F, {{22, 0.0, 1.1}}});
      registry.add("h_jet_zg", ";#it{z}_{g};Entries", {HistType::kTH1F, {{22, 0.0, 1.1}}});
    }

    if (doprocessChargedJetsMCP || doprocessChargedJetsMCPWeighted) {
      registry.add("h2_jet_pt_part_jet_zg_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{z}_{g}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
      registry.add("h2_jet_pt_part_jet_rg_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{R}_{g}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
      registry.add("h2_jet_pt_part_jet_thetag_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{theta}_{g}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
      registry.add("h2_jet_pt_part_jet_nsd_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{n}_{SD}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});
      registry.add("h_jet_thetag_MCP", ";#it{theta}_{g};Entries", {HistType::kTH1F, {{22, 0.0, 1.1}}});
      registry.add("h_jet_zg_MCP", ";#it{z}_{g};Entries", {HistType::kTH1F, {{22, 0.0, 1.1}}});
    }

    if (doprocessChargedJetsEventWiseSubData || doprocessChargedJetsEventWiseSubMCD || doprocessChargedJetsEventWiseSubMCDWeighted) { 
      registry.add("h2_jet_pt_jet_zg_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{z}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
      registry.add("h2_jet_pt_jet_rg_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
      registry.add("h2_jet_pt_jet_thetag_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{theta}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
      registry.add("h2_jet_pt_jet_nsd_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});
      registry.add("h_jet_thetag_eventwiseconstituentsubtracted", ";#it{theta}_{g};Entries", {HistType::kTH1F, {{22, 0.0, 1.1}}});
      registry.add("h_jet_zg_eventwiseconstituentsubtracted", ";#it{z}_{g};Entries", {HistType::kTH1F, {{22, 0.0, 1.1}}});
    }

    //Collisions
    if (doprocessCollisions || doprocessCollisionsWeighted) {
      registry.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h_fakecollisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h2_centrality_occupancy", "centrality vs occupancy; centrality; occupancy", {HistType::kTH2F, {centralityAxis, {60, 0, 30000}}});
      registry.add("h_collisions_Zvertex", "position of collision ;#it{Z} (cm)", {HistType::kTH1F, {{300, -15.0, 15.0}}});
      // registry.add("h_jets", ";Number of jets;Count", {HistType::kTH1F, {{1, 0.5, 1.5}}});
      // registry.add("h_tracks_per_collision", "Tracks per Collision;Collision Index;Counts", {HistType::kTH1F, {{1, 0.5, 1.5}}});
      // registry.add("h_collisionidex", "Collision Index;Collision Index;Counts", {HistType::kTH1F, {{1, 0.5, 1.5}}});
      if (doprocessCollisionsWeighted) {
        registry.add("h_collisions_weighted", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      }
    }
    //McCollisions
    if (doprocessMcCollisions || doprocessMcCollisionsWeighted) {
      registry.add("h_mcColl_counts", " number of mc events; event status; entries", {HistType::kTH1F, {{10, 0, 10}}});
      registry.get<TH1>(HIST("h_mcColl_counts"))->GetXaxis()->SetBinLabel(1, "allMcColl");
      registry.get<TH1>(HIST("h_mcColl_counts"))->GetXaxis()->SetBinLabel(2, "vertexZ");
      registry.get<TH1>(HIST("h_mcColl_counts"))->GetXaxis()->SetBinLabel(3, "noRecoColl");
      registry.get<TH1>(HIST("h_mcColl_counts"))->GetXaxis()->SetBinLabel(4, "recoEvtSel");
      registry.get<TH1>(HIST("h_mcColl_counts"))->GetXaxis()->SetBinLabel(5, "centralitycut");
      registry.get<TH1>(HIST("h_mcColl_counts"))->GetXaxis()->SetBinLabel(6, "occupancycut");
      registry.add("h_mc_zvertex", "position of collision ;#it{Z} (cm)", {HistType::kTH1F, {{300, -15.0, 15.0}}});
      if (doprocessMcCollisionsWeighted) {
        registry.add("h_mcColl_counts_weight", " number of mc events; event status; entries", {HistType::kTH1F, {{10, 0, 10}}});
        registry.get<TH1>(HIST("h_mcColl_counts_weight"))->GetXaxis()->SetBinLabel(1, "allMcColl");
        registry.get<TH1>(HIST("h_mcColl_counts_weight"))->GetXaxis()->SetBinLabel(2, "vertexZ");
        registry.get<TH1>(HIST("h_mcColl_counts_weight"))->GetXaxis()->SetBinLabel(3, "noRecoColl");
        registry.get<TH1>(HIST("h_mcColl_counts_weight"))->GetXaxis()->SetBinLabel(4, "recoEvtSel");
        registry.get<TH1>(HIST("h_mcColl_counts_weight"))->GetXaxis()->SetBinLabel(5, "centralitycut");
        registry.get<TH1>(HIST("h_mcColl_counts_weight"))->GetXaxis()->SetBinLabel(6, "occupancycut");
        registry.add("h_mc_zvertex_weight", "position of collision ;#it{Z} (cm)", {HistType::kTH1F, {{300, -15.0, 15.0}}});
      }
    }

    //data
    if (doprocessChargedJetsData) {
      registry.add("h_jet_pt_initial_data", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      registry.add("h_jet_pt_after_leadingtrackcut_data", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      // registry.add("h_jet_pt_after_grooming", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
    }

    //data + eventwise
    if (doprocessChargedJetsEventWiseSubData) {
      registry.add("h_jet_pt_initial_data_eventwise", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      registry.add("h_jet_pt_after_leadingtrackcut_data_eventwise", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      // registry.add("h_jet_pt_after_grooming", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
    }

    //MCD 
    if (doprocessChargedJetsMCD || doprocessChargedJetsMCDWeighted ) {
      registry.add("h_jet_pt_initial_mcd", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      registry.add("h_jet_pt_after_leadingtrackcut_mcd", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      //MCD + weighted
      if (doprocessChargedJetsMCDWeighted) {
        registry.add("h_jet_pthat_initial_mcd", "jet pT hat;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
        registry.add("h_jet_pthat_initial_mcd_weighted", "jet pT hat;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
        registry.add("h_jet_pt_after_leadingtrackcut_mcd_weighted", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
        registry.add("h_jet_pt_initial_mcd_weighted", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      }
    }

    //MCD + eventwise
    if (doprocessChargedJetsEventWiseSubMCD || doprocessChargedJetsEventWiseSubMCDWeighted ) {
      registry.add("h_jet_pt_initial_mcd_eventwise", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      registry.add("h_jet_pt_after_leadingtrackcut_mcd_eventwise", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      //MCD + eventwise + weighted
      if (doprocessChargedJetsEventWiseSubMCDWeighted) {
        registry.add("h_jet_pthat_initial_mcd_eventwise", "jet pT hat;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
        registry.add("h_jet_pthat_initial_mcd_eventwise_weighted", "jet pT hat;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
        registry.add("h_jet_pt_initial_mcd_eventwise_weighted", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
        registry.add("h_jet_pt_after_leadingtrackcut_mcd_eventwise_weighted", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      }
    }

    //MCP 
    if (doprocessChargedJetsMCP || doprocessChargedJetsMCPWeighted) {
      registry.add("h_jet_pt_initial_mcp", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      registry.add("h_jet_pt_after_leadingtrackcut_mcp", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      //MCP + weighted 
      if (doprocessChargedJetsMCPWeighted) {
        registry.add("h_jet_phat_initial_mcp", "jet pT hat;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
        registry.add("h_jet_pt_initial_mcp_weighted", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
        registry.add("h_jet_pthat_initial_mcp_weighted", "jet pT hat;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
        registry.add("h_jet_pt_after_leadingtrackcut_mcp_weighted", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      }
    }

    if (doprocessJetsMatched || doprocessJetsMatchedWeighted) {
      if (checkGeoMatched) {
        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint", "pT mcd vs. pT mcp;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcpetaconstraint", "pT mcd vs. pT mcp;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
        registry.add("h2_jet_eta_mcd_jet_eta_mcp_matchedgeo", "Eta mcd vs. Eta mcp;#eta_{jet}^{mcd};#eta_{jet}^{mcp}", {HistType::kTH2F, {jetEtaAxis, jetEtaAxis}});
        registry.add("h2_jet_phi_mcd_jet_phi_mcp_matchedgeo_mcdetaconstraint", "Phi mcd vs. Phi mcp;#varphi_{jet}^{mcd};#varphi_{jet}^{mcp}", {HistType::kTH2F, {phiAxis, phiAxis}});
        registry.add("h2_jet_phi_mcd_jet_phi_mcp_matchedgeo_mcpetaconstraint", "Phi mcd vs. Phi mcp;#varphi_{jet}^{mcd};#varphi_{jet}^{mcp}", {HistType::kTH2F, {phiAxis, phiAxis}});
        registry.add("h2_jet_ntracks_mcd_jet_ntracks_mcp_matchedgeo", "Ntracks mcd vs. Ntracks mcp;N_{jet tracks}^{mcd};N_{jet tracks}^{mcp}", {HistType::kTH2F, {{200, -0.5, 199.5}, {200, -0.5, 199.5}}});
        registry.add("h2_jet_pt_mcp_jet_pt_diff_matchedgeo", "jet mcp pT vs. delta pT / jet mcp pt;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); (#it{p}_{T,jet}^{mcp} (GeV/#it{c}) - #it{p}_{T,jet}^{mcd} (GeV/#it{c})) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_mcd_jet_pt_diff_matchedgeo", "jet mcd pT vs. delta pT / jet mcd pt;#it{p}_{T,jet}^{mcd} (GeV/#it{c}); (#it{p}_{T,jet}^{mcd} (GeV/#it{c}) - #it{p}_{T,jet}^{mcp} (GeV/#it{c})) / #it{p}_{T,jet}^{mcd} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_mcp_jet_pt_ratio_matchedgeo", "jet mcp pT vs. jet mcd pT / jet mcp pt;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); #it{p}_{T,jet}^{mcd} (GeV/#it{c}) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 5.0}}});
      }
      if (checkPtMatched) {
        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedpt_mcdetaconstraint", "pT mcd vs. pT mcp;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedpt_mcpetaconstraint", "pT mcd vs. pT mcp;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
        registry.add("h2_jet_eta_mcd_jet_eta_mcp_matchedpt", "Eta mcd vs. Eta mcp;#eta_{jet}^{mcd};#eta_{jet}^{mcp}", {HistType::kTH2F, {jetEtaAxis, jetEtaAxis}});
        registry.add("h2_jet_phi_mcd_jet_phi_mcp_matchedgpt_mcdetaconstraint", "Phi mcd vs. Phi mcp;#varphi_{jet}^{mcd};#varphi_{jet}^{mcp}", {HistType::kTH2F, {phiAxis, phiAxis}});
        registry.add("h2_jet_phi_mcd_jet_phi_mcp_matchedgpt_mcpetaconstraint", "Phi mcd vs. Phi mcp;#varphi_{jet}^{mcd};#varphi_{jet}^{mcp}", {HistType::kTH2F, {phiAxis, phiAxis}});
        registry.add("h2_jet_ntracks_mcd_jet_ntracks_mcp_matchedpt", "Ntracks mcd vs. Ntracks mcp;N_{jet tracks}^{mcd};N_{jet tracks}^{mcp}", {HistType::kTH2F, {{200, -0.5, 199.5}, {200, -0.5, 199.5}}});
        registry.add("h2_jet_pt_mcp_jet_pt_diff_matchedpt", "jet mcp pT vs. delta pT / jet mcp pt;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); (#it{p}_{T,jet}^{mcp} (GeV/#it{c}) - #it{p}_{T,jet}^{mcd} (GeV/#it{c})) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_mcd_jet_pt_diff_matchedpt", "jet mcd pT vs. delta pT / jet mcd pt;#it{p}_{T,jet}^{mcd} (GeV/#it{c}); (#it{p}_{T,jet}^{mcd} (GeV/#it{c}) - #it{p}_{T,jet}^{mcp} (GeV/#it{c})) / #it{p}_{T,jet}^{mcd} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_mcp_jet_pt_ratio_matchedpt", "jet mcp pT vs. jet mcd pT / jet mcp pt;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); #it{p}_{T,jet}^{mcd} (GeV/#it{c}) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 5.0}}});
      }
      if (checkGeoPtMatched) {
        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeopt_mcdetaconstraint", "pT mcd vs. pT mcp;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeopt_mcpetaconstraint", "pT mcd vs. pT mcp;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
        registry.add("h2_jet_eta_mcd_jet_eta_mcp_matchedgeopt", "Eta mcd vs. Eta mcp;#eta_{jet}^{mcd};#eta_{jet}^{mcp}", {HistType::kTH2F, {jetEtaAxis, jetEtaAxis}});
        registry.add("h2_jet_phi_mcd_jet_phi_mcp_matchedgeopt_mcdetaconstraint", "Phi mcd vs. Phi mcp;#varphi_{jet}^{mcd};#varphi_{jet}^{mcp}", {HistType::kTH2F, {phiAxis, phiAxis}});
        registry.add("h2_jet_phi_mcd_jet_phi_mcp_matchedgeopt_mcpetaconstraint", "Phi mcd vs. Phi mcp;#varphi_{jet}^{mcd};#varphi_{jet}^{mcp}", {HistType::kTH2F, {phiAxis, phiAxis}});
        registry.add("h2_jet_ntracks_mcd_jet_ntracks_mcp_matchedgeopt", "Ntracks mcd vs. Ntracks mcp;N_{jet tracks}^{mcd};N_{jet tracks}^{mcp}", {HistType::kTH2F, {{200, -0.5, 199.5}, {200, -0.5, 199.5}}});
        registry.add("h2_jet_pt_mcp_jet_pt_diff_matchedgeopt", "jet mcp pT vs. delta pT / jet mcp pt;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); (#it{p}_{T,jet}^{mcp} (GeV/#it{c}) - #it{p}_{T,jet}^{mcd} (GeV/#it{c})) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_mcd_jet_pt_diff_matchedgeopt", "jet mcd pT vs. delta pT / jet mcd pt;#it{p}_{T,jet}^{mcd} (GeV/#it{c}); (#it{p}_{T,jet}^{mcd} (GeV/#it{c}) - #it{p}_{T,jet}^{mcp} (GeV/#it{c})) / #it{p}_{T,jet}^{mcd} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_mcp_jet_pt_ratio_matchedgeopt", "jet mcp pT vs. jet mcd pT / jet mcp pt;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); #it{p}_{T,jet}^{mcd} (GeV/#it{c}) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 5.0}}});
      }
    }


    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    
  }

  // Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter collisionFilter = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centFT0M >= centralityMin && aod::jcollision::centFT0M < centralityMax);
  Filter mcCollisionFilter = nabs(aod::jmccollision::posZ)< vertexZCut;

  Preslice<aod::JetTracks> TracksPerCollision = aod::jtrack::collisionId;
  Preslice<aod::JetTracksSub> TracksPerCollisionDataSub = aod::bkgcharged::collisionId;
  Preslice<aod::JetParticles> ParticlesPerMcCollision = aod::jmcparticle::mcCollisionId;
  Preslice<ChargedMCDMatchedJets> mcdjetsPerJCollision = o2::aod::jet::collisionId;

  template <typename TTracks, typename TJets>
  bool isAcceptedJet(TJets const& jet, bool mcLevelIsParticleLevel = false)
  {
    if (jetAreaFractionMin > -98.0) {
      if (jet.area() < jetAreaFractionMin * o2::constants::math::PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    bool checkConstituentPt = true;
    bool checkConstituentMinPt = (ptLeadingTrackCut > -98.0);
    bool checkConstituentMaxPt = (ptLeadingTrackCutMax < 9998.0);
    if (!checkConstituentMinPt && !checkConstituentMaxPt) {
      checkConstituentPt = false;
    }
    if (mcLevelIsParticleLevel && !checkLeadConstituentPtForMcpJets) {
      checkConstituentPt = false;
    }

    if (checkConstituentPt) {
      bool isMinLeadingConstituent = !checkConstituentMinPt;
      bool isMaxLeadingConstituent = true;

      for (const auto& constituent : jet.template tracks_as<TTracks>()) {
        double pt = constituent.pt();

        if (checkConstituentMinPt && pt >= ptLeadingTrackCut) {
          isMinLeadingConstituent = true;
        }
        if (checkConstituentMaxPt && pt > ptLeadingTrackCutMax) {
          isMaxLeadingConstituent = false;
        }
      }
      return isMinLeadingConstituent && isMaxLeadingConstituent;
    }
    return true;
  }

  template <typename TBase, typename TTag>
  void fillMatchedHistograms(TBase const& jetMCD, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jetMCD.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin) {
      return;
    }
    // fill geometry matched histograms
    if (checkGeoMatched) {
      if (jetMCD.has_matchedJetGeo()) {
        for (const auto& jetMCP : jetMCD.template matchedJetGeo_as<std::decay_t<TTag>>()) {
          if (jetMCP.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin) {
            continue;
          }
          if (jetMCD.r() == round(selectedJetsRadius * 100.0f)) {
            double dpt = jetMCP.pt() - jetMCD.pt();
            if (jetfindingutilities::isInEtaAcceptance(jetMCD, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint"), jetMCD.pt(), jetMCP.pt(), weight);
              registry.fill(HIST("h2_jet_phi_mcd_jet_phi_mcp_matchedgeo_mcdetaconstraint"), jetMCD.phi(), jetMCP.phi(), weight);
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_diff_matchedgeo"), jetMCD.pt(), dpt / jetMCD.pt(), weight);
              registry.fill(HIST("h2_jet_ntracks_mcd_jet_ntracks_mcp_matchedgeo"), jetMCD.tracksIds().size(), jetMCP.tracksIds().size(), weight);
            }
            if (jetfindingutilities::isInEtaAcceptance(jetMCP, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcpetaconstraint"), jetMCD.pt(), jetMCP.pt(), weight);
              registry.fill(HIST("h2_jet_phi_mcd_jet_phi_mcp_matchedgeo_mcpetaconstraint"), jetMCD.phi(), jetMCP.phi(), weight);
              registry.fill(HIST("h2_jet_pt_mcp_jet_pt_diff_matchedgeo"), jetMCP.pt(), dpt / jetMCP.pt(), weight);
              registry.fill(HIST("h2_jet_pt_mcp_jet_pt_ratio_matchedgeo"), jetMCP.pt(), jetMCD.pt() / jetMCP.pt(), weight);
            }
            registry.fill(HIST("h2_jet_eta_mcd_jet_eta_mcp_matchedgeo"), jetMCD.eta(), jetMCP.eta(), weight);
          }
        }
      }
    }
    // fill pt matched histograms
    if (checkPtMatched) {
      if (jetMCD.has_matchedJetPt()) {
        for (const auto& jetMCP : jetMCD.template matchedJetPt_as<std::decay_t<TTag>>()) {
          if (jetMCP.pt() > pTHatMaxMCP * pTHat) {
            continue;
          }
          if (!jetfindingutilities::isInEtaAcceptance(jetMCP, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
            continue;
          }
          if (jetMCD.r() == round(selectedJetsRadius * 100.0f)) {
            double dpt = jetMCP.pt() - jetMCD.pt();
            if (jetfindingutilities::isInEtaAcceptance(jetMCD, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedpt_mcdetaconstraint"), jetMCD.pt(), jetMCP.pt(), weight);
              registry.fill(HIST("h2_jet_phi_mcd_jet_phi_mcp_matchedpt_mcdetaconstraint"), jetMCD.phi(), jetMCP.phi(), weight);
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_diff_matchedpt"), jetMCD.pt(), dpt / jetMCD.pt(), weight);
              registry.fill(HIST("h2_jet_ntracks_mcd_jet_ntracks_mcp_matchedpt"), jetMCD.tracksIds().size(), jetMCP.tracksIds().size(), weight);
            }
            if (jetfindingutilities::isInEtaAcceptance(jetMCP, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedpt_mcpetaconstraint"), jetMCD.pt(), jetMCP.pt(), weight);
              registry.fill(HIST("h2_jet_phi_mcd_jet_phi_mcp_matchedpt_mcpetaconstraint"), jetMCD.phi(), jetMCP.phi(), weight);
              registry.fill(HIST("h2_jet_pt_mcp_jet_pt_diff_matchedpt"), jetMCP.pt(), dpt / jetMCP.pt(), weight);
              registry.fill(HIST("h2_jet_pt_mcp_jet_pt_ratio_matchedpt"), jetMCP.pt(), jetMCD.pt() / jetMCP.pt(), weight);
            }
            registry.fill(HIST("h2_jet_eta_mcd_jet_eta_mcp_matchedpt"), jetMCD.eta(), jetMCP.eta(), weight);
          }
        }
      }
    }
    // fill geometry and pt histograms
    if (checkGeoPtMatched) {
      if (jetMCD.has_matchedJetGeo() && jetMCD.has_matchedJetPt()) {
        for (const auto& jetMCP : jetMCD.template matchedJetGeo_as<std::decay_t<TTag>>()) {
          if (jetMCP.pt() > pTHatMaxMCP * pTHat) {
            continue;
          }
          if (!jetfindingutilities::isInEtaAcceptance(jetMCP, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
            continue;
          }
          if (jetMCD.template matchedJetGeo_first_as<std::decay_t<TTag>>().globalIndex() == jetMCD.template matchedJetPt_first_as<std::decay_t<TTag>>().globalIndex()) { // not a good way to do this
            double dpt = jetMCP.pt() - jetMCD.pt();
            if (jetfindingutilities::isInEtaAcceptance(jetMCD, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeopt_mcdetaconstraint"), jetMCD.pt(), jetMCP.pt(), weight);
              registry.fill(HIST("h2_jet_phi_mcd_jet_phi_mcp_matchedgeopt_mcdetaconstraint"), jetMCD.phi(), jetMCP.phi(), weight);
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_diff_matchedgeopt"), jetMCD.pt(), dpt / jetMCD.pt(), weight);
              registry.fill(HIST("h2_jet_ntracks_mcd_jet_ntracks_mcp_matchedgeopt"), jetMCD.tracksIds().size(), jetMCP.tracksIds().size(), weight);
            }
            if (jetfindingutilities::isInEtaAcceptance(jetMCP, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
              registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeopt_mcpetaconstraint"), jetMCD.pt(), jetMCP.pt(), weight);
              registry.fill(HIST("h2_jet_phi_mcd_jet_phi_mcp_matchedgeopt_mcpetaconstraint"), jetMCD.phi(), jetMCP.phi(), weight);
              registry.fill(HIST("h2_jet_pt_mcp_jet_pt_diff_matchedgeopt"), jetMCP.pt(), dpt / jetMCP.pt(), weight);
              registry.fill(HIST("h2_jet_pt_mcp_jet_pt_ratio_matchedgeopt"), jetMCP.pt(), jetMCD.pt() / jetMCP.pt(), weight);
            }
            registry.fill(HIST("h2_jet_eta_mcd_jet_eta_mcp_matchedpt"), jetMCD.eta(), jetMCP.eta(), weight);
          }
        }
      }
    }
  }

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
        std::swap(parentSubJet1, parentSubJet2);//for 1 to be greather in pt than 2 (1: leading, 2:subleading)
      }
      registry.fill(HIST("h_leading_prong_pt"), parentSubJet1.pt(), weight);
      registry.fill(HIST("h_subleading_prong_pt"), parentSubJet2.pt(), weight);
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
          
          if constexpr (!isSubtracted && !isMCP) { //data & MCD no sub
            // LOGF(info, " Entering if statement for histograms :" );
            registry.fill(HIST("h2_jet_pt_jet_zg"), jet.pt(), zg, weight);
            registry.fill(HIST("h2_jet_pt_jet_rg"), jet.pt(), rg, weight);
            registry.fill(HIST("h2_jet_pt_jet_thetag"), jet.pt(), thetag, weight);
            registry.fill(HIST("h_jet_thetag"), thetag, weight);
            registry.fill(HIST("h_jet_zg"), zg, weight);
          }
          if constexpr (!isSubtracted && isMCP) { //MCP only no sub
            registry.fill(HIST("h2_jet_pt_part_jet_zg_part"), jet.pt(), zg, weight);
            registry.fill(HIST("h2_jet_pt_part_jet_rg_part"), jet.pt(), rg, weight);
            registry.fill(HIST("h2_jet_pt_part_jet_thetag_part"), jet.pt(), thetag, weight);
            registry.fill(HIST("h_jet_thetag_MCP"), thetag, weight);
            registry.fill(HIST("h_jet_zg_MCP"), zg, weight);
          }
          if constexpr (isSubtracted && !isMCP) { //data & MDD sub
            // LOGF(info, " Entering if statement for histograms :" );
            registry.fill(HIST("h2_jet_pt_jet_zg_eventwiseconstituentsubtracted"),jet.pt(), zg, weight);
            registry.fill(HIST("h2_jet_pt_jet_rg_eventwiseconstituentsubtracted"),jet.pt(), rg, weight);
            registry.fill(HIST("h2_jet_pt_jet_thetag_eventwiseconstituentsubtracted"),jet.pt(), thetag, weight);
            registry.fill(HIST("h_jet_thetag_eventwiseconstituentsubtracted"), thetag, weight);
            registry.fill(HIST("h_jet_zg_eventwiseconstituentsubtracted"), zg, weight);
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

  void processCollisions(soa::Filtered<aod::JetCollisions>::iterator const& collision)
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
    registry.fill(HIST("h2_centrality_occupancy"), collision.centFT0M(), collision.trackOccupancyInTimeRange());
    registry.fill(HIST("h_collisions_Zvertex"), collision.posZ());
  }
  PROCESS_SWITCH(JetSubstructureTask, processCollisions, "collisions Data and MCD", true);

  void processCollisionsWeighted(soa::Join<aod::JetCollisions, aod::JMcCollisionLbs>::iterator const& collision,
    aod::JetMcCollisions const&)
  {
  if (!collision.has_mcCollision()) {
    registry.fill(HIST("h_fakecollisions"), 0.5);
  }
  float eventWeight = collision.weight();
  registry.fill(HIST("h_collisions"), 0.5);
  registry.fill(HIST("h_collisions_weighted"), 0.5, eventWeight);
  if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
    return;
  }
  registry.fill(HIST("h_collisions"), 1.5);
  registry.fill(HIST("h_collisions_weighted"), 1.5, eventWeight);
  if (std::abs(collision.posZ()) > vertexZCut) {
    return;
  }
  registry.fill(HIST("h_collisions"), 2.5);
  registry.fill(HIST("h_collisions_weighted"), 2.5, eventWeight);
  registry.fill(HIST("h2_centrality_occupancy"), collision.centFT0M(), collision.trackOccupancyInTimeRange());
  registry.fill(HIST("h_collisions_Zvertex"), collision.posZ(), eventWeight);
  }
  PROCESS_SWITCH(JetSubstructureTask, processCollisionsWeighted, "weighted collsions for Data and MCD", false);

  void processMcCollisions(soa::Filtered<aod::JetMcCollisions>::iterator const& mcCollision,
                           soa::SmallGroups<aod::JetCollisionsMCD> const& collisions)
{ 
  float eventWeight = mcCollision.weight();
  registry.fill(HIST("h_mcColl_counts"), 0.5);

  if (std::abs(mcCollision.posZ()) > vertexZCut) {
    return;
  }
  registry.fill(HIST("h_mcColl_counts"), 1.5);

  if (collisions.size() < 1) {
    return;
  }
  registry.fill(HIST("h_mcColl_counts"), 2.5);

  bool hasSel8Coll = false;
  bool centralityIsGood = false;
  bool occupancyIsGood = false;
  for (auto const& collision : collisions) {
    if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
      hasSel8Coll = true;
    }
    if ((centralityMin < collision.centFT0M()) && (collision.centFT0M() < centralityMax)) {
      centralityIsGood = true;
    }
    if ((trackOccupancyInTimeRangeMin < collision.trackOccupancyInTimeRange()) && (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMax)) {
      occupancyIsGood = true;
    }
  }
  if (!hasSel8Coll) {
    return;
  }
  registry.fill(HIST("h_mcColl_counts"), 3.5);

  if (!centralityIsGood) {
    return;
  }
  registry.fill(HIST("h_mcColl_counts"), 4.5);

  if (!occupancyIsGood) {
    return;
  }
  registry.fill(HIST("h_mcColl_counts"), 5.5);
  registry.fill(HIST("h_mc_zvertex"), mcCollision.posZ());
}
PROCESS_SWITCH(JetSubstructureTask, processMcCollisions, "Mc collisions ", false);

  void processMcCollisionsWeighted(soa::Filtered<aod::JetMcCollisions>::iterator const& mcCollision,
                           soa::SmallGroups<aod::JetCollisionsMCD> const& collisions)
  { 
    float eventWeight = mcCollision.weight();
    registry.fill(HIST("h_mcColl_counts"), 0.5);
    registry.fill(HIST("h_mcColl_counts_weight"), 0.5, eventWeight);
    
    if (std::abs(mcCollision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("h_mcColl_counts"), 1.5);
    registry.fill(HIST("h_mcColl_counts_weight"), 1.5, eventWeight);
    
    if (collisions.size() < 1) {
      return;
    }
    registry.fill(HIST("h_mcColl_counts"), 2.5);
    registry.fill(HIST("h_mcColl_counts_weight"), 2.5, eventWeight);
  
    bool hasSel8Coll = false;
    bool centralityIsGood = false;
    bool occupancyIsGood = false;
    for (auto const& collision : collisions) {
      if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
        hasSel8Coll = true;
      }
      if ((centralityMin < collision.centFT0M()) && (collision.centFT0M() < centralityMax)) {
        centralityIsGood = true;
      }
      if ((trackOccupancyInTimeRangeMin < collision.trackOccupancyInTimeRange()) && (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMax)) {
        occupancyIsGood = true;
      }
    }
    if (!hasSel8Coll) {
      return;
    }
    registry.fill(HIST("h_mcColl_counts"), 3.5);
    registry.fill(HIST("h_mcColl_counts_weight"), 3.5, eventWeight);
    
    if (!centralityIsGood) {
      return;
    }
    registry.fill(HIST("h_mcColl_counts"), 4.5);
    registry.fill(HIST("h_mcColl_counts_weight"), 4.5, eventWeight);

    if (!occupancyIsGood) {
      return;
    }
    registry.fill(HIST("h_mcColl_counts"), 5.5);
    registry.fill(HIST("h_mcColl_counts_weight"), 5.5, eventWeight);
    registry.fill(HIST("h_mc_zvertex"), mcCollision.posZ());
    registry.fill(HIST("h_mc_zvertex_weight"), mcCollision.posZ(), eventWeight);
  }
  PROCESS_SWITCH(JetSubstructureTask, processMcCollisionsWeighted, "Mc collision weighted ", true);

  void processChargedJetsData(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                              aod::JetTracks const& tracksOfCollisions,
                              soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets)
  { 
    //Collision selection on : EventSelection, skipMBGapEvents, OccupancyCut
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    ///////////// leading track cut try : (because filter doesnt work)
    for (auto& jet : jets){
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      bool hasHighPtConstituent = false;
      registry.fill(HIST("h_jet_pt_initial_data"), jet.pt()); 
      for (auto& jetConstituent : jet.tracks_as<aod::JetTracks>()) {
        if (jetConstituent.pt() >= ptLeadingTrackCut) {
          hasHighPtConstituent = true;
          break; // Sortir de la boucle dès qu'un constituant valide est trouvé
        }
      }
      // Si un jet contient un constituant avec un pt > au critère, on l'analyse
      if (hasHighPtConstituent) {
        registry.fill(HIST("h_jet_pt_after_leadingtrackcut_data"), jet.pt()); 
        analyseCharged<false>(jet, tracksOfCollisions, jetSplittingsDataTable);//attention je donne TOUTES les traces pas juste les traces des jets comme dans jetSubstructure.cxx (jetConstituentS,jet.tracks_as<aod::JetTracks>()) en fait on s'en fou car c'est deja filtré dans analysisCharged() c'est meme mieux de tout donner si jamais il y manque
      }
    }
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsData, "charged jet substructure", false);

  void processChargedJetsEventWiseSubData(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                                          soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents> const& jets,
                                          aod::JetTracksSub const& tracksOfCollisions)
  {
    //rajouter les cuts de jetspectra

    //Collision selection on : EventSelection, skipMBGapEvents, OccupancyCut
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    // Leading track cut 
    for (auto& jet : jets){
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      bool hasHighPtConstituent = false;
      registry.fill(HIST("h_jet_pt_initial_data_eventwise"), jet.pt());
      // auto & jetConstituent0 = jet.tracks_as<aod::JetTracksSub>().iteratorAt(0)
      for (auto& jetConstituent : jet.tracks_as<aod::JetTracksSub>()) {
        if (jetConstituent.pt() >= ptLeadingTrackCut) {
          // LOGF(info, "Jet with leading constituent pt = %.2f found", jetConstituent.pt());
          // LOGF(info, "test1 ");
          hasHighPtConstituent = true;
          break; // Sortir de la boucle dès qu'un constituant valide est trouvé
        }
      }
      // Si un jet contient un constituant avec un pt > au critère, on l'analyse
      if (hasHighPtConstituent) {
        // LOGF(info, "test2 ");
        registry.fill(HIST("h_jet_pt_after_leadingtrackcut_data_eventwise"), jet.pt()); 
        analyseCharged<true>(jet, tracksOfCollisions, jetSplittingsDataSubTable);
        // registry.fill(HIST("h_jet_pt_after_grooming"), jet.pt()); //rajouter weight pour MC
      }
    }
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsEventWiseSubData, "eventwise-constituent subtracted charged jet substructure", false);

  void processChargedJetsMCD(soa::Filtered<aod::JetCollisionsMCD>::iterator const& collision,
                             soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                             aod::JetTracks const& tracks)
  { 
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }

    for (auto& jet : jets){
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      bool hasHighPtConstituent = false;
      registry.fill(HIST("h_jet_pt_initial_mcd"), jet.pt());  
      ///////////// leading track cut /////////////
      for (auto& jetConstituent : jet.tracks_as<aod::JetTracks>()) {
        if (jetConstituent.pt() >= ptLeadingTrackCut) {
          hasHighPtConstituent = true;
          break; // Sortir de la boucle dès qu'un constituant valide est trouvé
        }
      }
      // Si un jet contient un constituant avec un pt > au critère, on l'analyse
      if (hasHighPtConstituent) {
        registry.fill(HIST("h_jet_pt_after_leadingtrackcut_mcd"), jet.pt()); 
        analyseCharged<false>(jet, tracks, jetSplittingsMCDTable);
      }
    } 
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsMCD, "charged jet MCD substructure weighted", false);

  void processChargedJetsMCDWeighted(soa::Filtered<aod::JetCollisionsMCD>::iterator const& collision,
                                    //  aod::JetMcCollisions const&, //join the weight
                                     soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetEventWeights> const& jets,
                                     aod::JetTracks const& tracks)
  { 
    // LOGF(info, "processChargedJetsMCD: weight = %.4f", "1 : " ,jetweight);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    
    for (auto& jet : jets){
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      float jetweight = jet.eventWeight();
      LOGF(info, "jetweight = %.8f ",jetweight);
      float pTHat = 10. / (std::pow(jetweight, 1.0 / pTHatExponent));
      LOGF(info, "pTHat = %.8f ",pTHat);
      if (jet.pt() > pTHatMaxMCD * pTHat) {
        return;
      }
      bool hasHighPtConstituent = false;
      registry.fill(HIST("h_jet_pthat_initial_mcd"), pTHat); 
      registry.fill(HIST("h_jet_pthat_initial_mcd_weighted"), pTHat, jetweight); 
      registry.fill(HIST("h_jet_pt_initial_mcd"), jet.pt());  
      registry.fill(HIST("h_jet_pt_initial_mcd_weighted"), jet.pt(),jetweight);  

      ///////////// leading track cut /////////////
      for (auto& jetConstituent : jet.tracks_as<aod::JetTracks>()) {
        if (jetConstituent.pt() >= ptLeadingTrackCut) {
          hasHighPtConstituent = true;
          break; // Sortir de la boucle dès qu'un constituant valide est trouvé
        }
      }
      // Si un jet contient un constituant avec un pt > au critère, on l'analyse
      if (hasHighPtConstituent) {
        registry.fill(HIST("h_jet_pt_after_leadingtrackcut_mcd"), jet.pt()); 
        registry.fill(HIST("h_jet_pt_after_leadingtrackcut_mcd_weighted"), jet.pt(), jetweight); 
        LOGF(info, "MCD 7");
        LOGF(info, "jetweight = %.4f",jetweight);
        analyseCharged<false>(jet, tracks, jetSplittingsMCDTable, jetweight);
        // LOGF(info, "processChargedJetsMCD: weight = %.4f",jetweight);
        // LOGF(info, "processChargedJetsMCD: weight = %.4f", "1 : " ,jetweight, "2 : " , collision.mcCollision().weight());
      }
    }
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsMCDWeighted, "charged jet MCD substructure weighted", false);

    void processChargedJetsEventWiseSubMCD(soa::Filtered<aod::JetCollisionsMCD>::iterator const& collision,
                                           soa::Join<aod::ChargedMCDetectorLevelEventWiseSubtractedJets, aod::ChargedMCDetectorLevelEventWiseSubtractedJetConstituents> const& jets,
                                           aod::JetTracksSub const& tracks)
  { 
    //rajouter les cuts de jetspectra
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    ///////////// leading track cut /////////////
    for (auto& jet : jets){
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      registry.fill(HIST("h_jet_pt_initial_mcd_eventwise"), jet.pt()); 
      bool hasHighPtConstituent = false;
      for (auto& jetConstituent : jet.tracks_as<aod::JetTracksSub>()) {
        if (jetConstituent.pt() >= ptLeadingTrackCut) {
          hasHighPtConstituent = true;
          break; // Sortir de la boucle dès qu'un constituant valide est trouvé
        }
      }
      // Si un jet contient un constituant avec un pt > au critère, on l'analyse
      if (hasHighPtConstituent) {
        registry.fill(HIST("h_jet_pt_after_leadingtrackcut_mcd_eventwise"), jet.pt()); 
        analyseCharged<true>(jet, tracks, jetSplittingsMCDSubTable);
      }
    }
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsEventWiseSubMCD, "eventwise-constituent subtracted MCD charged jet substructure", false);

  void processChargedJetsEventWiseSubMCDWeighted(soa::Filtered<aod::JetCollisionsMCD>::iterator const& collision,
                                                soa::Join<aod::ChargedMCDetectorLevelEventWiseSubtractedJets, aod::ChargedMCDetectorLevelEventWiseSubtractedJetConstituents, aod::ChargedMCDetectorLevelEventWiseSubtractedJetEventWeights > const& jets,
                                                aod::JetTracksSub const& tracks)
{ 
  //rajouter les cuts de jetspectra
  if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
    return;
  }
  if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
    return;
  }
  ///////////// leading track cut /////////////
  for (auto& jet : jets){
    if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
      continue;
    }
    if (!isAcceptedJet<aod::JetTracks>(jet)) {
      continue;
    }
    float jetweight = jet.eventWeight();
    LOGF(info, "jetweight = %.8f ",jetweight);
    float pTHat = 10. / (std::pow(jetweight, 1.0 / pTHatExponent));
    LOGF(info, "pTHat = %.8f ",pTHat);
    if (jet.pt() > pTHatMaxMCD * pTHat) {
        return;
    }
    bool hasHighPtConstituent = false;
    registry.fill(HIST("h_jet_pthat_initial_mcd_eventwise"), pTHat); 
    registry.fill(HIST("h_jet_pthat_initial_mcd_eventwise_weighted"), pTHat,jetweight); 
    registry.fill(HIST("h_jet_pt_initial_mcd_eventwise"), jet.pt()); 
    registry.fill(HIST("h_jet_pt_initial_mcd_eventwise_weighted"), jet.pt(),jetweight); 
      
  
    ///////////// leading track cut /////////////
    for (auto& jetConstituent : jet.tracks_as<aod::JetTracksSub>()) {
      if (jetConstituent.pt() >= ptLeadingTrackCut) {
        hasHighPtConstituent = true;
      break; // Sortir de la boucle dès qu'un constituant valide est trouvé
      }
    }
    // Si un jet contient un constituant avec un pt > au critère, on l'analyse
    if (hasHighPtConstituent) {
      registry.fill(HIST("h_jet_pt_after_leadingtrackcut_mcd_eventwise"), jet.pt()); 
      registry.fill(HIST("h_jet_pt_after_leadingtrackcut_mcd_eventwise_weighted"), jet.pt(),jetweight); 
      analyseCharged<true>(jet, tracks, jetSplittingsMCDSubTable);
    }
  }
}
PROCESS_SWITCH(JetSubstructureTask, processChargedJetsEventWiseSubMCDWeighted, "Weighted eventwise-constituent subtracted MCD charged jet substructure ", false);

  void processChargedJetsMCP(soa::Filtered<aod::JetMcCollisions>::iterator const& mcCollision,
                             soa::SmallGroups<aod::JetCollisionsMCD> const& collisions,
                             soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                             aod::JetParticles const& particles)
  {
    // LOGF(info, " Entering processChargedJetsMCP " );

    // if (!jetderiveddatautilities::selectCollision(mcCollision, eventSelectionBits, skipMBGapEvents)) {
    //   return;
    // }
    // if (mcCollision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < mcCollision.trackOccupancyInTimeRange()) {
    //   return;
    // }


    //meme criteres que JetSpectra:
    bool mcLevelIsParticleLevel = true;
    if (std::abs(mcCollision.posZ()) > vertexZCut) {
      return;
    }
    if (collisions.size() < 1) {
      return;
    }
    bool hasSel8Coll = false;
    bool centralityIsGood = false;
    bool occupancyIsGood = false;
    for (auto const& collision : collisions) {
      if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
        hasSel8Coll = true;
      }
      if ((centralityMin < collision.centFT0M()) && (collision.centFT0M() < centralityMax)) {
        centralityIsGood = true;
      }
      if ((trackOccupancyInTimeRangeMin < collision.trackOccupancyInTimeRange()) && (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMax)) {
        occupancyIsGood = true;
      }
    }
    if (!hasSel8Coll) {
      return;
    }
    if (!centralityIsGood) {
      return;
    }
    if (!occupancyIsGood) {
      return;
    }
    for (auto& jet : jets){
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetParticles>(jet, mcLevelIsParticleLevel)) {
        continue;
      }
      // LOGF(info, " Entering boucle_jets " );
      bool hasHighPtConstituent = false;
      registry.fill(HIST("h_jet_pt_initial_mcp"), jet.pt()); 
      for (auto& jetConstituent : jet.tracks_as<aod::JetParticles>()) {
        // LOGF(info, " Entering boucle_jets_constituents " );
        if (jetConstituent.pt() >= ptLeadingTrackCut) {
          hasHighPtConstituent = true;
          break; // Sortir de la boucle dès qu'un constituant valide est trouvé
        }
      }
      if (hasHighPtConstituent) {
        // LOGF(info, " leading track cut applied " );
        registry.fill(HIST("h_jet_pt_after_leadingtrackcut_mcp"), jet.pt()); 
        //début de analyseCharged version MCP
        jetConstituents.clear();
        for (auto& jetConstituent : jet.tracks_as<aod::JetParticles>()) {
          // LOGF(info, " AnalyseCharged for MCP " );
          fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex(), static_cast<int>(JetConstituentStatus::track), pdg->Mass(jetConstituent.pdgCode()));
        }
        jetReclustering<true, false>(jet, jetSplittingsMCPTable , 1);
        //fin de analyseCharged version MCP
        // LOGF(info, "processChargedJetsMCP: weight = %.4f", mcCollision.weight());
      }
    }
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsMCP, "charged jet substructure on MC particle level", false);


  void processChargedJetsMCPWeighted(soa::Filtered<aod::JetMcCollisions>::iterator const& mcCollision,
                                      soa::SmallGroups<aod::JetCollisionsMCD> const& collisions,
                                      soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents , aod::ChargedMCParticleLevelJetEventWeights> const& jets,
                                      aod::JetParticles const& particles)
  {
  //meme criteres que JetSpectra:
  bool mcLevelIsParticleLevel = true;
  float eventWeight = mcCollision.weight();
  if (std::abs(mcCollision.posZ()) > vertexZCut) {
    return;
  }
  if (collisions.size() < 1) {
    return;
  }
  bool hasSel8Coll = false;
  bool centralityIsGood = false;
  bool occupancyIsGood = false;
  for (auto const& collision : collisions) {
    if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
      hasSel8Coll = true;
    }
    if ((centralityMin < collision.centFT0M()) && (collision.centFT0M() < centralityMax)) {
      centralityIsGood = true;
    }
    if ((trackOccupancyInTimeRangeMin < collision.trackOccupancyInTimeRange()) && (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMax)) {
      occupancyIsGood = true;
    }
  }
  if (!hasSel8Coll) {
    return;
  }
  if (!centralityIsGood) {
    return;
  }
  if (!occupancyIsGood) {
    return;
  }
  for (auto& jet : jets){
    if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
      continue;
    }
    if (!isAcceptedJet<aod::JetParticles>(jet, mcLevelIsParticleLevel)) {
      continue;
    }
    // LOGF(info, " Entering boucle_jets " );
    bool hasHighPtConstituent = false;
    float jetweight = jet.eventWeight();
    LOGF(info, "jetweight = %.8f",jetweight);
    double pTHat = 10. / (std::pow(jetweight, 1.0 / pTHatExponent));
    LOGF(info, "pTHat = %.8f ",pTHat);
    registry.fill(HIST("h_jet_pt_initial_mcp"), jet.pt());
    registry.fill(HIST("h_jet_pt_initial_mcp_weighted"), jet.pt(),jetweight);
    registry.fill(HIST("h_jet_phat_initial_mcp"), pTHat);
    registry.fill(HIST("h_jet_pthat_initial_mcp_weighted"), pTHat, jetweight); 
    for (auto& jetConstituent : jet.tracks_as<aod::JetParticles>()) {
      // LOGF(info, " Entering boucle_jets_constituents " );
      if (jetConstituent.pt() >= ptLeadingTrackCut) {
        hasHighPtConstituent = true;
        break; // Sortir de la boucle dès qu'un constituant valide est trouvé
      }
    }
    if (hasHighPtConstituent) {
      // LOGF(info, " leading track cut applied " );
      registry.fill(HIST("h_jet_pt_after_leadingtrackcut_mcp"), jet.pt());
      registry.fill(HIST("h_jet_pt_after_leadingtrackcut_mcp_weighted"), jet.pt(),jetweight);
      //début de analyseCharged version MCP
      jetConstituents.clear();
      for (auto& jetConstituent : jet.tracks_as<aod::JetParticles>()) {
        // LOGF(info, " AnalyseCharged for MCP " );
        fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex(), static_cast<int>(JetConstituentStatus::track), pdg->Mass(jetConstituent.pdgCode()));
      }
      jetReclustering<true, false>(jet, jetSplittingsMCPTable , jetweight);
      //fin de analyseCharged version MCP
      LOGF(info, "processChargedJetsMCP: weight = %.4f",jetweight);
    }
  }
}
PROCESS_SWITCH(JetSubstructureTask, processChargedJetsMCPWeighted, "charged jet substructure on MC particle level weighted", false);


};

void processJetsMCDMatchedMCP(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                        ChargedMCDMatchedJets const& mcdjets,
                        ChargedMCPMatchedJets const&,
                        aod::JetTracks const&, aod::JetParticles const&)
{
  if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
    return;
  }
  if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
    return;
  }

  for (const auto& mcdjet : mcdjets) {
    if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
      continue;
    }
    fillMatchedHistograms<ChargedMCDMatchedJets::iterator, ChargedMCPMatchedJets>(mcdjet);
  }
}
PROCESS_SWITCH(JetSubstructureTask, processJetsMCDMatchedMCP, "matched mcp and mcd jets", false);

void processJetsMCDMatchedMCPWeighted(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                                      ChargedMCDMatchedJetsWeighted const& mcdjets,
                                      ChargedMCPMatchedJetsWeighted const&,
                                      aod::JetTracks const&, aod::JetParticles const&)
{
  if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
    return;
  }
  if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
    return;
  }
  for (const auto& mcdjet : mcdjets) {
    if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
      continue;
    }
  fillMatchedHistograms<ChargedMCDMatchedJetsWeighted::iterator, ChargedMCPMatchedJetsWeighted>(mcdjet, mcdjet.eventWeight());
  }
}
PROCESS_SWITCH(JetSubstructureTask, processJetsMCDMatchedMCPWeighted, "matched mcp and mcd jets with weighted events", false);

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  return WorkflowSpec{adaptAnalysisTask<JetSubstructureTask>(
    cfgc, TaskName{"jet-substructure-softdrop"})};
}