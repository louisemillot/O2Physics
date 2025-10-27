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
#include <fstream>
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


#include <optional>

#include "EventFiltering/filterTables.h"


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct JetSubstructureTask {

  using ChargedMCDMatchedJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>; //1
  using ChargedMCPMatchedJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>; //2
  using ChargedMCDMatchedtoMCDEventWise = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCDetectorLevelEventWiseSubtractedJets>; //3 matching MCD to MCD eventwise
  using ChargedMCDEventWiseMatchedtoMCD = soa::Join<aod::ChargedMCDetectorLevelEventWiseSubtractedJets, aod::ChargedMCDetectorLevelEventWiseSubtractedJetConstituents, aod::ChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToChargedMCDetectorLevelJets>; //4 matching MCD eventwise to MCD


  using ChargedMCDMatchedJetsWeighted = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetEventWeights>;
  using ChargedMCPMatchedJetsWeighted = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetEventWeights>;
  using ChargedMCDEventWiseMatchedtoMCDWeighted = soa::Join<aod::ChargedMCDetectorLevelEventWiseSubtractedJets, aod::ChargedMCDetectorLevelEventWiseSubtractedJetConstituents, aod::ChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelEventWiseSubtractedJetEventWeights>; //4 matching MCD eventwise to MCD weighted


  using ChargedMCDEventWise = soa::Join<aod::ChargedMCDetectorLevelEventWiseSubtractedJets, aod::ChargedMCDetectorLevelEventWiseSubtractedJetConstituents>;

  Produces<aod::ChargedSPs> jetSplittingsDataTable;
  Produces<aod::ChargedEventWiseSubtractedSPs> jetSplittingsDataSubTable;
  Produces<aod::ChargedMCDetectorLevelSPs> jetSplittingsMCDTable;
  Produces<aod::ChargedMCParticleLevelSPs> jetSplittingsMCPTable;
  Produces<aod::ChargedMCDetectorLevelEventWiseSubtractedSPs> jetSplittingsMCDSubTable;
  // Produces<aod::ChargedEventWiseSubtractedSPs> jetSplittingsMCPSubTable;

  Configurable<float> zCut{"zCut", 0.1, "soft drop z cut"};
  Configurable<float> beta{"beta", 0.0, "soft drop beta"};
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
  Configurable<int> nBinsEta{"nBinsEta", 200, "number of bins for eta axes"};
  Configurable<float> selectedJetsRadius{"selectedJetsRadius", 0.2, "resolution parameter for histograms without radius"};







  Service<o2::framework::O2DatabasePDG> pdg;
  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  std::vector<float> energyMotherVec;
  std::vector<float> ptLeadingVec;
  std::vector<float> ptSubLeadingVec;
  std::vector<std::pair<float,float>> thetagMCDVec;
  std::vector<std::pair<float,float>> thetagMCPVec;
  std::vector<std::pair<float,float>> thetagMCDEventWiseVec;
  std::vector<float> thetaVec;
  std::vector<float> nSub;
  
  HistogramRegistry registry;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  void init(InitContext const&)
  {
    AxisSpec jetPtAxis = {200, 0., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec jetPtAxisMCD = {200, 0., 200., "#it{p}_{T}^{mcd} (GeV/#it{c})"};
    AxisSpec jetPtAxisMCP = {200, 0., 200., "#it{p}_{T}^{mcp} (GeV/#it{c})"};
    AxisSpec centralityAxis = {1200, -10., 110., "Centrality"};
    AxisSpec jetEtaAxis = {nBinsEta, -1.0, 1.0, "#eta"};
    AxisSpec phiAxis = {160, -1.0, 7.0, "#varphi"};
    AxisSpec thetagAxisMCD = {200, 0.0, 1.1, "{#theta}_{g}^{mcd}"}; 
    AxisSpec thetagAxisMCP = {200, 0.0, 1.1, "{#theta}_{g}^{mcp}"}; 
    AxisSpec thetagAxisMCDEventWise = {200, 0.0, 1.1, "{#theta}_{g}^{mcd_eventwise}"}; 




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
        registry.get<TH1>(HIST("h_mcColl_counts_weight"))->GetXaxis()->SetBinLabel(7, "centrality20to60");

        registry.add("h_mc_zvertex_weight", "position of collision ;#it{Z} (cm)", {HistType::kTH1F, {{300, -15.0, 15.0}}});
      }
    }

    //data
    if (doprocessChargedJetsData) {
      registry.add("h_jet_pt_initial_data", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      registry.add("h_jet_pt_after_leadingtrackcut_data", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      registry.add("h2_centrality_jet_pt", "centrality vs. jet pT;centrality; #it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH2F, {centralityAxis, jetPtAxis}});

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

    if (doprocessJetsMCDMatchedMCP || doprocessJetsMCDMatchedMCPWeighted || doprocessJetsMCDEventWiseMatchedMCP) {
      if (checkGeoMatched) {
        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint", "pT mcd vs. pT mcp;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisMCD, jetPtAxisMCP}});
        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcpetaconstraint", "pT mcd vs. pT mcp;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisMCD, jetPtAxisMCP}});
        registry.add("h2_jet_eta_mcd_jet_eta_mcp_matchedgeo", "Eta mcd vs. Eta mcp;#eta_{jet}^{mcd};#eta_{jet}^{mcp}", {HistType::kTH2F, {jetEtaAxis, jetEtaAxis}});
        registry.add("h2_jet_phi_mcd_jet_phi_mcp_matchedgeo_mcdetaconstraint", "Phi mcd vs. Phi mcp;#varphi_{jet}^{mcd};#varphi_{jet}^{mcp}", {HistType::kTH2F, {phiAxis, phiAxis}});
        registry.add("h2_jet_phi_mcd_jet_phi_mcp_matchedgeo_mcpetaconstraint", "Phi mcd vs. Phi mcp;#varphi_{jet}^{mcd};#varphi_{jet}^{mcp}", {HistType::kTH2F, {phiAxis, phiAxis}});
        registry.add("h2_jet_ntracks_mcd_jet_ntracks_mcp_matchedgeo", "Ntracks mcd vs. Ntracks mcp;N_{jet tracks}^{mcd};N_{jet tracks}^{mcp}", {HistType::kTH2F, {{200, -0.5, 199.5}, {200, -0.5, 199.5}}});
        registry.add("h2_jet_pt_mcp_jet_pt_diff_matchedgeo", "jet mcp pT vs. delta pT / jet mcp pt;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); (#it{p}_{T,jet}^{mcp} (GeV/#it{c}) - #it{p}_{T,jet}^{mcd} (GeV/#it{c})) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_mcd_jet_pt_diff_matchedgeo", "jet mcd pT vs. delta pT / jet mcd pt;#it{p}_{T,jet}^{mcd} (GeV/#it{c}); (#it{p}_{T,jet}^{mcd} (GeV/#it{c}) - #it{p}_{T,jet}^{mcp} (GeV/#it{c})) / #it{p}_{T,jet}^{mcd} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_mcp_jet_pt_ratio_matchedgeo", "jet mcp pT vs. jet mcd pT / jet mcp pt;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); #it{p}_{T,jet}^{mcd} (GeV/#it{c}) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 5.0}}});
        registry.add("h2_jet_thetag_mcd_jet_thetag_mcp_matchedgeo","#theta_{g}^{mcd} vs. #theta_{g}^{mcp};#theta_{g}^{mcd};#theta_{g}^{mcp}", {HistType::kTH2F, {thetagAxisMCD, thetagAxisMCP}});
        registry.add("h2_thetagMCD_vs_thetagMCP_pt_norange","#theta_{g}^{mcd} vs. #theta_{g}^{mcp};#theta_{g}^{mcd};#theta_{g}^{mcp}",{HistType::kTH2F, {thetagAxisMCD, thetagAxisMCP}});
        registry.add("h2_thetagMCD_vs_thetagMCP_pt_60_80","#theta_{g}^{mcd} vs. #theta_{g}^{mcp};#theta_{g}^{mcd};#theta_{g}^{mcp}",{HistType::kTH2F, {thetagAxisMCD, thetagAxisMCP}});
        registry.add("h4_ptMCD_ptMCP_thetagMCD_thetagMCP_norange","p_{T}^{MCD} vs p_{T}^{MCP} vs #theta_{g}^{MCD} vs #theta_{g}^{MCP};p_{T}^{MCD} (GeV/#it{c});p_{T}^{MCP} (GeV/#it{c});#theta_{g}^{MCD};#theta_{g}^{MCP}",{HistType::kTHnSparseF, {jetPtAxisMCD, jetPtAxisMCP, thetagAxisMCD, thetagAxisMCP}});


        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint_eventwise", "pT mcdeventwise vs. pT mcp;#it{p}_{T,jet}^{mcdeventwise} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisMCD, jetPtAxisMCP}});
        registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcpetaconstraint_eventwise", "pT mcdeventwise vs. pT mcp;#it{p}_{T,jet}^{mcdeventwise} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisMCD, jetPtAxisMCP}});
        registry.add("h2_jet_eta_mcd_jet_eta_mcp_matchedgeo_eventwise", "Eta mcdeventwise vs. Eta mcp;#eta_{jet}^{mcdeventwise};#eta_{jet}^{mcp}", {HistType::kTH2F, {jetEtaAxis, jetEtaAxis}});
        registry.add("h2_jet_phi_mcd_jet_phi_mcp_matchedgeo_mcdetaconstraint_eventwise", "Phi mcdeventwise vs. Phi mcp;#varphi_{jet}^{mcdeventwise};#varphi_{jet}^{mcp}", {HistType::kTH2F, {phiAxis, phiAxis}});
        registry.add("h2_jet_phi_mcd_jet_phi_mcp_matchedgeo_mcpetaconstraint_eventwise", "Phi mcdeventwise vs. Phi mcp;#varphi_{jet}^{mcdeventwise};#varphi_{jet}^{mcp}", {HistType::kTH2F, {phiAxis, phiAxis}});
        registry.add("h2_jet_ntracks_mcd_jet_ntracks_mcp_matchedgeo_eventwise", "Ntracks mcdeventwise vs. Ntracks mcp;N_{jet tracks}^{mcdeventwise};N_{jet tracks}^{mcp}", {HistType::kTH2F, {{200, -0.5, 199.5}, {200, -0.5, 199.5}}});
        registry.add("h2_jet_pt_mcp_jet_pt_diff_matchedgeo_eventwise", "jet mcp pT vs. delta pT / jet mcp pt;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); (#it{p}_{T,jet}^{mcp} (GeV/#it{c}) - #it{p}_{T,jet}^{mcdeventwise} (GeV/#it{c})) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_mcd_jet_pt_diff_matchedgeo_eventwise", "jet mcdeventwise pT vs. delta pT / jet mcdeventwise pt;#it{p}_{T,jet}^{mcdeventwise} (GeV/#it{c}); (#it{p}_{T,jet}^{mcdeventwise} (GeV/#it{c}) - #it{p}_{T,jet}^{mcp} (GeV/#it{c})) / #it{p}_{T,jet}^{mcdeventwise} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_mcp_jet_pt_ratio_matchedgeo_eventwise", "jet mcp pT vs. jet mcdeventwise pT / jet mcp pt;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); #it{p}_{T,jet}^{mcdeventwise} (GeV/#it{c}) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 5.0}}});
        registry.add("h2_jet_thetag_mcd_jet_thetag_mcp_matchedgeo_eventwise","#theta_{g}^{mcdeventwise} vs. #theta_{g}^{mcp};#theta_{g}^{mcdeventwise};#theta_{g}^{mcp}", {HistType::kTH2F, {thetagAxisMCDEventWise, thetagAxisMCP}});
        registry.add("h2_thetagMCD_vs_thetagMCP_pt_norange_eventwise","#theta_{g}^{mcdeventwise} vs. #theta_{g}^{mcp};#theta_{g}^{mcdeventwise};#theta_{g}^{mcp}",{HistType::kTH2F, {thetagAxisMCDEventWise, thetagAxisMCP}});
        registry.add("h2_thetagMCD_vs_thetagMCP_pt_60_80_eventwise","#theta_{g}^{mcdeventwise} vs. #theta_{g}^{mcp};#theta_{g}^{mcdeventwise};#theta_{g}^{mcp}",{HistType::kTH2F, {thetagAxisMCDEventWise, thetagAxisMCP}});


      }
      // if (checkPtMatched) {
      //   registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedpt_mcdetaconstraint", "pT mcd vs. pT mcp;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisMCD, jetPtAxisMCP}});
      //   registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedpt_mcpetaconstraint", "pT mcd vs. pT mcp;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisMCD, jetPtAxisMCP}});
      //   registry.add("h2_jet_eta_mcd_jet_eta_mcp_matchedpt", "Eta mcd vs. Eta mcp;#eta_{jet}^{mcd};#eta_{jet}^{mcp}", {HistType::kTH2F, {jetEtaAxis, jetEtaAxis}});
      //   registry.add("h2_jet_phi_mcd_jet_phi_mcp_matchedgpt_mcdetaconstraint", "Phi mcd vs. Phi mcp;#varphi_{jet}^{mcd};#varphi_{jet}^{mcp}", {HistType::kTH2F, {phiAxis, phiAxis}});
      //   registry.add("h2_jet_phi_mcd_jet_phi_mcp_matchedgpt_mcpetaconstraint", "Phi mcd vs. Phi mcp;#varphi_{jet}^{mcd};#varphi_{jet}^{mcp}", {HistType::kTH2F, {phiAxis, phiAxis}});
      //   registry.add("h2_jet_ntracks_mcd_jet_ntracks_mcp_matchedpt", "Ntracks mcd vs. Ntracks mcp;N_{jet tracks}^{mcd};N_{jet tracks}^{mcp}", {HistType::kTH2F, {{200, -0.5, 199.5}, {200, -0.5, 199.5}}});
      //   registry.add("h2_jet_pt_mcp_jet_pt_diff_matchedpt", "jet mcp pT vs. delta pT / jet mcp pt;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); (#it{p}_{T,jet}^{mcp} (GeV/#it{c}) - #it{p}_{T,jet}^{mcd} (GeV/#it{c})) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
      //   registry.add("h2_jet_pt_mcd_jet_pt_diff_matchedpt", "jet mcd pT vs. delta pT / jet mcd pt;#it{p}_{T,jet}^{mcd} (GeV/#it{c}); (#it{p}_{T,jet}^{mcd} (GeV/#it{c}) - #it{p}_{T,jet}^{mcp} (GeV/#it{c})) / #it{p}_{T,jet}^{mcd} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
      //   registry.add("h2_jet_pt_mcp_jet_pt_ratio_matchedpt", "jet mcp pT vs. jet mcd pT / jet mcp pt;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); #it{p}_{T,jet}^{mcd} (GeV/#it{c}) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 5.0}}});
      // }
      // if (checkGeoPtMatched) {
      //   registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeopt_mcdetaconstraint", "pT mcd vs. pT mcp;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisMCD, jetPtAxisMCP}});
      //   registry.add("h2_jet_pt_mcd_jet_pt_mcp_matchedgeopt_mcpetaconstraint", "pT mcd vs. pT mcp;#it{p}_{T,jet}^{mcd} (GeV/#it{c});#it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisMCD, jetPtAxisMCP}});
      //   registry.add("h2_jet_eta_mcd_jet_eta_mcp_matchedgeopt", "Eta mcd vs. Eta mcp;#eta_{jet}^{mcd};#eta_{jet}^{mcp}", {HistType::kTH2F, {jetEtaAxis, jetEtaAxis}});
      //   registry.add("h2_jet_phi_mcd_jet_phi_mcp_matchedgeopt_mcdetaconstraint", "Phi mcd vs. Phi mcp;#varphi_{jet}^{mcd};#varphi_{jet}^{mcp}", {HistType::kTH2F, {phiAxis, phiAxis}});
      //   registry.add("h2_jet_phi_mcd_jet_phi_mcp_matchedgeopt_mcpetaconstraint", "Phi mcd vs. Phi mcp;#varphi_{jet}^{mcd};#varphi_{jet}^{mcp}", {HistType::kTH2F, {phiAxis, phiAxis}});
      //   registry.add("h2_jet_ntracks_mcd_jet_ntracks_mcp_matchedgeopt", "Ntracks mcd vs. Ntracks mcp;N_{jet tracks}^{mcd};N_{jet tracks}^{mcp}", {HistType::kTH2F, {{200, -0.5, 199.5}, {200, -0.5, 199.5}}});
      //   registry.add("h2_jet_pt_mcp_jet_pt_diff_matchedgeopt", "jet mcp pT vs. delta pT / jet mcp pt;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); (#it{p}_{T,jet}^{mcp} (GeV/#it{c}) - #it{p}_{T,jet}^{mcd} (GeV/#it{c})) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
      //   registry.add("h2_jet_pt_mcd_jet_pt_diff_matchedgeopt", "jet mcd pT vs. delta pT / jet mcd pt;#it{p}_{T,jet}^{mcd} (GeV/#it{c}); (#it{p}_{T,jet}^{mcd} (GeV/#it{c}) - #it{p}_{T,jet}^{mcp} (GeV/#it{c})) / #it{p}_{T,jet}^{mcd} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
      //   registry.add("h2_jet_pt_mcp_jet_pt_ratio_matchedgeopt", "jet mcp pT vs. jet mcd pT / jet mcp pt;#it{p}_{T,jet}^{mcp} (GeV/#it{c}); #it{p}_{T,jet}^{mcd} (GeV/#it{c}) / #it{p}_{T,jet}^{mcp} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 5.0}}});
      // }
    }

    registry.add("h2_substructure_pt_mcd_substructure_pt_mcp_matchedgeo", 
      "pT substructure MCD vs pT MCP;#it{p}_{T,substructure}^{mcd} (GeV/#it{c});#it{p}_{T,substructure}^{mcp} (GeV/#it{c})", 
      {HistType::kTH2F, {jetPtAxisMCD, jetPtAxisMCP}});

    registry.add("h2_substructure_eta_mcd_substructure_eta_mcp_matchedgeo", 
      "Eta substructure MCD vs Eta MCP;#eta_{substructure}^{mcd};#eta_{substructure}^{mcp}", 
      {HistType::kTH2F, {jetEtaAxis, jetEtaAxis}});

    registry.add("h2_substructure_phi_mcd_substructure_phi_mcp_matchedgeo", 
      "Phi substructure MCD vs Phi MCP;#varphi_{substructure}^{mcd};#varphi_{substructure}^{mcp}", 
      {HistType::kTH2F, {phiAxis, phiAxis}});


    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    
  }

  // Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter collisionFilter = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centFT0C >= centralityMin && aod::jcollision::centFT0C < centralityMax);
  Filter mcCollisionFilter = (nabs(aod::jmccollision::posZ) < vertexZCut && aod::jmccollision::centFT0C >= centralityMin && aod::jmccollision::centFT0C < centralityMax);

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

  int count_surMCP = 0;
  int countthetagMCD_MCD_surMCP = 0;
  int countthetagMCP_MCD_surMCP = 0;
  int countthetagMCD_MCP_surMCP = 0;
  int countthetagMCP_MCP_surMCP = 0;
  std::vector<float> thetagMCDVecMatched;
  std::vector<float> thetagMCPVecMatched;
  // template <typename TBase, typename TTag, typename TableMCD, typename TableMCP>
  template <typename TBase, typename TTag >
  void fillMatchedHistograms(TBase const& jetMCD,
                            //  TableMCD& splittingTableMCD,
                            //  TableMCP& splittingTableMCP,
                            const std::vector<std::pair<float, float>>& thetagMCDVec,
                            const std::vector<std::pair<float, float>>& thetagMCPVec,
                            float weight = 1.0)
  { 
    // LOGF(info, "==== Jet numéro %d ====", jetMCD.globalIndex());
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
            count_surMCP++;
            double dpt = jetMCP.pt() - jetMCD.pt();
            if (jetfindingutilities::isInEtaAcceptance(jetMCD, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
              for (const auto& [thetagMCD, ptMCD] : thetagMCDVec) {
                if (ptMCD == jetMCD.pt()) {
                countthetagMCD_MCD_surMCP++;
                  for (const auto& [thetagMCP, ptMCP] : thetagMCPVec) {
                      if (ptMCP == jetMCP.pt()) {
                        countthetagMCP_MCD_surMCP++;
                        registry.fill(HIST("h2_thetagMCD_vs_thetagMCP_pt_norange"), thetagMCD, thetagMCP, weight);
                        // registry.fill(HIST("h4_ptMCD_ptMCP_thetagMCD_thetagMCP_norange"),jetMCD.pt(), jetMCP.pt(), thetagMCD, thetagMCP, weight);
                        // LOGF(info, "thetagMCD = %.4f, ptMCD = %.4f, thetagMCP = %.4f, ptMCP = %.4f", thetagMCD, ptMCD, thetagMCP, ptMCP);
                        // if (ptMCP >= 20.0 && ptMCP <= 80.0) {
                        //  registry.fill(HIST("h2_thetagMCD_vs_thetagMCP_pt_60_80"), thetagMCD, thetagMCP, weight);
                      } 
                  }
                }
              }
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
    std::cout << "nombre de MCD-MCP matchés - sur MCP : " << count_surMCP << std::endl;
    std::cout << "nombre de thetagMCD trouvés boucle for sur MCD - sur MCP : " << countthetagMCD_MCD_surMCP << std::endl;
    std::cout << "nombre de thetagMCP trouvés boucle for sur MCD - sur MCP : " << countthetagMCP_MCD_surMCP << std::endl; //nombre de thetagMCP trouvés parmis les thetagMCD trouves
    std::cout << "nombre de thetagMCP trouvés boucle for sur MCP - sur MCP : " << countthetagMCP_MCP_surMCP << std::endl;
    std::cout << "nombre de thetagMCD trouvés boucle for sur MCP - sur MCP : " << countthetagMCD_MCP_surMCP << std::endl; //nombre de thetagMCD trouvés parmis les thetagMCP trouves

    // std::cout << "Nombre de valeurs dans thetagMCDVec (colonne 1) = " << thetagMCDVec.size() << std::endl;
    // std::cout << "Nombre de valeurs dans thetagMCPVec (colonne 1) = " << thetagMCPVec.size() << std::endl;
    // fill pt matched histograms (a faire)
    // fill geometry and pt histograms (a faire)
  }

  int countMCDEW_MCD = 0;
  int countMCD_MCP = 0;
  template <typename TMCDEventWise, typename TMCDtoMCP, typename TMCP > // TMCDEventWise : ChargedMCDEventWiseMatchedtoMCD ; TMCDtoMCP : ChargedMCDMatchedJets ; TMCP : ChargedMCPMatchedJets
  void fillMatchedHistogramsEventWise(TMCDEventWise const& jetMCDEventWise, //le jetMCDEventWise est ici le TTag : ChargedMCDEventWiseMatchedtoMCD //iterator
                                      TMCP const&, //le jetMCP ici est le TBase : ChargedMCDMatchedtoMCDEventWise
                                      const std::vector<std::pair<float, float>>& thetagMCDEventWiseVec,
                                      const std::vector<std::pair<float, float>>& thetagMCPVec,
                                      float weight = 1.0) 
  {
    // LOGF(info, " fillMatchedHistogramsEventWise " );
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jetMCDEventWise.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin) {
      return;
    }
    // fill geometry matched histograms
    if (checkGeoMatched) { //true or false
      if (jetMCDEventWise.has_matchedJetGeo()) { //si il y a un match geometric entre MCD et MCDEventWise - 
        // for (const auto& jetMCP : (jetMCDEventWise.template matchedJetGeo_as<std::decay_t<TMCDtoMCP>>()).template matchedJetGeo_as<std::decay_t<TMCP>>()) { // - alors on boucle sur MCP qui ont un matching: MCDEventWise - MCD - MCP !!!! 2 ETAPES car sinon on ne peut pas build !!!!
        for (const auto& jetMCD : jetMCDEventWise.template matchedJetGeo_as<std::decay_t<TMCDtoMCP>>()){ // - alors on boucle sur les MCD qui ont un matching avec EventWiseMCD: MCDEventWise - MCD
          countMCDEW_MCD++;
          if (jetMCD.has_matchedJetGeo()) {
            for (const auto& jetMCP : jetMCD.template matchedJetGeo_as<std::decay_t<TMCP>>()){ // - puis on boucle sur MCD qui ont un matching avec MCP: MCD - MCP  !!!! 2 ETAPES car sinon on ne peut pas build !!!!
              if (jetMCP.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin) {
                continue;
              }
              if (jetMCDEventWise.r() == round(selectedJetsRadius * 100.0f)) { 
                countMCD_MCP++;
                double dpt = jetMCP.pt() - jetMCDEventWise.pt(); 
                if (jetfindingutilities::isInEtaAcceptance(jetMCDEventWise, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
                  for (const auto& [thetagMCD, ptMCD] : thetagMCDEventWiseVec) {
                    if (ptMCD == jetMCDEventWise.pt()) {
                      // LOGF(info, "thetagMCD = %.4f, ptMCD = %.4f, jetMCD.pt() = %.4f", thetagMCD, ptMCD, jetMCD.pt());
                      for (const auto& [thetagMCP, ptMCP] : thetagMCPVec) {
                        if (ptMCP == jetMCP.pt()) {
                          registry.fill(HIST("h2_thetagMCD_vs_thetagMCP_pt_norange_eventwise"), thetagMCD, thetagMCP, weight);
                          // LOGF(info, "thetagMCD = %.4f, ptMCD = %.4f, thetagMCP = %.4f, ptMCP = %.4f", thetagMCD, ptMCD, thetagMCP, ptMCP);
                          // if (ptMCP >= 20.0 && ptMCP <= 80.0) {
                            // registry.fill(HIST("h2_thetagMCD_vs_thetagMCP_pt_60_80_eventwise"), thetagMCD, thetagMCP, weight);
                          // } 
                        }
                      }
                    }
                  }
                  registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint_eventwise"), jetMCDEventWise.pt(), jetMCP.pt(), weight);
                  registry.fill(HIST("h2_jet_phi_mcd_jet_phi_mcp_matchedgeo_mcdetaconstraint_eventwise"), jetMCDEventWise.phi(), jetMCP.phi(), weight);
                  registry.fill(HIST("h2_jet_pt_mcd_jet_pt_diff_matchedgeo_eventwise"), jetMCDEventWise.pt(), dpt / jetMCDEventWise.pt(), weight);
                  registry.fill(HIST("h2_jet_ntracks_mcd_jet_ntracks_mcp_matchedgeo_eventwise"), jetMCDEventWise.tracksIds().size(), jetMCP.tracksIds().size(), weight);
                }
                if (jetfindingutilities::isInEtaAcceptance(jetMCP, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
                  registry.fill(HIST("h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcpetaconstraint_eventwise"), jetMCDEventWise.pt(), jetMCP.pt(), weight);
                  registry.fill(HIST("h2_jet_phi_mcd_jet_phi_mcp_matchedgeo_mcpetaconstraint_eventwise"), jetMCDEventWise.phi(), jetMCP.phi(), weight);
                  registry.fill(HIST("h2_jet_pt_mcp_jet_pt_diff_matchedgeo_eventwise"), jetMCP.pt(), dpt / jetMCP.pt(), weight);
                  registry.fill(HIST("h2_jet_pt_mcp_jet_pt_ratio_matchedgeo_eventwise"), jetMCP.pt(), jetMCDEventWise.pt() / jetMCP.pt(), weight);
                }
                registry.fill(HIST("h2_jet_eta_mcd_jet_eta_mcp_matchedgeo_eventwise"), jetMCDEventWise.eta(), jetMCP.eta(), weight);
              }
            }
          }
        }
      }
    }
    std::cout << "nombre de MCDEW-MCD matchés : " << countMCDEW_MCD << std::endl;
    std::cout << "nombre de MCD-MCP matchés aprês MCDEW-MCD: " << countMCD_MCP << std::endl;
    // fill pt matched histograms (a faire)
    // fill geometry and pt histograms (a faire)
  }
    

  template <bool isMCP, bool isSubtracted, typename T, typename U>
  std::tuple<std::vector<std::pair<float, float>>, std::vector<std::pair<float, float>>, std::vector<std::pair<float, float>> >
  jetReclustering(T const& jet, U& splittingTable, double weight)
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
    // LOGF(info, "==== Jet numéro %d ====", jet.globalIndex());

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
            thetagMCDVec.push_back({thetag, jet.pt()});
            // LOGF(info, "thetagMCD: %.4f et ptMCD: %.4f", thetag, jet.pt() );
            // LOGF(info, "thetagMCD = %.4f, ptJet = %.4f , Jet numéro %d ", thetag, jet.pt(), jet.globalIndex());


          }
          if constexpr (!isSubtracted && isMCP) { //MCP only no sub
            registry.fill(HIST("h2_jet_pt_part_jet_zg_part"), jet.pt(), zg, weight);
            registry.fill(HIST("h2_jet_pt_part_jet_rg_part"), jet.pt(), rg, weight);
            registry.fill(HIST("h2_jet_pt_part_jet_thetag_part"), jet.pt(), thetag, weight);
            registry.fill(HIST("h_jet_thetag_MCP"), thetag, weight);
            registry.fill(HIST("h_jet_zg_MCP"), zg, weight);
            thetagMCPVec.push_back({thetag, jet.pt()});
            // LOGF(info, "thetagMCP: %.4f et ptMCP: %.4f", thetag, jet.pt());

          }
          if constexpr (isSubtracted && !isMCP) { //data & MCD sub
            registry.fill(HIST("h2_jet_pt_jet_zg_eventwiseconstituentsubtracted"),jet.pt(), zg, weight);
            registry.fill(HIST("h2_jet_pt_jet_rg_eventwiseconstituentsubtracted"),jet.pt(), rg, weight);
            registry.fill(HIST("h2_jet_pt_jet_thetag_eventwiseconstituentsubtracted"),jet.pt(), thetag, weight);
            registry.fill(HIST("h_jet_thetag_eventwiseconstituentsubtracted"), thetag, weight);
            registry.fill(HIST("h_jet_zg_eventwiseconstituentsubtracted"), zg, weight);
            thetagMCDEventWiseVec.push_back({thetag, jet.pt()});

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
    // for (const auto& theta : thetaVec) { // boucle a changer pour acceder a l'element quand tu veux ptleading aussi 
    //   LOGF(info, "Theta: %.4f", theta);
    // }
    // for (const auto& [thetag, ptJet] : thetagMCDVec) {
    //   LOGF(info, "thetagMCD = %.4f, ptJet = %.4f , Jet numéro %d ", thetag, ptJet, jet.globalIndex());
    // }
    std::ofstream logFile1("thetagMCDVec.log"); // ouvre/crée le fichier
    logFile1 << "thetag\tpt\n"; // optionnel : écrire les en-têtes de colonnes
    for (const auto& [thetag, pt] : thetagMCDVec) {
      logFile1 << thetag << "\t" << pt << "\n"; // valeurs séparées par tabulation
    }
    logFile1.close();

    std::ofstream logFile2("thetagMCPVec.log"); // ouvre/crée le fichier
    logFile2 << "thetag\tpt\n"; // optionnel : écrire les en-têtes de colonnes
    for (const auto& [thetag, pt] : thetagMCPVec) {
      logFile2 << thetag << "\t" << pt << "\n"; // valeurs séparées par tabulation
    }
    logFile2.close();

    std::ofstream logFile3("thetagMCDEventWiseVec.log"); // ouvre/crée le fichier
    logFile3 << "thetag\tpt\n"; // optionnel : écrire les en-têtes de colonnes
    for (const auto& [thetag, pt] : thetagMCDEventWiseVec) {
      logFile3 << thetag << "\t" << pt << "\n"; // valeurs séparées par tabulation
    }
    logFile3.close();

    return {thetagMCDVec, thetagMCPVec, thetagMCDEventWiseVec};

  

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
    registry.fill(HIST("h2_centrality_occupancy"), collision.centFT0C(), collision.trackOccupancyInTimeRange());
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
  registry.fill(HIST("h2_centrality_occupancy"), collision.centFT0C(), collision.trackOccupancyInTimeRange());
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
    if ((centralityMin < collision.centFT0C()) && (collision.centFT0C() < centralityMax)) {
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
    bool centrality20to60 = false;

    for (auto const& collision : collisions) {
      if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
        hasSel8Coll = true;
      }
      if ((centralityMin < collision.centFT0C()) && (collision.centFT0C() < centralityMax)) {
        centralityIsGood = true;
      }
      if ((trackOccupancyInTimeRangeMin < collision.trackOccupancyInTimeRange()) && (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMax)) {
        occupancyIsGood = true;
      }
      if (collision.centFT0C() > 20 && collision.centFT0C() < 60) {
        centrality20to60 = true;
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

    if (!centrality20to60) {
      return;
    }
    registry.fill(HIST("h_mcColl_counts_weight"), 6.5, eventWeight);


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
    // LOGF(info, "test0 ");
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    // LOGF(info, "test1 ");
    LOGF(info, "Number of jets in event = %d", jets.size());
    ///////////// leading track cut try : (because filter doesnt work)
    for (auto& jet : jets){
      // LOGF(info, "test2 ");
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      // LOGF(info, "test3 ");
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      // LOGF(info, "test4 ");
      bool hasHighPtConstituent = false;
      registry.fill(HIST("h_jet_pt_initial_data"), jet.pt()); 
      for (auto& jetConstituent : jet.tracks_as<aod::JetTracks>()) {
        if (jetConstituent.pt() >= ptLeadingTrackCut) {
          hasHighPtConstituent = true;
          // LOGF(info, "test5 ");
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
    LOGF(info, "test0 ");
    // Leading track cut 
    for (auto& jet : jets){
      LOGF(info, "test1.2 ");
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      LOGF(info, "test1.5 ");
      bool hasHighPtConstituent = false;
      registry.fill(HIST("h_jet_pt_initial_data_eventwise"), jet.pt());
      // auto & jetConstituent0 = jet.tracks_as<aod::JetTracksSub>().iteratorAt(0)
      for (auto& jetConstituent : jet.tracks_as<aod::JetTracksSub>()) {
        if (jetConstituent.pt() >= ptLeadingTrackCut) {
          // LOGF(info, "Jet with leading constituent pt = %.2f found", jetConstituent.pt());
          LOGF(info, "test1 ");
          hasHighPtConstituent = true;
          break; // Sortir de la boucle dès qu'un constituant valide est trouvé
        }
      }
      // Si un jet contient un constituant avec un pt > au critère, on l'analyse
      if (hasHighPtConstituent) {
        LOGF(info, "test2 ");
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
    // LOGF(info, "processChargedJetsMCDWeighted ");
    // LOGF(info, "collision index = %d ", collision.globalIndex());
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
      // LOGF(info, "jetweight = %.8f ",jetweight);
      float pTHat = 10. / (std::pow(jetweight, 1.0 / pTHatExponent));
      // LOGF(info, "pTHat = %.8f ",pTHat);
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
        // auto thetagMCD = jetReclustering<false, false>(jet, jetSplittingsMCDTable, jetweight);
        // LOGF(info, "thetagMCD_process = %.4f", thetagMCD.value());
        registry.fill(HIST("h_jet_pt_after_leadingtrackcut_mcd"), jet.pt()); 
        registry.fill(HIST("h_jet_pt_after_leadingtrackcut_mcd_weighted"), jet.pt(), jetweight); 
        // LOGF(info, "jetweight = %.4f",jetweight);
        analyseCharged<false>(jet, tracks, jetSplittingsMCDTable, jetweight);
        
        // LOGF(info, "thetagMCD_process = %.4f", thetagMCD.value());
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
    // LOGF(info, "processChargedJetsEventWiseSubMCD");
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
      return;
    }
    // LOGF(info, "test1");
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    // LOGF(info, "test2");
    // LOGF(info, "collision index = %d ", collision.globalIndex());
    // LOGF(info, "Nombre de jets dans cet événement : %d", jets.size());
    ///////////// leading track cut /////////////
    for (auto& jet : jets){
      // LOGF(info, "test3");

      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracksSub>(jet)) {
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
  // //rajouter les cuts de jetspectra
  if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
    return;
  }
  if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
    return;
  }
  // ///////////// leading track cut /////////////
  for (auto& jet : jets){
    if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
      continue;
    }
    if (!isAcceptedJet<aod::JetTracksSub>(jet)) {
      continue;
    }
    float jetweight = jet.eventWeight();
    // LOGF(info, "jetweight = %.8f ",jetweight);
    float pTHat = 10. / (std::pow(jetweight, 1.0 / pTHatExponent));
    // LOGF(info, "pTHat = %.8f ",pTHat);
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
  //   // Si un jet contient un constituant avec un pt > au critère, on l'analyse
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
      if ((centralityMin < collision.centFT0C()) && (collision.centFT0C() < centralityMax)) {
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
      bool hasHighPtConstituent = false;
      registry.fill(HIST("h_jet_pt_initial_mcp"), jet.pt()); 
      for (auto& jetConstituent : jet.tracks_as<aod::JetParticles>()) {
        if (jetConstituent.pt() >= ptLeadingTrackCut) {
          hasHighPtConstituent = true;
          break; // Sortir de la boucle dès qu'un constituant valide est trouvé
        }
      }
      if (hasHighPtConstituent) {
        registry.fill(HIST("h_jet_pt_after_leadingtrackcut_mcp"), jet.pt()); 
        //début de analyseCharged version MCP
        jetConstituents.clear();
        for (auto& jetConstituent : jet.tracks_as<aod::JetParticles>()) {
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
  // LOGF(info, "processChargedJetsMCPWeighted ");
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
    if ((centralityMin < collision.centFT0C()) && (collision.centFT0C() < centralityMax)) {
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
    bool hasHighPtConstituent = false;
    float jetweight = jet.eventWeight();
    // LOGF(info, "jetweight = %.8f",jetweight);
    double pTHat = 10. / (std::pow(jetweight, 1.0 / pTHatExponent));
    // LOGF(info, "pTHat = %.8f ",pTHat);
    registry.fill(HIST("h_jet_pt_initial_mcp"), jet.pt());
    registry.fill(HIST("h_jet_pt_initial_mcp_weighted"), jet.pt(),jetweight);
    registry.fill(HIST("h_jet_phat_initial_mcp"), pTHat);
    registry.fill(HIST("h_jet_pthat_initial_mcp_weighted"), pTHat, jetweight); 
    for (auto& jetConstituent : jet.tracks_as<aod::JetParticles>()) {
      if (jetConstituent.pt() >= ptLeadingTrackCut) {
        hasHighPtConstituent = true;
        break; // Sortir de la boucle dès qu'un constituant valide est trouvé
      }
    }
    if (hasHighPtConstituent) {
      registry.fill(HIST("h_jet_pt_after_leadingtrackcut_mcp"), jet.pt());
      registry.fill(HIST("h_jet_pt_after_leadingtrackcut_mcp_weighted"), jet.pt(),jetweight);
      //début de analyseCharged version MCP
      jetConstituents.clear();
      for (auto& jetConstituent : jet.tracks_as<aod::JetParticles>()) {
        fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex(), static_cast<int>(JetConstituentStatus::track), pdg->Mass(jetConstituent.pdgCode()));
      }
      jetReclustering<true, false>(jet, jetSplittingsMCPTable , jetweight);
      // auto thetagMCP = jetReclustering<true, false>(jet, jetSplittingsMCPTable, jetweight);
      // LOGF(info, "thetagMCP_process = %.4f", thetagMCP.value());
      // LOGF(info, "jetMCP_process: pt = %.3f, eta = %.3f, phi = %.3f", jet.pt(), jet.eta(), jet.phi());
      //fin de analyseCharged version MCP
      // LOGF(info, "processChargedJetsMCP: weight = %.4f",jetweight);
    }
  }
}
PROCESS_SWITCH(JetSubstructureTask, processChargedJetsMCPWeighted, "charged jet substructure on MC particle level weighted", false);

void processJetsMCDMatchedMCP(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                              ChargedMCDMatchedJets const& mcdjets,
                              ChargedMCPMatchedJets const& mcpjets,
                              aod::JetTracks const& track, aod::JetParticles const&)
{
  if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
    return;
  }
  if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
    return;
  }

  for (const auto& mcdjet : mcdjets) {
    if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
      continue;
    }
    if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
      continue;
    }
    bool hasHighPtConstituent = false;
    ///////////// leading track cut /////////////
    for (auto& jetConstituent : mcdjet.tracks_as<aod::JetTracks>()) {
      if (jetConstituent.pt() >= ptLeadingTrackCut) {
        hasHighPtConstituent = true;
        break; // Sortir de la boucle dès qu'un constituant valide est trouvé
      }
    }
    if (hasHighPtConstituent) {
      fillMatchedHistograms<ChargedMCDMatchedJets::iterator, ChargedMCPMatchedJets>(mcdjet, thetagMCDVec, thetagMCPVec);
    }
  }
}
PROCESS_SWITCH(JetSubstructureTask, processJetsMCDMatchedMCP, "matched mcp and mcd jets", false);

int totalMCDjets = 0;
void processJetsMCDMatchedMCPWeighted(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                                      ChargedMCDMatchedJetsWeighted const& mcdjets,
                                      ChargedMCPMatchedJetsWeighted const&,
                                      aod::JetTracks const& tracks, aod::JetParticles const&)
{
  if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
    return;
  }
  if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
    return;
  }
  LOGF(info, "Total MCD jets in this timeframe = %d", mcdjets.size());
  totalMCDjets += mcdjets.size();
  for (const auto& mcdjet : mcdjets) {
    if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
      continue;
    }
    if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
      continue;
    }
    float jetweight = mcdjet.eventWeight();
    float pTHat = 10. / (std::pow(jetweight, 1.0 / pTHatExponent));
    if (mcdjet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    bool hasHighPtConstituent = false;
    ///////////// leading track cut /////////////
    for (auto& jetConstituent : mcdjet.tracks_as<aod::JetTracks>()) {
      if (jetConstituent.pt() >= ptLeadingTrackCut) {
        hasHighPtConstituent = true;
        break; // Sortir de la boucle dès qu'un constituant valide est trouvé
      }
    }
    if (hasHighPtConstituent) {
      // analyseCharged<false>(mcdjet, tracks, jetSplittingsMCDTable, jetweight);
      // auto thetagMCD = jetReclustering<false, false>(mcdjet, jetSplittingsMCDTable, jetweight);
      // LOGF(info, "thetagMCD = %.4f", thetagMCD.value());
      //if (doprocessChargedJetsMCD || doprocessChargedJetsMCDWeighted){ //doprocessChargedJetsEventWiseSubMCD
      // fillMatchedHistograms<ChargedMCDMatchedJetsWeighted::iterator, ChargedMCPMatchedJetsWeighted>(mcdjet,jetSplittingsMCDTable, jetSplittingsMCPTable, mcdjet.eventWeight());
      fillMatchedHistograms<ChargedMCDMatchedJetsWeighted::iterator, ChargedMCPMatchedJetsWeighted>(mcdjet, thetagMCDVec, thetagMCPVec, jetweight);


      //}
    }
  }
}
void endOfStream(o2::framework::EndOfStreamContext&) {
    LOGF(info, "====== Total MCD jets over the entire dataset = %d ======", totalMCDjets);
}
PROCESS_SWITCH(JetSubstructureTask, processJetsMCDMatchedMCPWeighted, "matched mcp and mcd jets with weighted events", false);

void processJetsMCDEventWiseMatchedMCP(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                                       ChargedMCDEventWiseMatchedtoMCD const& jetsMCDEventWise, 
                                       ChargedMCPMatchedJets const& jetsMCP,
                                       ChargedMCDMatchedJets const&,
                                       aod::JetTracksSub const& tracks)
{
  // LOGF(info, "entrering processJetsMCDEventWiseMatchedMCD");
  if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
    return;
  }
  // LOGF(info, "test3");

  if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
    return;
  }
  // LOGF(info, "test4");
  for (const auto& jetMCDEventWise : jetsMCDEventWise) {
    // LOGF(info, "test5");
    if (!jetfindingutilities::isInEtaAcceptance(jetMCDEventWise, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
      continue;
    }
    if (!isAcceptedJet<aod::JetTracksSub>(jetMCDEventWise)) {
      continue;
    }
    bool hasHighPtConstituent = false;
    ///////////// leading track cut /////////////
    for (auto& jetConstituent : jetMCDEventWise.tracks_as<aod::JetTracksSub>()) {
      if (jetConstituent.pt() >= ptLeadingTrackCut) {
        hasHighPtConstituent = true;
        break; // Sortir de la boucle dès qu'un constituant valide est trouvé
      }
    }
    if (hasHighPtConstituent) {
      // LOGF(info, "entrering fillMatchedHistogramsEventWise in processJetsMCDEventWiseMatchedMCD");
      fillMatchedHistogramsEventWise<ChargedMCDEventWiseMatchedtoMCD::iterator, ChargedMCDMatchedJets, ChargedMCPMatchedJets>(jetMCDEventWise, jetsMCP, thetagMCDEventWiseVec, thetagMCPVec);
    }
  }
}
PROCESS_SWITCH(JetSubstructureTask, processJetsMCDEventWiseMatchedMCP, "matched mcp and mcd jets eventwise", false);

void processJetsMCDEventWiseMatchedMCPWeighted(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                                       ChargedMCDEventWiseMatchedtoMCDWeighted const& jetsMCDEventWise, 
                                       ChargedMCPMatchedJetsWeighted const& jetsMCP,
                                       ChargedMCDMatchedJetsWeighted const&,
                                       aod::JetTracksSub const& tracks)
{
  // LOGF(info, "entrering processJetsMCDEventWiseMatchedMCDWeighted");
  if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
    return;
  }
  if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
    return;
  }
  for (const auto& jetMCDEventWise : jetsMCDEventWise) {
    if (!jetfindingutilities::isInEtaAcceptance(jetMCDEventWise, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
      continue;
    }
    if (!isAcceptedJet<aod::JetTracksSub>(jetMCDEventWise)) {
      continue;
    }
    float jetweight = jetMCDEventWise.eventWeight();
    float pTHat = 10. / (std::pow(jetweight, 1.0 / pTHatExponent));
    if (jetMCDEventWise.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    bool hasHighPtConstituent = false;
    ///////////// leading track cut /////////////
    for (auto& jetConstituent : jetMCDEventWise.tracks_as<aod::JetTracksSub>()) {
      if (jetConstituent.pt() >= ptLeadingTrackCut) {
        hasHighPtConstituent = true;
        break; // Sortir de la boucle dès qu'un constituant valide est trouvé
      }
    }
    if (hasHighPtConstituent) {
      // LOGF(info, "entrering fillMatchedHistogramsEventWise in processJetsMCDEventWiseMatchedMCD");
      fillMatchedHistogramsEventWise<ChargedMCDEventWiseMatchedtoMCDWeighted::iterator, ChargedMCDMatchedJetsWeighted, ChargedMCPMatchedJetsWeighted>(jetMCDEventWise, jetsMCP, thetagMCDEventWiseVec, thetagMCPVec, jetweight);
    }
  }
}
PROCESS_SWITCH(JetSubstructureTask, processJetsMCDEventWiseMatchedMCPWeighted, "matched mcp and mcd jets eventwise weighted", false);
};




WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  return WorkflowSpec{adaptAnalysisTask<JetSubstructureTask>(
    cfgc, TaskName{"jet-substructure-softdrop"})};
}
