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
  Produces<aod::CJetSSs> jetSubstructureDataTable;
  Produces<aod::CMCDJetSSs> jetSubstructureMCDTable;
  Produces<aod::CMCPJetSSs> jetSubstructureMCPTable;
  Produces<aod::CEWSJetSSs> jetSubstructureDataSubTable;

  Produces<aod::ChargedSPs> jetSplittingsDataTable;
  Produces<aod::ChargedMCDetectorLevelSPs> jetSplittingsMCDTable;
  Produces<aod::ChargedMCParticleLevelSPs> jetSplittingsMCPTable;
  Produces<aod::ChargedEventWiseSubtractedSPs> jetSplittingsDataSubTable;

  Produces<aod::ChargedPRs> jetPairsDataTable;
  Produces<aod::ChargedMCDetectorLevelPRs> jetPairsMCDTable;
  Produces<aod::ChargedMCParticleLevelPRs> jetPairsMCPTable;
  Produces<aod::ChargedEventWiseSubtractedPRs> jetPairsDataSubTable;

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

  void init(InitContext const&)
  {
    registry.add("h2_jet_pt_jet_zg", ";#it{p}_{T,jet} (GeV/#it{c});#it{z}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_rg", ";#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_nsd", ";#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});

    registry.add("h2_jet_pt_part_jet_zg_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{z}_{g}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_part_jet_rg_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{R}_{g}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_part_jet_nsd_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{n}_{SD}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});

    registry.add("h2_jet_pt_jet_zg_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{z}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_rg_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_nsd_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});

    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;
  }

  Preslice<aod::JetTracks> TracksPerCollision = aod::jtrack::collisionId;
  Preslice<aod::JetTracksSub> TracksPerCollisionDataSub = aod::bkgcharged::collisionId;
  Preslice<aod::JetParticles> ParticlesPerMcCollision = aod::jmcparticle::mcCollisionId;



  // Filter trackCuts = (aod::jtrack::pt >= trackQAPtMin && aod::jtrack::pt < trackQAPtMax && aod::jtrack::eta > trackQAEtaMin && aod::jtrack::eta < trackQAEtaMax);
  // Filter particleCuts = (aod::jmcparticle::pt >= trackQAPtMin && aod::jmcparticle::pt < trackQAPtMax && aod::jmcparticle::eta > trackQAEtaMin && aod::jmcparticle::eta < trackQAEtaMax);
  // Filter collisionFilter = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);

  template <bool isMCP, bool isSubtracted, typename T, typename U>
  void jetReclustering(T const& jet, U& splittingTable)
  {
    LOGF(info, " Entering jetReclustering " );
    energyMotherVec.clear(); //to be sure its empty before filling
    ptLeadingVec.clear();
    ptSubLeadingVec.clear();
    thetaVec.clear();
    jetReclustered.clear();
    fastjet::ClusterSequenceArea clusterSeq(jetReclusterer.findJets(jetConstituents, jetReclustered)); // jetConstituents is a vector of fastjet::PseudoJet containing (tracks/constituants)
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
        if (!softDropped) {//if the splitting hasent been already softdropped softdrop=false
          zg = z;
          rg = theta;
          if constexpr (!isSubtracted && !isMCP) {
            LOGF(info, " Entering if statement for histograms :" );
            registry.fill(HIST("h2_jet_pt_jet_zg"), jet.pt(), zg);
            registry.fill(HIST("h2_jet_pt_jet_rg"), jet.pt(), rg);
          }
          if constexpr (!isSubtracted && isMCP) {
            registry.fill(HIST("h2_jet_pt_part_jet_zg_part"), jet.pt(), zg);
            registry.fill(HIST("h2_jet_pt_part_jet_rg_part"), jet.pt(), rg);
          }
          if constexpr (isSubtracted && !isMCP) {
            registry.fill(HIST("h2_jet_pt_jet_zg_eventwiseconstituentsubtracted"), jet.pt(), zg);
            registry.fill(HIST("h2_jet_pt_jet_rg_eventwiseconstituentsubtracted"), jet.pt(), rg);
          }
          softDropped = true;//mark the splitting as true to avoid filling it again
        }
        nsd++;//step up the number of iterations
      }
      daughterSubJet = parentSubJet1; //following with the hardest branch
    }
    if constexpr (!isSubtracted && !isMCP) {
      LOGF(info, " Entering if statement for histograms: " );
      registry.fill(HIST("h2_jet_pt_jet_nsd"), jet.pt(), nsd);
    }
    if constexpr (!isSubtracted && isMCP) {
      registry.fill(HIST("h2_jet_pt_part_jet_nsd_part"), jet.pt(), nsd);
    }
    if constexpr (isSubtracted && !isMCP) {
      registry.fill(HIST("h2_jet_pt_jet_nsd_eventwiseconstituentsubtracted"), jet.pt(), nsd);
    }
  }

  template <typename T, typename U>
  void jetSubstructureSimple(T const& jet, U const& /*tracks*/)
  {
    angularity = 0.0;
    leadingConstituentPt = 0.0;
    for (auto& constituent : jet.template tracks_as<U>()) {
      if (constituent.pt() >= leadingConstituentPt) {
        leadingConstituentPt = constituent.pt();
      }
      angularity += std::pow(constituent.pt(), kappa) * std::pow(jetutilities::deltaR(jet, constituent), alpha);
    }
    angularity /= (jet.pt() * (jet.r() / 100.f));
  }

  template <bool isSubtracted, typename T, typename U, typename V, typename M, typename N>
  void analyseCharged(T const& jet, U const& tracks, V const& trackSlicer, M& outputTable, N& splittingTable)
  {
    LOGF(info, " Entering analyseCharged " );
    jetConstituents.clear();
    for (auto& jetConstituent : jet.template tracks_as<U>()) {
      fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
    }
    nSub = jetsubstructureutilities::getNSubjettiness(jet, tracks, tracks, tracks, 2, fastjet::contrib::CA_Axes(), true, zCut, beta);
    jetReclustering<false, isSubtracted>(jet, splittingTable);
    jetSubstructureSimple(jet, tracks);
    outputTable(energyMotherVec, ptLeadingVec, ptSubLeadingVec, thetaVec, nSub[0], nSub[1], nSub[2], pairJetPtVec, pairJetEnergyVec, pairJetThetaVec, pairJetPerpCone1PtVec, pairJetPerpCone1EnergyVec, pairJetPerpCone1ThetaVec, pairPerpCone1PerpCone1PtVec, pairPerpCone1PerpCone1EnergyVec, pairPerpCone1PerpCone1ThetaVec, pairPerpCone1PerpCone2PtVec, pairPerpCone1PerpCone2EnergyVec, pairPerpCone1PerpCone2ThetaVec, angularity, leadingConstituentPt, perpConeRho);
  }

  void processDummy(aod::JetTracks const&)
  {
  }
  PROCESS_SWITCH(JetSubstructureTask, processDummy, "Dummy process function turned on by default", true);

  void processChargedJetsData(soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>::iterator const& jet,
                              aod::JetCollisions const& collisions,
                              aod::JetTracks const& tracks)
  {
  //   LOGF(info, " Entering processChargedJetsData 1 " );
  //   bool hasHighPtConstituent = false;
  //   for (auto& jetConstituent : jet.tracks_as<aod::JetTracks>()) {
  //     if (jetConstituent.pt() >= 5.0f) {
  //       hasHighPtConstituent = true;
  //       break; 
  //     }
  //   }
  //  }
  // void processChargedJetsData(soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>::iterator const& jet,
  //                             soa::Filtered<aod::JetCollisions> const& collisions,
  //                             soa::Filtered<aod::JetTracks> const& tracks)
  // {
    
    LOGF(info, " Entering processChargedJetsData 1 " );

    ///////////// leading track cut try : (because filter doesnt work)

      bool hasHighPtConstituent = false;
      LOGF(info, " Entering processChargedJetsData 2 " );
      for (auto& jetConstituent : jet.tracks_as<aod::JetTracks>()) {
        LOGF(info, " Entering processChargedJetsData 3" );
        if (jetConstituent.pt() >= 5.0f) {
          LOGF(info, " Entering processChargedJetsData 4" );
          hasHighPtConstituent = true;
          break; // Sortir de la boucle dès qu'un constituant valide est trouvé
        }
      }
      LOGF(info, " Entering processChargedJetsData 5 " );
      // Si un jet contient un constituant avec un pt élevé, on l'analyse
      if (hasHighPtConstituent) {
        LOGF(info, " Entering if statement processChargedJetsData " );
        analyseCharged<false>(jet, tracks, TracksPerCollision, jetSubstructureDataTable, jetSplittingsDataTable);
      }
    
    /////////////// track selection try: (because filter doesnt work)

    // std::vector<int32_t> filteredTracks;
    // for (const auto& track : tracks) {
    //    if (track.pt() >= trackQAPtMin && track.pt() < trackQAPtMax &&
    //      track.eta() >= trackQAEtaMin && track.eta() < trackQAEtaMax) {
    //      filteredTracks.push_back(track);
    //    }
    //  }
    // if (filteredTracks.empty()) {
    //   return;
    // } //au lieu de mettre tracks dans analyseCharged on met filteredTracks 
    // analyseCharged<false>(jet, tracks, TracksPerCollision, jetSubstructureDataTable, jetSplittingsDataTable);
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsData, "charged jet substructure", false);

  void processChargedJetsEventWiseSubData(soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents>::iterator const& jet,
                                          aod::JetTracksSub const& tracks)
  {
    analyseCharged<true>(jet, tracks, TracksPerCollisionDataSub, jetSubstructureDataSubTable, jetSplittingsDataSubTable);
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsEventWiseSubData, "eventwise-constituent subtracted charged jet substructure", false);

  void processChargedJetsMCD(typename soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>::iterator const& jet,
                             aod::JetTracks const& tracks)
  {
    analyseCharged<false>(jet, tracks, TracksPerCollision, jetSubstructureMCDTable, jetSplittingsMCDTable);
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsMCD, "charged jet substructure", false);

  void processChargedJetsMCP(typename soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>::iterator const& jet,
                             aod::JetParticles const& particles)
  {
    jetConstituents.clear();
    for (auto& jetConstituent : jet.template tracks_as<aod::JetParticles>()) {
      fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex(), static_cast<int>(JetConstituentStatus::track), pdg->Mass(jetConstituent.pdgCode()));
    }
    nSub = jetsubstructureutilities::getNSubjettiness(jet, particles, particles, particles, 2, fastjet::contrib::CA_Axes(), true, zCut, beta);
    jetReclustering<true, false>(jet, jetSplittingsMCPTable);
    
    jetSubstructureSimple(jet, particles);
    jetSubstructureMCPTable(energyMotherVec, ptLeadingVec, ptSubLeadingVec, thetaVec, nSub[0], nSub[1], nSub[2], pairJetPtVec, pairJetEnergyVec, pairJetThetaVec, pairJetPerpCone1PtVec, pairJetPerpCone1EnergyVec, pairJetPerpCone1ThetaVec, pairPerpCone1PerpCone1PtVec, pairPerpCone1PerpCone1EnergyVec, pairPerpCone1PerpCone1ThetaVec, pairPerpCone1PerpCone2PtVec, pairPerpCone1PerpCone2EnergyVec, pairPerpCone1PerpCone2ThetaVec, angularity, leadingConstituentPt, perpConeRho);
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsMCP, "charged jet substructure on MC particle level", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  return WorkflowSpec{adaptAnalysisTask<JetSubstructureTask>(
    cfgc, TaskName{"jet-substructure-softdrop"})};
}