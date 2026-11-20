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

/// \file embeddingResponsePurity.cxx
/// \brief Response matrix (pT, theta_g, z_g), purity, efficiency, miss and fake for pp MC embedded in Pb-Pb data.
///        O2Physics version of the standalone macro fillResponsePurity.sh. It runs on the output of the standard
///        embedding chain : je-jet-deriveddata-embedding-producer, estimator-rho, subtractor-eventwiseconstituent,
///        jet finders (data / MCD / MCP charged) and the geometric jet matching tasks :
///          hybrid (event-wise subtracted) jet -> MCD jet : jet-matching-mcd-sub-ch
///          MCD jet -> MCP jet                          : jet-matching-mc-ch
///          hybrid jet -> unsubtracted hybrid jet       : jet-matching-sub-ch (used for the pp pT fraction)
///
///        pp fraction of a hybrid jet = pT of the embedded (MCD) constituents of its geometrically matched
///        unsubtracted jet / pT of the subtracted jet. The subtracted tracks (JetTracksSub) lose the embeddedTrack
///        bit, so the unsubtracted jet is the only place where the pp constituents can be identified.
///
///        Fakes are split by origin : pp fakes (pp fraction >= fakePPFractionMin) keep the MC weight and are
///        normalised like the signal, Pb-Pb fakes are filled with weight 1 and normalised per event (h_norm bin 2).
///
/// \author Louise Millot

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH3.h>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/contrib/SoftDrop.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
struct JetInfo {
  int64_t id;
  double pt;
  double eta;
  double phi;
  double thetaG;
  double zG;
  float weight;
  double ppFraction; // hybrid jets only
  int64_t matchedId; // hybrid : geometrically matched MCD jet, MCD : geometrically matched MCP jet, -1 if none
};

constexpr double PionMass = 0.13957;
constexpr double ZgUntaggedFill = 0.15;
constexpr double ThetaGUntaggedFill = -0.005;

enum Species { SpeciesMCP = 0,
               SpeciesMCD,
               SpeciesHybrid,
               SpeciesData,
               NSpecies };
enum Stage { BeforeVeto = 0,
             AfterVeto,
             NStages };

bool isEmbeddedTrack(uint8_t trackSel)
{
  return trackSel & (1ULL << jetderiveddatautilities::JTrackSel::embeddedTrack);
}
} // namespace

struct EmbeddingResponsePurityTask {
  using SubJets = soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents>;
  using UnsubJets = soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>;
  using McdJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>;
  using McpJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>;

  HistogramRegistry registry{"registry"};

  // event and track selection (keep them identical to the jet finders)
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", true, "skip MB gap events (jet finder default)"};
  Configurable<bool> applyRCTSelections{"applyRCTSelections", true, "apply RCT selections"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0, "|z| cut on the Pb-Pb collision"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum FT0M centrality of the Pb-Pb collision (embedding and data)"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum FT0M centrality of the Pb-Pb collision (embedding and data)"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT (Pb-Pb jets of the overlap veto, track QA)"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track |eta| (Pb-Pb jets of the overlap veto, track QA)"};
  Configurable<bool> useWeights{"useWeights", true, "fill with the MC collision weight"};

  // jets
  Configurable<float> jetR{"jetR", 0.2, "jet radius (must be produced by the jet finders)"};
  Configurable<float> jetEtaMax{"jetEtaMax", 0.7, "|eta| cut on MCP (efficiency) and hybrid jets"};
  Configurable<float> truthJetPtMin{"truthJetPtMin", 5.0, "pT cut on MCP and MCD jets"};
  Configurable<float> leadingTrackPtMax{"leadingTrackPtMax", 100.0, "reject jets with a constituent above this pT"};
  Configurable<float> zCut{"zCut", 0.2, "SoftDrop z_cut"};
  Configurable<float> beta{"beta", 0.0, "SoftDrop beta"};

  // matching (the DeltaR cut itself is the one of the matching tasks : maxMatchingDistancePerJetR, 0.12 for R = 0.2)
  Configurable<float> ppFractionMin{"ppFractionMin", 0.5, "min pp pT fraction of the hybrid jet for a hybrid-MCD match"};
  Configurable<bool> requireMcdMcpPtMatch{"requireMcdMcpPtMatch", false, "require the MCD-MCP pair to be pT-matched on top of the geometric match"};
  Configurable<float> fakePPFractionMin{"fakePPFractionMin", 0.5, "fake with pp fraction >= this value -> pp fake (weighted), otherwise Pb-Pb fake (weight 1)"};

  // dead zone in phi : the whole jet is rejected if one of its constituents falls in [phiVetoMin, phiVetoMax]
  Configurable<bool> applyPhiVeto{"applyPhiVeto", true, "apply the phi dead-zone veto"};
  Configurable<float> phiVetoMin{"phiVetoMin", 2.6, "phi veto lower edge"};
  Configurable<float> phiVetoMax{"phiVetoMax", 3.6, "phi veto upper edge"};

  // option 2 of the macro : reject the event if an MCP jet is close to the leading / subleading Pb-Pb jet.
  // The pairing is fixed by the BC, so a rejected event cannot be retried with another Pb-Pb event.
  Configurable<float> pbOverlapDR{"pbOverlapDR", -1.0, "reject the event if an MCP jet is within this DeltaR of the Pb-Pb leading/subleading jet, <= 0 disables"};
  Configurable<float> pbOverlapPtMin{"pbOverlapPtMin", 0.0, "min raw pT of the Pb-Pb jets used for the overlap veto"};

  // outliers : the whole event is rejected if its leading jet is above factor x pT hat
  Configurable<float> outlierFactorMCP{"outlierFactorMCP", 4.0, "MCP leading jet pT / pT hat above which the event is rejected"};
  Configurable<float> outlierFactorMCD{"outlierFactorMCD", 8.0, "MCD leading jet pT / pT hat above which the event is rejected"};
  Configurable<bool> useMcPtHard{"useMcPtHard", false, "take pT hat from the McCollision table instead of the weight formula"};
  Configurable<float> ptHatRef{"ptHatRef", 10.0, "pT hat = ptHatRef * w^(-1/ptHatBiasPow)"};
  Configurable<float> ptHatBiasPow{"ptHatBiasPow", 4.0, "pT hat = ptHatRef * w^(-1/ptHatBiasPow)"};

  Preslice<SubJets> subJetsPerCollision = aod::jet::collisionId;
  Preslice<McdJets> mcdJetsPerCollision = aod::jet::collisionId;
  Preslice<McpJets> mcpJetsPerMcCollision = aod::jet::mcCollisionId;
  Preslice<aod::JetTracks> tracksPerCollision = aod::jtrack::collisionId;
  Preslice<aod::JetParticles> particlesPerMcCollision = aod::jmcparticle::mcCollisionId;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;
  int jetRInt = 20;

  fastjet::JetDefinition jetDefPb;
  fastjet::JetDefinition reclusterDef;
  std::unique_ptr<fastjet::contrib::SoftDrop> softDrop;

  // phi control histograms (pT, eta, phi), before and after the phi veto, for jets and for their constituents
  std::array<std::array<std::shared_ptr<TH3>, NStages>, NSpecies> hJetPtEtaPhi;
  std::array<std::array<std::shared_ptr<TH3>, NStages>, NSpecies> hTrackPtEtaPhi;

  void init(InitContext const&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(eventSelections);
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(trackSelections);
    jetRInt = static_cast<int>(std::round(jetR.value * 100.));

    jetDefPb = fastjet::JetDefinition(fastjet::antikt_algorithm, static_cast<double>(jetR));
    reclusterDef = fastjet::JetDefinition(fastjet::cambridge_algorithm, 1.0); // all constituents of an R = 0.2 jet end in one C/A jet
    softDrop = std::make_unique<fastjet::contrib::SoftDrop>(static_cast<double>(beta), static_cast<double>(zCut));

    const AxisSpec ptAxis{200, 0., 200., "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thetaAxis{101, -0.01, 1., "#theta_{g}"};
    const AxisSpec zgAxis{std::vector<double>{0.10, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50}, "#it{z}_{g}"};
    const AxisSpec ptRangeAxis{std::vector<double>{0.15, 10., 20., 40., 60., 90., 115., 200.}, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec etaAxis{100, -1., 1., "#eta"};
    const AxisSpec phiAxis{100, 0., o2::constants::math::TwoPI, "#varphi"};
    const AxisSpec fracAxis{110, 0., 1.1, "pp #it{p}_{T} fraction"};

    registry.add("h_eventSelection", "event selection;;events", HistType::kTH1D, {{6, 0., 6.}});
    auto hEvSel = registry.get<TH1>(HIST("h_eventSelection"));
    const std::array<const char*, 6> evSelLabels{"all", "selected", "no MC collision", "Pb-Pb overlap veto", "outlier", "embedded"};
    for (std::size_t i = 0; i < evSelLabels.size(); i++) {
      hEvSel->GetXaxis()->SetBinLabel(i + 1, evSelLabels[i]);
    }
    // bin 1 = sum w of the embedded events (pp normalisation), bin 2 = N embedded events (Pb-Pb normalisation)
    registry.add("h_norm", "1: sum w (embedded) ; 2: N_{evt} hybrid", HistType::kTH1D, {{2, 0.5, 2.5}});
    // QA only : the MC is simulated at the vertex of the data collision, so this should be a peak at 0
    registry.add("h_collision_zdiff", ";z_{PbPb}-z_{MC} (cm);events", HistType::kTH1F, {{200, -10., 10.}});

    registry.add("h_hybrid_jets", "hybrid jets;#it{p}_{T};counts", HistType::kTH1F, {ptAxis}, true);
    registry.add("h_mcp_jets", "MCP jets;#it{p}_{T};counts", HistType::kTH1F, {ptAxis}, true);
    registry.add("h_nTracks", "N_{tracks} pp per collision;N_{tracks};counts", HistType::kTH1F, {{200, 0., 2000.}});
    registry.add("h_nParticles", "N_{particles} per MC collision;N_{particles};counts", HistType::kTH1F, {{200, 0., 5000.}});
    registry.add("h_track_pt", "hybrid track #it{p}_{T} (unsubtracted);#it{p}_{T};counts", HistType::kTH1F, {ptAxis}, true);
    registry.add("h_rho", "#rho (estimator-rho, hybrid event);#rho;events", HistType::kTH1F, {{100, 0., 500.}});
    registry.add("h2_rho_ntracks", "#rho vs N_{tracks} Pb-Pb;N_{tracks};#rho (GeV/area)", HistType::kTH2F, {{350, 0., 3500.}, {220, 0., 220.}});
    registry.add("h_pbpb_lead_pt", "Pb-Pb leading jet (raw #it{p}_{T});#it{p}_{T};events", HistType::kTH1F, {ptAxis});
    registry.add("h_pbpb_sublead_pt", "Pb-Pb subleading jet (raw #it{p}_{T});#it{p}_{T};events", HistType::kTH1F, {ptAxis});

    registry.add("h_response_pt", "pT MCP vs pT hybrid;#it{p}_{T} hybrid;#it{p}_{T} MCP", HistType::kTH2F, {ptAxis, ptAxis}, true);
    registry.add("h_response_theta", "#theta MCP vs #theta hybrid;#theta_{hybrid};#theta_{MCP}", HistType::kTH2F, {{100, 0., 1.}, {100, 0., 1.}}, true);
    registry.add("h_sd_theta_hybrid", "hybrid jets #theta_{g};#theta_{g};counts", HistType::kTH1F, {thetaAxis}, true);
    registry.add("h_sd_z_hybrid", "hybrid jets z_{g};z_{g};counts", HistType::kTH1F, {zgAxis}, true);
    registry.add("h_sd_theta_mcp", "MCP jets #theta_{g};#theta_{g};counts", HistType::kTH1F, {thetaAxis}, true);
    registry.add("h_sd_z_mcp", "MCP jets z_{g};z_{g};counts", HistType::kTH1F, {zgAxis}, true);

    // efficiency and purity (1D, ratios to be done offline)
    registry.add("h_jet_mcp_pt", "EFFICIENCY denominator : MCP jets (eta cut);#it{p}_{T} gen;counts", HistType::kTH1F, {ptAxis}, true);
    registry.add("h_jet_matched_pt", "EFFICIENCY numerator : MCP matched to hybrid;#it{p}_{T} gen;counts", HistType::kTH1F, {ptAxis}, true);
    registry.add("h_step0_hybrid_pt", "PURITY denominator : hybrid (eta cut);#it{p}_{T} hybrid;counts", HistType::kTH1F, {ptAxis}, true);
    registry.add("h_step1_hybmcd_pt", "hybrid geo-matched to MCD;#it{p}_{T} hybrid;counts", HistType::kTH1F, {ptAxis}, true);
    registry.add("h_step2_frac_pt", "hybrid matched to MCD passing the pp fraction cut;#it{p}_{T} hybrid;counts", HistType::kTH1F, {ptAxis}, true);
    registry.add("h_step3_mcdmcp_pt", "PURITY numerator : hybrid matched up to MCP;#it{p}_{T} hybrid;counts", HistType::kTH1F, {ptAxis}, true);
    registry.add("h_deltaPt", "#delta p_{T} = p_{T}^{hybrid} - p_{T}^{MCD};#delta p_{T};counts", HistType::kTH1F, {{100, -50., 50.}}, true);

    // response matrices
    registry.add("hResponse4D", "response;p_{T}^{truth};p_{T}^{reco};#theta^{truth};#theta^{reco}", HistType::kTHnSparseD, {ptAxis, ptAxis, thetaAxis, thetaAxis}, true);
    registry.add("hResponse6D", "response 6D;p_{T}^{truth};p_{T}^{reco};#theta^{truth};#theta^{reco};z_{g}^{truth};z_{g}^{reco}", HistType::kTHnSparseD, {ptAxis, ptAxis, thetaAxis, thetaAxis, zgAxis, zgAxis}, true);
    registry.add("hResponse6D_cnt", "response 6D (unweighted);p_{T}^{truth};p_{T}^{reco};#theta^{truth};#theta^{reco};z_{g}^{truth};z_{g}^{reco}", HistType::kTHnSparseD, {ptAxis, ptAxis, thetaAxis, thetaAxis, zgAxis, zgAxis}, true);

    // miss / fake and base histograms (pT, theta_g, z_g), weighted and _cnt (unweighted)
    //   purity axes = hybrid values, efficiency axes = MCP values
    //   hFake           : all fakes with the MC weight (as in the macro)
    //   hFake_pp        : fakes of pp origin, MC weight -> normalise like the signal (sum w)
    //   hFake_PbPb      : fakes of Pb-Pb origin, weight 1 -> normalise per event (N_evt)
    //   hFake_noJetColl : fakes of the events without MCD jet or without MCP jet
    for (const std::string name : {"hMiss", "hFake", "hFake_pp", "hFake_PbPb", "hPurNum", "hPurDen", "hEffNum", "hEffDen", "hFake_noJetColl"}) {
      registry.add(name, (name + ";p_{T};#theta_{g};z_{g}").c_str(), HistType::kTH3D, {ptAxis, thetaAxis, zgAxis}, true);
      registry.add(name + "_cnt", (name + " (unweighted);p_{T};#theta_{g};z_{g}").c_str(), HistType::kTH3D, {ptAxis, thetaAxis, zgAxis}, true);
    }

    // outliers
    registry.add("h_ptHatW", "#hat{p}_{T};#hat{p}_{T} (GeV);collisions", HistType::kTH1F, {{300, 0., 300.}});
    registry.add("h2_ptHatW_vs_leadMCP", "#hat{p}_{T} vs leading MCP jet;p_{T}^{lead,MCP};#hat{p}_{T}", HistType::kTH2F, {ptAxis, {300, 0., 300.}});
    registry.add("h_outlier_mcp_pt", "MCP leading jet of rejected events;p_{T};weighted events", HistType::kTH1F, {ptAxis}, true);
    registry.add("h_outlier_mcd_pt", "MCD leading jet of rejected events;p_{T};weighted events", HistType::kTH1F, {ptAxis}, true);
    for (const std::string level : {"MCP", "MCD", "Hybrid"}) {
      registry.add("h_cnt_" + level + "_before_outlier", (level + " jets before outlier cut;p_{T};counts (weighted)").c_str(), HistType::kTH1F, {ptAxis}, true);
      registry.add("h_cnt_" + level + "_after_outlier", (level + " jets after outlier cut;p_{T};counts (weighted)").c_str(), HistType::kTH1F, {ptAxis}, true);
    }

    // weight and pp fraction diagnostics, fakes vs matched
    registry.add("h2_w_pt_fake", "fakes : weight;p_{T}^{hybrid};log_{10}(w)", HistType::kTH2F, {ptAxis, {110, -9., 2.}});
    registry.add("h2_w_pt_matched", "matched : weight;p_{T}^{hybrid};log_{10}(w)", HistType::kTH2F, {ptAxis, {110, -9., 2.}});
    registry.add("h2_ptHat_pt_fake", "fakes (weighted) : #hat{p}_{T};p_{T}^{hybrid};#hat{p}_{T}", HistType::kTH2F, {ptAxis, {150, 0., 150.}}, true);
    registry.add("h2_ptHat_pt_matched", "matched (weighted) : #hat{p}_{T};p_{T}^{hybrid};#hat{p}_{T}", HistType::kTH2F, {ptAxis, {150, 0., 150.}}, true);
    registry.add("h2_ptHatCnt_pt_fake", "fakes (unweighted) : #hat{p}_{T};p_{T}^{hybrid};#hat{p}_{T}", HistType::kTH2F, {ptAxis, {150, 0., 150.}});
    registry.add("h2_ptHatCnt_pt_matched", "matched (unweighted) : #hat{p}_{T};p_{T}^{hybrid};#hat{p}_{T}", HistType::kTH2F, {ptAxis, {150, 0., 150.}});
    // to validate the pp / Pb-Pb fake threshold (fakePPFractionMin)
    registry.add("h2_fracPP_pt_fake", "fakes (unweighted) : pp fraction;p_{T}^{hybrid};pp #it{p}_{T} fraction", HistType::kTH2F, {ptAxis, fracAxis});
    registry.add("h2_fracPP_pt_matched", "matched (unweighted) : pp fraction;p_{T}^{hybrid};pp #it{p}_{T} fraction", HistType::kTH2F, {ptAxis, fracAxis});

    // phi control : the 7 pT slices of the macro are the bins of the pT axis (jet pT for jets, track pT for tracks)
    // Pb-Pb data spectrum to unfold (processData), same names as fillPbPbLeadCut.sh
    //   hPbPbPtThetaZg_all : all jets, untagged jets at (theta_g, z_g) = (-0.005, 0.15) -> input of the unfolding
    //   hPbPbPtThetaZg     : tagged jets only (theta_g > 0)
    registry.add("hPbPbPtThetaZg_all", "Pb-Pb jets;p_{T} (GeV/#it{c});#theta_{g};z_{g}", HistType::kTH3D, {ptAxis, thetaAxis, zgAxis}, true);
    registry.add("hPbPbPtThetaZg", "Pb-Pb jets (tagged);p_{T} (GeV/#it{c});#theta_{g};z_{g}", HistType::kTH3D, {ptAxis, {100, 0., 1., "#theta_{g}"}, zgAxis}, true);
    registry.add("hPbPbPtTheta_all", "Pb-Pb jets;p_{T} (GeV/#it{c});#theta_{g}", HistType::kTH2D, {ptAxis, thetaAxis}, true);
    registry.add("hPbPbPtTheta", "Pb-Pb jets (tagged);p_{T} (GeV/#it{c});#theta_{g}", HistType::kTH2D, {ptAxis, {100, 0., 1., "#theta_{g}"}}, true);
    registry.add("h_data_nEvents", "Pb-Pb data : selected events;;events", HistType::kTH1D, {{1, 0., 1.}});
    registry.add("h_data_rho", "Pb-Pb data : #rho (estimator-rho);#rho;events", HistType::kTH1F, {{100, 0., 500.}});
    registry.add("h_data_centrality", "Pb-Pb data : centrality;FT0M centrality (%);events", HistType::kTH1F, {{100, 0., 100.}});
    registry.add("h_data_leadcut_pt", "Pb-Pb data : jets rejected by the leading track cut;p_{T} jet (GeV/#it{c});counts", HistType::kTH1F, {ptAxis});
    // raw tracks (before subtraction), track pT slices of fillPbPbLeadCut.sh as bins
    const AxisSpec ptTrackRangeAxis{std::vector<double>{0.15, 2., 5., 7., 20., 30., 50., 80., 110., 150.}, "#it{p}_{T}^{track} (GeV/#it{c})"};
    registry.add("h_data_track_pt_eta_phi", "Pb-Pb data : raw tracks;p_{T}^{track};#eta;#varphi", HistType::kTH3F, {ptTrackRangeAxis, etaAxis, phiAxis});

    const std::array<std::string, NSpecies> speciesNames{"mcp", "mcd", "hybrid", "data"};
    const std::array<std::string, NStages> stageNames{"before", "after"};
    for (int isp = 0; isp < NSpecies; isp++) {
      for (int is = 0; is < NStages; is++) {
        const std::string suffix = speciesNames[isp] + "_" + stageNames[is];
        hJetPtEtaPhi[isp][is] = registry.add<TH3>("h_jet_pt_eta_phi_" + suffix, (speciesNames[isp] + " jets (" + stageNames[is] + " phi veto);p_{T}^{jet};#eta;#varphi").c_str(), HistType::kTH3F, {ptRangeAxis, etaAxis, phiAxis}, true);
        hTrackPtEtaPhi[isp][is] = registry.add<TH3>("h_track_pt_eta_phi_" + suffix, (speciesNames[isp] + " jet constituents (" + stageNames[is] + " phi veto);p_{T}^{track};#eta;#varphi").c_str(), HistType::kTH3F, {ptRangeAxis, etaAxis, phiAxis}, true);
      }
    }
  }

  double ptHatFromWeight(double w) const
  {
    return (w > 0.) ? ptHatRef.value * std::pow(w, -1.0 / ptHatBiasPow.value) : 0.;
  }

  template <typename T>
  double leadingConstituentPt(T const& constituents) const
  {
    double ptLeading = 0.;
    for (auto const& c : constituents) {
      ptLeading = std::max(ptLeading, static_cast<double>(c.pt()));
    }
    return ptLeading;
  }

  template <typename T>
  bool hasConstituentInPhiVeto(T const& constituents) const
  {
    if (!applyPhiVeto.value) {
      return false;
    }
    for (auto const& c : constituents) {
      if (c.phi() >= phiVetoMin.value && c.phi() <= phiVetoMax.value) {
        return true;
      }
    }
    return false;
  }

  template <typename J, typename T>
  void fillPhiControl(int species, J const& jet, T const& constituents, double w, int stage)
  {
    hJetPtEtaPhi[species][stage]->Fill(jet.pt(), jet.eta(), jet.phi(), w);
    for (auto const& c : constituents) {
      hTrackPtEtaPhi[species][stage]->Fill(c.pt(), c.eta(), c.phi(), w);
    }
  }

  // Pb-Pb collision selection, identical for the embedding and for the data
  template <typename C>
  bool selectPbPbCollision(C const& collision) const
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents.value, applyRCTSelections.value) || std::abs(collision.posZ()) > vertexZCut.value) {
      return false;
    }
    return collision.centFT0M() >= centralityMin.value && collision.centFT0M() < centralityMax.value;
  }

  // cuts on the reconstructed (event-wise subtracted) jets, identical for the hybrid jets of the embedding and
  // for the data jets : leading constituent, then phi dead-zone veto with control histograms before / after.
  // The |eta| cut is applied by the caller, after the inclusive histograms (as in the macros)
  template <typename J, typename T>
  bool passRecoJetCuts(int species, J const& jet, T const& constituents, double w)
  {
    if (leadingConstituentPt(constituents) >= leadingTrackPtMax) {
      return false;
    }
    fillPhiControl(species, jet, constituents, w, BeforeVeto);
    if (hasConstituentInPhiVeto(constituents)) {
      return false;
    }
    fillPhiControl(species, jet, constituents, w, AfterVeto);
    return true;
  }

  // theta_g = R_g / R, untagged jets -> (theta_g, z_g) = (-0.005, 0.15)
  template <bool IsParticle, typename T>
  bool computeSoftDrop(T const& constituents, double& thetaG, double& zG)
  {
    thetaG = ThetaGUntaggedFill;
    zG = ZgUntaggedFill;
    std::vector<fastjet::PseudoJet> parts;
    for (auto const& c : constituents) {
      double energy;
      if constexpr (IsParticle) {
        energy = c.e();
      } else {
        energy = std::sqrt(c.p() * c.p() + PionMass * PionMass);
      }
      parts.emplace_back(c.px(), c.py(), c.pz(), energy);
    }
    if (parts.empty()) {
      return false;
    }
    fastjet::ClusterSequence cs(parts, reclusterDef);
    auto jets = fastjet::sorted_by_pt(cs.inclusive_jets());
    if (jets.empty()) {
      return false;
    }
    fastjet::PseudoJet sdJet = (*softDrop)(jets[0]);
    if (sdJet.has_structure_of<fastjet::contrib::SoftDrop>() && sdJet.structure_of<fastjet::contrib::SoftDrop>().delta_R() > 0) {
      thetaG = sdJet.structure_of<fastjet::contrib::SoftDrop>().delta_R() / jetR.value;
      zG = std::min(sdJet.structure_of<fastjet::contrib::SoftDrop>().symmetry(), 0.5 - 1e-6);
      return true;
    }
    return false;
  }

  void fillFake(const JetInfo& hyb, double ptHat, bool fromNoJetEvent)
  {
    if (hyb.weight > 0.) {
      registry.fill(HIST("h2_w_pt_fake"), hyb.pt, std::log10(hyb.weight));
    }
    registry.fill(HIST("h2_ptHat_pt_fake"), hyb.pt, ptHat, hyb.weight);
    registry.fill(HIST("h2_ptHatCnt_pt_fake"), hyb.pt, ptHat);
    registry.fill(HIST("h2_fracPP_pt_fake"), hyb.pt, hyb.ppFraction);
    registry.fill(HIST("hFake"), hyb.pt, hyb.thetaG, hyb.zG, hyb.weight);
    registry.fill(HIST("hFake_cnt"), hyb.pt, hyb.thetaG, hyb.zG);
    if (hyb.ppFraction >= fakePPFractionMin) {
      registry.fill(HIST("hFake_pp"), hyb.pt, hyb.thetaG, hyb.zG, hyb.weight);
      registry.fill(HIST("hFake_pp_cnt"), hyb.pt, hyb.thetaG, hyb.zG);
    } else {
      registry.fill(HIST("hFake_PbPb"), hyb.pt, hyb.thetaG, hyb.zG, 1.);
      registry.fill(HIST("hFake_PbPb_cnt"), hyb.pt, hyb.thetaG, hyb.zG);
    }
    if (fromNoJetEvent) {
      registry.fill(HIST("hFake_noJetColl"), hyb.pt, hyb.thetaG, hyb.zG, hyb.weight);
      registry.fill(HIST("hFake_noJetColl_cnt"), hyb.pt, hyb.thetaG, hyb.zG);
    }
  }

  void fillMiss(const JetInfo& mcp)
  {
    registry.fill(HIST("hMiss"), mcp.pt, mcp.thetaG, mcp.zG, mcp.weight);
    registry.fill(HIST("hMiss_cnt"), mcp.pt, mcp.thetaG, mcp.zG);
  }

  // first geometrically matched jet of jet index i in a matching table, -1 if none
  template <typename M>
  int64_t firstGeoMatch(M const& matchingTable, int64_t i) const
  {
    auto const matching = matchingTable.iteratorAt(i);
    return matching.has_matchedJetGeo() ? matching.matchedJetGeoIds()[0] : -1;
  }

  void processEmbedding(soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos> const& collisions,
                        aod::JetMcCollisions const&,
                        SubJets const& subJets,
                        UnsubJets const& unsubJets,
                        McdJets const& mcdJets,
                        McpJets const& mcpJets,
                        aod::ChargedEventWiseSubtractedJetsMatchedToChargedMCDetectorLevelJets const& subToMcd,
                        aod::ChargedEventWiseSubtractedJetsMatchedToChargedJets const& subToUnsub,
                        aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets const& mcdToMcp,
                        aod::JetTracks const& tracks,
                        aod::JetTracksSub const&,
                        aod::JetParticles const& particles)
  {
    for (auto const& collision : collisions) {
      registry.fill(HIST("h_eventSelection"), 0.5);
      if (!selectPbPbCollision(collision)) {
        continue;
      }
      registry.fill(HIST("h_eventSelection"), 1.5);
      if (!collision.has_mcCollision()) { // no EMB collision in the same BC : nothing embedded
        registry.fill(HIST("h_eventSelection"), 2.5);
        continue;
      }
      auto const mcCollision = collision.mcCollision_as<aod::JetMcCollisions>();
      registry.fill(HIST("h_collision_zdiff"), collision.posZ() - mcCollision.posZ());
      const float weight = useWeights ? mcCollision.weight() : 1.f;
      const double ptHat = useMcPtHard ? static_cast<double>(mcCollision.ptHard()) : ptHatFromWeight(mcCollision.weight());

      auto const subJetsColl = subJets.sliceBy(subJetsPerCollision, collision.globalIndex());
      auto const mcdJetsColl = mcdJets.sliceBy(mcdJetsPerCollision, collision.globalIndex());
      auto const mcpJetsColl = mcpJets.sliceBy(mcpJetsPerMcCollision, mcCollision.globalIndex());
      auto const tracksColl = tracks.sliceBy(tracksPerCollision, collision.globalIndex());
      registry.fill(HIST("h_nParticles"), particles.sliceBy(particlesPerMcCollision, mcCollision.globalIndex()).size());

      // ---------------------------------------------------------------------------
      // outlier cut on the leading MCP / MCD jet of the event
      // ---------------------------------------------------------------------------
      double ptLeadMCP = 0., ptLeadMCD = 0.;
      for (auto const& jet : mcpJetsColl) {
        if (jet.r() == jetRInt) {
          ptLeadMCP = std::max(ptLeadMCP, static_cast<double>(jet.pt()));
        }
      }
      for (auto const& jet : mcdJetsColl) {
        if (jet.r() == jetRInt) {
          ptLeadMCD = std::max(ptLeadMCD, static_cast<double>(jet.pt()));
        }
      }
      registry.fill(HIST("h_ptHatW"), ptHat);
      registry.fill(HIST("h2_ptHatW_vs_leadMCP"), ptLeadMCP, ptHat);
      const bool rejMCP = ptLeadMCP > outlierFactorMCP * ptHat;
      const bool rejMCD = ptLeadMCD > outlierFactorMCD * ptHat;
      if (rejMCP) {
        registry.fill(HIST("h_outlier_mcp_pt"), ptLeadMCP, weight);
      }
      if (rejMCD) {
        registry.fill(HIST("h_outlier_mcd_pt"), ptLeadMCD, weight);
      }
      const bool rejected = rejMCP || rejMCD;
      for (auto const& jet : mcpJetsColl) {
        if (jet.r() != jetRInt) {
          continue;
        }
        registry.fill(HIST("h_cnt_MCP_before_outlier"), jet.pt(), weight);
        if (!rejMCP) {
          registry.fill(HIST("h_cnt_MCP_after_outlier"), jet.pt(), weight);
        }
      }
      for (auto const& jet : mcdJetsColl) {
        if (jet.r() != jetRInt) {
          continue;
        }
        registry.fill(HIST("h_cnt_MCD_before_outlier"), jet.pt(), weight);
        if (!rejMCD) {
          registry.fill(HIST("h_cnt_MCD_after_outlier"), jet.pt(), weight);
        }
      }
      for (auto const& jet : subJetsColl) {
        if (jet.r() != jetRInt) {
          continue;
        }
        registry.fill(HIST("h_cnt_Hybrid_before_outlier"), jet.pt(), weight);
        if (!rejected) {
          registry.fill(HIST("h_cnt_Hybrid_after_outlier"), jet.pt(), weight);
        }
      }

      // ---------------------------------------------------------------------------
      // Pb-Pb only jets (raw pT, non-embedded tracks) : control + overlap veto
      // ---------------------------------------------------------------------------
      std::vector<fastjet::PseudoJet> tracksPb;
      int nTracksPP = 0;
      for (auto const& track : tracksColl) {
        if (!jetderiveddatautilities::selectTrack(track, trackSelection) || track.pt() < trackPtMin || std::abs(track.eta()) > trackEtaMax) {
          continue;
        }
        if (isEmbeddedTrack(track.trackSel())) {
          nTracksPP++;
        } else {
          tracksPb.emplace_back(track.px(), track.py(), track.pz(), std::sqrt(track.p() * track.p() + PionMass * PionMass));
        }
      }
      fastjet::ClusterSequence csPb(tracksPb, jetDefPb);
      auto jetsPb = fastjet::sorted_by_pt(csPb.inclusive_jets());
      registry.fill(HIST("h_pbpb_lead_pt"), jetsPb.size() > 0 ? jetsPb[0].perp() : 0.);
      registry.fill(HIST("h_pbpb_sublead_pt"), jetsPb.size() > 1 ? jetsPb[1].perp() : 0.);
      if (pbOverlapDR > 0. && !rejected) {
        bool overlap = false;
        for (std::size_t ij = 0; ij < jetsPb.size() && ij < 2 && !overlap; ij++) {
          if (jetsPb[ij].perp() <= pbOverlapPtMin) {
            continue;
          }
          for (auto const& jet : mcpJetsColl) {
            if (jet.r() != jetRInt || jet.pt() < truthJetPtMin) {
              continue;
            }
            auto const constituents = jet.tracks_as<aod::JetParticles>();
            if (leadingConstituentPt(constituents) >= leadingTrackPtMax || hasConstituentInPhiVeto(constituents)) {
              continue;
            }
            fastjet::PseudoJet axis;
            axis.reset_PtYPhiM(jet.pt(), jet.eta(), jet.phi());
            if (axis.delta_R(jetsPb[ij]) < pbOverlapDR) {
              overlap = true;
              break;
            }
          }
        }
        if (overlap) {
          registry.fill(HIST("h_eventSelection"), 3.5);
          continue;
        }
      }

      // full rejection of the event if outlier at MCP or MCD level (no cut on the hybrid jet itself)
      if (rejected) {
        registry.fill(HIST("h_eventSelection"), 4.5);
        continue;
      }
      registry.fill(HIST("h_eventSelection"), 5.5);

      // ---------------------------------------------------------------------------
      // jet pools
      //   MCP with eta cut    -> efficiency / response / miss
      //   MCP without eta cut -> purity (1D)
      //   MCD never eta cut
      //   hybrid with eta cut -> purity, response, fake
      // ---------------------------------------------------------------------------
      std::map<int64_t, JetInfo> poolMCP, poolMCPNoEta, poolMCD;
      std::vector<JetInfo> poolHybrid;
      for (auto const& jet : mcpJetsColl) {
        if (jet.r() != jetRInt || jet.pt() < truthJetPtMin) {
          continue;
        }
        auto const constituents = jet.tracks_as<aod::JetParticles>();
        if (leadingConstituentPt(constituents) >= leadingTrackPtMax) {
          continue;
        }
        fillPhiControl(SpeciesMCP, jet, constituents, weight, BeforeVeto);
        if (hasConstituentInPhiVeto(constituents)) {
          continue;
        }
        fillPhiControl(SpeciesMCP, jet, constituents, weight, AfterVeto);
        double thetaG, zG;
        const bool tagged = computeSoftDrop<true>(constituents, thetaG, zG);
        const JetInfo info{jet.globalIndex(), jet.pt(), jet.eta(), jet.phi(), thetaG, zG, weight, 1., -1};
        poolMCPNoEta[jet.globalIndex()] = info;
        if (std::abs(jet.eta()) > jetEtaMax) {
          continue;
        }
        if (tagged) {
          registry.fill(HIST("h_sd_z_mcp"), zG, weight);
        }
        registry.fill(HIST("h_sd_theta_mcp"), thetaG, weight);
        registry.fill(HIST("h_mcp_jets"), jet.pt(), weight);
        poolMCP[jet.globalIndex()] = info;
      }
      for (auto const& jet : mcdJetsColl) {
        if (jet.r() != jetRInt || jet.pt() < truthJetPtMin) {
          continue;
        }
        auto const constituents = jet.tracks_as<aod::JetTracks>();
        if (leadingConstituentPt(constituents) >= leadingTrackPtMax) {
          continue;
        }
        fillPhiControl(SpeciesMCD, jet, constituents, weight, BeforeVeto);
        if (hasConstituentInPhiVeto(constituents)) {
          continue;
        }
        fillPhiControl(SpeciesMCD, jet, constituents, weight, AfterVeto);
        double thetaG, zG;
        computeSoftDrop<false>(constituents, thetaG, zG);
        int64_t mcpId = firstGeoMatch(mcdToMcp, jet.globalIndex());
        if (requireMcdMcpPtMatch && mcpId >= 0) {
          auto const matching = mcdToMcp.iteratorAt(jet.globalIndex());
          auto const ptIds = matching.matchedJetPtIds();
          if (std::find(ptIds.begin(), ptIds.end(), mcpId) == ptIds.end()) {
            mcpId = -1;
          }
        }
        poolMCD[jet.globalIndex()] = JetInfo{jet.globalIndex(), jet.pt(), jet.eta(), jet.phi(), thetaG, zG, weight, 1., mcpId};
      }

      for (auto const& info : poolMCP) {
        const JetInfo& m = info.second;
        registry.fill(HIST("h_jet_mcp_pt"), m.pt, m.weight);
        registry.fill(HIST("hEffDen"), m.pt, m.thetaG, m.zG, m.weight);
        registry.fill(HIST("hEffDen_cnt"), m.pt, m.thetaG, m.zG);
      }
      registry.fill(HIST("h_nTracks"), nTracksPP);
      registry.fill(HIST("h_rho"), collision.rho());
      registry.fill(HIST("h2_rho_ntracks"), tracksPb.size(), collision.rho());
      registry.fill(HIST("h_norm"), 1., weight);
      registry.fill(HIST("h_norm"), 2.);
      for (auto const& track : tracksColl) {
        if (jetderiveddatautilities::selectTrack(track, trackSelection) && track.pt() >= trackPtMin && std::abs(track.eta()) <= trackEtaMax) {
          registry.fill(HIST("h_track_pt"), track.pt(), weight);
        }
      }

      for (auto const& jet : subJetsColl) {
        if (jet.r() != jetRInt) {
          continue;
        }
        auto const constituents = jet.tracks_as<aod::JetTracksSub>();
        // no pT cut on the hybrid jet
        if (!passRecoJetCuts(SpeciesHybrid, jet, constituents, weight)) {
          continue;
        }
        double thetaG, zG;
        const bool tagged = computeSoftDrop<false>(constituents, thetaG, zG);
        if (tagged) {
          registry.fill(HIST("h_sd_z_hybrid"), zG, weight);
        }
        registry.fill(HIST("h_sd_theta_hybrid"), thetaG, weight);
        registry.fill(HIST("h_hybrid_jets"), jet.pt(), weight);
        if (std::abs(jet.eta()) > jetEtaMax) {
          continue;
        }
        // pp fraction from the embedded constituents of the matched unsubtracted jet
        double ppFraction = 0.;
        const int64_t unsubId = firstGeoMatch(subToUnsub, jet.globalIndex());
        if (unsubId >= 0) {
          double ptPP = 0.;
          for (auto const& track : unsubJets.iteratorAt(unsubId).tracks_as<aod::JetTracks>()) {
            if (isEmbeddedTrack(track.trackSel())) {
              ptPP += track.pt();
            }
          }
          ppFraction = ptPP / jet.pt();
        }
        poolHybrid.push_back(JetInfo{jet.globalIndex(), jet.pt(), jet.eta(), jet.phi(), thetaG, zG, weight, ppFraction, firstGeoMatch(subToMcd, jet.globalIndex())});
      }

      // ===========================
      // MATCHING (geometric matching of the O2Physics matching tasks + pp fraction)
      // ===========================
      for (const auto& hyb : poolHybrid) {
        registry.fill(HIST("hPurDen"), hyb.pt, hyb.thetaG, hyb.zG, hyb.weight);
        registry.fill(HIST("hPurDen_cnt"), hyb.pt, hyb.thetaG, hyb.zG);
      }
      if (poolMCD.empty() || poolMCP.empty()) {
        // no matching possible : all hybrid (|eta|<0.7) -> fake, all MCP (eta cut) -> miss
        for (const auto& hyb : poolHybrid) {
          fillFake(hyb, ptHat, true);
        }
        for (auto const& info : poolMCP) {
          fillMiss(info.second);
        }
        continue;
      }

      // PATH 1 : PURITY (1D), MCP without eta cut
      {
        std::set<int64_t> usedMCP;
        for (const auto& hyb : poolHybrid) {
          registry.fill(HIST("h_step0_hybrid_pt"), hyb.pt, hyb.weight);
          if (hyb.matchedId < 0 || !poolMCD.count(hyb.matchedId)) {
            continue;
          }
          registry.fill(HIST("h_step1_hybmcd_pt"), hyb.pt, hyb.weight);
          if (hyb.ppFraction < ppFractionMin) {
            continue;
          }
          registry.fill(HIST("h_step2_frac_pt"), hyb.pt, hyb.weight);
          const int64_t mcpId = poolMCD[hyb.matchedId].matchedId;
          if (mcpId < 0 || !poolMCPNoEta.count(mcpId) || usedMCP.count(mcpId)) {
            continue;
          }
          usedMCP.insert(mcpId);
          registry.fill(HIST("h_step3_mcdmcp_pt"), hyb.pt, hyb.weight);
        }
      }

      // PATH 2 : EFFICIENCY -> response matrices, purity / efficiency base histograms, miss and fake
      {
        std::set<int64_t> usedMCP;
        for (const auto& hyb : poolHybrid) {
          if (hyb.matchedId < 0 || !poolMCD.count(hyb.matchedId) || hyb.ppFraction < ppFractionMin) {
            fillFake(hyb, ptHat, false);
            continue;
          }
          const JetInfo& mcd = poolMCD[hyb.matchedId];
          if (mcd.matchedId < 0 || !poolMCP.count(mcd.matchedId) || usedMCP.count(mcd.matchedId)) {
            fillFake(hyb, ptHat, false);
            continue;
          }
          usedMCP.insert(mcd.matchedId);
          const JetInfo& mcp = poolMCP[mcd.matchedId];

          registry.fill(HIST("h_deltaPt"), hyb.pt - mcd.pt, hyb.weight);
          if (hyb.thetaG > 0. && mcp.thetaG > 0.) {
            registry.fill(HIST("h_response_pt"), hyb.pt, mcp.pt, hyb.weight);
            registry.fill(HIST("h_response_theta"), hyb.thetaG, mcp.thetaG, hyb.weight);
          }
          registry.fill(HIST("hResponse4D"), mcp.pt, hyb.pt, mcp.thetaG, hyb.thetaG, hyb.weight);
          registry.fill(HIST("hResponse6D"), mcp.pt, hyb.pt, mcp.thetaG, hyb.thetaG, mcp.zG, hyb.zG, hyb.weight);
          registry.fill(HIST("hResponse6D_cnt"), mcp.pt, hyb.pt, mcp.thetaG, hyb.thetaG, mcp.zG, hyb.zG);
          registry.fill(HIST("hPurNum"), hyb.pt, hyb.thetaG, hyb.zG, hyb.weight);
          registry.fill(HIST("hPurNum_cnt"), hyb.pt, hyb.thetaG, hyb.zG);
          registry.fill(HIST("hEffNum"), mcp.pt, mcp.thetaG, mcp.zG, mcp.weight);
          registry.fill(HIST("hEffNum_cnt"), mcp.pt, mcp.thetaG, mcp.zG);
          if (hyb.weight > 0.) {
            registry.fill(HIST("h2_w_pt_matched"), hyb.pt, std::log10(hyb.weight));
          }
          registry.fill(HIST("h2_ptHat_pt_matched"), hyb.pt, ptHat, hyb.weight);
          registry.fill(HIST("h2_ptHatCnt_pt_matched"), hyb.pt, ptHat);
          registry.fill(HIST("h2_fracPP_pt_matched"), hyb.pt, hyb.ppFraction);
          registry.fill(HIST("h_jet_matched_pt"), mcp.pt, mcp.weight);
        }
        for (auto const& info : poolMCP) {
          if (!usedMCP.count(info.first)) {
            fillMiss(info.second);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(EmbeddingResponsePurityTask, processEmbedding, "response, purity and efficiency for embedded pp in Pb-Pb", true);

  // Pb-Pb data spectrum to unfold (O2Physics version of fillPbPbLeadCut.sh). Runs on standard Pb-Pb data with
  // estimator-rho, subtractor-eventwiseconstituent and jet-finder-data-charged (processChargedEvtWiseSubJets) :
  // same jets, same cuts and same SoftDrop as the hybrid jets of processEmbedding.
  void processData(soa::Join<aod::JetCollisions, aod::BkgChargedRhos>::iterator const& collision,
                   SubJets const& jets,
                   aod::JetTracks const& tracks,
                   aod::JetTracksSub const&)
  {
    if (!selectPbPbCollision(collision)) {
      return;
    }
    registry.fill(HIST("h_data_nEvents"), 0.5);
    registry.fill(HIST("h_data_rho"), collision.rho());
    registry.fill(HIST("h_data_centrality"), collision.centFT0M());
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection) && track.pt() >= trackPtMin && std::abs(track.eta()) <= trackEtaMax) {
        registry.fill(HIST("h_data_track_pt_eta_phi"), track.pt(), track.eta(), track.phi());
      }
    }

    for (auto const& jet : jets) {
      if (jet.r() != jetRInt) {
        continue;
      }
      auto const constituents = jet.tracks_as<aod::JetTracksSub>();
      if (leadingConstituentPt(constituents) >= leadingTrackPtMax) {
        registry.fill(HIST("h_data_leadcut_pt"), jet.pt());
      }
      // no pT cut on the jet (the reco window is applied in the unfolding)
      if (!passRecoJetCuts(SpeciesData, jet, constituents, 1.)) {
        continue;
      }
      if (std::abs(jet.eta()) > jetEtaMax) {
        continue;
      }
      double thetaG, zG;
      computeSoftDrop<false>(constituents, thetaG, zG);
      registry.fill(HIST("hPbPbPtThetaZg_all"), jet.pt(), thetaG, zG);
      registry.fill(HIST("hPbPbPtTheta_all"), jet.pt(), thetaG);
      if (thetaG > 0.) {
        registry.fill(HIST("hPbPbPtThetaZg"), jet.pt(), thetaG, zG);
        registry.fill(HIST("hPbPbPtTheta"), jet.pt(), thetaG);
      }
    }
  }
  PROCESS_SWITCH(EmbeddingResponsePurityTask, processData, "Pb-Pb data spectrum (pT, theta_g, z_g) to unfold", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<EmbeddingResponsePurityTask>(cfgc, TaskName{"embedding-response-purity"})}; }
