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

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SlimTablesProducer {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<float> minPt{"minPt", 2.0, "min pT to save"};
  Configurable<float> maxDCA{"maxDCA", 0.1, "max DCA"};
  Configurable<float> etaWindow{"etaWindow", 0.8, "eta window"};
  Configurable<bool> skipUninterestingEvents{"skipUninterestingEvents", true, "skip collisions without particle of interest"};
  Configurable<int> minTPCNClsCrossedRows{"minTPCNClsCrossedRows", 80, "min TPC crossed rows"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "Event selection"};

  std::vector<int> eventSelectionBits;

  void init(InitContext&)
  {
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisPt{nBinsPt, 0, 10, "p_{T}"};
    histos.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
    histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt});
    histos.add("hTracksPerCollision", "Tracks per collision;collisionId;N tracks", {HistType::kTH1F, {{500, 0.0, 500.0}}});
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
  }

  Produces<o2::aod::SlimCollisions> slimCollisions;
  Produces<o2::aod::SlimTracks> slimTracks;

  // Look at primary tracks only
  Filter trackFilter = nabs(aod::track::dcaXY) < maxDCA && nabs(aod::track::eta) < etaWindow && aod::track::pt > minPt;
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut);

  using myCompleteTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;
  using myFilteredTracks = soa::Filtered<myCompleteTracks>;

  void process(soa::Filtered<aod::JetCollisions>::iterator const& collision, myFilteredTracks const& tracks)
  {
    histos.fill(HIST("h_collisions"), 0.5); // Compte tous les événements qui entrent dans la fonction, avant toute sélection
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, false)) {
      return;
    }
    histos.fill(HIST("h_collisions"), 1.5);

    if (tracks.size() < 1 && skipUninterestingEvents) // si l'event n'a aucune track ET j'ai demandé de skipper les événements inintéressants, on sort immédiatement.
      return;
    bool interestingEvent = false; // on suppose que l'événement n'est pas intéressant au depart
    for (const auto& track : tracks) {
      if (track.tpcNClsCrossedRows() < minTPCNClsCrossedRows) // On rejette les tracks avec pas assez de clusters TPC
        continue;                                             // on passe à la track suivante
      interestingEvent = true;                                // si une track a un NClsCrossedRows de 3 alors que j'ai demande 5 minimum, on l'ignore et on passe à la suivante et si la piste suivante est bonne alors interestingEvent devient true
    }
    if (!interestingEvent && skipUninterestingEvents) // si aucune track est de bonne qualité mais que skipUninterestingEvents est true alors on jet l'événement
      return;
    histos.fill(HIST("h_collisions"), 2.5);
    slimCollisions(collision.posZ());
    auto slimCollId = slimCollisions.lastIndex();
    int nTracksThisCollision = 0;
    for (const auto& track : tracks) {
      if (track.tpcNClsCrossedRows() < minTPCNClsCrossedRows)
        continue; // remove badly tracked
      nTracksThisCollision++;
      histos.get<TH1>(HIST("ptHistogram"))->Fill(track.pt());
      float mass = jetderiveddatautilities::mPion;
      float p = track.pt() * std::cosh(track.eta());
      float energy = std::sqrt(p * p + mass * mass);
      slimTracks(track.collisionId(), track.pt(), track.eta(), track.phi(), track.px(), track.py(), track.pz(), energy); // all that I need for posterior analysis!
      LOG(info) << "collision.globalIndex() = " << collision.globalIndex() << " track.collisionId() = " << track.collisionId() << " slimCollId = " << slimCollId;
      histos.get<TH1>(HIST("hTracksPerCollision"))->Fill(nTracksThisCollision);
    }
  }
  PROCESS_SWITCH(SlimTablesProducer, process, "Produce slim collision table", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SlimTablesProducer>(cfgc)};
}
