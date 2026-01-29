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

#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/SlimTables.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include <CommonConstants/MathConstants.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>

#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SlimTablesProducer {

  Configurable<float> minPt{"minPt", 0.15, "min pT to save"};
  Configurable<float> maxPt{"maxPt", 200.0, "max pT to save"};
  Configurable<float> minEta{"minEta", -0.9, "min eta to save"};
  Configurable<float> maxEta{"maxEta", 0.9, "max eta to save"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};

  Produces<o2::aod::SlimCollisions> slimCollisions;
  Produces<o2::aod::SlimTracks> slimTracks;

  Filter trackFilter = (aod::jtrack::pt >= minPt && aod::jtrack::pt < maxPt && aod::jtrack::eta > minEta && aod::jtrack::eta < maxEta);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut);

  void process(soa::Filtered<o2::aod::JetCollisions>::iterator const& collision,
               soa::Filtered<soa::Join<aod::JetTracks, aod::JTrackExtras, aod::JTrackPIs>> const& tracks)
  {

    int nTracksThisCollision = 0;
    int collisionId = collision.globalIndex();
    slimCollisions(collision.posZ());
    for (const auto& track : tracks) {
      nTracksThisCollision++;
      slimTracks(track.collisionId(), track.pt(), track.eta(), track.phi(), track.px(), track.py(), track.pz());
    }
    LOG(info) << "Number of tracks saved for collision " << collisionId << " : " << nTracksThisCollision;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SlimTablesProducer>(cfgc)};
}
