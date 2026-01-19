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

/// \file   trackEfficiency.cxx
/// \author Louise MILLOT <millot.louise@cern.ch>
/// \brief standalone procedure for embedding : this task creates slim tables for tracks, collisions, mcparticles and mccollisions.

#include "PWGJE/DataModel/SlimTables.h"

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

struct CountSlimTracks {

  // Pour pouvoir faire sliceBy(tracksPerCollision, collisionId)
  Preslice<aod::SlimTracks> tracksPerCollision = aod::slimtracks::collisionId;

  void processTracks(aod::SlimCollisions const& collisions,
                     aod::SlimTracks const& tracks)
  {
    for (auto const& coll : collisions) {

      // Filtre les tracks qui appartiennent à la collision coll
      auto tracksInColl = tracks.sliceBy(tracksPerCollision, coll.globalIndex());

      // Nombre de tracks dans cette collision
      int nTracks = tracksInColl.size();

      // Affiche le résultat
      LOG(info) << "Collision " << coll.globalIndex() << " has " << nTracks << " tracks";
    }
  }
  PROCESS_SWITCH(CountSlimTracks, processTracks, "number of tracks per collisions ", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CountSlimTracks>(cfgc)};
}
