// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details.

#ifndef PWGJE_DATAMODEL_SLIMTABLES_H_
#define PWGJE_DATAMODEL_SLIMTABLES_H_
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cstdint>


namespace o2::aod
{
namespace slimtracks
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);                    //! Id of collision
DECLARE_SOA_COLUMN(Pt, pt, float);

}

namespace slimparticles
{
// Track info
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);

}


DECLARE_SOA_TABLE(SlimTracks, "AOD", "SlimTracks",
                  o2::soa::Index<>,
                  slimtracks::Pt 
                  );

DECLARE_SOA_TABLE(SlimParticles, "AOD", "SlimParticles",
                  slimparticles::Px,
                  slimparticles::Py,
                  slimparticles::Pz
                  );

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_SLIMTABLES_H_
