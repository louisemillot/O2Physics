// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details.

#ifndef PWGJE_DATAMODEL_SLIMTABLES_H_
#define PWGJE_DATAMODEL_SLIMTABLES_H_

#include <Framework/ASoA.h>
#include <vector>

namespace o2::aod
{

// ---------------- SlimTracks ----------------
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);

DECLARE_SOA_TABLE(SlimTracks, "AOD", "SlimTracks",
                  o2::soa::Index<>,
                  Px, Py, Pz);

// ---------------- SlimParticles ----------------
DECLARE_SOA_COLUMN(Energy, energy, float);

DECLARE_SOA_TABLE(SlimParticles, "AOD", "SlimParticles",
                  o2::soa::Index<>,
                  Px, Py, Pz, Energy);

// ---------------- SlimCollision ----------------
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(Centrality, centrality, float);
DECLARE_SOA_COLUMN(EventSel, eventSel, uint8_t);
DECLARE_SOA_COLUMN(EventWeight, eventWeight, float);

DECLARE_SOA_TABLE(SlimCollision, "AOD", "SlimCollision",
                  PosZ, Centrality, EventSel, EventWeight);

// ---------------- SlimMcCollision ----------------
DECLARE_SOA_TABLE(SlimMcCollision, "AOD", "SlimMcCollision",
                  PosZ, Centrality, EventSel, EventWeight);

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_SLIMTABLES_H_