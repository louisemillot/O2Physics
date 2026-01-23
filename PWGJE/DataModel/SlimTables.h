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

///
/// \file   SlimTables.h
/// \author Millot Louise <louise.millot@cern.ch>
/// \since  2024-11-27
/// \brief  Header for the SlimTables task for the analysis of the reduced tables.
///

#ifndef PWGJE_DATAMODEL_SLIMTABLES_H_
#define PWGJE_DATAMODEL_SLIMTABLES_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <Rtypes.h>

namespace o2::aod
{
DECLARE_SOA_TABLE(SlimCollisions, "AOD", "SlimCollisions", o2::soa::Index<>,
                  o2::aod::collision::PosZ);
using SlimCollision = SlimCollisions::iterator;
// {
// DECLARE_SOA_INDEX_COLUMN(Collision, collision);
// DECLARE_SOA_COLUMN(PosZ, posZ, float);
// // DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);
// // DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);
// // DECLARE_SOA_COLUMN(Weight, weight, float);
// // DECLARE_SOA_COLUMN(EventSel, eventSel, uint16_t);
// // DECLARE_SOA_COLUMN(TrackOccupancyInTimeRange, trackOccupancyInTimeRange, int);
// } // namespace slimcollision

namespace slimtracks
{
DECLARE_SOA_INDEX_COLUMN(SlimCollision, slimcollision);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(E, energy, float);
// DECLARE_SOA_DYNAMIC_COLUMN(Px, px,
//                            [](float pt, float phi) -> float { return pt * std::cos(phi); });
// DECLARE_SOA_DYNAMIC_COLUMN(Py, py,
//                            [](float pt, float phi) -> float { return pt * std::sin(phi); });
// DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz,
//                            [](float pt, float eta) -> float { return pt * std::sinh(eta); });
// DECLARE_SOA_DYNAMIC_COLUMN(Energy, energy,
//                            [](float pt, float eta) -> float { return std::sqrt((pt * std::cosh(eta) * pt * std::cosh(eta)) + (jetderiveddatautilities::mPion * jetderiveddatautilities::mPion)); });
} // namespace slimtracks
DECLARE_SOA_TABLE(SlimTracks, "AOD", "SlimTracks",
                  o2::soa::Index<>,
                  slimtracks::SlimCollisionId,
                  slimtracks::Pt,
                  slimtracks::Eta,
                  slimtracks::Phi,
                  slimtracks::Px,
                  slimtracks::Py,
                  slimtracks::Pz,
                  slimtracks::E);
using SlimTrack = SlimTracks::iterator;
} // namespace o2::aod

#endif // PWGJE_DATAMODEL_SLIMTABLES_H_
