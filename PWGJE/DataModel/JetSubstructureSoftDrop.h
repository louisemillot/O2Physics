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

/// \file JetSubstructureSoftDrop.h
/// \brief Extension of jet splitting tables for SoftDrop analysis
/// \author Louise Millot <louise.millot@cern.ch>

#ifndef PWGJE_DATAMODEL_JETSUBSTRUCTURESOFTDROP_H_
#define PWGJE_DATAMODEL_JETSUBSTRUCTURESOFTDROP_H_

#include "PWGJE/DataModel/JetSubstructure.h"

namespace o2::aod
{

// ----------------------------------------------------------------------
// SoftDrop extension: MC Detector Level EventWiseSubtracted splitting
// ----------------------------------------------------------------------

JETSPLITTING_TABLE_DEF(ChargedMCDetectorLevelEventWiseSubtracted,
                       "CDEWS",
                       chargedmcdetectorleveleventwisesubtracted,
                       JTracks,
                       CMCDJetCOs)

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETSUBSTRUCTURESOFTDROP_H_
