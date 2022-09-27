// Copyright(C) 2022, ATA Engineering, Inc.
// 
// This program is free software; you can redistribute it and /or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, see
// < https://www.gnu.org/licenses/lgpl-3.0.html>.
// 
// Some of the source code for this program was generated with the 
// Pointwise Plugin SDK, and is subject to the Cadence Public License
// Version 1.0.


/****************************************************************************
 *
 * Pointwise Plugin utility functions
 *
 * (C) 2021 Cadence Design Systems, Inc. All rights reserved worldwide.
 *
 ***************************************************************************/

#ifndef _RTCAEPSUPPORTDATA_H_
#define _RTCAEPSUPPORTDATA_H_

/*! \cond */

/*------------------------------------*/
/* CaeUnsLociChem format item setup data */
/*------------------------------------*/
CAEP_BCINFO CaeUnsLociChemBCInfo[] = {
	{ "inflow", 100 },
	{ "supersonicInflow", 101 },
	{ "isentropicInflow", 102 },
	{ "fixedMass", 103},
	{ "fixedMassOutflow", 104 },
	{ "farfield", 105 },
	{ "outflow", 106 },
	{ "supersonicOutflow", 107 },
	{ "symmetry", 108 },
	{ "impermeable", 109 },
	{ "reflecting", 110 },
	{ "viscousWall", 111 },
	{ "wallLaw", 112 },
	{ "periodic", 113 },
	{ "interface", 114 },
	{ "turboInterface", 115 },
	{ "extrapolate", 116 },
	{ "reactive_wall", 117 },
	{ "inflowNRBC", 118 },
	{ "outflowNRBC", 119 },
};
/*------------------------------------*/
CAEP_VCINFO CaeUnsLociChemVCInfo[] = {
    { "volume tag", 200 },
};
/*------------------------------------*/
const char *CaeUnsLociChemFileExt[] = {
    "vars", "vog"
};
/*! \endcond */



/************************************************************************/
/*! \file
\brief Defines Support Data for the CAEP_RTITEM Array

The file \sf{%rtCaepSupportData.h} defines and initilizes the support data for
the global CAEP_RTITEM \ref caepRtItem[] array. The CAE Plugin SDK uses this data
to implement the functions and behaviors required by the \ref DOXGRP_APICAEP.
If you want to see the SDK implementation details, look in the
\sf{/shared/CAEP/apiCAEP.cxx} file.

The SDK file \sf{/shared/CAEP/apiCAEP.cxx} includes \sf{%rtCaepSupportData.h}
prior to the declaration of the \ref caepRtItem[] array as shown below.
\par
\dontinclude apiCAEP.cxx
\skip "rtCaepSupportData.h"
\until ARRAYSIZE(caepRtItem);

When copied from the \sf{src/plugins/templates/CAEP/} folder to your plugins
project folder, \sf{%rtCaepSupportData.h} will contain the support data needed
for the 3 example CAEP_RTITEM array items. This example support data must be
culled and edited as needed for your plugin's implementation.

The support data used by the CAE Plugin SDK includes the following:

\li Array of valid Boundary Conditions (CAEP_BCINFO) definitions.
\li Array of valid Volume Conditions (CAEP_VCINFO) definitions.
\li Array of valid File Extension definitions.

These support arrays are referenced in rtCaepInitItems.h to initialize the
corresponding data members.

\par Example Support Data Usage

The code segments below show how the example support data is implemented by
the SDK in the \sf{%rtCaepSupportData.h} and rtCaepInitItems.h template files.

The example support data declaration and initialization in
\sf{%rtCaepSupportData.h}:
\par
\dontinclude rtCaepSupportData.h
\skip CaeUnsLociChemBCInfo[]
\until "xml"
\skipline };

The support data referenced as initilizers in rtCaepInitItems.h:
\par
\dontinclude rtCaepInitItems.h
\skip == CAEP_BCINFO*
\until ARRAYSIZE(CaeUnsLociChemFileExt)
*/

#endif /* _RTCAEPSUPPORTDATA_H_ */
