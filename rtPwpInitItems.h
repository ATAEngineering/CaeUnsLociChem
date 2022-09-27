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

#ifndef _RTPWPINITITEMS_H_
#define _RTPWPINITITEMS_H_

/*............................*/
{
    {PWP_API_PLUGIN "/1.0", VERSION_PWP_INIT},
    0,
},
/*............................*/
{
    {PWP_API_CAE_EXPORT "/1.0", {1,0}},
    0,
},

/************************************************************************/
/*! \file
\brief Static Initialization Data for the PWU_RTITEM Array

The file \sf{%rtPwpInitItems.h} defines the static, compile-time initialization
of the global PWU_RTITEM \ref pwpRtItem[] array. The PWP Plugin SDK uses this
array to implement the functions and behaviors required by the
\ref DOXGRP_APIPWP. If you want to see the SDK implementation details, look in
the \sf{/shared/PWP/apiPWP.cxx} file.

The SDK file \sf{/shared/PWP/apiPWP.cxx} includes \sf{%rtPwpInitItems.h} as shown
below.
\par
\dontinclude apiPWP.cxx
\skip pwpRtItem[] =
\until };

A valid plugin will have a minimum of 2 PWU_RTITEM array items. One entry
identifies the supported \ref DOXGRP_APIPWP "PWP-API" version. The other
entry identifies an additional supported plugin API specification. See
\ref DOXGRP_APIPWP_BASENAMES for a list of the supported APIs.

The format of \sf{%rtPwpInitItems.h} must be valid for the static
initialization of an array of C-struct's. It is important to note that some of
\ref PWU_RTITEM's data members are also structs. This will require curly-braces
\p {} around these nested data members. If you are not familiar with static
initialization, see the \ref example_cstruct_init page.

When copied from the \sf{src/plugins/templates/PWP/} folder to your plugins
project folder, \sf{%rtPwpInitItems.h} will contain example initilization data
for 2 PWU_RTITEM array items. This example data must be culled and edited to
define the settings appropriate for your plugin's implementation.

The example data from the SDK \sf{%rtPwpInitItems.h} template file is shown
below. It specifies that the \em "Plugin-PWP/1.0" API specification is supported
along with the \em "Export-CAE/1.0" API specification.
\par
\dontinclude rtPwpInitItems.h
\skip .......
\until CAEP_API_EXPORT
\skipline },
*/

#endif /* _RTPWPINITITEMS_H_ */
