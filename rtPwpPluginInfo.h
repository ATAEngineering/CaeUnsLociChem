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

#ifndef _RTPWPPLUGININFO_H_
#define _RTPWPPLUGININFO_H_

/*! \cond */

    /* initialize the PWP_PLUGININFO data returned by
       PwpGetPluginInfo(PWP_PLUGININFO *pInfo)
    */
    VERSION_PWP_INIT,         // conforms to this PWP-API version
    VERSION_LIB_INIT,         // software library release version
    "ATA Engineering, Inc.",        // company/author description
    "https://www.ata-e.com",  // support description (phone, web-link).
    "(C) 2016-2022 ATA Engineering, Inc.", // copyright description
    0,                        // number of APIs (auto-set at runtime)
    0,                        // default msg callback (auto-set at runtime)
    0,                        // spy msg callback (auto-set at runtime)

/*! \endcond */

/************************************************************************/
/*! \file
\brief Static Initialization Data for the PWP_PLUGININFO structure

initialize the PWP_PLUGININFO data 

The file \sf{%rtPwpPluginInfo.h} defines the static, compile-time
initialization of the PWP_PLUGININFO data struct returned by the PWP-API
function PwpGetPluginInfo(). If you want to see the implementation details,
look in the \sf{/shared/PWP/apiPWP.cxx} file.

The SDK file \sf{/shared/PWP/apiPWP.cxx} includes \sf{%rtPwpPluginInfo.h} as
shown below.
\par
\dontinclude apiPWP.cxx
\skip PWP_PLUGININFO info =
\until };

The format of \sf{%rtPwpPluginInfo.h} must be valid for the static
initialization of a C-struct. If you are not familiar with static
initialization, see the \ref example_cstruct_init page.

When copied from the \sf{src/plugins/templates/PWP/} folder to your plugins
project folder, \sf{%rtPwpPluginInfo.h} will contain example initilization
data. This example data must be edited to define the values appropriate for
your plugin's implementation.

\par Example PWP_PLUGININFO data

The example data from the SDK \sf{%rtPwpPluginInfo.h} template file is shown
below.
\par
\dontinclude rtPwpPluginInfo.h
\skip VERSION_PWP_INIT
\until spy msg callback

Please notice that the last three data members are initilized to 0. These
values are set automatically at runtime by PwpGetPluginInfo() as shown in the
code below.
\par
\dontinclude apiPWP.cxx
\skip info.apiCount
\until info.spyCB
*/

#endif /* _RTPWPPLUGININFO_H_ */
