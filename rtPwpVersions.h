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


#ifndef _RTPWPVERSIONS_H_
#define _RTPWPVERSIONS_H_


/*------------------------------------------------
   PWP api version with which plugin conforms
------------------------------------------------*/
/*! \brief The PWP-API major version value.
*/
#define VERSION_PWP_MAJOR   1

/*! \brief The PWP-API minor version value.
*/
#define VERSION_PWP_MINOR   0

/*! \brief This macro is used for static initialization of a PWP_VERSION
struct to the current \sf{VERSION_PWP_MAJOR} and \sf{VERSION_PWP_MINOR} values.
*/
#define VERSION_PWP_INIT    {VERSION_PWP_MAJOR, VERSION_PWP_MINOR}


/*------------------------------------------------
        plugin software release version
------------------------------------------------*/
/*! \brief The software release major version value.
*/
#define VERSION_LIB_MAJOR   1

/*! \brief The software release minor version value.
*/
#define VERSION_LIB_MINOR   0

/*! \brief This macro is used for static initialization of a PWP_VERSION
struct to the current \sf{VERSION_LIB_MAJOR} and \sf{VERSION_LIB_MINOR} values.
*/
#define VERSION_LIB_INIT    {VERSION_LIB_MAJOR, VERSION_LIB_MINOR}


/************************************************************************/
/*! \file
\brief Defines Implementation Version Information

This file contains 2 version value macro pairs. These values are used
throughout the Plugin SDK. One pair (VERSION_PWP_ prefix) defines the PWP-API
version. The other pair (VERSION_LIB_ prefix) defines the plugin binary release
version.

The plugin developer should edit the VERSION_PWP_MAJOR, VERSION_PWP_MINOR,
VERSION_LIB_MAJOR, and VERSION_LIB_MINOR values as needed.

*/

#endif /* _RTPWPVERSIONS_H_ */
