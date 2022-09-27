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

#ifndef _RTCAEPINITITEMS_H_
#define _RTCAEPINITITEMS_H_

/*! \cond */
/*.................................................
    initialize caepRtItem[0]
*/
#define ID_CaeUnsLociChem  20
{
    //== CAEP_FORMATINFO FormatInfo
    {   PWP_SITE_GROUPNAME,     // const char *group
        "Loci-Chem",             // const char *name
        MAKEGUID(ID_CaeUnsLociChem),  // PWP_UINT32 id

        PWP_FILEDEST_BASENAME,  // PWP_ENUM_FILEDEST fileDest

        PWP_TRUE,               // PWP_BOOL allowedExportConditionsOnly
        PWP_FALSE,               // PWP_BOOL allowedVolumeConditions

        PWP_FALSE,               // PWP_BOOL allowedFileFormatASCII
        PWP_TRUE,              // PWP_BOOL allowedFileFormatBinary
        PWP_FALSE,              // PWP_BOOL allowedFileFormatUnformatted

        PWP_FALSE,               // PWP_BOOL allowedDataPrecisionSingle
        PWP_TRUE,               // PWP_BOOL allowedDataPrecisionDouble

        PWP_FALSE,               // PWP_BOOL allowedDimension2D
        PWP_TRUE                // PWP_BOOL allowedDimension3D
    },

    &pwpRtItem[1],  // PWU_RTITEM*

    CaeUnsLociChemBCInfo,            // CAEP_BCINFO* pBCInfo -- format BCs or NULL
    ARRAYSIZE(CaeUnsLociChemBCInfo), // PWP_UINT32 BCCnt -- # format BCs

    CaeUnsLociChemVCInfo,            // CAEP_VCINFO* pVCInfo -- format VCs or NULL
    ARRAYSIZE(CaeUnsLociChemVCInfo), // PWP_UINT32 VCCnt     -- # format VCs

    CaeUnsLociChemFileExt,            // const char **pFileExt -- file extensions
    ARRAYSIZE(CaeUnsLociChemFileExt), // PWP_UINT32 ExtCnt     -- # file extensions

    //== PWP_BOOL  elemType[PWGM_ELEMTYPE_SIZE]; -- un/supported elem
    {   PWP_TRUE,              // elemType[PWGM_ELEMTYPE_BAR]
        PWP_TRUE,              // elemType[PWGM_ELEMTYPE_HEX]
        PWP_TRUE,              // elemType[PWGM_ELEMTYPE_QUAD]
        PWP_TRUE,              // elemType[PWGM_ELEMTYPE_TRI]
        PWP_TRUE,              // elemType[PWGM_ELEMTYPE_TET]
        PWP_TRUE,              // elemType[PWGM_ELEMTYPE_WEDGE]
        PWP_TRUE,              // elemType[PWGM_ELEMTYPE_PYRAMID]
        PWP_TRUE },            // elemType[PWGM_ELEMTYPE_POINT]

    0,  // FILE *fp

    /*
    // PWU_UNFDATA UnfData
    {   0,                // PWP_UINT32 status
        0,                // FILE *fp
        0,                // sysFILEPOS fPos
        PWP_FALSE,        // PWP_BOOL inRec
        0,                // PWP_UINT32 recBytes
        0,                // PWP_UINT32 totRecBytes
        0,                // PWP_UINT32 recCnt
        PWP_ENDIAN_ERROR, // PWP_ENDIANNESS endianness
        0 },              // PWP_UINT32     fixedBytes
    */

    // this PWU_UNFDATA is from the old SDK, the new version above doesn't compile 
    // without error
    { 0,          // PWP_UINT32 status 
      0,          // FILE *fp 
      0,          // sysFILEPOS fPos 
      PWP_FALSE,  // PWP_BOOL hadError 
      PWP_FALSE,  // PWP_BOOL inRec 
      0,          // PWP_UINT32 recBytes 
      0,          // PWP_UINT32 totRecBytes 
      0 },        // PWP_UINT32 recCnt 

    0,  // PWGM_HGRIDMODEL model

    0,  // const CAEP_WRITEINFO *pWriteInfo

    0,    // PWP_UINT32 progTotal
    0,    // PWP_UINT32 progComplete
    {0},  // clock_t clocks[CAEPU_CLKS_SIZE];
    0,    // PWP_BOOL opAborted

    // if you added any custom data in rtCaepInstanceData.h, you need to
    // initialize it here. The init below matches the example CAEUNSLOCICHEM_DATA
    // struct given in rtCaepInstanceData.h
    /*
    {   0,
        0,
        0.0,
        "string" },
    */
},
/*! \endcond */


/************************************************************************/
/*! \file
\brief Static Initialization Data for the CAEP_RTITEM Array

The file \sf{%rtCaepInitItems.h} defines the static, compile-time initialization
of the global CAEP_RTITEM \ref caepRtItem[] array. The CAE Plugin SDK uses this
array to implement the functions and behaviors required by the
\ref DOXGRP_APICAEP. If you want to see the SDK implementation details,
look in the \sf{/shared/CAEP/apiCAEP.cxx} file.

The SDK file \sf{/shared/CAEP/apiCAEP.cxx} includes \sf{%rtCaepInitItems.h} as 
shown below.
\par
\dontinclude apiCAEP.cxx
\skip caepRtItem[] =
\until };

The format of \sf{%rtCaepInitItems.h} must be valid for the static
initialization of an array of C-struct's. It is important to note that some of
\ref CAEP_RTITEM's data members are also structs. This will require
curly-braces \p {} around these nested data members. If you are not familiar
with static initialization, see the \ref example_cstruct_init page.

\note
If you add custom data members to CAEP_RTITEM using rtCaepInstanceData.h, be
sure to add the additional static initializers when editing 
\sf{%rtCaepInitItems.h} to prevent compiler warnings or errors!
*/

#endif /* _RTCAEPINITITEMS_H_ */
