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

#ifndef _RTCAEPINSTANCEDATA_H_
#define _RTCAEPINSTANCEDATA_H_

/*----------------------------------------------------------------------
   You can uncomment and customize the section below to add your own
   CAEP_RTITEM data as a single struct.

    Remember to properly initialize any data you add in the
    rtCaepInitItems.h file.
----------------------------------------------------------------------*/
/*
struct CAEUNSLOCICHEM_DATA {
    int data1;
    int data2;
    float data3;
    char *pStr;
    // I assume you get the idea...
};

// This macro is included in the CAEP_RTITEM structure declaration:
#define CAEP_RUNTIME_INSTDATADECL   CAEUNSLOCICHEM_DATA myData;
*/




/************************************************************************/
/*! \file
\brief Customizes the typedef CAEP_RTITEM declaration

If your plugin needs custom, per-instance CAE data members added to the
CAEP_RTITEM typedef, you will need to modify the \sf{%rtCaepInstanceData.h} file.

The SDK file apiCAEPUtils.h includes \sf{%rtCaepInstanceData.h} before the
CAEP_RTITEM struct is typedef'ed. The \ref CAEP_RUNTIME_INSTDATADECL macro is
used at the end of the struct as shown below. By default, the
\ref CAEP_RUNTIME_INSTDATADECL macro resolves to nothing (no data) and does not
change the CAEP_RTITEM typedef.
\par
\dontinclude apiCAEPUtils.h
\skip CAEP_RTITEM_t {
\until FormatInfo;
...snip...
\skip opAborted;
\until CAEP_RTITEM;

\par Declaring the CAEP_RUNTIME_INSTDATADECL Macro

There are 2 similar approaches to adding data members to the CAEP_RTITEM
typedef.

\li Add a single struct containing all your custom data members.
\li Add each custom data member as a separate value.

The examples below show the same 4 data members being added to CAEP_RTITEM
using the two approaches. Which approach you choose is a matter of personal
preference.

\par Using a Single Struct to Extend CAEP_RTITEM

In \sf{%rtCaepInstanceData.h}, define your custom data struct and then define
the \ref CAEP_RUNTIME_INSTDATADECL Macro as if you were declaring an instance
of this struct. The code below shows this approach.
\par
\code
struct CAEUNSLOCICHEM_DATA {
    int data1;
    int data2;
    float data3;
    char *pStr;
};

// This macro is included in the CAEP_RTITEM structure declaration:
#define CAEP_RUNTIME_INSTDATADECL   CAEUNSLOCICHEM_DATA myData;
\endcode

Your single, struct data member is now part of CAEP_RTITEM. You can access
its members using the pRti pointer passed into the runtimeWrite() function as
shown in the code below.
\par
\code
PWP_BOOL runtimeWrite(CAEP_RTITEM *pRti, ...snip...)
{
    pRti->myData.data1;
    pRti->myData.data2;
    pRti->myData.data3;
    pRti->myData.pStr;
}
\endcode

\par Using a Separate Values to Extend CAEP_RTITEM

In \sf{%rtCaepInstanceData.h}, define the \ref CAEP_RUNTIME_INSTDATADECL Macro as if you
were declaring your data members as local variables. The code below shows this
approach.
\par
\code
// This macro is included in the CAEP_RTITEM structure declaration:
#define CAEP_RUNTIME_INSTDATADECL   int data1; \
                                    int data2; \
                                    float data3; \
                                    char *pStr;
\endcode

Your separate data members are now part of CAEP_RTITEM. You can access them
using the pRti pointer passed into the runtimeWrite() function as shown in the
code below.
\par
\code
PWP_BOOL runtimeWrite(CAEP_RTITEM *pRti, ...snip...)
{
    pRti->data1;
    pRti->data2;
    pRti->data3;
    pRti->pStr;
}
\endcode

\note
To prevent compiler warnings or errors, be sure to add the appropriate static
initializers to rtCaepInitItems.h for your custom data members!
*/
#endif /* _RTCAEPINSTANCEDATA_H_ */
