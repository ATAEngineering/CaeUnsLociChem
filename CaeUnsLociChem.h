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
 * class CaeUnsLociChem
 *
 * (C) 2021 Cadence Design Systems, Inc. All rights reserved worldwide.
 *
 ***************************************************************************/

#ifndef _CAEUNSLOCICHEM_H_
#define _CAEUNSLOCICHEM_H_

#include "apiGridModel.h"
#include "apiPWP.h"

#include "CaePlugin.h"
#include "CaeUnsGridModel.h"

#include <string>
#include <vector>
#include <array>
#include <utility>
#include <set>
#include <map>
#include <hdf5.h>

using std::string;
using std::vector;
using std::array;
using std::pair;
using std::set;
using std::map;

//***************************************************************************
//***************************************************************************
//***************************************************************************

// A struct to hold the position of nodes in the mesh.
// Having the data contiguous makes it easy to pass to HDF5.
struct node {
    double x;
    double y;
    double z;
};


class CaeUnsLociChem :
    public CaeUnsPlugin,
    public CaeFaceStreamHandler,
    public CaeUnsCustomIndexHandler {

public:
    CaeUnsLociChem(CAEP_RTITEM *pRti, PWGM_HGRIDMODEL model,
        const CAEP_WRITEINFO *pWriteInfo);

    ~CaeUnsLociChem() override;

    static bool create(CAEP_RTITEM &rti);
    static void destroy(CAEP_RTITEM &rti);

private:

    bool        beginExport() override;
    PWP_BOOL    write() override;
    bool        endExport() override;

	// added member functions
	bool writeVars() const;
	bool writeVog();
	bool writeVogSurfs(const hid_t&) const;
	bool writeVogNodes(const hid_t&) const;
	template<class T>
	void writeVectorToHDF5(const hid_t&, const vector<T>&, const string&, const hid_t&) const;
	vector<node> getAllNodes() const;
	bool writeVogFaces(const hid_t&);
	vector<unsigned char> encodeFaceCluster(set<int>&, set<int>&, set<int>&) const;
	set<int> faceCluster(set<int>&, vector<unsigned char>&, vector<unsigned short>&) const;
	int getFaceSize(const int&) const;
	PWGM_ENUM_ELEMTYPE getFaceType(const int&) const;
	bool writeBC(const CaeUnsPatch&, const bool&, PwpFile&) const;
	double getUnitsFactor() const;
	string getInputEnum(const string&, const string&) const;
	string getInputVal(const string&, const string&, bool = false) const;
	bool getInputBool(const string&) const;
	bool getTurbulenceModule(string&) const;
	bool usingNRBCs(string&) const;
	void swapOrientation(const int&, const PWGM_ENUM_ELEMTYPE&);
	void colorMatrix();
	bool isInterior(const int&) const;
	bool isInterior(const pair<int, int>&) const;

    // face streaming handlers
    PWP_UINT32 streamBegin(const PWGM_BEGINSTREAM_DATA &data) override;
    PWP_UINT32 streamFace(const PWGM_FACESTREAM_DATA &data) override;
    PWP_UINT32 streamEnd(const PWGM_ENDSTREAM_DATA &data) override;

    // custom indexing mapper
    PWP_UINT64 mapIndex(const PWP_UINT64 pwgmNdx) override;

private:
	PWGM_ENUM_FACEORDER faceOrder_;          // order in which faces will be streamed
	PWP_UINT32 totalNumFaces_;               // total number of faces in mesh
	PWP_UINT32 numInteriorFaces_;            // total number of interior faces in mesh
	vector<pair<int, int>> face2cell_;       // face to cell maps (cl, cr)
	vector<array<int, 4>> face2node_;        // face to node maps

};

// utility functions to write out data to HDF5
void writeUnsignedVal(vector<unsigned char>&, unsigned long long);
void writeSignedVal(vector<unsigned char>&, long long);
void writeTable(vector<unsigned char>&, vector<int>&);

string center(const string&, const int&);

void inverseMap(map<int, vector<int>>&, vector<pair<int, int>>&);

#endif // _CAEUNSLOCICHEM_H_
