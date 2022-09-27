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

#include "apiCAEP.h"
#include "apiCAEPUtils.h"
#include "apiGridModel.h"
#include "apiPWP.h"
#include "runtimeWrite.h"
#include "pwpPlatform.h"

#include "CaePlugin.h"
#include "CaeUnsGridModel.h"
#include "CaeUnsLociChem.h"

#include <string>
#include <vector>
#include <array>
#include <utility>
#include <set>
#include <map>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include <hdf5.h>

using std::string;
using std::ostringstream;
using std::vector;
using std::set;
using std::map;
using std::multimap;
using std::array;

/* This plugin writes out the grid in VOG file format for Loci/Chem.

It follows the process outlined in the FVMGridWriter.cc file packaged with Loci.
All of Loci/Chem's grid conversion utilities use this procedure. However,
before calling this procedure they call VOG::optimizeMesh to arrange the face and
node maps in a favorable order.

This exporter also mimics the behavior of VOG::colorMatrix to ensure that the face
maps are oriented in the way that the matrix solvers in Loci/Chem expect them to be.
The VOG::colorMatrix utilities are in vogtools.cc and distribute_container.cc

The VOG file format stores the mesh data in a face-centered way. To write a VOG
file it is required to have the coorinates of all the nodes, a map from each
face to the cell to the left, a map from each face to the cell to the right, a map
from each face to all the nodes making up that face, and the boundary surface IDs
and names. For this reason the data is accessed from Pointwise in the face-centered
view provided by the plugin API. Pointwise does not store face centered data, so the
face data is streamed one at a time. This plugin accesses each face as it is streamed
and constructs the necessary maps and stores them for writing out to HDF5 later. The
face streaming is only done once. STL vectors are used instead of maps because the
faces are accessed in order. For the face-to-cell maps, when a ghost cell is
encountered, the cell ID assigned is that of the boundary surface to which the face
is a part of, multiplied by -1. This is consistent with how the cl and cr maps are
constructed within Loci.

*/

//***************************************************************************
//***************************************************************************
//***************************************************************************

CaeUnsLociChem::CaeUnsLociChem(CAEP_RTITEM *pRti, PWGM_HGRIDMODEL
        model, const CAEP_WRITEINFO *pWriteInfo) :
    CaeUnsPlugin(pRti, model, pWriteInfo),
    faceOrder_(PWGM_FACEORDER_BOUNDARYFIRST),
    totalNumFaces_(0), numInteriorFaces_(0), face2cell_(0), face2node_(0)
{
}

CaeUnsLociChem::~CaeUnsLociChem()
{
}

bool
CaeUnsLociChem::beginExport()
{
    // bool ret = true; // all ok
    //PWP_BOOL doDump;
    // Attribute values set by:
    //   GUI using "CAE, Set Solver Attributes..."
    //   Glyph using "pw::Application setCAESolverAttribute"
    //model_.getAttribute("debugDump", doDump);
    //model_.getAttribute("quality", quality_);
    //setProgressMajorSteps(steps);

    // one step for vars, 2 more for vog
    setProgressMajorSteps(writeInfo_.conditionsOnly ? 1 : 3);
    return !aborted();
}

PWP_BOOL
CaeUnsLociChem::write()
{
    // bool ret = true; // all ok

    // export the grid here

    //-------------------------------------------------------------------------
    // REQUIRED call if MeshLinkSupported is true and IndexScheme is Custom.
    // Implement CaeUnsLociChem::mapIndex(const PWP_UINT64 pwgmNdx) below.
    //
    // For more information, see the "MESHLINK SUPPORT:" and "INDEXING SCHEME:"
    // sections in CaeUnsLociChem::create().
    //-------------------------------------------------------------------------
    //ret = ret && model_.customIndexSchemeReady(*this);

    writeVars();  // always write vars file
    if (!conditionsOnly()) {
        writeVog();
    }
    return !aborted();
}

bool
CaeUnsLociChem::endExport()
{
    return true;
}

// member function to write out vog file
bool CaeUnsLociChem::writeVog()
{
	// one step for streaming faces, one for writing data, one for matrix coloring
	progressBeginStep(3);

	// stream through faces to get face count, and maps face2cell, face2node	
	model_.streamFaces(faceOrder_, *this);
	progressEndStep();

	// color matrix
	colorMatrix();
	progressEndStep();

	// get file name for grid file
	auto dest = exportDestination();
	CAEP_FORMATINFO info;
	PWP_UINT32 ndx = 0;
	PwEnumCaeFormat(ndx++, &info);
	auto handle = PwCreateCaeByName(info.name);
	// use 1 for second extention listed in support files
	auto ext = PwCaeEnumFileExt(handle, 1);
	ostringstream name;
	name << dest << "." << ext;
	string fname = name.str();

	// open HDF5 vog file
	hid_t vogFile = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (vogFile) {
		// write surfaces
		sendInfoMsg("Writing boundary surfaces to vog...");
		writeVogSurfs(vogFile);
		progressIncrement();

		// write nodes
		sendInfoMsg("Writing nodal coordinates to vog...");
		writeVogNodes(vogFile);
		progressIncrement();

		// write faces
		sendInfoMsg("Writing face maps to vog...");
		writeVogFaces(vogFile);
		progressIncrement();
	}
	// close HDF5 vog file
	H5Fclose(vogFile);
	progressEndStep();

	return !aborted();
}

// member function to write out surface tags for vog file
bool CaeUnsLociChem::writeVogSurfs(const hid_t& vogFile) const
{
	CaeUnsPatch patches(model_);
	if (model_.patchCount() > 0) {
		hid_t group_id = H5Gcreate(vogFile, "surface_info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		for (auto ii = 0U; ii < model_.patchCount(); ii++) {

			// get name and id number of boundary condition
			string bcName = "undefined";
			auto bcIdNum = 0;
			bool ret = patches.isValid() && !aborted();
			if (ret) {
				PWGM_CONDDATA cData;
				if (patches.condition(cData)) {
					bcName = cData.name;
					bcIdNum = cData.id;
				}
			}

			// create group for boundary condition and write to file
			hid_t bc_id = H5Gcreate(group_id, bcName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

			hsize_t dims = 1;
			hid_t dataspace_id = H5Screate_simple(1, &dims, NULL);
			hid_t att_id = H5Acreate(bc_id, "Ident", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
			H5Awrite(att_id, H5T_NATIVE_INT, &bcIdNum);
			H5Aclose(att_id);
			H5Gclose(bc_id);

			patches.moveNext();
		}
		H5Gclose(group_id);
	}

	return !aborted();
}

// member function to write out nodal coordinates to vog file
bool CaeUnsLociChem::writeVogNodes(const hid_t& vogFile) const
{
	// open a group for the nodal coordinates
	hid_t group_id = H5Gcreate(vogFile, "node_info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// get a vector of all nodes in the mesh
	auto nodes = getAllNodes();

	// create an HDF5 group for the nodal coordinates and write out
	// create compound datatype for nodal coordinates
	hid_t datatype = H5Tcreate(H5T_COMPOUND, sizeof(node));
	H5Tinsert(datatype, "x", offsetof(struct node, x), H5T_NATIVE_DOUBLE);
	H5Tinsert(datatype, "y", offsetof(struct node, y), H5T_NATIVE_DOUBLE);
	H5Tinsert(datatype, "z", offsetof(struct node, z), H5T_NATIVE_DOUBLE);
	writeVectorToHDF5(group_id, nodes, "positions", datatype);
	H5Gclose(group_id);

	// write out file info
	group_id = H5Gcreate(vogFile, "file_info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hsize_t dims = 1;
	hid_t dataspace_id = H5Screate_simple(1, &dims, NULL);
	hid_t att_id = H5Acreate(group_id, "numNodes", H5T_STD_I64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
	long long nnodes = nodes.size();
	H5Awrite(att_id, H5T_NATIVE_LLONG, &nnodes);
	H5Aclose(att_id);

	H5Gclose(group_id);

	return !aborted();
}

// templated member function to write out a vector of some datatype to an HDF5 file
template<class T>
void CaeUnsLociChem::writeVectorToHDF5(const hid_t& group_id, const vector<T>& vec, const string& name,
	const hid_t& datatype) const
{
	// group_id -- HDF5 group to write data to
	// vec -- vector of some datatype to write out
	// name -- name of group for HDF5 file
	// datatype -- HDF5 datatype of data held in vec

	hsize_t count = vec.size();
	auto rank = 1;
	hid_t dataspace = H5Screate_simple(rank, &count, NULL);
	hsize_t start = 0;
	hsize_t stride = 1;
	hid_t dataset = H5Dcreate(group_id, name.c_str(), datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	if (count != 0) {
		H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL);
		hid_t memspace = H5Screate_simple(rank, &count, NULL);
		H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, &vec[0]);
		H5Sclose(memspace);
	}

	H5Dclose(dataset);
	H5Sclose(dataspace);
	H5Tclose(datatype);
}

// member function to get a copy of all nodal data organized in a vector for purposes
// of writing to HDF5
vector<node> CaeUnsLociChem::getAllNodes() const
{
	// get scale factor for units
	auto fac = getUnitsFactor();

	vector<node> nodes(model_.vertexCount());
	CaeUnsVertex vertex(model_);
	for (auto ii = 0U; ii < nodes.size(); ii++, vertex++) {
		if (vertex.isValid() && !aborted()) {
			PWGM_VERTDATA vdata;
			vertex.dataMod(vdata);
			nodes[ii].x = vdata.x * fac;
			nodes[ii].y = vdata.y * fac;
			nodes[ii].z = vdata.z * fac;
		}
	}

	return nodes;
}

// member funtion to write out face maps to vog file
bool CaeUnsLociChem::writeVogFaces(const hid_t& vogFile)
{
	// get number of faces and cells
	long long numFaces = totalNumFaces_;
	PWGM_ELEMCOUNTS elemCounts;
	long long numCells = model_.elementCount(&elemCounts);

	// write number of faces and number of cells
	hid_t group_id = H5Gopen(vogFile, "file_info", H5P_DEFAULT);
	hsize_t dims = 1;
	hid_t dataspace_id = H5Screate_simple(1, &dims, NULL);
	hid_t att_id = H5Acreate(group_id, "numFaces", H5T_STD_I64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(att_id, H5T_NATIVE_LLONG, &numFaces);
	H5Aclose(att_id);

	att_id = H5Acreate(group_id, "numCells", H5T_STD_I64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(att_id, H5T_NATIVE_LLONG, &numCells);
	H5Aclose(att_id);
	H5Gclose(group_id);

	// write face information
	group_id = H5Gcreate(vogFile, "face_info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// create a set of all face IDs
	set<int> faces;
	for (auto ii = 0U; ii < totalNumFaces_; ii++) {
		faces.insert(ii);
	}

	// write face clusters
	vector<unsigned char> clusterInfo;
	vector<unsigned short> clusterSizes;

	while (!faces.empty()) {
		auto fcluster = faceCluster(faces, clusterInfo, clusterSizes);

		// remove faces written to cluster from set of faces
		for (auto& val : fcluster) {
			faces.erase(val);
		}
	}

	// write face maps to HDF5
	writeVectorToHDF5(group_id, clusterSizes, "cluster_sizes", H5T_NATIVE_USHORT);
	writeVectorToHDF5(group_id, clusterInfo, "cluster_info", H5T_NATIVE_UCHAR);

	H5Gclose(group_id);

	return !aborted();
}

// member function to write vars file
bool CaeUnsLociChem::writeVars() const
{
	// get file name for input file
	auto dest = exportDestination();
	CAEP_FORMATINFO info;
	PWP_UINT32 ndx = 0;
	PwEnumCaeFormat(ndx++, &info);
	auto handle = PwCreateCaeByName(info.name);
	// use 0 for first extention listed in support files
	auto ext = PwCaeEnumFileExt(handle, 0);
	ostringstream name;
	name << dest << "." << ext;
	string fname = name.str();

	PwpFile inpFile(fname.c_str(), pwpWrite | pwpAscii);
	if (inpFile.isOpen()) {
		string comment = "//";
		string dashes = "--------------------";

		ostringstream mline;
		mline << comment << dashes << center("Load Modules", 21) <<
			dashes << "\n";
		auto smline = mline.str();
		inpFile.write(smline.c_str());

		string turbModule = "";
		if (getTurbulenceModule(turbModule)) {
			inpFile.write(turbModule.c_str());
		}
		// if using NRBCs, load module
		string nrbcs = "";
		if (usingNRBCs(nrbcs)) {
			inpFile.write(nrbcs.c_str());
		}

		// open bracket to start vars file
		inpFile.write("{\n");

		// write boundary conditions
		ostringstream bline;
		bline << comment << dashes << center("Boundary Conditions", 21) <<
			dashes << "\n";
		auto sbline = bline.str();
		inpFile.write(sbline.c_str());

		inpFile.write("boundary_conditions: <\n");

		// get patches and write bc for each patch
		if (progressBeginStep(model_.patchCount())) {
			CaeUnsPatch patches(model_);

			for (auto ii = 0U; ii < model_.patchCount(); ii++) {
				progressIncrement();
				const auto isLast = (ii == model_.patchCount() - 1) ? true : false;
				writeBC(patches, isLast, inpFile);
				patches.moveNext();
			}
		}
		inpFile.write(">\n\n");

		ostringstream iline;
		iline << comment << dashes << center("Initial Conditions", 21) <<
			dashes << "\n";
		auto siline = iline.str();
		inpFile.write(siline.c_str());

		inpFile.write(getInputVal("Initial Conditions", "initialConditions", true).c_str());
		inpFile.write(getInputEnum("Grid Coordinates", "gridCoordinates").c_str());
		inpFile.write("\n\n");

		ostringstream pline;
		pline << comment << dashes << center("Physics Models", 21) <<
			dashes << "\n";
		auto spline = pline.str();
		inpFile.write(spline.c_str());

		inpFile.write(getInputVal("Chemistry Model", "chemistry_model").c_str());
		inpFile.write(getInputEnum("Thermodynamic Model", "thermodynamic_model").c_str());
		inpFile.write(getInputEnum("Transport Model", "transport_model").c_str());
		inpFile.write(getInputEnum("Diffusion Model", "diffusion_model").c_str());
		inpFile.write(getInputEnum("Turbulence Model", "turbulence_model").c_str());
		inpFile.write("\n");
		inpFile.write(getInputEnum("Inviscid Flux", "inviscidFlux").c_str());
		inpFile.write(getInputEnum("Limiter", "limiter").c_str());
		inpFile.write(getInputEnum("Preconditioning Factor", "Minf").c_str());
		if (getInputBool("Hybrid RANS/LES")) {
			inpFile.write(getInputEnum("Multiscale", "multi_scale").c_str());
		}
		inpFile.write("\n\n");

		ostringstream tline;
		tline << comment << dashes << center("Time Incrementation", 21) <<
			dashes << "\n";
		auto stline = tline.str();
		inpFile.write(stline.c_str());

		inpFile.write(getInputEnum("Time Integration", "time_integration").c_str());
		inpFile.write(getInputVal("Newton Iterations", "newton_iter").c_str());
		inpFile.write(getInputEnum("Fluid Linear Solver", "fluidLinearSolver").c_str());
		inpFile.write(getInputVal("Gauss Seidel Iterations", "gauss_seidel_iter").c_str());
		inpFile.write(getInputVal("CFL Max", "cflmax").c_str());
		inpFile.write(getInputVal("Urelax", "urelax").c_str());
		inpFile.write(getInputVal("Global Time Step", "dtmax").c_str());
		inpFile.write(getInputVal("Iterations", "stop_iter").c_str());
		inpFile.write("\n\n");

		ostringstream oline;
		oline << comment << dashes << center("Output Control", 21) <<
			dashes << "\n";
		auto soline = oline.str();
		inpFile.write(soline.c_str());

		inpFile.write(getInputVal("Print Frequency", "print_freq").c_str());
		inpFile.write(getInputVal("Plot Output", "plot_output", false).c_str());
		inpFile.write(getInputVal("Plot Frequency", "plot_freq").c_str());
		inpFile.write(getInputVal("Plot Modulo", "plot_modulo").c_str());
		inpFile.write(getInputVal("Restart Frequency", "restart_freq").c_str());
		inpFile.write(getInputVal("Restart Modulo", "restart_modulo").c_str());
		inpFile.write("\n\n");


		// close bracket to end vars file
		inpFile.write("}\n");

		progressEndStep();
	}
	return !aborted();
}

// function to center a string within a space of a given size
string center(const string& str, const int& sz)
{
	auto needToAdd = sz - str.length();
	auto prefix = needToAdd / 2 + (needToAdd % 2);
	auto postfix = needToAdd / 2;

	ostringstream line;
	for (auto ii = 0U; ii < prefix; ii++) {
		line << " ";
	}
	line << str;
	for (auto ii = 0U; ii < postfix; ii++) {
		line << " ";
	}
	return line.str();
}

// function to out boundary conditions to vars file
bool CaeUnsLociChem::writeBC(const CaeUnsPatch& patch, const bool& isLast, PwpFile& inpFile) const
{
	bool ret = patch.isValid() && !aborted();
	if (ret) {
		PWGM_CONDDATA cData;
		if (patch.condition(cData)) {
			ostringstream bcline;
			bcline << "    " << cData.name << " = " << cData.type << "()";
			auto sbcline = bcline.str();
			inpFile.write(sbcline.c_str());

			if (isLast) {
				inpFile.write("\n");
			}
			else {
				inpFile.write(",\n");
			}
		}
	}
	return ret;
}

// utility routine for creating face cluster
// from FVMGridWriter.cc in Loci
void writeUnsignedVal(vector<unsigned char>& cluster, unsigned long long val)
{
	do {
		unsigned char byte = val & 0x7f;
		val = val >> 7;
		if (val != 0) {
			byte |= 0x80;
		}
		cluster.push_back(byte);
	} while (val != 0);
}

// utility routine for creating face cluster
// from FVMGridWriter.cc in Loci
void writeSignedVal(vector<unsigned char>& cluster, long long val)
{
	bool sign = false;
	if (val < 0) {
		sign = true;
		val = -val;
	}
	unsigned char byte = val & 0x3f;
	if (sign) {
		byte |= 0x40;
	}
	val = val >> 6;
	if (val != 0) {
		byte |= 0x80;
	}
	cluster.push_back(byte);
	if ((byte & 0x80) == 0x80) {
		writeUnsignedVal(cluster, val);
	}
}

// utility routine for creating face cluster
// from FVMGridWriter.cc in Loci
void writeTable(vector<unsigned char>& cluster, set<int>& dataset)
{
	auto iter = dataset.begin();
	unsigned char sz = dataset.size();
	cluster.push_back(sz);
	writeSignedVal(cluster, *iter);
	long long last = *iter;
	for (++iter; iter != dataset.end(); ++iter) {
		unsigned long diff = *iter - last;
		last = *iter;
		writeUnsignedVal(cluster, diff);
	}
}

// member function to encode face cluster for writing out to HDF5
// modified from FVMGridWriter.cc in Loci
vector<unsigned char>
CaeUnsLociChem::encodeFaceCluster(set<int>& fcluster, set<int>& nodeSet, set<int>& cellSet) const
{
	// create node and cell local maps
	map<int, int> node2local;
	map<int, int> cell2local;
	auto cnt = 0;
	for (auto& nid : nodeSet) {
		node2local[nid] = cnt++;
	}
	cnt = 0;
	for (auto cid : cellSet) {
		cell2local[cid] = cnt++;
	}


	// sort faces according to number of nodes
	vector<pair<int, int>> faceOrder(fcluster.size());
	cnt = 0;
	for (auto& fid : fcluster) {
		faceOrder[cnt++] = std::make_pair(getFaceSize(fid), fid);
	}
	std::sort(faceOrder.begin(), faceOrder.end());

	constexpr auto maxSize = 255;

	vector<pair<int, int>> rll;
	auto size = faceOrder[0].first;
	cnt = 0;
	for (auto ii = 0U; ii < faceOrder.size(); ++ii) {
		if (size != faceOrder[ii].first) {
			while (cnt > maxSize) {
				rll.push_back(std::make_pair(size, maxSize));
				cnt -= maxSize;
			}
			rll.push_back(std::make_pair(size, cnt));
			cnt = 0;
			size = faceOrder[ii].first;
		}
		cnt++;
	}

	while (cnt > maxSize) {
		rll.push_back(std::make_pair(size, maxSize));
		cnt -= maxSize;
	}
	rll.push_back(std::make_pair(size, cnt));

	// write out faces for each size category
	// face data consists of the node ids which make it up, and then the cl and cr ids
	vector<unsigned char> cluster;
	cnt = 0;
	for (auto ii = 0U; ii < rll.size(); ++ii) {
		cluster.push_back(rll[ii].first);
		cluster.push_back(rll[ii].second);
		auto nds = rll[ii].first;

		for (auto kk = 0; kk < rll[ii].second; ++kk) {
			auto fc = faceOrder[cnt].second;
			cnt++;
			for (auto jj = 0; jj < nds; ++jj) {
				cluster.push_back(node2local[face2node_[fc][jj]]);
			}
			cluster.push_back(cell2local[face2cell_[fc].first]);
			cluster.push_back(cell2local[face2cell_[fc].second]);
		}
	}

	// a zero face size marks end of cluster
	cluster.push_back(0);

	writeTable(cluster, nodeSet);
	writeTable(cluster, cellSet);

	return cluster;
}

// member function to encode face cluster for writing out to HDF5
// modified from FVMGridWriter.cc in Loci
set<int> CaeUnsLociChem::faceCluster(set<int>& faces, vector<unsigned char>& clusterInfo,
	vector<unsigned short>& clusterSizes) const
{
	set<int> fcluster;
	set<int> nodeSet;
	set<int> cellSet;

	auto nnodes = 0;
	auto ncells = 0;

	// loop over leftover faces here
	// construct cell and node sets and face cluster
	for (auto& ii : faces) {
		auto nnodesCluster = 0;
		auto ncellsCluster = 0;
		if (cellSet.find(face2cell_[ii].first) == cellSet.end()) {
			ncellsCluster++;
		}
		if (cellSet.find(face2cell_[ii].second) == cellSet.end()) {
			ncellsCluster++;
		}

		auto fsize = getFaceSize(ii);
		for (auto jj = 0; jj < fsize; ++jj) {
			if (nodeSet.find(face2node_[ii][jj]) == nodeSet.end()) {
				nnodesCluster++;
			}
		}

		if (nnodes + nnodesCluster > 256 || ncells + ncellsCluster > 256) {
			break;
		}
		cellSet.insert(face2cell_[ii].first);
		cellSet.insert(face2cell_[ii].second);

		for (auto jj = 0; jj < fsize; jj++) {
			nodeSet.insert(face2node_[ii][jj]);
		}

		fcluster.insert(ii);
		ncells = cellSet.size();
		nnodes = nodeSet.size();
	}

	auto cluster = encodeFaceCluster(fcluster, nodeSet, cellSet);

	clusterSizes.push_back(cluster.size());
	for (auto ii = 0U; ii < cluster.size(); ++ii) {
		clusterInfo.push_back(cluster[ii]);
	}

	return fcluster;
}

// member function to get the number of nodes on a face
// for a 3D mesh in Pointwise, cell faces are either triangles or quads
// when constructing face2node map, 4th node was set to -1 if it was a triangular face
int CaeUnsLociChem::getFaceSize(const int& index) const
{
	return (face2node_[index][3] == -1) ? 3 : 4;
}
PWGM_ENUM_ELEMTYPE CaeUnsLociChem::getFaceType(const int& index) const
{
	return (face2node_[index][3] == -1) ? PWGM_ELEMTYPE_TRI : PWGM_ELEMTYPE_QUAD;
}

// member function to get the units factor by which to scale the grid coordinates
// Loci/Chem wants the grid in meters, so if another unit system is specified, convert
double CaeUnsLociChem::getUnitsFactor() const
{
	const char* units;
	model_.getAttributeEnum("Grid Units", units, "fail");
	string sunits = units;

	auto lref = 0.0;
	if (sunits == "meters") {
		lref = 1.0;
	}
	else if (sunits == "centimeters") {
		lref = 0.01;
	}
	else if (sunits == "millimeters") {
		lref = 0.001;
	}
	else if (sunits == "feet") {
		lref = 0.3048;
	}
	else if (sunits == "inches") {
		lref = 0.0254;
	}

	ostringstream msgStream;
	msgStream << "Converting grid from " << sunits << " to meters with factor " << lref;
	auto msg = msgStream.str();
	sendInfoMsg(msg.c_str());

	return lref;
}

// member function to get the input specified from a solver attributes enum
string CaeUnsLociChem::getInputEnum(const string& param, const string& name) const
{
	const char* val;
	model_.getAttributeEnum(param.c_str(), val, "fail");

	ostringstream inputVar;
	inputVar << name << ": " << val << "\n";
	return inputVar.str();
}

// member function to get the input specified from a solver attributes integral value
string CaeUnsLociChem::getInputVal(const string& param, const string& name, bool bracket) const
{
	const char* val;
	model_.getAttribute(param.c_str(), val);

	ostringstream inputVar;
	inputVar << name << ": ";
	if (bracket) {
		inputVar << "<";
	}
	inputVar << val;
	if (bracket) {
		inputVar << ">";
	}
	inputVar << "\n";
	return inputVar.str();
}

// member function to get the input specified from a solver attributes bool
bool CaeUnsLociChem::getInputBool(const string& param) const
{
	const char* val;
	model_.getAttribute(param.c_str(), val);
	string inputBool = val;
	return (inputBool == "true") ? true : false;
}

// member function to determine if a module needs to be loaded for using a turbulence model
bool CaeUnsLociChem::getTurbulenceModule(string& line) const
{
	const char* val;
	model_.getAttributeEnum("Turbulence Model", val, "fail");
	string turbModel = val;

	auto needModule = false;

	if (turbModel == "SST" || turbModel == "BSL" || turbModel == "Wilcox08" ||
		turbModel == "Wilcox98" || turbModel == "KW") {
		line = "loadModule: KOmegaModel\n";
		needModule = true;
	}
	else if (turbModel == "RKE") {
		line = "loadModule: RKEModel\n";
		needModule = true;
	}
	else if (turbModel == "Spalart_Allmaras") {
		line = "loadModule: SAModel\n";
		needModule = true;
	}
	else if (turbModel == "wale" || turbModel == "smagorinsky") {
		line = "loadModule: subGridScale\n";
		needModule = true;
	}

	return needModule;
}

bool CaeUnsLociChem::usingNRBCs(string& line) const
{
	CaeUnsPatch patches(model_);
	auto usingNRBC = false;
	for (auto ii = 0U; ii < model_.patchCount(); ii++) {
		bool ret = patches.isValid() && !aborted();
		if (ret) {
			PWGM_CONDDATA cData;
			if (patches.condition(cData)) {
				ostringstream bcTypeStream;
				bcTypeStream << cData.type;
				auto bcType = bcTypeStream.str();
				if (bcType == "inflowNRBC" || bcType == "outflowNRBC") {
					line = "loadModule: nonreflecting\n";
					usingNRBC = true;
					break;
				}
			}
		}
		patches.moveNext();
	}
	return usingNRBC;
}

void CaeUnsLociChem::swapOrientation(const int& ind, const PWGM_ENUM_ELEMTYPE& etype)
{
	// swap cl, cr maps
	std::swap(face2cell_[ind].first, face2cell_[ind].second);

	// reverse node orientation
	if (etype == PWGM_ELEMTYPE_QUAD) {
		std::swap(face2node_[ind][0], face2node_[ind][3]);
		std::swap(face2node_[ind][1], face2node_[ind][2]);
	}
	else if (etype == PWGM_ELEMTYPE_TRI) {
		std::swap(face2node_[ind][0], face2node_[ind][2]);
	}
}

// This member function mimics the behavior of VOG::colorMatrix found in vogtools.cc.
// The purpose of this function is to orient the face maps in a way that the Loci/Chem
// matrix solvers expect them to be.
void CaeUnsLociChem::colorMatrix()
{
	sendInfoMsg("coloring matrix...");

	// VOG::colorMatrix is written for parallel processing. This plugin is a serial
	// program, so many steps can be skipped. For example the getDist() function is 
	// not needed because all entities are local

	// create cell map
	vector<pair<int, int>> cellmap(numInteriorFaces_ * 2);
	auto cnt = 0;
	for (auto& faceCellMap : face2cell_) {
		if (isInterior(faceCellMap)) {
			cellmap[cnt++] = std::make_pair(faceCellMap.first, faceCellMap.second);
			cellmap[cnt++] = std::make_pair(faceCellMap.second, faceCellMap.first);
		}
	}

	// create inverse map
	map<int, vector<int>> c2c;
	inverseMap(c2c, cellmap);

	// initialize color vector
	PWGM_ELEMCOUNTS elemCounts;
	auto numCells = model_.elementCount(&elemCounts);
	vector<int> color(numCells, -1);

	// create set of cells
	set<int> left_out;
	set<int> cells;  // geom_cells in vogtools.cc
	for (auto ii = 0U; ii < numCells; ii++) {
		left_out.insert(ii);
		cells.insert(ii);
	}

	// compared to vogtools.cc ctmp and color are the same because we are only using
	// one processor. There is no difference between a store and a distributed store
	auto col = 0;
	vector<int> visited;
	auto lo_p = 0;
	while (!left_out.empty()) {
		vector<int> work;
		work.push_back(*std::min_element(left_out.begin(), left_out.end()));

		while (work.size() != 0) {
			vector<int> working;
			for (auto ii = 0; ii < work.size(); ++ii) {
				auto cc = work[ii];
				if (color[cc] == -1) {
					color[cc] = col++;
					visited.push_back(cc);
					for (auto pi = c2c[cc].begin(); pi != c2c[cc].end(); ++pi) {
						if (cells.find(*pi) != cells.end() && color[*pi] == -1) {
							working.push_back(*pi);
						}
					}
				}
			}
			work.swap(working);
		}
		left_out.clear();
		for (auto ii = lo_p; ii < static_cast<int>(numCells); ii++) {
			if (color[ii] == -1) {
				left_out.insert(ii);
				break;
			}
		}
		if (!left_out.empty()) {
			lo_p = *std::min_element(left_out.begin(), left_out.end());
		}
	}

	// clone stuff here in vogtools.cc -- don't think this is needed

	// loop over interior faces and swap orientation if necessary
	for (auto ii = 0U; ii < totalNumFaces_; ii++) {
		if (isInterior(ii)) {
			auto color_l = color[face2cell_[ii].first];
			auto color_r = color[face2cell_[ii].second];
			if (color_l == color_r) {
				ostringstream err;
				err << "Left and right colors are both equal to " << color_l;
				auto serr = err.str();
				sendErrorMsg(serr.c_str());
			}
			if (color_l == -1 || color_r == -1) {
				sendErrorMsg("matrix coloring internal error");
			}

			if (color_l > color_r) {
				// change face orientation to match matrix coloring
				swapOrientation(ii, getFaceType(ii));
			}
		}
	}
}

// function to create an inverse map as outlined in distributed_inversMap in 
// distribute_container.cc this is a heavily modified version because many 
// things do not need to be done since we are only using one processor and the 
// inverse map is therefore not distributed
void inverseMap(map<int, vector<int>>& result, vector<pair<int, int>>& input) {

	// sort input according to second field
	std::sort(input.begin(), input.end(), []
	(const pair<int, int>& p1, const pair<int, int>& p2)
		{return p1.second < p2.second; });

	// adapted from distribute_container.cc
	// since this is not an MPI program, there is no sending or receiving
	// the send and recv variables are equal
	// count what to "send"
	auto size_recv = 0;
	for (auto ii = 0U; ii < input.size(); ++ii) {
		size_recv += 2;
	}

	// since only using one processor, send_displacement is always 0, as is current_p
	vector<int> recv_store(size_recv, 0);
	auto offset = 0;
	for (auto ii = 0U; ii < input.size(); ++ii) {
		auto to = input[ii].second;
		auto from = input[ii].first;
		recv_store[offset++] = to;
		recv_store[offset++] = from;
	}

	// local input image stuff -- don't think this is needed

	// populate sizes of vector<int> in map
	auto localInput = input;
	vector<int> sizes(localInput.size(), 0);
	for (auto ii = 0; ii < size_recv; ++ii) {
		auto indx = recv_store[ii++];
		sizes[indx]++;
	}

	// allocate vector<int>'s in result map
	for (auto ii = 0; ii < size_recv; ++ii) {
		auto indx = recv_store[ii++];
		result[indx].resize(sizes[indx], 0);
	}

	// create inverse map
	for (auto ii = 0; ii < size_recv; ++ii) {
		auto indx = recv_store[ii++];
		auto refer = recv_store[ii];
		sizes[indx]--;
		result[indx][sizes[indx]] = refer;
	}
}

// member functions to determine if face is interior or not
// for boundary faces, cr map will always be negative
bool CaeUnsLociChem::isInterior(const int& ind) const
{
	return (face2cell_[ind].second >= 0) ? true : false;
}
bool CaeUnsLociChem::isInterior(const pair<int, int>& facemap) const
{
	return (facemap.second >= 0) ? true : false;
}

//===========================================================================
// face streaming handlers
//===========================================================================

PWP_UINT32
CaeUnsLociChem::streamBegin(const PWGM_BEGINSTREAM_DATA &data)
{
	// get total number of faces and allocate face maps
	totalNumFaces_ = data.totalNumFaces;
	numInteriorFaces_ = 0;
	face2cell_.resize(totalNumFaces_);
	face2node_.resize(totalNumFaces_);
	return progressBeginStep(totalNumFaces_);
}

PWP_UINT32
CaeUnsLociChem::streamFace(const PWGM_FACESTREAM_DATA &data)
{
	// get face id and vertex information
	const auto fid = data.face;
	const auto nodes = data.elemData.vert;

	// construct face-to-node maps
	if (data.elemData.type == PWGM_ELEMTYPE_QUAD) {
		face2node_[fid][0] = nodes[0].id;
		face2node_[fid][1] = nodes[1].id;
		face2node_[fid][2] = nodes[2].id;
		face2node_[fid][3] = nodes[3].id;
	}
	else if (data.elemData.type == PWGM_ELEMTYPE_TRI) {
		face2node_[fid][0] = nodes[0].id;
		face2node_[fid][1] = nodes[1].id;
		face2node_[fid][2] = nodes[2].id;
		face2node_[fid][3] = -1;
	}
	else {
		sendErrorMsg("Error: Found non tri/quad face type");
	}

	// construct face-to-cell maps
	auto cl = data.owner.cellIndex;
	// possible that cr is a ghost cell, if so specify as negative of boundary id
	auto cr = 0;
	if (data.type == PWGM_FACETYPE_BOUNDARY) {
		PWGM_CONDDATA cData;
		PwDomCondition(data.owner.domain, &cData);
		cr = cData.id * -1;
	}
	else {
		cr = data.neighborCellIndex;
		numInteriorFaces_++;  // increment interior face counter
	}
	face2cell_[fid] = std::make_pair(cr, cl);

	// swap order on boundary faces so that area vector as calculated by Loci/Chem
	// will point from left cell center to right cell center
	if (data.type == PWGM_FACETYPE_BOUNDARY) {
		swapOrientation(fid, data.elemData.type);
	}

	return progressIncrement();
}

PWP_UINT32
CaeUnsLociChem::streamEnd(const PWGM_ENDSTREAM_DATA &data)
{
	(void)data.ok; // silence unused param warning
	return progressEndStep();
}

//===========================================================================
// custom index scheme mapper
//===========================================================================

PWP_UINT64
CaeUnsLociChem::mapIndex(const PWP_UINT64 pwgmNdx)
{
    // If MeshLinkSupported is true and IndexScheme is Custom, you must
    // implement this method. Invoked by calling customIndexSchemeReady(*this)
    // from CaeUnsLociChem::write().
    //
    // For more information, see the "MESHLINK SUPPORT:" and "INDEXING SCHEME:"
    // sections in CaeUnsLociChem::create().
    return pwgmNdx;
}


//===========================================================================
// called ONCE when plugin first loaded into memory
//===========================================================================

bool
CaeUnsLociChem::create(CAEP_RTITEM &rti)
{
    (void)rti.BCCnt; // silence unused arg warning
    bool ret = true;

    //-----------------------------------------------------------------------
    // BYTE ORDERING:
    //   Set the following flags to control the byte ordering options
    //   supported by the solver. If all flags are false, the plugin will use
    //   the platform's native byte ordering. Currently, Pointwise only runs on
    //   little endian, intel platforms. If the solver targeted by this plugin
    //   cannot import little endian files, you must set bigEndian to true and
    //   littleEndian to false.
    //-----------------------------------------------------------------------
    bool bigEndian = false;
    bool littleEndian = false;
    if (ret && (bigEndian || littleEndian)) {
        ret = allowByteOrders(rti, bigEndian, littleEndian);
    }

    //-----------------------------------------------------------------------
    // ELEMENT TOPOLOGY:
    //   Set the following flags to control the element topology options
    //   supported by the solver. If all flags are false, the allowed element
    //   topologies will be inferred from the supported element types. Unless
    //   this plugin has special needs, you should leave these all false.
    //-----------------------------------------------------------------------
    bool structured = false;
    bool unstructured = false;
    bool prismatic = false;
    if (ret && (structured || unstructured || prismatic)) {
        ret = allowElementTopologies(rti, structured, unstructured, prismatic);
    }

    //-----------------------------------------------------------------------
    // MESHLINK SUPPORT:
    //   Uncomment the line below if this plugin supports mesh link exporting.
    //   If "MeshLinkSupported" is not set, mesh link exporting is unsupported.
    //-----------------------------------------------------------------------
    //ret = ret && meshLinkSupported(rti, true);

    //-----------------------------------------------------------------------
    // INDEXING SCHEME:
    //   Uncomment the appropriate line below to specify the vertex index
    //   scheme used by the solver. If "IndexScheme" is not set, ZeroBased
    //   indexing is used. This setting is needed for proper mesh link export.
    //
    //   ZeroBased -
    //     Vertices are written to the mesh file in the same order as returned
    //     by PwModEnumVertices(). In the mesh file the vertices are indexed as
    //     0..NumVerts-1.
    //
    //   OneBased -
    //     Vertices are written to the mesh file in the same order as returned
    //     by PwModEnumVertices(). In the mesh file the vertices are indexed as
    //     1..NumVerts.
    //
    //   Custom - 
    //     Vertices are NOT written to the mesh file in the same order as
    //     returned by PwModEnumVertices(). In this case, the plugin must build
    //     an appropriate index map during export. After the mesh export is
    //     completed but before returning from runtimeWrite(), the plugin must
    //     call PwModCustomIndexSchemeReady(mapperCallback) to notify the
    //     Pointwise export framework that the vertex map is ready. Before
    //     returning from PwModCustomIndexSchemeReady(), the Pointwise export
    //     framework will invoke the mapperCallback as needed to map a PWGM
    //     index to its corresponding exported index.
    //     See PwModCustomIndexSchemeReady() for details.
    //-----------------------------------------------------------------------
    //ret = ret && useZeroBasedIndexing(rti);
    //ret = ret && useOneBasedIndexing(rti);
    //ret = ret && useCustomIndexing(rti);

    //-----------------------------------------------------------------------
    // NON-INFLATED BC TYPE NAMES:
    //   Specify which of the defined BCs types are non-inflated. Normally, if a
    //   BC is applied to a floating baffle domain, or block-to-block connection
    //   domain, the interior points are "inflated" to create a zero thickness
    //   wall. However, if a non-inflated BC is used, the points are NOT
    //   inflated. These BC types are often used for flow-through or porous
    //   conditions. Specify the names as a vbar delimited list of non-inflated
    //   BC type names (e.g. "bc_type_name1|bc_type_name2|bc_type_name3").
    //-----------------------------------------------------------------------
    const char * const shadowTypes = nullptr;
    if (ret && (nullptr != shadowTypes) && ('\0' != shadowTypes[0])) {
        ret = shadowBcTypes(rti, shadowTypes);
    }

    // These attributes are for example only. You can publish any attribute
    // needed for your solver.
    ret = ret &&
    //     publishUIntValueDef(rti, "iterations", 5, "Number of iterations", 0,
    //          2000) &&
    //     publishIntValueDef(rti, "magnitude", -5, "Signed int magnitude",
    //          -100, 100) &&
    //     publishRealValueDef(rti, "mach", 0.3, "Incoming flow velocity", 0.0,
    //          1000.0, 0.0, 50.0) &&
    //     publishRealValueDef(rti, "temperature", 77.5, "Ambient temperature",
    //          -5000.0, 5000.0, -100.0, 3000.0) &&
    //     publishEnumValueDef(rti, "temperature.units", "Fahrenheit",
    //          "Grid temperature units", "Fahrenheit|Celsius") &&
    //     publishEnumValueDef(rti, "units", "Inches", "Grid dimensional units",
    //          "Yards|Inches|Meters|Millimeters") &&
    //     publishStringValueDef(rti, "description", "", "Grid description") &&
    //     publishBoolValueDef(rti, "linear", false, "Grid is linear",
    //          "reject|accept");

	publishUIntValueDef(rti, "Iterations", 1, "Maximum number of iterations", 1,
		1000000) &&
	publishUIntValueDef(rti, "Plot Frequency", 0, "Plot frequency: How often to output data", 0,
		1000000) &&
	publishUIntValueDef(rti, "Plot Modulo", 0, "Plot modulo: When to overwrite plot data", 0,
		1000000) &&
	publishUIntValueDef(rti, "Restart Frequency", 0, "Restart frequency: How often to write restart data", 0,
		1000000) &&
	publishUIntValueDef(rti, "Restart Modulo", 0, "Restart modulo: When to overwrite restart data", 0,
		1000000) &&
	publishUIntValueDef(rti, "Print Frequency", 10, "Print frequency: How often to output data", 0,
		1000000) &&
	publishUIntValueDef(rti, "Newton Iterations", 1, "Number of nonlinear iterations", 1,
		1000) &&
	publishUIntValueDef(rti, "Gauss Seidel Iterations", 5, "Number of Gauss Seidel iterations", 1,
		1000) &&

	publishRealValueDef(rti, "Global Time Step", 1.0, "Global time step dtmax",
		1.0e-20, 1.0e20, 1.0e-8, 10.0) &&
	publishRealValueDef(rti, "CFL Max", 0.0, "Time step parameter",
		0.0, 1.0e20, 0.0, 1.0e5) &&
	publishRealValueDef(rti, "Urelax", 1.0, "Solution relaxation parameter",
		0.0, 1.0, 0.0, 1.0) &&
	publishRealValueDef(rti, "Preconditioning Factor", 1.0, "Preconditioning factor Minf",
		0.0, 100.0, 0.0, 1.0) &&

	publishEnumValueDef(rti, "Grid Units", "meters", "Grid dimensional units",
		"meters|centimeters|millimeters|feet|inches") &&
	publishEnumValueDef(rti, "Grid Coordinates", "cartesian", "Grid coordinate system",
		"cartesian|axisymmetric") &&
	publishEnumValueDef(rti, "Thermodynamic Model", "adaptive", "Thermodynamic model used in simulation",
		"adaptive|vibrational_equilibrium|curve_fit") &&
	publishEnumValueDef(rti, "Transport Model", "none", "Transport model used in simulation",
		"none|transportDB|chemkin|sutherland|powerLaw|const_viscosity") &&
	publishEnumValueDef(rti, "Diffusion Model", "none", "Diffusion model used in simulation",
		"none|transportDB|chemkin|laminarSchmidt|const_diffusivity") &&
	publishEnumValueDef(rti, "Turbulence Model", "none", "Turbulence model used in simulation",
		"none|SST|BSL|Wilcox08|Wilcox98|KW|RKE|Spalart_Allmaras|wale|smagorinsky") &&
	publishEnumValueDef(rti, "Inviscid Flux", "default", "Inviscid flux method",
		"default|roe|roeBTW|roePT|roeCT|hlle|hllc") &&
	publishEnumValueDef(rti, "Fluid Linear Solver", "sgs", "Matrix solver for implicit algorithm",
		"sgs|fsgs|lsgs|petsc") &&
	publishEnumValueDef(rti, "Time Integration", "euler", "Method of time incrementation",
		"euler|steady_state|time_accurate") &&
	publishEnumValueDef(rti, "Limiter", "venkatakrishnan", "Limiter used in simulation",
		"venkatakrishnan|none|barth|NB|V2|zero") &&
	publishEnumValueDef(rti, "Multiscale", "none", "multi_scale option",
		"none|LES|LES2DZ|LNS|DES|DDES|IDDES") &&

	publishStringValueDef(rti, "Plot Output", "put", "Additional vairables to output") &&
	publishStringValueDef(rti, "Initial Conditions", "", "Initial conditions for simulation") &&
	publishStringValueDef(rti, "Chemistry Model", "air_1s0r", "Chemistry model used in simulation") &&

	publishBoolValueDef(rti, "Hybrid RANS/LES", false, "Flag to use multi_scale option",
		"no|yes");

    return ret;
}


//===========================================================================
// called ONCE just before plugin unloaded from memory
//===========================================================================

void
CaeUnsLociChem::destroy(CAEP_RTITEM &rti)
{
    (void)rti.BCCnt; // silence unused arg warning
}
