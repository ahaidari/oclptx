/*  Copyright (C) 2014
 *    Afshin Haidari
 *    Steve Novakov
 *    Jeff Taylor
 */

/* samplemanager.h
 *
 *
 * Part of
 *    oclptx
 * OpenCL-based, GPU accelerated probtrackx algorithm module, to be used
 * with FSL - FMRIB's Software Library
 *
 * This file is part of oclptx.
 *
 * oclptx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * oclptx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with oclptx.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef  OCLPTX_SAMPLEMANAGER_H_
#define  OCLPTX_SAMPLEMANAGER_H_

#include "newimage/newimageall.h"
#include "oclptxOptions.h"
#include "customtypes.h"
#include <iostream>
#include <vector>
#include <string>
using namespace std;
using namespace NEWIMAGE;

#define __CL_ENABLE_EXCEPTIONS
// adds exception support from CL libraries
// define before CL headers inclusion


class SampleManager{
	public:
      static SampleManager& GetInstance();
      ~SampleManager();
      
	  //As of February 19th, you require 2 arguments to parse: a set of data, and a mask.
	  //Currently, you can load a preset mask (we are missing the initializing seedmask in our sample space, so for now,
	  //We just use a given mask, although the concept is the same. That volume is stored in a member variable here.
      void ParseCommandLine(int argc, char** argv);
      void LoadBedpostData(const std::string& aBasename);

      //Getters
      float const GetThetaData(int aFiberNum, int aSamp, int aX, int aY, int aZ);
      float const GetPhiData(int aFiberNum, int aSamp, int aX, int aY, int aZ);
      float const GetfData(int aFiberNum, int aSamp, int aX, int aY, int aZ);
      const unsigned short int* GetSeedMaskToArray();
      const volume<short int>* GetSeedMask();
      
      //WARNING: If you use these getters, you must access data from the BedpostXData vector as follows:
      //Ex Theta: thetaData.data.at(aFiberNum)[(aSamp)*(nx*ny*nz) + (aZ)*(nx*ny) + (aY)*nx + (aX)]
      //Where aFiberNum = the Fiber (0 or 1), aSamp = SampleNumber, nx ny nz = spacial dimensions stored in BedpostXData,
      //and aY aX, aZ = inputed spacial coordinates.
      //See definition of GetThetaData(...) above for example. 
      const BedpostXData* GetThetaDataPtr();
      const BedpostXData* GetPhiDataPtr();
      const BedpostXData* GetFDataPtr(); 

	private:
      SampleManager();

      void LoadBedpostDataHelper(
        const std::string& aThetaSampleName,
        const std::string& aPhiSampleName,
        const std::string& afSampleName,
        const volume<float>& aMask = volume<float>(),
        const int aFiberNum = 0);
        
      void PopulateMemberParameters(
        const volume4D<float> aLoadedData,
        BedpostXData& aTargetContainer,
        const volume<float> aMaskParams,
        const int aFiberNum);

      std::string IntTostring(const int& value);

      //Members
      static SampleManager* _manager;

      oclptxOptions& _oclptxOptions;
      volume<short int> _seedMask;
      BedpostXData _thetaData;
      BedpostXData _phiData;
      BedpostXData _fData;
};

#endif

//EOF
