/*  Copyright (C) 2014
 *    Afshin Haidari
 *    Steve Novakov
 *    Jeff Taylor
 */

/* samplemanager.cc
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



#include <fstream>
#include <sstream>


#define __CL_ENABLE_EXCEPTIONS
// adds exception support from CL libraries
// define before CL headers inclusion

#include "samplemanager.h"
#include "oclptxOptions.h"

//
// Assorted Functions Declerations
//

std::string SampleManager::IntTostring(const int& value)
{
   std::stringstream s;
   s << value;
   return s.str();
}

//Private method: Used for loading data directly into member containers
void SampleManager::LoadBedpostDataHelper(
  const std::string& aThetaSampleName,
  const std::string& aPhiSampleName,
  const std::string& afSampleName,
  const volume<float>& aMask,
  const int aFiberNum  )
{
   volume4D<float> loadedVolume4DTheta;
   volume4D<float> loadedVolume4DPhi;
   volume4D<float> loadedVolume4Df;

   //Load Theta/Phi/f samples
   NEWIMAGE::read_volume4D(loadedVolume4DTheta, aThetaSampleName);
   NEWIMAGE::read_volume4D(loadedVolume4DPhi, aPhiSampleName);
   NEWIMAGE::read_volume4D(loadedVolume4Df, afSampleName);

   if(aMask.xsize() > 0)
   {
      //_thetaSamples.push_back(loadedVolume4DTheta.matrix(aMask));
      //_phiSamples.push_back(loadedVolume4DPhi.matrix(aMask));
      //_fSamples.push_back(loadedVolume4Df.matrix(aMask));
   }
   else
   {
      PopulateMemberParameters(loadedVolume4DTheta, _thetaData, loadedVolume4DTheta[0], aFiberNum);
      PopulateMemberParameters(loadedVolume4DPhi, _phiData, loadedVolume4DPhi[0], aFiberNum);
      PopulateMemberParameters(loadedVolume4Df, _fData, loadedVolume4Df[0], aFiberNum);
   }
}

void SampleManager::PopulateMemberParameters(
  const volume4D<float> aLoadedData,
  BedpostXData& aTargetContainer,
  const volume<float> aMaskParams,
  const int aFiberNum)
{
  
   const int ns = aLoadedData.tsize();
   const int nx = aLoadedData.xsize();
   const int ny = aLoadedData.ysize();
   const int nz = aLoadedData.zsize();
   
   aTargetContainer.data.push_back( new float[ns*nx*ny*nz] );
   aTargetContainer.nx = nx;
   aTargetContainer.ny = ny;
   aTargetContainer.nz = nz;
   aTargetContainer.ns = ns;
   
   int xoff = aLoadedData[0].minx() - aMaskParams.minx();
   int yoff = aLoadedData[0].miny() - aMaskParams.miny();
   int zoff = aLoadedData[0].minz() - aMaskParams.minz();
   for (int z = aMaskParams.minz(); z <= aMaskParams.maxz(); z++)
   {
     for (int y = aMaskParams.miny(); y <= aMaskParams.maxy(); y++)
     {
        for (int x = aMaskParams.minx(); x <= aMaskParams.maxx(); x++)
        {
           if (aMaskParams(x,y,z) > 0)
           {
              for (int t = aLoadedData.mint(); t <= aLoadedData.maxt(); t++)
              {
                 aTargetContainer.data.at(aFiberNum)[t*nx*ny*nz + z*nx*ny + y*nx + x] = aLoadedData[t](x+xoff,y+yoff,z+zoff);
              }
           }
        }
     }
  }
}

void SampleManager::LoadBedpostData(const std::string& aBasename)
{
   std::cout<<"Loading Bedpost samples....."<<std::endl;
   if(aBasename == "")
   {
      std::cout<< "Bad File Name"<<std::endl;
      return;
   }

   //volume4D<float> loadedVolume4D;
   std::string thetaSampleNames;
   std::string phiSampleNames;
   std::string fSampleNames;

   //Single Fiber Case.
   if(NEWIMAGE::fsl_imageexists(aBasename+"_thsamples"))
   {
      thetaSampleNames = aBasename+"_thsamples";
      phiSampleNames = aBasename+"_phisamples";
      fSampleNames = aBasename+"_fsamples";
      LoadBedpostDataHelper(
        thetaSampleNames,phiSampleNames,fSampleNames);
   }
   //Multiple Fiber Case.
   else
   {
      int fiberNum = 1;
      std::string fiberNumAsstring = IntTostring(fiberNum);
      thetaSampleNames = aBasename+"_th"+fiberNumAsstring+"samples";
      bool doesFiberExist = NEWIMAGE::fsl_imageexists(thetaSampleNames);
      while(doesFiberExist)
      {
         phiSampleNames = aBasename+"_ph"+fiberNumAsstring+"samples";
         fSampleNames = aBasename+"_f"+fiberNumAsstring+"samples";

         LoadBedpostDataHelper(
          thetaSampleNames,phiSampleNames,fSampleNames);

         fiberNum++;
         fiberNumAsstring = IntTostring(fiberNum);
         thetaSampleNames = aBasename+"_th"+fiberNumAsstring+"samples";
         doesFiberExist = NEWIMAGE::fsl_imageexists(thetaSampleNames);
      }
      if(fiberNum == 1)
      {
         std::cout<<
          "Could not find samples. Exiting Program..."<<std::endl;
         exit(1);
      }
      std::cout<<"Finished Loading Samples from Bedpost"<<std::endl;
   }
}

void SampleManager::ParseCommandLine(int argc, char** argv)
{
   _oclptxOptions.parse_command_line(argc, argv);

   if (_oclptxOptions.verbose.value()>0)
   {
      _oclptxOptions.status();
   }

   if (_oclptxOptions.simple.value())
   {
      if (_oclptxOptions.matrix1out.value() ||
        _oclptxOptions.matrix3out.value())
      {
         std::cout<<
          "Error: cannot use matrix1 and matrix3 in simple mode"<<
            std::endl;
         exit(1);
      }
      std::cout<<"Running in simple mode"<<std::endl;
      this->LoadBedpostData(_oclptxOptions.basename.value());
      if(_oclptxOptions.seedref.value() == "")
      {
		  NEWIMAGE::read_volume(_seedMask, _oclptxOptions.maskfile.value());
	  }
	  else
	  {
		  NEWIMAGE::read_volume(_seedMask, _oclptxOptions.seedref.value());
	  }
   }
   else if (_oclptxOptions.network.value())
   {
      std::cout<<"Running in network mode"<<std::endl;
   }
   else
   {
      std::cout<<"Running in seedmask mode"<<std::endl;
   }
}

//AFSHIN TODO: Implement Masks

const volume<short int>* SampleManager::GetSeedMask()
{
	return &_seedMask;
}

const unsigned short int* SampleManager::GetSeedMaskToArray()
{
   const int maxZ = _seedMask.maxz();
   const int maxY = _seedMask.maxy();
   const int maxX = _seedMask.maxx();
   const int minZ = _seedMask.minz();
   const int minY = _seedMask.miny();
   const int minX = _seedMask.minx();
   const int sizeX = _seedMask.xsize();
   const int sizeY = _seedMask.ysize();
   const int sizeZ = _seedMask.zsize();
   
   unsigned short int* target = new unsigned short int[sizeX * sizeY * sizeZ];
   for (int z = minZ; z <= maxZ; z++)
   {
     for (int y = minY; y <= maxY; y++)
     {
        for (int x = minX; x <= maxX; x++)
        {
			target[z*sizeX*sizeY + y*sizeX + x] = _seedMask(x,y,z);
        }
     }
   }
   return target;
}

float const SampleManager::GetThetaData(int aFiberNum, int aSamp, int aX, int aY, int aZ)
{
   if(_thetaData.data.size()>0)
   {
      int nx = _thetaData.nx;
      int ny = _thetaData.ny;
      int nz = _thetaData.nz;
      return _thetaData.data.at(aFiberNum)[(aSamp)*(nx*ny*nz) + (aZ)*(nx*ny) + (aY)*nx + (aX)];
   }
   return 0.0f;
}

float const SampleManager::GetPhiData(int aFiberNum, int aSamp, int aX, int aY, int aZ)
{
   if(_phiData.data.size()>0)
   {
      int nx = _phiData.nx;
      int ny = _phiData.ny;
      int nz = _phiData.nz;
      return _phiData.data.at(aFiberNum)[(aSamp)*(nx*ny*nz) + (aZ)*(nx*ny) + (aY)*nx + (aX)];
   }
   return 0.0f;  
}

float const SampleManager::GetfData(int aFiberNum, int aSamp, int aX, int aY, int aZ)
{
   if(_fData.data.size()>0)
   {
      int nx = _fData.nx;
      int ny = _fData.ny;
      int nz = _fData.nz;
      return _fData.data.at(aFiberNum)[(aSamp)*(nx*ny*nz) + (aZ)*(nx*ny) + (aY)*nx + (aX)];
   }
   return 0.0f; 
}

const BedpostXData* SampleManager::GetThetaDataPtr()
{
   if(_thetaData.data.size() > 0)
   {
      return &_thetaData;
   }
   return NULL;
}

const BedpostXData* SampleManager::GetPhiDataPtr()
{
   if(_phiData.data.size() > 0)
   {
      return &_phiData;
   }
   return NULL;
}

const BedpostXData* SampleManager::GetFDataPtr()
{
   if(_fData.data.size() > 0)
   {
      return &_fData;
   }
   return NULL;
}
//*********************************************************************
// samplemanager Constructors/Destructors/Initializers
//*********************************************************************

//Use this to get the only instance of SampleManager in the program.

SampleManager& SampleManager::GetInstance()
{
   if(_manager == NULL)
   {
      _manager = new SampleManager();
   }
   return *_manager;
}SampleManager* SampleManager::_manager;

//Private Constructor.
SampleManager::SampleManager():_oclptxOptions(
  oclptxOptions::getInstance()){}

SampleManager::~SampleManager()
{
   for (unsigned int i = 0; i < _thetaData.data.size(); i++)
   {
      delete[] _thetaData.data.at(i);
   }
   
   for (unsigned int i = 0; i < _phiData.data.size(); i++)
   {
      delete[] _phiData.data.at(i);
   }
   
   for (unsigned int i = 0; i < _fData.data.size(); i++)
   {
      delete[] _fData.data.at(i);
   }
   
   delete _manager;
}


//EOF
