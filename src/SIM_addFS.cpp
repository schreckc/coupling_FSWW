/*
 * Copyright (c) 2019
 *      Side Effects Software Inc.  All rights reserved.
 *
 * Redistribution and use of Houdini Development Kit samples in source and
 * binary forms, with or without modification, are permitted provided that the
 * following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. The name of Side Effects Software may not be used to endorse or
 *    promote products derived from this software without specific prior
 *    written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SIDE EFFECTS SOFTWARE `AS IS' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN
 * NO EVENT SHALL SIDE EFFECTS SOFTWARE BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *----------------------------------------------------------------------------
 * addFS SIM (!experimental!)
 *----------------------------------------------------------------------------
 * change the velocity field according to a set of sources 
 */

#include "SIM_addFS.hpp"
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <PRM/PRM_Include.h>
#include <SIM/SIM_PRMShared.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_Object.h>
#include <GAS/GAS_SubSolver.h>
#include "definitions.hpp"

using namespace HDK_Sample;

int SIM_addFS::fr = 1;

void
initializeSIM(void *)
{
  IMPLEMENT_DATAFACTORY(SIM_addFS);
}
/// Standard constructor, note that BaseClass was crated by the
/// DECLARE_DATAFACTORY and provides an easy way to chain through
/// the class hierarchy.
SIM_addFS::SIM_addFS(const SIM_DataFactory *factory)
  : BaseClass(factory)
{
}
SIM_addFS::~SIM_addFS()
{
}
/// Used to automatically populate the node which will represent
/// this data type.
const SIM_DopDescription *
SIM_addFS::getDopDescription()
{
  static PRM_Name theGeometryName(GAS_NAME_GEOMETRY, "Geometry Sources");
  static PRM_Name surfaceFieldName(GAS_NAME_SURFACE, "Surface Field");
  static PRM_Name velocityFieldName(GAS_NAME_VELOCITY, "Velocity Field");
  static PRM_Name velxFieldName("velx", "Velocity X");
  static PRM_Name velyFieldName("vely", "Velocity Y");
  static PRM_Name velzFieldName("velz", "Velocity Z");

  static PRM_Template          theTemplates[] = {
						 PRM_Template(PRM_STRING, 1, &theGeometryName),
						 PRM_Template(PRM_STRING, 1, &surfaceFieldName),
						 PRM_Template(PRM_STRING, 1, &velocityFieldName),
						 PRM_Template(PRM_STRING, 1, &velxFieldName),
						 PRM_Template(PRM_STRING, 1, &velyFieldName),
						 PRM_Template(PRM_STRING, 1, &velzFieldName),
						 PRM_Template()
  };
  static SIM_DopDescription    theDopDescription(
						 true,               // Should we make a DOP?
						 "add_fs",   // Internal name of the DOP.
						 "add FS",          // Label of the DOP
						 "Solver",           // Default data name
						 classname(),        // The type of this DOP, usually the class.
						 theTemplates);      // Template list for generating the DOP
  // Make this a microsolver:
  setGasDescription(theDopDescription);
  return &theDopDescription;
}
bool
SIM_addFS::solveGasSubclass(SIM_Engine &engine,
			    SIM_Object *obj,
			    SIM_Time time,
			    SIM_Time timestep)
{
  SIM_VectorField *vel = getVectorField(obj, GAS_NAME_VELOCITY);
  SIM_VectorField *velx = getVectorField(obj, "velx");
  SIM_VectorField *vely = getVectorField(obj, "vely");
  SIM_VectorField *velz = getVectorField(obj, "velz");
  SIM_ScalarField *surface = getScalarField(obj, GAS_NAME_SURFACE);

  if (!velx || !vely || !velz) {
    addError(obj, SIM_MISSINGDATA,
	     "velx vely velz", UT_ERROR_MESSAGE);
    return false;
  }
  if (!vel) {
    addError(obj, SIM_MISSINGDATA,
	     "vel", UT_ERROR_MESSAGE);
    return false;
  }
   
  SIM_GeometryCopy    *geo = 0;
  geo = getGeometryCopy(obj, GAS_NAME_GEOMETRY);
  if (!geo)
    {
      std::cout<<"cannot find geo sources"<<std::endl;
      addError(obj, SIM_MISSINGDATA,
	       "Sources", UT_ERROR_MESSAGE);
      return false;
    }
  // Assume destination is an attribute we need to modify.
  SIM_GeometryAutoWriteLock lock(geo);
  GU_Detail &fs = lock.getGdp();
  UT_WorkBuffer  msg;

  GA_ROHandleF w_handle(fs.findAttribute(GA_ATTRIB_PRIMITIVE, "wavelengths"));
  GA_ROHandleI as_handle(fs.findAttribute(GA_ATTRIB_PRIMITIVE, "ampli_steps"));
  if (!w_handle.isValid()) {
    addError(obj, SIM_MESSAGE, "wavelenghts attribute not defined", UT_ERROR_WARNING);
    return false;
  }
  if (!as_handle.isValid()) {
    addError(obj, SIM_MESSAGE, "ampli_steps attribute not defined", UT_ERROR_WARNING);
    return false;
  }
  GA_ROHandleI bs_handle(fs.findAttribute(GA_ATTRIB_DETAIL, "buffer_size"));
  GA_ROHandleF damping_handle(fs.findAttribute(GA_ATTRIB_DETAIL, "damping"));
  if (!bs_handle.isValid()) {
    addError(obj, SIM_MESSAGE, "buffersize attribute not defined", UT_ERROR_WARNING);
    return false;// error();
  }
  if (!damping_handle.isValid()) {
    addError(obj, SIM_MESSAGE, "damping attribute not defined", UT_ERROR_WARNING);
    return false;// error();
  }
  buffer_size = bs_handle.get(0);
  damping_coef = damping_handle.get(0);
          
  const GA_Attribute *afs = fs.findFloatTuple(GA_ATTRIB_POINT, "ampli", buffer_size);
  if (!afs)	{
    addError(obj, SIM_MESSAGE, "ampli attribute not defined", UT_ERROR_WARNING);
    //      addError(SOP_ATTRIBUTE_INVALID, "ampli");
    return false;// error();
  }
  const GA_AIFTuple *tuple = afs->getAIFTuple();
  GA_Offset ptoff;
  if (!tuple) {
    addError(obj, SIM_MESSAGE, "tuple ampli not defined", UT_ERROR_WARNING);
    return false;
  }

    
    
  GA_Offset prim_off;
  GA_Range range_prim = fs.getPrimitiveRange();
  int w = 0;
  for (GA_Iterator itp = range_prim.begin(); itp != range_prim.end(); ++itp) {
    prim_off = *itp;
    addContribFromPrimitive(fs, prim_off, vel, surface);
  }
  splitVel(vel, velx, vely, velz);
  vel->pubHandleModification();
  velx->pubHandleModification();
  vely->pubHandleModification();
  velz->pubHandleModification();
  ++fr;
  return true;
}

 
void SIM_addFS::addContribFromPrimitivePartial(const GU_Detail &fs,
					       GA_Offset prim_off,
					       SIM_VectorField *vel,
					       SIM_ScalarField *surface,
					       const UT_JobInfo &info) {
  // void SOP_Deform_Surface_inter::addContribFromPrimitive(const GU_Detail *fs,
  // 								GA_Offset prim_off) {
  UT_Interrupt *boss = UTgetInterrupt();
  float dt = 0.1/3.0;
  float t = dt*fr;
  float amp = 1;//AMP(t);

  GA_ROHandleF w_handle(fs.findAttribute(GA_ATTRIB_PRIMITIVE, "wavelengths"));
  GA_ROHandleI as_handle(fs.findAttribute(GA_ATTRIB_PRIMITIVE, "ampli_steps"));
  const GA_Attribute *afs = fs.findFloatTuple(GA_ATTRIB_POINT, "ampli", buffer_size);
  const GA_AIFTuple *tuple;
  if (!w_handle.isValid()) {
    std::cerr<< "wavelengths attribute not defined"<<std::endl;
    return;
    //  return error();
  }
  if (!as_handle.isValid()) {
    std::cerr<< "ampli_steps attribute not defined"<<std::endl;
    return;
  }
  float wl = w_handle.get(prim_off);
  int as = as_handle.get(prim_off);
  float k = M_PI*2.0/wl;
  float om = omega(k);
  const GA_Primitive* prim = fs.getPrimitive(prim_off);
	
  GA_Range range = prim->getPointRange();
  GA_Offset ptoff;
  int i= 0, n;
  UT_VoxelArrayIteratorF vitx;
  UT_VoxelArrayIteratorF vity;
  UT_VoxelArrayIteratorF vitz;
  vitx.setArray(vel->getField(0)->fieldNC());
  vitx.setCompressOnExit(true);
  vitx.setPartialRange(info.job(), info.numJobs());
  UT_Vector3 pos;
  if (vel->getField(0)->isAligned(vel->getField(1)) && vel->getField(0)->isAligned(vel->getField(2)) ) {

    for (vitx.rewind(); !vitx.atEnd(); vitx.advance()) {
      if (vitx.isStartOfTile() && boss->opInterrupt()) {
        break;
      }
      vel->indexToPos(0,vitx.x(), vitx.y(), vitx.z(), pos);
      float d = surface->getField()->getCellValue(vitx.x(), vitx.y(), vitx.z());
      if (d > 0) {
	d = 0;
      } else {
	d = 1;
      }
      FLOAT hcoef = exp(k*pos.y())*d*om;
      COMPLEX ax = 0, ay = 0, az = 0;//exp(std::complex<float>(0, 1)*(k*Pvalue.x()));
      for(GA_Iterator it = range.begin(); it != range.end(); ++it) {
	UT_Vector3 P_fs = fs.getPos3((*it));
	float rx = pos.x() - P_fs.x();
	float rz = pos.z() - P_fs.z();
	float r = sqrt(pow(rx, 2) + pow(rz, 2));
	float kx  =rx/r;
	float kz  =rz/r;
	if (r == 0) {
	  kx = 0;
	  kz = 0;
	}
	float v = velocity(k, om);
	float ar = 0, ai = 0;
	    
	fpreal t_ret = t - r/v;
	int f_ret = floor(t_ret/(dt)) - 1;
	f_ret = floor((float)f_ret/(float)as); 
	if (t_ret < 0) {
	  --f_ret;
	}
	if (f_ret >= 0) {
	  afs = fs.findFloatTuple(GA_ATTRIB_POINT, "ampli", buffer_size);
	  tuple = afs->getAIFTuple();
	  tuple->get(afs, *it, ar, 2*f_ret);
	  tuple->get(afs, *it, ai, 2*f_ret+1);
	}
	std::complex<float> ampli(ar, ai);
	ax += hcoef*kx*ampli*fund_solution(k*r)*damping(damping_coef, r, k);
	az += hcoef*kz*ampli*fund_solution(k*r)*damping(damping_coef, r, k);
	ay += hcoef*ampli*fund_solution(k*r)*damping(damping_coef, r, k)*COMPLEX(0, 1);
      
      }
      vitx.setValue(real(amp*ax*exp(-COMPLEX(0, 1)*(om*(float)t))));
      vel->getField(1)->setCellValue(vitx.x(), vitx.y(), vitx.z(),
				     real(amp*ay*exp(-COMPLEX(0, 1)*(om*(float)t))));
      vel->getField(2)->setCellValue(vitx.x(), vitx.y(), vitx.z(),
				     real(amp*az*exp(-COMPLEX(0, 1)*(om*(float)t))));
    }
       
  }
}


void SIM_addFS::splitVelPartial(SIM_VectorField *vel,
				SIM_VectorField *velx,
				SIM_VectorField *vely,
				SIM_VectorField *velz,
				const UT_JobInfo &info) {
  UT_Interrupt *boss = UTgetInterrupt();
  UT_VoxelArrayIteratorF vitx;
  vitx.setArray(vel->getField(0)->fieldNC());
  vitx.setCompressOnExit(true);
  vitx.setPartialRange(info.job(), info.numJobs());
  UT_Vector3 pos;
  for (vitx.rewind(); !vitx.atEnd(); vitx.advance()) {
    if (vitx.isStartOfTile() && boss->opInterrupt()) {
      break;
    }
    velx->getField(0)->setCellValue(vitx.x(), vitx.y(), vitx.z(),
				    vel->getField(0)->getCellValue(vitx.x(), vitx.y(), vitx.z()));
    vely->getField(1)->setCellValue(vitx.x(), vitx.y(), vitx.z(),
				    vel->getField(1)->getCellValue(vitx.x(), vitx.y(), vitx.z()));
    velz->getField(2)->setCellValue(vitx.x(), vitx.y(), vitx.z(),
				    vel->getField(2)->getCellValue(vitx.x(), vitx.y(), vitx.z()));
  }
}
   
