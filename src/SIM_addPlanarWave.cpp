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
 * addPlanarWave SIM (!experimental!)
 *----------------------------------------------------------------------------
 * change the velocity field according to a planar wave
 */

#include "SIM_addPlanarWave.hpp"
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
///
/// This is the hook that Houdini grabs from the dll to link in
/// this.  As such, it merely has to implement the data factory
/// for this node.
///

int SIM_addPlanarWave::fr = 1;

void
initializeSIM(void *)
{
  IMPLEMENT_DATAFACTORY(SIM_addPlanarWave);
}
/// Standard constructor, note that BaseClass was crated by the
/// DECLARE_DATAFACTORY and provides an easy way to chain through
/// the class hierarchy.
SIM_addPlanarWave::SIM_addPlanarWave(const SIM_DataFactory *factory)
  : BaseClass(factory)
{
}
SIM_addPlanarWave::~SIM_addPlanarWave()
{
}

static PRM_Name Frequency("freq", "Frequency");
static PRM_Name Amplitude("ampli", "Amplitude");
static PRM_Name Wavy("wavy", "Wavy");
static PRM_Name Clamp("clamp", "Clamp");
static PRM_Name Aperiodic("aperiodic", "Aperiodic");


/// Used to automatically populate the node which will represent
/// this data type.
const SIM_DopDescription *
SIM_addPlanarWave::getDopDescription()
{
  static PRM_Name surfaceFieldName(GAS_NAME_SURFACE, "Surface Field");
  static PRM_Name velocityFieldName(GAS_NAME_VELOCITY, "Velocity Field");
  static PRM_Name velxFieldName("velx", "Velocity X");
  static PRM_Name velyFieldName("vely", "Velocity Y");
  static PRM_Name velzFieldName("velz", "Velocity Z");

  static PRM_Template          theTemplates[] = {
						 PRM_Template(PRM_STRING, 1, &surfaceFieldName),
						 PRM_Template(PRM_STRING, 1, &velocityFieldName),
						 PRM_Template(PRM_STRING, 1, &velxFieldName),
						 PRM_Template(PRM_STRING, 1, &velyFieldName),
						 PRM_Template(PRM_STRING, 1, &velzFieldName),
						 PRM_Template(PRM_FLT, 1, &Frequency),
						 PRM_Template(PRM_FLT, 1, &Amplitude),
						 PRM_Template(PRM_TOGGLE, 1, &Wavy),
						 PRM_Template(PRM_TOGGLE, 1, &Clamp),
						 PRM_Template(PRM_TOGGLE, 1, &Aperiodic),
						 PRM_Template()
  };
  static SIM_DopDescription    theDopDescription(
						 true,               // Should we make a DOP?
						 "add_planar",   // Internal name of the DOP.
						 "add Planar Wave",          // Label of the DOP
						 "Solver",           // Default data name
						 classname(),        // The type of this DOP, usually the class.
						 theTemplates);      // Template list for generating the DOP
  // Make this a microsolver:
  setGasDescription(theDopDescription);
  return &theDopDescription;
}
bool
SIM_addPlanarWave::solveGasSubclass(SIM_Engine &engine,
				    SIM_Object *obj,
				    SIM_Time time,
				    SIM_Time timestep)
{
  SIM_VectorField *vel = getVectorField(obj, GAS_NAME_VELOCITY);
  SIM_VectorField *velx = getVectorField(obj, "velx");
  SIM_VectorField *vely = getVectorField(obj, "vely");
  SIM_VectorField *velz = getVectorField(obj, "velz");
  SIM_ScalarField *surface = getScalarField(obj, GAS_NAME_SURFACE);
  FLOAT freq = getFreq();
  FLOAT ampli = getAmpli();

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
  if (!surface) {
    addError(obj, SIM_MISSINGDATA,
	     "surface", UT_ERROR_MESSAGE);
    return false;
  } 
  addWave(freq, ampli, vel, surface);
  splitVel(vel, velx, vely, velz);
  vel->pubHandleModification();
  velx->pubHandleModification();
  vely->pubHandleModification();
  velz->pubHandleModification();
  ++fr;
  return true;
}

 
void SIM_addPlanarWave::addWavePartial(float freq,
				       float ampli,
				       SIM_VectorField *vel,
				       SIM_ScalarField *surface,
				       const UT_JobInfo &info) {

  UT_Interrupt *boss = UTgetInterrupt();
  float dt = 0.1/3.0;
  float t = dt*fr;
  float om = 2*M_PI*freq/COEF_DISPERSION;
  float k = om*om/9.81;
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
      float d;
      float prof =  pos.y() - real(ampli*exp(COMPLEX(0, 1)*(om*(float)t - (float)M_PI/2.0f)));
      if (prof > 0) {
     	d = 0;
      } else {
     	d = 1;
      }
      float rx = pos.x();
      if (rx < 0 && getClamp()) {
	d = 0;
      }
       
      COMPLEX hcoef = d*om*exp(k*prof);
      if (getWavy()) {
	hcoef *= exp(COMPLEX(0, 1)*k*rx);
      }
      COMPLEX ax = 0, ay = 0, az = 0;
         
      ax += hcoef;
      az += 0;
      ay += hcoef;

      if (getAperiodic() && -rx > velocity(k, om)*t) {
	ax = 0;
	ay = 0;
	az = 0;
      }
      
      vel->getField(0)->setCellValue(vitx.x(), vitx.y(), vitx.z(),
				     real(ampli*ax*exp(COMPLEX(0, 1)*(om*(float)t - (float)M_PI/2.0f))));
      vel->getField(1)->setCellValue(vitx.x(), vitx.y(), vitx.z(),
       				     real(ampli*ay*exp(COMPLEX(0, 1)*(om*(float)t - (float)M_PI/2.0f))));
      vel->getField(2)->setCellValue(vitx.x(), vitx.y(), vitx.z(),
       				     imag(ampli*az*exp(COMPLEX(0, 1)*(om*(float)t - (float)M_PI/2.0f))));
    }
  }
}



void SIM_addPlanarWave::splitVelPartial(SIM_VectorField *vel,
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
    vely->getField(0)->setCellValue(vitx.x(), vitx.y(), vitx.z(),
				    vel->getField(1)->getCellValue(vitx.x(), vitx.y(), vitx.z()));
    velz->getField(0)->setCellValue(vitx.x(), vitx.y(), vitx.z(),
				    vel->getField(2)->getCellValue(vitx.x(), vitx.y(), vitx.z()));
  }
}
   
