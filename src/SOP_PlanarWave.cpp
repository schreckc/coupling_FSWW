/*
 * MIT License
 * 
 * Copyright (c) 2022 Camille Schreck
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 *----------------------------------------------------------------------------
 * PlanarWave SOP
 *----------------------------------------------------------------------------
 * At first frame, compute the spectrum a each point of the grid (input 0) by 
 * summing the contribution of all sources (each other input is a set of sources).
 * At each frame, sum the contribution of each wavelength at each point to get 
 * the height (according to the spectrum of the point).
 */


#include "SOP_PlanarWave.hpp"

#include <GU/GU_Detail.h>
#include <GA/GA_Handle.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_Director.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <UT/UT_DSOVersion.h>
#include <SYS/SYS_Math.h>

#include "definitions.hpp"

void newSopOperator(OP_OperatorTable *table) {
  table->addOperator(new OP_Operator("planar_wave",
				     "Planar Wave",
				     SOP_PlanarWave::myConstructor,
				     SOP_PlanarWave::myTemplateList,
				     1,
				     1,
				     0));
}

static PRM_Default dirDefaults[] = {PRM_Default(1),
				    PRM_Default(0),
				    PRM_Default(0)};

static PRM_Default phaseDefault((float)M_PI/2.0f);

static PRM_Name names[] = {
			   PRM_Name("amp",    "Amplitude"),
			   PRM_Name("freq",   "Frequency"),
			   PRM_Name("wl",   "WaveLenght"),
			   PRM_Name("use_wl",   "Use Wavelength"),
			   PRM_Name("wavy",   "Wavy"),
			   PRM_Name("clamp",   "Clamp"),
			   PRM_Name("aperiodic",   "Aperiodic"),
			   PRM_Name("sigma",   "Sigma"),
			   PRM_Name("dir",   "Direction"),
			   PRM_Name("use_time",   "Use time"),
			   PRM_Name("phase",   "Phase"),
};

PRM_Template SOP_PlanarWave::myTemplateList[] = {
						 PRM_Template(PRM_STRING,    1, &PRMgroupName, 0, &SOP_Node::pointGroupMenu,
							      0, 0, SOP_Node::getGroupSelectButton(GA_GROUP_POINT)),
						 PRM_Template(PRM_FLT_J,     1, &names[0], PRMoneDefaults, 0,
							      &PRMscaleRange),
						 PRM_Template(PRM_FLT_J,     1, &names[1], PRMoneDefaults, 0,
							      &PRMscaleRange),
						 PRM_Template(PRM_FLT_J,     1, &names[2], PRMoneDefaults, 0,
							      &PRMscaleRange),
						 PRM_Template(PRM_TOGGLE_J,  1, &names[3]),
						 PRM_Template(PRM_TOGGLE_J,  1, &names[4]),
						 PRM_Template(PRM_TOGGLE_J,  1, &names[5]),
						 PRM_Template(PRM_TOGGLE_J,  1, &names[6]),
						 PRM_Template(PRM_FLT_J,     1, &names[7], PRMoneDefaults, 0,
							      &PRMscaleRange),
						 PRM_Template(PRM_XYZ_J,     3, &names[8], dirDefaults),
						 PRM_Template(PRM_TOGGLE_J,     1, &names[9], PRMzeroDefaults),
						 PRM_Template(PRM_FLT_J,     1, &names[10], &phaseDefault, 0,
							      &PRMscaleRange),
						 PRM_Template(),
};


OP_Node *SOP_PlanarWave::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
  return new SOP_PlanarWave(net, name, op);
}

SOP_PlanarWave::SOP_PlanarWave(OP_Network *net, const char *name, OP_Operator *op)
  : SOP_Node(net, name, op) {
  mySopFlags.setManagesDataIDs(true);
}

SOP_PlanarWave::~SOP_PlanarWave()
{
}
OP_ERROR SOP_PlanarWave::cookInputGroups(OP_Context &context, int alone) {
  return cookInputPointGroups(context,
			      myGroup,
			      alone,
			      true,
			      0,
			      -1,
			      true,
			      false,
			      true,
			      0);
}


OP_ERROR SOP_PlanarWave::cookMySop(OP_Context &context) {
  OP_AutoLockInputs inputs(this);
  if (inputs.lock(context) >= UT_ERROR_ABORT)
    return error();
    
  flags().setTimeDep(1);
  fpreal frame = OPgetDirector()->getChannelManager()->getSample(context.getTime());
  frame *= 0.0333;
  float t = context.getTime();
  int fr = context.getFrame();
  float dt = 0.1/3.0;
  t = dt*fr;
  if (!T()) {
    t = 0;
  }

  float amp = AMP(t);
  float freq = FREQ(t);
  float sigma = SIGMA(t);
  float wl = WL(t);
  float ph = PHASE(t);
  UT_Vector3 dir(X(t), Y(t), Z(t));
  dir.normalize();
  
  duplicateSource(0, context); //grid

  GA_RWHandleV3 Phandle(gdp->findAttribute(GA_ATTRIB_POINT, "P"));

  float k, om;
  if (USE_WL(t)) {
    k = 2*M_PI/wl;
    om = COEF_DISPERSION*sqrt(k*9.81);
  } else {
    om = 2*M_PI*freq;
    k = pow(om/COEF_DISPERSION, 2)/9.81;
  }

  GA_Offset ptoff;
  GA_FOR_ALL_PTOFF(gdp, ptoff) {
    UT_Vector3 Pvalue = gdp->getPos3(ptoff);
    float rx = dot(Pvalue, dir);
    if (rx > 0 || !CLAMP(t)) {
      COMPLEX a = amp;
      if (APERIODIC(t) && -rx > velocity(k, om)*t) {
	a = 0;
      }
      if (WAVY(t)) {
	a *= exp(COMPLEX(0, 1)*(k*rx+ph));
      }
      Pvalue.y() = exp(rx*sigma)*real(a*exp(COMPLEX(0, 1)*(om*(float)t)));
      gdp->setPos3(ptoff, Pvalue);
    }
  }
  Phandle.bumpDataId();
  return error();
}
