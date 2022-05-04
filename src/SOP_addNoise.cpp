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
 * add noise SOP
 *----------------------------------------------------------------------------
 */

#include "SOP_addNoise.hpp"

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
  std::cout<<"Load add noise "<<std::endl;
  table->addOperator(new OP_Operator("add_noise",
				     "Add Noise FS",
				     SOP_addNoise::myConstructor,
				     SOP_addNoise::myTemplateList,
				     1,
				     1,
				     0));
}

static PRM_Name names[] = {
			   PRM_Name("amp",         "Amplitude"),
			   PRM_Name("dt_",         "Time Step"),
			   PRM_Name("win_size",    "Window Size"),
			   PRM_Name("shift_b",     "Shift B"),
			   PRM_Name("shift_u",     "Shift U"),
			   PRM_Name("ang_speed",   "Angular Speed"),
			   PRM_Name("freq",        "Frequency"),
			   PRM_Name("lacunarity",  "Lacunarity"),
			   PRM_Name("persistence", "Persistence"),
			   PRM_Name("ampli",       "Ampli Perlin"),
			   PRM_Name("octave",      "Octave")
};

PRM_Template
SOP_addNoise::myTemplateList[] = {
				  PRM_Template(PRM_STRING, 1, &PRMgroupName, 0, &SOP_Node::pointGroupMenu,
					       0, 0, SOP_Node::getGroupSelectButton(GA_GROUP_POINT)),
				  PRM_Template(PRM_FLT_J,  1, &names[0], new PRM_Default(0.01), 0,
					       &PRMscaleRange),
				  PRM_Template(PRM_FLT_J,  1, &names[1], PRMoneDefaults, 0,
					       &PRMscaleRange),
				  PRM_Template(PRM_INT_J,  1, &names[2], new PRM_Default(256)),
				  PRM_Template(PRM_INT_J,  1, &names[3], new PRM_Default(124)),
				  PRM_Template(PRM_INT_J,  1, &names[4], new PRM_Default(2)),
				  PRM_Template(PRM_FLT_J,  1, &names[5], PRMoneDefaults, 0,
					       &PRMscaleRange),
				  PRM_Template(PRM_FLT_J,  1, &names[6], new PRM_Default(1), 0,
					       &PRMscaleRange),
				  PRM_Template(PRM_FLT_J,  1, &names[7], new PRM_Default(2), 0,
					       &PRMscaleRange),
				  PRM_Template(PRM_FLT_J,  1, &names[8], new PRM_Default(0.5), 0,
					       &PRMscaleRange),
				  PRM_Template(PRM_FLT_J,  1, &names[9], new PRM_Default(1), 0,
					       &PRMscaleRange),
				  PRM_Template(PRM_INT_J,  1, &names[10], new PRM_Default(1)),
				  PRM_Template(),
};


OP_Node *SOP_addNoise::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
  return new SOP_addNoise(net, name, op);
}

SOP_addNoise::SOP_addNoise(OP_Network *net, const char *name, OP_Operator *op)
  : SOP_Node(net, name, op) {
  mySopFlags.setManagesDataIDs(true);
}

SOP_addNoise::~SOP_addNoise() {
}
OP_ERROR
SOP_addNoise::cookInputGroups(OP_Context &context, int alone) {
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


OP_ERROR SOP_addNoise::cookMySop(OP_Context &context) {
  OP_AutoLockInputs inputs(this);
  if (inputs.lock(context) >= UT_ERROR_ABORT)
    return error();
    
  flags().setTimeDep(1);
  float t = context.getTime();
  int fr = context.getFrame();
  float dt_ = 0.1/3.0;
  t = dt_*fr;
  
  float amp = AMP(t);
  int winsize = WIN_SIZE(t);
  uint shift_b = SHIFT_B(t);
  uint shift_u = SHIFT_U(t);
  FLOAT ang_speed = ANG_SPEED(t);

  duplicateSource(0, context); //grid
  GA_RWHandleV3 Phandle(gdp->findAttribute(GA_ATTRIB_POINT, "P"));

  float sample_rate = 1.0/(dt_);
  float freq_step = 1.0/(winsize*dt_);

  float f_min = freq_step;
  float f_max = sample_rate;
  float k_min = pow((2*M_PI*f_min), 2)/9.81;
  float k_max = pow((2*M_PI*f_max), 2)/9.81;
  float w_min = 2*M_PI/k_max;
  float w_max = 2*M_PI/k_min;
    
  int nb_wl = winsize/2;
  nb_wl -= shift_b + shift_u;
  if (nb_wl <= 0) {
    nb_wl = 1;
    shift_u = 0;
  }
  
  wave_lengths = std::vector<float>(winsize);
  float f_cur = f_min;
  for (int i = 0; i < winsize; ++i) {
    float k_cur = pow((2*M_PI*f_cur), 2)/9.81;
    float w_cur = 2*M_PI/k_cur;
    wave_lengths[/*winsize - 1 - */i] = w_cur;
    f_cur += freq_step;
  }
  
  GA_RWHandleF wl_attrib(gdp->findFloatTuple(GA_ATTRIB_DETAIL, "wavelength", nb_wl));
  if (!wl_attrib.isValid()) {
    wl_attrib = GA_RWHandleF(gdp->addFloatTuple(GA_ATTRIB_DETAIL, "wavelength", nb_wl));
  }
  if (!wl_attrib.isValid()) {
    addError(SOP_MESSAGE, "Failed to create attribute wavelength");
    return error();
  }
  for (int w = 0; w < nb_wl; ++w) {
    wl_attrib.set(0, w, wave_lengths[w+shift_u]);

  }
    
  
  std::vector<float>::iterator itw;
  srand(0);
  SimplexNoise sn(FREQ(t), AMPLI(t), LACUNARITY(t), PERSISTENCE(t));
  {
    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(gdp, ptoff) {
      UT_Vector3 Pvalue = gdp->getPos3(ptoff);
	
      FLOAT perl_r = sn.fractal(OCTAVE(t), Pvalue.x(), Pvalue.z());
      FLOAT perl_i = sn.fractal(OCTAVE(t), Pvalue.z(), Pvalue.x());
      COMPLEX perl(perl_r, perl_i);

      float om = ang_speed;
      Pvalue.y() += amp*real(perl*exp(-COMPLEX(0, 1)*(om*(float)t)));
      gdp->setPos3(ptoff, Pvalue);
    }

  }
  Phandle.bumpDataId();
  return error();
}


