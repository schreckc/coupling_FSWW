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
 * WaveParameters SOP
 *----------------------------------------------------------------------------
 * Define the parameters of the sources
 * Parameters:
 *    -Amplitude: multiply the amplitude of every sources by this number
 *    If defined by spectrum:
 *        - use the wavelengths obtained by doing an FFT on a time window of size <window size> (time step 0.3) using the dispertion relationship defined in definition.hpp. The parameters shiftB and shiftU remove respectively the <shifb> smaller wavelength and <shiftu> bigger ones.
 *   If not defined by spectrum  
 *        -Minimum/maximum wavelength, wavelength multiplicative step: used to compute the range of wl as a geometric sequence
 *    -Interactive sources
 *    -Size of the buffer: size of the buffer recording past amplitude
 *    -Damping
 */


#include "SOP_WaveParameters.hpp"

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
#include <vector>

void newSopOperator(OP_OperatorTable *table) {
  table->addOperator(new OP_Operator("wave_parameters_fs",
				     "Wave Parameters FS",
				     SOP_WaveParameters::myConstructor,
				     SOP_WaveParameters::myTemplateList,
				     0,
				     0,
				     nullptr,  
				     OP_FLAG_GENERATOR));
}
static PRM_Name names[] = {
			   PRM_Name("amp",     "Amplitude"),

			   PRM_Name("spectrum_def", "Defined by spectrum"),

			   PRM_Name("wl_min",  "Minimum Wavelength"),
			   PRM_Name("wl_max",  "Maximum Wavelength"),
			   PRM_Name("wl_step",  "Wavelength Multiplicative Step"),

			   PRM_Name("win_size",   "Window Size"),
			   PRM_Name("shift_b",   "Shift B"),
			   PRM_Name("shift_u",   "Shift U"),

			   PRM_Name("inter_src",   "Interactive sources"),
			   PRM_Name("buffer_size",   "Size of the buffer containing past amplitudes"),

			   PRM_Name("damping",   "Damping"),
			   PRM_Name("dt_",     "Time Step"),
};

PRM_Template
SOP_WaveParameters::myTemplateList[] = {
					PRM_Template(PRM_STRING,    1, &PRMgroupName, 0, &SOP_Node::pointGroupMenu,
						     0, 0, SOP_Node::getGroupSelectButton(GA_GROUP_POINT)),
					PRM_Template(PRM_FLT_J,     1, &names[0], PRMoneDefaults, 0, //ampli
						     &PRMscaleRange),
  
					PRM_Template(PRM_TOGGLE_E,  1, &names[1]), //spectrum def
  
					PRM_Template(PRM_FLT_E,     1, &names[2], PRMoneDefaults), //min wl
					PRM_Template(PRM_FLT_E,     1, &names[3], PRMoneDefaults), //max wl
					PRM_Template(PRM_FLT_E,     1, &names[4], PRMoneDefaults), //step wl

					PRM_Template(PRM_INT_E,  1, &names[5], new PRM_Default(256)), //win size
					PRM_Template(PRM_INT_E,  1, &names[6], new PRM_Default(124)), //shift b
					PRM_Template(PRM_INT_E,  1, &names[7], new PRM_Default(2)), //shift u

					PRM_Template(PRM_TOGGLE_E,  1, &names[8]), //interactive src
					PRM_Template(PRM_INT_J,     1, &names[9], new PRM_Default(500)), //buffer size
  
					PRM_Template(PRM_FLT_J,     1, &names[10], PRMzeroDefaults), //damping
					PRM_Template(PRM_FLT_E,     1, &names[11], new PRM_Default(0.03)), //time step

					PRM_Template(),
};


OP_Node *
SOP_WaveParameters::myConstructor(OP_Network *net, const char *name, OP_Operator *op) {
  return new SOP_WaveParameters(net, name, op);
}

SOP_WaveParameters::SOP_WaveParameters(OP_Network *net, const char *name, OP_Operator *op)
  : SOP_Node(net, name, op) {
  mySopFlags.setManagesDataIDs(true);
}

SOP_WaveParameters::~SOP_WaveParameters()
{
}
OP_ERROR
SOP_WaveParameters::cookInputGroups(OP_Context &context, int alone)
{
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


OP_ERROR SOP_WaveParameters::cookMySop(OP_Context &context) {
  OP_AutoLockInputs inputs(this);
  if (inputs.lock(context) >= UT_ERROR_ABORT)
    return error();
  flags().setTimeDep(0);

  float t = context.getTime();
  int fr = context.getFrame();

  float dt_ = DT();

  bool is_inter = INTER_SRC();
  uint buffer_size = 2;
  if (is_inter) {
    buffer_size = BUFFER_SIZE(0);
  }
  std::vector<float> wave_lengths;
  std::vector<int> ampli_steps;
  int nb_wl = 1;
  if (!SPECTRUM_DEF()) {
    // Compute the range of wavelength we want to use between WL_MIN and WL_MAX(t)
    // Note: multiplicative step between wl
    FLOAT wl = WL_MIN();
    wave_lengths.push_back(wl);

    float max_wl = WL_MAX();
    float step_wl = 1+WL_STEP();
    while (wl < max_wl) {
      wl *= step_wl;
      wave_lengths.push_back(wl);
      ++nb_wl;
    }
  
    // compute the number of time step between updates for each wavelength
    ampli_steps = std::vector<int>(nb_wl);
    wl = wave_lengths[0];
    FLOAT period = 0.5*wl/velocity(2*M_PI/wl); 
    int d_period = period/(dt_);
      
    if (d_period == 0) {
      ampli_steps[0] = 1;
    } else {
      ampli_steps[0] = d_period;
    }
    for (int w = 1; w < nb_wl; ++w) {
      wl = wave_lengths[w];
      period = 0.5*wl/velocity(2*M_PI/wl); 
      d_period = period/(dt_*ampli_steps[0]);
      if (d_period == 0) {
	ampli_steps[w] = 1;
      } else {
	ampli_steps[w] = d_period*ampli_steps[0];
      }
    }
  } else {
    /** compute a range a wavelength corresponding to the range of 
	(time)-frequencies obtained by running and FFT on a recorded 
	wave with a window of size winsize and a time step dt_ **/

    /** many of the smaller frequencies are irrelevant (depending 
	on the resolution of the grid and the resolution of the 
	simulation. The parameters shift_b and shift_u are used 
	to ignore too big or too small frequencies **/    
    int winsize = WIN_SIZE();
    uint shift_b = SHIFT_B();
    uint shift_u = SHIFT_U();
    float sample_rate = 1.0/(dt_);
    float freq_step = 1.0/(winsize*dt_);

    float f_min = freq_step;
    float f_max = sample_rate;
    float k_min = pow((2*M_PI*f_min/COEF_DISPERSION), 2)/9.81;
    float k_max = pow((2*M_PI*f_max/COEF_DISPERSION), 2)/9.81;
    float w_min = 2*M_PI/k_max;
    float w_max = 2*M_PI/k_min;
    // std::cout<<"wl min-max "<<w_min<<" - "<<w_max<<std::endl;
    // std::cout<<"freq min-max "<<f_min<<" - "<<f_max<<std::endl;
    
    nb_wl = winsize/2;
    nb_wl -= shift_b + shift_u;
    if (nb_wl <= 0) {
      nb_wl = 1;
      shift_u = 0;
    }
    // Compute the range of wavelength we want to use between WL_MIN and WL_MAX(t)
    // Note: multiplicative step between wl
    wave_lengths = std::vector<float>(winsize);
    float f_cur = f_min;
    wave_lengths[0] = 0; // TEST get wl 0
    for (int i = 1; i < winsize; ++i) { // TEST get wl 0
      float k_cur = pow((2*M_PI*f_cur/COEF_DISPERSION), 2)/9.81;
      float w_cur = 2*M_PI/k_cur;
      wave_lengths[/*winsize - 1 - */i] = w_cur;
      f_cur += freq_step;
      // std::cout<<i<<" "<<wave_lengths[i]<<std::endl;
    }
  
    // compute the number of time step between updates for each wavelength
    ampli_steps = std::vector<int>(winsize);
    float wl = wave_lengths[winsize - nb_wl - 1];
    FLOAT period = 0.125*wl/velocity(2*M_PI/wl); 
    int d_period = period/(dt_);

    std::cout<<"nb wl "<<nb_wl<<" " <<winsize - nb_wl - 1<<std::endl;
    if (d_period == 0) {
      ampli_steps[winsize - nb_wl - 1] = 1;
    } else {
      ampli_steps[winsize - nb_wl - 1] = d_period;
    }
    for (int w = winsize-2; w >= 0; --w) {
      wl = wave_lengths[w];
      if (wl == 0) { // TEST get wl 0
	wl = wave_lengths[w+1];
      }
      period = 0.125*wl/velocity(2*M_PI/wl); 
      d_period = period/(dt_*ampli_steps[winsize - nb_wl - 1]);
      if (d_period == 0) {
	ampli_steps[w] = 1;
      } else {
	ampli_steps[w] = d_period*ampli_steps[winsize - nb_wl - 1];
      }
    }

  }
  
  // begin creation of geometry
  gdp->clearAndDestroy();
  // creation of the primitive (one per wl), and link one source to each
  for (int w = 0; w < nb_wl; ++w) {
    gdp->appendPrimitive(GA_PRIMPOLY);
  }

  // creation of the primitve attibutes
  GA_RWHandleF wl_attrib(gdp->findFloatTuple(GA_ATTRIB_PRIMITIVE, "wavelengths", 1));
  GA_RWHandleI as_attrib(gdp->findIntTuple(GA_ATTRIB_PRIMITIVE, "ampli_steps", 1));
  if (!wl_attrib.isValid()) {
    wl_attrib = GA_RWHandleF(gdp->addFloatTuple(GA_ATTRIB_PRIMITIVE, "wavelengths", 1));
  }
  if (!wl_attrib.isValid()) {
    addError(SOP_MESSAGE, "Failed to create attribute wavelengths");
    return error();
  }
  if (!as_attrib.isValid()) {
    as_attrib = GA_RWHandleI(gdp->addIntTuple(GA_ATTRIB_PRIMITIVE, "ampli_steps", 1));
  }
  if (!as_attrib.isValid()) {
    addError(SOP_MESSAGE, "Failed to create attribute ampli_steps");
    return error();
  }
  // setting primitive attibutes
  GA_Offset prim_off;
  GA_Offset lcl_start, lcl_end;
  int w = 0;
  if (SPECTRUM_DEF()) {
    for (GA_Iterator lcl_it((gdp)->getPrimitiveRange()); lcl_it.blockAdvance(lcl_start, lcl_end); ) {
      for (prim_off = lcl_start; prim_off < lcl_end; ++prim_off) {
	wl_attrib.set(prim_off, wave_lengths[w+SHIFT_U()]);
	as_attrib.set(prim_off, ampli_steps[w+SHIFT_U()]);
	++w;
      }
    }
  } else {
    for (GA_Iterator lcl_it((gdp)->getPrimitiveRange()); lcl_it.blockAdvance(lcl_start, lcl_end); ) {
      for (prim_off = lcl_start; prim_off < lcl_end; ++prim_off) {
	wl_attrib.set(prim_off, wave_lengths[w]);
	as_attrib.set(prim_off, ampli_steps[w]);
	++w;
      }
    }
  }

  // creation of the detail attibutes
  GA_RWHandleI bs_attrib(gdp->findIntTuple(GA_ATTRIB_DETAIL, "buffer_size", 1));
  GA_RWHandleF damping_attrib(gdp->findFloatTuple(GA_ATTRIB_DETAIL, "damping", 1));
  GA_RWHandleI inter_attrib(gdp->findIntTuple(GA_ATTRIB_DETAIL, "inter", 1));
  GA_RWHandleF ampli_attrib(gdp->findFloatTuple(GA_ATTRIB_DETAIL, "ampli", 1));
  GA_RWHandleF dt_attrib(gdp->findFloatTuple(GA_ATTRIB_DETAIL, "time_step", 1));
  GA_RWHandleI w_attrib(gdp->findIntTuple(GA_ATTRIB_DETAIL, "winsize", 1));
  GA_RWHandleI sb_attrib(gdp->findIntTuple(GA_ATTRIB_DETAIL, "shift_b", 1));
  GA_RWHandleI su_attrib(gdp->findIntTuple(GA_ATTRIB_DETAIL, "shift_u", 1));
  if (!bs_attrib.isValid()) {
    bs_attrib = GA_RWHandleI(gdp->addIntTuple(GA_ATTRIB_DETAIL, "buffer_size", 1));
  }
  if (!bs_attrib.isValid()) {
    addError(SOP_MESSAGE, "Failed to create attribute buffer_size");
    return error();
  }
  if (!damping_attrib.isValid()) {
    damping_attrib = GA_RWHandleF(gdp->addFloatTuple(GA_ATTRIB_DETAIL, "damping", 1));
  }
  if (!damping_attrib.isValid()) {
    addError(SOP_MESSAGE, "Failed to create attribute damping");
    return error();
  }
  if (!inter_attrib.isValid()) {
    inter_attrib = GA_RWHandleI(gdp->addIntTuple(GA_ATTRIB_DETAIL, "inter", 1));
  }
  if (!inter_attrib.isValid()) {
    addError(SOP_MESSAGE, "Failed to create attribute inter");
    return error();
  }
  if (!ampli_attrib.isValid()) {
    ampli_attrib = GA_RWHandleF(gdp->addFloatTuple(GA_ATTRIB_DETAIL, "ampli", 1));
  }
  if (!ampli_attrib.isValid()) {
    addError(SOP_MESSAGE, "Failed to create attribute ampli");
    return error();
  }
  if (!dt_attrib.isValid()) {
    dt_attrib = GA_RWHandleF(gdp->addFloatTuple(GA_ATTRIB_DETAIL, "time_step", 1));
  }
  if (!dt_attrib.isValid()) {
    addError(SOP_MESSAGE, "Failed to create attribute time_step");
    return error();
  }
  if (!w_attrib.isValid()) {
    w_attrib = GA_RWHandleI(gdp->addIntTuple(GA_ATTRIB_DETAIL, "winsize", 1));
  }
  if (!sb_attrib.isValid()) {
    sb_attrib = GA_RWHandleI(gdp->addIntTuple(GA_ATTRIB_DETAIL, "shift_b", 1));
  }
  if (!su_attrib.isValid()) {
    su_attrib = GA_RWHandleI(gdp->addIntTuple(GA_ATTRIB_DETAIL, "shift_u", 1));
  }
   
  // setting detail attibutes
  bs_attrib.set(0,buffer_size);
  damping_attrib.set(0,DAMPING(t));
  inter_attrib.set(0,is_inter);
  ampli_attrib.set(0,AMP(t));
  dt_attrib.set(0,dt_);
  w_attrib.set(0,WIN_SIZE());
  su_attrib.set(0,SHIFT_U());
  sb_attrib.set(0,SHIFT_B());
  
 
  wl_attrib->bumpDataId();
  as_attrib->bumpDataId();
  bs_attrib->bumpDataId();
  damping_attrib->bumpDataId();
  inter_attrib->bumpDataId();
  ampli_attrib->bumpDataId();
  dt_attrib->bumpDataId();
  w_attrib->bumpDataId();
  su_attrib->bumpDataId();
  sb_attrib->bumpDataId();
  return error();
}
