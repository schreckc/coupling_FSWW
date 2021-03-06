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


#ifndef __SOP_WaveParameters_h__
#define __SOP_WaveParameters_h__

#include <SOP/SOP_Node.h>

class SOP_WaveParameters : public SOP_Node {
  
public:
  SOP_WaveParameters(OP_Network *net, const char *name, OP_Operator *op);
  virtual ~SOP_WaveParameters();

  static PRM_Template myTemplateList[];
  static OP_Node *myConstructor(OP_Network*, const char *, OP_Operator *);

  virtual OP_ERROR             cookInputGroups(OP_Context &context, 
					       int alone = 0);

protected:
  virtual OP_ERROR cookMySop(OP_Context &context);
private:
  void        getGroups(UT_String &str){ evalString(str, "group", 0, 0); }
  fpreal      AMP(fpreal t)            { return evalFloat("amp", 0, t); }
  fpreal      WL_MIN()             { return evalFloat("wl_min", 0, 0); }
  fpreal      WL_MAX()             { return evalFloat("wl_max", 0, 0); }
  fpreal      WL_STEP()             { return evalFloat("wl_step", 0, 0); }
  int         BUFFER_SIZE(fpreal t) {return evalInt("buffer_size", 0, t);}
  fpreal      DAMPING(fpreal t) { return evalFloat("damping", 0, t); }
  int        INTER_SRC() {return evalInt("inter_src", 0, 0);}
  int        SPECTRUM_DEF() {return evalInt("spectrum_def", 0, 0);}
  fpreal      DT()             { return evalFloat("dt_", 0, 0); }
  int         WIN_SIZE() {return evalInt("win_size", 0, 0);}
  int         SHIFT_B() {return evalInt("shift_b", 0, 0);}
  int         SHIFT_U() {return evalInt("shift_u", 0, 0);}
  
  const GA_PointGroup *myGroup;
};

#endif
