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


#ifndef __SOP_PlanarWave_h__
#define __SOP_PlanarWave_h__

#include <SOP/SOP_Node.h>

class SOP_PlanarWave : public SOP_Node {
public:
  SOP_PlanarWave(OP_Network *net, const char *name, OP_Operator *op);
  virtual ~SOP_PlanarWave();

  static PRM_Template myTemplateList[];
  static OP_Node *myConstructor(OP_Network*, const char *, OP_Operator *);

  virtual OP_ERROR             cookInputGroups(OP_Context &context, 
					       int alone = 0);

protected:
  virtual OP_ERROR cookMySop(OP_Context &context);
private:
  void        getGroups(UT_String &str){ evalString(str, "group", 0, 0); }
  fpreal      AMP(fpreal t)           { return evalFloat("amp", 0, t); }
  fpreal      FREQ(fpreal t)           { return evalFloat("freq", 0, t); }
  fpreal      WL(fpreal t)           { return evalFloat("wl", 0, t); }
  fpreal      USE_WL(fpreal t)           { return evalFloat("use_wl", 0, t); }
  fpreal      WAVY(fpreal t)           { return evalFloat("wavy", 0, t); }
  fpreal      CLAMP(fpreal t)           { return evalFloat("clamp", 0, t); }
  fpreal      APERIODIC(fpreal t)           { return evalFloat("aperiodic", 0, t); }
  fpreal      SIGMA(fpreal t)           { return evalFloat("sigma", 0, t); }
  fpreal      X(fpreal t) { return evalFloat("dir", 0, t); }
  fpreal      Y(fpreal t) { return evalFloat("dir", 1, t); }
  fpreal      Z(fpreal t) { return evalFloat("dir", 2, t); }
  int        T() { return evalInt("use_time", 0, 0); }
  fpreal      PHASE(fpreal t) { return evalFloat("phase", 0, t); }
  const GA_PointGroup *myGroup;
};


#endif
