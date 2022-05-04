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
 *
 *----------------------------------------------------------------------------
 * Circle Obstalce SOP
 *---------------------------------------------------------------------------
 * SOP_CircleObstacle_Src (Square/Texture -> need update):
 * Create set of point sources along an offset surface of the obstacle.
 * Need the waveParameters as input.
 * Create a subset of sources for each wl.
 * Spacing depends on the wavelength and the parameter "density".
 * Note: the number of boundary points should be at least twice the number of points for the biggest subset of sources (the one corresponding to the smallest wavelength).
 */


#ifndef __SOP_Circle_Obstacle_Src_h__
#define __SOP_Circle_Obstacle_Src_h__

#include <SOP/SOP_Node.h>
#include "definitions.hpp"
#include <Eigen/SVD>

class SOP_Circle_Obstacle_Src : public SOP_Node {
  
public:
  SOP_Circle_Obstacle_Src(OP_Network *net, const char *name, OP_Operator *op);
  virtual ~SOP_Circle_Obstacle_Src();

  static PRM_Template myTemplateList[];
  static OP_Node *myConstructor(OP_Network*, const char *, OP_Operator *);

  virtual OP_ERROR             cookInputGroups(OP_Context &context, 
					       int alone = 0);

protected:
  virtual OP_ERROR cookMySop(OP_Context &context);
private:
  void        getGroups(UT_String &str){ evalString(str, "group", 0, 0); }
  fpreal      OFF(fpreal t)            { return evalFloat("off", 0, t); }
  fpreal      DENSITY(fpreal t)          { return evalFloat("density", 0, t); }
  fpreal      RADIUS(fpreal t)          { return evalFloat("radius", 0, t); }
  fpreal      CX(fpreal t) { return evalFloat("center", 0, t); }
  fpreal      CY(fpreal t) { return evalFloat("center", 1, t); }
  fpreal      CZ(fpreal t) { return evalFloat("center", 2, t); }
  int        INTER_SRC(fpreal t) {return evalInt("inter_src", 0, t);}

  const GA_PointGroup *myGroup;
  
};

#endif
