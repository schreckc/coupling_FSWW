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
 * Boundary_Points SOP
 *----------------------------------------------------------------------------
 * Record the height at the position of the boundary points (input 0) to compute the amplitude of each frequency component using the wave parameters (input 1).
 */


#ifndef __SOP_Boundary_Points_h__
#define __SOP_Boundary_Points_h__

#include <SOP/SOP_Node.h>
#include "InputPoint.hpp"

class SOP_Boundary_Points : public SOP_Node {
  
public:
  SOP_Boundary_Points(OP_Network *net, const char *name, OP_Operator *op);
  virtual ~SOP_Boundary_Points();

  static PRM_Template myTemplateList[];
  static OP_Node *myConstructor(OP_Network*, const char *, OP_Operator *);

  virtual OP_ERROR             cookInputGroups(OP_Context &context, 
					       int alone = 0);

protected:
  std::ofstream & record(std::ofstream & file);
  std::ifstream & read(std::ifstream & file);
  virtual OP_ERROR cookMySop(OP_Context &context);
private:
  void        getGroups(UT_String &str){ evalString(str, "group", 0, 0); }
  
  std::list<InputPoint> inputPoints;
  std::vector<float> wave_lengths;
  std::vector<int> ampli_steps;
  int nb_pts;
  
  const GA_PointGroup *myGroup;
};

#endif
