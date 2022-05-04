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
 * Deform_Surface interactive (aperiodic) SOP
 *----------------------------------------------------------------------------
 * At each frame, sum the contribution of all the inputs for all wavelengths
 * at each point to get the height.
 */


#ifndef __SOP_Deform_Surface_inter_h__
#define __SOP_Deform_Surface_inter_h__

#include <SOP/SOP_Node.h>
#include <UT/UT_ThreadedAlgorithm.h>
#include "Time.hpp"


class SOP_Deform_Surface_inter : public SOP_Node {
public:
    SOP_Deform_Surface_inter(OP_Network *net, const char *name, OP_Operator *op);
    virtual ~SOP_Deform_Surface_inter();

    static PRM_Template myTemplateList[];
    static OP_Node *myConstructor(OP_Network*, const char *, OP_Operator *);

    virtual OP_ERROR             cookInputGroups(OP_Context &context, 
						 int alone = 0);

    int n_grid;
  THREADED_METHOD2(                   // Construct two parameter threaded method
                    SOP_Deform_Surface_inter,                // Name of class
                    n_grid > 100,     // Evaluated to see if we should multithread.
                    addContribFromPrimitive,                // Name of function
                    const GU_Detail *, fs,            // An integer parameter named p1
                    GA_Offset, prim_off)          // A float parameter named p2
      void addContribFromPrimitivePartial(const GU_Detail *fs, GA_Offset prim_off, const UT_JobInfo &info);
protected:
  //  void addContribFromPrimitive(const GU_Detail *fs, GA_Offset prim_off);
    virtual OP_ERROR cookMySop(OP_Context &context);
private:
    void        getGroups(UT_String &str){ evalString(str, "group", 0, 0); }
    fpreal      AMP(fpreal t)           { return evalFloat("amp", 0, t); }
  int        PARALLEL(fpreal t) {return evalInt("parallel", 0, t);}

  const GA_PointGroup *myGroup;
  int fr;
  int buffer_size;
  float damping_coef;
  Times time;
};

#endif
