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
 * Normalize height SOP (make the average y-coord = 0)
 *----------------------------------------------------------------------------
 * At each frame, sum the contribution of all the inputs for all wavelengths
 * at each point to get the height.
 */


#include "SOP_NormalizeHeight.hpp"

#include <GU/GU_Detail.h>
#include <GA/GA_Handle.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_Director.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <UT/UT_DSOVersion.h>
#include <SYS/SYS_Math.h>
#include <iostream>

void newSopOperator(OP_OperatorTable *table) {
  table->addOperator(new OP_Operator("normalize_height_fs",
				     "Normalize Height FS",
				     SOP_NormalizeHeight::myConstructor,
				     SOP_NormalizeHeight::myTemplateList,
				     1,
				     1,
				     0));
}

static PRM_Name names[] = {
  PRM_Name("amp",     "Amplitude"),
  PRM_Name("parallel", "Parallel")
};

PRM_Template
SOP_NormalizeHeight::myTemplateList[] = {
  PRM_Template(PRM_STRING,    1, &PRMgroupName, 0, &SOP_Node::pointGroupMenu,
	       0, 0, SOP_Node::getGroupSelectButton(
						    GA_GROUP_POINT)),
  PRM_Template(PRM_FLT_J,     1, &names[0], PRMoneDefaults, 0,
	       &PRMscaleRange),
  PRM_Template(PRM_TOGGLE_J,  1, &names[1]),
  PRM_Template(),
};


OP_Node *SOP_NormalizeHeight::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
  return new SOP_NormalizeHeight(net, name, op);
}

SOP_NormalizeHeight::SOP_NormalizeHeight(OP_Network *net, const char *name, OP_Operator *op)
  : SOP_Node(net, name, op) {
  mySopFlags.setManagesDataIDs(true);
}

SOP_NormalizeHeight::~SOP_NormalizeHeight() {
}
OP_ERROR
SOP_NormalizeHeight::cookInputGroups(OP_Context &context, int alone) {
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



OP_ERROR SOP_NormalizeHeight::cookMySop(OP_Context &context) {
  OP_AutoLockInputs inputs(this);
  if (inputs.lock(context) >= UT_ERROR_ABORT)
    return error();
    
  flags().setTimeDep(1);
  
  duplicateSource(0, context); 
  GA_RWHandleV3 Phandle(gdp->findAttribute(GA_ATTRIB_POINT, "P"));

  GA_Offset ptoff;
  float sum = 0;
  int nbp = 0;
  { GA_FOR_ALL_PTOFF(gdp, ptoff) {
    UT_Vector3 Pvalue = gdp->getPos3(ptoff);
    sum += Pvalue.y();
    ++nbp;
    }}
  float avg = sum / (float)nbp;
  
  { GA_FOR_ALL_PTOFF(gdp, ptoff) {
    UT_Vector3 Pvalue = gdp->getPos3(ptoff);
    Pvalue.y() -= avg;
    gdp->setPos3(ptoff, Pvalue);
    }}

  return error();
}


