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


#include "SOP_Boundary_Points.hpp"

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
  table->addOperator(new OP_Operator("boundary_points_fs",
				     "Boundary Points FS",
				     SOP_Boundary_Points::myConstructor,
				     SOP_Boundary_Points::myTemplateList,
				     2,
				     2,
				     nullptr,  
				     OP_FLAG_GENERATOR));
}

PRM_Template
SOP_Boundary_Points::myTemplateList[] = {
					 PRM_Template(PRM_STRING,    1, &PRMgroupName, 0, &SOP_Node::pointGroupMenu,
						      0, 0, SOP_Node::getGroupSelectButton(GA_GROUP_POINT)),
					 PRM_Template(),
};


OP_Node *
SOP_Boundary_Points::myConstructor(OP_Network *net, const char *name, OP_Operator *op) {
  return new SOP_Boundary_Points(net, name, op);
}

SOP_Boundary_Points::SOP_Boundary_Points(OP_Network *net, const char *name, OP_Operator *op)
  : SOP_Node(net, name, op) {
  mySopFlags.setManagesDataIDs(true);
}

SOP_Boundary_Points::~SOP_Boundary_Points()
{
  std::list<InputPoint>::iterator it;
  for (it = inputPoints.begin(); it !=inputPoints.end(); ++it) {  
    (*it).clear();
  }
}
OP_ERROR
SOP_Boundary_Points::cookInputGroups(OP_Context &context, int alone)
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

std::ofstream & SOP_Boundary_Points::record(std::ofstream & file) {
  std::list<InputPoint>::iterator it;
  for (it = inputPoints.begin(); it !=inputPoints.end(); ++it) {
    (*it).record_samples(file);
  }
  return file;
}

std::ifstream & SOP_Boundary_Points::read(std::ifstream & file) {
  std::list<InputPoint>::iterator it;
  for (it = inputPoints.begin(); it !=inputPoints.end(); ++it) {
    FLOAT next_sample;
    file >> next_sample;
    (*it).update(next_sample);
  }
  return file;
}

OP_ERROR SOP_Boundary_Points::cookMySop(OP_Context &context) {
  OP_AutoLockInputs inputs(this);
  if (inputs.lock(context) >= UT_ERROR_ABORT)
    return error();

  flags().setTimeDep(1);

  float t = context.getTime();
  int fr = context.getFrame();
  
  // TODO: user def time step
  float dt_ = 0.1/3.0;
  t = dt_*fr;

  const GU_Detail *fs = inputGeo(0);

  GA_ROHandleI ws_handle(inputGeo(1)->findAttribute(GA_ATTRIB_DETAIL, "winsize"));
  if (!ws_handle.isValid()) {
    addError(SOP_ATTRIBUTE_INVALID, "winsize");
    return error();
  }
  GA_ROHandleI su_handle(inputGeo(1)->findAttribute(GA_ATTRIB_DETAIL, "shift_u"));
  if (!su_handle.isValid()) {
    addError(SOP_ATTRIBUTE_INVALID, "shift_u");
    return error();
  }
  GA_ROHandleI sb_handle(inputGeo(1)->findAttribute(GA_ATTRIB_DETAIL, "shift_b"));
  if (!sb_handle.isValid()) {
    addError(SOP_ATTRIBUTE_INVALID, "shift_b");
    return error();
  }
  int winsize = ws_handle.get(0);
  int shift_u = su_handle.get(0);
  int shift_b = sb_handle.get(0);
  
  if (fr == 1) {
  
    nb_pts = fs->getPointRange().getEntries();
    GA_Range range = fs->getPointRange();

    // create InputPoint that record height of input simu
    uint i = 0;
    inputPoints.clear();
    uint nb_ip = range.getEntries();
    for(GA_Iterator itfs = range.begin(); itfs != range.end(); ++itfs) {
      UT_Vector3 pos_fs = fs->getPos3(*itfs);
      InputPoint ip(winsize, dt_);
      ip.setPos(pos_fs(0), pos_fs(2));
      inputPoints.push_back(ip);
      ++i;
    }
  
    //  begin creation of geometry
    gdp->clearAndDestroy();
    // get details and primitives attibutes from the input sources
    duplicateSource(1, context);
    int nb_wl = gdp->getPrimitiveRange().getEntries();
    GA_RWHandleF w_handle(gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "wavelengths"));
    if (!w_handle.isValid()) {
      addError(SOP_ATTRIBUTE_INVALID, "wavelengths");
      return error();
    }
    GA_RWHandleI as_handle(gdp->findIntTuple(GA_ATTRIB_PRIMITIVE, "ampli_steps", 1));
    if (!as_handle.isValid()) {
      addError(SOP_MESSAGE, "Failed to create attribute ampli_steps");
      return error();
    }
    
    GA_Range prange = gdp->getPrimitiveRange();
    for(GA_Iterator it = prange.begin(); it != prange.end(); ++it) {
      GA_Offset prim_off = *it;
      float wl = w_handle.get(prim_off);
      int as = as_handle.get(prim_off);
	
      GA_Offset ptoff = gdp->appendPointBlock(nb_pts);
      GA_Offset vtxoff;
      GA_Offset prim_off_2 = gdp->appendPrimitivesAndVertices(GA_PRIMPOLY, 1, nb_pts, vtxoff, true);
      GA_Range range = fs->getPointRange();
      uint i = 0;
      for(GA_Iterator itfs = range.begin(); itfs != range.end(); ++itfs) {
	UT_Vector3 pos_fs = fs->getPos3(*itfs);
	pos_fs(1) = 0;
	gdp->getTopology().wireVertexPoint(vtxoff+i,ptoff+i);
	gdp->setPos3(ptoff+i, pos_fs);
	++i;
      }

      w_handle.set(prim_off_2, wl);
      as_handle.set(prim_off_2, as);
      gdp->destroyPrimitiveOffset(prim_off);
    }
    //create amplitude buffer for each source and fill it
    GA_RWHandleF ampli_attrib(gdp->findFloatTuple(GA_ATTRIB_POINT, "ampli", 2));
    if (!ampli_attrib.isValid()) {
      ampli_attrib = GA_RWHandleF(gdp->addFloatTuple(GA_ATTRIB_POINT, "ampli", 2));
    }
    if (!ampli_attrib.isValid()) {
      addError(SOP_MESSAGE, "Failed to create attribute ampli");
      return error();
    }

    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(gdp, ptoff) {
      ampli_attrib.set(ptoff, 0, 0);
      ampli_attrib.set(ptoff, 1, 0);
    }

  }
  if (nb_pts != fs->getPointRange().getEntries()) {
    addWarning(SOP_ATTRIBUTE_INVALID, "not the right number of boundary points");
    return error();
  }

  GA_Range rangefs = fs->getPointRange();
  GA_Range rangegdp = gdp->getPointRange();
  GA_Iterator itfs = rangefs.begin();
  std::list<InputPoint>::iterator it = inputPoints.begin();
  
  if (rangefs.getEntries() != inputPoints.size()) {
    addWarning(SOP_ATTRIBUTE_INVALID, "not the right number of inputpoints");
  }
  for(; itfs != rangefs.end() && it != inputPoints.end(); ++itfs, ++it) {
    UT_Vector3 pos_fs = fs->getPos3(*itfs);
    (*it).update(pos_fs(1));
  }
 
  GA_RWHandleF w_handle(gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "wavelengths"));
  if (!w_handle.isValid()) {
    addError(SOP_ATTRIBUTE_INVALID, "wavelengths");
    return error();
  }
  GA_RWHandleI as_handle(gdp->findIntTuple(GA_ATTRIB_PRIMITIVE, "ampli_steps", 1));
  if (!as_handle.isValid()) {
    addError(SOP_MESSAGE, "Failed to create attribute ampli_steps");
    return error();
  }
    
  GA_RWHandleF ampli_attrib(gdp->findFloatTuple(GA_ATTRIB_POINT, "ampli", 2));
  if (!ampli_attrib.isValid()) {
    addError(SOP_MESSAGE, "cannot find attribute ampli");
    return error();
  }
	     
  GA_Offset ptoff;
  GA_Range range_i = gdp->getPrimitiveRange();
  int w = 0;
  for(GA_Iterator itp = range_i.begin(); itp != range_i.end(); ++itp, ++w) {
    GA_Offset prim_off = *itp;
    float wl = w_handle.get(prim_off);
    int as = as_handle.get(prim_off);
     
    float k = M_PI*2.0/wl;
    float om = omega(k);

    if (wl == 0) {
      k = 0;
      om = 0;
    }
      
    const GA_Primitive* prim = gdp->getPrimitive(*itp);
    GA_Range range_bp = prim->getPointRange();
    it = inputPoints.begin();
    for(GA_Iterator itbp = range_bp.begin();
	itbp != range_bp.end() && it != inputPoints.end(); ++itbp, ++it) {
      COMPLEX a = COMPLEX((*it).spectrum_re[w+shift_u+1], (*it).spectrum_im[w+shift_u+1]); 
      a /= ((FLOAT)winsize/2.0);
      a *= -(float)((w+shift_u+1)%2*2-1)*exp(COMPLEX(0, 1)*(om*(float)fr*dt_));
      ampli_attrib.set(*itbp, 0, real(a));
      ampli_attrib.set(*itbp, 1, imag(a));
    }
  }

  return error();
}
