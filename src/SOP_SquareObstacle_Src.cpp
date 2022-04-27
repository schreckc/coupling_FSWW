/*
 * MIT License
 * 
 * Copyright (c) 2019 Camille Schreck
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
 * Square Obstacle SOP
 *---------------------------------------------------------------------------
 * Create set of point sources along an offset surface of the obstacle.
 * Range of wavelength (and time step per wl), and is copied from input geometry, as well as
 *    the detail attibute.
 * Create on primitve and subset of sources for each wl.
 * Spacing depends on the wavelength and the parameter "density".
 * I use this node also to create a set of point sampling the border of the obstacle (offset=0)
 *    for the boundary conditions.
 * Note: the density of the boundary points should be at least twice the density of the
 *    sources.
 * Note 2: for the aperiodic version, do not forget to check the "interactive sources" box in
 *   the parameter of the node creating the sources of the obstacle.
 */

#include "SOP_SquareObstacle_Src.hpp"

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
  table->addOperator(new OP_Operator("square_obstacle_src_fs",
				     "Square Obstacle_Src Sources FS",
				     SOP_Square_Obstacle_Src::myConstructor,
				     SOP_Square_Obstacle_Src::myTemplateList,
				     1,
				     1,
				     nullptr,  
				     OP_FLAG_GENERATOR));
}
static PRM_Name names[] = {
  PRM_Name("center",  "Center"),
  PRM_Name("off",     "Offset distance"),
  PRM_Name("density",   "Density"),
  PRM_Name("length",   "Length"),
  PRM_Name("width",   "Width"),
  PRM_Name("inter_src",   "Interactive sources"),
};

PRM_Default* off_default = new PRM_Default(0.3);
PRM_Default* dens_default = new PRM_Default(3);

PRM_Template
SOP_Square_Obstacle_Src::myTemplateList[] = {
  PRM_Template(PRM_STRING,    1, &PRMgroupName, 0, &SOP_Node::pointGroupMenu,
	       0, 0, SOP_Node::getGroupSelectButton(GA_GROUP_POINT)),
  PRM_Template(PRM_XYZ_J,     3, &names[0], PRMzeroDefaults),
  PRM_Template(PRM_FLT_J,     1, &names[1], off_default),
  PRM_Template(PRM_FLT_J,     1, &names[2], dens_default),
  PRM_Template(PRM_FLT_J,     1, &names[3], PRMoneDefaults),
  PRM_Template(PRM_FLT_J,     1, &names[4], PRMoneDefaults),
  PRM_Template(PRM_TOGGLE_J,    1, &names[5]),
  PRM_Template(),
};


OP_Node *SOP_Square_Obstacle_Src::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
  return new SOP_Square_Obstacle_Src(net, name, op);
}

SOP_Square_Obstacle_Src::SOP_Square_Obstacle_Src(OP_Network *net, const char *name, OP_Operator *op)
  : SOP_Node(net, name, op) {
  mySopFlags.setManagesDataIDs(true);
}

SOP_Square_Obstacle_Src::~SOP_Square_Obstacle_Src() {
}

OP_ERROR SOP_Square_Obstacle_Src::cookInputGroups(OP_Context &context, int alone) {
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


/*OP_ERROR SOP_Square_Obstacle_Src::cookMySop(OP_Context &context)
  {
  flags().setTimeDep(0);
  float t = context.getTime();
  bool is_inter = INTER_SRC(t);
  
  OP_AutoLockInputs inputs(this);
  if (inputs.lock(context) >= UT_ERROR_ABORT)
  return error();
  
  gdp->clearAndDestroy();

  const GU_Detail *is = inputGeo(0); //input sources (used to get wavelength)
  int nb_wl = is->getPrimitiveRange().getEntries();

  std::vector<float> wave_lengths(nb_wl);
  std::vector<int> ampli_steps(nb_wl);

  GA_ROHandleF w_handle(is->findAttribute(GA_ATTRIB_PRIMITIVE, "wavelengths"));
  GA_ROHandleF as_handle(is->findAttribute(GA_ATTRIB_PRIMITIVE, "ampli_steps"));
  if (!w_handle.isValid()) {
  addError(SOP_ATTRIBUTE_INVALID, "wavelengths input sources");
  return error();
  }
  if (!as_handle.isValid()) {
  addError(SOP_ATTRIBUTE_INVALID, "amplis_steps input sources");
  return error();
  }

  GA_ROHandleI bs_handle(is->findAttribute(GA_ATTRIB_DETAIL, "buffer_size"));
  GA_ROHandleF damping_handle(is->findAttribute(GA_ATTRIB_DETAIL, "damping"));
  if (!bs_handle.isValid()) {
  addError(SOP_ATTRIBUTE_INVALID, "buffer sizes input sources");
  return error();
  }
  if (!damping_handle.isValid()) {
  addError(SOP_ATTRIBUTE_INVALID, "damping input sources");
  return error();
  }
  int buffer_size = 2;
  if (is_inter) {
  buffer_size= bs_handle.get(0);
  }
  float damping = damping_handle.get(0);
  
  int w = 0;
  GA_Range range_is = is->getPrimitiveRange();
  for(GA_Iterator itis = range_is.begin(); itis != range_is.end(); ++itis, ++w) {
  GA_Offset prim_off = *itis;
  float wl = w_handle.get(prim_off);
  wave_lengths[w] = wl;
  ampli_steps[w]= as_handle.get(prim_off);
  }

  //create set of points for each wavelength, and liked them to their corresponding primitve
  for (int w = 0; w < nb_wl; ++w) {
  float wl = wave_lengths[w];
  float density = DENSITY(t)/wl;
  float off = OFF(t);
  float width = WIDTH(t) - 2*off;;
  float length = LENGTH(t) - 2*off;;
  int nb_points_w = width*density + 1;
  int nb_points_l = length*density + 1;
 
  VEC3 center(CX(t), CY(t), CZ(t));
  std::vector<VEC3> corners(4);
  corners[0] = center + VEC3(width/2.0, 0, length/2.0);
  corners[1] = center + VEC3(-width/2.0, 0, length/2.0);
  corners[2] = center + VEC3(-width/2.0, 0, -length/2.0);
  corners[3] = center + VEC3(width/2.0, 0, -length/2.0);
  GA_Offset ptoff = gdp->appendPointBlock(2*nb_points_w + 2*nb_points_l);
  GA_Offset vtxoff;
  GA_Offset prim_off = gdp->appendPrimitivesAndVertices(GA_PRIMPOLY, 1, 2*nb_points_w + 2*nb_points_l, vtxoff, true);
  int i = 0;
  for (int c = 0; c < 4; ++c) {
      
  VEC3 dir = corners[(c+1)%4] - corners[c];
  float l = dir.norm();
  int nb_points = l*density + 1;
  float step = 1.0/(float)nb_points;
  float d = 0;
  for (int j = 0; j < nb_points;  ++i, ++j) {
  VEC3 pos =  corners[c] + d*dir;
  gdp->getTopology().wireVertexPoint(vtxoff+i,ptoff+i);
  gdp->setPos3(ptoff+i, UT_Vector3(pos(0), pos(1), pos(2)));
  d += step; 
  }
  }

  }
  
  // create and set attibutes for primitves and detail
  GA_RWHandleF wl_attrib(gdp->findFloatTuple(GA_ATTRIB_PRIMITIVE, "wavelengths", 1));
  GA_RWHandleF as_attrib(gdp->findFloatTuple(GA_ATTRIB_PRIMITIVE, "ampli_steps", 1));
  GA_RWHandleF bs_attrib(gdp->findFloatTuple(GA_ATTRIB_DETAIL, "buffer_size", 1));
  GA_RWHandleF damping_attrib(gdp->findFloatTuple(GA_ATTRIB_DETAIL, "damping", 1));
  if (!wl_attrib.isValid()) {
  wl_attrib = GA_RWHandleF(gdp->addFloatTuple(GA_ATTRIB_PRIMITIVE, "wavelengths", 1));
  }
  if (!wl_attrib.isValid()) {
  addError(SOP_MESSAGE, "Failed to create attribute wavelengths");
  return error();
  }
  if (!as_attrib.isValid()) {
  as_attrib = GA_RWHandleF(gdp->addFloatTuple(GA_ATTRIB_PRIMITIVE, "ampli_steps", 1));
  }
  if (!as_attrib.isValid()) {
  addError(SOP_MESSAGE, "Failed to create attribute ampli_steps");
  return error();
  }
  if (!bs_attrib.isValid()) {
  bs_attrib = GA_RWHandleF(gdp->addFloatTuple(GA_ATTRIB_DETAIL, "buffer_size", 1));
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

  bs_attrib.set(0, buffer_size);
  damping_attrib.set(0, damping);

  w = 0;
  GA_Offset prim_off, lcl_start, lcl_end;
  for (GA_Iterator lcl_it((gdp)->getPrimitiveRange()); lcl_it.blockAdvance(lcl_start, lcl_end); ) {
  for (prim_off = lcl_start; prim_off < lcl_end; ++prim_off) {
  wl_attrib.set(prim_off, wave_lengths[w]);
  as_attrib.set(prim_off, ampli_steps[w]);
  ++w;
  }
  }

  //create amplitude buffer attributes for the point sources and set them to zero  
  GA_RWHandleF ampli_attrib(gdp->findFloatTuple(GA_ATTRIB_POINT, "ampli", buffer_size));
  if (!ampli_attrib.isValid()) {
  ampli_attrib = GA_RWHandleF(gdp->addFloatTuple(GA_ATTRIB_POINT, "ampli", buffer_size));
  }
  if (!ampli_attrib.isValid()) {
  addError(SOP_MESSAGE, "Failed to create attribute ampli");
  return error();
  }
  {
  GA_Offset ptoff; 
  GA_FOR_ALL_PTOFF(gdp, ptoff) {
  for (uint i = 0; i < buffer_size; ++i) {
  ampli_attrib.set(ptoff, i, 0);
  }
  }
  }
  gdp->bumpDataIdsForAddOrRemove(true, true, true);

  return error();
  }
*/


OP_ERROR SOP_Square_Obstacle_Src::cookMySop(OP_Context &context) {
  OP_AutoLockInputs inputs(this);
  if (inputs.lock(context) >= UT_ERROR_ABORT)
    return error();

  flags().setTimeDep(0);
  float t = context.getTime();
  //  bool is_inter = INTER_SRC(t);

  gdp->clearAndDestroy();

  // get details and primitives attibutes from the input sources
  duplicateSource(0, context);

  int nb_inputs = getInputsArraySize();
  int nb_wl = gdp->getPrimitiveRange().getEntries();
    
  GA_ROHandleI bs_handle(gdp->findAttribute(GA_ATTRIB_DETAIL, "buffer_size"));
  // GA_ROHandleF damping_handle(is->findAttribute(GA_ATTRIB_DETAIL, "damping"));
  if (!bs_handle.isValid()) {
    addError(SOP_ATTRIBUTE_INVALID, "buffer sizes input sources");
    return error();
  }
  int buffer_size= bs_handle.get(0);
 
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
  //create set of points for each wavelength, and link them to their corresponding primitve
  //for (int w = 0; w < nb_wl; ++w) {
  GA_Range prange = gdp->getPrimitiveRange();
  for(GA_Iterator it = prange.begin(); it != prange.end(); ++it) {
    //float wl = wave_lengths[w];
    GA_Offset prim_off = *it;
    float wl = w_handle.get(prim_off);
    const GA_Primitive* prim = gdp->getPrimitive(prim_off);
    float density = DENSITY(t)/wl;
    float off = OFF(t);
    //
    // float radius = RADIUS(t) - wl*off;;
    // int nb_points = 2*M_PI*radius*density;
    // float angle_step = 2*M_PI/nb_points;
      float width = WIDTH(t) - 2*off;;
  float length = LENGTH(t) - 2*off;;
  int nb_points_w = width*density + 1;
  int nb_points_l = length*density + 1;

    // VEC3 dir(1, 0, 0);
    // VEC3 center(CX(t), CY(t), CZ(t));
    // MAT3 rotation;
    // rotation <<
    //   cos(angle_step), 0, -sin(angle_step),
    //   0, 0, 0, 
    //   sin(angle_step), 0, cos(angle_step);

    // GA_Offset ptoff = gdp->appendPointBlock(nb_points);

    // GA_Offset vtxoff;// = gdp->appendVertexBlock(nb_points);
	
    // GA_Offset prim_off_2 = gdp->appendPrimitivesAndVertices(GA_PRIMPOLY, 1, nb_points, vtxoff, true);
    int as = as_handle.get(prim_off);
    GA_Offset ptoff = gdp->appendPointBlock(2*nb_points_w + 2*nb_points_l);
  GA_Offset vtxoff;
  GA_Offset prim_off_2 = gdp->appendPrimitivesAndVertices(GA_PRIMPOLY, 1, 2*nb_points_w + 2*nb_points_l, vtxoff, true);
    w_handle.set(prim_off_2, wl);
    as_handle.set(prim_off_2, as);
    gdp->destroyPrimitiveOffset(prim_off);
    // for (int i = 0; i < nb_points;  ++i) {
    //   VEC3 pos = center + radius*dir;
    //    gdp->getTopology().wireVertexPoint(vtxoff+i,ptoff+i);
    //   gdp->setPos3(ptoff+i, UT_Vector3(pos(0), pos(1), pos(2)));
    //   dir = rotation*dir;
    //   dir.normalize();
    // }
      VEC3 center(CX(t), CY(t), CZ(t));
  std::vector<VEC3> corners(4);
  corners[0] = center + VEC3(width/2.0, 0, length/2.0);
  corners[1] = center + VEC3(-width/2.0, 0, length/2.0);
  corners[2] = center + VEC3(-width/2.0, 0, -length/2.0);
  corners[3] = center + VEC3(width/2.0, 0, -length/2.0);
    int i = 0;
    for (int c = 0; c < 4; ++c) {
	  
      VEC3 dir = corners[(c+1)%4] - corners[c];
      float l = dir.norm();
      int nb_points = l*density + 1;
      float step = 1.0/(float)nb_points;
      float d = 0;
      for (int j = 0; j < nb_points;  ++i, ++j) {
	VEC3 pos =  corners[c] + d*dir;
	gdp->getTopology().wireVertexPoint(vtxoff+i,ptoff+i);
	gdp->setPos3(ptoff+i, UT_Vector3(pos(0), pos(1), pos(2)));
	d += step; 
      }
    }
	
  }


  //create amplitude buffer attributes for the point sources and set them to zero  
  GA_RWHandleF ampli_attrib(gdp->findFloatTuple(GA_ATTRIB_POINT, "ampli", buffer_size));
  if (!ampli_attrib.isValid()) {
    ampli_attrib = GA_RWHandleF(gdp->addFloatTuple(GA_ATTRIB_POINT, "ampli", buffer_size));
  }
  if (!ampli_attrib.isValid()) {
    addError(SOP_MESSAGE, "Failed to create attribute ampli");
    return error();
  }
  {
    GA_Offset ptoff; 
    GA_FOR_ALL_PTOFF(gdp, ptoff) {
      for (uint i = 0; i < buffer_size; ++i) {
  	ampli_attrib.set(ptoff, i, 0);
      }
    }
  }
  gdp->bumpDataIdsForAddOrRemove(true, true, true);
  gdp->bumpAllDataIds();
  return error();
}
