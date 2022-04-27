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
 * Deform_Surface interactive (aperiodic) SOP
 *----------------------------------------------------------------------------
 * At each frame, sum the contribution of all the inputs for all wavelengths
 * at each point to get the height.
 */


#include "SOP_Deform_Surface_inter.hpp"

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
  table->addOperator(new OP_Operator("def_surf_fs_inter",
				     "Deform Surface_inter FS",
				     SOP_Deform_Surface_inter::myConstructor,
				     SOP_Deform_Surface_inter::myTemplateList,
				     2,
				     OP_MULTI_INPUT_MAX,
				     0));
}

static PRM_Name names[] = {
  PRM_Name("amp",     "Amplitude"),
  PRM_Name("parallel", "Parallel")
};

PRM_Template
SOP_Deform_Surface_inter::myTemplateList[] = {
  PRM_Template(PRM_STRING,    1, &PRMgroupName, 0, &SOP_Node::pointGroupMenu,
	       0, 0, SOP_Node::getGroupSelectButton(
						    GA_GROUP_POINT)),
  PRM_Template(PRM_FLT_J,     1, &names[0], PRMoneDefaults, 0,
	       &PRMscaleRange),
  PRM_Template(PRM_TOGGLE_J,  1, &names[1]),
  PRM_Template(),
};


OP_Node *SOP_Deform_Surface_inter::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
  return new SOP_Deform_Surface_inter(net, name, op);
}

SOP_Deform_Surface_inter::SOP_Deform_Surface_inter(OP_Network *net, const char *name, OP_Operator *op)
  : SOP_Node(net, name, op) {
  mySopFlags.setManagesDataIDs(true);
  time.init();
}

SOP_Deform_Surface_inter::~SOP_Deform_Surface_inter() {
}
OP_ERROR
SOP_Deform_Surface_inter::cookInputGroups(OP_Context &context, int alone) {
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


void SOP_Deform_Surface_inter::addContribFromPrimitivePartial(const GU_Detail *fs,
 							      GA_Offset prim_off,
 							      const UT_JobInfo &info) {
  // void SOP_Deform_Surface_inter::addContribFromPrimitive(const GU_Detail *fs,
  // 								GA_Offset prim_off) {
  float dt = 0.1/3.0;
  float t = dt*fr;
  float amp = AMP(t);
  GA_ROHandleF w_handle(fs->findAttribute(GA_ATTRIB_PRIMITIVE, "wavelengths"));
  GA_ROHandleI as_handle(fs->findAttribute(GA_ATTRIB_PRIMITIVE, "ampli_steps"));
  const GA_Attribute *afs = fs->findFloatTuple(GA_ATTRIB_POINT, "ampli", buffer_size);
  const GA_AIFTuple *tuple;
  if (!w_handle.isValid()) {
    addError(SOP_ATTRIBUTE_INVALID, "wavelengths...!");
    return;
    //  return error();
  }
  if (!as_handle.isValid()) {
    addError(SOP_ATTRIBUTE_INVALID, "ampli_steps");
    return;
    // return error();
  }
  float wl = w_handle.get(prim_off);
  //  std::cout<<"WL  "<<wl<<std::endl;
  int as = as_handle.get(prim_off);
  float k = M_PI*2.0/wl;
  float om = omega(k);
  if (wl == 0) {// TEST get wl 0
    k = 0;
    om = 0;
  }
  const GA_Primitive* prim = fs->getPrimitive(prim_off);
	
  GA_Range range = prim->getPointRange();
  GA_Offset ptoff;
  GA_Range range_grid = gdp->getPointRange();
  int i= 0, n;
  //GA_Iterator itgr = range_grid.begin();
  //for(GA_Iterator itgr = range_grid.begin(); itgr != range_grid.end(); ++itgr) {
  for (info.divideWork(n_grid, i, n); i < n; i++) {
    ptoff = gdp->pointOffset(i);
    // std::cout<<"i "<<i<<" "<<amp<<" "<<as<<" "<<k<<" "<<om<<std::endl;
    //  	GA_FOR_ALL_PTOFF(gdp, ptoff) {
    UT_Vector3 Pvalue = gdp->getPos3(ptoff);
    //std::cout<<"Pval "<<Pvalue<<std::endl;
    COMPLEX a = 0;//exp(std::complex<float>(0, 1)*(k*Pvalue.x()));
    for(GA_Iterator it = range.begin(); it != range.end(); ++it) {
      UT_Vector3 P_fs = fs->getPos3((*it));
      float r = sqrt(pow(Pvalue.x() - P_fs.x(), 2) + pow(Pvalue.z() - P_fs.z(), 2));
      float v = velocity(k, om);//0.5*omega/k;
      float ar = 0, ai = 0;
      float arn = 0, ain = 0;
      //  float arn = 0, ain = 0;
	    
      fpreal t_ret = t - r/v;
      OP_Context c_ret(t_ret);
      // float ret = t_ret/dt;
      // int f_ret = floor(ret) - 1;
      int f_ret = floor(t_ret/(dt)) - 1;
      float coef = t_ret/(dt) - floor(t_ret/(dt));
      // float q =  ret - f_ret;
      //	    std::cout<<"q  "<<q<<std::endl;
      f_ret = floor((float)f_ret/(float)as); 
      if (t_ret < 0) {
	--f_ret;
      }
      if (f_ret >= 0) {
	afs = fs->findFloatTuple(GA_ATTRIB_POINT, "ampli", buffer_size);
	tuple = afs->getAIFTuple();
	tuple->get(afs, *it, ar, 2*f_ret);
	tuple->get(afs, *it, ai, 2*f_ret+1);
	tuple->get(afs, *it, arn, 2*f_ret+2);
	tuple->get(afs, *it, ain, 2*f_ret+3);
      }
      // if (itgr == range_grid.begin()) {
      //   std::cout<<"Pval "<<Pvalue<<std::endl;
      //   std::cout<<"P_fs "<<P_fs<<" "<<t<<" "<<r<<" "<<t_ret<<" "<<f_ret<<std::endl;
      // }

      // if (f_ret >= 0) {
      //   afs = fs->findFloatTuple(GA_ATTRIB_POINT, "ampli", buffer_size);
      //   tuple = afs->getAIFTuple();
      //   tuple->get(afs, *it, arn, 2*f_ret);
      //   tuple->get(afs, *it, ain, 2*f_ret+1);
      // }
      std::complex<float> ampli((1-coef)*ar + coef*arn, (1-coef)*ai + coef*ain);
      //	    std::complex<float> amplin(arn, ain);
      a += ampli*fund_solution(k*r)*damping(damping_coef, r, k);
    }
    // if (real(a) != 0 || imag(a) != 0) {
    //    std::cout<<"i "<<i<<" "<<a<<std::endl;
    //  }
    Pvalue.y() += real(amp*a*exp(-COMPLEX(0, 1)*(om*(float)t)));
    {UT_AutoJobInfoLock alock(info);
      gdp->setPos3(ptoff, Pvalue);
    }
    //    ++itgr;

  }
}


OP_ERROR SOP_Deform_Surface_inter::cookMySop(OP_Context &context) {
  OP_AutoLockInputs inputs(this);
  if (inputs.lock(context) >= UT_ERROR_ABORT)
    return error();
    
  flags().setTimeDep(1);
  float t = context.getTime();
  // int
  fr = context.getFrame();
  float dt = 0.1/3.0;
  t = dt*fr;

  // if (fr == 1) {
  //   _cuda_surface = CudaWaterSurface();
  //   // _cuda_surface.init(1);
  //   _cuda_surface.clear();
  // }
  time.tick(Times::deform_time_);
  float amp = AMP(t);
  int nb_inputs = getInputsArraySize();
  
  duplicateSource(0, context); //grid
  GA_RWHandleV3 Phandle(gdp->findAttribute(GA_ATTRIB_POINT, "P"));
  
  int nb_wl = inputGeo(1)->getPrimitiveRange().getEntries();
  GA_Offset ptoff;
  GA_FOR_ALL_PTOFF(gdp, ptoff) {
    UT_Vector3 Pvalue = gdp->getPos3(ptoff);
    Pvalue.y() = 0;
    gdp->setPos3(ptoff, Pvalue);
  }


  for (int input = 1; input < nb_inputs; ++input) { //all these inputs should be sets of sources
    const GU_Detail *fs = inputGeo(input);
    std::cout<<"Handling "<<input<<std::endl; 
    GA_ROHandleF w_handle(fs->findAttribute(GA_ATTRIB_PRIMITIVE, "wavelengths"));
    GA_ROHandleI as_handle(fs->findAttribute(GA_ATTRIB_PRIMITIVE, "ampli_steps"));
    if (!w_handle.isValid()) {
      addWarning(SOP_ATTRIBUTE_INVALID, "wavelengths...");
      std::cout<<"break "<<input<<std::endl; 
      break;
      //      return error();
    }
    if (!as_handle.isValid()) {
      addError(SOP_ATTRIBUTE_INVALID, "ampli_steps");
      return error();
    }
    GA_ROHandleI bs_handle(fs->findAttribute(GA_ATTRIB_DETAIL, "buffer_size"));
    GA_ROHandleF damping_handle(fs->findAttribute(GA_ATTRIB_DETAIL, "damping"));
    if (!bs_handle.isValid()) {
      addError(SOP_ATTRIBUTE_INVALID, "buffer sizes input sources");
      return error();
    }
    if (!damping_handle.isValid()) {
      addError(SOP_ATTRIBUTE_INVALID, "damping input sources");
      return error();
    }
    buffer_size = bs_handle.get(0);
    damping_coef = damping_handle.get(0);
          
    const GA_Attribute *afs = fs->findFloatTuple(GA_ATTRIB_POINT, "ampli", buffer_size);
    if (!afs)	{
      addError(SOP_ATTRIBUTE_INVALID, "ampli");
      return error();
    }
    const GA_AIFTuple *tuple = afs->getAIFTuple();
    GA_Offset ptoff;
    if (!tuple) {
      addError(SOP_ATTRIBUTE_INVALID, "ampli tuple");
      return error();
    }
      
    GA_Offset prim_off;
    GA_Range range_prim = fs->getPrimitiveRange();
    int w = 0;
    // for (GA_Iterator lcl_it((fs)->getPrimitiveRange()); lcl_it.blockAdvance(lcl_start, lcl_end); ) {
    //   for (prim_off = lcl_start; prim_off < lcl_end; ++prim_off) {
    n_grid = gdp->getPointRange().getEntries();
    //    std::cout<<"n_grid "<<n_grid<<std::endl;
    if (PARALLEL(t)) {
      for (GA_Iterator itp = range_prim.begin(); itp != range_prim.end(); ++itp) {
	prim_off = *itp;
	addContribFromPrimitive(fs, prim_off);
      }
    } else {
    for (GA_Iterator itp = range_prim.begin(); itp != range_prim.end(); ++itp) {
      prim_off = *itp;
      //      addContribFromPrimitive(fs, prim_off);
      float wl = w_handle.get(prim_off);
      //      std::cout<<"WL  "<<wl<<std::endl;
      int as = as_handle.get(prim_off);
      float k = M_PI*2.0/wl;
      float om = omega(k);//sqrtf(9.81*k + 0.074/1000*pow(k, 3));
      const GA_Primitive* prim = fs->getPrimitive(prim_off);
      GA_Range range = prim->getPointRange();
	
      GA_Offset ptoff;
      GA_Range range_grid = gdp->getPointRange();
      for(GA_Iterator itgr = range_grid.begin(); itgr != range_grid.end(); ++itgr) {
      	//  	GA_FOR_ALL_PTOFF(gdp, ptoff) {
      	UT_Vector3 Pvalue = gdp->getPos3((*itgr));
      	COMPLEX a = 0;//exp(std::complex<float>(0, 1)*(k*Pvalue.x()));
      	for(GA_Iterator it = range.begin(); it != range.end(); ++it) {
      	  UT_Vector3 P_fs = fs->getPos3((*it));
      	  float r = sqrt(pow(Pvalue.x() - P_fs.x(), 2) + pow(Pvalue.z() - P_fs.z(), 2));
      	  float v = velocity(k, om);//0.5*omega/k;
      	  float ar = 0, ai = 0;
      	  //  float arn = 0, ain = 0;
	    
      	  fpreal t_ret = t - r/v;
      	  OP_Context c_ret(t_ret);
      	  // float ret = t_ret/dt;
      	  // int f_ret = floor(ret) - 1;
      	  int f_ret = floor(t_ret/(dt)) - 1;
      	  // float q =  ret - f_ret;
      	  //	    std::cout<<"q  "<<q<<std::endl;
      	  f_ret = floor((float)f_ret/(float)as); 
      	  if (t_ret < 0) {
      	    --f_ret;
      	  }
      	  if (f_ret >= 0) {
      	    afs = fs->findFloatTuple(GA_ATTRIB_POINT, "ampli", buffer_size);
      	    tuple = afs->getAIFTuple();
      	    tuple->get(afs, *it, ar, 2*f_ret);
      	    tuple->get(afs, *it, ai, 2*f_ret+1);
      	  }
      	  // if (f_ret >= 0) {
      	  //   afs = fs->findFloatTuple(GA_ATTRIB_POINT, "ampli", buffer_size);
      	  //   tuple = afs->getAIFTuple();
      	  //   tuple->get(afs, *it, arn, 2*f_ret);
      	  //   tuple->get(afs, *it, ain, 2*f_ret+1);
      	  // }
      	  std::complex<float> ampli(ar, ai);
      	  //	    std::complex<float> amplin(arn, ain);
      	  //std::cout<<"damping_coef "<<damping_coef<<"  damping "<<damping(damping_coef, r, k)<<std::endl;
      	  a += ampli*fund_solution(k*r)*damping(damping_coef, r, k);
      	}
      	Pvalue.y() += real(amp*a*exp(-COMPLEX(0, 1)*(om*(float)t)));
      	gdp->setPos3((*itgr), Pvalue);
      }
      ++w;
    }
    }
  }
  Phandle.bumpDataId();
  time.tock(Times::deform_time_);
  std::cout<<"Deform time : "<<time.getTime(Times::deform_time_)<<std::endl;
  time.next_loop();
  return error();
}


