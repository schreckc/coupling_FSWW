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
 * Solve FS interactive (aperiodic) SOP
 *----------------------------------------------------------------------------
 * Compute amplitude of the obstacle sources (input 0) such that the boundary conditions computed at the boundary points (input 1) given the wave parameters (input 2)
 */

#include "SOP_Solve_FS_inter.hpp"

#include <GU/GU_Detail.h>
#include <GA/GA_Handle.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_Director.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <UT/UT_DSOVersion.h>
#include <SYS/SYS_Math.h>


void newSopOperator(OP_OperatorTable *table) {
  table->addOperator(new OP_Operator("solve_fs_inter",
				     "Solve FS Interactif",
				     SOP_Solve_FS_inter::myConstructor,
				     SOP_Solve_FS_inter::myTemplateList,
				     3,
				     4,
				     nullptr,  
				     OP_FLAG_GENERATOR));
}

static PRM_Name names[] = {
			   PRM_Name("obstacle",  "Obstacle"),
};

PRM_Template SOP_Solve_FS_inter::myTemplateList[] = {
						     PRM_Template(PRM_STRING,    1, &PRMgroupName, 0, &SOP_Node::pointGroupMenu,
								  0, 0, SOP_Node::getGroupSelectButton(
												       GA_GROUP_POINT)),
						     PRM_Template(PRM_TOGGLE_J,  1, &names[0]),
						     PRM_Template(),
};


OP_Node *SOP_Solve_FS_inter::myConstructor(OP_Network *net, const char *name, OP_Operator *op) {
  return new SOP_Solve_FS_inter(net, name, op);
}

SOP_Solve_FS_inter::SOP_Solve_FS_inter(OP_Network *net, const char *name, OP_Operator *op)
  : SOP_Node(net, name, op) {
  mySopFlags.setManagesDataIDs(true);
  time.init();
  nb_wl = 0;
}

SOP_Solve_FS_inter::~SOP_Solve_FS_inter() {

}
OP_ERROR
SOP_Solve_FS_inter::cookInputGroups(OP_Context &context, int alone) {
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

OP_ERROR SOP_Solve_FS_inter::cookMySop(OP_Context &context) {
  OP_AutoLockInputs inputs(this);
  if (inputs.lock(context) >= UT_ERROR_ABORT)
    return error();

  
  flags().setTimeDep(1);
  float t = context.getTime();
  int fr = context.getFrame();
  float dt = 0.1/3.0;
  t = dt*fr;
   
  time.tick(Times::solve_time_);
  gdp->clearAndDestroy();  
  duplicateSource(0, context);
  

  const GU_Detail *bp = inputGeo(1); //boundary points
  const GU_Detail *is = inputGeo(2); //input sources if any
  const GU_Detail *s = inputGeo(0);  // sources
    
  std::cout <<"Number of wavelengths ---- "<<s->getPrimitiveRange().getEntries()<<std::endl; 
  if (fr == 1 || nb_wl != s->getPrimitiveRange().getEntries()) {
    std::cout<<"Initiation ..."<<std::endl;
    nb_wl = s->getPrimitiveRange().getEntries();
   
    // get wavelenght range and time step for each wavelength from input sources
    wave_lengths = std::vector<float>(nb_wl);
    ampli_steps = std::vector<int>(nb_wl);
    p_in = std::vector<VectorXcf>(nb_wl);
    svd = std::vector<Eigen::BDCSVD<MatrixXcf> >(nb_wl);
    matrices = std::vector<MatrixXcf>(nb_wl);
    GA_ROHandleF w_handle(s->findAttribute(GA_ATTRIB_PRIMITIVE, "wavelengths"));
    GA_ROHandleI as_handle(s->findAttribute(GA_ATTRIB_PRIMITIVE, "ampli_steps"));
    if (!w_handle.isValid()) {
      addError(SOP_ATTRIBUTE_INVALID, "wavelengths input sources");
      return error();
    }
    if (!as_handle.isValid()) {
      addError(SOP_ATTRIBUTE_INVALID, "ampli_steps input sources");
      return error();
    }
            
    int w = 0;
    GA_Range range_is = s->getPrimitiveRange();
    for(GA_Iterator itis = range_is.begin(); itis != range_is.end(); ++itis, ++w) {
      GA_Offset prim_off = *itis;
      float wl = w_handle.get(prim_off);
      wave_lengths[w] = wl;
      ampli_steps[w]= as_handle.get(prim_off);
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
    buffer_size = bs_handle.get(0);
    damping_coef = damping_handle.get(0);

    shifts = std::vector<std::vector<int> >(nb_wl);
    
    // create transfer matrix
    for (int w = 0; w < nb_wl; ++w) {
      float wl = wave_lengths[w];
      int as = ampli_steps[w];
      float k = M_PI*2.0/wl;
      float v = velocity(k);
      const GA_Primitive* prim_fs = gdp->getPrimitiveByIndex(w);
      GA_Range range_fs = prim_fs->getPointRange();
      int nb_fs = range_fs.getEntries();
      const GA_Primitive* prim_bp = bp->getPrimitiveByIndex(w);
      GA_Range range_bp = prim_bp->getPointRange();
      int nb_bp = range_bp.getEntries();
      if (nb_fs <= 0) {
  	addError(SOP_ATTRIBUTE_INVALID, "no fundamental sources defined");
  	return error();
      }
      if (nb_bp <= 0) {
	nb_wl = 0;
  	return error();
      }
      shifts[w] = std::vector<int>(nb_fs);
      int ifs = 0;
      for(GA_Iterator itfs = range_fs.begin(); itfs != range_fs.end(); ++itfs, ++ifs) {
      	shifts[w][ifs] = 0;
      	UT_Vector3 pos_fs = gdp->getPos3(*itfs);
      	FLOAT d_min = distance3d(pos_fs, bp->getPos3(*(range_bp.begin())));
      	for(GA_Iterator itbp = range_bp.begin(); itbp != range_bp.end(); ++itbp) {
      	  UT_Vector3 pos_b = bp->getPos3(*itbp);
      	  FLOAT d_cur = distance3d(pos_fs, pos_b);
      	  if (d_cur < d_min) {
      	    d_min = d_cur;
      	  }
	
      	}
      	FLOAT ret = d_min/velocity(k);
      	int shift = floor(ret/(dt*as));
      	shifts[w][ifs] = shift;
      	float q = interpolation(ret, shift*(dt*as), (shift+1)*dt*as);
      }
	  
      MatrixXcf T = MatrixXcf(nb_bp, nb_fs);
      uint i = 0, j = 0;
      for(GA_Iterator itfs = range_fs.begin(); itfs != range_fs.end(); ++itfs) {
  	UT_Vector3 pos_fs = gdp->getPos3(*itfs);
  	i = 0;
  	uint nb = 0;
  	for(GA_Iterator itbp = range_bp.begin(); itbp != range_bp.end(); ++itbp) {
  	  UT_Vector3 pos_b = bp->getPos3(*itbp);
  	  float r = sqrt(pow(pos_b.x() - pos_fs.x(), 2) + pow(pos_b.z() - pos_fs.z(), 2));
  	  float ret = r/v;
  	  float q = interpolation(ret, shifts[w][j]*(dt*as), (shifts[w][j]+1)*dt*as);
  	  if (q > 0 && r != 0) {
  	    q = 1;
	    T(i, j) = q*fund_solution(k*r)*damping(damping_coef, r, k);
  	    ++nb;
  	  } else {
  	    T(i, j) = 0;
  	  }
  	  ++i;
  	}
  	++j;
      }
      svd[w] =  BDCSVD<MatrixXcf>(T,ComputeFullU | ComputeFullV);
      matrices[w] = T;
    }
  }

  GA_ROHandleI bs_es_handle(gdp->findAttribute(GA_ATTRIB_DETAIL, "buffer_size"));
  if (!bs_es_handle.isValid()) {
    addError(SOP_ATTRIBUTE_INVALID, "buffer sizes obsacle sources");
    return error();
  }
  if (bs_es_handle.get(0) != buffer_size) {
    addError(SOP_ATTRIBUTE_INVALID, "Buffer size of obstacle sources do not match buffer size of input sources. Did you forget to check the interactive sources box ?");
    return error();
  }

  const GA_Attribute *afs;
  const GA_AIFTuple *tuple; 

  // fill the p_in for each wavelength
  GA_ROHandleF abp_attrib(bp->findFloatTuple(GA_ATTRIB_POINT, "ampli", 2));
  if (abp_attrib.isValid()) {
    int winsize = bp->getPrimitiveRange().getEntries();
  }
  GA_ROHandleF a_handle(is->findFloatTuple(GA_ATTRIB_POINT, "ampli", buffer_size));
  for (int w = 0; w < nb_wl; ++w) {
    int as = ampli_steps[w];
    if (fr%as == 0) {
      float wl = wave_lengths[w];
      float k = M_PI*2.0/wl;
      float om = omega(k);
      if (wl == 0) {
	k = 0;
	om = 0;
      }
      int i = 0;
      const GA_Primitive* prim_bp = bp->getPrimitiveByIndex(w);
      GA_Range range_bp = prim_bp->getPointRange();
      int nb_bp = range_bp.getEntries();
      p_in[w] = VectorXcf(nb_bp);
      for(GA_Iterator itbp = range_bp.begin(); itbp != range_bp.end(); ++itbp) {
  	UT_Vector3 pos_b = bp->getPos3(*itbp);
  	p_in[w](i) = 0;
  	// add contribution from the spectrum computed at bp
  	if (!OBS(0) && abp_attrib.isValid()) {
  	  COMPLEX a(abp_attrib.get(*itbp, 0), abp_attrib.get(*itbp, 1));
  	  p_in[w](i) += a;
  	}
  	// add contribution from input sources
  	if (OBS(0)) {
	  const GA_Primitive* prim_is = is->getPrimitiveByIndex(w);
	  GA_Range range_is = prim_is->getPointRange();
	  for(GA_Iterator it = range_is.begin(); it != range_is.end(); ++it) {
	    UT_Vector3 pos_is = is->getPos3(*it);
	    float r = sqrt(pow(pos_b.x() - pos_is.x(), 2) + pow(pos_b.z() - pos_is.z(), 2));
	    float v = velocity(k);
	    fpreal t_ret = t - r/v;
	    OP_Context c_ret(t_ret);
	    int f_ret = floor((t_ret)/(dt*as))+1;
	    float ar = 0, ai = 0;
	    if (f_ret >= 1) {
	      afs = is->findFloatTuple(GA_ATTRIB_POINT, "ampli", buffer_size);
	      tuple = afs->getAIFTuple();
	      tuple->get(afs, *it, ar, 2*f_ret);
	      tuple->get(afs, *it, ai, 2*f_ret+1);
	    }
	    std::complex<float> ampli(ar, ai);
	    p_in[w](i) -= ampli*fund_solution(k*r)*damping(damping_coef, r, k);
	  }
  	}

  	if (getInputsArraySize() == 4) {
  	  const GU_Detail *obs = inputGeo(3);
  	  if (obs->getPrimitiveRange().getEntries() > w) {
  	    const GA_Primitive* prim_obs = obs->getPrimitiveByIndex(w);
  	    GA_Range range_obs = prim_obs->getPointRange();
  	    for(GA_Iterator it = range_obs.begin(); it != range_obs.end(); ++it) {
  	      UT_Vector3 pos_obs = obs->getPos3(*it);
  	      float r = sqrt(pow(pos_b.x() - pos_obs.x(), 2) + pow(pos_b.z() - pos_obs.z(), 2));
  	      float v = velocity(k);
  	      fpreal t_ret = t - r/v;
  	      OP_Context c_ret(t_ret);
  	      int f_ret = floor((t_ret)/(dt*as))+1;
  	      float ar = 0, ai = 0;
  	      if (f_ret >= 1) {
  		afs = obs->findFloatTuple(GA_ATTRIB_POINT, "ampli", buffer_size);
  		tuple = afs->getAIFTuple();
  		tuple->get(afs, *it, ar, 2*f_ret);
  		tuple->get(afs, *it, ai, 2*f_ret+1);
  	      }
  	      std::complex<float> ampli(ar, ai);
  	      p_in[w](i) -= ampli*fund_solution(k*r)*damping(damping_coef, r, k);
  	    }
  	  } else{
  	    addError(SOP_ATTRIBUTE_INVALID, "not the right nb of wavelength for entry 3");
  	  }
  	}
	
  	// add contribution from other sources of the obstacle
  	const GA_Primitive* prim = gdp->getPrimitiveByIndex(w);
  	GA_Range range = prim->getPointRange();
  	for(GA_Iterator it = range.begin(); it != range.end(); ++it) {
  	  UT_Vector3 pos_fs = gdp->getPos3(*it);
  	  float r = sqrt(pow(pos_b.x() - pos_fs.x(), 2) + pow(pos_b.z() - pos_fs.z(), 2));
  	  float ret = r/velocity(k);
  	  int f_ret = floor((t - ret)/(dt*as))+1;
  	  if (ret > t) {
  	    --f_ret;
  	  }
  	  float q = interpolation(ret,  shifts[w][i]*(dt*as), (shifts[w][i]+1)*dt*as);//0, (dt*as));
  	  float ar = 0, ai = 0;
  	  if (q <= 0) {
  	    if (f_ret >= 0) {
      	      afs = gdp->findFloatTuple(GA_ATTRIB_POINT, "ampli", buffer_size);
      	      tuple = afs->getAIFTuple();
      	      tuple->get(afs, *it, ar, 2*f_ret);
      	      tuple->get(afs, *it, ai, 2*f_ret+1);
      	    }
  	    COMPLEX a(ar, ai);
  	    p_in[w](i) -= a*fund_solution(k*r)*damping(damping_coef, r, k);
  	  } else {
  	  }
  	}
  	++i;
      }
    }
  }
  
  // solve for each wavelength and set new anpli
  GA_RWHandleF ampli_attrib(gdp->findFloatTuple(GA_ATTRIB_POINT, "ampli", buffer_size));
  if (!ampli_attrib.isValid()) {
    addError(SOP_ATTRIBUTE_INVALID, "ampli");
    return error();
  }
  for (int w = 0; w < nb_wl; ++w) {
    int as = ampli_steps[w];
    if (fr%as == 0) {
      const GA_Primitive* prim = gdp->getPrimitiveByIndex(w);
      GA_Range range = prim->getPointRange();
      int nb_es = range.getEntries();
      VectorXcf c(nb_es);
      c = svd[w].solve(p_in[w]);
      int i = 0;
      for(GA_Iterator it = range.begin(); it != range.end(); ++it, ++i) {
  	FLOAT prevr = 0, previ = 0;
  	if ((fr/as - shifts[w][i]) >= 0) {
	  FLOAT ar = real(c[i]);
	  FLOAT ai = imag(c[i]);
	  ampli_attrib.set(*it, 2*(fr/as - shifts[w][i]), ar);
	  ampli_attrib.set(*it, 2*(fr/as - shifts[w][i])+1, ai);
  	}
      }
      ampli_attrib.bumpDataId();
    }
  }
  time.tock(Times::solve_time_);
  std::cout<<"Solve time : "<<time.getTime(Times::solve_time_)<<std::endl; 
  time.next_loop();
  return error();
}

