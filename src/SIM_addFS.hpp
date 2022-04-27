/*
 * Copyright (c) 2019
 *      Side Effects Software Inc.  All rights reserved.
 *
 * Redistribution and use of Houdini Development Kit samples in source and
 * binary forms, with or without modification, are permitted provided that the
 * following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. The name of Side Effects Software may not be used to endorse or
 *    promote products derived from this software without specific prior
 *    written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SIDE EFFECTS SOFTWARE `AS IS' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN
 * NO EVENT SHALL SIDE EFFECTS SOFTWARE BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *----------------------------------------------------------------------------
 */
#ifndef ADD_FS_HPP
#define ADD_FS_HPP
#include <GAS/GAS_SubSolver.h>
#include <GAS/GAS_Utils.h>
#include <SIM/SIM_VectorField.h>
#include <SIM/SIM_ScalarField.h>

namespace HDK_Sample {
/// A simple field manipulation class that will add fields
/// together.  
class SIM_addFS : public GAS_SubSolver
{
public:
    /// These macros are used to create the accessors
    /// getFieldDstName and getFieldSrcName functions we'll use
    /// to access our data options.
    GET_DATA_FUNC_S(GAS_NAME_GEOMETRY, Geometry);
    static int fr;

protected:
    explicit             SIM_addFS(const SIM_DataFactory *factory);
    virtual             ~SIM_addFS();
    /// Used to determine if the field is complicated enough to justify
    /// the overhead of multithreading.
    bool                 shouldMultiThread(const SIM_RawField *field) const 
                         { return field->field()->numTiles() > 10; }
    /// The overloaded callback that GAS_SubSolver will invoke to
    /// perform our actual computation.  We are giving a single object
    /// at a time to work on.
    virtual bool         solveGasSubclass(SIM_Engine &engine,
                                SIM_Object *obj,
                                SIM_Time time,
                                SIM_Time timestep);
  THREADED_METHOD4(SIM_addFS, /*shouldMultiThread(vel->getField(0))*/true,
		   addContribFromPrimitive,
		   const GU_Detail &, fs,
		   GA_Offset, prim_off,
		   SIM_VectorField *, vel,
		   SIM_ScalarField *, surface);
  THREADED_METHOD4(SIM_addFS, /*shouldMultiThread(vel->getField(0))*/true,
		   splitVel,
		   SIM_VectorField *, vel,
		   SIM_VectorField *, velx,
		   SIM_VectorField *, vely,
		   SIM_VectorField *, velz);
  void addContribFromPrimitivePartial(const GU_Detail &fs,
				      GA_Offset prim_off,
				      SIM_VectorField *vel,
				      //SIM_RawField *vel,
				      SIM_ScalarField *surface,
				      const UT_JobInfo &info);
void splitVelPartial(SIM_VectorField *vel,
		     SIM_VectorField *velx,
		     SIM_VectorField *vely,
		     SIM_VectorField *velz,
		     const UT_JobInfo &info);

private:
    /// We define this to be a DOP_Auto node which means we do not
    /// need to implement a DOP_Node derivative for this data.  Instead,
    /// this description is used to define the interface.
    static const SIM_DopDescription     *getDopDescription();
    /// These macros are necessary to bind our node to the factory and
    /// ensure useful constants like BaseClass are defined.
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(SIM_addFS,
                        GAS_SubSolver,
                        "add FS to field",
                        getDopDescription());
  int buffer_size;
  float damping_coef;
};
  
} // End HDK_Sample namespace
#endif
