Coupling 3D Liquid Simulation with 2D Wave Propagation for Large Scale Water Surface Animation Using the Equivalent Sources Method
Camille Schreck 1,2 and Chris Wojtan 1
2
1 IST Austria
Université de Lorraine, CNRS, Inria, LORIA, France
-----------------------------------------------------------------------------

This code implements a simplified version the method described in

    "Coupling 3D Liquid Simulation with 2D Wave Propagation for Large Scale Water Surface Animation Using the Equivalent Sources Method"
    by Camille Schreck and Chris Wojtan
    published in Computer Graphics Forum (Eurographics 2022)

This code is not particularly well-tested or optimized, but it can be used to
reproduce (some) of the examples in the paper. This is not a GPU or multithreaded version.
The code depends on Houdini. We have only tested the code using Houdini FX 18.5 under
Ubuntu 20.10.

To use the code, you first have to compile the plug-in. A Makefile is provided
in the root directory. Before you can compile the code, you need to
load in the HDK environment variables.

    $ pushd /opt/hfs18.5
    $ source houdini_setup
    $ popd

Then you can compile the code.

    $ make

Once the plug-in is installed in your home directory (there should be files like "SOP_Create_Source.so" in
~/houdini18.5/dso/) you can load up the included Houdini scene
file.

    $ houdini extended_waves.hipnc

------------------------------------------------------------------------------
Camille Schreck <camille.schreck@inria.fr
Last updated 2022-05-03

--------------------------------------------------------------------------------
Implementation
--------------------------------------------------------------------------------
If you want to understand or modify the code, first thing to do is to read the article.
Here are the details specific to the houdini implementation.

NOTE: the houdini name of any the node created by this code contains "FS" so you can find them
  more easily in Houdini.

A set of sources is encoded as a geometry detail:
  * Detail:
      -Attibutes:
         buffer_size (int): (aperiodic version) size of the buffer recording past amplitude
	                          for each source (number of float, as amplitudes are complexes
				  buffer_size/2 amplitudes are recorded)
			    (periodic version)	should be 2.  
	 damping (float)
	 ampli
	 inter (1 if interactive sources, 0 if not)
	 shiftb, shiftu, windowsize
  *Primitive: subset of source of same wavelength
      -Attibutes:
         wavelengths (float)
	 ampli_step (int): amplitudes of the sources updated every <ampli_step> time step
  *Point: one source of one wavelength
      -Attibutes:
        P (UT_Vector3): position, (x, z) should be the position of the surface, y=0
	ampli (FloatTuple of size <buffer_size>): complex amplitudes of the source
	           ampli[2n], ampli[2n+1]: real and imaginay part of the amplitude after n time step
(Note that the wavelength of one source point is defined be their corresponding primitive)

Representation of the surface:
We use a regular or projected grid, but any planar surface mesh should work.
For the periodic method, each node of the grid store a spectrum (a complex amplitude for each 
  wavelength).

SOP_WaveParameters:
Define the parameters of the sources
  Parameters:
     -Amplitude: multiply the amplitude of every sources by this number
     If defined by spectrum:
         - use the wavelengths obtained by doing an FFT on a time window of size <window size> (time step 0.3) using the dispertion relationship defined in definition.hpp. The parameters shiftB and shiftU remove respectively the <shifb> smaller wavelength and <shiftu> bigger ones.
    If not defined by spectrum  
         -Minimum/maximum wavelength, wavelength multiplicative step: used to compute the range of wl as a geometric sequence
     -Interactive sources
     -Size of the buffer: size of the buffer recording past amplitude
     -Damping

definition.hpp:
cointains notably the definition of the damping, dispersion relationship and fundamental solution

SOP_CircleObstacle_Src (Square/Texture -> need update):
Create set of point sources along an offset surface of the obstacle.
Need the waveParameters in entry.
Create a subset of sources for each wl.
Spacing depends on the wavelength and the parameter "density".
Note: the number of boundary points should be at least twice the number of points for the biggest subset of sources (the one corresponding to the smallest wavelength).

SOP_Boundary_Points:
Record the height at the position of the boundary points (input 0) to compute the amplitude of each frequency component using the wave parameters (input 1).

SOP_Solve_FS_inter:
Compute amplitude of the obstacle sources (input 0) such that the boundary conditions computed at the boundary points (input 1) given the wave parameters (input 2)

SOP_Deform_Surface_inter:
Compute height at each point of the grid (input 0) by summing all the sources.
(each other input is a set of sources)