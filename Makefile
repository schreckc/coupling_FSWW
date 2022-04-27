SRC=src/
HOUDINI_DIR=/opt/hfs17.0.459/
SHELL := /bin/bash
EIGEN=/usr/include/eigen3

all: inter freq obstacles boundary_points merge read_grid write_grid wave_param addplanar addfs

env:
	pushd $(HOUDINI_DIR);source houdini_setup;popd

inter: create_source solve_FS_inter deform_surface_inter
freq: create_source solve_FS deform_surface

obstacles: circle square point #texture

create_source: $(SRC)SOP_Create_Source.cpp $(SRC)SOP_Create_Source.hpp $(SRC)definitions.hpp
	hcustom $(SRC)SOP_Create_Source.cpp -I$(EIGEN) -g

create_sources: $(SRC)SOP_Create_Sources.cpp $(SRC)SOP_Create_Sources.hpp $(SRC)definitions.hpp
	hcustom $(SRC)SOP_Create_Sources.cpp -I$(EIGEN) -g

solve_FS_inter: $(SRC)SOP_Solve_FS_inter.cpp $(SRC)SOP_Solve_FS_inter.hpp $(SRC)definitions.hpp
	hcustom $(SRC)SOP_Solve_FS_inter.cpp -I$(EIGEN) -g

deform_surface_inter: $(SRC)SOP_Deform_Surface_inter.cpp $(SRC)definitions.hpp
	hcustom $(SRC)SOP_Deform_Surface_inter.cpp -I$(EIGEN) -g

solve_FS: $(SRC)SOP_Solve_FS.cpp $(SRC)SOP_Solve_FS.hpp $(SRC)definitions.hpp
	hcustom $(SRC)SOP_Solve_FS.cpp -I$(EIGEN) -g

deform_surface: $(SRC)SOP_Deform_Surface.cpp $(SRC)SOP_Deform_Surface.cpp $(SRC)definitions.hpp
	hcustom $(SRC)SOP_Deform_Surface.cpp -I$(EIGEN) -g

circle: $(SRC)SOP_CircleObstacle_Src.cpp $(SRC)SOP_CircleObstacle_Src.hpp $(SRC)definitions.hpp
	hcustom $(SRC)SOP_CircleObstacle_Src.cpp -I$(EIGEN) -g

square: $(SRC)SOP_SquareObstacle_Src.cpp $(SRC)SOP_SquareObstacle_Src.hpp $(SRC)definitions.hpp
	hcustom $(SRC)SOP_SquareObstacle_Src.cpp -I$(EIGEN) -g

point: $(SRC)SOP_PointObstacle_Src.cpp $(SRC)SOP_PointObstacle_Src.hpp $(SRC)definitions.hpp
	hcustom $(SRC)SOP_PointObstacle_Src.cpp -I$(EIGEN) -g

merge: $(SRC)SOP_Merge_Sources.cpp $(SRC)SOP_Merge_Sources.hpp $(SRC)definitions.hpp
	hcustom $(SRC)SOP_Merge_Sources.cpp -I$(EIGEN) -g

boundary_points: $(SRC)SOP_Boundary_Points.cpp $(SRC)SOP_Boundary_Points.hpp $(SRC)InputPoint.hpp $(SRC)FFT.hpp $(SRC)definitions.hpp
	hcustom $(SRC)SOP_Boundary_Points.cpp -I$(EIGEN) -g

boundary_points_test: $(SRC)SOP_Boundary_Points_test.cpp $(SRC)SOP_Boundary_Points_test.hpp $(SRC)InputPoint.hpp $(SRC)FFT.hpp $(SRC)definitions.hpp
	hcustom $(SRC)SOP_Boundary_Points_test.cpp -I$(EIGEN) -g

recording: $(SRC)SOP_Recording.cpp $(SRC)SOP_Recording.hpp $(SRC)SOP_Recording_one.cpp $(SRC)SOP_Recording_one.hpp $(SRC)InputPoint.hpp $(SRC)FFT.hpp $(SRC)definitions.hpp
	hcustom $(SRC)SOP_Recording.cpp -I$(EIGEN) -g
	hcustom $(SRC)SOP_Recording_one.cpp -I$(EIGEN) -g

write_grid: $(SRC)SOP_WriteGrid.cpp $(SRC)SOP_WriteGrid.hpp $(SRC)definitions.hpp
	hcustom $(SRC)SOP_WriteGrid.cpp -I$(EIGEN) -g

read_grid: $(SRC)SOP_ReadGrid.cpp $(SRC)SOP_ReadGrid.hpp $(SRC)definitions.hpp
	hcustom $(SRC)SOP_ReadGrid.cpp -I$(EIGEN) -g

noise: $(SRC)SOP_addNoise.cpp $(SRC)SOP_addNoise.hpp $(SRC)definitions.hpp $(SRC)SimplexNoise.h
	hcustom $(SRC)SOP_addNoise.cpp -I$(EIGEN) -g

addfs: $(SRC)SIM_addFS.cpp $(SRC)SIM_addFS.hpp $(SRC)definitions.hpp
	hcustom $(SRC)SIM_addFS.cpp -I$(EIGEN) -g

addplanar: $(SRC)SIM_addPlanarWave.cpp $(SRC)SIM_addPlanarWave.hpp $(SRC)SOP_PlanarWave.cpp $(SRC)SOP_PlanarWave.hpp $(SRC)definitions.hpp
	hcustom $(SRC)SIM_addPlanarWave.cpp -I$(EIGEN) -g
	hcustom $(SRC)SOP_PlanarWave.cpp -I$(EIGEN) -g	
wave_param: $(SRC)SOP_WaveParameters.cpp $(SRC)SOP_WaveParameters.hpp $(SRC)definitions.hpp
	hcustom $(SRC)SOP_WaveParameters.cpp -I$(EIGEN) -g

clean:
	rm ($SRC)*.o
