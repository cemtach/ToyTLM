//
// ToyTLM.c
//
// Copyright (C) 1999 Paul Hayes, Matthew O'Keefe 
//
// This   program is free   software;  you can redistribute it  and/or
// modify it  under the terms   of the GNU  General  Public License as
// published by the Free Software Foundation;  either version 2 of the
// License,  or any   later version,  with   the  following conditions
// attached in addition to any and all conditions of the GNU
//
// General Public License: When reporting or displaying any results or
// animations created using this code or modification of this code,
// make the appropriate citation referencing ToyFDTD by name and
// including the version number.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
//
// Contacting the authors:
//
// Paul Hayes, Matthew O'Keefe
// Department of Electrical and Computer Engineering
//      200 Union Street S. E.
//      Minneapolis, MN 55455
//
// info@cemtach.com
// 
// http://cemtach.com
//
// This code is here for everyone, but not everyone will need
// something so simple, and not everyone will need to read all the
// comments.
//
// This ToyTLM is a stripped-down, minimalist 3D TLM code.  It
// illustrates minimum factors that must be considered to create a
// simple TLM simulation.
//
#include <math.h> 
#include <stdio.h>         
#include <float.h>
//
// malloc may be defined locally in one of these locations.. Should
// either cause compile difficulties, simply comment out.  Try [man
// malloc] which should indicate the proper header file location.
#include <stdlib.h>
#include <malloc.h>

//
// Total number of timesteps to be computed
//
#define MAXIMUM_ITERATION 1003
//
// The program will output 3D data every PLOT_MODULUS timesteps,
// except the last iteration computed is always output.  So if
// MAXIMUM_ITERATION is not an integer multiple of PLOT_MODULUS, the
// last timestep output will come after a shorter interval than that
// separating previous outputs.
//
#define PLOT_MODULUS 5
//
// Frequency of the stimulus in Hertz
//
#define FREQUENCY 10.0e9     
//
// Waveguide width in meters
//
#define GUIDE_WIDTH 0.0229
//
// Waveguide height in meters
//
#define GUIDE_HEIGHT 0.0102
//
// Length of the waveguide in wavelengths of the stimulus wave
//
#define LENGTH_IN_WAVELENGTHS 5.0
//
// Minimum number of grid cells per wavelength in the x, y, and z directions
//
#define CELLS_PER_WAVELENGTH 25.0


//
// Speed of light in a vacuum in meters/second
//
#define LIGHT_SPEED 299792458.0       
//
// Permeability of free space in henry/meter
//
#define MU_0 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6
//
// Permittivity of free space in farad/meter
//
#define EPSILON_0 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12
//
// Impedance of freespace or sqrt(MU_0/EPS_0)
//
#define Z_0 376.734309182110149436898630439319


typedef struct tlmVoxelStruct TlmVoxel;

struct tlmVoxelStruct
{
  //
  // First 12 voltages essentially model the electric and magnetic
  // fields within the voxel.
  //
  float v1;
  float v2;
  float v3;
  float v4;
  float v5;
  float v6;
  float v7;
  float v8;
  float v9;
  float v10;
  float v11;
  float v12;
  //
  // Last 6 voltages essentially model spatial distortions of the cell
  // or material parameters within the cell.
  //
  float v13;
  float v14;
  float v15;
  float v16;
  float v17;
  float v18;
};

// Function Prototypes;
int main(void);


int
main()	
{
  //
  // Commonly utilized integer access into arrays
  //
  int i,j,k;
  //
  // Total number of cells along the x, y, and z axes, respectively
  //
  int nx, ny, nz;
  //
  // Counter to keep track of dynamically allocated memory, ordinarily this
  // is the largest memory requirement and not the variables declared here.
  //
  int allocatedBytes = 0;
  //
  // The current iteration along in time signifying a delta in time
  //
  int iteration = 0;
  //
  // Value of the excitation at a given time
  //
  float stimulus = 0.0;
  //
  // The current simulation time in seconds
  //
  float currentSimulatedTime = 0.0;
  //
  // The total simulation time the simulation should proceed towards in seconds
  //
  float totalSimulatedTime = 0.0;
  //
  // The radian frequency and wavelength in freespace of the excitation
  //
  float omega;
  float lambda0;
  //
  // Spatial deltas in each cartesian coordinate direction
  //
  float dx, dy, dz; 
  //
  // The time delta
  //
  float dt;
  //
  // Voltage array pointer to contain all 18 base voltages which represent
  // the electric and magnetic fields within the voxel.
  //
  TlmVoxel ***voxel;
  //
  // A filename variable for flexible naming with 1024 chosen as the 
  // maximum path possible under POSIX.
  //
  char filename[1024];
  //
  // Autoscaling values for the min and max during each output
  //
  float min, max;
  //
  // Norm is the largest absolute value of either min or max.  This yields
  // a scaling centered with the largest values at the extreme of the scaling.
  // in other words, zero should always be at the center of the scale.
  //
  float norm;
  //
  // The actual scaling multiplier once the norm, min and max have been addressed.
  //
  float scalingValue;
  //
  // File access variables.
  //
  FILE *openFilePointer;
  FILE *vizFilePointer;
  //
  // Stub quantities which may be commonly used in this simple program
  // as simple common variables.
  //
  float yx,yy,yz,zx,zy,zz,g0;
  float zx_by_2,zy_by_2,zz_by_2;
  float yx_by_2,yy_by_2,yz_by_2;
  //
  // Commonly utilized temporary variables
  //
  float rtmp1, rtmp2; // , rtmp3, rtmp4; unused.
  float ex,ey,ez,hx,hy,hz;
  float v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18;
  float mur, epsr;
  float delt1, delt2, delt3, delt4, delt5, delt6;
  //
  // Mesh scaling values calculated from dx, dy and dz
  //
  float u,v,w; 

  //
  //
  // David K. Cheng, Field and Wave Electromagnetics, 2nd ed., pages
  // 554-555.  Rectangular waveguide, interior width = 2.29cm,
  // interior height = 1.02cm.  This is a WG-16 waveguide useful for
  // X-band applications.
  //
  // There should be at least 20 cells per wavelength in each
  // direction, but we'll go with 25 so the animation will look
  // prettier.  (CELLS_PER_WAVELENGTH was set to 25.0 in the global
  // constants at the beginning of the code.)
  //
  // The number of cells along the width of the guide and the width of
  // those cells should fit the guide width exactly, so that ny*dy =
  // GUIDE_WIDTH meters.  The same should be true for nz*dz =
  // GUIDE_HEIGHT meters.
  //
  // dx is chosen to be dy or dz -- whichever is smaller
  //
  // nx is chosen to make the guide at least LENGTH_IN_WAVELENGTHS
  // wavelengths long.
  // 
  // dt is chosen for Courant stability; the time step must be kept
  // small enough so that the plane wave only travels one cell length
  // (one dx) in a single timestep.  Otherwise FDTD cannot keep up
  // with the signal propagation, since FDTD computes a cell only from
  // it's immediate neighbors.
  //
  // Calculate the freespace wavelength and radian frequency of the excitation.
  //
  lambda0 = LIGHT_SPEED/FREQUENCY;
  omega = 2.0*M_PI*FREQUENCY; 
  //
  // Starting with a ridiculously small ny
  //
  ny = 3;
  //
  // Estimate dy from the guide width and ny.
  //
  dy = GUIDE_WIDTH/ny;
  //
  // Continue in this vein until the dy spatial delta is less than or equal to
  // lambda0/CELLS_PER_WAVELENGTH.
  //
  while(dy >= lambda0/CELLS_PER_WAVELENGTH)
    {
    ny++;
    dy = GUIDE_WIDTH/ny;
    }
  //
  // Starting with a ridiculously small nz
  //
  nz = 3;  
  //
  // Estimate dz from the guide width and nz.
  //
  dz = GUIDE_HEIGHT/nz;
  //
  // Continue in this vein until the dz spatial delta is less than or equal to
  // lambda0/CELLS_PER_WAVELENGTH.
  //
  while(dz >= lambda0/CELLS_PER_WAVELENGTH)
    {
    nz++;
    dz = GUIDE_HEIGHT/nz;
    }
  //
  // Now that dy and dz have been established, the dx will be established
  // to maintain at least as good spatial properties as the two already chosen.
  // One possible way is to force the dx to simply be the smallest of dy and dz.
  //
  dx = (dy < dz) ? dy : dz;
  //
  // The chosen dx effectively sets the nx or number of cells in the x direction.
  //
  nx = (int)(LENGTH_IN_WAVELENGTHS*lambda0/dx);


  //
  // Since dx is the smallest spatial delta, the grading value in the x
  // direction, u, becomes 1.0.  All others are then based off those scalings
  // of the spatial deltas.  These are relative delta values which scale the 
  // stub quantities and essentially distort the voxel to the requested size.
  //
  u = 1.0;
  v = dy/dx;
  w = dz/dx;

  //
  // Allocate memory for the voltage arrays which represent the
  // electric and magnetic fields within a voxel.  Clear the arrays
  // to force an empty field system.
  //
  voxel = (TlmVoxel ***)malloc(nx*sizeof(TlmVoxel **));
  for(i=0; i<nx; i++)
    {
    voxel[i] = (TlmVoxel **)malloc(ny*sizeof(TlmVoxel *));
    for(j=0; j<ny; j++)
      {
      voxel[i][j] = (TlmVoxel *)malloc(nz*sizeof(TlmVoxel));
      for(k=0; k<nz; k++)
	{
	voxel[i][j][k].v1 = 0.0;
	voxel[i][j][k].v2 = 0.0;
	voxel[i][j][k].v3 = 0.0;
	//
	voxel[i][j][k].v4 = 0.0;
	voxel[i][j][k].v5 = 0.0;
	voxel[i][j][k].v6 = 0.0;
	//
	voxel[i][j][k].v7 = 0.0;
	voxel[i][j][k].v8 = 0.0;
	voxel[i][j][k].v9 = 0.0;
	//
	voxel[i][j][k].v10 = 0.0;
	voxel[i][j][k].v11 = 0.0;
	voxel[i][j][k].v12 = 0.0;
	//
	voxel[i][j][k].v13 = 0.0;
	voxel[i][j][k].v14 = 0.0;
	voxel[i][j][k].v15 = 0.0;
	//
	voxel[i][j][k].v16 = 0.0;
	voxel[i][j][k].v17 = 0.0;
	voxel[i][j][k].v18 = 0.0;
	}
      }
    }
  //
  // This is a rough estimate and does not include the actual pointer
  // values, yet, most of the memory is included in the TlmVoxel
  // mallocs.
  //
  allocatedBytes += nx*ny*nz*sizeof(TlmVoxel);
  //
  // Calculate the time delta based on the size of the actual voxels
  // throughout the entire mesh which may be distorted in size or have
  // various materials.
  //
  rtmp1 = FLT_MAX;
  for(k=0; k<nz; k++)
    {
    for(j=0; j<ny; j++)
      {
      for(i=0; i<nx; i++)
	{
	mur = 1.0;
	epsr = 1.0;
	delt1=v*w*epsr/(u*2.0*LIGHT_SPEED);
	delt2=u*w*epsr/(v*2.0*LIGHT_SPEED);
	delt3=v*u*epsr/(w*2.0*LIGHT_SPEED);
	delt4=v*w*mur/(u*2.0*LIGHT_SPEED);
	delt5=u*w*mur/(v*2.0*LIGHT_SPEED);
	delt6=v*u*mur/(w*2.0*LIGHT_SPEED);
	rtmp1 = (rtmp1 < delt1) ? rtmp1:delt1;   // pick smallest of rtmp1 and delt1
	rtmp1 = (rtmp1 < delt2) ? rtmp1:delt2;   // pick smallest of rtmp1 and delt2
	rtmp1 = (rtmp1 < delt3) ? rtmp1:delt3;   // pick smallest of rtmp1 and delt3
	rtmp1 = (rtmp1 < delt4) ? rtmp1:delt4;   // pick smallest of rtmp1 and delt4
	rtmp1 = (rtmp1 < delt5) ? rtmp1:delt5;   // pick smallest of rtmp1 and delt5
	rtmp1 = (rtmp1 < delt6) ? rtmp1:delt6;   // pick smallest of rtmp1 and delt6
	}
      }
    }
  //
  // Establish the time delta based on the previous estimate code
  //
  dt = rtmp1;
  //
  // Base the total simulation time off the requested number of
  // iterations and time delta
  //
  totalSimulatedTime = MAXIMUM_ITERATION*dt*dx;




  //
  // Calculate the actual stub values based on freespace and
  // precalculate constants which will be used in the update
  // equations.
  //
  mur = 1.0;
  epsr = 1.0;
  rtmp1 = epsr/(dt*LIGHT_SPEED);
  rtmp2 = mur/(dt*LIGHT_SPEED);
  yx = 2.0*(v*w*rtmp1/u - 2.0);
  yy = 2.0*(u*w*rtmp1/v - 2.0);
  yz = 2.0*(v*u*rtmp1/w - 2.0);
  zx = 2.0*(v*w*rtmp2/u - 2.0);
  zy = 2.0*(u*w*rtmp2/v - 2.0);
  zz = 2.0*(v*u*rtmp2/w - 2.0);
  //
  yx_by_2 = yx/2.0;
  yy_by_2 = yy/2.0;
  yz_by_2 = yz/2.0;
  zx_by_2 = zx/2.0;
  zy_by_2 = zy/2.0;
  zz_by_2 = zz/2.0;
  g0 = 0.0; 
  

  //
  // Output some of the simulation information to the user which is
  // particularly useful after the simulation is done.
  //
  fprintf(stdout, "\n"); 
  fprintf(stdout, "bob -cmap gbry.cmap -s %dx%dx%d *.bob\n", nx, ny, nz); 
  fprintf(stdout, "\n"); 
  fprintf(stdout, "viz ToyTLMc.viz\n"); 
  fprintf(stdout, "\n"); 
  fprintf(stdout, "Dynamically allocated %d bytes\n", allocatedBytes); 
  fprintf(stdout, "\n"); 
  fprintf(stdout, "Meshing parameters:\n"); 
  fprintf(stdout, "%dx%dx%d cells\n", nx, ny, nz); 
  fprintf(stdout, "dx=%lg, dy=%lg, dz=%lg meters\n", dx, dy, dz); 
  fprintf(stdout, "u=%lg, v=%lg, w=%lg\n", u, v, w); 
  fprintf(stdout, "%lg x %lg x %lg meter^3 simulation region\n",  
	  GUIDE_WIDTH, GUIDE_HEIGHT, LENGTH_IN_WAVELENGTHS*lambda0); 
  fprintf(stdout,"dt=%g \n",dt*dx);
  fprintf(stdout, "\n"); 




  //
  // Open and start writing the .viz file with header information
  //
  while ((vizFilePointer = fopen("ToyTLMc.viz", "w")) == NULL)
    {
    fprintf(stderr, "Difficulty opening ToyTLMc.viz");
    perror(" ");
    }
  fprintf(vizFilePointer, "#Viz V1.0\n");
  fprintf(vizFilePointer, "time: %lg %lg\n", currentSimulatedTime, dt*dx);
  fprintf(vizFilePointer, "color: gbry.cmap\n");
  fprintf(vizFilePointer, "\n");


  for(iteration = 0; iteration < MAXIMUM_ITERATION; iteration++)
    {
    //
    // Time in simulated seconds that the simulation has progressed.
    //
    currentSimulatedTime = dt*dx*(float)iteration;  
    //
    // Print to standard output the iteration number and current simulated time.
    //
    fprintf(stdout, "#%d %lgsec", iteration, currentSimulatedTime);
    //
    // 3D data output every PLOT_MODULUS timesteps:
    //     The first time through the main loop all the data written to 
    //     file will be zeros.  If anything is nonzero, there's a bug.  :>
    //
    if ( (iteration % PLOT_MODULUS) == 0)
      {
      //
      // Create the filename for this iteration, which includes the iteration number.
      //
      sprintf(filename, "c_%06d.bob", iteration);
      //
      // open a new data file for this iteration:
      //
      while ((openFilePointer = fopen(filename, "wb")) == NULL)
	{
	//
	// If unable to open the file, exit with a descriptive failure.
	//
	fprintf(stderr, "Difficulty opening c_%06d.bob", iteration);
	perror(" ");
	}
      //
      // Locate the min and max values in order to autoscale
      //
      min = FLT_MAX;
      max = -FLT_MAX;
      for(k=0;k<nz;k++)
	{
	for(j=0;j<ny;j++)
	  {
	  for(i=0;i<nx;i++)
	    {
	    //
//  	    ex=2.0*(voxel[i][j][k].v1 + voxel[i][j][k].v2 +
// 	 	    voxel[i][j][k].v9 + voxel[i][j][k].v12 +
// 		    yx*voxel[i][j][k].v13) / (u*(4.0+yx));
// 	    //
//  	    ey=2.0*(voxel[i][j][k].v3 + voxel[i][j][k].v4 +
// 		    voxel[i][j][k].v8 + voxel[i][j][k].v11 +
// 		    yy*voxel[i][j][k].v14) / (v*(4.0+yy));
	    //
	    ez=2.0*(voxel[i][j][k].v5 + voxel[i][j][k].v6 +
		    voxel[i][j][k].v7 + voxel[i][j][k].v10 +
		    yz*voxel[i][j][k].v15) / (w*(4.0+yz));
	    //
//  	    hx=-2.0*(voxel[i][j][k].v4 - voxel[i][j][k].v5 +
// 		     voxel[i][j][k].v7 - voxel[i][j][k].v8 -
// 		     voxel[i][j][k].v16) / (Z_0*u*(4.0+zx));
// 	    //
//  	    hy=-2.0*(-voxel[i][j][k].v2 + voxel[i][j][k].v6 +
// 		     voxel[i][j][k].v9 - voxel[i][j][k].v10 -
// 		     voxel[i][j][k].v17)/(Z_0*v*(4.0+zy));
// 	    //
//  	    hz=-2.0*(-voxel[i][j][k].v3 + voxel[i][j][k].v1 +
// 		     voxel[i][j][k].v11 - voxel[i][j][k].v12 -
// 		     voxel[i][j][k].v18)/(Z_0*w*(4.0+zz));
	    //
	    //
	    //
	    if (ez < min)
	      {
	      min = ez;
	      }
	    if (ez > max)
	      {
	      max = ez;
	      }
	    }
	  }
	}
      //
      // Set norm to be max or min, whichever is greater in magnitude.
      //
      norm = (fabs(max) > fabs(min)) ? fabs(max) : fabs(min);
      if (norm == 0.0)
	{
	//
	// If everything is zero, give norm a tiny value to avoid division by zero.
	//
	norm = DBL_EPSILON;
	}
      scalingValue = 127.0/norm;
      //
      // Write to standard output the minimum and maximum values from this iteration
      // and the minimum and maximum values that will be written to the bob file this iteration.
      //
      fprintf(stdout, "\t%lg(%d) < ez BoB < %lg(%d)",
	      min, (int)(128.0 + scalingValue*min),
	      max, (int)(128.0 + scalingValue*max));
      //
      // Scale each ez value in the mesh to the range of integers from zero through 254
      // and write them to the output file for this iteration.
      //
      for(k=0;k<nz;k++)
	{
	for(j=0;j<ny;j++)
	  {
	  for(i=0;i<nx;i++)
	    {
	    //
//  	    ex=2.0*(voxel[i][j][k].v1 + voxel[i][j][k].v2 +
// 	 	    voxel[i][j][k].v9 + voxel[i][j][k].v12 +
// 		    yx*voxel[i][j][k].v13) / (u*(4.0+yx));
// 	    //
//  	    ey=2.0*(voxel[i][j][k].v3 + voxel[i][j][k].v4 +
// 		    voxel[i][j][k].v8 + voxel[i][j][k].v11 +
// 		    yy*voxel[i][j][k].v14) / (v*(4.0+yy));
	    //
	    ez=2.0*(voxel[i][j][k].v5 + voxel[i][j][k].v6 +
		    voxel[i][j][k].v7 + voxel[i][j][k].v10 +
		    yz*voxel[i][j][k].v15) / (w*(4.0+yz));
	    //
//  	    hx=-2.0*(voxel[i][j][k].v4 - voxel[i][j][k].v5 +
// 		     voxel[i][j][k].v7 - voxel[i][j][k].v8 -
// 		     voxel[i][j][k].v16) / (Z_0*u*(4.0+zx));
// 	    //
//  	    hy=-2.0*(-voxel[i][j][k].v2 + voxel[i][j][k].v6 +
// 		     voxel[i][j][k].v9 - voxel[i][j][k].v10 -
// 		     voxel[i][j][k].v17)/(Z_0*v*(4.0+zy));
// 	    //
//  	    hz=-2.0*(-voxel[i][j][k].v3 + voxel[i][j][k].v1 +
// 		     voxel[i][j][k].v11 - voxel[i][j][k].v12 -
// 		     voxel[i][j][k].v18)/(Z_0*w*(4.0+zz));
	    //
	    // Put the value out to the file
	    //
	    putc((int)(128.0 + scalingValue*ez), openFilePointer);
	    }
	  }
	}
 
      //
      // Close the output file for this iteration.
      //
      fclose(openFilePointer);

      //
      // Write the dimensions and name of the output file for this
      // iteration to the viz control.
      //
      fprintf(vizFilePointer, "%dx%dx%d %s\n", nx+1, ny+1, nz, filename);
      //
      // Write identification of the corners of the mesh and the max and
      // min values for this iteration to the viz control file.
      //
      fprintf(vizFilePointer, "bbox: 0.0 0.0 0.0 %lg %lg %lg %lg %lg\n",
	      dx*(double)nx, dy*(double)ny, dz*(double)nz, min, max);

      }
   

    //
    // Compute the stimulus: a plane wave emanates from the x=0 face:
    // The length of the guide lies in the x-direction, the width of
    // the guide lies in the y-direction, and the height of the guide
    // lies in the z-direction.  So the guide is sourced by all the ez
    // components on the stimulus face.
    //
    stimulus = sin(omega*currentSimulatedTime);
    for (i=0; i<1; i++)
      { 
      for(j=0; j<ny; j++)
	{
	for(k=0; k<nz; k++)
	  {
	  //
	  // Calculate the actual fields to excite
	  //
	  ex = 0.0;
	  ey = 0.0;
 	  ez = stimulus;
	  hx = 0.0;
	  hy = 0.0;
	  hz = 0.0;

	  //
	  // Map those excited fields into the voltages utilized in
	  // the TLM method.
	  //
	  voxel[i][j][k].v1  = +(u*ex + w*Z_0*hz)/2.0;
	  voxel[i][j][k].v2  = +(u*ex + v*Z_0*hy)/2.0;
	  voxel[i][j][k].v3  = +(v*ey - w*Z_0*hz)/2.0;
	  //
	  voxel[i][j][k].v4  = +(v*ey + u*Z_0*hx)/2.0;
	  voxel[i][j][k].v5  = +(w*ez - u*Z_0*hx)/2.0;
	  voxel[i][j][k].v6  = +(w*ez + v*Z_0*hy)/2.0;
	  //
	  voxel[i][j][k].v7  = +(w*ez + u*Z_0*hx)/2.0;
	  voxel[i][j][k].v8  = +(v*ey - u*Z_0*hx)/2.0;
	  voxel[i][j][k].v9  = +(u*ex + v*Z_0*hy)/2.0;
	  //
	  voxel[i][j][k].v10 = +(w*ez - v*Z_0*hy)/2.0;
	  voxel[i][j][k].v11 = +(v*ey + w*Z_0*hz)/2.0;
	  voxel[i][j][k].v12 = +(u*ex - w*Z_0*hz)/2.0;
	  //
	  voxel[i][j][k].v13 = +u*ex/2.0;
	  voxel[i][j][k].v14 = +v*ey/2.0;
	  voxel[i][j][k].v15 = +w*ez/2.0;
	  //
	  voxel[i][j][k].v16 = -zx*Z_0*u*hx/2.0;
	  voxel[i][j][k].v17 = -zy*Z_0*v*hy/2.0;
	  voxel[i][j][k].v18 = -zz*Z_0*w*hz/2.0;
	  }
	}
      }



    //
    // Scattering which accounts for the propagation of the electric
    // and magnetic fields within a voxel.  This section would roughly
    // correspond to the field update equations in the FDTD method.
    //
    for (k=0; k<nz; k++)
      {
      for (j=0; j<ny; j++)
	{
	for (i=0; i<nx; i++)
	  {
	  //
	  // Store the voltages in temporary variables since the
	  // update equations will modify the values and thus skew
	  // results.
	  //
	  v1 = voxel[i][j][k].v1;
	  v2 = voxel[i][j][k].v2;
	  v3 = voxel[i][j][k].v3;
	  //                
	  v4 = voxel[i][j][k].v4;
	  v5 = voxel[i][j][k].v5;
	  v6 = voxel[i][j][k].v6;
	  //                
	  v7 = voxel[i][j][k].v7;
	  v8 = voxel[i][j][k].v8;
	  v9 = voxel[i][j][k].v9;
	  //                
	  v10 = voxel[i][j][k].v10;
	  v11 = voxel[i][j][k].v11;
	  v12 = voxel[i][j][k].v12;
	  //                
	  v13 = voxel[i][j][k].v13;
	  v14 = voxel[i][j][k].v14;
	  v15 = voxel[i][j][k].v15;
	  //                
	  v16 = voxel[i][j][k].v16;
	  v17 = voxel[i][j][k].v17;
	  v18 = voxel[i][j][k].v18;

	  //
	  // Upgrade the voltages representing the electric and magnetic fields
	  //
	  voxel[i][j][k].v1 = (2.0*(v3-v11+v18)+(v1-v12)*zz_by_2)
	    /(4.0+zz)+(2.0*(v2+v9+v13*yx)
		      -(v1+v12)*yx_by_2)/(yx+4.0+g0);

	  voxel[i][j][k].v2 = (2.0*(v6-v10-v17)+(v2-v9)*zy_by_2)
	    /(4.0+zy)+(2.0*(v1+v12+v13*yx)
		      -(v2+v9)*yx_by_2)/(4.0+yx+g0);

	  voxel[i][j][k].v3 = (2.0*(v1-v12-v18)+(v3-v11)*zz_by_2)
	    /(4.0+zz)+(2.0*(v4+v8+v14*yy)
		      -(v3+v11)*yy_by_2)/(4.0+yy+g0);

	  voxel[i][j][k].v4 = (2.0*(v5-v7+v16)+(v4-v8)*zx_by_2)
	    /(4.0+zx)+(2.0*(v3+v11+v14*yy)
		      -(v4+v8)*yy_by_2)/(4.0+yy+g0);

	  voxel[i][j][k].v5 = (2.0*(v4-v8-v16)+(v5-v7)*zx_by_2)
	    /(4.0+zx)+(2.0*(v6+v10+v15*yz)
		      -(v5+v7)*yz_by_2)/(4.0+yz+g0);

	  voxel[i][j][k].v6 = (2.0*(v2-v9+v17)+(v6-v10)*zy_by_2)
	    /(4.0+zy)+(2.0*(v5+v7+v15*yz)
		      -(v6+v10)*yz_by_2)/(4.0+yz+g0);

	  voxel[i][j][k].v7 = (2.0*(-v4+v8+v16)-(v5-v7)*zx_by_2)
	    /(4.0+zx)+(2.0*(v6+v10+v15*yz)
		      -(v5+v7)*yx_by_2)/(4.0+yz+g0);

	  voxel[i][j][k].v8 = (2.0*(-v5+v7-v16)-(v4-v8)*zx_by_2)
	    /(4.0+zx)+(2.0*(v3+v11+v14*yy)
		      -(v4+v8)*yy_by_2)/(4.0+yy+g0);

	  voxel[i][j][k].v9 = (2.0*(-v6+v10+v17)-(v2-v9)*zy_by_2)
	    /(4.0+zy)+(2.0*(v1+v12+v13*yx)
		      -(v2+v9)*yx_by_2)/(4.0+yx+g0);

	  voxel[i][j][k].v10 = (2.0*(-v2+v9-v17)-(v6-v10)*zy_by_2)
	    /(4.0+zy)+(2.0*(v5+v7+v15*yz)
		      -(v6+v10)*yz_by_2)/(4.0+yz+g0);

	  voxel[i][j][k].v11 = (2.0*(-v1+v12+v18)-(v3-v11)*zz_by_2)
	    /(4.0+zz)+(2.0*(v4+v8+v14*yy)
		      -(v3+v11)*yy_by_2)/(4.0+yy+g0);

	  voxel[i][j][k].v12 = (2.0*(-v3+v11-v18)-(v1-v12)*zz_by_2)
	    /(4.0+zz)+(2.0*(v2+v9+v13*yx)
		      -(v1+v12)*yx_by_2)/(4.0+yx+g0);

	  //
	  // Upgrade the voltages representing the materials and voxel size
	  //
	  voxel[i][j][k].v13 = (2.0*(v1+v2+v9+v12)+(yx-4.0)*v13)/(4.0+yx+g0);

	  voxel[i][j][k].v14 = (2.0*(v3+v4+v8+v11)+(yy-4.0)*v14)/(4.0+yy+g0);

	  voxel[i][j][k].v15 = (2.0*(v5+v6+v7+v10)+(yz-4.0)*v15)/(4.0+yz+g0);

	  voxel[i][j][k].v16 = -(2.0*zx*(v4-v5+v7-v8)+(4.0-zx)*v16)/(4.0+zx);

	  voxel[i][j][k].v17 = -(2.0*zy*(-v2+v6+v9-v10)+(4.0-zy)*v17)/(4.0+zy);

	  voxel[i][j][k].v18 = -(2.0*zz*(v1-v3+v11-v12)+(4.0-zz)*v18)/(4.0+zz);

	

	  }
	}
      }









    //
    // Apply reflections along the outer faces of each voxel that
    // would require.  In the simple waveguide problem here, only
    // reflections on the voxels along the outer faces of the
    // simulation volume will require reflections applied.  The
    // boundary conditions are specified here to create a PEC
    // waveguide with a short termination.  By changing the
    // orientation of the reflections here, the waveguide could easily
    // be `turned' in the cartesian mesh.
    //
    //


    for(i=0;i<nx;i++)
      {
      for(j=0;j<ny;j++)
	{
	for(k=0;k<nz;k++)
	  {
	  //
	  // PEC reflections at the most positive x face
	  //
	  if(i==nx-1)
	    {
	    voxel[i][j][k].v10 *= -1.0;
	    voxel[i][j][k].v11 *= -1.0;
 	    }
	  //
	  // PEC reflections at the most positive y face
	  //
	  if(j==ny-1)
	    {
	    voxel[i][j][k].v7 *= -1.0;
	    voxel[i][j][k].v12 *= -1.0;
	    }
	  //
	  // PEC reflections at the most positive z face
	  //
	  if(k==nz-1)
	    {
	    voxel[i][j][k].v8 *= -1.0;
	    voxel[i][j][k].v9 *= -1.0;
	    }
	  //
	  // Simple absorbing boundary at the most negative x face
	  //
	  if(i==0)
	    {
	    voxel[i][j][k].v3 *= 0.0;
	    voxel[i][j][k].v6 *= 0.0;
	    }
	  //
	  // PEC reflections at the most negative y face
	  //
	  if(j==0)
	    {
	    voxel[i][j][k].v1 *= -1.0;
	    voxel[i][j][k].v5 *= -1.0;
	    }
	  //
	  // PEC reflections at the most negative z face
	  //
	  if(k==0)
	    {
	    voxel[i][j][k].v2 *= -1.0;
	    voxel[i][j][k].v4 *= -1.0;
	    }



	  }
	}
      }





    //
    // Connect the voxels in the x direction.  Note that temporary
    // variables must be used since this is in essence a memory swap.
    // Without the tempoaries, the originals would be overwritten
    // before the swap was finished.
    //
    for(i=0;i<nx-1;i++)
      {
      for(j=0;j<ny;j++)
	{
	for(k=0;k<nz;k++)
	  {
	  //
	  v3 = voxel[i+1][j][k].v3;
	  voxel[i+1][j][k].v3 = voxel[i][j][k].v11;
	  voxel[i][j][k].v11 = v3;
	  //
	  v6 = voxel[i+1][j][k].v6;
	  voxel[i+1][j][k].v6 = voxel[i][j][k].v10;
	  voxel[i][j][k].v10 = v6;
	  }
	}
      }

    //
    // Connect the voxels in the y direction
    //
    for(j=0;j<ny-1;j++)
      {
      for(i=0;i<nx;i++)
	{
	for(k=0;k<nz;k++)
	  {
	  //
	  v5 = voxel[i][j+1][k].v5;
	  voxel[i][j+1][k].v5 = voxel[i][j][k].v7;
	  voxel[i][j][k].v7 = v5;
	  //
	  v1 = voxel[i][j+1][k].v1;
	  voxel[i][j+1][k].v1 = voxel[i][j][k].v12;
	  voxel[i][j][k].v12 = v1;
	  }
	}
      }

    //
    // Connect the voxels in the z direction
    //
    for(k=0;k<nz-1;k++)
      {
      for(i=0;i<nx;i++)
	{
	for(j=0;j<ny;j++)
	  {
 	  //
	  v4 = voxel[i][j][k+1].v4;
	  voxel[i][j][k+1].v4 = voxel[i][j][k].v8;
	  voxel[i][j][k].v8 = v4;
	  //
	  v2 = voxel[i][j][k+1].v2;
	  voxel[i][j][k+1].v2 = voxel[i][j][k].v9;
	  voxel[i][j][k].v9 = v2;
	  }
	}
      }




    fprintf(stdout, "\n");



    }







  //
  // Time in simulated seconds that the simulation has progressed.
  //
  currentSimulatedTime = dt*dx*(float)iteration;  
  //
  // Print to standard output the iteration number and current simulated time.
  //
  fprintf(stdout, "#%d %lgsec", iteration, currentSimulatedTime);


  //
  // Create the filename for this iteration, which includes the iteration number.
  //
  sprintf(filename, "c_%06d.bob", iteration);
  //
  // open a new data file for this iteration:
  //
  while ((openFilePointer = fopen(filename, "wb")) == NULL)
    {
    //
    // If unable to open the file, exit with a descriptive failure.
    //
    fprintf(stderr, "Difficulty opening c_%06d.bob", iteration);
    perror(" ");
    }
  //
  // Locate the min and max values in order to autoscale
  //
  min = FLT_MAX;
  max = -FLT_MAX;
  for(k=0;k<nz;k++)
    {
    for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
	{
	//
//  	ex=2.0*(voxel[i][j][k].v1 + voxel[i][j][k].v2 +
// 		voxel[i][j][k].v9 + voxel[i][j][k].v12 +
// 		yx*voxel[i][j][k].v13) / (u*(4.0+yx));
 	//
// 	ey=2.0*(voxel[i][j][k].v3 + voxel[i][j][k].v4 +
// 		voxel[i][j][k].v8 + voxel[i][j][k].v11 +
// 		yy*voxel[i][j][k].v14) / (v*(4.0+yy));
	//
	ez=2.0*(voxel[i][j][k].v5 + voxel[i][j][k].v6 +
		voxel[i][j][k].v7 + voxel[i][j][k].v10 +
		yz*voxel[i][j][k].v15) / (w*(4.0+yz));
	//
// 	hx=-2.0*(voxel[i][j][k].v4 - voxel[i][j][k].v5 +
// 		 voxel[i][j][k].v7 - voxel[i][j][k].v8 -
// 		 voxel[i][j][k].v16) / (Z_0*u*(4.0+zx));
 	//
// 	hy=-2.0*(-voxel[i][j][k].v2 + voxel[i][j][k].v6 +
// 		 voxel[i][j][k].v9 - voxel[i][j][k].v10 -
// 		 voxel[i][j][k].v17)/(Z_0*v*(4.0+zy));
 	//
// 	hz=-2.0*(-voxel[i][j][k].v3 + voxel[i][j][k].v1 +
// 		 voxel[i][j][k].v11 - voxel[i][j][k].v12 -
// 		 voxel[i][j][k].v18)/(Z_0*w*(4.0+zz));
	//
	//
	//
	if (ez < min)
	  {
	  min = ez;
	  }
	if (ez > max)
	  {
	  max = ez;
	  }
	}
      }
    }
  //
  // Set norm to be max or min, whichever is greater in magnitude.
  //
  norm = (fabs(max) > fabs(min)) ? fabs(max) : fabs(min);
  if (norm == 0.0)
    {
    //
    // If everything is zero, give norm a tiny value to avoid division by zero.
    //
    norm = DBL_EPSILON;
    }
  scalingValue = 127.0/norm;
  //
  // Write to standard output the minimum and maximum values from this iteration
  // and the minimum and maximum values that will be written to the bob file this iteration.
  //
  fprintf(stdout, "\t%lg(%d) < ez BoB < %lg(%d)",
	  min, (int)(128.0 + scalingValue*min),
	  max, (int)(128.0 + scalingValue*max));
  //
  // Scale each ez value in the mesh to the range of integers from zero through 254
  // and write them to the output file for this iteration.
  //
  for(k=0;k<nz;k++)
    {
    for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
	{
	//
// 	ex=2.0*(voxel[i][j][k].v1 + voxel[i][j][k].v2 +
// 		voxel[i][j][k].v9 + voxel[i][j][k].v12 +
// 		yx*voxel[i][j][k].v13) / (u*(4.0+yx));
 	//
// 	ey=2.0*(voxel[i][j][k].v3 + voxel[i][j][k].v4 +
// 		voxel[i][j][k].v8 + voxel[i][j][k].v11 +
// 		yy*voxel[i][j][k].v14) / (v*(4.0+yy));
	//
	ez=2.0*(voxel[i][j][k].v5 + voxel[i][j][k].v6 +
		voxel[i][j][k].v7 + voxel[i][j][k].v10 +
		yz*voxel[i][j][k].v15) / (w*(4.0+yz));
	//
// 	hx=-2.0*(voxel[i][j][k].v4 - voxel[i][j][k].v5 +
// 		 voxel[i][j][k].v7 - voxel[i][j][k].v8 -
// 		 voxel[i][j][k].v16) / (Z_0*u*(4.0+zx));
 	//
// 	hy=-2.0*(-voxel[i][j][k].v2 + voxel[i][j][k].v6 +
// 		 voxel[i][j][k].v9 - voxel[i][j][k].v10 -
// 		 voxel[i][j][k].v17)/(Z_0*v*(4.0+zy));
 	//
// 	hz=-2.0*(-voxel[i][j][k].v3 + voxel[i][j][k].v1 +
// 		 voxel[i][j][k].v11 - voxel[i][j][k].v12 -
// 		 voxel[i][j][k].v18)/(Z_0*w*(4.0+zz));
 	//
	// Put the value out to the file
	//
	putc((int)(128.0 + scalingValue*ez), openFilePointer);
	}
      }
    }
 
  //
  // Close the output file for this iteration.
  //
  fclose(openFilePointer);
  //
  // Write the dimensions and name of the output file for this
  // iteration to the viz control.
  //
  fprintf(vizFilePointer, "%dx%dx%d %s\n", nx+1, ny+1, nz, filename);
  //
  // Write identification of the corners of the mesh and the max and
  // min values for this iteration to the viz control file.
  //
  fprintf(vizFilePointer, "bbox: 0.0 0.0 0.0 %lg %lg %lg %lg %lg\n",
	  dx*(double)nx, dy*(double)ny, dz*(double)nz, min, max);
  //
  // Close the viz control file
  //
  fclose(vizFilePointer);



  fprintf(stdout, "\n");

  //
  // Output some of the simulation information to the user which is
  // particularly useful after the simulation is done.
  //
  fprintf(stdout, "\n"); 
  fprintf(stdout, "bob -cmap gbry.cmap -s %dx%dx%d *.bob\n", nx, ny, nz); 
  fprintf(stdout, "\n"); 
  fprintf(stdout, "viz ToyTLMc.viz\n"); 
  fprintf(stdout, "\n"); 
  fprintf(stdout, "Dynamically allocated %d bytes\n", allocatedBytes); 
  fprintf(stdout, "\n"); 
  fprintf(stdout, "Meshing parameters:\n"); 
  fprintf(stdout, "%dx%dx%d cells\n", nx, ny, nz); 
  fprintf(stdout, "dx=%lg, dy=%lg, dz=%lg meters\n", dx, dy, dz); 
  fprintf(stdout, "u=%lg, v=%lg, w=%lg\n", u, v, w); 
  fprintf(stdout, "%lg x %lg x %lg meter^3 simulation region\n",  
	  GUIDE_WIDTH, GUIDE_HEIGHT, LENGTH_IN_WAVELENGTHS*lambda0); 
  fprintf(stdout,"dt=%g \n",dt*dx);
  fprintf(stdout, "\n"); 

  // Return exit status to the environment.
  return(0);
}
