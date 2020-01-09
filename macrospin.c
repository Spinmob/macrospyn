#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


////////////////////////////////
// CONSTANTS
////////////////////////////////
#define PI 3.1415926535897
#define u0 1.25663706e-6

// File handle for log
FILE *log_file;

// Main data structure, allowing user to specify parameters and 
// external drives.
struct DATA;
typedef struct DATA
{

  // Integer for how much information to log.
  int log_level;

  //////////////////////////////
  // Settings and parameters
  //////////////////////////////

  double Byx1;           // in-plane hard axis demag field (T), element 1
  double Bzx1;           // out-of-plane hard axis demag field (T), element 1
  double Byx2;           // in-plane hard axis demag field (T), element 2
  double Bzx2;           // out-of-plane hard axis demag field (T), element 2

  double Bex1;           // Exchange-like field (T) experienced by element 1 (from element 2)
  double Bex2;           // Exchange-like field (T) experienced by element 2 (from element 1)

  double damping1;       // Gradient damping (unitless), element 1
  double damping2;       // Gradient damping (unitless), element 2

  double g1;             // magnitude of the gyromagnetic ratio (rad / s T) for element 1
  double g2;             // magnitude of the gyromagnetic ratio (rad / s T) for element 2

  //////////////////////////////////////////
  // USER SUPPLIED ARRAYS
  //////////////////////////////////////////
  
  // Array of fields (T) externally applied to by each element
  double *applied_Bx1, *applied_By1, *applied_Bz1;
  double *applied_Bx2, *applied_By2, *applied_Bz2; 

  // Array of additional torques (inverse seconds) applied to each element's magnetization unit vector
  double *applied_tx1, *applied_ty1, *applied_tz1;
  double *applied_tx2, *applied_ty2, *applied_tz2;

  //////////////////////////////////////////
  // SOLVER STUFF
  //////////////////////////////////////////

  int    steps; // total number of steps
  double dt;    // time step (s)

  // Solution arrays
  double *Bx1, *By1, *Bz1, *Bx2, *By2, *Bz2; // Total effective fields.
  double *tx1, *ty1, *tz1, *tx2, *ty2, *tz2; // Total additional torques.
  double *mx1, *my1, *mz1, *mx2, *my2, *mz2; // Magnetization unit vectors.

  // Intermediate values (differential steps returned by Dm)
  double dmx1, dmy1, dmz1;
  double dmx2, dmy2, dmz2;

} DATA;


///////////////////////////////////
// ACTUAL SOLVER CODE
///////////////////////////////////

// Definitions of time differentials for the free and fixed layers
double calculate_dms(DATA *d, int n)
{
  // Calculate the instantaneous total effective magnetic field, including:
  //          applied                 demag                         exchange
  d->Bx1[n] = d->applied_Bx1[n]   -   0                      +   d->Bex1 * d->mx2[n];
  d->By1[n] = d->applied_By1[n]   -   d->Byx1 * d->my1[n]    +   d->Bex1 * d->my2[n];  
  d->Bz1[n] = d->applied_Bz1[n]   -   d->Bzx1 * d->mz1[n]    +   d->Bex1 * d->mz2[n];  
  
  // Calculate the instantaneous total additional torques, including:
  //          applied                 spin transfer
  d->tx1[n] = d->applied_tx1[n];
  d->ty1[n] = d->applied_ty1[n];
  d->tz1[n] = d->applied_tz1[n];

  // Calculate the magnetization unit vector step, including:
  
  //          Precession around the total effective field                     Gradient damping: steps down the energy gradient from the total magnetic field.                                                                                       Other torques
  d->dmx1 = ( d->g1 * ( d->By1[n] * d->mz1[n] - d->Bz1[n] * d->my1[n] )   +   d->g1 * d->damping1 * ( ( d->my1[n] * d->my1[n] + d->mz1[n] * d->mz1[n] ) * d->Bx1[n] - d->mx1[n] * d->my1[n] * d->By1[n] - d->mx1[n] * d->mz1[n] * d->Bz1[n] )   +   d->tx1[n]     )*d->dt;
  d->dmy1 = ( d->g1 * ( d->Bz1[n] * d->mx1[n] - d->Bx1[n] * d->mz1[n] )   +   d->g1 * d->damping1 * ( ( d->mz1[n] * d->mz1[n] + d->mx1[n] * d->mx1[n] ) * d->By1[n] - d->my1[n] * d->mz1[n] * d->Bz1[n] - d->my1[n] * d->mx1[n] * d->Bx1[n] )   +   d->ty1[n]     )*d->dt;
  d->dmz1 = ( d->g1 * ( d->Bx1[n] * d->my1[n] - d->By1[n] * d->mx1[n] )   +   d->g1 * d->damping1 * ( ( d->mx1[n] * d->mx1[n] + d->my1[n] * d->my1[n] ) * d->Bz1[n] - d->mz1[n] * d->mx1[n] * d->Bx1[n] - d->mz1[n] * d->my1[n] * d->By1[n] )   +   d->tz1[n]     )*d->dt;
}


void Solve_Huen(DATA *d)
{
  // This just runs through the time-domain once, populating all the solution arrays.

          // Log output indented for easy viewing of "important" code.
          if(d->log_level) log_file = fopen("macrospin.log", "w");
          if(d->log_level) fprintf(log_file, "Solve_Huen... FIGHT!\n");

  // changes in the unit vectors per dt
  double dmx1,dmy1,dmz1,dmx2,dmy2,dmz2;

  // Factor by which to normalize the components at each step
  // This variable will be adjusted after each step.
  double norminator = 1.0/sqrt(d->mx1[0] * d->mx1[0] 
                             + d->my1[0] * d->my1[0] 
                             + d->mz1[0] * d->mz1[0]);
  d->mx1[0] = d->mx1[0]*norminator;
  d->my1[0] = d->my1[0]*norminator;
  d->mz1[0] = d->mz1[0]*norminator;

  // solve the time dependence
  int    n = 0;                      // step number
  double dmx1i, dmy1i, dmz1i;        // initial Euler step
  for(int n=0; n <= d->steps-2; n++) // This loop calculates the n+1 index, and doesn't need to do anything when n=steps-1
  {
    /* 
    
    Heun method: with our derivative step D(d,n), we calculate intermediate value
    
      yi[n+1] = y[n] + D(y,n)

    then get a better estimate
      
      y[n+1] = y[n] + 0.5*( D(y,n) + D(yi,n+1) )

    Importantly, Dm(d,n) involves the current magnetization, field, etc, 
    whereas Dm(di, n+1) involves the intnermediate magnetization, field, etc at the next step.

    */

    // Calculate the Euler-method step. This populates the step values d->dmx1, d->dmy1, ...
    calculate_dms(d, n);
    
    // Keep these steps for the Heun result below
    dmx1i = d->dmx1;
    dmy1i = d->dmy1;
    dmz1i = d->dmz1;

    // Step the magnetization unit vector to the intermediate value using Euler.
    // This will be overwritten (corrected) at the end of the Heun method.
    d->mx1[n+1] = d->mx1[n] + dmx1i; 
    d->my1[n+1] = d->my1[n] + dmy1i;
    d->mz1[n+1] = d->mz1[n] + dmz1i;

    // Normalize the new magnetization using the Taylor expansion of sqrt() near 1 to speed up the calculation.
    norminator = 1.0/(1.0 + 0.5 * (d->mx1[n+1]*d->mx1[n+1] + d->my1[n+1]*d->my1[n+1] + d->mz1[n+1]*d->mz1[n+1] - 1) );
    d->mx1[n+1] = d->mx1[n+1]*norminator;
    d->my1[n+1] = d->my1[n+1]*norminator;
    d->mz1[n+1] = d->mz1[n+1]*norminator;

    // Calculate the magnetization change for the *next* step using the intermediate Euler result.
    calculate_dms(d, n+1);

    // Overwrite the Euler result with the Heun estimate
    d->mx1[n+1] = d->mx1[n] + 0.5*(dmx1i + d->dmx1);
    d->my1[n+1] = d->my1[n] + 0.5*(dmx1i + d->dmy1);
    d->mz1[n+1] = d->mz1[n] + 0.5*(dmx1i + d->dmz1);

    // Normalize the new magnetization using the Taylor expansion of sqrt() near 1 to speed up the calculation.
    norminator = 1.0/(1.0 + 0.5 * (d->mx1[n+1]*d->mx1[n+1] + d->my1[n+1]*d->my1[n+1] + d->mz1[n+1]*d->mz1[n+1] - 1) );
    d->mx1[n+1] = d->mx1[n+1]*norminator;
    d->my1[n+1] = d->my1[n+1]*norminator;
    d->mz1[n+1] = d->mz1[n+1]*norminator;

  } // end of for loop

  // Update the total field & torque values by calculating the next step.
  calculate_dms(d,n);

  if(d->log_level) fprintf(log_file, "Ada!");
  if(d->log_level) fclose(log_file);
}
