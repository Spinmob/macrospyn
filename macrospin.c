#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


////////////////////////////////
// CONSTANTS
////////////////////////////////
#define PI 3.1415926535897
#define ME 9.10938356e-31
#define KB 1.38064852e-23
#define EC 1.60217662e-19
#define u0 1.25663706e-6
#define uB 9.274009994e-24
#define HBAR 1.0545718e-34
#define TO_RADIANS 0.0174532925                     // for speed.

FILE *log_file;

struct DATA;
typedef struct DATA
{
  double time_scale;                        // converts to nanoseconds.

  //////////////////////////////
  // Settings and parameters
  //////////////////////////////

  // Physical system
  double T;              // Temperature in K
  double thickness_fm;   // thickness of the ferromagnetic layer (nm)
  double thickness_nm;   // thickness of the normal metal layer (nm)
  double resistivity_fm; // uOhm-cm (only the ratio matters, though)
  double resistivity_nm; // uOhm-cm
  double length;         // in-plane length of structure (nm)
  double width;          // in-plane width of structure (nm)
  double ms;             // saturation magnetization of free layer (T) used only for spin transfer efficiency and Langevin force
  double Byx;            // in-plane hard axis demag field (T)
  double Bzx;            // out-of-plane demag field (T)
  double damping;        // Gilbert damping (unitless)

  // We're doing a hacky modification of a spin-valve code to model our nanowires.
  // We run our current through a Py/Pt nanowire. To estimate the vertical spin current
  // density, we multiply the charge current density in the Pt layer by the spin Hall angle,
  // then use the "Sinusoidal" spin torque. The fixed layer in this case simply defines
  // the polarization (though it doesn't exist in the real system).
  //
  // This uses the sinusoidal torque and efficiency after calculating the current density.

  int    torque_type; // 1 for spin Hall, 2 for sinusoid
  double hall_angle;  // Conversion from charge to spin current density (also having units of "Amps/m^2"
  double g_factor;    // Electron g-factor (2, usually)
  double efficiency;  // number of hbar/2's each electron imparts (at most)
  double gyromagnetic_magnitude; // 8.794033700300592e10*d->g_factor gyromagnetic ratio magnitude (rad / s T)  
  
  // Current drive I(t) = I0 + I1*cos(2*PI*f1*t)
  double I0;        // DC current (mA)
  double I1;        // RF amplitude (mA)
  double f1;        // RF frequency (GHz)
  double Bx_per_mA; // B-field along x per mA of current
  double By_per_mA; // B-field along y per mA of current
  double Bz_per_mA; // B-field along z per mA of current

  // Applied field
  double B0;          // Applied field (T)
  double Bx_hat;      // Components of B_hat unit vector
  double By_hat;
  double Bz_hat;

  //////////////////////////////////////////
  // SOLVER STUFF
  //////////////////////////////////////////

  int    steps; // total number of steps
  double t0;    // initial time (used only for RF current calculation)
  double dt;    // time step (ns). 1 period at 50 GHz = 0.02

  // Instantaneous values
  int    n;     // current step
  double mx;
  double my;
  double mz;
  double Mx;
  double My;
  double Mz;
  double Bx;
  double By;
  double Bz;
  double I;
  double Js;

  // Solution arrays
  double *solution_mx;
  double *solution_my;
  double *solution_mz;

  // Bureaucracy
  int log_level;

} DATA;




////////////////////////////////////////
// Standalone FUNCTIONS
////////////////////////////////////////
double Random_Gaussian()
{
  // this is from numerical recipes (sort of) I got it from some website, and it $
  // book's algorythm for finding a random number with a gaussian distribution
  // of width 1

  static int noExtra = 1;  // the algorythm is a cute trick
                           // that produces 2 numbers, if noExtra, make two new ones
  static double gset;
  double fac, r, v1, v2;

  if (noExtra)
  {
    do // get two random numbers between -1 and 1, and make sure they are within the unit circle
    {
      v1 = 2.0*rand()/RAND_MAX - 1;
      v2 = 2.0*rand()/RAND_MAX - 1;
      r = v1*v1 + v2*v2;           // need two random numbers in the unit circle
    } while (r >= 1.0 || r == 0);

    // now do the transformation
    fac = sqrt(-2.0*log(r)/r);
    gset = v1*fac;
    noExtra = 0; // we have an extra now
    return v2*fac;
  }
  else
  {
    noExtra = 1;
    return gset;
  }
}



//////////////////////////////
// Actuall solver functions
//////////////////////////////

double Beta(DATA *d) // spin transfer coefficient (depends on m's)
{
  // Prefactor on the spin transfer torque.
  if(d->torque_type==1 || d->torque_type==2)
    return(
           // sin torque (also for Hall); the spin current is calculated differently for Hall.
           // Note ms is in Tesla, so there is that extra u0
	   d->efficiency*d->g_factor*uB*u0/(2*EC*d->thickness_fm*1e-9*d->ms)
	   );
  // For other torque types (if I re-add them), this can depend on the magnetizations
}



// Definitions of time differentials for the free and fixed layers
double D_mx_Dt(DATA *d)
{
  return(
	  (
       // Precession around field
  	   d->gyromagnetic_magnitude*( d->By*d->mz-d->Bz*d->my )

       // Gradient damping
       + d->gyromagnetic_magnitude*d->damping*( (d->my*d->my+d->mz*d->mz)*d->Bx
                                - d->mx*d->my             *d->By
                                - d->mx*d->mz             *d->Bz )
       // Spin torque
       + Beta(d) * d->Js * ( d->mx*d->my*d->My + d->mx*d->mz*d->Mz - (d->my*d->my+d->mz*d->mz)*d->Mx )
	  ) * d->time_scale
	);

}

double D_my_Dt(DATA *d)	// returns the current value of dmydt at t, mx and my
{
  return(
	 (
	    d->gyromagnetic_magnitude*( d->Bz*d->mx-d->Bx*d->mz )
	  + d->gyromagnetic_magnitude*d->damping*( (d->mz*d->mz+d->mx*d->mx)*d->By - d->my*d->mz*d->Bz - d->my*d->mx*d->Bx )
	  + Beta(d)*d->Js *        ( d->my*d->mz*d->Mz + d->my*d->mx*d->Mx - (d->mz*d->mz+d->mx*d->mx)*d->My )
   ) * d->time_scale

	);
}
double D_mz_Dt(DATA *d)	// returns the current value of dmydt at t, mx and my
{
  // if(d->log_level >= 2)
  //   fprintf(log_file, " D_mz_Dt: Bx=%.1G By=%.1G mx=%.1G my=%.1G p=%.1G",
  //   d->Bx, d->By, d->mx, d->my, d->gyromagnetic_magnitude*( d->Bx * d->my - d->By * d->mx ));

  return(
	 (
	    d->gyromagnetic_magnitude*( d->Bx * d->my - d->By * d->mx )
	  + d->gyromagnetic_magnitude*d->damping*( (d->mx*d->mx+d->my*d->my)*d->Bz - d->mz*d->mx*d->Bx - d->mz*d->my*d->By )
	  + Beta(d)*d->Js*        ( d->mz*d->mx*d->Mx + d->mz*d->my*d->My - (d->mx*d->mx+d->my*d->my)*d->Mz )
   ) * d->time_scale
	);
}


void Solve_Huen(DATA *d)
{
  // This just runs through the time-domain once, populating all the solution_* variables.

          // Log output indented for easy viewing of "important" code.
          if(d->log_level) log_file = fopen("macrospin.log", "w");
          if(d->log_level) fprintf(log_file, "Solve_Huen... FIGHT!\n");

  // Seed the randomnumber generator
  srand(time(NULL));




  // Langevin field
  double Bx_langevin;
  double By_langevin;
  double Bz_langevin;

  // changes in the unit vectors per dt
  double dmx1,dmy1,dmz1;
  double dmx2,dmy2,dmz2;

  // factor by which to normalize the components at each step
  double norminator  = 1;




  // single-calculation of langevin pre-factor (speed boost)
  double Bl = sqrt(2.0*d->T*d->damping*KB/(d->gyromagnetic_magnitude*d->ms/u0*d->length*d->width*d->thickness_fm*1e-27)/(d->dt*d->time_scale));

          if(d->log_level) fprintf(log_file, "%.2f %.2f %.2f %.2f %.2f %.2f %.2G %.2G\n",
                        d->T, d->damping, d->ms, d->length, d->width, d->thickness_fm, d->dt, d->time_scale);

  // Conversion from mA of current to A/m^2 charge current
  double I_to_Js;
  if(d->torque_type==1)
    //   hall angle * current density through spin hall layer
    I_to_Js = d->hall_angle * 0.001 / (d->width*d->thickness_nm*1e-18)
              / (1.0 + d->resistivity_nm*d->thickness_fm/(d->resistivity_fm*d->thickness_nm));
  else
    I_to_Js = 0.001/(d->width*d->length*1e-18);



  
  // Initial normalization of all unit vectors
  norminator = 1.0/sqrt(d->mx*d->mx 
                      + d->my*d->my 
                      + d->mz*d->mz);
  d->mx = d->mx*norminator;
  d->my = d->my*norminator;
  d->mz = d->mz*norminator;
  
  norminator = 1.0/sqrt(d->Mx*d->Mx 
                      + d->My*d->My 
                      + d->Mz*d->Mz);
  d->Mx = d->Mx*norminator;
  d->My = d->My*norminator;
  d->Mz = d->Mz*norminator;  
    
  norminator = 1.0/sqrt(d->Bx_hat*d->Bx_hat 
                      + d->By_hat*d->By_hat 
                      + d->Bz_hat*d->Bz_hat);
  d->Bx_hat = d->Bx_hat*norminator;
  d->By_hat = d->By_hat*norminator;
  d->Bz_hat = d->Bz_hat*norminator;




  // Calculate the n=0 values to start the loop
  d->I  = d->I0 + d->I1*cos(2*PI*d->f1*d->t0); // Current through wire at t=t0 (mA)
  d->Js = I_to_Js*d->I;                        // Spin current density at t=t0 (A/m^2)
  Bx_langevin = Bl*Random_Gaussian();
  By_langevin = Bl*Random_Gaussian();
  Bz_langevin = Bl*Random_Gaussian();

          if(d->log_level) fprintf(log_file, "Bl=(%.1G,%.1G,%.1G) %.2f\n", Bx_langevin, By_langevin, Bz_langevin, Bl);



  // solve the time dependence
  d->n = 0;
  while(d->n < d->steps)
  {
    // Store the current step in the solution arrays
    d->solution_mx[d->n] = d->mx;
    d->solution_my[d->n] = d->my;
    d->solution_mz[d->n] = d->mz;

            if(d->log_level >= 2)
              fprintf(log_file, " %i m=(%.1G,%0.1G,%0.1G)",
              d->n, d->solution_mx[d->n], d->solution_my[d->n], d->solution_mz[d->n]);

    // Re-calculate the instantaneous field using the current values of
    // mx, my, mz, and this step's Langevin field.
    d->Bx = d->B0*d->Bx_hat                + Bx_langevin + d->Bx_per_mA*d->I;
    d->By = d->B0*d->By_hat - d->Byx*d->my + By_langevin + d->By_per_mA*d->I;
    d->Bz = d->B0*d->Bz_hat - d->Bzx*d->mz + Bz_langevin + d->Bz_per_mA*d->I;

            if(d->log_level >= 2)
              fprintf(log_file, " B=(%.1G,%.1G,%.1G) I=%.1G",
              d->Bx, d->By, d->Bz, d->I);

    // Euler method aproximate step, using this iteration's values
    dmx1 = d->dt * D_mx_Dt(d);
    dmy1 = d->dt * D_my_Dt(d);
    dmz1 = d->dt * D_mz_Dt(d);
    d->mx = d->mx + dmx1;
    d->my = d->my + dmy1;
    d->mz = d->mz + dmz1;

    // Normalize
    // Use a Taylor expansion of sqrt() near 1 to speed up the calculation.
    norminator = 1.0/(1.0 + 0.5 * (d->mx*d->mx + d->my*d->my + d->mz*d->mz - 1) );
    d->mx = d->mx*norminator;
    d->my = d->my*norminator;
    d->mz = d->mz*norminator;

    // Here we have set the actual values of mx, my, and mz to this intermediate
    // approximation, so that we can make the second step for the Huen approximation
    // For this we also want to increment the time and all other values
    d->n++;

    // Update the next step's charge current (mA) and spin current density (A/m^2)
    // SPEEDUP: Could use an approximate form here
    //          Could evaluate cos once, and have another solution array
    //          calculating the overlap with this function.
    if(d->I1 != 0) d->I  = d->I0 + d->I1*cos(2*PI*d->f1*(d->n*d->dt + d->t0));
    else           d->I  = d->I0;  // saves a cos() call for I1=0 simulations
    d->Js = d->I*I_to_Js;

    // Next step's Langevin field (T).
    // SPEEDUP: Could use rand() straight up (approximating Gaussian as flat)
    //          It will converge to the same result from the central limit theorem
    if(d->T>0)
    {
      Bx_langevin = Bl*Random_Gaussian();
      By_langevin = Bl*Random_Gaussian();
      Bz_langevin = Bl*Random_Gaussian();
    }
    else // Saves some processor time.
    {
      Bx_langevin = 0;
      By_langevin = 0;
      Bz_langevin = 0;
    }

    // Next step's *approximate* field based on Euler mx, my, and mz
    // (do NOT use this for the next iteration!)
    d->Bx = d->B0*d->Bx_hat                + Bx_langevin + d->Bx_per_mA*d->I;
    d->By = d->B0*d->By_hat - d->Byx*d->my + By_langevin + d->By_per_mA*d->I;
    d->Bz = d->B0*d->Bz_hat - d->Bzx*d->mz + Bz_langevin + d->Bz_per_mA*d->I;

    // Second value of dm/dt
    dmx2 = d->dt * D_mx_Dt(d);
    dmy2 = d->dt * D_my_Dt(d);
    dmz2 = d->dt * D_mz_Dt(d);

    // Now add the averaged increment to the original magnetization.
    // This is the Heun estimate for the next step.
    d->mx = d->solution_mx[d->n-1] + 0.5*(dmx1+dmx2);
    d->my = d->solution_my[d->n-1] + 0.5*(dmy1+dmy2);
    d->mz = d->solution_mz[d->n-1] + 0.5*(dmz1+dmz2);

    // normalize the final vector
    // Use a Taylor expansion of sqrt() near 1 to speed up the calculation.
    norminator = 1.0/(1.0 + 0.5 * (d->mx*d->mx + d->my*d->my + d->mz*d->mz - 1) );
    d->mx = d->mx*norminator;
    d->my = d->my*norminator;
    d->mz = d->mz*norminator;

            if(d->log_level>=2) fprintf(log_file, " dN=%.1G", norminator-1);

    // Note the field should now be re-calculated at the beginning of the loop using
    // these more accurate values of mx, my, and mz.

            if(d->log_level>=2) fprintf(log_file, "\n");

  } // end of for loop

  // Update the field and time t0. mx, my, and mz should already be updated.
  d->Bx = d->B0*d->Bx_hat                + Bx_langevin + d->Bx_per_mA*d->I;
  d->By = d->B0*d->By_hat - d->Byx*d->my + By_langevin + d->By_per_mA*d->I;
  d->Bz = d->B0*d->Bz_hat - d->Bzx*d->mz + Bz_langevin + d->Bz_per_mA*d->I;
  d->t0 += d->dt * d->steps;

  if(d->log_level) fprintf(log_file, "Ada!");
  if(d->log_level) fclose(log_file);
}
