/**
 * This file is part of the Macrospyn distribution 
 * (https://github.com/Spinmob/macrospyn).
 * Copyright (c) 2002-2020 Jack Childress (Sankey).
 * 
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.

 * Building on Linux requires gcc.
 * Building on Windows is easiest with the dev-C++ IDE.
**/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int log_level=0;
FILE *log_file;

///////////////////////////////
// DOMAIN class
///////////////////////////////

// Data structure for each domain. 
// We use a struct because there is no cross-platform-compiler means of
// having Python interact with a C++ class. Blame the C++ compilers.
typedef struct domain {

    //////////////////////////////
    // Model Inputs 
    //////////////////////////////

    // Magnitude of the gyromagnetic ratio [radians / (sec T)]
    double gamma, *gammas;

    // Saturation magnetization u0*Ms [T]
    double M, *Ms;

    // Gilbert damping parameter [unitless]
    double alpha, *alphas;         

    // Exchange-like field strength [T], applied in the direction of the other domain's unit vector
    double X, *Xs; 

    // Spin transfer torque (rate) parallel to other domain [rad / s]
    double s, *ss;

    // Other torque (rate) unrelated to either domain [rad / s]
    double tx, *txs;
    double ty, *tys; 
    double tz, *tzs; 

    // Applied field components [T]
    double Bx, *Bxs;
    double By, *Bys; 
    double Bz, *Bzs; 

    // Anisotropy tensor elements [unitless], defined such that Nxx+Nyy+Nzz=1 for an aligned ellipsoid
    double Nxx, *Nxxs, Nxy, *Nxys, Nxz, *Nxzs; 
    double Nyx, *Nyxs, Nyy, *Nyys, Nyz, *Nyzs;
    double Nzx, *Nzxs, Nzy, *Nzys, Nzz, *Nzzs;

    // Dipole tensor [unitless], representing the fraction of the other layer's saturation magnetization
    double Dxx, *Dxxs, Dxy, *Dxys, Dxz, *Dxzs;       
    double Dyx, *Dyxs, Dyy, *Dyys, Dyz, *Dyzs; 
    double Dzx, *Dzxs, Dzy, *Dzys, Dzz, *Dzzs;

    //////////////////////////////
    // Settings
    //////////////////////////////

    // 0 = disabled, 1 = LLG
    int mode;    

    // Initial conditions
    double x0, y0, z0;

    //////////////////////////////////////////
    // SOLVER STUFF
    //////////////////////////////////////////

    // Solution arrays
    double *x, *y, *z;  // Magnetization unit vector

} domain; 

// Sets all the instantaneous model inputs for step n.
// We do this to ease the user's ability to input parameters
// vs arrays.
void get_input_parameters(domain *a, int n) {

    // Always check that the array exists first, then assume it's
    // of sufficient length.
    if(a->gammas != NULL) a->gamma = a->gammas[n]; // gyro
    if(a->Ms     != NULL) a->M     = a->Ms[n];     // magnetization
    if(a->alphas != NULL) a->alpha = a->alphas[n]; // damping
    if(a->Xs     != NULL) a->X     = a->Xs[n];     // exchange
    if(a->ss     != NULL) a->s     = a->ss[n];     // spin torque
    
    if(a->Bxs != NULL) a->Bx = a->Bxs[n]; // applied field
    if(a->Bys != NULL) a->By = a->Bys[n];
    if(a->Bzs != NULL) a->Bz = a->Bzs[n];

    if(a->txs != NULL) a->tx = a->txs[n]; // applied torque
    if(a->tys != NULL) a->ty = a->tys[n];
    if(a->tzs != NULL) a->tz = a->tzs[n];

    if(a->Nxxs != NULL) a->Nxx = a->Nxxs[n]; // anisotropy
    if(a->Nxys != NULL) a->Nxy = a->Nxys[n];
    if(a->Nxzs != NULL) a->Nxz = a->Nxzs[n];
    if(a->Nyxs != NULL) a->Nyx = a->Nyxs[n];
    if(a->Nyys != NULL) a->Nyy = a->Nyys[n];
    if(a->Nyzs != NULL) a->Nyz = a->Nyzs[n];
    if(a->Nzxs != NULL) a->Nzx = a->Nzxs[n];
    if(a->Nzys != NULL) a->Nzy = a->Nzys[n];
    if(a->Nzzs != NULL) a->Nzz = a->Nzzs[n];
    
    if(a->Dxxs != NULL) a->Dxx = a->Dxxs[n]; // dipole
    if(a->Dxys != NULL) a->Dxy = a->Dxys[n];
    if(a->Dxzs != NULL) a->Dxz = a->Dxzs[n];
    if(a->Dyxs != NULL) a->Dyx = a->Dyxs[n];
    if(a->Dyys != NULL) a->Dyy = a->Dyys[n];
    if(a->Dyzs != NULL) a->Dyz = a->Dyzs[n];
    if(a->Dzxs != NULL) a->Dzx = a->Dzxs[n];
    if(a->Dzys != NULL) a->Dzy = a->Dzys[n];
    if(a->Dzzs != NULL) a->Dzz = a->Dzzs[n];
};

// Calculate a single step for this domain, if enabled.
// Parameters
//   domain *a    The domain whose step we wish to calculate.
//   domain *b    The "other" domain that exerts exchange fields, dipolar fields, and spin transfer.
//   int n        The step at which to calculate.
void D(domain *a, domain *b, int n, double dt, double *dx, double *dy, double *dz) {

    // At each step (including intermediate steps), make sure to get 
    // the most current model input values.
    get_input_parameters(a, n);
    get_input_parameters(b, n);
    
    // If our domain's dynamics are not enabled, no step
    if(a->mode == 0) {
        *dx = *dy = *dz = 0;
        return;
    }

    // Intermediate values
    double Nx, Ny, Nz; // aNisotropy field [T] from this domain
    double Dx, Dy, Dz; // Dipolar field [T] from "other" domain
    double Xx, Xy, Xz; // Exchange field [T] from "other" domain
    double Tx, Ty, Tz; // Non-damping effective field [T]
    double vx, vy, vz; // Total non-damping torque [rad/sec]

    // First let's calculate the aNisotropy field
    Nx = a->M*(a->Nxx*a->x[n] + a->Nxy*a->y[n] + a->Nxz*a->z[n]);
    Ny = a->M*(a->Nyy*a->y[n] + a->Nyz*a->z[n] + a->Nyx*a->x[n]);
    Nz = a->M*(a->Nzz*a->z[n] + a->Nzx*a->x[n] + a->Nzy*a->y[n]);
    
    // Now the Dipolar field from b
    Dx = b->M*(a->Dxx*b->x[n] + a->Dxy*b->y[n] + a->Dxz*b->z[n]);
    Dy = b->M*(a->Dyy*b->y[n] + a->Dyz*b->z[n] + a->Dyx*b->x[n]);
    Dz = b->M*(a->Dzz*b->z[n] + a->Dzx*b->x[n] + a->Dzy*b->y[n]);
    
    // Now the eXchange field from b
    Xx = a->X*b->x[n];
    Xy = a->X*b->y[n];
    Xz = a->X*b->z[n];

    // Now we can get the components of T
    Tx = -a->gamma*(a->Bx+Nx+Dx+Xx) + a->s*(a->z[n]*b->y[n] - a->y[n]*b->z[n]) + a->z[n]*a->ty - a->y[n]*a->tz;
    Ty = -a->gamma*(a->By+Ny+Dy+Xy) + a->s*(a->x[n]*b->z[n] - a->z[n]*b->x[n]) + a->x[n]*a->tz - a->z[n]*a->tx;
    Tz = -a->gamma*(a->Bz+Nz+Dz+Xz) + a->s*(a->y[n]*b->x[n] - a->x[n]*b->y[n]) + a->y[n]*a->tx - a->x[n]*a->ty;
    
    // Now we can compute the total non-damping torque for this step.
    vx = a->y[n]*Tz - a->z[n]*Ty; 
    vy = a->z[n]*Tx - a->x[n]*Tz; 
    vz = a->x[n]*Ty - a->y[n]*Tx; 

    // We store the step magnitude to help with Heun method.
    double scale = dt/(1.0+a->alpha*a->alpha);
    *dx = ( vx + a->alpha*(a->y[n]*vz-a->z[n]*vy) + a->alpha*a->alpha*a->x[n]*(a->x[n]*vx+a->y[n]*vy+a->z[n]*vz) ) * scale;
    *dy = ( vy + a->alpha*(a->z[n]*vx-a->x[n]*vz) + a->alpha*a->alpha*a->y[n]*(a->y[n]*vy+a->z[n]*vz+a->x[n]*vx) ) * scale;
    *dz = ( vz + a->alpha*(a->x[n]*vy-a->y[n]*vx) + a->alpha*a->alpha*a->z[n]*(a->z[n]*vz+a->x[n]*vx+a->y[n]*vy) ) * scale;
}


///////////////////////////////////
// LOG STUFF
///////////////////////////////////
void log_step(domain *a, domain *b, int n) {
    fprintf(log_file, "n=%i --------------------------%p\n", n, a->x);
    fprintf(log_file, "  a->gamma=%f, pointer=%p\n",     a->gamma, a->gammas);
    fprintf(log_file, "  a->M=%f,     pointer=%p\n",     a->M,     a->Ms);
    fprintf(log_file, "  a->alpha=%f, pointer=%p\n",     a->alpha, a->alphas);
    fprintf(log_file, "  a->X=%f,     pointer=%p\n",     a->X,     a->Xs);
    fprintf(log_file, "  a->s=%f,     pointer=%p\n",     a->s,     a->ss);

    fprintf(log_file, "\n");
}


///////////////////////////////////
// SOLVER
///////////////////////////////////

int solve_heun(domain *a, domain *b, double dt, int N) {
 
  long t0 = time(0);

  // Log file
  if(log_level > 0) {
    log_file = fopen("engine.log", "w");
    fprintf(log_file, "solve_heun() beings %li\n------------------------------------------------\n\n", t0);
    log_step(a,b,0);
  }

  // The initial condition of the magnetization is assumed to be the 
  // first element of the array, but we should make sure it's length is 1!
  double scale;
  
  // Normalize a
  scale = 1.0/sqrt(a->x[0]*a->x[0] + a->y[0]*a->y[0] + a->z[0]*a->z[0]);
  a->x[0] *= scale;
  a->y[0] *= scale;
  a->z[0] *= scale;

  // Normalize b
  scale = 1.0/sqrt(b->x[0]*b->x[0] + b->y[0]*b->y[0] + b->z[0]*b->z[0]);
  b->x[0] *= scale;
  b->y[0] *= scale;
  b->z[0] *= scale;
  
  // These will hold the step values calculated by D();
  double adx1, ady1, adz1, bdx1, bdy1, bdz1;
  double adx2, ady2, adz2, bdx2, bdy2, bdz2;
  
  if(log_level >=1) fprintf(log_file, "STARTING LOOP: N=%i steps\n", N);

  // Now do the Heun loop
  // We don't go to the end because we don't want to overwrite the first step.
  for(int n=0; n<=N-2; n++) {
   
    //  Heun method: with our derivative step dy(y,n), we calculate intermediate value
    //
    //    yi[n+1] = y[n] + dy(y,n)
    //
    //  then get a better estimate
    //    
    //    y[n+1] = y[n] + 0.5*( dy(y,n) + dy(yi,n+1) )
    //
    //  Importantly, dy(y,n) involves the current magnetization, field, etc, 
    //  whereas dy(yi, n+1) involves the intnermediate magnetization, field, etc at the next step.
    
    // Calculate dy(y,n)
    D(a, b, n, dt, &adx1, &ady1, &adz1);
    D(b, a, n, dt, &bdx1, &bdy1, &bdz1);

    if(log_level >= 3 && (n % (int)(N/5) == 0 || n<5 || N-n<7)) log_step(a, b, n);
    
    // Store the intermediate value yi at n+1
    a->x[n+1] = a->x[n] + adx1; 
    a->y[n+1] = a->y[n] + ady1; 
    a->z[n+1] = a->z[n] + adz1;
    b->x[n+1] = b->x[n] + bdx1; 
    b->y[n+1] = b->y[n] + bdy1; 
    b->z[n+1] = b->z[n] + bdz1;
    
    // Calculate dy(yi,n+1)
    D(a, b, n+1, dt, &adx2, &ady2, &adz2);
    D(b, a, n+1, dt, &bdx2, &bdy2, &bdz2);

    // Get the Heun step
    a->x[n+1] = a->x[n] + 0.5*(adx1 + adx2);
    a->y[n+1] = a->y[n] + 0.5*(ady1 + ady2);
    a->z[n+1] = a->z[n] + 0.5*(adz1 + adz2);
    b->x[n+1] = b->x[n] + 0.5*(bdx1 + bdx2);
    b->y[n+1] = b->y[n] + 0.5*(bdy1 + bdy2);
    b->z[n+1] = b->z[n] + 0.5*(bdz1 + bdz2);

    // Normalize the new magnetization using the Taylor expansion of sqrt() near 1 to speed up the calculation.
    double norminator;
    
    // domain a
    norminator = 1.0/(1.0 + 0.5 * (a->x[n+1]*a->x[n+1] + a->y[n+1]*a->y[n+1] + a->z[n+1]*a->z[n+1] - 1.0) );
    a->x[n+1] *= norminator;
    a->y[n+1] *= norminator;
    a->z[n+1] *= norminator;

    // domain b
    norminator = 1.0/(1.0 + 0.5 * (b->x[n+1]*b->x[n+1] + b->y[n+1]*b->y[n+1] + b->z[n+1]*b->z[n+1] - 1.0) );
    b->x[n+1] *= norminator;
    b->y[n+1] *= norminator;
    b->z[n+1] *= norminator;

  } // End of for loop.

  // At this point, the whole solution arrays should be populated.
  if(log_level>0) {
    fprintf(log_file, "\n\n------------------------------------------------\nDone after %li", time(0)-t0);
    fclose(log_file);
  }
}

