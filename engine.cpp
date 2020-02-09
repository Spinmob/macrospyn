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

Building on linux requires gcc and g++ installed.
*/


///////////////////////////////
// Domain class
///////////////////////////////

// Data structure for each domain. 
// We use a struct because there is no cross-platform-compiler means of
// having Python interact with a C++ class.
struct Domain {

    //////////////////////////////
    // Model Inputs 
    //////////////////////////////

    // Magnitude of the gyromagnetic ratio [radians / (sec T)]
    double g, *gs;

    // Saturation magnetization u0*Ms [T]
    double M, *Ms;

    // Gilbert damping parameter (alpha) [unitless]
    double d, *ds;         

    // Exchange-like field strength [T], applied in the direction of the other domain's unit vector
    double X, *Xs; 

    // Spin transfer torque (rate) parallel to other domain [rad / s]
    double s, *ss;

    // Applied field components [T]
    double Bx, *Bxs, By, *Bys, Bz, *Bzs; 

    // Other torque (rate) unrelated to either domain [rad / s]
    double tx, *txs, ty, *tys, tz, *tzs; 

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

    bool mode;    // 0 = disabled, 1 = Gilbert damping

    //////////////////////////////////////////
    // SOLVER STUFF
    //////////////////////////////////////////

    // Solution arrays
    double *x, *y, *z;  // Magnetization unit vector

    // Temporary step values for Heun method
    double dx, dy, dz;
}; 

// Sets all the instantaneous model inputs for step n.
// We do this to ease the user's ability to input parameters
// vs arrays.
void get_input_parameters(Domain *a, int n) {

    // Always check that the array exists first, then assume it's
    // of sufficient length.
    if(a->gs != nullptr) a->g = a->gs[n]; // gyro
    if(a->Ms != nullptr) a->M = a->Ms[n]; // magnetization
    if(a->ds != nullptr) a->d = a->ds[n]; // damping
    if(a->Xs != nullptr) a->X = a->Xs[n]; // exchange
    if(a->ss != nullptr) a->s = a->ss[n]; // spin torque
    
    if(a->Bxs != nullptr) a->Bx = a->Bxs[n]; // applied field
    if(a->Bys != nullptr) a->By = a->Bys[n];
    if(a->Bzs != nullptr) a->Bz = a->Bzs[n];

    if(a->txs != nullptr) a->tx = a->txs[n]; // applied torque
    if(a->tys != nullptr) a->ty = a->tys[n];
    if(a->tzs != nullptr) a->tz = a->tzs[n];

    if(a->Nxxs != nullptr) a->Nxx = a->Nxxs[n]; // anisotropy
    if(a->Nxys != nullptr) a->Nxy = a->Nxys[n];
    if(a->Nxzs != nullptr) a->Nxz = a->Nxzs[n];
    if(a->Nyxs != nullptr) a->Nyx = a->Nyxs[n];
    if(a->Nyys != nullptr) a->Nyy = a->Nyys[n];
    if(a->Nyzs != nullptr) a->Nyz = a->Nyzs[n];
    if(a->Nzxs != nullptr) a->Nzx = a->Nzxs[n];
    if(a->Nzys != nullptr) a->Nzy = a->Nzys[n];
    if(a->Nzzs != nullptr) a->Nzz = a->Nzzs[n];
    
    if(a->Dxxs != nullptr) a->Dxx = a->Dxxs[n]; // dipole
    if(a->Dxys != nullptr) a->Dxy = a->Dxys[n];
    if(a->Dxzs != nullptr) a->Dxz = a->Dxzs[n];
    if(a->Dyxs != nullptr) a->Dyx = a->Dyxs[n];
    if(a->Dyys != nullptr) a->Dyy = a->Dyys[n];
    if(a->Dyzs != nullptr) a->Dyz = a->Dyzs[n];
    if(a->Dzxs != nullptr) a->Dzx = a->Dzxs[n];
    if(a->Dzys != nullptr) a->Dzy = a->Dzys[n];
    if(a->Dzzs != nullptr) a->Dzz = a->Dzzs[n];
}

// Calculate a single step for this domain, if enabled.
// Parameters
//   Domain *a    The domain whose step we wish to calculate.
//   Domain *b    The "other" domain that exerts exchange fields, dipolar fields, and spin transfer.
//   int n        The step at which to calculate.
void D(Domain *a, Domain *b, int n, double dt) {

    // At each step (including intermediate steps), make sure to get 
    // the most current model input values.
    get_input_parameters(a, n);
    get_input_parameters(b, n);
    
    // If our domain's dynamics are not enabled, no step
    if(a->mode == 0) {
    a->dx = a->dy = a->dz = 0;
    return;
    }

    // Intermediate values
    double Nx, Ny, Nz; // aNisotropy field [T] from this domain
    double Dx, Dy, Dz; // Dipolar field [T] from "other" domain
    double Xx, Xy, Xz; // Exchange field [T] from "other" domain
    double Tx, Ty, Tz; // Non-damping effective field [T]
    double vx, vy, vz; // Total non-damping torque [rad/sec]

    // We assume the instantaneous model input values have already
    // been calculated for both domains.

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
    Tx = -a->g*(a->Bx+Nx+Dx+Xx) + a->s*(a->z[n]*b->y[n] - a->y[n]*b->z[n]) + a->z[n]*a->ty - a->y[n]*a->tz;
    Ty = -a->g*(a->By+Ny+Dy+Xy) + a->s*(a->x[n]*b->z[n] - a->z[n]*b->x[n]) + a->x[n]*a->tz - a->z[n]*a->tx;
    Tz = -a->g*(a->Bz+Nz+Dz+Xz) + a->s*(a->y[n]*b->x[n] - a->x[n]*b->y[n]) + a->y[n]*a->tx - a->x[n]*a->ty;
    
    // Now we can compute the total non-damping torque for this step.
    vx = a->y[n]*Tz - a->z[n]*Ty; 
    vy = a->z[n]*Tx - a->x[n]*Tz; 
    vz = a->x[n]*Ty - a->y[n]*Tx; 

    // ... which can be used to get the next value for the magnetization
    
    // We store the d's to help with Heun method.
    double scale = dt/(1+a->d*a->d);
    a->dx = ( vx + a->d*(a->y[n]*vz-a->z[n]*vy) + a->d*a->d*a->x[n]*(a->x[n]*vx+a->y[n]*vy+a->z[n]*vz) ) * scale;
    a->dy = ( vy + a->d*(a->z[n]*vx-a->x[n]*vz) + a->d*a->d*a->y[n]*(a->y[n]*vy+a->z[n]*vz+a->x[n]*vx) ) * scale;
    a->dz = ( vz + a->d*(a->x[n]*vy-a->y[n]*vx) + a->d*a->d*a->z[n]*(a->z[n]*vz+a->x[n]*vx+a->y[n]*vy) ) * scale;
}



///////////////////////////////////
// SOLVER
///////////////////////////////////

int solve_heun(Domain *a, Domain *b, double dt) {

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
  
  // Get the number of steps from the supplied array length
  int N = sizeof(a->x)/sizeof(double);
  double adx, ady, adz, bdx, bdy, bdz;

  // Now do the Heun loop
  // We don't go to the end because we don't want to overwrite the first step.
  for(int n=0; n<N-2; n++) {
    
    /*
      Heun method: with our derivative step D(d,n), we calculate intermediate value
      
        yi[n+1] = y[n] + dy(y,n)

      then get a better estimate
        
        y[n+1] = y[n] + 0.5*( dy(y,n) + dy(yi,n+1) )

      Importantly, Dm(d,n) involves the current magnetization, field, etc, 
      whereas Dm(di, n+1) involves the intnermediate magnetization, field, etc at the next step.
    */

    // Calculate dy(y,n)
    D(a, b, n, dt);
    D(b, a, n, dt);

    // Save the values of dy(y,n)
    adx = a->dx; ady = a->dy; adz = a->dz;
    bdx = b->dx; bdy = b->dy; bdz = b->dz;

    // Store the intermediate value yi at n+1
    a->x[n+1] = a->x[n] + adx; a->y[n+1] = a->y[n] + ady; a->z[n+1] = a->z[n] + adz;
    b->x[n+1] = b->x[n] + bdx; b->y[n+1] = b->y[n] + bdy; b->z[n+1] = b->z[n] + bdz;
    
    // Calculate dy(yi,n+1)
    D(a, b, n+1, dt);
    D(b, a, n+1, dt);

    // Get the Heun step
    a->x[n+1] = a->x[n] + 0.5*(adx + a->dx);
    a->y[n+1] = a->y[n] + 0.5*(ady + a->dy);
    a->z[n+1] = a->z[n] + 0.5*(adz + a->dz);
    b->x[n+1] = b->x[n] + 0.5*(bdx + b->dx);
    b->y[n+1] = b->y[n] + 0.5*(bdy + b->dy);
    b->z[n+1] = b->z[n] + 0.5*(bdz + b->dz);

    // Normalize the new magnetization using the Taylor expansion of sqrt() near 1 to speed up the calculation.
    double norminator;
    
    // Domain a
    norminator = 1.0/(1.0 + 0.5 * (a->x[n+1]*a->x[n+1] + a->y[n+1]*a->y[n+1] + a->z[n+1]*a->z[n+1] - 1.0) );
    a->x[n+1] *= a->x[n+1]*norminator;
    a->y[n+1] *= a->y[n+1]*norminator;
    a->z[n+1] *= a->z[n+1]*norminator;

    // Domain b
    norminator = 1.0/(1.0 + 0.5 * (b->x[n+1]*b->x[n+1] + b->y[n+1]*b->y[n+1] + b->z[n+1]*b->z[n+1] - 1.0) );
    b->x[n+1] *= b->x[n+1]*norminator;
    b->y[n+1] *= b->y[n+1]*norminator;
    b->z[n+1] *= b->z[n+1]*norminator;

  } // End of for loop.

  // At this point, the whole solution arrays should be populated.
}
