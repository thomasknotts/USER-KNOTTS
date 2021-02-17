/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/*
CHANGES THAT HAVE BEEN MADE
-changed keyword in FixStyle to be wall/hydro
-Changed Defined the variable in all instances to be FixWallHydro
-Changed the ifndef and define to be _HYDRO_H

-Changed all the FixWall to FixHwall
-Changed all WALL_HYDRO to HWALL_HYDRO
-Changed include to fix_hwall

-Added in private  **coefficients and **offset, similar to the pair_eten
   -thus also created int allocated
-Changed the old coefficient variables so we don't have naming issues
   (implementation algorithm is still lj12-6)
-Added the deconstructor virtual function

-Removed any un-necessary lines
*/

#ifdef FIX_CLASS

FixStyle(wall/hydro,FixHwallHydro)

#else

#ifndef LMP_FIX_HWALL_HYDRO_H
#define LMP_FIX_HWALL_HYDRO_H

#include "fix_hwall.h"

namespace LAMMPS_NS {

class FixHwallHydro : public FixHwall {
 public:
  FixHwallHydro(class LAMMPS *, int, char **);
  virtual ~FixHwallHydro();
  void precompute(int);
  void wall_particle(int, int, double);
  int allocated; 					    // 0/1 = whether arrays are allocated

 private:
  double **coeff1,**coeff2,**coeff3,**coeff4,**coeff5,**coeff6,**offset;  //coeff 1,2,3 -> force, coeff 4,5,6 -> energy, cutoff offset 
                                                            //The array of 6 indicates all possible surface faces:
                                                            //   xlo, xhi, ylo, yhi, zlo, zhi

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Particle on or inside fix wall surface

Particles must be "exterior" to the wall in order for energy/force to
be calculated.

*/



/*
 ADDISON THINKINHG!!!

I think the actual values that are pulled (face, coord, epsilon, sigma, and cutoff) are assigned
in a different routine. To add in additional function variables (H-Index for surface and H-Index for AA)
I will need to modify that one. See atom.h
-Look into adding in a new keyword???

*/


