/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <math.h>
#include "fix_hwall_hydro.h"
#include "atom.h"
#include "error.h"
#include "memory.h"
//#include "math_const.h"  //Not working for some reason...

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixHwallHydro::FixHwallHydro(LAMMPS *lmp, int narg, char **arg) :
  FixHwall(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

FixHwallHydro::~FixHwallHydro()
{
  if (allocated) {
    memory->destroy(coeff1);
    memory->destroy(coeff2);
    memory->destroy(coeff3);
    memory->destroy(coeff4);
    memory->destroy(coeff5);
    memory->destroy(coeff6);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void FixHwallHydro::precompute(int m) 
{
  int nlocal = atom->nlocal;
  double roh = 1.0;
  double MY_PI = 3.14159265358979323846;
  double theta1 = 0.2340000000;
  double theta2 = 0.4962857143; 
  double theta3 = 0.1333333333;
  double thetas = 0.01;  
  double thetap = 0.05; 

  allocated = 1;
  memory->create(coeff1,6,nlocal,"fix:hwall/coeff1");
  memory->create(coeff2,6,nlocal,"fix:hwall/coeff2");
  memory->create(coeff3,6,nlocal,"fix:hwall/coeff3");
  memory->create(coeff4,6,nlocal,"fix:hwall/coeff4");
  memory->create(coeff5,6,nlocal,"fix:hwall/coeff5");
  memory->create(coeff6,6,nlocal,"fix:hwall/coeff6");
  memory->create(offset,6,nlocal,"fix:hwall/offset");

  int flag = 1;
  int aaindex = atom->find_custom("AAcid", flag);
  int *aa = atom->ivector[aaindex];
  int aaeindex = atom->find_custom("wc_epsilon", flag);
  double *aaepsilon = atom->dvector[aaeindex];
  int aasindex = atom->find_custom("wc_sigma", flag);
  double *aasigma = atom->dvector[aasindex];

  double eprime, sprime, aachi;
  double coeff0;
  double rinv = 1.0/cutoff[m];
  double r3inv = rinv * rinv * rinv;
  double r4inv = r3inv * rinv;
  double r6inv = r3inv * r3inv;
  for (int i = 0; i < nlocal; i++) {
    eprime = sqrt(epsilon[m] * aaepsilon[i]);
    sprime = 0.5 * (sigma[m] + aasigma[i]);
    if      (aa[i] == 1)  aachi =  1.8;  //Alanine
    else if (aa[i] == 18) aachi = -4.5;  //Arginine
    else if (aa[i] == 14) aachi = -3.5;  //Asparagine
    else if (aa[i] == 4)  aachi = -3.5;  //Aspartic Acid
    else if (aa[i] == 3)  aachi =  2.5;  //Cysteine
    else if (aa[i] == 5)  aachi = -3.5;  //Glutamic Acid
    else if (aa[i] == 17) aachi = -3.5;  //Glutamine
    else if (aa[i] == 7)  aachi = -0.4;  //Glycine
    else if (aa[i] == 8)  aachi = -3.2;  //Histidine
    else if (aa[i] == 9)  aachi =  4.5;  //Isoleucine
    else if (aa[i] == 12) aachi =  3.8;  //Leucine
    else if (aa[i] == 11) aachi = -3.9;  //Lysine
    else if (aa[i] == 13) aachi =  1.9;  //Methionine
    else if (aa[i] == 6)  aachi =  2.8;  //Phenylalanine
    else if (aa[i] == 16) aachi =  1.6;  //Proline
    else if (aa[i] == 19) aachi = -0.8;  //Serine
    else if (aa[i] == 20) aachi = -0.7;  //Threonine
    else if (aa[i] == 23) aachi = -0.9;  //Tryptophan
    else if (aa[i] == 25) aachi = -1.3;  //Tyrosine
    else if (aa[i] == 22) aachi =  4.2;  //Valine
    else error->one(FLERR, "Improper amino acid identifier");
    coeff0 = MY_PI * roh * sprime * sprime * eprime;
    coeff1[m][i] = coeff0 * 9 * theta1 * pow(sprime,10);
    coeff2[m][i] = coeff0 * 7 * theta2 * pow(sprime,8);
    coeff3[m][i] = coeff0 * 3 * (theta3 - thetas*(chi[m]-4.5) - thetap*aachi) * pow(sprime,4);
    coeff0 = coeff0 * sprime;
    coeff4[m][i] = coeff0 * theta1 * pow(sprime,9);
    coeff5[m][i] = coeff0 * theta2 * pow(sprime,7);
    coeff6[m][i] = coeff0 * (theta3 - thetas*(chi[m]-4.5) - thetap*aachi) * pow(sprime,3);
    offset[m][i] = r3inv * (coeff4[m][i]*r6inv - coeff5[m][i]*r4inv + coeff6[m][i]);
  }
}

/* ----------------------------------------------------------------------
   interaction of all particles in group with a wall
   m = index of wall coeffs
   which = xlo,xhi,ylo,yhi,zlo,zhi
   error if any particle is on or behind wall
------------------------------------------------------------------------- */

void FixHwallHydro::wall_particle(int m, int which, double coord)
{
  double delta,rinv,r3inv,r4inv,r6inv,fwall;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  int dim = which / 2;
  int side = which % 2;
  if (side == 0) side = -1;

  int onflag = 0;
 
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (side < 0) delta = x[i][dim] - coord;
      else delta = coord - x[i][dim];
      if (delta >= cutoff[m]) continue;
      if (delta <= 0.0) {
        onflag = 1;
        continue;
      }
      rinv = 1.0/delta;
      r3inv = rinv*rinv*rinv;
      r4inv = r3inv*rinv;
      r6inv = r3inv*r3inv;
      fwall = side * r4inv * (coeff1[m][i]*r6inv - coeff2[m][i]*r4inv + coeff3[m][i]);
      f[i][dim] -= fwall;
      ewall[0] += r3inv * (coeff4[m][i]*r6inv - coeff5[m][i]*r4inv + coeff6[m][i]) - offset[m][i];
      ewall[m+1] += fwall;
    }

  if (onflag) error->one(FLERR,"Particle on or inside fix wall surface");
}
