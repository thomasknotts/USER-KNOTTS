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

/*
Changes
-Changed all FixWall the FixHwall
-Changed the include to be fix_hwall.h
---WORKING up to this point---WORKING up to this point

-Added all the surface hydropoty index "chi" functionality
    -See comments below functions changed: the constructor, the destructor, init(), post_force()
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_hwall.h"
#include "atom.h"
#include "input.h"
#include "variable.h"
#include "domain.h"
#include "lattice.h"
#include "update.h"
#include "modify.h"
#include "respa.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{XLO=0,XHI=1,YLO=2,YHI=3,ZLO=4,ZHI=5};
enum{NONE=0,EDGE,CONSTANT,VARIABLE};

/* ---------------------------------------------------------------------- */

FixHwall::FixHwall(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nwall(0)
{
  scalar_flag = 1;
  vector_flag = 1;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  // parse args

  int scaleflag = 1;
  fldflag = 0;
  int pbcflag = 0;

  for (int i = 0; i < 6; i++) xstr[i] = estr[i] = sstr[i] = cstr[i] = NULL; 	//Added cstr[i]

  int iarg = 3;
  while (iarg < narg) {
    if ((strcmp(arg[iarg],"xlo") == 0) || (strcmp(arg[iarg],"xhi") == 0) ||
        (strcmp(arg[iarg],"ylo") == 0) || (strcmp(arg[iarg],"yhi") == 0) ||
        (strcmp(arg[iarg],"zlo") == 0) || (strcmp(arg[iarg],"zhi") == 0)) {
      if (iarg+6 > narg) error->all(FLERR,"Illegal fix wall command");		//increased the argument size to 6

      int newwall;
      if      (strcmp(arg[iarg],"xlo") == 0) newwall = XLO;
      else if (strcmp(arg[iarg],"xhi") == 0) newwall = XHI;
      else if (strcmp(arg[iarg],"ylo") == 0) newwall = YLO;
      else if (strcmp(arg[iarg],"yhi") == 0) newwall = YHI;
      else if (strcmp(arg[iarg],"zlo") == 0) newwall = ZLO;
      else if (strcmp(arg[iarg],"zhi") == 0) newwall = ZHI;

      for (int m = 0; (m < nwall) && (m < 6); m++)
        if (newwall == wallwhich[m])
          error->all(FLERR,"Wall defined twice in fix wall command");

      wallwhich[nwall] = newwall;
      if (strcmp(arg[iarg+1],"EDGE") == 0) {
        xstyle[nwall] = EDGE;
        int dim = wallwhich[nwall] / 2;
        int side = wallwhich[nwall] % 2;
        if (side == 0) coord0[nwall] = domain->boxlo[dim];
        else coord0[nwall] = domain->boxhi[dim];
      } else if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        xstyle[nwall] = VARIABLE;
        int n = strlen(&arg[iarg+1][2]) + 1;
        xstr[nwall] = new char[n];
        strcpy(xstr[nwall],&arg[iarg+1][2]);
      } else {
        xstyle[nwall] = CONSTANT;
        coord0[nwall] = force->numeric(FLERR,arg[iarg+1]);
      }

      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
        int n = strlen(&arg[iarg+2][2]) + 1;
        estr[nwall] = new char[n];
        strcpy(estr[nwall],&arg[iarg+2][2]);
        estyle[nwall] = VARIABLE;
      } else {
        epsilon[nwall] = force->numeric(FLERR,arg[iarg+2]);
        estyle[nwall] = CONSTANT;
      }

      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) {
        int n = strlen(&arg[iarg+3][2]) + 1;
        sstr[nwall] = new char[n];
        strcpy(sstr[nwall],&arg[iarg+3][2]);
        sstyle[nwall] = VARIABLE;
      } else {
        sigma[nwall] = force->numeric(FLERR,arg[iarg+3]);
        sstyle[nwall] = CONSTANT;
      }

      //New Chi addition
      if (strstr(arg[iarg+4],"v_") == arg[iarg+4]) { 	//If you pass-in a variable chi
        int n = strlen(&arg[iarg+4][2]) + 1;
        cstr[nwall] = new char[n];
        strcpy(cstr[nwall],&arg[iarg+4][2]);
        cstyle[nwall] = VARIABLE;			//actual chi definition in function post-force (locallly defined)
      } else {
        chi[nwall] = force->numeric(FLERR,arg[iarg+4]);	//if you pass-in a constant chi (most common)
        cstyle[nwall] = CONSTANT;
      }
      ////

      cutoff[nwall] = force->numeric(FLERR,arg[iarg+5]);		//changed to iarg+5
      nwall++;
      iarg += 6;							//Changed to iarg+=6

    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix wall command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"fld") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall command");
      if (strcmp(arg[iarg+1],"no") == 0) fldflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) fldflag = 1;
      else error->all(FLERR,"Illegal fix wall command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"pbc") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall command");
      if (strcmp(arg[iarg+1],"yes") == 0) pbcflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) pbcflag = 0;
      else error->all(FLERR,"Illegal fix wall command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix wall command");
  }

  size_vector = nwall;

  // error checks

  if (nwall == 0) error->all(FLERR,"Illegal fix wall command");
  for (int m = 0; m < nwall; m++)
    if (cutoff[m] <= 0.0)
      error->all(FLERR,"Fix wall cutoff <= 0.0");

  for (int m = 0; m < nwall; m++)
    if ((wallwhich[m] == ZLO || wallwhich[m] == ZHI) && domain->dimension == 2)
      error->all(FLERR,"Cannot use fix wall zlo/zhi for a 2d simulation");

  if (!pbcflag) {
    for (int m = 0; m < nwall; m++) {
      if ((wallwhich[m] == XLO || wallwhich[m] == XHI) && domain->xperiodic)
        error->all(FLERR,"Cannot use fix wall in periodic dimension");
      if ((wallwhich[m] == YLO || wallwhich[m] == YHI) && domain->yperiodic)
        error->all(FLERR,"Cannot use fix wall in periodic dimension");
      if ((wallwhich[m] == ZLO || wallwhich[m] == ZHI) && domain->zperiodic)
        error->all(FLERR,"Cannot use fix wall in periodic dimension");
    }
  }

  // scale factors for wall position for CONSTANT and VARIABLE walls

  int flag = 0;
  for (int m = 0; m < nwall; m++)
    if (xstyle[m] != EDGE) flag = 1;

  if (flag) {
    if (scaleflag) {
      xscale = domain->lattice->xlattice;
      yscale = domain->lattice->ylattice;
      zscale = domain->lattice->zlattice;
    }
    else xscale = yscale = zscale = 1.0;

    for (int m = 0; m < nwall; m++) {
      if (xstyle[m] != CONSTANT) continue;
      if (wallwhich[m] < YLO) coord0[m] *= xscale;
      else if (wallwhich[m] < ZLO) coord0[m] *= yscale;
      else coord0[m] *= zscale;
    }
  }

  // set xflag if any wall positions are variable
  // set varflag if any wall positions or parameters are variable
  // set wstyle to VARIABLE if either epsilon, sigma or chi is a variable

  varflag = xflag = 0;
  for (int m = 0; m < nwall; m++) {
    if (xstyle[m] == VARIABLE) xflag = 1;
    if (xflag || estyle[m] == VARIABLE || sstyle[m] == VARIABLE || cstyle[m] == VARIABLE) varflag = 1;	//Added cstyle to this check
    if (estyle[m] == VARIABLE || sstyle[m] == VARIABLE || cstyle[m] == VARIABLE) wstyle[m] = VARIABLE;	//Added cstyle to this check
    else wstyle[m] = CONSTANT;
  }

  eflag = 0;
  for (int m = 0; m <= nwall; m++) ewall[m] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixHwall::~FixHwall()
{
  for (int m = 0; m < nwall; m++) {
    delete [] xstr[m];
    delete [] estr[m];
    delete [] sstr[m];
    delete [] cstr[m]; 		//Added cstr to this release
  }
}

/* ---------------------------------------------------------------------- */

int FixHwall::setmask()
{
  int mask = 0;

  // FLD implicit needs to invoke wall forces before pair style

  if (fldflag) mask |= PRE_FORCE;
  else mask |= POST_FORCE;

  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHwall::init()
{
  for (int m = 0; m < nwall; m++) {
    if (xstyle[m] == VARIABLE) {
      xindex[m] = input->variable->find(xstr[m]);
      if (xindex[m] < 0)
        error->all(FLERR,"Variable name for fix wall does not exist");
      if (!input->variable->equalstyle(xindex[m]))
        error->all(FLERR,"Variable for fix wall is invalid style");
    }
    if (estyle[m] == VARIABLE) {
      eindex[m] = input->variable->find(estr[m]);
      if (eindex[m] < 0)
        error->all(FLERR,"Variable name for fix wall does not exist");
      if (!input->variable->equalstyle(eindex[m]))
        error->all(FLERR,"Variable for fix wall is invalid style");
    }
    if (sstyle[m] == VARIABLE) {
      sindex[m] = input->variable->find(sstr[m]);
      if (sindex[m] < 0)
        error->all(FLERR,"Variable name for fix wall does not exist");
      if (!input->variable->equalstyle(sindex[m]))
        error->all(FLERR,"Variable for fix wall is invalid style");
    }
    if (cstyle[m] == VARIABLE) {						//Added this VARIABLE cstyle set and check
      cindex[m] = input->variable->find(cstr[m]);
      if (cindex[m] < 0)
        error->all(FLERR,"Variable name for fix wall does not exist");
      if (!input->variable->equalstyle(cindex[m]))
        error->all(FLERR,"Variable for fix wall is invalid style");
    }

  }

  // setup coefficients

  for (int m = 0; m < nwall; m++) precompute(m);				//This is where my custom Hwall_hydro.cpp is called

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixHwall::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet")) {
    if (!fldflag) post_force(vflag);
  } else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixHwall::min_setup(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   only called if fldflag set, in place of post_force
------------------------------------------------------------------------- */

void FixHwall::pre_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixHwall::post_force(int vflag)
{
  eflag = 0;
  for (int m = 0; m <= nwall; m++) ewall[m] = 0.0;

  // coord = current position of wall
  // evaluate variables if necessary, wrap with clear/add
  // for epsilon/sigma/chi variables need to re-invoke precompute()
  // NOTE defines eps/sig/chi values if any are VARIABLE

  if (varflag) modify->clearstep_compute();

  double coord;
  for (int m = 0; m < nwall; m++) {
    if (xstyle[m] == VARIABLE) {
      coord = input->variable->compute_equal(xindex[m]);
      if (wallwhich[m] < YLO) coord *= xscale;
      else if (wallwhich[m] < ZLO) coord *= yscale;
      else coord *= zscale;
    } else coord = coord0[m];
    if (wstyle[m] == VARIABLE) {
      if (estyle[m] == VARIABLE) {
        epsilon[m] = input->variable->compute_equal(eindex[m]);
        if (epsilon[m] < 0.0)
          error->all(FLERR,"Variable evaluation in fix wall gave bad value");
      }
      if (sstyle[m] == VARIABLE) {
        sigma[m] = input->variable->compute_equal(sindex[m]);
        if (sigma[m] < 0.0)
          error->all(FLERR,"Variable evaluation in fix wall gave bad value");
      }
      if (cstyle[m] == VARIABLE) { 						//Added chi functionality
        chi[m] = input->variable->compute_equal(cindex[m]);
        if (chi[m] < 0.0)
          error->all(FLERR,"Variable evaluation in fix wall gave bad value");
      }
      precompute(m); 								//Note: Re-call precompute
    }

    wall_particle(m,wallwhich[m],coord);
  }

  if (varflag) modify->addstep_compute(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

void FixHwall::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixHwall::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of wall interaction
------------------------------------------------------------------------- */

double FixHwall::compute_scalar()
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(ewall,ewall_all,nwall+1,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return ewall_all[0];
}

/* ----------------------------------------------------------------------
   components of force on wall
------------------------------------------------------------------------- */

double FixHwall::compute_vector(int n)
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(ewall,ewall_all,nwall+1,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return ewall_all[n+1];
}
