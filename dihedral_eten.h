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

/* ----------------------------------------------------------------------
   Written by Thomas Allen Knotts IV                                     
   Based on dihedral_charmm.h                                          
------------------------------------------------------------------------- */

#ifdef DIHEDRAL_CLASS

DihedralStyle(eten,DihedralEten)

#else

#ifndef LMP_DIHEDRAL_ETEN_H
#define LMP_DIHEDRAL_ETEN_H

#include <stdio.h>
#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralEten : public Dihedral {
 public:
  DihedralEten(class LAMMPS *);
  virtual ~DihedralEten();
  virtual void compute(int, int);
  virtual void coeff(int, char **);
  virtual void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);

 protected:
  double *k,*weight,*cos_shift,*sin_shift, *shift;
  int *multiplicity;
  double **lj14_1,**lj14_2,**lj14_3,**lj14_4;
  int implicit,weightflag;

  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

W: Dihedral problem: %d %ld %d %d %d %d

Conformation of the 4 listed dihedral atoms is extreme; you may want
to check your simulation geometry.

E: Incorrect args for dihedral coefficients

Self-explanatory.  Check the input script or data file.

E: Incorrect multiplicity arg for dihedral coefficients

Self-explanatory.  Check the input script or data file.

E: Incorrect weight arg for dihedral coefficients

Self-explanatory.  Check the input script or data file.

E: Dihedral eten is incompatible with Pair style

Dihedral style eten must be used with a pair style eten
in order for the 1-4 epsilon/sigma parameters to be defined.

*/