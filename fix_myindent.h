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

#ifdef FIX_CLASS

FixStyle(myindent,FixMyindent)

#else

#ifndef LMP_FIX_MYINDENT_H
#define LMP_FIX_MYINDENT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMyindent : public Fix {
 public:
  FixMyindent(class LAMMPS *, int, char **);
  ~FixMyindent();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

 private:
  int istyle,scaleflag,side;
  double k,k3;
  char *xstr,*ystr,*zstr,*rstr,*pstr;
  int xvar,yvar,zvar,rvar,pvar;
  
  // Beibei Wang, cone shape indenter
  double xvalue,yvalue,zvalue,rvalue,theta,htip;

  //Yuzi He, Berkovich indenter
  //h_sharp is defined as the tangent of the angle A
  //angle A is between horizontal and edge of indenter
  double h_sharp;

  //at the tip of the indenter, it is treated as a flat head
  //h_head is the distance from tip to the flat head.
  double h_head;
  

  int indenter_flag;
  double indenter[4],indenter_all[4];
  int cdim,varflag;
  int ilevel_respa;

  void options(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Variable name for fix indent does not exist

Self-explanatory.

E: Variable for fix indent is invalid style

Only equal-style variables can be used.

E: Variable for fix indent is not equal style

Only equal-style variables can be used.

*/
