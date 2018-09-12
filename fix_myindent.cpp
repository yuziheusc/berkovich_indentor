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

/* ----------------------------------------------------------------------
   Contributing author: Ravi Agrawal (Northwestern U)
   ------------------------------------------------------------------------- */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fix_myindent.h"
#include "atom.h"
#include "input.h"
#include "variable.h"
#include "domain.h"
#include "lattice.h"
#include "update.h"
#include "modify.h"
#include "output.h"
#include "respa.h"
#include "error.h"
#include "force.h"

double frad = 3.14159265/180.0;

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,SPHERE,CONE,BKVCH};
enum{INSIDE,OUTSIDE};

/* ---------------------------------------------------------------------- */

FixMyindent::FixMyindent(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  xstr(NULL), ystr(NULL), zstr(NULL), rstr(NULL), pstr(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal fix indent command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  k = force->numeric(FLERR,arg[3]);
  k3 = k/3.0;

  // read options from end of input line

  options(narg-4,&arg[4]);

  // setup scaling

  double xscale,yscale,zscale;
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // apply scaling factors to geometry

  if (istyle == SPHERE || istyle == CONE || istyle == BKVCH) {
    if (!xstr) xvalue *= xscale;
    if (!ystr) yvalue *= yscale;
    if (!zstr) zvalue *= zscale;
    if (!rstr) rvalue *= xscale;
  } else error->all(FLERR,"Illegal fix indent command");

  varflag = 0;
  if (xstr || ystr || zstr || rstr || pstr) varflag = 1;

  indenter_flag = 0;
  indenter[0] = indenter[1] = indenter[2] = indenter[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixMyindent::~FixMyindent()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] rstr;
  delete [] pstr;
}

/* ---------------------------------------------------------------------- */

int FixMyindent::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMyindent::init()
{
  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix indent does not exist");
    if (!input->variable->equalstyle(xvar))
      error->all(FLERR,"Variable for fix indent is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix indent does not exist");
    if (!input->variable->equalstyle(yvar))
      error->all(FLERR,"Variable for fix indent is not equal style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix indent does not exist");
    if (!input->variable->equalstyle(zvar))
      error->all(FLERR,"Variable for fix indent is not equal style");
  }
  if (rstr) {
    rvar = input->variable->find(rstr);
    if (rvar < 0)
      error->all(FLERR,"Variable name for fix indent does not exist");
    if (!input->variable->equalstyle(rvar))
      error->all(FLERR,"Variable for fix indent is not equal style");
  }
  if (pstr) {
    pvar = input->variable->find(pstr);
    if (pvar < 0)
      error->all(FLERR,"Variable name for fix indent does not exist");
    if (!input->variable->equalstyle(pvar))
      error->all(FLERR,"Variable for fix indent is not equal style");
  }

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixMyindent::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixMyindent::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMyindent::post_force(int vflag)
{
  // indenter values, 0 = energy, 1-3 = force components
  // wrap variable evaluations with clear/add

  if (varflag) modify->clearstep_compute();

  indenter_flag = 0;
  indenter[0] = indenter[1] = indenter[2] = indenter[3] = 0.0;

  // spherical indenter

  if (istyle == SPHERE) {

    // ctr = current indenter center
    // remap into periodic box

    double ctr[3];
    if (xstr) ctr[0] = input->variable->compute_equal(xvar);
    else ctr[0] = xvalue;
    if (ystr) ctr[1] = input->variable->compute_equal(yvar);
    else ctr[1] = yvalue;
    if (zstr) ctr[2] = input->variable->compute_equal(zvar);
    else ctr[2] = zvalue;
    domain->remap(ctr);

    double radius;
    if (rstr) radius = input->variable->compute_equal(rvar);
    else radius = rvalue;

    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double delx,dely,delz,r,dr,fmag,fx,fy,fz;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        delx = x[i][0] - ctr[0];
        dely = x[i][1] - ctr[1];
        delz = x[i][2] - ctr[2];
        domain->minimum_image(delx,dely,delz);
        r = sqrt(delx*delx + dely*dely + delz*delz);
        if (side == OUTSIDE) {
          dr = r - radius;
          fmag = k*dr*dr;
        } else {
          dr = radius - r;
          fmag = -k*dr*dr;
        }
        if (dr >= 0.0) continue;
        fx = delx*fmag/r;
        fy = dely*fmag/r;
        fz = delz*fmag/r;
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        indenter[0] -= k3 * dr*dr*dr;
        indenter[1] -= fx;
        indenter[2] -= fy;
        indenter[3] -= fz;
      }

    // conical indenter (with round tip, cone angle cannot be changed once initilized)

  } else if (istyle == CONE) {

    // ctr = current indenter center
    // remap into periodic box

    double ctr[3];
    if (xstr) ctr[0] = input->variable->compute_equal(xvar);
    else ctr[0] = xvalue;
    if (ystr) ctr[1] = input->variable->compute_equal(yvar);
    else ctr[1] = yvalue;
    if (zstr) ctr[2] = input->variable->compute_equal(zvar);
    else ctr[2] = zvalue;
    domain->remap(ctr);

    double radius;
    if (rstr) radius = input->variable->compute_equal(rvar);
    else radius = rvalue;

    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double delx,dely,delz,r,raxis,rcone,dr,fmag,fx,fy,fz;
    double rinv, r2inv, r6inv, r12inv, r13inv;
    double fcos = cos(theta*0.5*frad), fsin = sin(theta*0.5*frad), ftan = tan(theta*0.5*frad), vcos;
    // the interaction betwen the indenter and MoS2 is lj repulsive part:
    // type list: 1: Mo, 2: Sulfur, energy unit is eV, distance is Angstrom
    double eps = 0.007355, sig= 3.219, rcut = 10.0, rcut_cone = rcut/fcos;
    double ffac = 48*eps*pow(sig, 12);
    double efac = ffac/12.0;
    double E_corr = -1*efac/pow(rcut,12);
    double ztb = ctr[2] - radius/fsin, zt0 = ctr[2] - radius - rcut, zt1 = ctr[2] - radius;
    double zt2 = ctr[2] - radius*fsin, zt3 = ctr[2] - radius + htip;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (type[i] == 1) {
          eps = 0.003325; sig = 2.818;
          ffac = 48*eps*pow(sig,12); efac = ffac/12.0; E_corr = -1*efac/pow(rcut,12);
        }
        delx = x[i][0] - ctr[0];
        dely = x[i][1] - ctr[1];
        delz = x[i][2] - ctr[2];
        domain->minimum_image(delx,dely,delz);

	// bellow the indenter tip or above the whole indenter:
        if (x[i][2] <= zt0 or x[i][2] >= zt3) continue; 
	
	else {
          r = sqrt(delx*delx + dely*dely + delz*delz);
          vcos = -1*delz/r;
          if (vcos<=1 && vcos>= 0.5) {
	    // in the spherical region:
            dr = r - radius - rcut;
	    if (dr >= 0.0) continue; // out of the sphere tip
            else {
	      rinv = 1.0/(rcut+dr);
	      r2inv = rinv*rinv; r6inv = r2inv*r2inv*r2inv; r12inv = r6inv*r6inv; r13inv = r12inv*rinv;
	      fmag = ffac*r13inv;
              fx = delx*fmag/r;
              fy = dely*fmag/r;
              fz = delz*fmag/r;
              f[i][0] += fx;
              f[i][1] += fy;
              f[i][2] += fz;
              indenter[0] += efac*r12inv + E_corr;
              indenter[1] -= fx;
              indenter[2] -= fy;
              indenter[3] -= fz;
            }
	  } else {
            // in the conical region:
	    raxis = sqrt(delx*delx + dely*dely);
	    rcone = (x[i][2]-ztb)*ftan;
            if(rcone < 0) continue;
            dr = raxis - rcone - rcut_cone;
	    if(dr >= 0) continue;
	    else {
	      rinv = 1.0/(rcut+dr*fcos);
	      r2inv = rinv*rinv; r6inv = r2inv*r2inv*r2inv; r12inv = r6inv*r6inv; r13inv = r12inv*rinv;
	      fmag = ffac*r13inv;
              fx = delx*fmag*fcos/raxis;
              fy = dely*fmag*fcos/raxis;
              fz = -1*fmag*fsin;
              f[i][0] += fx;
              f[i][1] += fy;
              f[i][2] += fz;
              indenter[0] += efac * r12inv + E_corr;
              indenter[1] -= fx;
              indenter[2] -= fy;
              indenter[3] -= fz;
	    }
          }
	}
      }
  }else if(istyle == BKVCH) {
    // force calculation here
    const double PI = 4.*atan(1.0);
    const double ETA = 1e-8;
      
    double ctr[3];
    if (xstr) ctr[0] = input->variable->compute_equal(xvar);
    else ctr[0] = xvalue;
    if (ystr) ctr[1] = input->variable->compute_equal(yvar);
    else ctr[1] = yvalue;
    if (zstr) ctr[2] = input->variable->compute_equal(zvar);
    else ctr[2] = zvalue;
    domain->remap(ctr);

    //printf("center %.2lf %.2lf %.2lf\n", ctr[0], ctr[1], ctr[2]);
    
    //first define three vectors of the berkovich indenter
    double vec1[3],vec2[3],vec3[3], lvec;
    vec1[0] = 0; vec1[1] = -1; vec1[2] = h_sharp;
    vec2[0] =  sqrt(3.0)/2; vec2[1] = 0.5; vec2[2] = h_sharp;
    vec3[0] = -sqrt(3.0)/2; vec3[1] = 0.5; vec3[2] = h_sharp;
    lvec = sqrt(1.0 + h_sharp*h_sharp);
    
    //displace the three vectors to the center of the indenter
    //here center measn the tip
    //vec1[0] += ctr[0]; vec1[1] += ctr[1]; vec1[2] += ctr[2];
    //vec2[0] += ctr[0]; vec2[1] += ctr[1]; vec2[2] += ctr[2];
    //vec3[0] += ctr[0]; vec3[1] += ctr[1]; vec3[2] += ctr[2];
        
    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    
    double delx,dely,delz,fmag,fx,fy,fz;
    
    // mark the compress plane of indenter for a specific particle 
    double *vec_a,*vec_b;
    
    for(int i=0; i < nlocal; ++i)
      if (mask[i] & groupbit){
        // calculate the displacement of a atom to the tip
        delx = x[i][0] - ctr[0];
        dely = x[i][1] - ctr[1];
        delz = x[i][2] - ctr[2];
        domain->minimum_image(delx,dely,delz);

        //calculate the domian of the particle
        //decide which plan will interact with the particle
        double z_ang = atan(dely/(delx+ETA));
        if(delx < 0.0) z_ang = PI + z_ang;

        //get the vector which shows the plane
        if(z_ang>=-PI/2 && z_ang<PI/6){vec_a = vec2;vec_b = vec1;}
        else if(z_ang>=PI/6 && z_ang<PI*5/6){vec_a = vec3;vec_b = vec2;}
        else {vec_a = vec1;vec_b = vec3;}

        //calculate the distance into the plane

        //1. calculate the norm
        double vec_n[3],l_vec_n;
        vec_n[0] =   vec_a[1]*vec_b[2] - vec_a[2]*vec_b[1];
        vec_n[1] = -(vec_a[0]*vec_b[2] - vec_a[2]*vec_b[0]);
        vec_n[2] =   vec_a[0]*vec_b[1] - vec_a[1]*vec_b[0];
        l_vec_n = sqrt(vec_n[0]*vec_n[0]+vec_n[1]*vec_n[1]+vec_n[2]*vec_n[2]);

        //2. calculate the distance of particle into the plane
        double d_into = (delx*vec_n[0]+dely*vec_n[1]+delz*vec_n[2])/l_vec_n;

        //3. use F(r) = - K d^2 to calculate force 
        if(d_into > 0.0) d_into = 0.0;

        fmag = k * d_into*d_into; 
                
        fx = fmag*vec_n[0]; fy = fmag*vec_n[1]; fz = fmag*vec_n[2];
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        indenter[0] -= k3 * d_into * d_into * d_into;
        indenter[1] -= fx;
        indenter[2] -= fy;
        indenter[3] -= fz;
    
      }// end if masgroupbit
  }

  if (varflag) modify->addstep_compute(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

void FixMyindent::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMyindent::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of indenter interaction
   ------------------------------------------------------------------------- */

double FixMyindent::compute_scalar()
{
  // only sum across procs one time

  if (indenter_flag == 0) {
    MPI_Allreduce(indenter,indenter_all,4,MPI_DOUBLE,MPI_SUM,world);
    indenter_flag = 1;
  }
  return indenter_all[0];
}

/* ----------------------------------------------------------------------
   components of force on indenter
   ------------------------------------------------------------------------- */

double FixMyindent::compute_vector(int n)
{
  // only sum across procs one time

  if (indenter_flag == 0) {
    MPI_Allreduce(indenter,indenter_all,4,MPI_DOUBLE,MPI_SUM,world);
    indenter_flag = 1;
  }
  return indenter_all[n+1];
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
   ------------------------------------------------------------------------- */

void FixMyindent::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal fix indent command");

  istyle = NONE;
  xstr = ystr = zstr = rstr = pstr = NULL;
  xvalue = yvalue = zvalue = rvalue = 0.0;
  scaleflag = 1;
  side = OUTSIDE;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"sphere") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix indent command");

      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        xstr = new char[n];
        strcpy(xstr,&arg[iarg+1][2]);
      } else xvalue = force->numeric(FLERR,arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
        int n = strlen(&arg[iarg+2][2]) + 1;
        ystr = new char[n];
        strcpy(ystr,&arg[iarg+2][2]);
      } else yvalue = force->numeric(FLERR,arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) {
        int n = strlen(&arg[iarg+3][2]) + 1;
        zstr = new char[n];
        strcpy(zstr,&arg[iarg+3][2]);
      } else zvalue = force->numeric(FLERR,arg[iarg+3]);
      if (strstr(arg[iarg+4],"v_") == arg[iarg+4]) {
        int n = strlen(&arg[iarg+4][2]) + 1;
        rstr = new char[n];
        strcpy(rstr,&arg[iarg+4][2]);
      } else rvalue = force->numeric(FLERR,arg[iarg+4]);

      istyle = SPHERE;
      iarg += 5;

    } else if (strcmp(arg[iarg],"cone") == 0){
      if (iarg+7 > narg) error->all(FLERR,"Illegal fix indent command");

      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        xstr = new char[n];
        strcpy(xstr,&arg[iarg+1][2]);
      } else xvalue = force->numeric(FLERR,arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
        int n = strlen(&arg[iarg+2][2]) + 1;
        ystr = new char[n];
        strcpy(ystr,&arg[iarg+2][2]);
      } else yvalue = force->numeric(FLERR,arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) {
        int n = strlen(&arg[iarg+3][2]) + 1;
        zstr = new char[n];
        strcpy(zstr,&arg[iarg+3][2]);
      } else zvalue = force->numeric(FLERR,arg[iarg+3]);

      rvalue = force->numeric(FLERR,arg[iarg+4]);
      theta = force->numeric(FLERR,arg[iarg+5]);
      htip = force->numeric(FLERR,arg[iarg+6]);

      istyle = CONE;
      iarg += 7;

    
    }else if (strcmp(arg[iarg],"bkvch") == 0){
      if (iarg+6 > narg) error->all(FLERR,"Illegal fix indent command");
      if(strstr(arg[iarg+1],"v_") == arg[iarg+1]){
        int n = strlen(&arg[iarg+1][2]) + 1;
        xstr = new char[n];
        strcpy(xstr,&arg[iarg+1][2]);
      }else xvalue = force->numeric(FLERR,arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
        int n = strlen(&arg[iarg+2][2]) + 1;
        ystr = new char[n];
        strcpy(ystr,&arg[iarg+2][2]);
      } else yvalue = force->numeric(FLERR,arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) {
        int n = strlen(&arg[iarg+3][2]) + 1;
        zstr = new char[n];
        strcpy(zstr,&arg[iarg+3][2]);
      } else zvalue = force->numeric(FLERR,arg[iarg+3]);
      h_sharp = force->numeric(FLERR, arg[iarg+4]);
      h_head = force->numeric(FLERR, arg[iarg+5]);
      istyle = BKVCH;
      iarg += 6;
      
    }
    else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix indent command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix indent command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"side") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix indent command");
      if (strcmp(arg[iarg+1],"in") == 0) side = INSIDE;
      else if (strcmp(arg[iarg+1],"out") == 0) side = OUTSIDE;
      else error->all(FLERR,"Illegal fix indent command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix indent command");
  }
}
