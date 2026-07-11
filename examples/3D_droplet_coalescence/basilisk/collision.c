/**
# 3D binary droplet collision — Qian & Law case p

Basilisk translation of `examples/3D_droplet_coalescence` case **p**
(tetradecane droplets in nitrogen). This is the near-grazing, high-Weber
collision that sits in the stretching-separation regime:

    We = rho_l*Ur^2*D/sigma = 64.9
    Re = rho_l*Ur*D/mu_l    = 312.8
    B  = impact parameter    = 0.71   (offset / D)
    R  = 177 um,  D = 2R

MFC solves this with an artificial-compressible five-equation model (liquid
sound speed reduced to c_l = 100 m/s so the acoustic dt is affordable). The
true Mach number is ~0.002, so the physics is incompressible. Here we solve it
the way Basilisk is built for it: incompressible two-phase Navier-Stokes with a
VOF interface, geometric surface tension, and octree AMR. Units are SI, so the
time axis (0 .. 3 ms) and lengths line up directly with the MFC output for
cross-comparison.

The momentum-conserving scheme + filtered properties are used because the
density ratio is large (rho_l/rho_g = 666).
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#define FILTERED 1
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "tag.h"

/**
## Parameters

Physical properties (SI), from `case.py` / `case_p_v100.py` for case p. */

#define RHO_L  763.0
#define RHO_G  (763.0/666.0)
#define MU_L   2.183033e-3
#define MU_G   1.834482e-5
#define SIGMA  0.0266
#define RAD    177e-6
#define DIAM   (2.0*RAD)

/** Collision kinematics. Each drop moves at UD toward the other, so the
relative velocity is Ur = 2*UD. Centres are offset by SEP along x and by the
impact parameter B*D across it (split half-and-half in y). */

#define UR     2.528129            // relative velocity  sqrt(We*sigma/(rho_l*D))
#define UD     (UR/2.0)            // per-drop speed
#define SEP    (0.68*DIAM)         // x half-separation of the centres
#define BIMP   (0.71*DIAM)         // impact-parameter offset (across the line of approach)

#define X1  (-SEP)
#define Y1  (-BIMP/2.0)
#define X2  ( SEP)
#define Y2  ( BIMP/2.0)

/** Numerics. Domain is a cube of side L0 centred on the origin; the pair sits
near the middle and the far field is filled cheaply by coarse cells. `MAXLEVEL`
sets the finest resolution: with L0 = 6D, level 10 gives dx_min ~ D/170
(comparable to the MFC production grid, D/147-155) and level 11 gives ~D/341.
Drop to 7-8 for a quick smoke test. Override at compile time with
`-DMAXLEVEL=...`. */

#define L0_D   6.0                 // domain side, in diameters
#ifndef MAXLEVEL
  #define MAXLEVEL 10
#endif
#define MINLEVEL 4
#define FEMAX   1e-3               // VOF adaptivity tolerance
#define UEMAX   (1e-2*UR)          // velocity adaptivity tolerance

/** Time. t* = D/Ur = 0.140 ms; run to 3 ms (~21 t*) and save 300 frames at a
10 us cadence, matching the MFC output. */

#define TEND   3e-3
#define TSAVE  (TEND/300.0)

#define V0     (2.0*4.0/3.0*pi*cube(RAD))   // initial liquid volume (two drops)

scalar f0[];                       // helper for the initial per-drop velocity

int main (int argc, char * argv[])
{
  size (L0_D*DIAM);
  origin (-L0/2., -L0/2., -L0/2.);
  init_grid (1 << MINLEVEL);

  rho1 = RHO_L; rho2 = RHO_G;
  mu1  = MU_L;  mu2  = MU_G;
  f.sigma = SIGMA;

  /* Tighten the projection a bit for the high density ratio. */
  TOLERANCE = 1e-4;

  run();
}

/**
## Boundary conditions

Open far field on all six faces (non-reflecting outflow), the incompressible
analogue of MFC's `bc = -8`. The pair is far from the walls and the net
momentum is zero, so this mostly just pins a reference pressure. */

u.n[left]   = neumann(0.);  p[left]   = dirichlet(0.);  pf[left]   = dirichlet(0.);
u.n[right]  = neumann(0.);  p[right]  = dirichlet(0.);  pf[right]  = dirichlet(0.);
u.n[bottom] = neumann(0.);  p[bottom] = dirichlet(0.);  pf[bottom] = dirichlet(0.);
u.n[top]    = neumann(0.);  p[top]    = dirichlet(0.);  pf[top]    = dirichlet(0.);
u.n[back]   = neumann(0.);  p[back]   = dirichlet(0.);  pf[back]   = dirichlet(0.);
u.n[front]  = neumann(0.);  p[front]  = dirichlet(0.);  pf[front]  = dirichlet(0.);

/**
## Initial condition

Two spheres of radius R, union'd into the VOF field, each carrying its own
approach velocity. The interface is refined to MAXLEVEL before the fraction is
computed so the initial VOF is sharp. */

#define sphere(xc,yc,zc) (RAD - sqrt(sq(x-(xc)) + sq(y-(yc)) + sq(z-(zc))))

event init (t = 0)
{
  if (!restore (file = "restart")) {
    int it = 0;
    do {
      fraction (f, max (sphere(X1,Y1,0.), sphere(X2,Y2,0.)));
    } while (++it < 10 &&
             adapt_wavelet ({f}, (double[]){FEMAX}, MAXLEVEL, MINLEVEL).nf);

    fraction (f, max (sphere(X1,Y1,0.), sphere(X2,Y2,0.)));

    /* Per-drop velocity: whichever sphere a cell is closer to sets its sign;
       weighted by f so the gas stays quiescent. */
    foreach() {
      double s1 = sphere(X1,Y1,0.), s2 = sphere(X2,Y2,0.);
      u.x[] = f[]*(s1 > s2 ? UD : -UD);
      u.y[] = 0.;
      u.z[] = 0.;
    }
  }
}

/**
## Adaptivity */

event adapt (i++) {
  adapt_wavelet ({f, u.x, u.y, u.z},
                 (double[]){FEMAX, UEMAX, UEMAX, UEMAX}, MAXLEVEL, MINLEVEL);
}

/**
## Diagnostics

Per-step scalar log: time, dt, cell count, kinetic energy, interface area,
liquid volume (mass-conservation check) and the interface's x/y bounding box
(how far the ligament stretches). */

event logfile (i += 10)
{
  double ke = 0., vol = 0.;
  double xmin = HUGE, xmax = -HUGE, ymin = HUGE, ymax = -HUGE;
  foreach (reduction(+:ke) reduction(+:vol)
           reduction(min:xmin) reduction(max:xmax)
           reduction(min:ymin) reduction(max:ymax)) {
    double dv = f[]*dv();
    vol += dv;
    foreach_dimension()
      ke += 0.5*rho(f[])*sq(u.x[])*dv();
    if (f[] > 1e-3 && f[] < 1. - 1e-3) {   // interface-containing cells
      if (x < xmin) xmin = x;  if (x > xmax) xmax = x;
      if (y < ymin) ymin = y;  if (y > ymax) ymax = y;
    }
  }
  double area = interface_area (f);

  static FILE * fp = NULL;
  if (!fp) {
    fp = fopen ("stats.dat", "w");
    fprintf (fp, "# t dt cells ke area volume xmin xmax ymin ymax\n");
  }
  fprintf (fp, "%g %g %ld %g %g %g %g %g %g %g\n",
           t, dt, grid->tn, ke, area, vol,
           xmin/DIAM, xmax/DIAM, ymin/DIAM, ymax/DIAM);
  fflush (fp);

  fprintf (stderr, "t=%.4e ms  dt=%.2e  cells=%ld  ke=%.4e  V/V0=%.6f\n",
           t*1e3, dt, grid->tn, ke, vol/V0);
}

/**
## Interface output

Interface polygons every frame (feeds the existing Blender / isosurface
pipeline), plus a full-field `dump` every 0.1 ms for restart and post. */

event interface (t += TSAVE)
{
  char name[80];
  sprintf (name, "facets-%04d.dat", (int)(t/TSAVE + 0.5));
  FILE * fp = fopen (name, "w");
  output_facets (f, fp);
  fclose (fp);
}

event snapshot (t += 1e-4)
{
  char name[80];
  sprintf (name, "dump-%04d", (int)(t/1e-4 + 0.5));
  dump (name);
}

/**
## Movies

A `z = 0` slice of the VOF field and the velocity magnitude. */

event movie (t += TSAVE)
{
  scalar umag[];
  foreach()
    umag[] = norm(u);
  output_ppm (f, file = "f.mp4", n = 800,
              box = {{-2.*DIAM,-2.*DIAM},{2.*DIAM,2.*DIAM}},
              min = 0, max = 1, linear = false);
  output_ppm (umag, file = "umag.mp4", n = 800,
              box = {{-2.*DIAM,-2.*DIAM},{2.*DIAM,2.*DIAM}},
              min = 0, max = UR, linear = true);
}

event end (t = TEND)
{
  dump ("dump-final");
}
