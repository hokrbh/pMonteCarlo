#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../vector.h"
#include "../hybridTaus.h"
#include "mcmlSingle.h"

PHOTON_DATA mcmlSingle(PROP cfg, TAUS_SEED *globalTausSeed, double *absData, int *error)
{
  PHOTON p;
  PHOTON_DATA dat;
  // Initialize photon seed
  p.seed = *globalTausSeed;
  dat.initialSeed = p.seed;
  // Initialize photon position
  p.r.x = 0.0;
  p.r.y = 0.0;
  p.r.z = 0.0;
  // Initialize photon directions
  p.v.x = 0.0;
  p.v.y = 0.0;
  p.v.z = 1.0;
  p.rid = 0;
  p.scat = false;
  p.det = 0;
  p.t = 0.0;
  p.w = 1.0;
  p.s = 0.0;
  
  // Propagate photon
  while( p.det == 0 )
  {
    if( (p.v.x*p.v.x + p.v.y*p.v.y + p.v.z*p.v.z - 1.0) > 1.0e-6 )
    {
      printf("Velocity not a unit vector: |(%g, %g, %g)| = %g\n", p.v.x, p.v.y, p.v.z, p.v.x*p.v.x + p.v.y*p.v.y + p.v.z*p.v.z);
    }
    if( p.s == 0.0 )
    {
      p.s = -log( hybridTaus( &(p.seed) ) );
    }
    double dr;
    if( cfg.layer[p.rid].us_permm == 0.0 ) // No scattering
    {
      dr = cfg.maxStep_mm;
    }
    else
    {
      dr = p.s/cfg.layer[p.rid].us_permm;
    }
    // Determine if the photon will reach a boundary before traveling a distance dr
    bool boundary = false;
    unsigned int newRegion;
    // Photon must be moving left or right to hit a boundary
    if( p.v.z != 0.0 ) 
    {
      double drBoundary;
      // Photon is moving left
      if( p.v.z < 0.0 )
      {
        drBoundary = (cfg.layer[p.rid].leftZ_mm-p.r.z)/p.v.z;
        newRegion = p.rid-1;
      }
      // Photon is moving right
      else
      {
        drBoundary = (cfg.layer[p.rid].rightZ_mm-p.r.z)/p.v.z;
        newRegion = p.rid+1;
        // Region to the right will exit the sample
        if( newRegion > cfg.numLayers )
        {
          newRegion = 0;
        }
      }
      if( dr >= drBoundary )
      {
        dr = drBoundary;
        boundary = true; // Photon will reach interface
      }
    }
    // Compute the deposited absoprtion while moving photon if enabled
    // We must add cfg.epsilon to distances for photons moving left to avoid
    // death by round off errors
    if( cfg.logAbsProfile != 0 )
    {
      double d_rem = dr;
      REAL3 temp = {0.0,0.0,-500.0};
      REAL3 pos = p.r;
      while( d_rem > 0.0 )
      {
        if( pos.x != temp.x || pos.y != temp.y || pos.z != temp.z )
        {
          temp = pos;
        }
        else
        {
          printf("Error\n");
          printf("%g\t%g\t%g\n", pos.x, pos.y, pos.z);
          printf("%g\t%g\t%g\n", temp.x, temp.y, temp.z);
          exit(0);
        }

        //printf("Position: %g\t%g\t%g\n", pos.x, pos.y, pos.z);
        //printf("Velocity: %g\t%g\t%g\n", p.v.x, p.v.y, p.v.z);
        
        // Make sure photon is in the ROI
        if( pos.x >= cfg.grid.x.min && pos.x < cfg.grid.x.max &&
            pos.y >= cfg.grid.y.min && pos.y < cfg.grid.y.max &&
            pos.z >= cfg.grid.z.min && pos.z < cfg.grid.z.max )
        {
          // Find current voxel
          unsigned int idx = (unsigned int)floor((pos.x-cfg.grid.x.min)/cfg.grid.x.d);
          unsigned int idy = (unsigned int)floor((pos.y-cfg.grid.y.min)/cfg.grid.y.d);
          unsigned int idz = (unsigned int)floor((pos.z-cfg.grid.z.min)/cfg.grid.z.d);

          //printf("In region: %u\t%u\t%u\n", idx, idy, idz);

          // Find distance to sides in each coordinate
          double dx, dy, dz;

          if( p.v.x > 0.0 ) // Use left + d
          {
            dx = (cfg.grid.x.lefts[idx]+cfg.grid.x.d - pos.x)/p.v.x + cfg.epsilon;
          }
          else if( p.v.x < 0.0 ) // Use left
          {
            dx = (cfg.grid.x.lefts[idx] - pos.x)/p.v.x + cfg.epsilon;
          }
          else // Not moving in the x-direction, set to HUGE_VAL to avoid NaN's
          {
            dx = HUGE_VAL;
          }

          if( p.v.y > 0.0 ) // Use left + d
          {
            dy = (cfg.grid.y.lefts[idy]+cfg.grid.y.d - pos.y)/p.v.y + cfg.epsilon;
          }
          else if( p.v.y < 0.0 ) // Use left
          {
            dy = (cfg.grid.y.lefts[idy] - pos.y)/p.v.y + cfg.epsilon;
          }
          else // Not moving in the y-direction, set to HUGE_VAL to avoid NaN's
          {
            dy = HUGE_VAL;
          }

          if( p.v.z > 0.0 ) // Use left + d
          {
            dz = (cfg.grid.z.lefts[idz]+cfg.grid.z.d - pos.z)/p.v.z + cfg.epsilon;
          }
          else if( p.v.z < 0.0 ) // Use left
          {
            dz = (cfg.grid.z.lefts[idz] - pos.z)/p.v.z + cfg.epsilon;
          }
          else // Not moving in the z-direction, set to HUGE_VAL to avoid NaN's
          {
            dz = HUGE_VAL;
          }

          //printf("%g\t%g\t%g\n", dx, dy, dz);

          // Whichever d is the smallest tells us which side will hit first
          double d = dx;
          if (d > dy) d = dy;
          if (d > dz) d = dz;
          
          //printf("%g\t%g\n", d, d_rem);

          if (d > d_rem) d = d_rem; // Makes sure we don't overshoot

          // Add the absorbed energy to the voxel
          double dw = p.w*(1.0-exp(-cfg.layer[p.rid].ua_permm*d)); // Amount of light absorbed
          absData[voxIndex(idx,idy,idz,cfg)] += dw;
          p.w -= dw; // Update the photons weight
          // Move the photon a distance d
          pos.x += p.v.x*d;
          pos.y += p.v.y*d;
          pos.z += p.v.z*d;
          d_rem -= d;
          
          //printf("%g\t%g\t%g\n\n", pos.x, pos.y, pos.z);
        }
        else // we need to check if the photon will enter the ROI before d_rem
        {
          double dx, dy, dz;
          if( pos.x < cfg.grid.x.min && p.v.x > 0.0 ) // Photon must be moving towards ROI
          {
            dx = (cfg.grid.x.min - pos.x)/p.v.x;
          }
          else if( pos.x > cfg.grid.x.max && p.v.x < 0.0 )
          {
            dx = (cfg.grid.x.max - pos.x)/p.v.x + cfg.epsilon; // so it trips the < 
          }
          else
          {
            dx = HUGE_VAL;
          }
          if( pos.y < cfg.grid.y.min && p.v.y > 0.0 ) // Photon must be moving towards ROI
          {
            dy = (cfg.grid.y.min - pos.y)/p.v.y;
          }
          else if( pos.y > cfg.grid.y.max && p.v.y < 0.0 )
          {
            dy = (cfg.grid.y.max - pos.y)/p.v.y + cfg.epsilon; // so it trips the < 
          }
          else
          {
            dy = HUGE_VAL;
          }
          if( pos.z < cfg.grid.z.min && p.v.z > 0.0 ) // Photon must be moving towards ROI
          {
            dz = (cfg.grid.z.min - pos.z)/p.v.z;
          }
          else if( pos.z > cfg.grid.z.max && p.v.z < 0.0 )
          {
            dz = (cfg.grid.z.max - pos.z)/p.v.z + cfg.epsilon; // so it trips the < 
          }
          else
          {
            dz = HUGE_VAL;
          }
          double d = dx;
          if( d > dy ) d = dy;
          if( d > dz ) d = dz;
          if( d < d_rem ) // Photon will make it to ROI, move to edge of ROI 
          {
            pos.x += p.v.x*d;
            pos.y += p.v.y*d;
            pos.z += p.v.z*d;
            p.w *= exp(-cfg.layer[p.rid].ua_permm*d);
            d_rem -= d;
          }
          else // Photon did not make it into ROI, move a full distance d_rem
          {
            pos.x += p.v.x*d_rem;
            pos.y += p.v.y*d_rem;
            pos.z += p.v.z*d_rem;
            p.w *= exp(-cfg.layer[p.rid].ua_permm*d_rem);
            d_rem = 0.0;
          }
        }
      }
    }
    else
    {
      p.w *= exp(-cfg.layer[p.rid].ua_permm*dr); // Decrease weight to account for absorption
    }
    // Move photon a distance dr
    p.r.x += p.v.x*dr;
    p.r.y += p.v.y*dr;
    p.r.z += p.v.z*dr;
    p.s -= cfg.layer[p.rid].us_permm*dr;
    p.t += cfg.layer[p.rid].n/C0*dr;
    // Make sure time hasn't expired
    if( p.t > cfg.maxTime_ps )
    {
      p.det = 5;
      storePhotonData( &dat, p );
    }
    // Compute what happens at the boundary
    if( boundary == true )
    {
      // Incident angle to the boundary
      double alpha = acos(fabs(p.v.z));
      double beta;
      double ref;
      // If moving from a high index to a low index total internal reflection is possible
      if( cfg.layer[p.rid].n > cfg.layer[newRegion].n && alpha >= asin(cfg.layer[newRegion].n/cfg.layer[p.rid].n) )
      {
          ref = 1.0;
          beta = 0.0; // This doesn't matter, but the compiler won't shut up about it
      }
      // Compute Fresnel reflections
      else
      {
          // Transmission angle (Snell's law)
          beta = asin(cfg.layer[p.rid].n/cfg.layer[newRegion].n*sin(alpha));
          double cosAlpha = cos(alpha);
          double cosBeta = cos(beta);
          double temp1 = (cfg.layer[p.rid].n*cosBeta - cfg.layer[newRegion].n*cosAlpha);
          double temp2 = (cfg.layer[p.rid].n*cosBeta + cfg.layer[newRegion].n*cosAlpha);
          double temp3 = (cfg.layer[p.rid].n*cosAlpha - cfg.layer[newRegion].n*cosBeta);
          double temp4 = (cfg.layer[p.rid].n*cosAlpha + cfg.layer[newRegion].n*cosBeta);
          ref = 0.5*( (temp1*temp1)/(temp2*temp2) + (temp3*temp3)/(temp4*temp4) );
      }
      if( hybridTaus( &(p.seed) ) > ref )
      {
        // Photon transmits through the boundary
        p.v.x = p.v.x*cfg.layer[p.rid].n/cfg.layer[newRegion].n;
        p.v.y = p.v.y*cfg.layer[p.rid].n/cfg.layer[newRegion].n;
        p.v.z = copysign(1.0,p.v.z)*cos(beta);
        /*if( (p.v.x*p.v.x + p.v.y*p.v.y + p.v.z*p.v.z - 1.0) > 1.0e-6 )
        {
          printf("%g\n", copysign(1.0,p.v.z));
          printf("%g\n", cos(alpha));
          printf("Old layer: %u, new layer: %u\n", p.rid, newRegion);
          printf("Old layer index: %g, new layer index: %g\n", cfg.layer[p.rid].n, cfg.layer[newRegion].n);
          printf("Refraction: |(%g, %g, %g)| = %g\n", p.v.x, p.v.y, p.v.z, p.v.x*p.v.x + p.v.y*p.v.y + p.v.z*p.v.z);
        }*/
        p.rid = newRegion;
        if( p.rid == 0 )
        {
          //Photon has exited the region
          if( p.v.z < 0.0 )
          {
            // Photon exits in reflection
            p.det = 1;
            storePhotonData( &dat, p );
          }
          if( p.v.z > 0.0 )
          {
            // Photon exits in transmission
            p.det = 2;
            storePhotonData( &dat, p );
          }
        }
      }
      else
      {
        // Photon reflects at the boundary
        p.v.z = -p.v.z;
        /*if( (p.v.x*p.v.x + p.v.y*p.v.y + p.v.z*p.v.z - 1.0) > 1.0e-6 )
        {
          printf("Reflection: |(%g, %g, %g)| = %g\n", p.v.x, p.v.y, p.v.z, p.v.x*p.v.x + p.v.y*p.v.y + p.v.z*p.v.z);
        }*/
        if( p.rid == 0 )
        {
          // Photon was reflected specularly
          p.det = 3;
          storePhotonData( &dat, p );
        }
      }
    }
    
    // Scatter photon if it has reached a scattering event
    if( p.s == 0.0 )
    {
      double phi = 2.0*PI*hybridTaus( &(p.seed) );
      double xi = hybridTaus( &(p.seed) );
      double cTheta = cosTheta( xi, cfg.layer[p.rid].g );
      double sTheta = sinTheta( xi, cfg.layer[p.rid].g );
      double cPhi = cos( phi );
      double sPhi = sin( phi );
      if( fabs(p.v.z) > 1.0-1.0E-10 )
      {
        p.v.x = sTheta*cPhi;
        p.v.y = sTheta*sPhi;
        p.v.z = copysign(1.0,p.v.z)*cTheta;
        /*if( (p.v.x*p.v.x + p.v.y*p.v.y + p.v.z*p.v.z - 1.0) > 1.0e-6 )
        {
          printf("Scatter along z: |(%g, %g, %g)| = %g\n", p.v.x, p.v.y, p.v.z, p.v.x*p.v.x + p.v.y*p.v.y + p.v.z*p.v.z);
        }*/
      }
      else
      {
        double temp = sqrt( 1.0 - p.v.z*p.v.z );
        REAL3 vTemp;
        vTemp.x = sTheta/temp*(p.v.x*p.v.z*cPhi - p.v.y*sPhi) + p.v.x*cTheta;
        vTemp.y = sTheta/temp*(p.v.y*p.v.z*cPhi + p.v.x*sPhi) + p.v.y*cTheta;
        vTemp.z = -sTheta*cPhi*temp + p.v.z*cTheta;
        /*if( fabs(vTemp.x*vTemp.x + vTemp.y*vTemp.y + vTemp.z*vTemp.z - 1.0) > 1.0e-6 )
        {
          printf("xi = %g\n", xi);
          printf("p.rid = %d\n", p.rid);
          printf("g = %g\n", cfg.layer[p.rid].g);
          printf("sin(theta) = %g, cos(theta) = %g, mag = %g\n", sTheta, cTheta, sTheta*sTheta + cTheta*cTheta);
          printf("sin(phi) = %g, cos(phi) = %g, mag = %g\n", sPhi, cPhi, sPhi*sPhi + cPhi*cPhi);
          printf("Initial Vel: |(%g, %g, %g)| = %g\n", p.v.x, p.v.y, p.v.z, p.v.x*p.v.x + p.v.y*p.v.y + p.v.z*p.v.z);
          printf("Scatter: |(%g, %g, %g)| = %g\n", vTemp.x, vTemp.y, vTemp.z, vTemp.x*vTemp.x + vTemp.y*vTemp.y + vTemp.z*vTemp.z);
        }*/
        p.v = vTemp;
      }
      p.scat = true;
    }
    
    if( p.w < cfg.weightThreshold )
    {
      if( cfg.useRoulette == 1 )
      {
        // Play russian roulette
        if( hybridTaus( &(p.seed) ) < cfg.rouletteProb )
        {
          p.w = p.w/cfg.rouletteProb;
        }
        else
        {
          p.w = 0.0;
          p.det = 4;
          storePhotonData( &dat, p );
        }
      }
      else
      {
        // Quit tracking photons below a certain weight
        p.w = 0.0;
        p.det = 4;
        storePhotonData( &dat, p );
      }
    }

    // Store photon data to array
  }
  return(dat);
}

unsigned int voxIndex( unsigned int idx, unsigned int idy, unsigned int idz, PROP cfg )
{
  return( cfg.grid.x.n*cfg.grid.y.n*idz + cfg.grid.x.n*idy + idx );
}

// Computes the cosine of the scattering angle for a given g using the Heyney-Greenstien distribution for a given uniform random number xi
double cosTheta( double xi, double g )
{
  double ans;
  if( g == 0.0 )
  {
    ans = 2*xi - 1;
  }
  else
  {
    double temp = (1.0-g*g)/(1.0-g+2.0*g*xi);
    ans = 1.0/(2.0*g)*( 1.0 + g*g - temp*temp );
  }
  return(ans);
}

// Computes the sine of the scattering angle for a given g using the Heyney-Greenstien distribution for a given uniform random number xi
double sinTheta( double xi, double g )
{
  double cTheta = cosTheta( xi, g );
  return( sqrt( 1.0 - cTheta*cTheta ) );
}

void storePhotonData( PHOTON_DATA *dat, PHOTON p )
{
  dat->r.x = p.r.x;
  dat->r.y = p.r.y;
  dat->v = p.v;
  dat->w = p.w;
  dat->t = p.t;
  dat->det = p.det;
}
