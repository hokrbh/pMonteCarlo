#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "vector.h"
#include "hybridTaus.h"
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
        int newRegion;
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
        // Move photon a distance dr
        p.r.x += p.v.x*dr;
        p.r.y += p.v.y*dr;
        p.r.z += p.v.z*dr;
        // Decrease weight to account for absorption
        p.w *= exp(-cfg.layer[p.rid].ua_permm*dr);
        p.s -= cfg.layer[p.rid].us_permm*dr;
        p.t += cfg.layer[p.rid].n/C0*dr;
        // Compute what happens at the boundary
        if( boundary == true )
        {
            // Incident angle to the boundary
            double alpha = acos(fabs(p.v.z));
            double ref;
            // If moving from a high index to a low index total internal reflection is possible
            if( cfg.layer[p.rid].n > cfg.layer[newRegion].n && alpha >= asin(cfg.layer[newRegion].n/cfg.layer[p.rid].n) )
            {
                ref = 1.0;
            }
            // Compute Fresnel reflections
            else
            {
                // Transmission angle (Snell's law)
                double beta = asin(cfg.layer[p.rid].n/cfg.layer[newRegion].n*sin(alpha));
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
                p.v.z = copysign(1.0,p.v.z)*cos(alpha);
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
            }
            else
            {
                double temp = sqrt( 1.0 - p.v.z*p.v.z );
                REAL3 vTemp;
                vTemp.x = sTheta/temp*(p.v.y*sPhi-p.v.z*p.v.x*cPhi) + p.v.x*cTheta;
                vTemp.y = sTheta/temp*(p.v.x*sPhi-p.v.z*p.v.y*cPhi) + p.v.y*cTheta;
                vTemp.z = sTheta*temp*cPhi + p.v.z*cTheta;
                p.v = vTemp;
            }
            p.scat = true;
        }
        
        if( p.w < cfg.weightThreshold )
        {
            // Quit tracking photons below a certain weight
            p.w = 0.0;
            p.det = 4;
            // Play russian roulette
            /*if( hybridTaus( &(p.seed) ) < cfg.rouletteProb )
            {
                p.w = p.w/cfg.rouletteProb;
            }
            else
            {
                p.w = 0.0;
                p.det = 4;
            }*/
            storePhotonData( &dat, p );
        }
    }
    return(dat);
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
		ans = 1.0/(2*g)*( 1.0 + g*g - temp*temp );
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
