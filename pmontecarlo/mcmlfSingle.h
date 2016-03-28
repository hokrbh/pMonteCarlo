#ifndef MCMLfSINGLE_H
#define MCMLfSINGLE_H

#ifndef STR_SIZE
#define STR_SIZE            1024
#endif

#ifndef PI
#define PI					3.141592653589793
#endif

#ifndef C0
#define C0					0.299792458 // Speed of light in vacuum (mm/ps)
#endif

typedef struct
{
  /* Photon Data */
  REAL3 r;				// Position vector of photon
  REAL3 v;				// Unit vector giving the direction of the photon
  double w;				// Photon weight
  double t;				// Time
  double s;				// Demensionless step size (see Wang1995)
  TAUS_SEED seed;			// Seed for the random number generator
  /* Flags */
  unsigned int rid;		// Region photon is in
  bool scat;				// Flag for whether or not the photon has scattered
  int det;				// Flag indicating if the photon has been detected ( 0-no; 1-backwards; 2-forwards; 3-specular; 4-killed by roulette )	
} PHOTON;

typedef struct
{
  REAL2 r;
  REAL3 v;
  double w;
  double t;
  int det;
  TAUS_SEED initialSeed;
  // det is 0-not detected, 1-reflection, 2-transmission, 3-specular, 4-killed by absorption threshold
} PHOTON_DATA;

typedef struct
{
  double leftZ_mm;                // Left z-value of the layer
  double rightZ_mm;               // Right z-value of the layer
  double n;                       // Index of refraction of the layer
  double g;                       // Anisotropy (<cos(theta)>) of the layer
  double us_permm;                // Scattering coefficient in mm^(-1)
  double ua_permm;                // Absorption coefficient in mm^(-1)
} LAYER_PROP;

typedef struct
{
  double min;			    // Minimum of the grid
  double max;			    // Maximum of the grid
  unsigned int n;		    // Number of grid boxes
  double d;                       // Length of a grid box
  double *lefts;                  // Value of the left side of a grid box
  double *centers;                // Value of the center of a grid box
} GRID_PARAM;

typedef struct
{
  GRID_PARAM x;		    // Parameters for the x-direction
  GRID_PARAM y;		    // Parameters for the y-direction
  GRID_PARAM z;		    // Parameters for the z-direction
} VOXEL_PARAM;

typedef struct
{
  unsigned int numPhotons;	    // Number of photons to simulate
  unsigned int globalSeed;	    // Global seed for the simulation
  double backgroundIndex;         // Index of refraction of the background medium
  LAYER_PROP *layer;              // Struct to hold the layer properties
  unsigned int numLayers;         // Number of layers used
  TAUS_SEED globalTausSeed;       // Global 4 element Taus seed
  unsigned int writeDetData;	    // Flag to determine if detection data is to be written
  char detDataFilename[STR_SIZE]; // Filename for the detection data
  double maxStep_mm;              // Max distance photon is allowed to travel
  double weightThreshold;         // Threshold weight
  double rouletteProb;            // Roulette probability
  unsigned int logAbsProfile;	    // Flag to track absorption data and write file
  char absDataFilename[STR_SIZE]; // Filename for the absorption profile data
  VOXEL_PARAM grid;		    // Grid properties for tracking deposited energy
  double epsilon;                 // Fudge factor to trip less than for floats
} PROP;

PHOTON_DATA mcmlfSingle(PROP cfg, TAUS_SEED *globalTausSeed, double *absData, int *error);
unsigned int voxIndex( unsigned int idx, unsigned int idy, unsigned int idz, PROP cfg );
double cosTheta( double xi, double g );
double sinTheta( double xi, double g );
void storePhotonData( PHOTON_DATA *dat, PHOTON p );

#endif
