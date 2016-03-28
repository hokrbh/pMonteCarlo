#include <Python.h>
#include <stdbool.h>
#include <float.h>
#include "../vector.h"
#include "../hybridTaus.h"
#include "../allocate.h"
#include "mcmlSingle.h"
#include "mcmlModule.h"

/* Anything called run is the function name while mcml_run is the module name */

static PyObject *_run(PyObject* self, PyObject* args)
{
  int error = 0;
  // Read in paramters as a string from Python
  const char* cfgString;
  if (!PyArg_ParseTuple(args, "s", &cfgString))
    return NULL;
  
  PROP cfg;
  if(parseInputString(cfgString, &cfg) < 0)
  {
    fprintf(stderr,"Error parsing inputs\n");
    error = -2; // Currently doesn't do anything
    return NULL;
  }
  // Set background layer properties
  cfg.layer[0].leftZ_mm = -1.0; // Unused
  cfg.layer[0].rightZ_mm = 0.0;
  cfg.layer[0].n = cfg.backgroundIndex;
  cfg.layer[0].g = 0.0;
  cfg.layer[0].us_permm = 0.0;
  cfg.layer[0].ua_permm = 0.0;
  
  // Initialize MC simulation
  PHOTON_DATA *data = (PHOTON_DATA*)allocate(cfg.numPhotons*sizeof(PHOTON_DATA));
  srand(cfg.globalSeed);
  cfg.globalTausSeed.z1 = rand();
  cfg.globalTausSeed.z2 = rand();
  cfg.globalTausSeed.z3 = rand();
  cfg.globalTausSeed.z4 = rand();
  cfg.maxStep_mm = 1.0E6;
  cfg.weightThreshold = 1.0E-9;
  cfg.rouletteProb = 0.1;

  // Initialize absorption grid
  cfg.epsilon = 1.0E-9;
  double *absData; // Declare so there isn't an error in passing it, but don't allocate
  if(cfg.logAbsProfile != 0)
  {
    cfg.grid.x.d = (cfg.grid.x.max-cfg.grid.x.min)/(double)cfg.grid.x.n;
    cfg.grid.x.lefts = (double*)allocate(cfg.grid.x.n*sizeof(double));
    cfg.grid.x.centers = (double*)allocate(cfg.grid.x.n*sizeof(double));
    for( unsigned int i=0; i<cfg.grid.x.n; i++ )
    {
      cfg.grid.x.lefts[i] = cfg.grid.x.min + (double)i*cfg.grid.x.d;
      cfg.grid.x.centers[i] = cfg.grid.x.lefts[i]+0.5*cfg.grid.x.d;
    }
    cfg.grid.y.d = (cfg.grid.y.max-cfg.grid.y.min)/(double)cfg.grid.y.n;
    cfg.grid.y.lefts = (double*)allocate(cfg.grid.y.n*sizeof(double));
    cfg.grid.y.centers = (double*)allocate(cfg.grid.y.n*sizeof(double));
    for( unsigned int i=0; i<cfg.grid.y.n; i++ )
    {
      cfg.grid.y.lefts[i] = cfg.grid.y.min + (double)i*cfg.grid.y.d;
      cfg.grid.y.centers[i] = cfg.grid.y.lefts[i]+0.5*cfg.grid.y.d;
    }
    cfg.grid.z.d = (cfg.grid.z.max-cfg.grid.z.min)/(double)cfg.grid.z.n;
    cfg.grid.z.lefts = (double*)allocate(cfg.grid.z.n*sizeof(double));
    cfg.grid.z.centers = (double*)allocate(cfg.grid.z.n*sizeof(double));
    for( unsigned int i=0; i<cfg.grid.z.n; i++ )
    {
      cfg.grid.z.lefts[i] = cfg.grid.z.min + (double)i*cfg.grid.z.d;
      cfg.grid.z.centers[i] = cfg.grid.z.lefts[i]+0.5*cfg.grid.z.d;
    }
    absData = (double*)callocate(cfg.grid.x.n*cfg.grid.y.n*cfg.grid.z.n*sizeof(double));
  }
  else
  {
    // Go ahead and initialize absData to avoid annoying compiler warning,
    // but make it small so we don't waste RAM needlessly
    absData = (double*)callocate(1*sizeof(double));
  }

  double reflectionCoeff = 0.0;
  double transmissionCoeff = 0.0;
  double specularCoeff = 0.0;
  // Run Monte Carlo simulations
  for( unsigned int i = 0; i < cfg.numPhotons; i++ )
  {
    TAUS_SEED photonSeed;
    photonSeed.z1 = hybridTausInt(&(cfg.globalTausSeed));
    photonSeed.z2 = hybridTausInt(&(cfg.globalTausSeed));
    photonSeed.z3 = hybridTausInt(&(cfg.globalTausSeed));
    photonSeed.z4 = hybridTausInt(&(cfg.globalTausSeed));
    data[i] = mcmlSingle(cfg, &photonSeed, absData, &error);
    if(error != 0)
    {
        fprintf(stderr, "Error propagating photon %d\n", i);
        return NULL;
    }
    // Calculate reflection and transmission coefficients
    if( data[i].det == 1 )
    {
        // Reflection
        reflectionCoeff += data[i].w;
    }
    else if( data[i].det == 2 )
    {
        // Transmission
        transmissionCoeff += data[i].w;
    }
    else if( data[i].det == 3 )
    {
        // Specular Reflection
        specularCoeff += data[i].w;
    }
  }
  reflectionCoeff /= (double)cfg.numPhotons;
  transmissionCoeff /= (double)cfg.numPhotons;
  specularCoeff /= (double)cfg.numPhotons;
  double totalRefCoeff = reflectionCoeff+specularCoeff;
  
  // If flag is true write out detection data
  if( cfg.writeDetData != 0 )
  {
    FILE *detData = fopen(cfg.detDataFilename, "w");
    if(detData == NULL)
    {
      fprintf(stderr, "Error opening detection data file %s\n", cfg.detDataFilename);
      return NULL;
    }
    if( cfg.writeDetData == 1 ) // Write out binary data
    {
      fwrite(data, sizeof(PHOTON_DATA), cfg.numPhotons, detData);
    }
    else if( cfg.writeDetData == 2 ) // Write out ASCII data
    {
      for( unsigned int i = 0; i<cfg.numPhotons; i++ )
      {
          fprintf(detData, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%u\t%u\t%u\t%u\n", data[i].r.x, data[i].r.y, data[i].v.x, data[i].v.y, data[i].v.z, data[i].w, data[i].t, data[i].det, data[i].initialSeed.z1, data[i].initialSeed.z2, data[i].initialSeed.z3, data[i].initialSeed.z4);
      }
    }
    else
    {
      fprintf(stderr, "Error: Unrecognized writeDetData flag: %u\n", cfg.writeDetData);
      error = -3; // Currently this does nothing
      return(NULL);
    }
    fclose(detData); 
  }
  else
  {
    printf("Not writing detection data\n");
  }

  // If flag is true write out absorption data
  if( cfg.logAbsProfile != 0 )
  {
    FILE *absFile = fopen(cfg.absDataFilename, "w");
    if(absFile == NULL)
    {
      fprintf(stderr, "Error opening absorption data file %s\n", cfg.absDataFilename);
      return NULL;
    }
    if( cfg.logAbsProfile == 1 ) // Write out binary data
    {
      fwrite(absData, sizeof(double), cfg.grid.x.n*cfg.grid.y.n*cfg.grid.z.n, absFile);
    }
    else if( cfg.logAbsProfile == 2 ) // Write out ascii data
    {
      for( unsigned int k=0; k<cfg.grid.z.n; k++ )
      {
        for( unsigned int j=0; j<cfg.grid.y.n; j++ )
        {
          for( unsigned int i=0; i<cfg.grid.x.n; i++ )
          {
            fprintf(absFile, "%g\t", absData[ voxIndex( i, j, k, cfg ) ] );
          }
          fprintf(absFile, "\n");
        }
        fprintf(absFile, "\n");
      }
    }
    else
    {
      fprintf(stderr, "Error: Unrecognized writeDetData flag: %u\n", cfg.writeDetData);
      error = -3; // Currently this does nothing
      return(NULL);
    }
    fclose(absFile);
  }
  else
  {
    printf("Not writing absorption data\n");
  }
  
  //printf("%lu\t%lu\t%lu\n", sizeof(PHOTON_DATA), sizeof(TAUS_SEED), sizeof(int));

  //printf("R = %lg\tT = %lg\n", reflectionCoeff, transmissionCoeff);
  
  // Build the return values
  PyObject *ret = Py_BuildValue("dddd", totalRefCoeff, transmissionCoeff, reflectionCoeff, specularCoeff);
  
  free(cfg.layer);
  free(data);

  //Py_RETURN_NONE;
  return ret;
}

static PyMethodDef runMethods[] =
{
  {"run", _run, METH_VARARGS, "Monte Carlo Multi-Layer (Wang1995)."},
  {NULL, NULL, 0, NULL}
};

#ifdef IS_PY3K
static struct PyModuleDef run = {
  PyModuleDef_HEAD_INIT,
  "mcml_run",   /* name of module */
  NULL, /* module documentation, may be NULL */
  -1,       /* size of per-interpreter state of the module,
              or -1 if the module keeps state in global variables. */
  runMethods
};

PyMODINIT_FUNC PyInit_mcml_run(void)
{
  return PyModule_Create(&run);
}
#else
PyMODINIT_FUNC initmcml_run(void)
{
    (void) Py_InitModule3("run", runMethods, "A Monte Carlo library for python.");
}
#endif


// Parses the input string, cfgString, and stores it in the properties struct, cfg.
int parseInputString(const char *cfgString, PROP *cfg)
{
  // Parse string for values
  char leftZString[STR_SIZE];
  char rightZString[STR_SIZE];
  char indexString[STR_SIZE];
  char anisotropyString[STR_SIZE];
  char usString[STR_SIZE];
  char uaString[STR_SIZE];

  //printf("%s\n", cfgString);
  
  if( sscanf(cfgString, "%u,%u,%lg,[%[^]]],[%[^]]],[%[^]]],[%[^]]],[%[^]]],[%[^]]],%u,%[^,],%u,%[^,],%lg,%lg,%u,%lg,%lg,%u,%lg,%lg,%u", &cfg->numPhotons, &cfg->globalSeed, &cfg->backgroundIndex, leftZString, rightZString, indexString, anisotropyString, usString, uaString, &cfg->writeDetData, cfg->detDataFilename, &cfg->logAbsProfile, cfg->absDataFilename, &cfg->grid.x.min, &cfg->grid.x.max, &cfg->grid.x.n, &cfg->grid.y.min, &cfg->grid.y.max, &cfg->grid.y.n, &cfg->grid.z.min, &cfg->grid.z.max, &cfg->grid.z.n ) < 22 )
  {
    printf("Error reading input string\n");
    return(-1);
  }
  printf("Running %d photons\n", cfg->numPhotons);
  // printf("%s\n", cfg->absDataFilename);
  // Find the number of layers used
  unsigned int i = 0;
  char* tempString;
  char tempString2[STR_SIZE];
  strcpy(tempString2, leftZString);
  tempString = strtok(tempString2,", ");
  while( tempString != NULL )
  {
    i++;
    tempString = strtok(NULL, ", ");
    if( i > 100 )
    {
        printf("Error reading in layer parameters\n");
        break;
    }
  }
  cfg->numLayers = i;
  //printf("Number of layers: %u\n", cfg->numLayers);
  cfg->layer = (LAYER_PROP*)allocate( (cfg->numLayers+1)*sizeof(LAYER_PROP) );
  tempString = strtok(leftZString,", ");
  for( unsigned int i = 0; i < cfg->numLayers; i++ )
  {
    cfg->layer[i+1].leftZ_mm = atof(tempString);
    tempString = strtok(NULL, ", ");
  }
  tempString = strtok(rightZString,", ");
  for( unsigned int i = 0; i < cfg->numLayers; i++ )
  {
    cfg->layer[i+1].rightZ_mm = atof(tempString);
    tempString = strtok(NULL, ", ");
  }
  tempString = strtok(indexString,", ");
  for( unsigned int i = 0; i < cfg->numLayers; i++ )
  {
    cfg->layer[i+1].n = atof(tempString);
    tempString = strtok(NULL, ", ");
  }
  tempString = strtok(anisotropyString,", ");
  for( unsigned int i = 0; i < cfg->numLayers; i++ )
  {
    cfg->layer[i+1].g = atof(tempString);
    tempString = strtok(NULL, ", ");
  }
  tempString = strtok(usString,", ");
  for( unsigned int i = 0; i < cfg->numLayers; i++ )
  {
    cfg->layer[i+1].us_permm = atof(tempString);
    tempString = strtok(NULL, ", ");
  }
  tempString = strtok(uaString,", ");
  for( unsigned int i = 0; i < cfg->numLayers; i++ )
  {
    cfg->layer[i+1].ua_permm = atof(tempString);
    tempString = strtok(NULL, ", ");
  }
  
  return(0);
}
