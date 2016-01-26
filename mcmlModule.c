#include <Python.h>
#include <stdbool.h>
#include <float.h>
#include "vector.h"
#include "hybridTaus.h"
#include "allocate.h"
#include "mcmlSingle.h"
#include "mcmlModule.h"

static PyObject *_mcml(PyObject* self, PyObject* args)
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
    cfg.weightThreshold = 1.0E-25;
    cfg.rouletteProb = 0.1;
    
    double reflectionCoeff = 0.0;
    double transmissionCoeff = 0.0;
    double specularCoeff = 0.0;
    // Run Monte Carlo simulations
    for( int i = 0; i < cfg.numPhotons; i++ )
    {
        TAUS_SEED photonSeed;
        photonSeed.z1 = hybridTausInt(&(cfg.globalTausSeed));
        photonSeed.z2 = hybridTausInt(&(cfg.globalTausSeed));
        photonSeed.z3 = hybridTausInt(&(cfg.globalTausSeed));
        photonSeed.z4 = hybridTausInt(&(cfg.globalTausSeed));
        data[i] = mcmlSingle(cfg, &photonSeed, &error);
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
    
    //printf("R = %lg\tT = %lg\n", reflectionCoeff, transmissionCoeff);
    
    // Build the return values
    PyObject *ret = Py_BuildValue("ddd", reflectionCoeff, transmissionCoeff, specularCoeff);
    
    free(cfg.layer);
    free(data);

    //Py_RETURN_NONE;
    return ret;
}

static PyMethodDef mcmlMethods[] =
{
     {"mcml", _mcml, METH_VARARGS, "Monte Carlo Multi-Layer (Wang1995)."},
     {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initpMonteCarlo(void)
{
     (void) Py_InitModule3("pMonteCarlo", mcmlMethods, "A Monte Carlo library for python.");
}

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
    if( sscanf(cfgString, "%u,%u,%lg,[%[^]]],[%[^]]],[%[^]]],[%[^]]],[%[^]]],[%[^]]]", &cfg->numPhotons, &cfg->globalSeed, &cfg->backgroundIndex, leftZString, rightZString, indexString, anisotropyString, usString, uaString) < 9 )
    {
        printf("Error reading input string\n");
        return(-1);
    }
    printf("Running %d photons\n", cfg->numPhotons);
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
    printf("Number of layers: %u\n", cfg->numLayers);
    cfg->layer = (LAYER_PROP*)allocate( (cfg->numLayers+1)*sizeof(LAYER_PROP) );
    tempString = strtok(leftZString,", ");
    for( int i = 0; i < cfg->numLayers; i++ )
    {
        cfg->layer[i+1].leftZ_mm = atof(tempString);
        tempString = strtok(NULL, ", ");
    }
    tempString = strtok(rightZString,", ");
    for( int i = 0; i < cfg->numLayers; i++ )
    {
        cfg->layer[i+1].rightZ_mm = atof(tempString);
        tempString = strtok(NULL, ", ");
    }
    tempString = strtok(indexString,", ");
    for( int i = 0; i < cfg->numLayers; i++ )
    {
        cfg->layer[i+1].n = atof(tempString);
        tempString = strtok(NULL, ", ");
    }
    tempString = strtok(anisotropyString,", ");
    for( int i = 0; i < cfg->numLayers; i++ )
    {
        cfg->layer[i+1].g = atof(tempString);
        tempString = strtok(NULL, ", ");
    }
    tempString = strtok(usString,", ");
    for( int i = 0; i < cfg->numLayers; i++ )
    {
        cfg->layer[i+1].us_permm = atof(tempString);
        tempString = strtok(NULL, ", ");
    }
    tempString = strtok(uaString,", ");
    for( int i = 0; i < cfg->numLayers; i++ )
    {
        cfg->layer[i+1].ua_permm = atof(tempString);
        tempString = strtok(NULL, ", ");
    }
    return(0);
}
