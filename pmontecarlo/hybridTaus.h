#ifndef HYBRID_TAUS_H
#define HYBRID_TAUS_H

typedef struct
{
	unsigned int z1, z2, z3, z4;
} TAUS_SEED;

// Function definitions
unsigned int tausStep( unsigned int *z, int S1, int S2, int S3, unsigned int M );
unsigned int LCGStep( unsigned int *z, unsigned int a, unsigned int c );
double hybridTaus( TAUS_SEED *seed );
unsigned int hybridTausInt( TAUS_SEED *seed );

#endif
