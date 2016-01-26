#include <float.h>
#include <math.h>

#include "hybridTaus.h"

unsigned int tausStep( unsigned int *z, int S1, int S2, int S3, unsigned int M )
{
	return *z = ( ( ( *z & M ) << S3 ) ^ ( ( ( *z << S1 ) ^ *z ) >> S2 ) );
}

unsigned int LCGStep( unsigned int *z, unsigned int a, unsigned int c )
{
	return *z = ( a * *z + c );
}

double hybridTaus( TAUS_SEED *seed )
{
	return 2.3283064365387e-10*( tausStep( &seed->z1, 13, 19, 12, 4294967294UL )
		^ tausStep( &seed->z2, 2, 25, 4, 4294967288UL )
		^ tausStep( &seed->z3, 3, 11, 17, 4294967280UL )
		^ LCGStep( &seed->z4, 1664525, 1013904223UL ) );
}

unsigned int hybridTausInt( TAUS_SEED *seed )
{
	return tausStep( &seed->z1, 13, 19, 12, 4294967294UL )
		^ tausStep( &seed->z2, 2, 25, 4, 4294967288UL )
		^ tausStep( &seed->z3, 3, 11, 17, 4294967280UL )
		^ LCGStep( &seed->z4, 1664525, 1013904223UL );
}
