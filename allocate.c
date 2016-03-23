// Include standard C libraries
#include <stdlib.h>
#include <stdio.h>

// Include header file
#include "allocate.h"

void *allocate( unsigned int size )
{
	void *ptr = malloc(size);
	if(ptr == NULL)
	{
		fprintf(stderr,"Memory allocation error, terminating\n");
		exit(EXIT_FAILURE);
	}
	return(ptr);
}

void *callocate( unsigned int size )
{
	void *ptr = calloc(size,1);
	if(ptr == NULL)
	{
		fprintf(stderr,"Memory allocation error, terminating\n");
		exit(EXIT_FAILURE);
	}
	return(ptr);
}
