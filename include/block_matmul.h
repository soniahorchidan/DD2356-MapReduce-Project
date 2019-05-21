#ifndef __CANNON_H__
#define __CANNON_H__

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include <cblas.h>

void init_matmul(char*, char*, char*);
void cleanup_matmul();
void compute_fox();
void compute_cannon();

#endif
