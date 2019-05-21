#ifndef __FRAMWORK_H__
#define __FRAMWORK_H__

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>

void read_file(char *);

void flat_map();

void reduce();

void write_file(char *);
#endif
