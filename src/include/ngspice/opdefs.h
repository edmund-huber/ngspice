/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

#ifndef DCOP
#define DCOP


#include "ngspice/cktdefs.h"
#include "ngspice/jobdefs.h"
#include "ngspice/tskdefs.h"
    /*
     * structures used to describe D.C. operationg point analyses to
     * be performed.
     */

typedef struct {
    int JOBtype;
    JOB *JOBnextJob;
    char *JOBname;
} OP;

#endif /*DCOP*/
