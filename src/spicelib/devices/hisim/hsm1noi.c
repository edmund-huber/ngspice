/***********************************************************************
 HiSIM v1.1.0
 File: hsm1noi.c of HiSIM v1.1.0

 Copyright (C) 2002 STARC

 June 30, 2002: developed by Hiroshima University and STARC
 June 30, 2002: posted by Keiichi MORIKAWA, STARC Physical Design Group
***********************************************************************/

/*
 * Modified by Paolo Nenzi 2002
 * ngspice integration
 */

#include "ngspice.h"
#include "hsm1def.h"
#include "cktdefs.h"
#include "iferrmsg.h"
#include "noisedef.h"
#include "suffix.h"
#include "const.h"  /* jwan */

/*
 * HSM1noise (mode, operation, firstModel, ckt, data, OnDens)
 *    This routine names and evaluates all of the noise sources
 *    associated with MOSFET's.  It starts with the model *firstModel and
 *    traverses all of its insts.  It then proceeds to any other models
 *    on the linked list.  The total output noise density generated by
 *    all of the MOSFET's is summed with the variable "OnDens".
 */

/*
 Channel thermal and flicker noises are calculated based on the value
 of model->HSM1_noise.
 If model->HSM1_noise = 1,
    Channel thermal noise = SPICE2 model
    Flicker noise         = SPICE2 model
 If model->HSM1_noise = 2,
    Channel thermal noise = HiSIM1 model corresponding to BSIM3 model
    Flicker noise         = HiSIM1 model
 If model->HSM1_noise = 3,
    Channel thermal noise = SPICE2 model 
    Flicker noise         = HiSIM1 model
 If model->HSM1_noise = 4,
    Channel thermal noise = HiSIM1 model corresponding to BSIM3 model
    Flicker noise         = SPICE2 model
 If model->HSM1_noise = 5,
    Channel thermal noise = NONE
    Flicker noise         = HiSIM1 model
 */

extern void   NevalSrc();
extern double Nintegrate();

int HSM1noise (int mode, int operation, GENmodel *inModel, CKTcircuit *ckt, 
               Ndata *data, double *OnDens)
{
  HSM1model *model = (HSM1model *)inModel;
  HSM1instance *here;
  char name[N_MXVLNTH];
  double tempOnoise;
  double tempInoise;
  double noizDens[HSM1NSRCS];
  double lnNdens[HSM1NSRCS];
  register int error, i;

  /* define the names of the noise sources */
  static char * HSM1nNames[HSM1NSRCS] = {
    /* Note that we have to keep the order
       consistent with the index definitions 
       in hsm1defs.h */
    ".rd",              /* noise due to rd */
    ".rs",              /* noise due to rs */
    ".id",              /* noise due to id */
    ".1ovf",            /* flicker (1/f) noise */
    ""                  /* total transistor noise */
  };
  
  for ( ;model != NULL; model = model->HSM1nextModel ) {
    for ( here = model->HSM1instances; here != NULL;
	  here = here->HSM1nextInstance ) {
	  
      if (here->HSM1owner != ARCHme)
	      continue;	  
	  
      switch (operation) {
      case N_OPEN:
	/* see if we have to to produce a summary report */
	/* if so, name all the noise generators */
	  
	if (((NOISEAN*)ckt->CKTcurJob)->NStpsSm != 0) {
	  switch (mode) {
	  case N_DENS:
	    for ( i = 0; i < HSM1NSRCS; i++ ) { 
	      (void) sprintf(name, "onoise.%s%s", 
			     (char *)here->HSM1name, HSM1nNames[i]);
	      data->namelist = 
		(IFuid *) trealloc((char *) data->namelist,
				   (data->numPlots + 1) * sizeof(IFuid));
	      if (!data->namelist)
		return(E_NOMEM);
	      (*(SPfrontEnd->IFnewUid)) 
		(ckt, &(data->namelist[data->numPlots++]),
		 (IFuid) NULL, name, UID_OTHER, (void **) NULL);
	    }
	    break;
	  case INT_NOIZ:
	    for ( i = 0; i < HSM1NSRCS; i++ ) {
	      (void) sprintf(name, "onoise_total.%s%s", 
			     (char *)here->HSM1name, HSM1nNames[i]);
	      data->namelist = 
		(IFuid *) trealloc((char *) data->namelist,
				   (data->numPlots + 1) * sizeof(IFuid));
	      if (!data->namelist)
		return(E_NOMEM);
	      (*(SPfrontEnd->IFnewUid)) 
		(ckt, &(data->namelist[data->numPlots++]),
		 (IFuid) NULL, name, UID_OTHER, (void **) NULL);
	      
	      (void) sprintf(name, "inoise_total.%s%s", 
			     (char *)here->HSM1name, HSM1nNames[i]);
	      data->namelist = 
		(IFuid *) trealloc((char *) data->namelist,
				   (data->numPlots + 1) * sizeof(IFuid));
	      if (!data->namelist)
		return(E_NOMEM);
	      (*(SPfrontEnd->IFnewUid)) 
		(ckt, &(data->namelist[data->numPlots++]),
		 (IFuid) NULL, name, UID_OTHER, (void **)NULL);
	    }
	    break;
	  }
	}
	break;
      case N_CALC:
	switch (mode) {
	case N_DENS:
	  NevalSrc(&noizDens[HSM1RDNOIZ], &lnNdens[HSM1RDNOIZ], 
		   ckt, THERMNOISE,
		   here->HSM1dNodePrime, here->HSM1dNode,
		   here->HSM1drainConductance * here->HSM1_m);
	  
	  NevalSrc(&noizDens[HSM1RSNOIZ], &lnNdens[HSM1RSNOIZ], 
		   ckt, THERMNOISE,
		   here->HSM1sNodePrime, here->HSM1sNode,
		   here->HSM1sourceConductance * here->HSM1_m);

	  switch( model->HSM1_noise ) {
	    double I;
	  case 1:
	  case 3:
	    I = here->HSM1_gm + here->HSM1_gds + here->HSM1_gmbs;
	    I *= (I < 0.0) ? -1.0 : 1.0;
	    I *= 2.0/3.0;
	    I *=  here->HSM1_m; /* PN */
	    NevalSrc(&noizDens[HSM1IDNOIZ], &lnNdens[HSM1IDNOIZ], 
		     ckt, THERMNOISE, 
		     here->HSM1dNodePrime, here->HSM1sNodePrime, I);
	    break;
	  case 2:
	  case 4:
	    I = -1.0 * (here->HSM1_qg + here->HSM1_qb)
	      / (here->HSM1_weff * here->HSM1_leff);
	    I *= (I < 0.0) ? -1.0 : 1.0;
	    I *= here->HSM1_mu;
	    I *= here->HSM1_m; /* PN */
	    NevalSrc(&noizDens[HSM1IDNOIZ], &lnNdens[HSM1IDNOIZ], 
		     ckt, THERMNOISE, 
		     here->HSM1dNodePrime, here->HSM1sNodePrime, I);
	    break;
	  case 5:
	    NevalSrc(&noizDens[HSM1IDNOIZ], &lnNdens[HSM1IDNOIZ], 
		     ckt, THERMNOISE, 
		     here->HSM1dNodePrime, here->HSM1sNodePrime, 0.0);
	    break;
	  }
	  NevalSrc(&noizDens[HSM1FLNOIZ], (double*) NULL,
		   ckt, N_GAIN, 
		   here->HSM1dNodePrime, here->HSM1sNodePrime, 
		   (double) 0.0);
	  
	  /* flicker noise */
	  switch ( model->HSM1_noise ) {
	  case 1:
	  case 4: /* SPICE2 model */
	    noizDens[HSM1FLNOIZ] *= model->HSM1_kf
	      * exp(model->HSM1_af * log(MAX(fabs(here->HSM1_ids * here->HSM1_m), N_MINLOG)))
	      / (pow(data->freq, model->HSM1_ef) * here->HSM1_leff
		 * here->HSM1_leff * (3.453133e-11 / model->HSM1_tox));
	    /*                        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cox  */
	    break;
	  case 2:
	  case 3:
	  case 5:
	    /* from HiSIM */
	    noizDens[HSM1FLNOIZ] *= here->HSM1_nfc / data->freq; 
	    break;
	  }
	  
	  lnNdens[HSM1FLNOIZ] = log(MAX(noizDens[HSM1FLNOIZ], N_MINLOG));
	  
	  noizDens[HSM1TOTNOIZ] = noizDens[HSM1RDNOIZ] + noizDens[HSM1RSNOIZ]
	    + noizDens[HSM1IDNOIZ] + noizDens[HSM1FLNOIZ];
	  lnNdens[HSM1TOTNOIZ] = log(MAX(noizDens[HSM1TOTNOIZ], N_MINLOG));
	  
	  *OnDens += noizDens[HSM1TOTNOIZ];
	  
	  if ( data->delFreq == 0.0 ) {
	    /* if we haven't done any previous 
	       integration, we need to initialize our
	       "history" variables.
	    */
	    
	    for ( i = 0; i < HSM1NSRCS; i++ ) 
	      here->HSM1nVar[LNLSTDENS][i] = lnNdens[i];
	    
	    /* clear out our integration variables
	       if it's the first pass
	    */
	    if (data->freq == ((NOISEAN*) ckt->CKTcurJob)->NstartFreq) {
	      for (i = 0; i < HSM1NSRCS; i++) {
		here->HSM1nVar[OUTNOIZ][i] = 0.0;
		here->HSM1nVar[INNOIZ][i] = 0.0;
	      }
	    }
	  }
	  else {
	    /* data->delFreq != 0.0,
	       we have to integrate.
	    */
	    for ( i = 0; i < HSM1NSRCS; i++ ) {
	      if ( i != HSM1TOTNOIZ ) {
		tempOnoise = 
		  Nintegrate(noizDens[i], lnNdens[i],
			     here->HSM1nVar[LNLSTDENS][i], data);
		tempInoise = 
		  Nintegrate(noizDens[i] * data->GainSqInv, 
			     lnNdens[i] + data->lnGainInv,
			     here->HSM1nVar[LNLSTDENS][i] + data->lnGainInv,
			     data);
		here->HSM1nVar[LNLSTDENS][i] = lnNdens[i];
		data->outNoiz += tempOnoise;
		data->inNoise += tempInoise;
		if ( ((NOISEAN*)ckt->CKTcurJob)->NStpsSm != 0 ) {
		  here->HSM1nVar[OUTNOIZ][i] += tempOnoise;
		  here->HSM1nVar[OUTNOIZ][HSM1TOTNOIZ] += tempOnoise;
		  here->HSM1nVar[INNOIZ][i] += tempInoise;
		  here->HSM1nVar[INNOIZ][HSM1TOTNOIZ] += tempInoise;
		}
	      }
	    }
	  }
	  if ( data->prtSummary ) {
	    for (i = 0; i < HSM1NSRCS; i++) {
	      /* print a summary report */
	      data->outpVector[data->outNumber++] = noizDens[i];
	    }
	  }
	  break;
	case INT_NOIZ:
	  /* already calculated, just output */
	  if ( ((NOISEAN*)ckt->CKTcurJob)->NStpsSm != 0 ) {
	    for ( i = 0; i < HSM1NSRCS; i++ ) {
	      data->outpVector[data->outNumber++] = here->HSM1nVar[OUTNOIZ][i];
	      data->outpVector[data->outNumber++] = here->HSM1nVar[INNOIZ][i];
	    }
	  }
	  break;
	}
	break;
      case N_CLOSE:
	/* do nothing, the main calling routine will close */
	return (OK);
	break;   /* the plots */
      }       /* switch (operation) */
    }    /* for here */
  }    /* for model */
  
  return(OK);
}



