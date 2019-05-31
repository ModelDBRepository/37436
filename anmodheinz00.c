/* anmodheinz00.C : Code for AN Model from Heinz et al. (2000) */

/* Cleaned up from Manuscript code 9/10/99 by M.G. Heinz
   (mgheinz@mit.edu) */

/* Basic multiple-fiber LINEAR, Human auditory-nerve model with all
   model-stage outputs saved to ascii files in the following format: 

 	  DATA FILES set up as follows: 
 	    line1: MAXCHS [time_vector] 
 	    line2: Cf1 [ifr_Cf1] 
 	    line3: Cf2 [ifr_Cf2] 
 	    ... 

   The variables: savestim, savebm, saveihc, saveifr can be set to 0
   to avoid saving the stimulus, basilar membrane, inner hair cell. or
   instantaneous firing rate, respectively.  */
	  
/* The model and physiological responses of the model are described in
   the file: model_descript.ps, which is an excerpt from the
   manuscript.  Filter bandwidths are based on human psychophysical
   tuning curves, and DO NOT change with level (some code is left in
   this version that allows the filters to vary with level [set
   health=1], but the NONLINEAR version has not been tested.  Only the
   LINEAR version of the model was used in Heinz et al., 2000.

   Basic model response: Rsat=210 sp/sec, Rate Threshold=0 dB SPL,
   DR=20 dB, Max Onset rate = 850 sp/sec, max synch at LFs = 0.8,
   rolloff matches Johnson (1980). Synch Threshold = -13 dB. PSTs look
   good.  */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>

#define MAXTIMPTS 30000L   /* Maximum number of samples in time */
#define MAXCHS 60          /* Maximum number of frequency channels */

#define PIe   3.1415926536  

#define DOFOR(i, to) for(i = 0; i < to; i++)

#define CMULTR(X,Y) ((X).x*(Y).x-(X).y*(Y).y)
#define CMULTI(X,Y) ((X).y*(Y).x+(X).x*(Y).y)
#define CTREAL(z,X,re) {(z).x=(X).x*(re);(z).y=(X).y*(re);}
#define CMULT(z,X,Y) {(z).x=CMULTR((X),(Y));(z).y=CMULTI((X),(Y));}
#define CADD(z,X,Y) {(z).x=(X).x+(Y).x;(z).y=(X).y+(Y).y;}

long stimtype,numchs,health,nstimpts;
long savestim,savebm,saveihc,saveifr;
long ichan,i,resp;
double xdum,Cf,delx,x,SPER;
double lowcf,xlowcf,highcf,xhighcf;

double stim[MAXTIMPTS],gtf[MAXTIMPTS],ihcL[MAXTIMPTS],ihc2[MAXTIMPTS];
double ifr[MAXTIMPTS];   /* channel x time */
double anfCfs[MAXCHS]; /* chan: anf: ANFS */
double freq,phase,stimrft,levdbS,stimdur,Textra;

/* file pointers for writing to files */
FILE *fpstim,*fpbm,*fpihc,*fpifr;

void error(char *fmt, ...);
void stimulus();
void gamma4();
void synapse();
double erbGM(double);
double cmaph_x2f(double);
double cmaph_f2x(double);

struct complex { double x; double y;} ;
struct complex compexp(double);
struct complex compmult(double, struct complex);
struct complex compprod(struct complex, struct complex);
struct complex comp2sum(struct complex, struct complex);
double REAL(struct complex);

long round(double);

void main(void)
{

  /* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
  printf("\n Heinz, Colburn, and Carney (2000) Auditory-Nerve Model (Linear, Human) - \n\tVersion 1.0 (9/10/99)\n");
  printf("This version reads in files with <= 30,000 pts. \n");

  savestim=1;     /* Boolean save variables: 0=NO, 1=YES */
  savebm=1;
  saveihc=1;
  saveifr=1;
  /* set these flag to write out the output at different stage of the
     model to files for plotting (slows down program). */
	
  /* initialize variables */
  printf("Filterbank parameters: Set Center-filter freq [0], or Range of freqs [1]?: "); scanf("%ld",&resp); printf(" %d\n",resp);
  if(resp == 0)    {
    printf("CF of center filter in bank (Hz): ");scanf("%lf", &Cf); printf(" %f\n",Cf);
    printf("# of Filters in Bank: ");scanf("%ld", &numchs);  printf(" %d\n",numchs);
    printf("Distance between filters (mm along BM)- ");scanf("%lf", &delx);  printf(" %f\n",delx);  
    /* NOTE: spacing between filters of 0.01mm ~= 1 HC spacing along BM */
    xlowcf=cmaph_f2x(Cf)-(numchs-1.)/2.*delx;
    xhighcf=xlowcf+(numchs-1)*delx;
    lowcf=round(cmaph_x2f(xlowcf));
    highcf=round(cmaph_x2f(xhighcf));
    printf("Resulting Filter Bank Center Frequency Range:  %.2f - %.2f Hz\n",lowcf,highcf);
  }
  else  {
    printf("Lowest frequency filter(Hz): ");scanf("%lf", &lowcf); printf(" %f\n",lowcf);
    printf("Highest frequency filter(Hz): ");scanf("%lf", &highcf); printf(" %f\n",highcf);
    printf("# of Filters in Bank: ");scanf("%ld", &numchs);  printf(" %d\n",numchs);
    if(numchs>MAXCHS) {
      printf("\n*** MAXCHS=%d in code. Need to increase!\n\n",MAXCHS);
      error("STOPPED");
    }
    xlowcf=cmaph_f2x(lowcf);
    xhighcf=cmaph_f2x(highcf);
    delx = (xhighcf-xlowcf)/(numchs-1);
    printf("Resulting delx = %f [Distance between filters (mm along BM)]\n",delx);  
  }
  
  
  /* ***** CHANNEL CFs SETUP (anf) ***** */
  /* Set up ANF channels */
  xdum=xlowcf;
  anfCfs[0]=round(cmaph_x2f(xlowcf));
/*    printf("Channel %d = %4.2f Hz\n",ichan+1,anfCfs[0]);     */
  DOFOR(ichan,numchs-1) {
    xdum=xdum+delx;
    anfCfs[ichan+1]=round(cmaph_x2f(xdum));
/*      printf("Channel %d = %4.2f Hz\n",ichan+2,anfCfs[ichan+1]);     */
  }

  stimulus();

  health=-1;
  /* (1=NL(healthy);
     0=impaired(LIN,broad,high threshold);
     -1=Previous Analytical Models(LIN,narrow,low threshold)) 
     -2=To isolate BW effect (LIN,broad,low threshold)) */
  printf("Health [-1=LIN(sharp,Low Thresh); 1=NL; 0=LIN(broad,HT); -2:LIN(broad,LT)]:  %d\n",health);
  
  /* Save2file -- Stimulus */
  if(savestim) {
    fpstim = fopen("stim.dat","w");
    fprintf(fpstim, "%.4f",freq);
    DOFOR(i,nstimpts) fprintf(fpstim, " %e",i*SPER*1e3);
    fprintf(fpstim,"\n");
    fprintf(fpstim, "%.4f",levdbS);
    DOFOR(i,nstimpts) fprintf(fpstim," %e",stim[i]);   
    fprintf(fpstim,"\n");
    fclose(fpstim);
  } 

  /* Save2file -- BM Filter Outputs */
  if(savebm) {
    fpbm = fopen("bm.dat","w"); 
    fprintf(fpbm,"%d",numchs);
    DOFOR(i,nstimpts) fprintf(fpbm, " %f",i*SPER*1e3); 
    fprintf(fpbm,"\n"); 
  }
  
  /* Save2file -- IHC Outputs */
  if(saveihc) {
    fpihc = fopen("ihc.dat","w"); 
    fprintf(fpihc,"%d",numchs);  
    DOFOR(i,nstimpts) fprintf(fpihc, " %f",i*SPER*1e3); 
    fprintf(fpihc,"\n"); 
  }
  
  /* Save2file -- IFR */ 
  if(saveifr) {
    fpifr = fopen("ifr.dat","w");
    fprintf(fpifr,"%d",numchs);
    DOFOR(i,nstimpts) fprintf(fpifr, " %e",i*SPER*1e3);  
    fprintf(fpifr,"\n");   
  }
  
  /* Channel Loop */
  
  DOFOR(ichan,numchs) {
    Cf = anfCfs[ichan];
/*      x=cmaph_f2x(Cf); */
    x = 11.9 * log10(.8 + (Cf/456.));
    printf("chan = %d out of %d  Cf= %f\n",ichan+1,numchs,Cf); 

    /* Calculate IFR responses for stimulus */
    gamma4();
    synapse(); 
    
    /* Save2file -- BM Filter Outputs */
    if(savebm) {
      fprintf(fpbm,"%f",anfCfs[ichan]);
      DOFOR(i,nstimpts) fprintf(fpbm," %e",gtf[i]);
      fprintf(fpbm,"\n");
    }
    
    /* Save2file -- IHC Outputs */
    if(saveihc) {
      fprintf(fpihc,"%f",anfCfs[ichan]);
      DOFOR(i,nstimpts) fprintf(fpihc," %e",ihcL[i]);
      fprintf(fpihc,"\n");
    }
    
    /* Save2file -- IFR */
    if(saveifr) {
      fprintf(fpifr,"%f",anfCfs[ichan]);   
      DOFOR(i,nstimpts) fprintf(fpifr," %e",ifr[i]);   
      fprintf(fpifr,"\n");   
    }
    
  } /* end ichan loop */
  
  if(savebm)      fclose(fpbm); 
  if(saveihc)      fclose(fpihc); 
  if(saveifr)     fclose(fpifr); 
  
  printf("\n AN filter bank simulation is complete. \n");
  
}  /* end of main() */

void stimulus(void)
{
/* stimulus amplitudes are in "pascals" 
     Conversion: 90dB SPL re 20uPa = .632 Pa */
  
char waveform_file[40];
  int i;
  
  printf("Stimulus Waveform Filename: ");
  scanf("%s",waveform_file); printf(" %s\n",waveform_file);
  printf("Time Step Size in input file (secs): ");
  scanf("%lf", &SPER); printf(" %f\n",SPER);
  printf("Input Duration of Simulation (in msec): ");
  scanf("%lf", &stimdur); printf(" %f\n",stimdur);
  printf("  [Zeroes will be added to end of input waveform out to this duration.]\n ");
  /* C convert to seconds */
  stimdur=stimdur*1e-3;
  
  nstimpts=round(stimdur/SPER);
  DOFOR(i,nstimpts) stim[i] = 0.;     /* zero out buffer */
  
  printf("\nReady to read in from the waveform file: %s\n  First ten values are:\n",waveform_file);
  
  fpstim = fopen(waveform_file,"r");
  DOFOR(i,nstimpts) {
    /* For single column, ascii, floating point input waveform, use next line */
		fscanf(fpstim,"%lf",&stim[i]);
		if(i < 10) printf("i = %d  stim[i] = %f\n",i,stim[i]);
    }
  printf("\n");
  fclose(fpstim);
  return;
}  /* End stimulus() */


/* NEW GAMMA4 5/22/98 */
void gamma4(void)
{
  long i,j,idelay;
  double A0,A1,ss0,cc0,ss1,cc1;
  double c, Fc, c1LPihc, c2LPihc, c1LPfb=0, c2LPfb=0;
  double fb = 0., fbl = 0., fbtemp = 0., fbtempl = 0.;
  double wavenow = 0.,wavel = 0.;
  double x = 0., tauLL = 0.;
  double delay = 0., tau0 = 0., tau = 0., Kihc=0, Kfb=0.;
  double asymihc=0., taurange=0, betafb=0;
  double asymfb=0, betaihc=0;
  double Fcfb=0, DC=0.;
  struct complex gtf2[5],gtf2l[5];  


/* Initialization */
  for(i = 0; i < 5; i++)
    {
      gtf2[i].x = 0.; gtf2[i].y = 0.; gtf2l[i].x = 0.; gtf2l[i].y = 0.;
    }
  
  /* parameters for Tau0 vs. CF */
  ss0 = 6.;
  cc0 = 1.1;
  ss1 = 2.2;
  cc1 = 1.1;
  
  c = 2. / SPER;  /* for Bilinear transformation */
/* IHC LP filter parameters */
  Fc = 4800.; /* Fc is nominal cutoff freq (i.e., -3*n dB down point,
		 n=order) for IHC low-pass filters (Hz) to obtain a
		 3-dB cutoff frequency of 2500 Hz according to Weiss
		 and Rose (1988), and to match Johnson (1980) by
		 eye. */
  c1LPihc = ( c - 2 * PIe * Fc ) / ( c + 2 * PIe * Fc );
  c2LPihc = 2 * PIe * Fc / (2 * PIe * Fc + c );

/* FB LPF parameters */  
  Fcfb = 1./(2 * PIe * .002); /* for 2 msec time constant in FB */
  c1LPfb = ( c - 2 * PIe * Fcfb ) / ( c + 2 * PIe * Fcfb ); 
  c2LPfb = 2 * PIe * Fcfb / (2 * PIe * Fcfb + c );
 
/* Find tau0 for this CF */
  /* x = 11.9 * log10(0.80 + Cf / 456.); */ /* position of cf unit;
      from Liberman's map */
  /* Cat filters */ 
  /* tau0 = ( cc0 * exp( -x / ss0) + cc1 * exp( -x /ss1) ) * 1e-3; */ 
  /* in sec */
  
  /* This is the setting for HUMAN filters. Bandwidths are based on
     Glasberg and Moore's (1990) ERB=f(CF) equation.  These values of
     ERB are used at low levels, and then the "high-level" ERB
     corresponding to tau0 is set to be twice as wide. */
  /* tau is (2pi*(1.019*ERB)) */
  /* tau0 in sec */
  /* Gammatone bandwidth is 1/(2*pi*tau) */
  tauLL = 1./(2*PIe*1.019*erbGM(Cf));

/* Set parameters for IHC and FB nonlinearities */
  asymihc = 3; /* asymmetry - this sets positive:negative asymmetry of
                  ihc NL*/
  betaihc = tan(PIe * (-0.5 + 1./(asymihc + 1.))); /* used below to bias NL */
  Kihc = 1225.;
  /* Gain on input of NL- effectively determines threshold & dynamic
     range of IHC - set this using anrhode - saturates at ~60 dB SPL
     (see Dallos) - also influences alpha for rate-level function */

  asymfb = 3; /* asymmetry  - this sets positive:negative asymmetry of ihc NL*/
  betafb = tan(PIe * (-0.5 + 1./(asymfb + 1.)));
  Kfb = 3000.; /* set using anrhode - this determines threshold and
                  range of compression*/
  taurange = 0.5 * tauLL; /*range of tau variation - determines
			  'strength' & influences threshold of
			  compression. Adjust using anrhode and anra
			  (phase) */
  tau0=tauLL-taurange;

  DC = (1.-1./asymfb)/2.; /* - this is asymptotic DC when max of NL is normed to 1 */ 

  /* i=0 */  
  gtf2l[0] = compmult( stim[0], compexp( -2*PIe * Cf * SPER));  /* init */

  for(i = 1; i < nstimpts; i++)   /* Time Loop */
    {	
      /* FREQUENCY SHIFT THE ARRAY BUF  */
      gtf2[0] = compmult(stim[i], compexp( -2*PIe * Cf * SPER * i));
      if(health == 1) tau = tau0 + taurange * (DC - fbl)/DC ; 
      if(health == 0) tau = tau0;
      if(health == -1) tau = tau0 + taurange;
      if(health ==-2) tau = tau0;

      for(j = 1; j < 5; j++)  /* IIR Bilinear transformation LPF */
	gtf2[j] = comp2sum(compmult(1./(tau*c+1.),comp2sum(gtf2[j-1],gtf2l[j-1]) ), 
			  compmult((tau*c-1.)/(tau*c+1.),gtf2l[j])); 
      
      /* FREQUENCY SHIFT BACK UP   */
      /* Factor of tau  put in front of filter 11/26/97  - normalization by tau0 included for now */
      gtf[i] = tau*tau*tau*tau/(tau0*tau0*tau0*tau0) 
	* REAL(compprod(compexp(2*PIe * Cf * SPER * i), gtf2[4]));
      
      if(health==-2) gtf[i]=gtf[i]*(tau0+taurange)*(tau0+taurange)*(tau0+taurange)*(tau0+taurange)/(tau*tau*tau*tau);
      
      fbtemp = gtf[i]; /* filter output used for feedback */

      /* IHC NL */
      wavenow = (atan(Kihc * gtf[i] + betaihc) 
		 - atan(betaihc))/(PIe/2. - atan(betaihc));
      
      /* FB NL */
      fbtemp = (atan(Kfb * fbtemp + betafb) 
		 - atan(betafb))/(PIe/2. - atan(betafb));
      
      /* The following LPFs are IIR Bilinear transformation filters */
      ihcL[i] = c1LPihc * ihcL[i-1] 
	+ c2LPihc * (wavenow + wavel);  /* lp filter the IHC*/

      fb = c1LPfb * fbl 
	+ c2LPfb * (fbtemp + fbtempl);  /* lp filter the fb tau signal*/

      /* save all loop parameters */
      for(j = 0; j < 5; j++) gtf2l[j] = gtf2[j];
      wavel = wavenow;	
      fbl = fb;
      fbtempl = fbtemp;
    }    /* END of TIME LOOP */ 
  

  /* lowpass filter the IHC voltage  more (these could be merged
     with loop above....) */
  /* There's no need to keep entire arrays for the intermediate IHC
     signals, once everything's debugged */

  for( i = 1; i < nstimpts; i++) 
    ihc2[i] = c1LPihc * ihc2[i-1] + c2LPihc * (ihcL[i] + ihcL[i-1]); 

  for( i = 1; i < nstimpts; i++) 
    ihcL[i] = c1LPihc * ihcL[i-1] + c2LPihc * (ihc2[i] + ihc2[i-1]); 

  for( i = 1; i < nstimpts; i++) 
    ihc2[i] = c1LPihc * ihc2[i-1] + c2LPihc * (ihcL[i] + ihcL[i-1]); 

  for( i = 1; i < nstimpts; i++) 
    ihcL[i] = c1LPihc * ihcL[i-1] + c2LPihc * (ihc2[i] + ihc2[i-1]); 

  for( i = 1; i < nstimpts; i++) 
    ihc2[i] = c1LPihc * ihc2[i-1] + c2LPihc * (ihcL[i] + ihcL[i-1]); 

  for( i = 1; i < nstimpts; i++) 
    ihcL[i] = c1LPihc * ihcL[i-1] + c2LPihc * (ihc2[i] + ihc2[i-1]); 
  
  /* DELAY THE WAVEFORM (delay gtf and ihcL for display purposes)  */
  /* Note: Latency vs. CF for click responses is available for Cat only (not human) */ 
  /* Use original fit for Tl (latency vs. CF in msec) from Carney & Yin '88 
     and then correct by .75 cycles to go from PEAK delay to ONSET delay */
/*    A0 = 8.13;  */ /* from Carney and Yin '88 */
/*    A1 = 6.49; */
/*    delay = A0 * exp( -x/A1 ) * 1e-3 - 1./Cf;    */
/*    printf("delay=%e\n",delay); */
/*    idelay = delay / SPER; */
/*    for(i = nstimpts; i > idelay; i--) { */
/*      gtf[i] = gtf[i-idelay];  */
/*      ihcL[i] = ihcL[i-idelay];  */
/*    } */
/*    for(i = 1; i < (idelay + 1); i++) { */
/*      gtf[i] = 0.;  */
/*      ihcL[i] = 0.;  */
/*    } */
  /*   printf("No delay\n");  */ /*reminder message, if delay is commented out*/
  
  return;
}  /* End of gamma4() */


/* NEW SYNAPSE 5/22/99 */
void synapse(void)
{
  long i,j,isp;
  double PIrest,PPI,PImax,PL,PG,CI,CL,CG,VI,VL,Pfactor,Vfactor;
  double c0,s0,c1,s1,dead,rtime,rsptime,rint,prob;
  double g,spont,ftemp,Rsat;
  double p1,p3;
  long option; /* Option for PPI: 1: half-wave rectify, 2: NL Smoother */

  spont = 50.;  /* "spont" rate, before adaptation, refractoriness */
/*    Rsat=165.; */

  Pfactor = .03;  /* ** Controllable ** This might be scaled later to
		     go from inst. rate mode to spikes mode*/
  PIrest = Pfactor/2.5;  /* .012; this will be lower for spikes version */
  PG = Pfactor;        /* .03 */
  PL = Pfactor * 2.;   /*.06 roughly from W&S 1988 figs */

  Vfactor = 0.0005;   /* **Controllable**  */
  VI = Vfactor;
  VL = 10. * Vfactor; 

  CI = spont / PIrest;
  CL = CI * (PIrest + PL)/PL;  /* for stability in steady-state */
  CG = CL * (1. + PL/PG) - CI * PL / PG;  /* so that system is in steady
					     state at spont */
/*    PImax=PL*PG/(CG*PL*PG/Rsat-PG-PL); */   /* 0.18 newsyn5 */
  PImax=0.6;

  ifr[0] = spont;

  option=2;

  for( i = 1; i < nstimpts ; i ++) 
    {

      if (option==1) {
      /* Option 1: Linear equation between ihcL (in range [-1/3,1]) and
	 PPI, the half-wave rectify PPI */

	PPI = (PImax-PIrest) * ihcL[i] + PIrest;
	if (PPI<0.0) PPI=0.0;
      }
      
      if (option==2) {
      /* Option 2: NonLinear relation between ihcL (in range [-1/3,1])
	 and PPI, such that PImax and PIrest are achieved and the PPI
	 goes to 0 smoothly for negative ihcL */

	p1=log(exp(log(2)*PImax/PIrest)-1);
	p3=p1*PIrest/log(2);

	PPI = p3/p1*log(1+exp(p1*ihcL[i]));
      }

      CI = CI + (SPER/VI)*(-PPI*CI + PL*(CL - CI));
      CL = CL + (SPER/VL)*(-PL*(CL - CI) + PG*(CG - CL));
      ifr[i] = CI * PPI;
    }
  
  /* Now, ifr[i] contains instantaneous discharge rate vs. time */
  
  return;
}  /* End of synapse() */
  

long round(double value)
{
  if((value-floor(value))>=0.5) return(ceil(value));
  else return(floor(value));
}

double erbGM(double CF)
{
  double erbCf;

  erbCf=24.7*(4.37*Cf/1000+1);

  return(erbCf);
}

double cmaph_f2x(double f)
{
  double x;
  
  if((f>20677)||(f<20)) error("frequency out of human range, [in cmaph_f2x(f)]");
  x=(1.0/0.06)*log10((f/165.4)+0.88);
  return(x);
}

double cmaph_x2f(double x)
{
  double f;
  
  if((x>35)||(x<0)) error("BM distance out of human range, [in cmaph_x2f(x)]");
  f=165.4*(pow(10,(0.06*x))-0.88);
  return(f);
}

struct complex compexp(double theta)
     /* this returns a complex number equal to exp(i*theta)  */
{
  struct complex answer;
  
  answer.x = cos(theta);
  answer.y = sin(theta);
  return answer;
}

struct complex compmult(double scalar, struct complex compnum)
     /* Multiply a complex number by a scalar */
{
  struct complex answer;
  
  CTREAL(answer,compnum,scalar);
  return answer;
}

struct complex compprod(struct complex compnum1, struct complex compnum2)
     /* Find the product of 2 complex numbers */
{
  struct complex answer;

  CMULT(answer,compnum1,compnum2);
  return answer;
}

struct complex comp2sum(struct complex summand1, struct complex summand2)
     /* add 2 complex numbers  */
{
  struct complex answer;

  CADD(answer,summand1,summand2);
  return answer;
}

double REAL(struct complex compnum)
{
  return compnum.x;
}

/* error: print an error message and die gracefully */
/* Takes arguments like printf */
/* Copied from Kernighan and Ritchie, p 174 */
void error(char *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  fprintf(stderr, "error: ");
  vfprintf(stderr, fmt, args);
  fprintf(stderr, "\n");
  va_end(args);
  exit(1);  /* closes all open file */
}
