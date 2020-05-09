/****************************************************************/
/*								*/
/*		aok4.c	         (Version 4.0) 		        */
/*								*/
/*	        Adaptive Optimal-Kernel (AOK)                   */
/*              Time-Frequency Representation                   */
/*								*/
/*		Douglas L. Jones (author)			*/
/*		  University of Illinois			*/
/*		  E-mail: jones@uicsl.csl.uiuc.edu		*/	
/*								*/
/*		Richard G. Baraniuk (facilitator)		*/	
/*	  	  Rice University				*/
/*		  E-mail: richb@rice.edu			*/
/*								*/
/*	Written:	December 26, 1991  (version 1.0) 	*/
/* 	Modified:	June 10, 1992	   (version 2.0)	*/
/*		  	November 20, 1992  (version 3.0) 	*/
/*			February 8, 1992   (version 4.0)	*/
/*                      January 28, 1996   (renamed aok)        */
/*								*/
/****************************************************************/
/*								*/
/*	This version interpolates the polar STAF from		*/
/*	the rectangular STAF.  It implicitly applies a		*/
/*	rectangular window to the data to create the		*/
/*	short-time ambiguity function, and includes		*/
/*	all non-zero time lags in the STAF.			*/
/*
*-----------------------------------------------------------------
* Copyright (C) 1992, 1993, 1994, 1995, 1996 the Board of Trustees of
* the University of Illinois.  All Rights Reserved.  Permission is
* hereby given to use, copy, modify, and distribute this software
* provided that (1) the headers, copyright and proprietary notices are
* retained in each copy and (2) any files that are modified are
* identified as such (see below).  The University of Illinois makes no
* representations or warranties of any kind concerning this software or
* its use.
*
* Any modifications made to this file must be commented and dated 
* in the following style:
*
*  Source file:   	aok4.c
*  Modifications:	Richard Baraniuk, November 25, 1992
*       		  Inserted this sample edit history
*			Douglas Jones, February 8, 1993
*			  Implemented gradient-project algorithm
*			  in terms of sigma rather than sigma^2
*                       Richard Baraniuk, January 29, 1996
*                         Renamed runrgk --> aok
*
*       Please log any further modifications made to this file:
*
*---------------------------------------------------------------*/


#include <stdio.h>
#include <math.h>	      /* Compile with " cc aok4.c -lm " */

main()
{
	FILE	*fopen(), *ofp, *setoutfile();
	FILE	*ifp, *setinfile();


	int	i,j,k,ii;
	int	itemp;
	int	xlen, tlen, nraf,tlag,nlag,mfft,nrad, nphi, nits;
	int	slen,fftlen;
	int	tstep,fstep,outct;

	int	maxrad[1024];	/* max radius for each phi s.t. theta < pi */

	
	char	name[10];


	double	pi, rtemp, rtemp1, rtemp2, mu, forget;
	double	vol;			/* total kernel volume	*/
	double	outdelay;		/* delay in samples of output time from current sample time	*/
	double	xr[1024],xi[1024];	/* recent data samples	*/
	double	rectafr[70000];		/* real part of rect running AF	*/
	double	rectafi[70000];		/* imag part of rect running AF	*/
	double  rectafm2[70000];	/* rect |AF|^2	*/
	double  polafm2[70000];		/* polar |AF|^2	*/
	double	rectrotr[70000];	/* real part of rect AF phase shift	*/
	double	rectroti[70000];	/* imag part of rect AF phase shift	*/
	double	req[70000];		/* rad corresp. to rect coord	*/
	double  pheq[70000];		/* phi corresp. to rect coord	*/
	double	plag[70000];		/* lag index at polar AF sample points	*/
	double	ptheta[70000];		/* theta index at polar AF sample points	*/
	double  sigma[1024];		/* optimal kernel spreads	*/
	double  rar[1024],rai[1024];	/* poles for running FFTs for rect AF	*/
	double  rarN[70000];		/* poles for running FFTs for rect AF	*/
	double  raiN[70000];
	double  tfslicer[1024];		/* freq slice at current time	*/
	double  tfslicei[1024];

extern	double	ccmr(),ccmi(),cmr(),cmi(),rectkern();

/********************************/
/*   initialization and setup   */
/********************************/

#define	ALPHA	0.01;

        printf("\n");           /* Print header */
        printf("ADAPTIVE OPTIMAL-KERNEL (AOK) TIME-FREQUENCY REPRESENTATION\n");

        printf("   \n");
        printf("   Version 4.0\n");
        printf("   \n");

	ifp =  setinfile();
	ofp =  setoutfile();

	printf("Number of signal samples in input file?\n");    /* get input parameters */
	scanf("%d", &xlen);

        printf("Length of sliding analysis window?  (power of two, no larger than 256)\n");
	printf("  (Number of samples along each dimension of the STAF)\n");
	scanf("%d", &tlag);

	printf("Number of output frequency samples per time-slice?  (power of two)\n", (2*tlag));
	scanf("%d", &fftlen);

	if ( fftlen < (2*tlag) )
	   {
	    fstep = 2*tlag/fftlen;
	    fftlen = 2*tlag;
	   }
	  else
	   {
	    fstep = 1;
	   }

	printf("Time increment in samples between time-slice outputs?\n");
	scanf("%d", &tstep);


        printf("Normalized volume of optimal kernel?\n");
        printf("  (Typically between 1 and 5)\n");
	scanf("%lf", &vol);
        printf("  \n");			    /* Skip a line */
/*	fprintf(ofp, "vol = %g ; \n", vol); */


/*	tlag = 64; */		/* total number of rectangular AF lags	*/
	nits = (int) log2((double) tstep+2);	/* number of gradient steps to take each time	*/
/*	nits = 2; */
/*	vol = 2.0; */		/* kernel volume (1.0=Heisenberg limit)	*/
	mu = 0.5;		/* gradient descent factor	*/

	forget = 0.001;		/* set no. samples to 0.5 weight on running AF	*/
	nraf = tlag;		/* theta size of rectangular AF	*/
	nrad = tlag;		/* number of radial samples in polar AF	*/
	nphi = tlag;		/* number of angular samples in polar AF */
	outdelay = tlag/2;	/* delay in effective output time in samples	*/
				/* nlag-1 < outdelay < nraf to prevent "echo" effect */

	nlag = tlag + 1;	/* one-sided number of AF lags	*/
	mfft = po2(fftlen);
	slen = 1.42*(nlag-1) + nraf + 3;	/* number of delayed samples to maintain	*/

	pi = 3.141592654;
	vol = (2.0*vol*nphi*nrad*nrad)/(pi*tlag);	/* normalize volume	*/


	kfill((nrad*nphi),0.0,polafm2);
	kfill((nraf*nlag),0.0,rectafr);
	kfill((nraf*nlag),0.0,rectafi);
	kfill(slen,0.0,xr);
	kfill(slen,0.0,xi);
	kfill(nphi,1.0,sigma);


/*****************/
/*   main loop   */
/*****************/

	tlen = xlen + nraf + 2;
	rectamake(nlag,nraf,forget,rar,rai,rarN,raiN);	/* make running rect AF parms	*/
	plagmake(nrad,nphi,nlag,plag);
	pthetamake(nrad,nphi,nraf,ptheta,maxrad);	/* make running polar AF parms	*/
	rectrotmake(nraf,nlag,outdelay,rectrotr,rectroti);
	rectopol(nraf,nlag,nrad,nphi,req,pheq);


	outct = 0;

	for (ii=0; ii < tlen; ii++)
	 {
	  cshift(slen,xr);
	  cshift(slen,xi);
	  if ( ii < xlen )				/* get data	*/
	    fscanf(ifp,"%lf %lf", &xr[0], &xi[0]);
	  else
	   {
	    xr[0] = 0.0;
	    xi[0] = 0.0;
	   }

	  rectaf(xr,xi,nlag,nraf,rar,rai,rarN,raiN,rectafr,rectafi);

	 if ( (ii - (ii/tstep)*tstep) == 0 )	/* output t-f slice	*/
	 {
	  outct = outct + 1;

	  mkmag2((nlag*nraf),rectafr,rectafi,rectafm2);
	  polafint(nrad,nphi,nraf,maxrad,nlag,plag,ptheta,rectafm2,polafm2);
	  sigupdate(nrad,nphi,nits,vol,mu,maxrad,polafm2,sigma);

	  for (i=0; i < nlag-1; i++)	/* for each tau	*/
	   {
	    tfslicer[i] = 0.0;
	    tfslicei[i] = 0.0;


	    for (j = 0; j < nraf; j++)	/* integrate over theta	*/
	     {
	      rtemp = ccmr(rectafr[i*nraf+j],rectafi[i*nraf+j],rectrotr[i*nraf+j],rectroti[i*nraf+j]);
	      rtemp1 =ccmi(rectafr[i*nraf+j],rectafi[i*nraf+j],rectrotr[i*nraf+j],rectroti[i*nraf+j]);

	      rtemp2 = rectkern(i,j,nraf,nphi,req,pheq,sigma);
	      tfslicer[i] = tfslicer[i] + rtemp*rtemp2;
	      tfslicei[i] = tfslicei[i] + rtemp1*rtemp2;
/*	fprintf(ofp," %d , %d , %g, %g, %g , %g , %g \n", i,j,rectafr[i*nraf+j],rectafi[i*nraf+j],rtemp,rtemp1,rtemp2); */
	     }
	   }
	  for (i=nlag-1; i < (fftlen-nlag+2); i++)	/* zero pad for FFT	*/
	   {
	    tfslicer[i] = 0.0;
	    tfslicei[i] = 0.0;
	   }

	  for (i=(fftlen-nlag+2); i < fftlen; i++)	/* fill in c.c. symmetric half of array	*/
	   {
	    tfslicer[i] = tfslicer[fftlen - i];
	    tfslicei[i] = - tfslicei[fftlen - i];
	   }



	  fft(fftlen,mfft,tfslicer,tfslicei);

	  itemp = fftlen/2 + fstep;
	  j = 1;				/* print output slice	*/
	  for (i=itemp; i < fftlen; i=i+fstep)
	   {
/*	    fprintf(ofp, "x( %d , %d )= %g ; \n", outct, j, tfslicer[i]); */
	    fprintf(ofp, "%g ", tfslicer[i]);
	    j = j + 1;
	   }
	  for (i=0; i < itemp; i=i+fstep)
	   {
/*	    fprintf(ofp, "x( %d , %d )= %g ; \n", outct, j, tfslicer[i]); */
	    fprintf(ofp, "%g ", tfslicer[i]);
	    j = j + 1;
	   }
	    fprintf(ofp, "\n");
/*	    printf("outct = %d \n", outct); */
	 }
	 }


/*************************/
/*   close out program   */
/*************************/

	fclose(ifp);
	fclose(ofp);
}
/*								*/
/*   setinfile: opens an input file				*/
/*								*/
FILE	*setinfile()
{
	FILE	*fopen(), *ifp;
	char	filename[40];


	printf ( "Name of input file?\n");
	scanf  ( "%s", filename );
	ifp = fopen( filename, "r" );

	return(ifp);
}
/*								*/
/*   setoutfile: opens an output file				*/
/*								*/
FILE	*setoutfile()
{
	FILE	*fopen(), *ofp;
	char	filename[40];


	printf ( "Name of output file?\n");
	scanf  ( "%s", filename );
	ofp = fopen( filename, "w" );

	return(ofp);
}
/****************************************************************/
/*		fft.c						*/
/*		Douglas L. Jones				*/
/*		University of Illinois at Urbana-Champaign	*/
/*		January 19, 1992				*/
/*								*/
/*   fft: in-place radix-2 DIT DFT of a complex input		*/
/*								*/
/*   input:							*/
/*	n:	length of FFT: must be a power of two		*/
/*	m:	n = 2**m					*/
/*   input/output						*/
/*	x:	double array of length n with real part of data	*/
/*	y:	double array of length n with imag part of data	*/
/*								*/
/*   Permission to copy and use this program is granted		*/
/*   as long as this header is included.			*/
/****************************************************************/
fft(n,m,x,y)
int	n,m;
double	x[],y[];
{
	int	i,j,k,n1,n2;
	double	c,s,e,a,t1,t2;


	j = 0;				/* bit-reverse	*/
	n2 = n/2;
	for (i=1; i < n - 1; i++)
	 {
	  n1 = n2;
	  while ( j >= n1 )
	   {
	    j = j - n1;
	    n1 = n1/2;
	   }
	  j = j + n1;

	  if (i < j)
	   {
	    t1 = x[i];
	    x[i] = x[j];
	    x[j] = t1;
	    t1 = y[i];
	    y[i] = y[j];
	    y[j] = t1;
	   }
	 }


	n1 = 0;				/* FFT	*/
	n2 = 1;

	for (i=0; i < m; i++)
	 {
	  n1 = n2;
	  n2 = n2 + n2;
	  e = -6.283185307179586/n2;
	  a = 0.0;

	  for (j=0; j < n1; j++)
	   {
	    c = cos(a);
	    s = sin(a);
	    a = a + e;

	    for (k=j; k < n; k=k+n2)
	     {
	      t1 = c*x[k+n1] - s*y[k+n1];
	      t2 = s*x[k+n1] + c*y[k+n1];
	      x[k+n1] = x[k] - t1;
	      y[k+n1] = y[k] - t2;
	      x[k] = x[k] + t1;
	      y[k] = y[k] + t2;
	     }
	   }
	 }
	
	    
	return;
}
/*								*/
/*   po2: find the smallest power of two >= input value		*/
/*								*/
int	po2(n)
int	n;
{
	int	m, mm;

	mm = 1;
	m = 0;
	while (mm < n) {
	   ++m;
	   mm = 2*mm;
	}

	return(m);
}
/*								*/
/*   power: compute x^n, x and n positve integers		*/
/*								*/
power(x, n)
int  x,n;
{
	int	i,p;

	p = 1;
	for (i=1; i<=n; ++i)
	     p = p*x;
	return(p);
}
/*								*/
/*   zerofill: set array elements to constant			*/
/*								*/
kfill(len,k,x)
int	len;
double	k,x[];
{
	int	i;

	for (i=0; i < len; i++)
	  x[i] = k;

	return;
}
/*								*/
/*   cshift: circularly shift an array				*/
/*								*/
cshift(len,x)
int	len;
double	x[];
{
	int	i;
	double	rtemp;


	rtemp = x[len-1];

	for (i=len-1; i > 0; i--)
	  x[i] = x[i-1];

	x[0] = rtemp;


	return;
}
/*								*/
/*   cmr: computes real part of x times y			*/
/*								*/
double	cmr(xr,xi,yr,yi)
double	xr,xi,yr,yi;
{
	double	rtemp;

	rtemp = xr*yr - xi*yi;

	return(rtemp);
}
/*								*/
/*   cmi: computes imaginary part of x times y			*/
/*								*/
double	cmi(xr,xi,yr,yi)
double	xr,xi,yr,yi;
{
	double	rtemp;

	rtemp = xi*yr + xr*yi;

	return(rtemp);
}
/*								*/
/*   ccmr: computes real part of x times y*			*/
/*								*/
double	ccmr(xr,xi,yr,yi)
double	xr,xi,yr,yi;
{
	double	rtemp;

	rtemp = xr*yr + xi*yi;

	return(rtemp);
}
/*								*/
/*   ccmi: computes imaginary part of x times y*		*/
/*								*/
double	ccmi(xr,xi,yr,yi)
double	xr,xi,yr,yi;
{
	double	rtemp;

	rtemp = xi*yr - xr*yi;

	return(rtemp);
}
/*								*/
/*   rectamake: make vector of poles for rect running AF	*/
/*								*/
rectamake(nlag,n,forget,ar,ai,arN,aiN)
int	nlag,n;
double	forget,ar[],ai[],arN[],aiN[];
{
	int	i,j;
	double	trig,decay;
	double	trigN,decayN;


	trig = 6.283185307/n;
	decay = exp(-forget);

	for (j=0; j < n; j++)
	 {
	  ar[j] = decay*cos(trig*j);
	  ai[j] = decay*sin(trig*j);
	 }

	for (i=0; i < nlag; i++)
	 {
	  trigN = 6.283185307*(n-i);
	  trigN = trigN/n;
	  decayN = exp(-forget*(n-i));

	  for (j=0; j < n; j++)
	   {
	    arN[i*n+j] = decayN*cos(trigN*j);
	    aiN[i*n+j] = decayN*sin(trigN*j);
	   }
	 }


	return;
}
/*								*/
/*   pthetamake: make matrix of theta indices for polar samples	*/
/*								*/
pthetamake(nrad,nphi,ntheta,ptheta,maxrad)
int	nrad,nphi,ntheta,maxrad[];
double	ptheta[];
{
	int	i,j;
	double	theta,rtemp,deltheta;


	deltheta = 6.283185307/ntheta;

	for (i=0; i < nphi; i++)	/* for all phi ...	*/
	 {
	  maxrad[i] = nrad;

	  for (j = 0; j < nrad; j++)	/* and all radii	*/
	   {
	    theta = - ((4.442882938/nrad)*j)*cos((3.141592654*i)/nphi);
	    if ( theta >= 0.0 )
	       {
	        rtemp = theta/deltheta;
		if ( rtemp > (ntheta/2 - 1) )
		 {
		  rtemp = -1.0;
	          if (j < maxrad[i])  maxrad[i] = j;
		 }
	       }
	      else
	       {
	        rtemp = (theta + 6.283185307)/deltheta;
		if ( rtemp < (ntheta/2 + 1) )
		 {
		  rtemp = -1.0;
	          if (j < maxrad[i])  maxrad[i] = j;
		 }
	       }
		
	    ptheta[i*nrad+j] = rtemp;
	   }
	 }


	return;
}
/*								*/
/*   plagmake: make matrix of lags for polar running AF		*/
/*								*/
plagmake(nrad,nphi,nlag,plag)
int	nrad,nphi,nlag;
double	plag[];
{
	int	i,j;

extern	double	mklag();


	for (i=0; i < nphi; i++)        /* for all phi ...      */
	 {
	  for (j = 0; j < nrad; j++)    /* and all radii        */
	   {
	    plag[i*nrad+j] = mklag(nrad,nphi,nlag,i,j);
	   }
	 }


	return;
}
/*								*/
/*   rectopol: find polar indices corresponding to rect samples	*/
/*								*/
rectopol(nraf,nlag,nrad,nphi,req,pheq)
int	nraf,nlag,nrad,nphi;
double	req[],pheq[];
{
	int	i,j,jt;
	double	pi,deltau,deltheta,delrad,delphi;


	pi = 3.141592654;

	deltau = sqrt(pi/(nlag-1));
	deltheta = 2.0*sqrt((nlag-1)*pi)/nraf;
	delrad = sqrt(2.0*pi*(nlag-1))/nrad;
	delphi = pi/nphi;

	for (i=0; i < nlag; i++)
	 {
	  for (j=0; j < nraf/2; j++)
	   {
	    req[i*nraf +j] = sqrt(i*i*deltau*deltau + j*j*deltheta*deltheta)/delrad;
	    if ( i == 0 )
	      pheq[i*nraf +j] = 0.0;
	    else pheq[i*nraf +j] = (atan((j*deltheta)/(i*deltau)) + 1.570796327)/delphi;
	   }

	  for (j=0; j < nraf/2; j++)
	   {
	    jt = j - nraf/2;
	    req[i*nraf + nraf/2 + j]  = sqrt(i*i*deltau*deltau + jt*jt*deltheta*deltheta)/delrad;
	    if ( i == 0 )
	      pheq[i*nraf + nraf/2 + j] = 0.0;
	    else pheq[i*nraf + nraf/2 + j] = (atan((jt*deltheta)/(i*deltau)) + 1.570796327)/delphi;
	   }
	 }


	return;
}
/*								*/
/*   rectrotmake: make array of rect AF phase shifts		*/
/*								*/
rectrotmake(nraf,nlag,outdelay,rectrotr,rectroti)
int	nraf,nlag;
double	outdelay,rectrotr[],rectroti[];
{
	int	i,j;
	double	twopin;

	twopin = 6.283185307/nraf;


	for (i=0; i < nlag; i++)
	 {
	  for (j=0; j < nraf/2; j++)
	   {
	    rectrotr[i*nraf+j] = cos( (twopin*j)*(outdelay - ((double) i)/2.0 ) );
	    rectroti[i*nraf+j] = sin( (twopin*j)*(outdelay - ((double) i)/2.0 ) );
	   }
	  for (j=nraf/2; j < nraf; j++)
	   {
	    rectrotr[i*nraf+j] = cos( (twopin*(j-nraf))*(outdelay - ((double) i)/2.0 ) );
	    rectroti[i*nraf+j] = sin( (twopin*(j-nraf))*(outdelay - ((double) i)/2.0 ) );
	   }
	 }


	return;
}
/*								*/
/*   rectaf: generate running AF on rectangular grid;		*/
/*	     negative lags, all DFT frequencies			*/
/*								*/
rectaf(xr,xi,laglen,freqlen,alphar,alphai,alpharN,alphaiN,afr,afi)
int	laglen,freqlen;
double	xr[],xi[],alphar[],alphai[],alpharN[],alphaiN[],afr[],afi[];
{
	int	i,j;
	double	rtemp,rr,ri,rrN,riN;
extern	double	cmr(),cmi();
extern	double	ccmr(),ccmi();

	for (i=0; i < laglen; i++)
	 {
	  rr = ccmr(xr[0],xi[0],xr[i],xi[i]);
	  ri = ccmi(xr[0],xi[0],xr[i],xi[i]);

	  rrN = ccmr(xr[freqlen-i],xi[freqlen-i],xr[freqlen],xi[freqlen]);
	  riN = ccmi(xr[freqlen-i],xi[freqlen-i],xr[freqlen],xi[freqlen]);

	  for (j = 0; j < freqlen; j++)
	   {
	    rtemp = cmr(afr[i*freqlen+j],afi[i*freqlen+j],alphar[j],alphai[j]) - cmr(rrN,riN,alpharN[i*freqlen+j],alphaiN[i*freqlen+j]) + rr;
	    afi[i*freqlen + j] = cmi(afr[i*freqlen+j],afi[i*freqlen+j],alphar[j],alphai[j]) - cmi(rrN,riN,alpharN[i*freqlen+j],alphaiN[i*freqlen+j]) + ri;
	    afr[i*freqlen + j] = rtemp;
	   }
	 }


	return;
}
/*								*/
/*   polafint: interpolate AF on polar grid;			*/
/*								*/
polafint(nrad,nphi,ntheta,maxrad,nlag,plag,ptheta,rectafm2,polafm2)
int	nrad,nphi,ntheta,nlag,maxrad[];
double	plag[],ptheta[],rectafm2[],polafm2[];
{
	int	i,j;
	int	ilag,itheta,itheta1;
	double	rlag,rtheta,rtemp,rtemp1;


	for (i=0; i < nphi/2; i++)	/* for all phi ...	*/
	 {
	  for (j = 0; j < maxrad[i]; j++)	/* and all radii with |theta| < pi */
	   {
	    ilag = (int) plag[i*nrad+j];
	    rlag = plag[i*nrad+j] - ilag;

	    if ( ilag >= nlag )
	       {
		polafm2[i*nrad+j] = 0.0;
	       }
	      else
	       {
		itheta = (int) ptheta[i*nrad+j];
		rtheta = ptheta[i*nrad+j] - itheta;

		itheta1 = itheta + 1;
		if ( itheta1 >= ntheta )  itheta1 = 0;

		rtemp =  (rectafm2[ilag*ntheta+itheta1] - rectafm2[ilag*ntheta+itheta])*rtheta + rectafm2[ilag*ntheta+itheta];
		rtemp1 =  (rectafm2[(ilag+1)*ntheta+itheta1] - rectafm2[(ilag+1)*ntheta+itheta])*rtheta + rectafm2[(ilag+1)*ntheta+itheta];
		polafm2[i*nrad+j] = (rtemp1-rtemp)*rlag + rtemp;
	       }
	   }
	 }


	for (i=nphi/2; i < nphi; i++)	/* for all phi ...	*/
	 {
	  for (j = 0; j < maxrad[i]; j++)	/* and all radii with |theta| < pi */
	   {
	    ilag = (int) plag[i*nrad+j];
	    rlag = plag[i*nrad+j] - ilag;

	    if ( ilag >= nlag )
	       {
		polafm2[i*nrad+j] = 0.0;
	       }
	      else
	       {
		itheta = (int) ptheta[i*nrad+j];
		rtheta = ptheta[i*nrad+j] - itheta;

		rtemp =  (rectafm2[ilag*ntheta+itheta+1] - rectafm2[ilag*ntheta+itheta])*rtheta + rectafm2[ilag*ntheta+itheta];
		rtemp1 =  (rectafm2[(ilag+1)*ntheta+itheta+1] - rectafm2[(ilag+1)*ntheta+itheta])*rtheta + rectafm2[(ilag+1)*ntheta+itheta];
		polafm2[i*nrad+j] = (rtemp1-rtemp)*rlag + rtemp;
	       }
	   }
	 }


	return;
}
/*								*/
/*   mklag: compute radial sample lag				*/
/*								*/
double	mklag(nrad,nphi,nlag,iphi,jrad)
int	nrad,nphi,nlag,iphi,jrad;
{
	double	delay;

	delay = ((1.414213562*(nlag-1)*jrad)/nrad)*sin((3.141592654*iphi)/nphi);


	return(delay);
}
/*								*/
/*   rectkern: generate kernel samples on rectangular grid	*/
/*								*/
double	rectkern(itau,itheta,ntheta,nphi,req,pheq,sigma)
int	itau,itheta,ntheta,nphi;
double	req[],pheq[],sigma[];
{
	int	iphi,iphi1;
	double	kern,tsigma;


	iphi = (int) pheq[itau*ntheta + itheta];
	iphi1 = iphi + 1;
	if (iphi1 > (nphi-1))  iphi1 = 0;
	tsigma = sigma[iphi] + (pheq[itau*ntheta + itheta] - iphi)*(sigma[iphi1]-sigma[iphi]);

/*  Tom Polver says on his machine, exp screws up when the argument of */
/*  the exp function is too large */
	kern = exp(-req[itau*ntheta+itheta]*req[itau*ntheta+itheta]/(tsigma*tsigma));

	return(kern);
}
/*								*/
/*   sigupdate: update RG kernel parameters			*/
/*								*/
sigupdate(nrad,nphi,nits,vol,mu0,maxrad,polafm2,sigma)
int	nrad,nphi,nits,maxrad[];
double	vol,mu0,polafm2[],sigma[];
{
	int	ii,i,j;
	double	grad[1024],gradsum,gradsum1,tvol,volfac,eec,ee1,ee2,mu;


	for (ii=0; ii < nits; ii++)
	 {
	  gradsum = 0.0;
	  gradsum1 = 0.0;

	  for (i=0; i < nphi; i++)
	   {
	    grad[i] = 0.0;

	    ee1 = exp( - 1.0/(sigma[i]*sigma[i]) );	/* use Kaiser's efficient method */
	    ee2 = 1.0;
	    eec = ee1*ee1;

	    for (j=1; j < maxrad[i]; j++)
	     {
	      ee2 = ee1*ee2;
	      ee1 = eec*ee1;

	      grad[i] = grad[i] + j*j*j*ee2*polafm2[i*nrad+j];
	     }
	    grad[i] = grad[i]/(sigma[i]*sigma[i]*sigma[i]);

	    gradsum = gradsum + grad[i]*grad[i];
	    gradsum1 = gradsum1 + sigma[i]*grad[i];
	   }

	  gradsum1 = 2.0*gradsum1;
	  if ( gradsum < 0.0000001 )  gradsum = 0.0000001;
	  if ( gradsum1 < 0.0000001 )  gradsum1 = 0.0000001;

	  mu = ( sqrt(gradsum1*gradsum1 + 4.0*gradsum*vol*mu0) - gradsum1 ) / ( 2.0*gradsum );


	  tvol = 0.0;

	  for (i=0; i < nphi; i++)
	   {
	    sigma[i] = sigma[i] + mu*grad[i];
	    if (sigma[i] < 0.5)  sigma[i] = 0.5;
/*	    printf("sigma[%d] = %g\n", i,sigma[i]); */
	    tvol = tvol + sigma[i]*sigma[i];
	   }

	  volfac = sqrt(vol/tvol);
	  for (i=0; i < nphi; i++)  sigma[i] = volfac*sigma[i];
	 }


	return;
}
/*								*/
/*   mkmag2: compute squared magnitude of an array		*/
/*								*/
mkmag2(tlen,xr,xi,xm2)
int	tlen;
double	xr[],xi[],xm2[];
{
	int	i;


	for (i=0; i < tlen; i++)
	 {
	  xm2[i] = xr[i]*xr[i] + xi[i]*xi[i];
	 }

	return;
}
/*								*/
/*   printvec: print out a vector				*/
/*								*/
printvec(ofp,len,vec)
int	len;
double	vec[];
FILE	*ofp;
{
	int	i;


	for (i=0; i < len; i++)
	 {
	  fprintf(ofp, " %d , %g \n", i, vec[i]);
	 }

	return;
}
/*								*/
/*   printmat: print matrix in Matlab format			*/
/*								*/
printmat(ofp,rlen,clen,name,mat)
int	rlen,clen;
char	name[];
double	mat[];
FILE	*ofp;
{
	int	i,j;


	for (i=0; i < rlen; i++)
	 {
	  for (j=0; j < clen; j++)
	   {
	    fprintf(ofp, "%s ( %d , %d )=  %g ; \n", name, (i+1),(j+1), mat[i*clen +j]);
	   }
	 }

	return;
}




