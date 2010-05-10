/*
copied from http://www.gnu.org/software/gsl/manual/html_node/Basis-Splines.html
*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_sf_erf.h>

void geterf_(double *xanderfx) {
  xanderfx[1] = gsl_sf_erf(xanderfx[0]);
}

//spline locations held fixed in Mpc^-1; CMB basically fixed P(k) in these units
void dopksmoothbspline_(double *kvals, double *lnpklinear, double *lnpksmooth, int *npts)  {

	double kmaxsuppress = 0.01*0.7;
	size_t n, ncoeffs, nbreak;
	gsl_bspline_workspace *bw;
	gsl_vector *B;
	gsl_vector *c, *w, *x, *y;
	gsl_matrix *X, *cov;
	gsl_multifit_linear_workspace *mw;
	double deltak,lastk;
	int i,j,countkeep;

	nbreak = 9;
	gsl_vector *mybreaks = gsl_vector_alloc(nbreak);
	gsl_vector_set(mybreaks,0,(0.001*0.7));
	gsl_vector_set(mybreaks,1,(0.025*0.7));
	gsl_vector_set(mybreaks,2,(0.075*0.7));
	gsl_vector_set(mybreaks,3,(0.125*0.7));
	gsl_vector_set(mybreaks,4,(0.175*0.7));
	gsl_vector_set(mybreaks,5,(0.225*0.7));
	gsl_vector_set(mybreaks,6,(0.275*0.7));
	gsl_vector_set(mybreaks,7,(0.325*0.7));
	gsl_vector_set(mybreaks,8,(0.375*0.7));

	countkeep = 0;
	for(i=0;i<(*npts);i++)  {
		if((kvals[i]) >= gsl_vector_get(mybreaks,0) && (kvals[i]) <= gsl_vector_get(mybreaks,nbreak-1)) {
			countkeep += 1;
			}
		}
	n = countkeep;
	ncoeffs = nbreak + 2;

	/* allocate a cubic bspline workspace (k = 4) */
	bw = gsl_bspline_alloc(4, nbreak);
	B = gsl_vector_alloc(ncoeffs);     
	x = gsl_vector_alloc(n);
	y = gsl_vector_alloc(n);
	X = gsl_matrix_alloc(n, ncoeffs);
	c = gsl_vector_alloc(ncoeffs);
	w = gsl_vector_alloc(n);
	cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
	mw = gsl_multifit_linear_alloc(n, ncoeffs);
	i=0;
	for(j=0;j<(*npts);j++)  {
		if((kvals[j]) >= gsl_vector_get(mybreaks,0) && (kvals[j]) <= gsl_vector_get(mybreaks,nbreak-1)) {
			gsl_vector_set(x,i,(kvals[j]));
			gsl_vector_set(y,i,exp(lnpklinear[j])*pow(kvals[j],1.5));
			if(j>0)  {
				deltak = kvals[j] - kvals[j-1];
				}
			else {
				deltak = kvals[0];
				if(kvals[1] - kvals[0] < deltak)  {
					deltak = kvals[1]-kvals[0];
					}
				}
			gsl_vector_set(w,i,deltak);
			i+=1;
			}
		}
	gsl_bspline_knots(mybreaks,bw);
	for(i=0;i<n;i++)  {
		double xi = gsl_vector_get(x,i);
		gsl_bspline_eval(xi,B,bw);
		for(j=0;j<ncoeffs;j++)  {
			double Bj = gsl_vector_get(B,j);
			gsl_matrix_set(X,i,j,Bj);
			}
		}
	//do fit
	double yi,yierr,chisq;
	gsl_multifit_wlinear(X,w,y,c,cov,&chisq,mw);
	i = 0;
	for(j=0;j<(*npts);j++)  {
		if((kvals[j]) >= gsl_vector_get(mybreaks,0) && (kvals[j]) <= gsl_vector_get(mybreaks,nbreak-1)) {
			gsl_bspline_eval(gsl_vector_get(x,i),B,bw);
			gsl_multifit_linear_est(B,c,cov,&yi,&yierr);
			lnpksmooth[j] = log(yi*pow(kvals[j],-1.5));
			i += 1;
			}
		else {
			lnpksmooth[j] = lnpklinear[j];
			}
		//spline is wacky at small k -- suppress difference at k < 0.01
		if(kvals[j] < kmaxsuppress)  {
			lnpksmooth[j] = lnpklinear[j];
			}
		}
	assert(i==n);
	gsl_bspline_free(bw);
	gsl_vector_free(B);
	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_vector_free(mybreaks);
	gsl_matrix_free(X);
	gsl_vector_free(c);
	gsl_vector_free(w);
	gsl_matrix_free(cov);
	gsl_multifit_linear_free(mw);
	}
/*

FILE *open_file_read(char *filename) {
	FILE *ifp;

	if(!(ifp=fopen(filename,"r")))  {
		fprintf(stderr,"Can't open %s\n",filename);
		return NULL;
		}
	else {
		return ifp;
		}
	}	

int get_file_length(FILE *ifp, int *headercount, int *linecount) {
	int i;
	char line[256];
	
	*headercount = 0;
	*linecount = 0;
			
    while(fgets(line,256,ifp))  {
        if(!feof(ifp))  {
            (*linecount)++;
            }
		else {
			break;
			}
		if(strstr(line,"#"))  {
			(*headercount)++;
			}
        }
	(void) rewind(ifp);
	}



int main()  {
	double *lnpksmooth, *lnpklinear, *kvals;
	FILE *ifp;
	int headercount, linecount,i,j;
	double kval,pnw,plin;
	ifp = open_file_read("/Users/breid/cosmomcmnu/test1.cut");
	get_file_length(ifp,&headercount,&linecount);
	lnpksmooth = (double *) malloc(sizeof(double)*linecount);
	lnpklinear = (double *) malloc(sizeof(double)*linecount);
	kvals = (double *) malloc(sizeof(double)*linecount);
	for(i=0;i<linecount;i++)  {
		fscanf(ifp,"%lf %lf %lf\n",&kval,&pnw,&plin);
		lnpklinear[i] = plin;
		kvals[i] = kval*0.701;
		//kvals[i] = kval;
		}
	dopksmoothbspline_(kvals,lnpklinear,lnpksmooth,&linecount);
	for(i=0;i<linecount;i++)  {
		printf("%f %f %f\n",kvals[i],lnpklinear[i],lnpksmooth[i]);
		}
	return 0;
	}
*/
