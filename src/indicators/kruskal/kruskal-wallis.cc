
/* kruskal-wallis.cc  (C) Joshua Knowles, 2005

Implements a nonparametric test for differences between multiple independent samples,
as described in W.J.Conover (1999) "Practical Nonparametric Statistics (3rd Edition)", Wiley.

Compile and link with the attached Makefile:
   make kruskal

Run:
   ./kruskal <indicator_file> <param_file> <output_file>
   
   where:

   <indicator_file> is the name of a file containing a single column of
     indicator values. Blank lines in the file divide the separate sample
     populations;
   <param_file> is the name of a file with the following one-line format

       alpha 0.05

     where the alpha value specifies the significance level and should be in the range (0,0.1];
   <output_file> is a filename to write to.

Output:

   If and only if a first test for significance of any differences between the
   samples is passed, at the given alpha value, then the output will be the one-taileed p-values
   for rejecting a null hypothesis of no significant difference between two samples,
   for each pair-wise combination.

   In the case that the first test fails, the output is simply, "H0".

   With VERBOSE set to true, some output to stdout is also given.

*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "dcdflib.h"

// using namespace std;

#define RN rand()/(RAND_MAX+1.0)
#define MAX_STR_LENGTH 100
#define MAX_LINE_LENGTH 100
#define MAX_DISTS 30 
#define VERBOSE true

typedef struct data
{
  double value;
  int label;
  double rank;
}D;

D *d;
int N; // the total number of samples from all distributions
int *Nsamp; // the number of samples in each distribution
int ndist; // the number of distributions

FILE *fp;

double myt(double t, double df);
double mychi(double x, double df);
double myabs(double v);
int compare(const void *, const void *);
int assign_ranks(D *d, int N);
double sum_of_ranks(D *d, int index, int N);
double sum_squared_ranks(D *d, int N);
double Tvalue(D *d, int N, int ndist, int *Nsamp);
double S_squared(D *d, int N, int ndist);
void  check_file(FILE  *fp, int  *no_runsp, int  *max_pointsp, int *Nsamp);
void  read_file(FILE  *fp, int  *no_pointsp, D *d);
double pairwise(int a, int b, D *d, int N, int *Nsamp, double T);

int main(int argc, char **argv)
{
  int i, j;
  double alpha;
  
  if(argc!=4)
    {
      fprintf(stderr,"./kruskal <indicator_file> <param_file> <output_file>\n");
      exit(1);
    }

  if((fp=fopen(argv[2],"rb")))
    {
      if(fscanf(fp, "%*s %lg\n", &alpha)==EOF)	
	fprintf(stderr, "Error occurred in parameter file.\n"), exit(1);

      fclose(fp);
    }
  else
    {
      fprintf(stderr, "Couldn't open param file %s for reading.\n", argv[2]);
      exit(1);
    }

  if((alpha>0.1)||(alpha<=0))
    {
      fprintf(stderr, "The significance, alpha, must be in the range (0,0.1]\n");
      exit(1);
    }

  Nsamp = (int *)malloc(MAX_DISTS*sizeof(int));
  for( j=0;j<MAX_DISTS;j++) Nsamp[j] = 0;

  if((fp=fopen(argv[1],"rb")))
    {
      check_file(fp, &ndist, &N, Nsamp);
      if(VERBOSE)
	fprintf(stdout,"Number of sample populations = %d. Total number of values in the input = %d\n", ndist, N);
      rewind(fp);
      d = (D *)malloc(N *sizeof(D));
      read_file(fp, &N, d);
      fclose(fp);
    }
  else
    {
      fprintf(stderr,"Couldn't open %s for reading\n", argv[1]);
      exit(1);
    }
  
  for( i=0;i<ndist;i++)
    if(Nsamp[i]<20)
      {
	fprintf(stderr, "Warning: Sample population %d is of size %d. This software is not using a correction for small samples. Your samples should contain at least 20 values: the p-values returned for tests with this sample will be approximate.\n", i+1, Nsamp[i]);
      }

  qsort(d, N, sizeof(D), compare);
  int t = assign_ranks(d, N);

  if(VERBOSE)
    {
      for( i=0;i<N;i++)
	{
	  fprintf(stdout, "%g %d %.2g\n", d[i].value, d[i].label, d[i].rank);
	}
      fprintf(stdout,"Total number of ties =%d\n", t);
    }

  if(VERBOSE)
    for( j=0;j<ndist;j++)    
      fprintf(stdout, "Number of samples = %d; sum = %g\n", Nsamp[j], sum_of_ranks(d, j, N));

  double T;
  T=Tvalue(d, N, ndist, Nsamp);

  if(VERBOSE)
    fprintf(stdout, "Corrected T value =%g\n", T );

  

  
  double allsame=mychi(T, ndist-1);
  
  if(VERBOSE)
    fprintf(stdout, "p-value to accept the null hypothesis that all distribution functions are identical = %.9g\n", allsame);
  
  
  if((fp=fopen(argv[3],"w")))
    {
      if(allsame<=alpha)
	{
	  // fprintf(fp, "Overall p-value = %g. Null hypothesis rejected (alpha %g)\n", allsame, alpha);
	  for( i=0;i<ndist;i++)
	    for( j=0;j<ndist;j++)
	      {		
		if(i==j)
		  continue;
		fprintf(fp, "%d better than %d with a p-value of %g\n", j+1,i+1, myt(pairwise(i, j, d, N, Nsamp,T),N-ndist));
	      }

	}
      else
	fprintf(fp, "H0");
      fclose(fp);
    }
  else
    {
      fprintf(stderr, "Couldn't open output file for writing\n");
      exit(1);
    }
  
  return(0);

}

double pairwise(int a, int b, D *d, int N, int *Nsamp, double T)
{
  // Implements Equation 6, page 290 of Conover (1999).
  double value;

  value = (sum_of_ranks(d, a, N)/double(Nsamp[a])) - (sum_of_ranks(d, b, N)/double(Nsamp[b]));

  double denom;
  denom = sqrt(S_squared(d, N, ndist)*(N-1.0-T)/(N-ndist)) * sqrt(1.0/Nsamp[a]+1.0/Nsamp[b]);

  return(value/denom);
  
}

double myabs(double v)
{
  if(v>=0)
    return v;
  else
    return -v;
}

double mychi(double x, double df)
{
  // returns the probability that the chi-square distribtion of degree df will have a value >= x
  double bound;
  double p;
  double q;
  int status;
  int which=1;

  cdfchi( &which, &p, &q, &x, &df, &status, &bound );  // library function for the cdf of the chi-square dist.
  return(q);
  
}

double myt(double t, double df)
{
  // returns the probability that the t distribtion of degree df will have a value >= t
  double bound;
  double p;
  double q;
  int status;
  int which=1;

  cdft ( &which, &p, &q, &t, &df, &status, &bound );  // library function for the cdf of the t dist.
  return(q);

}


double Tvalue(D *d, int N, int ndist, int *Nsamp)
{
  // Equation 3, page 289 Conover (1999)

  double T;
  double *R;
  double S2;
  int i;

  R = (double *)malloc(ndist*sizeof(double));  

  for(i=0;i<ndist;i++)
    R[i] = sum_of_ranks(d,i,N);
  
  S2 = S_squared(d, N, ndist);

  double sum=0.0;
  for(i=0;i<ndist;i++)
    {
      sum += (R[i]*R[i])/double(Nsamp[i]);
    }

  T = (1.0/S2)*(sum - ((N*(N+1.0)*(N+1.0))/4.0));
  
  return(T);
  
}


double S_squared(D *d, int N, int ndist)
{ 
  // Equation 4, page 289 of Conover (1999)
  return ( (1.0/(N-1.0))*(sum_squared_ranks(d, N) - ((N*(N+1.0)*(N+1.0))/4.0)) );
}

double sum_squared_ranks(D *d, int N)
{
  double sum=0.0;
  int i;
  
  for(i=0;i<N;i++)
    sum += pow(d[i].rank,2.0);
  return(sum);
}


double sum_of_ranks(D *d, int index, int N)
{
  double sum=0.0;
  int i;
  
  for(i=0;i<N;i++)
    {
      if(d[i].label == index)
	sum+=d[i].rank;
    }

  return(sum);
}

int assign_ranks(D *d, int N)
{
  // assign ranks to the N values, giving the same (averaged) rank to any tied values
  // NOTE: the N values in d must be in sorted order
  int i,j;
  int crank=1;
  int totalrank;
  int count;
  int total_ties=0;
  D *ahead;

  i=0;
  while(i<N)
    {
      ahead = &(d[i+1]);
      if(ahead->value == d[i].value)
	{	
	  totalrank=crank;
	  count=0;
	  do
	    {
	      ahead++;
	      count++;
	      totalrank+=(crank+count);
	      //  printf("i+count=%d\n", i+count);
	    }while((ahead->value == d[i].value)&&(i+count<N-1));
	  // set all the ranks to the average value
	  for(j=0;j<=count;j++)
	    d[i+j].rank = (double(totalrank)/double(count+1));
	  i+=count+1;
	  crank+=count+1;
	  total_ties+=count;
	}
      else
	{
	  d[i].rank=crank;
	  crank++;
	  i++;
	}
    }  
  return(total_ties);

}


int compare(const void *i, const void *j)
{
  double x;
  x = (*(D *)i).value - (*(D *)j).value;
  
  if(x<0)
    return(-1);

  else if (x>0)
    return(1);
  
  else
    return(0);

}


void  check_file(FILE  *fp, int  *no_runsp, int  *totalp, int *Nsamp)
{
  char  line[MAX_STR_LENGTH];
  double  number;
  int new_run;
  
  *totalp = 0;
  *no_runsp = 0;
  new_run = 1;

  while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) {
      if (sscanf(line, "%lf", &number) != 1) {
			new_run = 1;
	  } else {
			if (new_run == 1) (*no_runsp)++;
			new_run = 0;
			(*totalp)++;
			if(*no_runsp<=MAX_DISTS) {
				(Nsamp[*no_runsp-1])++;
			} else {
				fprintf(stderr,"Please edit MAX_DISTS. Number of sample distributions exceeded the current setting.\n");
				exit(1);
			}
	  }
	
  }
  
}

void  read_file(FILE  *fp, int  *no_pointsp, D *d)
{
  char  line[MAX_STR_LENGTH];
  double  number;
  int clabel=0;
  
  *no_pointsp = 0;
  while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) 
    {
      if (sscanf(line, "%lf", &number) != 1)
	{
	  clabel++;
	}
      else 
	{
	  d[*no_pointsp].value = number;
	  d[*no_pointsp].label = clabel;
	  (*no_pointsp)++;
	}
    } 
}
