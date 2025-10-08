/* mann-whit.cc  (C) Joshua Knowles, 2005

Implements a nonparametric test for differences between precisely two independent samples,
as described in W.J.Conover (1999) "Practical Nonparametric Statistics (3rd Edition)", Wiley.
The version here uses T1 - the equation to correct for ties in the inputs. However, it does
not compute corrections for small samples. Instead, a warning is issued if either
sample is smaller than 20.

This version accepts multiple (more than two) samples. In that case, the test is
carried out for each pair and a warning, which advises the p-values are not
accurate because the samples are no longer independent random samples, is issued.



Compile and link with the attached Makefile:
   make mann-whit
   OU
   g++ mann-whit.cc dcdflib.o -o mann-whit // by Felipe

Run:
   ./mann-whit <indicator_file> <param_file> <output_file>
   
   where:

   <indicator_file> is the name of a file containing a single column of
     indicator values. Blank lines in the file divide the separate sample
     populations;
   <param_file> is the name of a file that, at present, may be empty;
   <output_file> is a filename to write to.

Output:

   The output is the p-value of the two-tailed test for rejecting 
     the null hypothesis of no difference between each pair of sample populations.

   A warning message is output if there are multiple (more than two) samples.
     
   With VERBOSE set to true, some output to stdout is also given.

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../utils/dcdflib/dcdflib.h"
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
D *pair;
int N; // the total number of values in the input
int *Nsamp; // the number of values in each sample population
int ndist; // the number of sample populations
FILE *fp;

double myabs(double v);
double myZ(double x);
int compare(const void *, const void *);
int assign_ranks(D *d, int N);
double sum_of_ranks(D *d, int index, int N);
double sum_squared_ranks(D *d, int N);
double corrected_Tvalue(D *d, int n, int m, int N, int idx);
void  check_file(FILE  *fp, int  *no_runsp, int  *max_pointsp, int *Nsamp);
void  read_file(FILE  *fp, int  *no_pointsp, D *d);

int main(int argc, char **argv)
{
    int i, j, k;
  int test=1;
  
  if(argc!=4)
    {
      fprintf(stderr,"./mann <indicator_file> <param_file> <output_file> \n");
      exit(1);
    }
  
  if((fp=fopen(argv[2],"rb")))
    {
      // Parameter file not used yet.

      fclose(fp);
    }


  Nsamp = (int *)malloc(MAX_DISTS*sizeof(int));
  for(j=0;j<MAX_DISTS;j++) Nsamp[j] = 0;

  if((fp=fopen(argv[1],"rb")))
    {
      check_file(fp, &ndist, &N, Nsamp);
      if(VERBOSE)
    fprintf(stdout, "Number of samples (populations) = %d. Total number of values in the input = %d\n", ndist, N);
      rewind(fp);
      if(VERBOSE)
    {
      fprintf(stdout,"Numbers of values in each sample = ");
      for( j=0;j<ndist;j++)
        fprintf(stdout,"%d ", Nsamp[j]);
      fprintf(stdout,"\n");
    }      
      d = (D *)malloc(N *sizeof(D));
      pair = (D *)malloc(N *sizeof(D));
      read_file(fp, &N, d);
      fclose(fp);
    }
  else
    {
      fprintf(stderr,"Couldn't open %s for reading\n", argv[1]);
      exit(1);
    }
  

  for(i=0;i<ndist;i++)
    if(Nsamp[i]<20)
      {
    fprintf(stderr, "Warning: Sample population %d is of size %d. This software is not using a correction for small samples. Your samples should contain at least 20 values: the p-values returned for tests with this sample will be approximate.\n", i+1, Nsamp[i]);
      }

  if((fp=fopen(argv[3],"w")))
    fclose(fp);
  else
    {
      fprintf(stderr,"Couldn't open %s for writing.\n", argv[3]);
      exit(1);
    }
  

  for(j=0;j<ndist;j++)
    {
      for(k=0;k<ndist;k++)
    {
      if(j==k)
        continue;
      if(VERBOSE)
        fprintf(stdout, "\n\n**** Test %d between %d and %d ****\n", test++, j+1, k+1);
      int start=0;
      for(i=0;i<j;i++)
        start+=Nsamp[i];
      for(i=0;i<Nsamp[j];i++)
        {
          pair[i].value = d[i+start].value;
          pair[i].label = j;
        }
      start=0;
      for(i=0;i<k;i++)
        start+=Nsamp[i];
      for(i=0;i<Nsamp[k];i++)
        {
          pair[i+Nsamp[j]].value = d[i+start].value;
          pair[i+Nsamp[j]].label = k;
        }

      qsort(pair, Nsamp[j]+Nsamp[k], sizeof(D), compare);
      int t = assign_ranks(pair, Nsamp[j]+Nsamp[k]);

      if(VERBOSE)
        {
          for(i=0;i<Nsamp[j]+Nsamp[k];i++)
        {
          fprintf(stdout, "%g %d %.2g\n", pair[i].value, pair[i].label, pair[i].rank);
        }
          fprintf(stdout,"Total number of ties =%d\n", t);
        }
      
      if(VERBOSE)
        {
          fprintf(stdout, "Number of samples = %d; sum = %g\n", Nsamp[j], sum_of_ranks(pair, j, Nsamp[j]));
          fprintf(stdout, "Number of samples = %d; sum = %g\n", Nsamp[k], sum_of_ranks(pair, k, Nsamp[k]));
        }
      
      double T;
      double p_value;
      T=corrected_Tvalue(pair, Nsamp[j], Nsamp[k], Nsamp[j]+Nsamp[k], j);
      p_value= (1.0-myZ(T));
      if(VERBOSE)
        fprintf(stdout, "Corrected T value =%g\n", T );
      if(VERBOSE)
        fprintf(stdout, "One-tailed p-value = %.9g\n", p_value);
      if((fp=fopen(argv[3],"a")))
        {
          fprintf(fp, "%d better than %d with a p-value of  %.9g\n", k+1, j+1, p_value);
          fclose(fp);
        }
      else
        {
          fprintf(stderr, "Couldn't open output file for writing\n");
          exit(1);
        }
    }
    }
  if(ndist>2)
    fprintf(stderr, "Warning: the p-values for accepting the null hypothesis that these are two samples from the same underlying distribution are not correct because multiple tests have been carried out using the same sample. Therefore, these values should only be used in preliminary (explorative) tests, and do not indicate true probabilities. Consider collecting new, independent random samples for each statistical test to be performed. Alternatively, use the Kruskal-Wallis test.\n");
 
  return(0);
}

double myZ(double x)
{ 
  // returns the p-value that a normal random variate with mean 0 and SD 1 would have a value <= x
  int which=1;
  double p;
  double q;
  double mean=0.0;
  double sd=1.0;
  int status;
  double bound;
  cdfnor(&which, &p, &q, &x, &mean, &sd, &status, &bound); // library function for computing the cdf of the normal Z distribution
  return(p);
}

double corrected_Tvalue(D *d, int n, int m, int N, int idx)
{
  // Equation 2, page 273 of Conover (1999)
  double T;
  double T1;
  double denom;
  double ssR; // sum of squared Ranks

  T = sum_of_ranks(d,idx,N);
  ssR = sum_squared_ranks(d, N);

  
  T1 = (T - ((n*(N+1))/2.0));

  denom = (double(n*m)/double(N*(N-1))*ssR) - double((m*n*(N+1)*(N+1))/double(4*(N-1)));
  if (T1 == 0.0) return 0.0;
  
  return(T1/sqrt(denom));
  
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
  
  double dtr;
  double dcount;

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
      dtr = (double)totalrank;
      dcount = (double)(count+1);
      //      printf("%g\n", dtr/dcount);
      for(j=0;j<=count;j++)
        {
          d[i+j].rank = dtr/dcount; //(double(totalrank)/double(count+1));
          //  printf("%g\n", d[i+j].rank);
        }
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

double myabs(double v)
{
  if(v>=0)
    return v;
  else
    return -v;
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
