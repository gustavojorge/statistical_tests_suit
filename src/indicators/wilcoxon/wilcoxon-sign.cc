
/* wilcoxon-sign.cc  (C) Joshua Knowles, 2005

Implements a nonparametric test for differences between two paired (or matched) samples,
as described in W.J.Conover (1999) "Practical Nonparametric Statisticsn (3rd Edition)", Wiley.
If the sample size n>=50 then the normal approximation is used. If n<=50, the exact quantiles
from a table are used. In this case, only the lowest from a set of critical p-values that
exceeds the ranksum is returned: no interpolation of values in the table is used.

This version accepts multiple (more than two) samples. In that case, the test is
carried out for each pair and a warning, which advises the p-values are not
accurate because the samples are no longer independent random samples, is issued.



Compile and link with the attached Makefile:
   make wilcoxon
    g++ wilcoxon-sign.cc dcdflib.o -o wilcoxon-sign \\ by Felipe

Run:
   ./wilcoxon <indicator_file> <param_file> <output_file>
   
   where:

   <indicator_file> is the name of a file containing a single column of
     indicator values. Blank lines in the file divide the sample
     populations. All the sample populations must be the same size, and ordered
     so that the entries are matched pairs;
   <param_file> is the name of a file that, at present, may be empty;
   <output_file> is a filename to write to.

Output:

   The output is the one-tailed p-value that the null hypothesis is true for each pair.
   If the normal approximation has been used then this is indicated in brackets.

   A warning message is output if there are multiple (more than two) samples.

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

typedef struct pairs
{
  double diff; // signed difference
  double rank; // always positive
  double sigrank; // the signed rank
}P;

P *p;
D *d;

int N; // the total number of values in the input
int *Nsamp; // the number of values in each sample population
int ndist; // the number of sample populations

double W[51][9];
double Wsig[9];

FILE *fp;

void init_W();
double myZ(double x);
double myt(double t, double df);
double mychi(double x, double df);
double myabs(double v);
int comparediff(const void *, const void *);
int assign_ranks(P *p, int N);
double sum_of_ranks(D *d, int index, int N);
double sum_squared_ranks(D *d, int N);
double Tvalue(D *d, int N, int ndist, int *Nsamp);
double S_squared(D *d, int N, int ndist);
void  check_file(FILE  *fp, int  *no_runsp, int  *max_pointsp, int *Nsamp);
void  read_file(FILE  *fp, int  *no_pointsp, D *d);
double pairwise(int a, int b, D *d, int N, int *Nsamp, double T);

int main(int argc, char **argv)
{
    int a, b, i, j;
    
  init_W();
  
  if(argc!=4)
    {
      fprintf(stderr,"./wilcoxon <indicator_file> <param_file> <output_file>\n");
      exit(1);
    }

  if((fp=fopen(argv[2],"rb")))
    {
      // Parameter file not used yet.
      
      fclose(fp);
    }
  

  Nsamp = (int *)malloc(MAX_DISTS*sizeof(int));
  for( j=0;j<MAX_DISTS;j++) Nsamp[j] = 0;

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
      p = (P *)malloc(Nsamp[0] *sizeof(P));
      read_file(fp, &N, d);
      fclose(fp);
    }
  else
    {
      fprintf(stderr,"Couldn't open %s for reading\n", argv[1]);
      exit(1);
    }

  if((fp=fopen(argv[3],"w")))
    fclose(fp);
  else
    {
      fprintf(stderr,"Couldn't open %s for writing.\n", argv[3]);
      exit(1);
    }

  for( a=0;a<ndist;a++)
    {
      for( b=0;b<ndist;b++)
	{
	  if(a==b)
	    continue;
	  // check for pairs that do not differ (ties) and remove them
	  int ties=0;
	  for( i=0;i<Nsamp[0]; i++)
	    {
	      p[i].diff = d[Nsamp[0]*b+i].value-d[Nsamp[0]*a+i].value;
	      if(VERBOSE)
		fprintf(stdout, "%g %g ", d[Nsamp[0]*a+i].value,d[Nsamp[0]*b+i].value);
	      if(VERBOSE)
		fprintf(stdout, "%g\n", p[i].diff);
	      if(p[i].diff==0)
		{
		  p[i].diff=10e99;
		  ties++;
		}
	    }
	  
  
	  qsort(p, Nsamp[0], sizeof(P), comparediff);
	  int n=Nsamp[0]-ties;

	  assign_ranks(p,n);

	  for(i=0;i<n;i++)
	    {
	      if(p[i].diff<0)
		p[i].sigrank=-p[i].rank;
	      else
		p[i].sigrank=p[i].rank;
	    }
	      
	  if(VERBOSE)
	    {
	      fprintf(stdout, "__diff__\tabs_diff\t__rank__\tsign_rnk:\n");
	      for( i=0;i<n;i++)
		{
		  fprintf(stdout, "%8g\t%8g\t", p[i].diff, myabs(p[i].diff));
		  fprintf(stdout, "%8g\t%8g\n", p[i].rank, p[i].sigrank);
		}
	    }
	  if(n<4)
	    {
	      fprintf(stderr,"Need at least 4 values in a sample to perform signed-rank test.");
	      exit(1);
	    }
	  
	      
	  if((n>50)||(double(ties)/double(Nsamp[0])>0.5))
	    {
	    
	      if(VERBOSE)
		fprintf(stdout,"Using the standard normal approximation because n>50 or there are many ties\n");
	  
	      double sum_of_ranks=0.0; 
	      double sum_of_sq_ranks=0.0;
	      
	      for( i=0;i<n;i++)
		{
		  sum_of_ranks+=p[i].sigrank;
		  sum_of_sq_ranks+=pow(p[i].sigrank,2.0);
		}
	      
	      if(VERBOSE)
		fprintf(stdout, "sum of signed ranks = %g\n", sum_of_ranks);
	      if(VERBOSE)
		fprintf(stdout, "sum of squared ranks = %g\n", sum_of_sq_ranks);
	      
	      double lowerp = 1.0-myZ((sum_of_ranks+1.0)/sqrt(sum_of_sq_ranks)); // Equation 7, page 354 of Conover, 1999.
	      
	      double upperp = myZ((sum_of_ranks-1.0)/sqrt(sum_of_sq_ranks)); // Equation 8, page 354 of Conover, 1999.
	      if(VERBOSE)
		fprintf(stdout, "upper p = %g\n", upperp);
	      if(VERBOSE)
		fprintf(stdout, "lower p = %g\n", lowerp);
	      
	      
	      double pvalue;
	      if(lowerp<upperp)
		pvalue=2*lowerp;
	      else
		pvalue=2*upperp;
	      
	      pvalue=upperp; // for a 1-tailed test
	      
	      if(VERBOSE)
		fprintf(stdout, "The one-tailed p-value for accepting the null hypothesis that the expected value of the difference is zero is p=%g\n", pvalue);
	      if((fp=fopen(argv[3],"a")))
		{
		  fprintf(fp, "%d better than %d with a p-value of %g\n", b+1, a+1, pvalue); // Normal approximation
		  fclose(fp);     
		}
	      else
		{
		  fprintf(stderr,"Couldn't open output file for writing.\n");
		  exit(1);
		}
	    }
	  else
	    {
	      double upperp, lowerp, pvalue;
	      if(VERBOSE)
		fprintf(stdout,"Using exact critical values of the W distribution from a lookup table\n");
	      double Tplus=0.0;
	      for( i=0;i<n;i++)
		{	       
		  if(p[i].diff>0)
		    Tplus += p[i].sigrank;  // Equation 3, page 353 of Conover (1999)
		}
	          
	      if(VERBOSE)
		fprintf(stdout, "Tplus =%g\n",Tplus);


	      int j=0;
	      while(W[n][j]<Tplus && j <8)
		{
		  j++;

		}
	      upperp = Wsig[j];
	      j=0;
	      while( (n*(n-1))/2 - W[n][j] > Tplus && j < 8)
		j++;
	      lowerp = Wsig[j];
	      
	      
	      
	      if(VERBOSE)
		fprintf(stdout, "upper p = %g\n", upperp);
	      if(VERBOSE)
		fprintf(stdout, "lower p = %g\n", lowerp);
	      
	      
	      if(lowerp<upperp)
		pvalue=2*lowerp;
	      else
		pvalue=2*upperp;
	      
	      pvalue=upperp; // for a 1-tailed test
	      
	      if(VERBOSE)
		fprintf(stdout, "The one-tailed p-value for accepting the null hypothesis that the expected value of the difference is zero is p=%g\n", pvalue);
	      if((fp=fopen(argv[3],"a")))
		{
		  fprintf(fp, "%d better than %d with a p-value of %g\n", b+1, a+1, pvalue); // Critical value of exact W table
		  fclose(fp);     
		}
	      else
		{
		  fprintf(stderr,"Couldn't open output file for writing.\n");
		  exit(1);
		}
	    }
	  
	}
    }
  
  if(ndist>2)
    fprintf(stderr, "Warning: the p-values for accepting the null hypothesis that the expected differences are zero, are not correct because multiple tests have been carried out using the same sample. Therefore, these values should only be used in preliminary (explorative) tests, and do not indicate true probabilities. Consider collecting new, independent random samples for each statistical test to be performed.\n");
  return(0);
  
  
}


double myZ(double x)
{
  int which=1;
  double p;
  double q;
  double mean=0.0;
  double sd=1.0;
  int status;
  double bound;
  cdfnor(&which, &p, &q, &x, &mean, &sd, &status, &bound);
  return(p);
}


double pairwise(int a, int b, D *d, int N, int *Nsamp, double T)
{
  double value;

  value =myabs( (sum_of_ranks(d, a, N)/double(Nsamp[a])) - (sum_of_ranks(d, b, N)/double(Nsamp[b])));

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
  double bound;
  double p;
  double q;
  int status;
  int which=1;

  cdfchi( &which, &p, &q, &x, &df, &status, &bound );
  return(q);
  
}

double myt(double t, double df)
{
  double bound;
  double p;
  double q;
  int status;
  int which=1;

  cdft ( &which, &p, &q, &t, &df, &status, &bound );
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

int assign_ranks(P *p, int N)
{
  int i,j;
  int crank=1;
  int totalrank;
  int count;
  int total_ties=0;
  P *ahead;

  i=0;
  while(i<N)
    {
      ahead = &(p[i+1]);
      if(myabs(ahead->diff) == myabs(p[i].diff))
	{	
	  totalrank=crank;
	  count=0;
	  do
	    {
	      ahead++;
	      count++;
	      totalrank+=(crank+count);
	      //  printf("i+count=%d\n", i+count);
	    }while((myabs(ahead->diff) == myabs(p[i].diff))&&(i+count<N-1));
	  // set all the ranks to the average value
	  for(j=0;j<=count;j++)
	    p[i+j].rank = (double(totalrank)/double(count+1));
	  i+=count+1;
	  crank+=count+1;
	  total_ties+=count;
	}
      else
	{
	  p[i].rank=crank;
	  crank++;
	  i++;
	}
    }  
  return(total_ties);

}


int comparediff(const void *i, const void *j)
{
  double x;
  x = myabs((*(P *)i).diff) - myabs((*(P *)j).diff);
  
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
  int j, new_run;
  
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
  for(j=1;j<*no_runsp;j++)
  if(Nsamp[j]!=Nsamp[0])
    {
      fprintf(stderr,"Two samples of indicator values are not of the same size. This program computes statistics for paired samples only. Exiting.\n");
      exit(1);
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

void init_W()
{
  // W values and significance levels from Table A12 of Conover (1999)
  Wsig[0]=0.005;
  Wsig[1]=0.01;
  Wsig[2]=0.025;
  Wsig[3]=0.05;
  Wsig[4]=0.1;
  Wsig[5]=0.2;
  Wsig[6]=0.3;
  Wsig[7]=0.4;
  Wsig[8]=0.5;

  W[4][0]=0;
  W[4][1]=0;
  W[4][2]=0;
  W[4][3]=0;
  W[4][4]=1;
  W[4][5]=3;
  W[4][6]=3;
  W[4][7]=4;
  W[4][8]=5;

  W[5][0]=0;
  W[5][1]=0;
  W[5][2]=0;
  W[5][3]=1;
  W[5][4]=3;
  W[5][5]=4;
  W[5][6]=5;
  W[5][7]=6;
  W[5][8]=7.5;
  
  W[6][0]=0;
  W[6][1]=0;
  W[6][2]=1;
  W[6][3]=3;
  W[6][4]=4;
  W[6][5]=6;
  W[6][6]=8;
  W[6][7]=9;
  W[6][8]=10.5;

  W[7][0]=0;
  W[7][1]=1;
  W[7][2]=3;
  W[7][3]=4;
  W[7][4]=6;
  W[7][5]=9;
  W[7][6]=11;
  W[7][7]=12;
  W[7][8]=14;


  W[8][0]=1;
  W[8][1]=2;
  W[8][2]=4;
  W[8][3]=6;
  W[8][4]=9;
  W[8][5]=12;
  W[8][6]=14;
  W[8][7]=16;
  W[8][8]=18;

W[9][0]=2;
W[9][1]=4;
W[9][2]=6;
W[9][3]=9;
W[9][4]=11;
W[9][5]=15;
W[9][6]=18;
W[9][7]=20;
W[9][8]=22.5;
W[10][0]=4;
W[10][1]=6;
W[10][2]=9;
W[10][3]=11;
W[10][4]=15;
W[10][5]=19;
W[10][6]=22;
W[10][7]=25;
W[10][8]=27.5;
W[11][0]=6;
W[11][1]=8;
W[11][2]=11;
W[11][3]=14;
W[11][4]=18;
W[11][5]=23;
W[11][6]=27;
W[11][7]=30;
W[11][8]=33;
W[12][0]=8;
W[12][1]=10;
W[12][2]=14;
W[12][3]=18;
W[12][4]=22;
W[12][5]=28;
W[12][6]=32;
W[12][7]=36;
W[12][8]=39;
W[13][0]=10;
W[13][1]=13;
W[13][2]=18;
W[13][3]=22;
W[13][4]=27;
W[13][5]=33;
W[13][6]=38;
W[13][7]=42;
W[13][8]=45.5;
W[14][0]=13;
W[14][1]=16;
W[14][2]=22;
W[14][3]=26;
W[14][4]=32;
W[14][5]=39;
W[14][6]=44;
W[14][7]=48;
W[14][8]=52.5;
W[15][0]=16;
W[15][1]=20;
W[15][2]=26;
W[15][3]=31;
W[15][4]=37;
W[15][5]=45;
W[15][6]=51;
W[15][7]=55;
W[15][8]=60;
W[16][0]=20;
W[16][1]=24;
W[16][2]=30;
W[16][3]=36;
W[16][4]=43;
W[16][5]=51;
W[16][6]=58;
W[16][7]=63;
W[16][8]=68;
W[17][0]=24;
W[17][1]=28;
W[17][2]=35;
W[17][3]=42;
W[17][4]=49;
W[17][5]=58;
W[17][6]=65;
W[17][7]=71;
W[17][8]=76.5;
W[18][0]=28;
W[18][1]=33;
W[18][2]=41;
W[18][3]=48;
W[18][4]=56;
W[18][5]=66;
W[18][6]=73;
W[18][7]=80;
W[18][8]=85.5;
W[19][0]=33;
W[19][1]=38;
W[19][2]=47;
W[19][3]=54;
W[19][4]=63;
W[19][5]=74;
W[19][6]=82;
W[19][7]=89;
W[19][8]=95;
W[20][0]=38;
W[20][1]=44;
W[20][2]=53;
W[20][3]=61;
W[20][4]=70;
W[20][5]=83;
W[20][6]=91;
W[20][7]=98;
W[20][8]=105;
W[21][0]=44;
W[21][1]=50;
W[21][2]=59;
W[21][3]=68;
W[21][4]=78;
W[21][5]=91;
W[21][6]=100;
W[21][7]=108;
W[21][8]=115.5;
W[22][0]=49;
W[22][1]=56;
W[22][2]=67;
W[22][3]=76;
W[22][4]=87;
W[22][5]=100;
W[22][6]=110;
W[22][7]=119;
W[22][8]=126.5;
W[23][0]=55;
W[23][1]=63;
W[23][2]=74;
W[23][3]=84;
W[23][4]=95;
W[23][5]=110;
W[23][6]=120;
W[23][7]=130;
W[23][8]=138;
W[24][0]=62;
W[24][1]=70;
W[24][2]=82;
W[24][3]=92;
W[24][4]=105;
W[24][5]=120;
W[24][6]=131;
W[24][7]=141;
W[24][8]=150;
W[25][0]=69;
W[25][1]=77;
W[25][2]=90;
W[25][3]=101;
W[25][4]=114;
W[25][5]=131;
W[25][6]=143;
W[25][7]=153;
W[25][8]=162.5;
W[26][0]=76;
W[26][1]=85;
W[26][2]=99;
W[26][3]=111;
W[26][4]=125;
W[26][5]=142;
W[26][6]=155;
W[26][7]=165;
W[26][8]=175.5;
W[27][0]=84;
W[27][1]=94;
W[27][2]=108;
W[27][3]=120;
W[27][4]=135;
W[27][5]=154;
W[27][6]=167;
W[27][7]=178;
W[27][8]=189;
W[28][0]=92;
W[28][1]=102;
W[28][2]=117;
W[28][3]=131;
W[28][4]=146;
W[28][5]=166;
W[28][6]=180;
W[28][7]=192;
W[28][8]=203;
W[29][0]=101;
W[29][1]=111;
W[29][2]=127;
W[29][3]=141;
W[29][4]=158;
W[29][5]=178;
W[29][6]=193;
W[29][7]=206;
W[29][8]=217.5;
W[30][0]=110;
W[30][1]=121;
W[30][2]=138;
W[30][3]=152;
W[30][4]=170;
W[30][5]=191;
W[30][6]=207;
W[30][7]=220;
W[30][8]=232.5;
W[31][0]=119;
W[31][1]=131;
W[31][2]=148;
W[31][3]=164;
W[31][4]=182;
W[31][5]=205;
W[31][6]=221;
W[31][7]=235;
W[31][8]=248;
W[32][0]=129;
W[32][1]=141;
W[32][2]=160;
W[32][3]=176;
W[32][4]=195;
W[32][5]=219;
W[32][6]=236;
W[32][7]=250;
W[32][8]=264;
W[33][0]=139;
W[33][1]=152;
W[33][2]=171;
W[33][3]=188;
W[33][4]=208;
W[33][5]=233;
W[33][6]=251;
W[33][7]=266;
W[33][8]=280.5;
W[34][0]=149;
W[34][1]=163;
W[34][2]=183;
W[34][3]=201;
W[34][4]=222;
W[34][5]=248;
W[34][6]=266;
W[34][7]=282;
W[34][8]=297.5;
W[35][0]=160;
W[35][1]=175;
W[35][2]=196;
W[35][3]=214;
W[35][4]=236;
W[35][5]=263;
W[35][6]=283;
W[35][7]=299;
W[35][8]=315;
W[36][0]=172;
W[36][1]=187;
W[36][2]=209;
W[36][3]=228;
W[36][4]=251;
W[36][5]=279;
W[36][6]=299;
W[36][7]=317;
W[36][8]=333;
W[37][0]=184;
W[37][1]=199;
W[37][2]=222;
W[37][3]=242;
W[37][4]=266;
W[37][5]=296;
W[37][6]=316;
W[37][7]=335;
W[37][8]=351.5;
W[38][0]=196;
W[38][1]=212;
W[38][2]=236;
W[38][3]=257;
W[38][4]=282;
W[38][5]=312;
W[38][6]=334;
W[38][7]=353;
W[38][8]=370.5;
W[39][0]=208;
W[39][1]=225;
W[39][2]=250;
W[39][3]=272;
W[39][4]=298;
W[39][5]=329;
W[39][6]=352;
W[39][7]=372;
W[39][8]=390;
W[40][0]=221;
W[40][1]=239;
W[40][2]=265;
W[40][3]=287;
W[40][4]=314;
W[40][5]=347;
W[40][6]=371;
W[40][7]=391;
W[40][8]=410;
W[41][0]=235;
W[41][1]=253;
W[41][2]=280;
W[41][3]=303;
W[41][4]=331;
W[41][5]=365;
W[41][6]=390;
W[41][7]=411;
W[41][8]=430.5;
W[42][0]=248;
W[42][1]=267;
W[42][2]=295;
W[42][3]=320;
W[42][4]=349;
W[42][5]=384;
W[42][6]=409;
W[42][7]=431;
W[42][8]=451.5;
W[43][0]=263;
W[43][1]=282;
W[43][2]=311;
W[43][3]=337;
W[43][4]=366;
W[43][5]=403;
W[43][6]=429;
W[43][7]=452;
W[43][8]=473;
W[44][0]=227;
W[44][1]=297;
W[44][2]=328;
W[44][3]=354;
W[44][4]=385;
W[44][5]=422;
W[44][6]=450;
W[44][7]=473;
W[44][8]=495;
W[45][0]=292;
W[45][1]=313;
W[45][2]=344;
W[45][3]=372;
W[45][4]=403;
W[45][5]=442;
W[45][6]=471;
W[45][7]=495;
W[45][8]=517.5;
W[46][0]=308;
W[46][1]=329;
W[46][2]=362;
W[46][3]=390;
W[46][4]=423;
W[46][5]=463;
W[46][6]=492;
W[46][7]=517;
W[46][8]=540.5;
W[47][0]=324;
W[47][1]=346;
W[47][2]=379;
W[47][3]=408;
W[47][4]=442;
W[47][5]=484;
W[47][6]=514;
W[47][7]=540;
W[47][8]=564;
W[48][0]=340;
W[48][1]=363;
W[48][2]=397;
W[48][3]=428;
W[48][4]=463;
W[48][5]=505;
W[48][6]=536;
W[48][7]=563;
W[48][8]=588;
W[49][0]=357;
W[49][1]=381;
W[49][2]=416;
W[49][3]=447;
W[49][4]=483;
W[49][5]=527;
W[49][6]=559;
W[49][7]=587;
W[49][8]=612.5;
W[50][0]=374;
W[50][1]=398;
W[50][2]=435;
W[50][3]=467;
W[50][4]=504;
W[50][5]=550;
W[50][6]=583;
W[50][7]=611;
W[50][8]=637.5;

}
