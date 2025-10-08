/* bound.cc  (C) Joshua Knowles, January 2005

A program that reads in a file <datafile> consisting of a collection
of approximation sets, assumed to be from all the optimizers
under comparison, and not necessarily internally nondominated sets. 
A separate parameter file <param> gives 
the dimension of the data and whether each objective should
be minimized or maximized. It also gives a parameter phi.
The output is the upper and lower bound in each objective, the ideal 
point and the nadir point, and the latter when shifted
by phi times the range (in each objective).


For more information see sections below.

** Please contact me - Joshua Knowles - if you have any comments, suggestions
or questions regarding this program or metrics for measuring nondominated sets. 
My email address is j.knowles@manchester.ac.uk

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version. 

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details. 

   The GNU General Public License is available at:
      http://www.gnu.org/copyleft/gpl.html
   or by writing to: 
        The Free Software Foundation, Inc., 
        675 Mass Ave, Cambridge, MA 02139, USA.  

   

   COMPILE:
      g++ bound.cc -o bound -lm -Wall -pedantic

   RUN:
      ./bound [<param>] <datafile> <outfile>


    The format of the parameter file <param> is

      dim <number>
      obj <+|-> <+|-> ...
      phi <number>

      where dim specifies the number of objective dimensions
      obj specifies whether each objective should be minimized (-) or maximized (+)
      phi specifies the amount that the ideal and nadir point should be shifted. Any positive real is allowed.

    If the parameter file is omitted, default parameters are taken and
    the number of objectives is determined from the data file. 


   The format of the collection of approximation sets in <datafile> is
      <number> <number> ...
      <number> <number> ...
      [blank line]
      <number> <number> ...
      <number> <number> ...
      [blank line]
      .
      .
      <number> <number> ...

Note: blank lines have no effect - since all points
are considered together.

    The output of bound to <outfile> is:

     lower_bound <number> <number>... <number>
     upper_bound <number> <number>... <number>
     ideal <number> <number>... <number>
     nadir <number> <number>... <number>
     utopia_bound <number> <number>... <number>
     nadir_bound <number> <number> ...<number> 

*******************************************************************/



#include <ctime>
#include <iostream>       
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace std;

#define LARGE 10e50
#define MAX_LINE_LENGTH 1024
#define MAX_STR_LENGTH 200
#define VERBOSE 1
#define error(X,Y)  if (X) fprintf(stderr, Y "\n"), exit(1)

FILE *fp;

struct dnode{
  double *o;
  bool dominated;
  struct dnode *next;
  struct dnode *prev;
};


dnode **po;
dnode *ref;

int *minmax1;
double *best;
double *worst;
double phi;
int nobjs;


double myabs(double a);
void d_append (struct dnode **s, double *vec);
void d_display ( struct dnode *q );
void  check_file(FILE  *fp, int  *no_runsp, int  *max_pointsp);
int  determine_dim(FILE  *fp);
void  read_file(FILE  *fp, int  *no_pointsp, dnode **po);


void d_append (struct dnode **s, double *vec) 
/* adds a new node at the end of the doubly linked list */ 
{
  struct dnode *r, *q = *s ; 
  /* if the linked list is empty */ 
  
  if ( *s == NULL ) 
    { 
      
      /*create a new node */ 
      
      *s = (dnode *) malloc ( sizeof ( struct dnode ) ) ;
      ( *s ) -> prev = NULL ; 
      (*s)->o = (double *) malloc(nobjs * sizeof (double) );
      for(int j=0;j<nobjs;j++)
	( *s ) -> o[j] = vec[j]; 
      (*s)->dominated = false;
      ( *s ) -> next = NULL ; 
    } 
  
  else 
    { 
      
      /* traverse the linked list till the last node is reached */ 
      while ( q -> next != NULL ) 
	q = q -> next; 

      /* add a new node at the end */ 

      r = (dnode *)malloc ( sizeof ( struct dnode ) ) ; 
      r->o = (double *)malloc(nobjs * sizeof (double) );
      for(int j=0;j<nobjs;j++)
	r -> o[j] = vec[j]; 
      r->dominated = false;
      r -> next = NULL; 
      r -> prev = q; 
      q -> next = r; 
    } 
  
} 

/* displays the contents of the linked list */ 
void d_display ( struct dnode *q ) 
{ 
  /* traverse the entire linked list */ 
  while ( q != NULL )     
    { 
      if(1) // (q->dominated==false)
	{
	  for(int j=0;j<nobjs;j++)
	    printf("%.9e ", q->o[j]);
	  printf("\n");
	}
      q = q -> next;       
    } 
  
} 


int main(int argc, char **argv)
{
  int num;
  int nruns;
  int count;
  struct dnode *p;
  p=NULL;  
  int i;  
  char str[MAX_STR_LENGTH];
  
  error(argc!=4 && argc != 3,"./bound [<paramfile>] <datafile> <outfile>");
  
   
  /* read in the parameter file */
  if (argc == 4) {
      if((fp = fopen(argv[1], "rb")))
      {
	  fscanf(fp, "%s", str);
	  error(strcmp(str, "dim") != 0, "error in parameter file");
	  fscanf(fp, "%d", &nobjs);
	  fscanf(fp, "%s", str);
	  error(strcmp(str, "obj") != 0, "error in parameter file");
	  minmax1 = (int *)malloc(nobjs*sizeof(int));
	  for (i = 0; i < nobjs; i++) 
	  {
	      fscanf(fp, "%s", str);
	      error(str[0] != '-' && str[0] != '+', "error in parameter file");
	      if (str[0] == '-')
		  minmax1[i] = -1;
	      else
		  minmax1[i] = 1;
	  }
	  best = (double *)malloc(nobjs*sizeof(double));
	  worst = (double *)malloc(nobjs*sizeof(double));
	  fscanf(fp, "%s", str);error(strcmp(str, "phi") != 0, "error in parameter file");
	  fscanf(fp, "%lf", &phi);
	  error((phi<0),"phi should be a positive real number");
	  fclose(fp);
      }
      else
      {
	  fprintf(stderr,"Couldn't open param file\n");
	  exit(1);
      }
  }
  else {
	fp = fopen(argv[1], "r");
	error(fp == NULL, "data file not found");
	nobjs = determine_dim(fp);
	error(nobjs < 1, "error in data file");
	fclose(fp);
	minmax1 = (int *)malloc(nobjs*sizeof(int));
	for (i = 0; i < nobjs; i++) 
	    minmax1[i] = -1;
	best = (double *)malloc(nobjs*sizeof(double));
	worst = (double *)malloc(nobjs*sizeof(double));
	phi = 0.1;
  }
  

  /* read in each of the approximation sets */
  if((fp=fopen(argv[(argc == 4 ? 2 : 1)], "rb")))
    {
      check_file(fp, &nruns, &num);
      rewind(fp);
      while(!feof(fp))
	read_file(fp, &num, &p);      
    }
  else
    {
      fprintf(stderr,"Couldn't open %s", argv[(argc == 4 ? 2 : 1)]);
      exit(1);
    }
  
  for(i=0;i<nobjs;i++)
    {
      best[i] = p->o[i];
      worst[i] = p->o[i];
    }
  

  struct dnode *di = p;
  count=0;
  while( di !=NULL)
    {
      for(int i=0;i<nobjs;i++)
	{
	  if(di->o[i]*minmax1[i] > best[i]*minmax1[i])
	    best[i] = di->o[i];
	  if(di->o[i]*minmax1[i] < worst[i]*minmax1[i])
	    worst[i] = di->o[i];
	}
      count++;
      di = di->next;
    }
  //  fprintf(stderr, "count =%d\n",count);


  if(!(fp=fopen(argv[(argc == 4 ? 3 : 2)],"wb")))
    {
      fprintf(stderr, "Couldn't open %s for writing. Exiting\n",
	      argv[(argc == 4 ? 3 : 2)]);
      exit(1);
    }

  // fprintf(fp, "nadir "); // by Felipe 09/07/2018
  // for(i=0;i<nobjs;i++)
  //   fprintf(fp, "%.9e ", worst[i]);
  // fprintf(fp, "\n");
  // fprintf(fp, "ideal "); // by Felipe 09/07/2018
  // for(i=0;i<nobjs;i++)
  //   fprintf(fp, "%.9e ", best[i]);
  // fprintf(fp, "\n");
  fprintf(fp, "lower_bound ");
  for(i=0;i<nobjs;i++)
    {
      if(minmax1[i]==-1)
	fprintf(fp, "%.9e ", best[i]);
      else
	fprintf(fp, "%.9e ", worst[i]);
    }
  fprintf(fp, "\n");
  fprintf(fp, "upper_bound ");
  for(i=0;i<nobjs;i++)
    {
      if(minmax1[i]==-1)
	fprintf(fp, "%.9e ", worst[i]);
      else
	fprintf(fp, "%.9e ", best[i]);
    }
  fprintf(fp, "\n");  
  // fprintf(fp, "utopia_bound "); // by Felipe 09/07/2018
  // for(i=0;i<nobjs;i++)
  //   fprintf(fp, "%.9e ", best[i]+myabs(best[i]-worst[i])*minmax1[i]*phi);
  // fprintf(fp, "\n");

  // fprintf(fp, "nadir_bound "); // by Felipe 09/07/2018
  // for(i=0;i<nobjs;i++)
  //   fprintf(fp, "%.9e ", worst[i]-myabs(best[i]-worst[i])*minmax1[i]*phi);
  // fprintf(fp, "\n");

  fclose(fp);
  exit(0);
  return(0);

}

void  check_file(FILE  *fp, int  *no_runsp, int  *max_pointsp)
    /* check_file() function by Eckart Zitzler
       determines the maximum number of points and the number of runs
       for the data resp. the reference set file; if the array v is
       specified, the data read in will be stored in v
    */
{
    char  line[MAX_STR_LENGTH];
    int  i, j;
    int  new_run;
    int  no_points;
    double  number;

    no_points = 0;
    *max_pointsp = 0;
    *no_runsp = 0;
    new_run = 1;
    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) {
	if (sscanf(line, "%lf", &number) != 1)
	    new_run = 1;
	else {
	    if (new_run == 1)
	    {
		(*no_runsp)++;
		if (*max_pointsp < no_points)
		    *max_pointsp = no_points;
		no_points = 0;
	    }
	    new_run = 0;
	    i = 0;
	    for (j = 1; j < nobjs; j++) {
		while (line[i] != ' ' && line[i] != '\n' && line[i] != '\0')
		    i++;
		if((sscanf(&(line[i]), "%lf", &number)) <= 0)
		  {
		    fprintf(stderr,"error in data or reference set file");
		    exit(0);
		  }
		
		while (line[i] == ' ' && line[i] != '\0')
		    i++;
	    }
	    no_points++;
	}
    }
    if (*max_pointsp < no_points)
	*max_pointsp = no_points;
}


int  determine_dim(FILE  *fp)
{
    char  line[MAX_STR_LENGTH];
    int  i, no_obj;
    int  line_found, number_found;
    double  number;
    
    no_obj = 0;
    line_found = 0;
    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL && !line_found)
        line_found = sscanf(line, "%lf", &number);
    if (line_found) {
	i = 0;
	do {
	    no_obj++;
	    while (line[i] != ' ' && line[i] != '\n' && line[i] != '\0')
		i++;
	    number_found = sscanf(&(line[i]), "%lf", &number);
	    while (line[i] == ' ' && line[i] != '\0')
		i++;
	} while (number_found == 1);
    }
    
    return no_obj;
}

void  read_file(FILE  *fp, int  *no_pointsp, dnode **po)
{
  /* read_file() function by Eckart Zitzler */
  char  line[MAX_STR_LENGTH];
  int  i, j, k;
  int  reading;
  double  number;
  double *vector;

  vector=(double *)malloc(nobjs*sizeof(double));
  
  reading = 0;
  *no_pointsp = 0;
  //  printf("read_file\n");
  while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) 
    {
      //      printf("OK\n");
      k=0;
      if (sscanf(line, "%lf", &number) != 1)
	{
	  //	  printf("break\n");
	  if (reading) break;
	}
      else 
	{
	  reading = 1;
	  vector[k++] = number;
	  i = 0;
	  for (j = 1; j < nobjs; j++) 
	    {
	      while (line[i] != ' ' && line[i] != '\n' && line[i] != '\0')
		i++;
	      if((sscanf(&(line[i]), "%lf", &number)) <= 0)
		{
		    fprintf(stderr,"error in data or reference set file");
		    exit(0);
		}
	      
	      vector[k++] = number;
	      while (line[i] == ' ' && line[i] != '\0')
		i++;
	    }
	  (*no_pointsp)++;
	  // printf("%.9e", vector[0]);
	  d_append(po, vector);
	}
    } 
}

double myabs(double a)
{
  if(a>=0)
    return(a);
  else
    return(-a);
}



