/* filter.cc  (C) Joshua Knowles 2005

A program that reads in a file <datafile> consisting of a collection
of approximation sets and outputs a collection
of internally nondominated approximation sets. 
A separate parameter file <param> gives 
the dimension of the data and whether each objective should
be minimized or maximized.

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
      g++ filter.cc -o filter -lm -Wall -pedantic

   RUN:
      ./filter [<param>] <datafile> <outfile>


    The format of the parameter file <param> is

      dim <number>
      obj <+|-> <+|-> ...
      method <0|1>

    where dim specifies the number of objective dimensions
    obj specifies whether each objective should be minimized (-) or maximized (+)
    method specifies whether each approximation set is treated separately (0)
    or the nondominated set among all approximation sets is computed (1)

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

      
   The output of filter is a file <outfile> in
      the same format, but all points should 
      be internally nondominated within each set
      separated by the blank lines.

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

int *minmax1;
int nobjs;
int method;


void d_append (struct dnode **s, double *vec);
void d_display ( struct dnode *q );
void  check_file(FILE  *fp, int  *no_runsp, int  *max_pointsp);
int  determine_dim(FILE  *fp);
void  read_file(FILE  *fp, int  *no_pointsp, dnode **po);
int dominates(double *a, double *b, int *minmax1, int n);
bool are_identical(double *first, double *second, int *minmax1, int n);


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
  int i;  
  char str[MAX_STR_LENGTH];
  
  error(argc!=4 && argc != 3,"./filter [<paramfile>] <datafile> <outfile>");
  
  
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
      fscanf(fp, "%s", str);
      error(strcmp(str, "method") != 0, "error in parameter file");
      fscanf(fp, "%d", &method);
      error(method != 0 && method != 1, "error in parameter file");
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
      method = 1;
  }
  

  /* read in each of the approximation sets */
  if((fp=fopen(argv[(argc == 4 ? 2 : 1)], "rb")))
  {
      check_file(fp, &nruns, &num);
      po = (dnode **)malloc(nruns*sizeof(dnode *));
      rewind(fp);
      for(int j=0;j<nruns;j++)
      {
	  po[j]=NULL;
	  read_file(fp, &num, &(po[j]));      
      }
  }
  else
  {
      fprintf(stderr,"Couldn't open %s", argv[(argc == 4 ? 2 : 1)]);
      exit(1);
  }
  if(!(fp=fopen(argv[(argc == 4 ? 3 : 2)],"wb")))
  {
      fprintf(stderr,"Couldn't open %s for writing\n",
	      argv[(argc == 4 ? 3 : 2)]);
      exit(0);
  }
  
  for(i=0; i<nruns;i++)
  {
      struct dnode *di = po[i];
      while( di !=NULL)
      {
	  struct dnode *dj = po[i];
	  while( dj !=NULL)
	  {
	      if( di != dj )
	      {		  
		  if(dominates( di->o, dj->o, minmax1, nobjs)==-1)
		  {
		      di->dominated=true;
		      break;
		  }		
	      }
	      dj = dj->next;
	  }
	  di = di->next;
      }
      
      di = po[i];
      while( di !=NULL)
      {
	  struct dnode *dj = di->next;
	  while( dj !=NULL)
	  {
	      if(are_identical( di->o, dj->o, minmax1, nobjs))
	      {
		  dj->dominated=true;
		  break;
	      }
	      dj = dj->next;
	  }
	  di = di->next;
      }      
      
      struct dnode *q = po[i];
      while ( q != NULL )     
      { 
	  if(q->dominated==false)
	  {
	      for(int j=0;j<nobjs;j++)
		  fprintf(fp, "%.9e ", q->o[j]);
	      fprintf(fp, "\n");
	  }
	  q = q -> next;       
      }       
      fprintf(fp, "\n");
  }
  
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
        if (sscanf(line, "%lf", &number) != 1) {
	    if (method == 0)
		new_run = 1;
        }
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
      if (sscanf(line, "%lf", &number) != 1) {
	  if (reading)
	  {
	      //	  printf("break\n");
	      break;
	  }
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
	  if (method != 0)  reading = 0;
	}
    } 
}

int dominates(double *a, double *b, int *minmax1, int n)
{
  // Returns 1 if a dominates b
  // Returns -1 if b dominates a
  // Returns 0 if a==b or a and b are mutually nondominating (a~b)

  // This dominance relationship is measured using the dimensions 
  // defined in the array minmax1[]. If minmax1[i]=1 then objective i
  // is compared assuming maximization of i. If minmax1[i]=-1 then
  // objective i is compared assuming minimization of i. If
  // minmax1[i]=0 then no comparison is made on objective i. The 
  // total number of dimensions which could potentially be compared
  // is given by the argument n.

  double diff;
  int abb=0; // counts number of dimensions where a beats b
  int bba=0;
  
  for(int i=0;i<n;i++)
    {
      if(!(minmax1[i]==0))
	{
	  diff=a[i]-b[i];
	  if(diff>0)
	    {
	      if(minmax1[i]==1)
		abb++;
	      else if (minmax1[i]==-1)
		bba++;
	      else
		{
		  fprintf(stderr, "minmax1 out of range\n");
		  exit(0);
		}
	    }
	  else if(diff<0)
	    {
	      if(minmax1[i]==1)
		bba++;
	      else if (minmax1[i]==-1)
		abb++;
	      else
		{
		  fprintf(stderr, "minmax1 out of range\n");
		  exit(0);
		}
	    }
	}
      if((bba>0)&&(abb>0))
	return(0);
    }
  if(abb>0)
    return(1);
  else if(bba>0)
    return(-1);
  else
    return(0);

}


bool are_identical(double *first, double *second, int *minmax1, int n)
{
  // Returns true if two vectors are the same in all of 
  // the nonzero components of the minmax1 array, false otherwise.
  // Returns false if all components of minmax1 array are zero.

  int i;
  bool onenonzero=false;
  
  for(i=0;i<n;i++)
    {
      if(minmax1[i]!=0)
	{
	  if(first[i]!=second[i])
	    return(false);      
	  onenonzero=true;
	}
    }
  if(onenonzero)
    return(true);
  else
    return(false);
}


