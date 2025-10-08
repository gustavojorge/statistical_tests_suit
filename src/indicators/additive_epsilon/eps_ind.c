/*===========================================================================*
 * eps_ind.c: implements the unary epsilon indicator as proposed in
 *            Zitzler, E., Thiele, L., Laumanns, M., Fonseca, C., and
 *            Grunert da Fonseca, V (2003): Performance Assessment of
 *            Multiobjective Optimizers: An Analysis and Review. IEEE
 *            Transactions on Evolutionary Computation, 7(2), 117-132.
 *
 * Compile:
 *   gcc -lm -o eps_ind eps_ind.c
 *
 * Usage:
 *   eps_ind [<param_file>] <data_file> <reference_set> <output_file>
 *
 *   <param_file> specifies the name of the parameter file for eps_ind; the
 *     file has the following format:
 *
 *       dim <integer>
 *       obj <+|-> <+|-> ...
 *       method <0|1>
 *
 *     The first line defines the number of objectives, the second for each
 *     objective whether it is minimized (-) or maximized, and the third
 *     line determines whether the additive (0) or the multiplicate (1)
 *     version of the epsilon indicator is used.
 *     If the parameter file is omitted, default parameters are taken (all
 *     objective are to be minimized, method = 0) and the number of objectives
 *     is determined from the data file.
 *
 *   <data_file> specifies a file that contains the output of one or
 *     several runs of a selector/variator pair; the format corresponds to
 *     the one defined in the specification of the PISA monitor program.
 *
 *   <reference_set> is the name of a file that contains the reference set
 *     according to which the indicator values are calculated; the file
 *     format is the same as for the data file.
 *
 *   <output_file> defines the name of the file to which the computed
 *     indicator values are written to.
 *
 * IMPORTANT:
 *   To make the epsilon indicator work for mixed optimization problems
 *   where some objectives are to be maximized while others are to be
 *   minimized, in the case of minimization the value -epsilon (for the
 *   additive case) resp. 1/epsilon (for the multiplicative version) is
 *   considered and returned. Example: suppose f1 is to be minimized and
 *   f2 to be maximized, and the multiplicative epsilon value computed by
 *   this program is 3.0; this means that the considered nondominated front
 *   needs to be multiplied by 1/3 for all f1 values and by 3 for all
 *   f2 values. Thus, independently of which type of problem one considers
 *   (minimization, maximization, mixed minimization/maximization), a lower
 *   indicator value corresponds to a better approximation set.
 *
 * Author:
 *   Eckart Zitzler, February 3, 2005 / last update August 9, 2005
 */

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define error(X,Y)  if (X) fprintf(stderr, Y "\n"), exit(1)

#define MAX_LINE_LENGTH  2048 /* maximal length of lines in the files */
#define MAX_STR_LENGTH  256 /* maximal length of strings in the files */

int  dim;  /* number of objectives */
int  *obj;  /* obj[i] = 0 means objective i is to be minimized */
int  method;  /* 0 = additive, 1 = multiplicative */


double  calc_ind_value(double  *a, int  size_a, double  *b, int  size_b)
{
    int  i, j, k;
    double  eps, eps_j, eps_k, eps_temp;

    if (method == 0)
	eps = DBL_MIN;
    else
	eps= 0;
    
    for (i = 0; i < size_a; i++) {
	for (j = 0; j < size_b; j++) {
	    for (k = 0; k < dim; k++) {
		switch (method) {
		case 0:
		    if (obj[k] == 0)
			eps_temp = b[j * dim + k] - a[i * dim + k];
		    else
			eps_temp = a[i * dim + k] - b[j * dim + k];
		    break;
		default:
		    error((a[i * dim + k] < 0 && b[j * dim + k] > 0) ||
			  (a[i * dim + k] > 0 && b[j * dim + k] < 0) ||
			  a[i * dim + k] == 0 || b[j * dim + k] == 0,
			  "error in data file");
		    if (obj[k] == 0)
			eps_temp = b[j * dim + k] / a[i * dim + k];
		    else
			eps_temp = a[i * dim + k] / b[j * dim + k];
		    break;
		}
		if (k == 0)
		    eps_k = eps_temp;
		else if (eps_k < eps_temp)
		    eps_k = eps_temp;
	    }
	    if (j == 0)
		eps_j = eps_k;
	    else if (eps_j > eps_k)
		eps_j = eps_k;
	}
	if (i == 0)
	    eps = eps_j;
	else if (eps < eps_j)
	    eps = eps_j;
    }
    
    return eps;
}

void  read_params(FILE  *fp)
{
    char str[MAX_STR_LENGTH];
    int  i;
    
    fscanf(fp, "%s", str);
    error(strcmp(str, "dim") != 0, "error in parameter file");
    fscanf(fp, "%d", &dim);
    error(dim <= 0, "error in parameter file");
    obj = malloc(dim * sizeof(int));
    error(obj == NULL, "memory overflow");

    fscanf(fp, "%s", str);
    error(strcmp(str, "obj") != 0, "error in parameter file");
    for (i = 0; i < dim; i++) {
	fscanf(fp, "%s", str);
	error(str[0] != '-' && str[0] != '+', "error in parameter file");
	if (str[0] == '-')
	    obj[i] = 0;
	else
	    obj[i] = 1;
    }

    fscanf(fp, "%s", str);

    error(strcmp(str, "method") != 0, "error in parameter file");
    fscanf(fp, "%d", &method);
    error(method != 0 && method != 1, "error in parameter file");
}

void  set_params(void)
{
    int  i;
    
    obj = malloc(dim * sizeof(int));
    error(obj == NULL, "memory overflow");
    for (i = 0; i < dim; i++)
	obj[i] = 0;
    method = 0;
}

void  check_file(FILE  *fp, int  *no_runsp, int  *max_pointsp)
    /* determines the maximum number of points and the number of runs
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
	    for (j = 1; j < dim; j++) {
		while (line[i] != ' ' && line[i] != '\n' && line[i] != '\0')
		    i++;
		error(sscanf(&(line[i]), "%lf", &number) <= 0,
		      "error in data or reference set file");
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

void  read_file(FILE  *fp, int  *no_pointsp, double  *points)
{
    char  line[MAX_STR_LENGTH];
    int  i, j, k;
    int  reading;
    double  number;

    k = 0;
    reading = 0;
    *no_pointsp = 0;
    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) {
        if (sscanf(line, "%lf", &number) != 1) {
	    if (reading)
	        break;
	}
	else {
	    reading = 1;
	    points[k++] = number;
	    i = 0;
	    for (j = 1; j < dim; j++) {
		while (line[i] != ' ' && line[i] != '\n' && line[i] != '\0')
		    i++;
		error(sscanf(&(line[i]), "%lf", &number) <= 0,
		      "error in data or reference set file");
		points[k++] = number;
		while (line[i] == ' ' && line[i] != '\0')
		    i++;
	    }
	    (*no_pointsp)++;
	}
    } 
}

int  main(int  argc, char  *argv[])
{
    int  i;
    int  no_runs;  /* number of runs */
    int  max_points;  /* maximum number of points per run */
    int  ref_set_size;  /* number of points in the reference set */
    int  curr_run_size;  /* number of points associated with the current run */
    double  *ref_set;  /* reference set */
    double  *curr_run; /* objective vectors fur current run */
    double  ind_value;
    FILE  *fp, *out_fp;
    
    error(argc != 4 && argc != 5,
	  "Epsilon indicator - wrong number of arguments:\neps_ind [parFile] datFile refSet outFile");

    /* set parameters */
    if (argc == 5) {
	fp = fopen(argv[1], "r");
	error(fp == NULL, "parameter file not found");
	read_params(fp);
	fclose(fp);
    }
    else {
	fp = fopen(argv[1], "r");
	error(fp == NULL, "data file not found");
	dim = determine_dim(fp);
	error(dim < 1, "error in data file");
	fclose(fp);
	obj = malloc(dim * sizeof(int));
	error(obj == NULL, "memory overflow");
	for (i = 0; i < dim; i++)
	    obj[i] = 0;
	method = 0;	
    }

    /* read reference set */
    if (argc == 5)
	fp = fopen(argv[3], "r");
    else
	fp = fopen(argv[2], "r");
    error(fp == NULL, "reference set file not found");
    check_file(fp, &no_runs, &max_points);
    error(no_runs != 1 || max_points < 1, "error in reference set file");
    ref_set = malloc(dim * max_points * sizeof(double));
    error(ref_set == NULL, "memory overflow");
    rewind(fp);
    read_file(fp, &ref_set_size, ref_set);
    fclose(fp);
    no_runs = 0;
    max_points = 0;
    
    /* check data file */
    if (argc == 5)
	fp = fopen(argv[2], "r");
    else
	fp = fopen(argv[1], "r");
    error(fp == NULL, "data file not found");
    check_file(fp, &no_runs, &max_points);
    error(no_runs < 1 || max_points < 1, "error in data file");
    curr_run = malloc(dim * max_points * sizeof(double));
    rewind(fp);

    /* process data */
    if (argc == 5)
	out_fp = fopen(argv[4], "w");
    else
	out_fp = fopen(argv[3], "w");
    error(out_fp == NULL, "output file could not be generated");
    while (no_runs > 0) {
	read_file(fp, &curr_run_size, curr_run);
	ind_value = calc_ind_value(ref_set, ref_set_size,
				   curr_run, curr_run_size);
	fprintf(out_fp, "%.9e\n", ind_value);
	no_runs--;
    }
    fclose(out_fp);
    fclose(fp);
}
