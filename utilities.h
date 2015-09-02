#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <math.h>
#include "csparse.h"


#ifndef _utilities_h
#define _utilities_h


struct term_2 { //lista gia stoixeia 2 akrodektwn
	char type;
	char *name;
	char *pos_term;
	char *neg_term;
	double value;
	struct term_2 *nxt;	
};

struct term_2 *root_2;

struct term_3 { //lista gia stoixeia 3 akrodektwn
	char type;
	char *name;
	char *collector;
	char *base;
	char *emitter;
	double value;
	struct term_3 *nxt;	
};

struct term_4 { //lista gia stoixeia 4 akrodektwn
	char type;
	char *name;
	char *drain;
	char *gate;
	char *source;
	char *body;
	double Lvalue;
	double Wvalue;
	struct term_4 *nxt;	
};

typedef struct {      // struct pou krataei ta dedomena gia analysh DC Sweep
	int DC_flag;
	char *source;
	double start_value;
	double end_value;
	double inc;
} DC_sweep;


typedef struct {      // To struct plot krataei ta onomata twn komvwn pros ektypwsh
	int plot_flag;
	char **source;
	int count;    // o arithmos twn komvwn  pou prokeitai na ektypwthoun
} Plot;


struct term_2 *root_2; // root komvos stin lista me ta duo termatika
void init_2(); // arxikopoiisi listas me ta duo termatika
void insert_2(char type, char* name,char* pos_term,char* neg_term,double value);//eisagwgh stoixeiwn sth lista
struct term_3 *root_3; 
void init_3();
void insert_3(char type, char* name,char* collector,char* base,char* emitter,double value);
struct term_4 *root_4; 
void init_4();
void insert_4(char type, char* name,char* drain,char* gate,char* source,char* body,double Lvalue,double Wvalue);

//void prnt();
int read_netlist(char *filename, int *n, int *m2, char **method, double *itol, int *nz, int *cnz, DC_sweep *dc_sweep, Plot *plot);

void computeSparseMna(char **method, int *n, int *m2, double *itol, int *nz, int *cnz, DC_sweep *dc_sweep, Plot *plot);
void computeMna(char **method, int *n, int *m2, double *itol, DC_sweep *dc_sweep, Plot *plot);

int BiCG_solve(gsl_matrix *A_matrix, gsl_vector *b_vector, gsl_vector *x_vector, int DIM, double *itol_n);
int CG_solve(gsl_matrix *A_matrix, gsl_vector * b_vector,gsl_vector *x_vector, int DIM, double *itol_n);
void precond_create(gsl_matrix *A_matrix, gsl_vector *M, int DIM);
void precond_solve(gsl_vector *z, const gsl_vector *r, const gsl_vector *M, int DIM);
double norm_gslcalc(const gsl_vector *x, int DIM);
double vv_dot_gslcalc(const gsl_vector *x, const gsl_vector *y, int DIM);
void mv_prod_calc(gsl_vector *y, const gsl_matrix *A, const gsl_vector *x, int transpose, int DIM);

int BiCG_solve_Sparse(cs *C, double *RHS, double *x_vect, int DIM, double *itol_n); 
int CG_solve_Sparse(cs *C, double *RHS, double *x_vect, int DIM, double *itol_n);
void precond_spcreate(cs *C, double *M, int DIM);
void precond_spsolve(double *z, const double *r, const double *M, int DIM);
double norm_calc(const double *x, int DIM); 
double vv_dot_calc(const double *x, const double *y, int DIM);
void mv_spprod_calc(double *y, const cs *A, cs *C, const double *x, int transpose, int DIM);

#endif