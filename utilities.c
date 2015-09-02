#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "utilities.h"
#include "csparse.h"
#include "hash_table.h"


bool flag = false;   // boolean metavlhth gia anakalypsh tou an yparxei komvos geiwshs sto netlist

int read_netlist(char *filename, int *n, int *m2, char **method, double *itol, int *nz, int *cnz, DC_sweep *dc_sweep, Plot *plot) { //diavazei to netlist apo to arxeio
   
   FILE *fp;
   int choice=0;
   double value, Lvalue, Wvalue;
   char *line, *options, *temp, *temp2;
   char type, t;
   char *type_name, *name, *pos_node, *neg_node;
   char *drain, *gate, *source, *body;
   char *collector, *base, *emitter;
   char *t_str, *i_str;	// metavlhtes gia anagnwsh parametrwn pwl
   char *chartemp;
   double t1,i1;

   char **nodes;
   int i;
 

   /* arxikopoihsh listwn */
   init_2();
   init_3();
   init_4();
     
   line = (char *)malloc(100*sizeof(char));

   fp = fopen(filename, "r"); 
   if( fp == NULL ) {
	printf("cannot open file\n");
	exit(0);
   }    
  
   while( fgets(line, 50, fp) != NULL) {	
	type_name = strtok(line, " ");
	type      = *type_name;
	
	if (type == '*') continue;

	else if (type == '.') {
		if( strcmp((type_name + 1), "OPTIONS") == 0 ) {
			temp2 = strdup(strtok(NULL, "\n"));
			if(temp2==NULL)
			{
			  printf("NO TEMP2");
			}
			//printf("%s %d\n", temp2, strlen(temp2));
			options = strdup(temp2); 
			//printf("options = %s, %d\n", options, strlen(options)); 
			temp = strtok(temp2, " ");
			printf("%s", options);
 			//printf("%d %d", strlen(options), strlen(temp));
			if ( strcmp(temp, "SPARSE") == 0 ) { 
				if ( strcmp(options, "SPARSE") == 0 ){
				  strcpy(*method, options);
				  choice=1;
				}
				else if( strcmp(options, "SPARSE SPD") == 0 ){
				  strcpy(*method, options);
				  choice=2;
				}
				else if( strcmp(options, "SPARSE ITER") == 0 ){
				  strcpy(*method, options);
				  choice=3;
				}
				else if( strcmp(options, "SPARSE SPD ITER") == 0 ){
				  //printf("????");
				  strcpy(*method, options);
				  choice=4;
				}
				else {
					*method = strtok(options, " ");
					*method = strtok(NULL, " ");
					*method = strtok(NULL, " ");
					if( strcmp(*method, "ITOL") == 0 ) {
						*method = strdup("SPARSE ITER");
						*itol = atof(strtok(NULL, "\n"));
						continue;
					}
					else {
					  printf("ginetai malakia?\n");
					//	*method = strtok(NULL, " ");
					//	if( strcmp(*method, "ITOL") == 0 ) {
					//		*method = strdup("SPARSE SPD ITER");
					//		*itol = atof(strtok(NULL, "\n"));
					//	}
					}
				}
			}
			
			else {  
				if( strcmp(options, "SPD") == 0 ){
				  strcpy(*method, options); 
				  choice=5;
				}
				else if( strcmp(options, "ITER") == 0 ){
				  strcpy(*method, options);
				  choice=6;
				}
				else if( strcmp(options, "SPD ITER") == 0 ){
				  printf("mpainei");
				  strcpy(*method, options);
				  choice=7;
				}
				else {
				  printf(">>>>>");
					*method = strtok(options, " ");
			 		*method = strtok(NULL, " ");
					if( strcmp(*method, "ITOL") == 0 ) {
						*method = strdup("ITER");
						*itol = atof(strtok(NULL, "\n"));
						continue;
					}
					else {
						*method = strtok(NULL, " ");
						if( strcmp(*method, "ITOL") == 0 ) {
							*method = strdup("ITER SPD");
							*itol = atof(strtok(NULL, "\n"));
						}
					}
				}
			}
		}
		else if( strcmp((type_name + 1), "DC") == 0 ) {
			dc_sweep->DC_flag     = 1;
			dc_sweep->source      = strdup(strtok(NULL, " "));
			dc_sweep->start_value = atof(strtok(NULL, " "));
			dc_sweep->end_value   = atof(strtok(NULL, " "));
			dc_sweep->inc         = atof(strtok(NULL, "\n"));
		}
		//else if( strcmp((type_name + 1), "TRAN") == 0 ) {
		//	trans_analysis->trans_flag = 1;
		//	trans_analysis->time_step  = atof(strtok(NULL, " "));
		//	trans_analysis->fin_time   = atof(strtok(NULL, " "));
		//}
		else if( strcmp((type_name + 1), "PLOT") == 0 ) {
			plot->plot_flag = 1;
			plot->source    = (char **) malloc(sizeof(char *));

			chartemp      = strtok(NULL, " ");
			*plot->source = (char *)strdup(strtok(NULL, " "));

			i = 0;
			while(1) {
				i++;
				chartemp = strtok(NULL, " "); 
				if(chartemp == NULL)   break;
				*(plot->source+i) = (char *)strdup(strtok(NULL, " \n"));	
			}
			plot->count = i;
		}
	}
	
	else if( (type == 'v') || (type == 'V') || (type == 'r') || (type == 'R') || (type == 'i') || (type == 'I') || (type == 'c') || (type == 'C')  || (type == 'l') || (type == 'L') ) {
		name     = type_name + 1 ;
        	pos_node = strtok(NULL, " ");
        	neg_node = strtok(NULL, " ");
        	value    = atof(strtok(NULL, " ")); 

		
		if( (strcmp(pos_node, "0") == 0) || (strcmp(neg_node, "0") == 0) ) {
			flag = true;  
		}
		
		if( (type == 'R') || (type == 'r') ) {
			// an o enas komvos antistashs einai syndedemenos 
			// ayksanoume ta non-zeros kata 1
			//if( (strcmp(pos_node, "0") == 0) || (strcmp(neg_node, "0") == 0) )   *nz += 1;
			// diaforetika ayksanoume ta non-zeros kata 4
			//else *nz += 4; 
			*nz += 4;
			//printf("eimai sto r kai exw %d\n", *nz);
		}

		if( (type == 'V') || (type == 'v') ) {
			 *m2 = *m2 + 1;

			// an o enas komvos ths phghs tashs einai syndedemenos 
			// ayksanoume ta non-zeros kata 2
			//if( (strcmp(pos_node, "0") == 0) || (strcmp(neg_node, "0") == 0) )   *nz += 2; 
			// diaforetika ayksanoume ta non-zeros kata 4
			//else *nz += 4;
			*nz += 2;
			printf("einai sto v kai exw %d\n", *nz);
		}

		if( (type == 'L') || (type == 'l') ) {
			 //*cnz += 1;
			 *m2 = *m2 + 1;
			 *nz += 4;
			 printf("eimai sto l kai exw %d\n", *nz);
		}		
		
		if( (type == 'C') || (type == 'c') ) {
			// an o enas komvos ths phghs tashs einai syndedemenos 
			// ayksanoume ta non-zeros kata 2
			//if( (strcmp(pos_node, "0") == 0) || (strcmp(neg_node, "0") == 0) )   *cnz += 1; 
			// diaforetika ayksanoume ta non-zeros kata 4
			//else *cnz += 4;
		}

		//*n = insert_table(pos_node, *n);
		*n = addEntry(pos_node, *n);
		//*n = insert_table(neg_node, *n);
		*n = addEntry(neg_node, *n);

		insert_2(type, name, pos_node, neg_node, value);
	}

	else if( (type == 'm') || (type == 'M') ) {
		name   = type_name + 1 ;
		drain  = strtok(NULL, " ");
		gate   = strtok(NULL, " ");
		source = strtok(NULL, " ");
		body   = strtok(NULL, " ");
		Lvalue = atof(strtok(NULL, " "));
		Wvalue = atof(strtok(NULL, " "));

		insert_4(type, name, drain, gate, source, body, Lvalue, Wvalue);
	}

	else if( (type == 'q') || (type == 'Q') ) {
		name      = type_name + 1;
		collector = strtok(NULL, " ");
		base      = strtok(NULL, " ");
		emitter   = strtok(NULL, " ");
		value     = atof(strtok(NULL, " "));

		insert_3(type, name, collector, base, emitter, value);
	}
   }

   if (flag == false) {
	printf("Error. No ground node exists.\nProgram exit\n");
	exit(0);
   }
   
   free(line);
   return choice;
}


void init_2() {
	root_2 = (struct term_2*)malloc(sizeof(struct term_2));
	root_2->nxt = NULL;

}

void insert_2(char type, char *name, char *pos_term, char *neg_term, double value) {
	struct term_2 *new_node_2 = (struct term_2*)malloc(sizeof(struct term_2)); //new node stin lista me ta duo termatika

	new_node_2->type = type;
	new_node_2->name = strdup(name);
	
	new_node_2->pos_term   = strdup(pos_term);
	new_node_2->neg_term   = strdup(neg_term);
	new_node_2->value      = value;

	new_node_2->nxt = root_2;
	root_2          = new_node_2;
}

void init_3() {
	root_3 = (struct term_3*)malloc(sizeof(struct term_3));
	root_3->nxt = NULL;
}

void insert_3(char type, char* name, char* collector, char* base, char* emitter, double value) {
	struct term_3 *new_node_3 = (struct term_3*)malloc(sizeof(struct term_3));; //new node stin lista me ta 3 termatika

	new_node_3->type = type;
	new_node_3->name = strdup(name);
	
	new_node_3->collector = strdup(collector);
	new_node_3->base      = strdup(base);
	new_node_3->emitter   = strdup(emitter);
	new_node_3->value     = value;

	new_node_3->nxt = root_3;
	root_3          = new_node_3;
}

void init_4() {
	root_4 = (struct term_4*)malloc(sizeof(struct term_4));
	root_4->nxt = NULL;
}

void insert_4(char type, char* name, char* drain, char* gate, char* source, char* body, double Lvalue, double Wvalue){
	struct term_4 *new_node_4 = (struct term_4*)malloc(sizeof(struct term_4));; //new node stin lista me ta 4 termatika

	new_node_4->type   = type;
	new_node_4->name   = strdup(name);
	new_node_4->drain  = strdup(drain);
	new_node_4->gate   = strdup(gate);
	new_node_4->source = strdup(source);
	new_node_4->body   = strdup(body);
	new_node_4->Lvalue = Lvalue;
	new_node_4->Wvalue = Wvalue;

	new_node_4->nxt = root_4;
	root_4          = new_node_4;
}


void prnt() {//ektupwsh listwn
	struct term_2 *look = root_2;
	
	for(look=root_2; look->nxt != NULL; look=look->nxt) {
              // printf("%c%s %s %s %.6f %s %f %f %f %f %f %f\n\n", look->type, look->name, look->pos_term, look->neg_term, look->value, look->trans.spec, *look->trans.PWL.i1, *look->trans.PWL.i2, look->trans.EXP.td1, look->trans.EXP.td2, look->trans.EXP.tc1, look->trans.EXP.tc2); 
	}

	struct term_3 *look3;
	for(look3=root_3; look3->nxt != NULL; look3=look3->nxt) { 
               printf("%c%s\t%s\t%s\t%s\t%.6f\n\n", look3->type, look3->name, look3->collector, look3->base, look3->emitter, 	look3->value);
	}

	struct term_4 *look4;
	for(look4=root_4; look4->nxt != NULL; look4=look4->nxt) {
               printf("%c%s\t%s\t%s\t%s\t%s\t%.6f\t%.2f\n\n", look4->type, look4->name, look4->drain, look4->gate, 		  			look4->source, look4->body, look4->Lvalue, look4->Wvalue);
	}
}




void computeSparseMna(char **method, int *n, int *m2, double *itol, int *nz, int *cnz, DC_sweep *dc_sweep, Plot *plot)
{
  //Dc_Sweep dc_sweep;
  printf("\n*******SPARSE WITH %s METHOD, %d NODES, %d M2, %f ITOL, %d NZ, %d CNZ*******\n", *method, *n , *m2 , *itol, *nz, *cnz);
  struct term_2 *look;
	
	int k = 0, DIM, n1;
	int inc = 0;
	int i, pos, neg, index;
        int Vindex_in_RHS = 0;
	
	double *RHS;
	double *RHS_copy;
	cs *CSA;
	cs *CSC;
	cs *C;
	
	double *x, *x_temp;   // vohthhtiko dianysma
	double *x_vector;
	double val;
	
	css *S;	     
	csn *N;


	FILE *fp;
    FILE *fp1;
	 //S = cs_calloc (1, sizeof (css));
	 // N = cs_calloc (1, sizeof (csn));
	
	n1=*n;
	x = (double *)malloc(n1*sizeof(double));


	//b_vector = gsl_vector_alloc(DIM);
	//x_vector = gsl_vector_alloc(DIM);
	
	DIM = *n + *m2;
	printf("\nDIM : %d\n", DIM);
	CSA      = cs_spalloc(DIM, DIM, *nz, 1, 1);
	CSA->nz  = *nz;
	CSC      = cs_spalloc(DIM, DIM, *nz, 1, 1);
	CSC->nz  = *cnz;

	RHS      = (double *)malloc((DIM)*sizeof(double));
	RHS_copy = (double *)malloc((DIM)*sizeof(double));
	x_vector = (double *)malloc((DIM)*sizeof(double));

	
	for(look=root_2; look->nxt != NULL; look=look->nxt) {
	  pos = -1;
	  neg = -1;
	  if( (look->type == 'R') || (look->type == 'r') )
	  {
	    pos = searchHashTable(look->pos_term);
            neg = searchHashTable(look->neg_term);
            
            //printf("\n Sparse computation (R) \n");
            if( (pos != -1) && (neg != -1) ){
                CSA->i[k] = pos;
                CSA->p[k] = pos;
                CSA->x[k] = 1.0/look->value;
                k++;
                CSA->i[k] = pos;
                CSA->p[k] = neg;
                CSA->x[k] = -1.0/look->value;
                k++;
                CSA->i[k] = neg;
                CSA->p[k] = pos;
                CSA->x[k] = -1.0/look->value;
                k++;
                CSA->i[k] = neg;
                CSA->p[k] = neg;
                CSA->x[k] = 1.0/look->value;
                k++;
                
            }
            else if(pos != -1){
                CSA->i[k] = pos;
                CSA->p[k] = pos;
                CSA->x[k] = 1.0/look->value;
                k++;
            }
            else if(neg != -1){
                CSA->i[k] = neg;
                CSA->p[k] = neg;
                CSA->x[k] = 1.0/look->value;
                k++;
            }
            
        }
	  else if( (look->type == 'I') || (look->type == 'i') )
	  {
            
            pos = searchHashTable(look->pos_term);
            neg = searchHashTable(look->neg_term);
            //printf("\n Value1 : %lf \n",look->value);
            //printf("\n pos : %d \n",pos);
            //printf("\n neg : %d \n",neg);
            if ( ((look->type == 'V') || (look->type == 'v')) && strcmp(dc_sweep->source+1, look->name) == 0 ) 
				Vindex_in_RHS = (*n + inc);
            //printf("\n Sparse computation (I) \n");
            
            if( (pos != -1) && (neg != -1) )
            {
                RHS[pos] -= look->value;
                RHS[neg] +=  look->value;
                
            }
            else if(pos != -1)
            {
                RHS[pos] -= look->value;
                
            }
            else if(neg != -1)
            {
                RHS[neg] +=  look->value;
                
            }
            
        }
	  else if( (look->type == 'V') || (look->type == 'v') || (look->type == 'l') || (look->type == 'L'))
	  {
            
            pos = searchHashTable(look->pos_term);
            neg = searchHashTable(look->neg_term);
	    
            if ( ((look->type == 'V') || (look->type == 'v')) && strcmp(dc_sweep->source+1, look->name) == 0 ) 
				Vindex_in_RHS = (*n + inc);
            //printf(" Value1 : %lf \n",look->value);
            //printf(" pos : %d \n",pos);
            //printf(" neg : %d \n",neg);
            //printf("\n Sparse computation (V) \n");
            
            if( (pos != -1) && (neg != -1) )
            {
                CSA->i[k] = pos;
                CSA->p[k] = *n+inc;
                CSA->x[k] = 1;
                k++;
                CSA->i[k] = *n+inc;
                CSA->p[k] = pos;
                CSA->x[k] = 1;
                k++;
                CSA->i[k] = neg;
                CSA->p[k] = *n+inc;
                CSA->x[k] = -1;
                k++;
                CSA->i[k] = *n+inc;
                CSA->p[k] = neg;
                CSA->x[k] = -1;
                k++;
            
                
            }
            else if(pos != -1)
            {
                CSA->i[k] = pos;
                CSA->p[k] = *n+inc;
                CSA->x[k] = 1;
                k++;
                
                CSA->i[k] = *n+inc;
                CSA->p[k] = pos;
                CSA->x[k] = 1;
                k++;
                
                
            }
            else if(neg != -1)
            {
                CSA->i[k] = neg;
                CSA->p[k] = *n+inc;
                CSA->x[k] = -1;
                k++;
                CSA->i[k] = *n+inc;
                CSA->p[k] = neg;
                CSA->x[k] = -1;
                k++;

                
            }
            if((look->type == 'V') || (look->type == 'v'))
            {
                    RHS[*n+inc] = look->value;
		    //printf("%d\n", *n+inc-1);
		    //printf("type:%c\t b[%d]:%f\n",look->type,*n+inc-1,look->value);
            }
            inc++;
            
        }
           

	}	
	fp = fopen("MNAmatrixSparse", "w");
	//cs_print(CSA, "MNAmatrixSparse", 0);
	//exit(0);
    
    /*COPY OF RHS************************/
	for(i=0; i<DIM; i++) 
	{
	  RHS_copy[i] = RHS[i];
	  //printf("b[%d] = %f\n", i, RHS[i]);
	}
	//cs_print(C, "MNAmatrixSparse", 0);
    
    
	 C=cs_compress(CSA);
	//cs_spfree(CSA);
	cs_dupl(C);
	if ( strcmp(*method, "SPARSE") == 0 ) {
	  printf("\n------LU Sparse------\n\n");
	  S = cs_sqr(2, C, 0);
	  if (S == NULL) { 
	    printf("EMPTY S\n");
	    exit(0);
	  }
	  N = cs_lu(C, S, 1);
	 if (N == NULL) { 
	    printf("EMPTY N\n");
	    exit(0);
	  }
	  //cs_spfree(C);
	  cs_ipvec(N->pinv, RHS_copy, x, DIM);
	  cs_lsolve(N->L, x);
	  cs_usolve(N->U, x);
	  cs_ipvec(S->q, x, RHS_copy, DIM);
		
	  fprintf(fp,"\n#-------LU_SPARSE-------\n");
	  fprintf(fp,"\n#-------Vector: x_vector------- \n");
	  for(i=0; i<DIM; i++){
	    fprintf(fp, "#x[%d] = %f\n",i, RHS_copy[i]);   // h lysh apothhkeyetai sto dianysma deksiou melous b
	    printf("x[%d] = %#9.6f\n",i, RHS_copy[i]); 
	  }
	}
	else if ( strcmp(*method, "SPARSE SPD") == 0 ) {
	  printf("\n-------CHOLESKY SPARSE--------\n");
	  S = cs_schol(1, C);
	  N = cs_chol(C, S);
	  
	  cs_ipvec(S->pinv, RHS_copy, x, DIM);
	  cs_lsolve(N->L, x);
	  cs_ltsolve(N->L, x);
	  cs_pvec(S->pinv, x, RHS_copy, DIM);
	  
	  fprintf(fp,"\n-------CHOLESKY_SPARSE-------\n");
	  fprintf(fp,"\n-------Vector: x_vector ------- \n");
	  for(i=0; i<DIM; i++){
	    fprintf(fp, "x[%d] = %#9.6f\n",i, RHS_copy[i]);   // h lysh apothhkeyetai sto dianysma deksiou melous b
	    printf("x[%d] = %#9.6f\n",i, RHS_copy[i]); 
	  }
	  
	}

	else if ( strcmp(*method, "SPARSE ITER") == 0 ) {
	  
	  fprintf(fp,"\n------BiCG SPARSE--------\n");
	  BiCG_solve_Sparse(C, RHS, x_vector, DIM, itol);
	  fprintf(fp,"\n -------Vector: x_vector ------- \n");
	  for(i=0; i<DIM; i++)
	  {
	    fprintf(fp,"x[%d] = %#9.6f\n",i, x_vector[i]);
	    printf("x[%d] = %#9.6f\n", i, x_vector[i]);
	  }
	  
	}

	else if ( strcmp(*method, "SPARSE SPD ITER") == 0 ) { 
	  
	  fprintf(fp,"\n--------CG SPARSE--------\n");
	  CG_solve_Sparse(C, RHS, x_vector, DIM, itol);
	  fprintf(fp,"\n -------Vector: x_vector ------- \n");
	  for(i=0; i<DIM; i++)
	  {
	    fprintf(fp,"x[%d] = %#9.6f\n",i, x_vector[i]);
	    printf("x[%d] = %#9.6f\n", i, x_vector[i]);
	  }
	}
	
	printf("\nDCSWEEP : %d flag, %d id, %s Source, %f Step, %f Start Value\n", dc_sweep->DC_flag, Vindex_in_RHS, dc_sweep->source, dc_sweep->inc, dc_sweep->start_value);
	if(dc_sweep->DC_flag == 1)
	{
	    if(strcmp(*method, "SPARSE SPD ITER") == 0)
	    {
          fprintf(fp,"\n-------DC_Sweep:CG SPARSE---------\n");
	      for(val = dc_sweep->start_value; val <= dc_sweep->end_value; val = val + dc_sweep->inc)
	      {
            
              for (i = 0; i < DIM; i++) {
                  RHS_copy[i] = RHS[i];
                  //printf("b[%d] = %#9.6f\n", i , RHS[i]);
                  x_vector[i] = 0.0;
              }
              RHS_copy[Vindex_in_RHS] = val;
              CG_solve_Sparse(C, RHS_copy, x_vector, DIM, itol);
              //fprintf(fp,"\n -------Vector: x_vector ------- \n");
              fprintf(fp,"Sweep source current at %#9.6f :   Node %s value: %.6f \n", val, dc_sweep->source, x_vector[Vindex_in_RHS]);
              // printf("x[%d] = %#9.6f\n", i, x_vector[i]);
		
	      }
	      cs_spfree(CSA);
	      cs_spfree(C);
	      return;
	     
	    }
	    else if(strcmp(*method, "SPARSE SPD") == 0)
	    {
            S = cs_schol(1, C);
            N = cs_chol(C, S);
            fprintf(fp,"\n-------DC_Sweep:CHOLESKY SPARSE---------\n");
            for(val = dc_sweep->start_value; val <= dc_sweep->end_value; val = val + dc_sweep->inc)
            {
		
                for (i = 0; i < DIM; i++) {
                    RHS_copy[i] = RHS[i];
                    //printf("b[%d] = %#9.6f\n", i , RHS[i]);
                    x[i] = 0.0;
                }
                RHS_copy[Vindex_in_RHS] = val;
                cs_ipvec(S->pinv, RHS_copy, x, DIM);
                cs_lsolve(N->L, x);
                cs_ltsolve(N->L, x);
                cs_pvec(S->pinv, x, RHS_copy, DIM);
	  
                //fprintf(fp,"\n-------CHOLESKY_SPARSE-------\n");
                //fprintf(fp,"\n-------Vector: x_vector ------- \n");
              //  for(i=0; i<DIM; i++){
                    fprintf(fp,"Sweep source current at %#9.6f :   Node %s value: %.6f \n", val, dc_sweep->source, RHS_copy[Vindex_in_RHS]);
                    //printf("x[%d] = %#9.6f\n",i, RHS_copy[i]);
                //}
            }
	    }
	    else if( strcmp(*method, "SPARSE") == 0 )
	    {
            cs_sfree(S);
            cs_nfree(N);
            S = cs_sqr(2, C, 0);
            N = cs_lu(C, S, 1);
            
            fprintf(fp,"\n-------DC_Sweep:LU SPARSE---------\n");
            for(val = dc_sweep->start_value; val <= dc_sweep->end_value; val = val + dc_sweep->inc)
            {
                
                for (i = 0; i < DIM; i++) {
                    RHS_copy[i] = RHS[i];
                    //printf("b[%d] = %#9.6f\n", i , RHS[i]);
                    x[i] = 0.0;
                }
                RHS_copy[Vindex_in_RHS] = val;
                printf("%d %f\n", Vindex_in_RHS, val);;
                cs_ipvec(N->pinv, RHS_copy, x, DIM);
                cs_lsolve(N->L, x);
                cs_usolve(N->U, x);
                cs_ipvec(S->q, x, RHS_copy, DIM);
	
                /* the solution is stored in temp vector b */
                //fprintf(fp,"\n-------LU_SPARSE-------\n");
                //fprintf(fp,"\n-------Vector: x_vector ------- \n");
                //for(i=0; i<DIM; i++){
                fprintf(fp,"Sweep source current at %#9.6f :   Node %s value: %.6f \n", val, dc_sweep->source, RHS_copy[Vindex_in_RHS]);  // h lysh apothhkeyetai sto dianysma deksiou melous b
                //printf("x[%d] = %#9.6f\n",i, RHS_copy[i]);
                //}
            }
	    }
	    else if( strcmp(*method, "SPARSE ITER") == 0 )
	    {
            fprintf(fp,"\n-------DC_Sweep:BiCG SPARSE---------\n");
            for(val = dc_sweep->start_value; val <= dc_sweep->end_value; val = val + dc_sweep->inc)
            {
                
                for (i = 0; i < DIM; i++) {
                    RHS_copy[i] = RHS[i];
                    //printf("b[%d] = %#9.6f\n", i , RHS[i]);
                    x_vector[i] = 0.0;
                }
                RHS_copy[Vindex_in_RHS] = val;
                BiCG_solve_Sparse(C, RHS_copy, x_vector, DIM, itol);
                //fprintf(fp,"\n -------Vector: x_vector ------- \n");
                //for(i=0; i<DIM; i++)
                //{
                    fprintf(fp,"Sweep source current at %#9.6f :   Node %s value: %.6f \n", val, dc_sweep->source, x_vector[Vindex_in_RHS]);
                    //printf("x[%d] = %#9.6f\n", i, x_vector[i]);
                //}
            }
            cs_spfree(CSA);
            cs_spfree(C);
            return;
        }
	}
  
	cs_sfree(S);
	cs_nfree(N);
	cs_spfree(CSA);
	cs_spfree(C);
	return;
}

void computeMna(char **method, int *n, int *m2, double *itol, DC_sweep *dc_sweep, Plot *plot)
{
    printf("\n****** MNA with %s Method, %d N, %d M2, %f ITOL******\n", *method, *n, *m2, *itol);
    struct term_2 *look;
    int i, j, index, DIM, s;
	int pos, neg;
	int inc = 0;  // xrhsimopoieitai gia arithmish twn kladwn tou group 2.
	//double itol_new;
	
	int Vindex_in_RHS;   // thesh sto deksi melos opou vrisketai h timh ths phghs tashs pou prepei na tropopoioume se kathe
                             // vhma ths DC_sweep	
	
	double **A;	     // pinakas pou dhmiourgeitai kata thn kataskeysh tou MNA systhmatos (G)
	double *RHS;
	double *x;
	double val;
    double result;
	
    gsl_matrix *A_matrix;		// matrix A pou xrhsimopoieitai kai gia MNA
    gsl_vector *b_vector;		// vector b
    gsl_vector *x_vector;		//pinakas agnwstw
    
    /* Element for solving the system with the LU direct method */
    gsl_permutation *perm_vector;	/* permutation vector */
    
    FILE *fp;
	
    DIM=*n+*m2;
    //itol_new = itol;
    A_matrix = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_set_zero(A_matrix);
    
    
    b_vector = gsl_vector_alloc(DIM);
    gsl_vector_set_zero(b_vector);
    
    
    x_vector = gsl_vector_alloc(DIM);
    gsl_vector_set_zero(x_vector);
	
	
	
	
	/*A = (double **)malloc((DIM)*sizeof(double *));
	for(i=0; i<(DIM); i++) {
	  A[i]  = (double *)malloc((DIM)*sizeof(double));
	}
	RHS = (double *)malloc((DIM)*sizeof(double));
	x = (double *)malloc((DIM)*sizeof(double));
	
	for(i=0; i<(DIM); i++) {
	  for(j=0; j<(DIM); j++) {
	    A[i][j] = 0;  
	  }
	  RHS[i] = 0;
	  x[i] = 0;
	}*/
    A = (double **)calloc(DIM,sizeof(double*));
    for(i=0;i<DIM;i++){
        A[i] = (double*)calloc(DIM,sizeof(double));
    }
    RHS = (double*)calloc(DIM,sizeof(double));
    x = (double*)calloc(DIM,sizeof(double));
	/* kataskeyh tou MNA systhmatos */
	for(look=root_2; look->nxt != NULL; look=look->nxt) {

		pos = -1;
		neg = -1;
		
        if( (look->type == 'R') || (look->type == 'r') ) {
            pos = searchHashTable(look->pos_term);
            neg = searchHashTable(look->neg_term);

            if( (pos != -1) && (neg != -1) ) {
				A[pos][pos] += 1/look->value;
				A[pos][neg] -= 1/look->value;
				A[neg][pos] -= 1/look->value;
				A[neg][neg] += 1/look->value;
			}
			else if(pos != -1) {
				A[pos][pos] += 1/look->value;
			}
			else if(neg != -1) {
				A[neg][neg] += 1/look->value;
			}
		}

		else if( (look->type == 'I') || (look->type == 'i') ) {

		  
			pos = searchHashTable(look->pos_term);
			neg = searchHashTable(look->neg_term);
			
			if ( ((look->type == 'I') || (look->type == 'i')) && strcmp(dc_sweep->source+1, look->name) == 0 ) 
				 Vindex_in_RHS = pos;

			if( (pos != -1) && (neg != -1) ) {
				RHS[pos] -= look->value;
				RHS[neg] +=  look->value;
			}
			else if(pos != -1) {
				RHS[pos] -= look->value;
			}
			else if(neg != -1) {
				RHS[neg] +=  look->value;
			}
		}

		else if( (look->type == 'V') || (look->type == 'v') || (look->type == 'L') || (look->type == 'l') ) {

			/* Apothhkeyoume sto Vindex_in_RHS th thesh sto dianysma tou deksiou melous opou
			   apothhkeyetai h timh ths phgh tashs gia thn opoia tha kanoume DC_sweep analysh.
			   To kanoume ayto wste na mh xreiazetai na kanoume compute_MNA se kathe vhma ths 
			   DC_sweep, alla na allazoume mono th sygkekrimenh thesh sto deksi melos.
			*/
			if ( ((look->type == 'V') || (look->type == 'v')) && strcmp(dc_sweep->source+1, look->name) == 0 ) 
				 Vindex_in_RHS = (*n + inc);

			pos = searchHashTable(look->pos_term);
			neg = searchHashTable(look->neg_term);
		
			if( (pos != -1) && (neg != -1) ) {
				A[pos][*n+inc] = 1;
				A[inc+*n][pos] = 1;
				A[neg][*n+inc] = -1;
				A[inc+*n][neg] = -1;
			}
			else if(pos != -1) {
				A[pos][*n+inc] = 1;
				A[inc+*n][pos] = 1;
			}
			else if(neg != -1) {
				A[neg][*n+inc] = -1;
				A[inc+*n][neg] = -1;
			}

			if( (look->type == 'V') || (look->type == 'v') ) {
				RHS[*n+inc] = look->value;
			}

			inc++;
		}
	}

	//printf("%d %f\n",DIM, A[0][0]);
	fp=fopen("MNAmatrix.txt", "w");
	fprintf(fp," -------Matrix MNA-------\n");
	
	for(i=0; i<DIM; i++)
	{
	  for(j=0; j<DIM; j++)
	  {
	    fprintf(fp,"%.3f\t",A[i][j]);
	  }
	  fprintf(fp,"\n");
	}
	
	for(i=0; i<DIM; i++)
	{
	  for(j=0; j<DIM; j++)
	  {
            gsl_matrix_set(A_matrix, i, j, A[i][j]);
	  }
	}
	///To deksi melos dhladh to dianusma b
	for(i=0; i<DIM; i++)
	{
	  gsl_vector_set(b_vector, i, RHS[i]);
	}
	
	
	//printf("%s\n", *method);
	if( strcmp(*method, "SPD") == 0){
        
        printf("\n-------CHOLESKY METHOD------\n");
        gsl_linalg_cholesky_decomp(A_matrix);
        gsl_linalg_cholesky_solve(A_matrix,b_vector, x_vector);
        //fprintf(fp,"\n -------Matrix: CHOLESKY ------- \n");
        //for (i = 0; i < DIM; i++) {
          //  for (j = 0; j < DIM; j++) {
            //    fprintf(fp, "%f\t", gsl_matrix_get(A_matrix, i, j));
           // }
           // fprintf(fp, "\n");
        //}
        //fprintf(fp," -------Matrix: CHOLESKY ------- \n\n");
	  
        fprintf(fp,"\n -------Vector: x_vector ------- \n");
        gsl_vector_fprintf(fp, x_vector, "%f");;
	  
	}
	else if( strcmp(*method, "ITER") == 0){
        
        printf("\n------- BI-CG METHOD -------\n");
        BiCG_solve(A_matrix, b_vector, x_vector, DIM, itol);
        fprintf(fp,"\n -------Vector: x_vector ------- \n");
        gsl_vector_fprintf(fp, x_vector, "%f");
	}
	else if( strcmp(*method, "LU") == 0 )
	{
	  printf("\n -------LU METHOD------ \n");
	  perm_vector = gsl_permutation_calloc(DIM);
	  
	 // fprintf(fp, "\nperm_old: \n");
	 // gsl_permutation_fprintf (fp, perm_vector, " %2u");
	  gsl_linalg_LU_decomp(A_matrix, perm_vector, &s);
	  gsl_linalg_LU_solve(A_matrix, perm_vector, b_vector, x_vector);
	  fprintf(fp, "\npermutation vector: \n");
	  gsl_permutation_fprintf (fp, perm_vector, " %2u");
	  //fprintf(fp,"\n -------Matrix: LU ------- \n");
	  /*for (i = 0; i < DIM; i++) {
		for (j = 0; j < DIM; j++) {
			fprintf(fp, "%f\t", gsl_matrix_get(A_matrix, i, j));
		}
		fprintf(fp, "\n");
	  }*/
	  //fprintf(fp," -------Matrix: LU ------- \n\n");
	  
	  fprintf(fp,"\n -------Vector: x_vector ------- \n");
	  gsl_vector_fprintf(fp, x_vector, "%f");
	  
	}
	else if( strcmp(*method, "SPD ITER") == 0)
	{
	  printf("\n-------CG------\n");
	  CG_solve(A_matrix, b_vector, x_vector, DIM, itol);
	  fprintf(fp,"\n -------Vector: x_vector ------- \n");
	  gsl_vector_fprintf(fp, x_vector, "%f");
	
	}
  
    printf("\nDCSWEEP : %d flag, %d id, %s Source, %f Step, %f Start Value\n\n", dc_sweep->DC_flag, Vindex_in_RHS, dc_sweep->source, dc_sweep->inc, dc_sweep->start_value);
	
    if(dc_sweep->DC_flag == 1)
	{
        for(val = dc_sweep->start_value; val <= dc_sweep->end_value; val = val + dc_sweep->inc)
        {
            if(strcmp(*method, "SPD ITER") == 0)
            {
                fprintf(fp,"\n-------DC_Sweep:CG---------\n");
                gsl_vector_set(b_vector, Vindex_in_RHS, val);
                CG_solve(A_matrix, b_vector, x_vector, DIM, itol);
                //result = gsl_vector_get(b_vector, Vindex_in_RHS);
                fprintf(fp,"-------Vector: b_vector ------- \n");
                gsl_vector_fprintf(fp, b_vector, "%f");
                fprintf(fp,"-------Vector: x_vector ------- \n");
                gsl_vector_fprintf(fp, x_vector, "%f");
                //fprintf(fp,"SWEEP SOURCE CURRENT AT %#9.6f : NODE %s value : %.6f",val, dc_sweep->source, result);
            }
            else if(strcmp(*method, "SPD") == 0)
            {
                fprintf(fp,"\n-------DC_Sweep:CHOLESKY---------\n");
                gsl_vector_set(b_vector, Vindex_in_RHS, val);
                gsl_linalg_cholesky_solve(A_matrix,b_vector, x_vector);
                fprintf(fp,"-------Vector: b_vector ------- \n");
                gsl_vector_fprintf(fp, b_vector, "%f");
                fprintf(fp,"-------Vector: x_vector ------- \n");
                gsl_vector_fprintf(fp, x_vector, "%f");
            }
            else if( strcmp(*method, "LU") == 0 )
            {
                fprintf(fp,"\n-------DC_Sweep:LU---------\n");
                gsl_vector_set(b_vector, Vindex_in_RHS, val);
                gsl_linalg_LU_solve(A_matrix, perm_vector, b_vector, x_vector);
                fprintf(fp,"-------Vector: b_vector ------- \n");
                gsl_vector_fprintf(fp, b_vector, "%f");
                fprintf(fp,"-------Vector: x_vector ------- \n");
                gsl_vector_fprintf(fp, x_vector, "%f");
            }
            else if( strcmp(*method, "ITER") == 0 )
            {
                fprintf(fp,"\n-------DC_Sweep:BiCG---------\n");
                gsl_vector_set(b_vector, Vindex_in_RHS, val);
                BiCG_solve(A_matrix, b_vector, x_vector, DIM, itol);
                fprintf(fp,"-------Vector: b_vector ------- \n");
                gsl_vector_fprintf(fp, b_vector, "%f");
                fprintf(fp,"-------Vector: x_vector ------- \n");
                gsl_vector_fprintf(fp, x_vector, "%f");
            }
        }
    }
  free(A);
  free(RHS);
  gsl_matrix_free(A_matrix);
  gsl_vector_free(b_vector);
  gsl_vector_free(x_vector);
  if( strcmp(*method, "LU") == 0)
  {
    gsl_permutation_free(perm_vector);
  }
  fclose(fp);
}

int BiCG_solve(gsl_matrix *A_matrix, gsl_vector *b_vector, gsl_vector *x_vector, int DIM, double *itol_n) {
    int i, iter, MAX_ITER, ret = 1;
	double norm_r, norm_b;  /* Euclidean norms of vectors_r and vector_b */
	double rho, rho_1=1, alpha, beta, omega, temp;
	//const int DIM = node_cnt+volt_src_cnt; /* dimension of matrices/vectors */
	const double EPS = 1e-14;  /* zero value for doubles */
	double itol;
	gsl_vector *A_mul_x = gsl_vector_calloc(DIM); /* A*x result */
	
	gsl_vector *M = gsl_vector_calloc(DIM); /* preconditioner */
	gsl_vector *r_vector = gsl_vector_calloc(DIM); /* residual */
	gsl_vector *z_vector = gsl_vector_calloc(DIM);
	gsl_vector *p_vector = gsl_vector_calloc(DIM); /* search dir. */
	gsl_vector *q_vector = gsl_vector_calloc(DIM);
	gsl_vector *r_tilde_vector = gsl_vector_calloc(DIM); /* 2nd residual */
	gsl_vector *z_tilde_vector = gsl_vector_calloc(DIM);
	gsl_vector *p_tilde_vector = gsl_vector_calloc(DIM); /* 2nd search dir.*/
	gsl_vector *q_tilde_vector = gsl_vector_calloc(DIM);

	itol = *itol_n;
	/* set max number of iterations */
	MAX_ITER = DIM*DIM; 
	
	/* x_vector has already been initialized to zero-values */

	/* create the preconditioner */
	precond_create(A_matrix, M, DIM);
	
	/* A dot x */
	gsl_blas_dgemv(CblasNoTrans, 1.0, A_matrix, x_vector, 0.0, A_mul_x);
	
	/* compute initial residuals */
	for (i = 0; i < DIM; i++) {
		temp = gsl_vector_get(b_vector, i) - gsl_vector_get(A_mul_x, i);
		gsl_vector_set(r_vector, i, temp);
		gsl_vector_set(r_tilde_vector, i, temp);
	}
	gsl_vector_free(A_mul_x);
	
	/* compute norms */
	norm_r = gsl_blas_dnrm2(r_vector);
	norm_b = gsl_blas_dnrm2(b_vector);
	
	if (norm_b == 0) { 
		norm_b = 1;
	}
	
	iter = 0;
	while ((norm_r/norm_b > itol) && (iter < MAX_ITER)) {
		iter++;
		
		/* solve M * z = r */
		precond_solve(z_vector, r_vector, M, DIM);
		
		/* solve M^T* z~ = r~ */
		precond_solve(z_tilde_vector, r_tilde_vector, M, DIM);

		/* rho = z dot r~*/
		gsl_blas_ddot(z_vector, r_tilde_vector, &rho);

		if (fabs(rho) < EPS) {
			/* failure */
			printf("Bi-CG failed!\n");
			ret = 0;
			break;
		}
		if (iter == 1) {
			for(i = 0; i< DIM; i++) {
				gsl_vector_set(p_vector, i, gsl_vector_get(z_vector, i));
				gsl_vector_set(p_tilde_vector, i, gsl_vector_get(z_tilde_vector, i));
			}
		}
		else {
			beta = rho / rho_1;
			/* p = z  + beta * p */
			for (i = 0; i < DIM; i++) {
				temp = gsl_vector_get(p_vector, i);
				gsl_vector_set(p_vector, i, gsl_vector_get(z_vector, i) + beta*temp);
			}
			
			/* p~ = z~  + beta * p~ */
			for (i = 0; i < DIM; i++) {
				temp = gsl_vector_get(p_tilde_vector, i);
				gsl_vector_set(p_tilde_vector, i, gsl_vector_get(z_tilde_vector, i) + beta*temp);
			}
		}
		rho_1 = rho;
		
		/* q = A * p */
		mv_prod_calc(q_vector, A_matrix, p_vector, 0, DIM);
		
		/* q~ = A^T * p~  */
		mv_prod_calc(q_tilde_vector, A_matrix, p_tilde_vector, 1, DIM);
		
		/* omega = p~ dot q */
		gsl_blas_ddot(p_tilde_vector, q_vector, &omega);
		
		if (fabs(omega) < EPS) {
			/* failure */
			printf("Bi-CG failed!\n");
			ret = 0;
			break;
		} 
		
		alpha = rho / omega;
		
		/* x = x + alpha * p */
		for (i = 0; i < DIM; i++) {
			temp = gsl_vector_get(x_vector, i);
			gsl_vector_set(x_vector, i, temp + alpha * gsl_vector_get(p_vector, i));
		}
		/* r = r - alpha * q */
		for (i = 0; i < DIM; i++) {
			temp = gsl_vector_get(r_vector, i);
			gsl_vector_set(r_vector, i,  temp - alpha * gsl_vector_get(q_vector, i));
		}
		
		/* r~ = r~ - alpha * q~ */
		for (i = 0; i < DIM; i++) {
			temp = gsl_vector_get(r_tilde_vector, i);
			gsl_vector_set(r_tilde_vector, i,  
			temp - alpha * gsl_vector_get(q_tilde_vector, i));
		}
		
		
		
		/* re-compute norm */
		norm_r = gsl_blas_dnrm2(r_vector);

	}
	
	/* free the allocated memory */
	gsl_vector_free(M);
	gsl_vector_free(r_vector);
	gsl_vector_free(z_vector);
	gsl_vector_free(p_vector);
	gsl_vector_free(q_vector);
	gsl_vector_free(r_tilde_vector);
	gsl_vector_free(z_tilde_vector);
	gsl_vector_free(p_tilde_vector);
	gsl_vector_free(q_tilde_vector);
	 
	return (ret);
}

void precond_create(gsl_matrix *A_matrix, gsl_vector *M, int DIM) {
	int i;
	double a_ii;
	
	for (i = 0; i < DIM; i++) {
		a_ii = gsl_matrix_get(A_matrix, i, i);
		if (a_ii == 0) { a_ii = 1; }
		gsl_vector_set(M, i, a_ii);
	}
// 	print_gsl_vector("M", M);
	
	return;
}

void precond_solve(gsl_vector *z, const gsl_vector *r, const gsl_vector *M, int DIM) {
	int i;
	
	for (i = 0; i < DIM; i++) {
		gsl_vector_set(z, i, gsl_vector_get(r, i) / gsl_vector_get(M, i));
	}
	
	return;
}

double norm_gslcalc(const gsl_vector *x, int DIM) {
	int i;
	double sum = 0, temp;
	
	for (i = 0; i < DIM; i++) { 
		temp = gsl_vector_get(x, i);
		sum += temp*temp;
	}
	return (sqrt(sum));
}

double vv_dot_gslcalc(const gsl_vector *x, const gsl_vector *y, int DIM) {
	int i;
	double sum = 0;
	
	for (i = 0; i < DIM; i++) { 
			sum += gsl_vector_get(x, i) * gsl_vector_get(y, i);
	}
	
	return (sum);
}

void mv_prod_calc(gsl_vector *y, const gsl_matrix *A, const gsl_vector *x, int transpose, int DIM) {
	int i, j;
	double sum;

	if (!transpose) {
		/* y = A * x  */
		for (i = 0; i < DIM; i++) {
			sum = 0;
			for (j = 0;  j < DIM; j++) {
				sum += gsl_matrix_get(A, i, j)*gsl_vector_get(x, j);
			}
			gsl_vector_set(y, i, sum);
		}
// 		gsl_blas_dgemv(CblasNoTrans, 1.0, A, x, 0.0, y);
	}
	else {
		/* y = A^T * x  */
		for (i = 0; i < DIM; i++) {
			sum = 0;
			for (j = 0;  j < DIM; j++) {
				sum += gsl_matrix_get(A, j, i) * gsl_vector_get(x, j);
			}
			gsl_vector_set(y, i, sum);
		}
	}
	return;
	
}


int CG_solve(gsl_matrix *A_matrix, gsl_vector * b_vector,gsl_vector *x_vector, int DIM, double *itol_n)
{
  int i, iter, MAX_ITER;
	double norm_r, norm_b;	/* Euclidean norms of r_vector and b_vector */
	double rho, rho_1=1, alpha, beta, p_dot_q, temp;
	double itol;
	gsl_vector *A_mul_x = gsl_vector_calloc(DIM); /* A * x result */
	
	gsl_vector *M = gsl_vector_calloc(DIM); /* preconditioner */
	gsl_vector *r_vector = gsl_vector_calloc(DIM); /* residual */
	gsl_vector *z_vector = gsl_vector_calloc(DIM);
	gsl_vector *p_vector = gsl_vector_calloc(DIM); /* search dir. */
	gsl_vector *q_vector = gsl_vector_calloc(DIM);
	
	/* set max number of iterations */
	MAX_ITER = DIM;
	itol = *itol_n;
	/* x_vector has already been initialized to zero-values */
	
	/* create the preconditioner */
	precond_create(A_matrix, M, DIM);
	
	/* A dot x */
	gsl_blas_dgemv(CblasNoTrans, 1.0, A_matrix, x_vector, 0.0, A_mul_x);
	
	/* compute initial residual */
	for (i = 0; i < DIM; i++) {
		temp = gsl_vector_get(b_vector, i) - gsl_vector_get(A_mul_x, i);
		gsl_vector_set(r_vector, i, temp);
	}
	gsl_vector_free(A_mul_x);

	/* compute norms */
	norm_r = gsl_blas_dnrm2(r_vector);
	norm_b = gsl_blas_dnrm2(b_vector);
// 	norm_r = norm_gslcalc(r_vector);
// 	norm_b = norm_gslcalc(b_vector);
	
	if (norm_b == 0) { 
		norm_b = 1;
	}
	
	iter = 0;
	while ((norm_r/norm_b > itol) && (iter < MAX_ITER)) {
		iter++;
		
		/* solve M * z = r */
		precond_solve(z_vector, r_vector, M, DIM);
		
		/* rho = r dot z */
		gsl_blas_ddot(r_vector, z_vector, &rho);
// 		rho = vv_dot_gslcalc(r_vector, z_vector);
		
		if (iter == 1) {
			/* p = z */
			for(i = 0; i< DIM; i++) {
				gsl_vector_set(p_vector, i, gsl_vector_get(z_vector, i));
			}
		}
		else {
			beta = rho / rho_1;
			/* p = z  + beta * p */
			for (i = 0; i < DIM; i++) {
				temp = gsl_vector_get(p_vector, i);
				gsl_vector_set(p_vector, i, gsl_vector_get(z_vector, i) + beta*temp);
			}
		}
		rho_1 = rho;
		
		/* q = A * p */;
		mv_prod_calc(q_vector, A_matrix, p_vector, 0, DIM);
		
		/* p dot q */
		gsl_blas_ddot(p_vector, q_vector, &p_dot_q);
// 		p_dot_q = vv_dot_gslcalc(p_vector, q_vector);
		
		alpha = rho / p_dot_q; 
		
		/* x = x + alpha * p */
		for (i = 0; i < DIM; i++) {
			temp = gsl_vector_get(x_vector, i);
			gsl_vector_set(x_vector, i, temp + alpha*gsl_vector_get(p_vector, i));
		}
		/* r = r - alpha * q */
		for (i = 0; i < DIM; i++) {
			temp = gsl_vector_get(r_vector, i);
			gsl_vector_set(r_vector, i,  temp - alpha*gsl_vector_get(q_vector, i));
		}
		
		/* re-compute norm */
		norm_r = gsl_blas_dnrm2(r_vector);
// 		norm_r = norm_gslcalc(r_vector);
	}
	
	/* free the allocated memory */
	gsl_vector_free(M);
	gsl_vector_free(r_vector);
	gsl_vector_free(z_vector);
	gsl_vector_free(p_vector);
	gsl_vector_free(q_vector);
	
	return (1);
}


int BiCG_solve_Sparse(cs *C, double *RHS, double *x_vect, int DIM, double *itol_n) {
	int i, iter, MAX_ITER, ret = 1;
	double norm_r, norm_b, itol;	/* Euclidean norms of vectors_r and vector_b */
	double rho, rho_1=1, alpha, beta, omega;
	const double EPS = 1e-14; 	/* zero value for doubles */
	double *A_mul_x;	/* A * x result */	
	double *M;		/* preconditioner, the diagonal/Jacobian is used */
	double *r_vec;	/* residual */
	double *z_vec;	
	double *p_vec;	/* search dir. */
	double *q_vec;
	double *r_tilde_vec;	/* 2nd residual */
	double *z_tilde_vec;	
	double *p_tilde_vec;	/* 2nd search dir.*/
	double *q_tilde_vec;
	
	
	M = (double*)calloc(DIM, sizeof(*M));
	if (M == NULL) { 
	  printf("Memory Error\n");
	  exit(0); 
	}
	
	r_vec = (double*)calloc(DIM, sizeof(*r_vec));
	if (r_vec == NULL) { 
	  printf("Memory Error\n");
	  exit(0);
	}
	
	z_vec = (double*)calloc(DIM, sizeof(*z_vec));
	if (z_vec == NULL) { 
	  printf("Memory Error\n");
	  exit(0);
	}
	
	p_vec = (double*)calloc(DIM, sizeof(*p_vec));
	if (p_vec == NULL) { 
	  printf("Memory Error\n");
	  exit(0);
	}
	
	q_vec = (double*)calloc(DIM, sizeof(*q_vec));
	if (q_vec == NULL) { 
	  printf("Memory Error\n");
	  exit(0);
	}
	
	r_tilde_vec = (double*)calloc(DIM, sizeof(*r_tilde_vec));
	if (r_tilde_vec == NULL) { 
	  printf("Memory Error\n");
	  exit(0);
	}
	
	z_tilde_vec = (double*)calloc(DIM, sizeof(*z_tilde_vec));
	if (z_tilde_vec == NULL) { 
	  printf("Memory Error\n");
	  exit(0);
	}
	
	p_tilde_vec = (double*)calloc(DIM, sizeof(*p_tilde_vec));
	if (p_tilde_vec == NULL) { 
	  printf("Memory Error\n");
	  exit(0);
	}
	
	q_tilde_vec = (double*)calloc(DIM, sizeof(*q_tilde_vec));
	if (q_tilde_vec == NULL) { 
	  printf("Memory Error\n");
	  exit(0);
	}
	
	/* set max number of iterations */
	MAX_ITER = DIM*DIM;
	itol = *itol_n;
	/* x_vector has already been initialized to zero-values */

	/* create the preconditioner M */
	precond_spcreate(C, M, DIM);
	
	/* A_matrix dot x_vector */
	A_mul_x = (double*)calloc(DIM, sizeof(*q_tilde_vec));
	if (q_tilde_vec == NULL) { 
	  printf("Memory Error\n");
	  exit(0);
	}
	cs_gaxpy(C, x_vect, A_mul_x);
	
	/* compute initial residuals */
	for (i = 0; i < DIM; i++) {
		r_vec[i] = RHS[i] - A_mul_x[i];
		r_tilde_vec[i] = RHS[i] - A_mul_x[i];
	}
	free(A_mul_x);
	
	/* compute norms */
	norm_r = norm_calc(r_vec, DIM);
	norm_b = norm_calc(RHS, DIM);
	
	if (norm_b == 0) { 
		norm_b = 1;
	}

	iter = 0;
	while ((norm_r/norm_b > itol) && (iter < MAX_ITER)) {
		iter++;
		
		/* solve M * z = r */
		precond_spsolve(z_vec, r_vec, M, DIM);
		
		/* solve M^T* z~ = r~ */
		precond_spsolve(z_tilde_vec, r_tilde_vec, M, DIM);

		/* rho = z dot r~*/
		rho = vv_dot_calc(r_tilde_vec, z_vec, DIM);
		
		if (fabs(rho) < EPS) {
			/* failure */
			printf("Bi-CG failed!\n");
			exit(0);
		}
		if (iter == 1) {
			for(i = 0; i < DIM; i++) {
				p_vec[i] = z_vec[i];
				p_tilde_vec[i] = z_tilde_vec[i];
			}
		}
		else {
			beta = rho / rho_1;
			/* p = z  + beta * p  and  p~ = z~  + beta * p~ */
			for (i = 0; i < DIM; i++) {
				p_vec[i] = z_vec[i] + beta * p_vec[i];
				p_tilde_vec[i] = z_tilde_vec[i] + beta * p_tilde_vec[i];
			}
		}
		rho_1 = rho;
		
		/* q = A * p */
		mv_spprod_calc(q_vec, C, C, p_vec, 0, DIM);
		
		/* q~ = A^T * p~ */
		mv_spprod_calc(q_tilde_vec, C, C, p_tilde_vec, 1, DIM);
		
		/* omega = p~ dot q */
		omega = vv_dot_calc(p_tilde_vec, q_vec, DIM);
		if (fabs(omega) < EPS) {
			/* failure */
			printf("Bi-CG failed!\n");
			exit(0);
		} 
		
		alpha = rho / omega;
		
		/* x = x + alpha * p */
		for (i = 0; i < DIM; i++) {
			x_vect[i] +=  alpha * p_vec[i];
		}
		/* r = r - alpha * q */
		for (i = 0; i < DIM; i++) {
			r_vec[i] -=  alpha * q_vec[i];
		}
		
		/* r~ = r~ - alpha * q~ */
		for (i = 0; i < DIM; i++) {
			r_tilde_vec[i] -=  alpha * q_tilde_vec[i];
		}
		
		
		/* re-compute norm */
		norm_r = norm_calc(r_vec, DIM);
	}
	
	free(M);
	free(r_vec);
	free(z_vec);
	free(p_vec);
	free(q_vec);
	free(r_tilde_vec);
	free(z_tilde_vec);
	free(p_tilde_vec);
	free(q_tilde_vec);
	
	return (ret);
}

void precond_spcreate(cs *C, double *M, int DIM) {
	int p, j, n, *Ap, *Ai;
	double a_ii;
	double *Ax;
	
	n = C->n;
	Ap = C->p;
	Ai = C->i;
	Ax = C->x;
	
	for (j = 0; j < n; j++) {
		a_ii = 1;
		for (p = Ap[j]; p < Ap[j + 1]; p++) {
			if (j == Ai[p]) { a_ii = Ax[p]; }
		}
		M[j] = a_ii;
	}
// 	print_vector("M", M);
	
	return;
}

void precond_spsolve(double *z, const double *r, const double *M, int DIM) {
	int i;
	
	for (i = 0; i < DIM; i++) {
		z[i] = r[i] / M[i];
	}
	
	return;
}

double norm_calc(const double *x, int DIM) {
	int i;
	double sum = 0;
	
	for (i = 0; i < DIM; i++) { 
			sum += x[i] * x[i];
	}
	return (sqrt(sum));
}

double vv_dot_calc(const double *x, const double *y, int DIM) {
	int i;
	double sum = 0;
	
	for (i = 0; i < DIM; i++) { 
			sum += x[i] * y[i];
	}
	
	return (sum);
}

void mv_spprod_calc(double *y, const cs *A, cs *C, const double *x, int transpose, int DIM) {
	int j, p;
	
	if (!transpose) {
		/* y = A * x  */
		for (j = 0; j < DIM; j++) {
			y[j] = 0.0;
		}
		
		/* returns 1 or 0 */
		//cs_gaxpy(A, p_temp, q_temp);
		for (j = 0; j < DIM; j++) {
			for (p = C->p[j]; p < C->p[j+1]; p++) {
				y[A->i[p]] = y[A->i[p]] + A->x[p] * x[j];
			}
		}
	}
	else {
		/* y = A^T * x  */
		for (j = 0; j < DIM; j++) {
			y[j] = 0.0;
			for (p = A->p[j]; p < A->p[j+1]; p++) {
				y[j] = y[j] + A->x[p] * x[A->i[p]];
			}
		}
	}
	return;
}

int CG_solve_Sparse(cs *C, double *RHS, double *x_vect, int DIM, double *itol_n) {
	int i, iter, MAX_ITER;
	double norm_r, norm_b, itol;	/* Euclidean norms of r_vector and b_vector */
	double rho, rho_1=1, alpha, beta, p_dot_q;
	double *A_mul_x;	/* A * x result */
	
	double *M;		/* preconditioner, the diagonal/Jacobian is used */
	double *r_vec;	/* residual */
	double *z_vec;	
	double *p_vec;	/* search dir. */
	double *q_vec;
	
	M = (double*)calloc(DIM, sizeof(*M));
	if (M == NULL) { 
	  printf("Memory Error\n");
	  exit(0);
	}
	
	r_vec = (double*)calloc(DIM, sizeof(*r_vec));
	if (r_vec == NULL) { 
	  printf("Memory Error\n");
	  exit(0);
	}
	
	z_vec = (double*)calloc(DIM, sizeof(*z_vec));
	if (z_vec == NULL) { 
	  printf("Memory Error\n");
	  exit(0);
	}
	
	p_vec = (double*)calloc(DIM, sizeof(*p_vec));
	if (p_vec == NULL) { 
	  printf("Memory Error\n");
	  exit(0);
	}
	
	q_vec = (double*)calloc(DIM, sizeof(*q_vec));
	if (q_vec == NULL) { 
	  printf("Memory Error\n");
	  exit(0); 
	}
	
	
	/* set max number of iterations */
	MAX_ITER = DIM;
	itol = *itol_n;
	/* x_vector has already been initialized to zero-values */
	
	/* create the preconditioner M */
	precond_spcreate(C, M, DIM);
	
	/* A_matrix dot x_vector */
	A_mul_x = malloc(DIM*sizeof(*A_mul_x));
	cs_gaxpy(C, x_vect, A_mul_x);
	
	/* compute initial residual */
	for (i = 0; i < DIM; i++) {
		r_vec[i] = RHS[i] - A_mul_x[i];
	}
	free(A_mul_x);
	
	/* compute norms */
	norm_r = norm_calc(r_vec, DIM);
	norm_b = norm_calc(RHS, DIM);
	
	if (norm_b == 0) { 
		norm_b = 1;
	}
	
	iter = 0;
	while ((norm_r/norm_b > itol) && (iter < MAX_ITER)) {
		iter++;
		
		/* solve M * z = r */
		precond_spsolve(z_vec, r_vec, M, DIM);
		
		/* rho = r dot z */
		rho = vv_dot_calc(r_vec, z_vec, DIM);
		
		if (iter == 1) {
			for (i = 0; i < DIM; i++) {
				p_vec[i] = z_vec[i];
			}
		}
		else {
			beta = rho / rho_1;
			/* p = z  + beta * p */
			for (i = 0; i < DIM; i++) {
				p_vec[i] = z_vec[i] + beta * p_vec[i];
			}
		}
		rho_1 = rho;
		
		/* q = A * p */
		mv_spprod_calc(q_vec, C, C, p_vec, 0, DIM);
// 		cs_gaxpy(A_comp, p_vec, q_vec);

		/* p dot q */
		p_dot_q = vv_dot_calc(p_vec, q_vec, DIM);

		alpha = rho / p_dot_q; 
		
		/* x = x + alpha * p */
		for (i = 0; i < DIM; i++) {
			x_vect[i] += alpha * p_vec[i];
		}
		/* r = r - alpha * q */
		for (i = 0; i < DIM; i++) {
			r_vec[i] -= alpha * q_vec[i];
		}	
		
		/* re-compute norm */
		norm_r = norm_calc(r_vec, DIM);
	}
	
// 	print_vector("x (operating point)", x_vec);
	free(M);
	free(r_vec);
	free(z_vec);
	free(p_vec);
	free(q_vec);
	
	return (1);
}