///etc/resolv.conf


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utilities.h"
#include "hash_table.h"


/******************************************
 MAIN FUNCTION
 *******************************************/
int main(int argc, char * argv[])
{
    int i;
    int n=0;
    int m2= 0;
    int choice;
    int *pointern=&n;
    int *pointerm=&m2;
    //int *nz=&nz1;
    char *method = strdup("LU");
    double itol = 1e-3;  // itol gia methodous CG, Bi-CG
    int nz = 0;	     // arithmos non-zeros
    int cnz = 0;	     // arithmos non-zeros pinaka C
    //struct hashTable *runner;
    //runner=hashEntry;
    DC_sweep dc_sweep;
    Plot plot;
    //printf("to runner einai %s", runner->nodeName);
    //initHashTable();
    init_hash();
    
    dc_sweep.DC_flag     = 0;
    dc_sweep.source      = "";
    dc_sweep.start_value = 0;
    dc_sweep.end_value   = 0;
    dc_sweep.inc         = 0;
    
    plot.plot_flag = 0;
    //anoigoume to arxeio pou periexei to kuklwma
    //FILE * input_file_ptr;
    //input_file_ptr=fopen( argv[1], "r" );
    //printf( "**** opening file %s\n",argv[1] );
    
    //if( !input_file_ptr )
    //{
      //  printf( "error: can't find the specified file\n " );
    //}
    
    choice = read_netlist(argv[1], &n, &m2, &method, &itol,  &nz, &cnz, &dc_sweep, &plot);
    //printf("method is %s and choice is %d %d\n", method, choice, n+m2);
    //exit(0);
    
    if(choice == 1 || choice == 2 || choice == 3 || choice == 4)
    {
      computeSparseMna(&method, &n, &m2, &itol, &nz, &cnz, &dc_sweep, &plot);
    }
    else if(choice == 0 || choice == 5|| choice == 6 || choice == 7)
    {
      computeMna(&method, &n, &m2, &itol, &dc_sweep, &plot);
    }
    else{
      printf("ERROR 1\n");
    }
    
    freeHashTable();
    
    return 0;
}