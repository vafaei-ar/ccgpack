# include <stdlib.h>
# include <stdio.h>
# include <math.h>
//# include "Array.h"

using namespace std;
//using namespace Array;

double** reader(int rows, int columns, double** field){
    int i,j;

    double** dfield = new double*[rows];
    for (i = 0; i < rows; i++){
      dfield[i] = new double[columns]; 
      for (j = 0; j < columns; j++){dfield[i][j] =0;}}

    return dfield;
