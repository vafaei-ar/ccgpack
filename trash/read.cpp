#include <stdio.h>
#include <stdlib.h>
# include <fstream>
# include <iostream>

using namespace std;

int main(void) {
    float x;
		float data[1024][1024];


    FILE *file;  /* declare a FILE pointer  */
    file = fopen("1", "r");  /* open a text file for reading */

		int n=1024,i,j;


		for (int i = 1; i < n; i++) {
		 for (int j = 1; j < n; j++) {
				fscanf(file, "%f", &data[i][j]);
		}}

    fclose(file);
    return 0;
}
