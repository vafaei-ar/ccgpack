# include "fdct_wrapping.hpp"
# include "fdct_wrapping_inline.hpp"
# include <stdio.h>
# include <fstream>
# include <iostream>
# include <cstdio>

using namespace std;
using namespace fdct_wrapping_ns;

#ifdef __cplusplus
extern "C" {  // only need to export C interface if
              // used by C++ source code
#endif
void curvelet(double* array, int m, int n, int nbscales, int sreq, int nbangles_coarse, int ac ) {
    int i,j;
  clock_t ck0, ck1;
		
		CpxNumMat x(m,n);

    for(i=0; i<m; i++){
        for(j=0; j<n; j++){
        x(i,j) = array[i*n + j];
    }}

//  ck0 = clock();
  //fdct_wrapping_
  vector< vector<CpxNumMat> > c;  //vector<int> extra;
  fdct_wrapping(m, n, nbscales, nbangles_coarse, ac, x, c);

//  ck1 = clock();  cout<<"FDCT_WRAPPING_  takes "<<double(ck1-ck0)/CLOCKS_PER_SEC<<" seconds"<<endl;  
//ck0 = ck1;

  vector< vector<double> > sx, sy;
  vector< vector<double> > fx, fy;
  vector< vector<int> > nx, ny;
  fdct_wrapping_param(m, n, nbscales, nbangles_coarse, ac, sx, sy, fx, fy, nx, ny);

int s,w,nbangles;
  for(s=0; s<nbscales; s++){
		nbangles = nbangles_coarse * pow(2,(int)(s/2.));
		if (s==0 or s==nbscales-1) nbangles=1;
		if (s==sreq) nbangles=0;
  	for(w=0; w<nbangles; w++){
			for (i = 0; i < nx[s][w]; i++) {
		 		for (j = 0; j < ny[s][w]; j++) {
					c[s][w](i,j)=(0,0);
}}}}

	clear(x);
  ifdct_wrapping(m, n, nbscales, nbangles_coarse, ac, c, x);

//  ck1 = clock();  cout<<"IFDCT_WRAPPING_ takes "<<double(ck1-ck0)/CLOCKS_PER_SEC<<" seconds"<<endl;  
//ck0 = ck1;

    for(i=0; i<m; i++){
        for(j=0; j<n; j++){
        array[i*n + j] = real(x(i,j));
    }}
  

}
#ifdef __cplusplus
}
#endif


//#ifdef __cplusplus
//extern "C" {  // only need to export C interface if
//              // used by C++ source code
//#endif

//int prnt(double* ptr, int ny, int nx)
//{
//    int i, j;
//    for(i=0; i<ny; i++)
//    {
//        for(j=0; j<nx; j++)
//        {
//            printf("%.3f \t", ptr[i*nx + j]);
//        }
//        printf("\n");
//    }

//    return 0;
//}

//#ifdef __cplusplus
//}
//#endif

//extern "C" double **cur(double arr[1024][1024])
//{
//int i,j;

//for (i = 0; i < 1024; i++) {
//	for (j = 0; j < 1024; j++) {
//		arr[i][j]+=1;
//}}

//return arr;
//}

//extern "C" int square(int x)
//{
//  return x*x;
//}


