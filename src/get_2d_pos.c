#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>

#define LIMIT_REF_POINTS 20
#define LIMIT_SAMPLE 20

double REF_PixcelPos[LIMIT_REF_POINTS][2];
double REF_RealPos[LIMIT_REF_POINTS][2];

int NumOfSample;
int NumOfRef;


double PixcelPos[LIMIT_SAMPLE][2];
double RealPos[LIMIT_SAMPLE][2];


#define NUM_OF_DLT_PARAM 8 // DLT変数の数


double DLT_x[NUM_OF_DLT_PARAM]; // DLT Parameters
double DLT_b[LIMIT_REF_POINTS*2];
double DLT_A[LIMIT_REF_POINTS*2][NUM_OF_DLT_PARAM];
double DLT_AtA[NUM_OF_DLT_PARAM][NUM_OF_DLT_PARAM];




void gauss_solve(double A[NUM_OF_DLT_PARAM][NUM_OF_DLT_PARAM],
		 double b[NUM_OF_DLT_PARAM],
		 int num)
{
  int i, j, k;
  double sigma;
    
  /* forward elimination */
  for(k=0 ; k<num-1 ; k++){
    for(i=k+1 ; i<num ; i++){
      for(j=k+1 ; j<num ; j++){
	      A[i][j] = A[i][j]-A[k][j]*A[i][k]/A[k][k];
      }
      b[i] = b[i]-A[i][k]*b[k]/A[k][k];
    }
  }

  /* back substitution */
  b[num-1] = b[num-1] / A[num-1][num-1];
  for(i=num-2 ; i>=0 ; i--){
    sigma = 0.0;
    for(j=i+1 ; j<num ; j++){
      sigma += A[i][j]*b[j];
    }
    b[i] = 1/A[i][i]*(b[i]-sigma);
  }
}


void calc_DLT_param()
{
  int i, j, k;
  char str[256];
  FILE *fp;
  double k1, k2, b1, b2;

  // Read Pixcel Reference Points
  fp = fopen("image-ref-points.csv", "r");
  if(fp == NULL){
    fprintf(stderr, "Can't Open File  \"image-ref-points.csv\"  \n");
    fflush(stderr);
    exit(1);
  }

  int ret;
  int nnn = 0;
  while((ret = fscanf(fp, "%[^,],%lf,%lf", 
	   str,
	   &REF_PixcelPos[nnn][0],
	   &REF_PixcelPos[nnn][1])) != EOF){

      if(ret != 3){ break; }
      /*
      printf("P%d:  (%e, %e) \n",
	     nnn, REF_PixcelPos[nnn][0], REF_PixcelPos[nnn][1]);
      */

      nnn++;
  }
  NumOfRef = nnn;
  
  fclose(fp);


  // Read Real Reference Points
  fp = fopen("real-ref-points.csv", "r");
  if(fp == NULL){
    fprintf(stderr, "Can't Open File  \"real-ref-points.csv\"  \n");
    fflush(stderr);
    exit(1);
  }

  nnn = 0;
  while((ret = fscanf(fp, "%[^,],%lf,%lf", 
	   str,
	   &REF_RealPos[nnn][0],
	   &REF_RealPos[nnn][1])) != EOF){

      if(ret != 3){ break;}
      /* 
      printf("P%d:  (%e, %e) \n",
	     nnn, REF_RealPos[nnn][0], REF_RealPos[nnn][1]);
      */

      nnn++;
  }

  if(NumOfRef != nnn){
    fprintf(stderr, "Error:  NumOfRef = %d, nnn %d \n", NumOfRef, nnn);
    fflush(stderr);
    exit(1);    
  }

  fclose(fp);


  for(i=0 ; i<NumOfRef ; i++){
    DLT_A[i*2][0] = REF_RealPos[i][0]; 
    DLT_A[i*2+1][0] = 0.0;

    DLT_A[i*2][1] = 0.0;
    DLT_A[i*2+1][1] = REF_RealPos[i][0]; 

    DLT_A[i*2][2] = REF_RealPos[i][1]; 
    DLT_A[i*2+1][2] = 0.0;

    DLT_A[i*2][3] = 0.0;
    DLT_A[i*2+1][3] = REF_RealPos[i][1]; 

    DLT_A[i*2][4] = 1.0;
    DLT_A[i*2+1][4] = 0.0;

    DLT_A[i*2][5] = 0.0;
    DLT_A[i*2+1][5] = 1.0;

    DLT_A[i*2][6] = - REF_RealPos[i][0]*REF_PixcelPos[i][0]; 
    DLT_A[i*2+1][6] = - REF_RealPos[i][0]*REF_PixcelPos[i][1]; 

    DLT_A[i*2][7] = - REF_RealPos[i][1]*REF_PixcelPos[i][0]; 
    DLT_A[i*2+1][7] = - REF_RealPos[i][1]*REF_PixcelPos[i][1]; 

    DLT_b[i*2] = REF_PixcelPos[i][0];
    DLT_b[i*2+1] = REF_PixcelPos[i][1];
  }


  for(i=0 ; i<NUM_OF_DLT_PARAM ; i++){
    for(j=0 ; j<NUM_OF_DLT_PARAM ; j++){
      DLT_AtA[i][j] = 0.0;
      for(k=0 ; k<NumOfRef ; k++){
	      DLT_AtA[i][j] += DLT_A[k][i]*DLT_A[k][j];
      }
    }
    DLT_x[i] = 0.0;
    for(k=0 ; k<NumOfRef ; k++){
      DLT_x[i] += DLT_A[k][i]*DLT_b[k];
    }
  }

  gauss_solve(DLT_AtA, DLT_x, NUM_OF_DLT_PARAM);

}


void PxPos_to_RealPos(double rpos[2], double pxpos[2])
{
  double mat[2][2];
  double bb[2];
  double detj=0.0;
  int i, j;

  mat[0][0] = DLT_x[3] - DLT_x[7]*pxpos[1];
  mat[1][1] = DLT_x[0] - DLT_x[6]*pxpos[0];
  mat[1][0] = -(DLT_x[1] - DLT_x[6]*pxpos[1]);
  mat[0][1] = -(DLT_x[2] - DLT_x[7]*pxpos[0]);
		
  detj = mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];

  for(i=0 ; i<2 ; i++){
    for(j=0 ; j<2 ; j++){
      mat[i][j] /= detj;
    }
  }

  bb[0] = pxpos[0] - DLT_x[4];
  bb[1] = pxpos[1] - DLT_x[5];

  rpos[0] = mat[0][0]*bb[0] + mat[0][1]*bb[1];
  rpos[1] = mat[1][0]*bb[0] + mat[1][1]*bb[1];
}


void calc_2d_pos()
{
  FILE *fp;
  char str1[256], str2[256], str3[256];
  double real_pos[2];
  double px_pos[2];


  fp = fopen("sample_image_points.csv", "r");
  if(fp == NULL){
    fprintf(stderr, "Can't Open File  \"sample_image_points.csv\"  \n");
    fflush(stderr);
    exit(1);
  }


  if(fscanf(fp, "%[^,],%[^,],%s", str1,str2,str3) != EOF)
  {
    /*
    printf("sample_header:   %s, %s, %s\n", str1,str2,str3);
    fflush(stdout);
    */
  } else{
    fprintf(stderr, "Read Error  \"sample_image_points.csv\"  \n");
    fflush(stderr);
    exit(1);
  }

  int ret;
  NumOfSample = 0;
  while((ret = fscanf(fp, "%[^,],%lf,%lf", 
	   str1,
	   &PixcelPos[NumOfSample][0],
	   &PixcelPos[NumOfSample][1])) != EOF){

      if(ret != 3){ break;  }
      /*
      printf("Sample%d:  (%e, %e) \n",
	    NumOfSample+1, PixcelPos[NumOfSample][0], PixcelPos[NumOfSample][1]);
      fflush(stdout);
      */

      NumOfSample++;
  }
  

  for(int i=0 ; i<NumOfSample ; i++){
    px_pos[0] = PixcelPos[i][0];
    px_pos[1] = PixcelPos[i][1];

    PxPos_to_RealPos(real_pos, px_pos);

    RealPos[i][0] = real_pos[0];
    RealPos[i][1] = real_pos[1];
  }

  fp = fopen("output_2d_points.csv", "w");
  if(fp == NULL){
    fprintf(stderr, "Can't Open File  \"output_2d_points.csv\"  \n");
    fflush(stderr);
    exit(1);
  }


  for(int i=0; i<NumOfSample; i++){
    fprintf(fp, "%d,%e,%e\n",
    i+1, RealPos[i][0], RealPos[i][1]);
  }


  fclose(fp);
}



int main()
{
  calc_DLT_param();
  calc_2d_pos();
  return 0;
}
