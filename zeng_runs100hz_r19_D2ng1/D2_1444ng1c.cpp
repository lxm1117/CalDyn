//CaM is reduecd to 1/voxel in layer1.
//include the Ca pump, Ca influx locate in the middle 16 voxel.
//increase the rate constant from bound state to trapped state to 10 /S
//include no pp1

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "zeng_canng.h"

using namespace std;

double const D_CaMKII_C1 = 0.002;
double const D_CaMKII_C2 = 0;
double const D_Ca_C = 0.223; 
double const D_CaM00_C = 0.0223;
double const D_CaM01_C = 0.0223;
double const D_CaM10_C = 0.0223;
double const D_CaM11_C = 0.0223;
double const D_CaM02_C = 0.0223;
double const D_CaM20_C = 0.0223;
double const D_CaM12_C = 0.0223;
double const D_CaM21_C = 0.0223;
double const D_CaM22_C = 0.0223;
double const D_CB00 = 0.02;
double const D_CB10 = 0.02; 
double const D_CB01 = 0.02; 
double const D_CB11 = 0.02; 
double const D_CB20 = 0.02; 
double const D_CB02 = 0.02; 
double const D_CB21 = 0.02;
double const D_CB12 = 0.02;
double const D_CB22 = 0.02;

struct voxel layer1[10][10];
struct voxel layer2[5][5][5];
struct voxel layer3[8];
struct voxel layer4[4][2][2];

void initial();
struct CaMKII *creat(int n);
double transition_lamda();
double layer1_diffusion_lamda(int i, int j);
double layer2_diffusion_lamda(int i, int j, int k);
double layer3_diffusion_lamda(int i);
double layer4_diffusion_lamda(int i, int j, int k);
double voxel_chemical_lamda(struct voxel grid);
void reaction_layer1(double choose2, int p, int q);
void reaction_layer2(double choose2, int p, int q, int r);
void reaction_layer3(double choose2, int p);
void reaction_layer4(double choose2, int p, int q, int r);
void chemical_reaction(struct voxel *p1, double choose2, int index);
int CaMKII_number(struct voxel grid);
double ran2(long *idum);
void reaction();
struct CaMKII * delete_CaMKII(struct voxel grid, struct CaMKII *point);
struct CaMKII * insert_CaMKII(struct voxel grid, struct CaMKII *point);

double reaction_lamda;
double random_num;
long *idum;
int x;
int y;
int z;
int total_Ca_influx;
int total_Ca_pumped;
int total_Ca_pumped1;
int total_Ca_pumped2;
int total_Ca_pumped3;
int total_Ca_pumped4;

int main(){
  int i, j, k, m, n, num, index, flag;
  struct CaMKII *p;
  long l;
  double time, influx[100000], input, old_Ca_influx, Ca_lamda, Ca_lamda1, total_Ca, total_Ca1;
  FILE *fp1;
  fp1=fopen("D2_1444ng1c.dat","wb"); 

  int Ca, Ca1, Ca2, Ca4, CaMKII, CaMKIICaM00, CaMKIICaM10, CaMKIICaM01, CaMKIICaM11, CaMKIICaM20, CaMKIICaM02, CaMKIICaM12, CaMKIICaM21, CaMKIICaM22;
  int Trapped00, Trapped01, Trapped10, Trapped11, Trapped02, Trapped20, Trapped12, Trapped21, Trapped22, Auton, Capped;
  int CaM00, CaM00_1, CaM01, CaM10, CaM11, CaM20, CaM02, CaM12, CaM21, CaM22;
  int CaN, CaNCaM00, CaNCaM10, CaNCaM01, CaNCaM11, CaNCaM02, CaNCaM20, CaNCaM12, CaNCaM21, CaNCaM22;
  int Ng, NgCaM00, NgCaM10, NgCaM01, NgCaM11, NgCaM02, NgCaM20, NgCaM12, NgCaM21, NgCaM22;
  int Trapped0, Trapped1, Trapped2, Trapped3;
  int CB0, CB1, CB2, CB3, CB4;
  int bound, trapped;
  
  FILE *fp;

  if((fp = fopen("../gr19ca.100.1444", "r")) == NULL){
    printf("Cannot open file\n");
    exit(0);
  }
 
  for(i=0; i<100000; i++){
    for(j=0; j<10; j++){
      fscanf(fp, "%lf", &input);
      if(j == 5) influx[i] = input;
    }
  }
  
  l = 1000;
  idum = &l;
  time = 0;
  num = 0;
  Ca = 0;
  Ca1 = 0;
  Ca2 = 0;
  Ca4 = 0;
  CaMKII = 0;
  CaMKIICaM00 = 0;
  CaMKIICaM01 = 0;
  CaMKIICaM10 = 0;
  CaMKIICaM11 = 0;
  CaMKIICaM20 = 0;
  CaMKIICaM02 = 0;
  CaMKIICaM12 = 0;
  CaMKIICaM21 = 0;
  CaMKIICaM22 = 0;
  Trapped00 = 0;
  Trapped01 = 0;
  Trapped10 = 0;
  Trapped11 = 0;
  Trapped02 = 0;
  Trapped20 = 0;
  Trapped12 = 0;
  Trapped21 = 0;
  Trapped22 = 0;
  Auton = 0;
  Capped = 0;
  CaM00 = 0;
  CaM00_1 = 0;
  CaM01 = 0;
  CaM10 = 0;
  CaM11 = 0;
  CaM02 = 0;
  CaM20 = 0;
  CaM12 = 0;
  CaM21 = 0;
  CaM22 = 0;
  Ng = 0;
  NgCaM00 = 0;
  NgCaM01 = 0;
  NgCaM10 = 0;
  NgCaM11 = 0;
  NgCaM02 = 0;
  NgCaM20 = 0;
  NgCaM12 = 0;
  NgCaM21 = 0;
  NgCaM22 = 0;
  CaN = 0;
  CaNCaM00 = 0;
  CaNCaM01 = 0;
  CaNCaM10 = 0;
  CaNCaM11 = 0;
  CaNCaM02 = 0;
  CaNCaM20 = 0;
  CaNCaM12 = 0;
  CaNCaM21 = 0;
  CaNCaM22 = 0;
  Trapped0 = 0;
  Trapped1 = 0;
  Trapped2 = 0;
  Trapped3 = 0;
  CB0 = 0;
  CB1 = 0;
  CB2 = 0;
  CB3 = 0;
  CB4 = 0;
  bound = 0;
  trapped = 0;
  total_Ca_influx=0;
  total_Ca_pumped=0;
  total_Ca_pumped1=0;
  total_Ca_pumped2=0;
  total_Ca_pumped3=0;
  total_Ca_pumped4=0;
  flag = 1;
  
  initial();

  reaction_lamda = transition_lamda();
  old_Ca_influx = 0;
  Ca_lamda1 = 0;
  total_Ca = 0;
  while(time < 1000){
    total_Ca1 = total_Ca;
    if(time < 1000) {
       index = (int) (time / 0.01);
       total_Ca = (influx[index] * 6.02 * pow(10, 4) / (2 * 9.648456));
      }else {
	total_Ca = 0; 
    }
    if(total_Ca != total_Ca1){
        for(i=3; i<7; i++)
	for(j=3; j<7; j++)
	  layer1[i][j].Ca_influx = total_Ca / 16;
   }

    if(influx[index] != old_Ca_influx){
      old_Ca_influx = influx[index];
      Ca_lamda = 0;
      for(i=0; i<10; i++)
	for(j=0; j<10; j++){
	  layer1[i][j].diffusion_lamda += layer1[i][j].Ca_influx - total_Ca1 / 100;
	  Ca_lamda += layer1[i][j].Ca_influx;
	}
      reaction_lamda += Ca_lamda - Ca_lamda1;
      Ca_lamda1 = Ca_lamda;
    }

    time += log(1.0/ran2(idum)) / reaction_lamda;
     
    if(time > (num * 1)){
      for(i=0; i<10; i++)
	for(j=0; j<10; j++){
	  Ca += layer1[i][j].Ca;
	  Ca1 += layer1[i][j].Ca;
	  Ca2 += layer1[i][j].Ca;
	  CaM00 += layer1[i][j].CaM00;
	  CaM00_1 += layer1[i][j].CaM00;
	  CaM01 += layer1[i][j].CaM01;
          CaM10 += layer1[i][j].CaM10;
	  CaM11 += layer1[i][j].CaM11;
	  CaM02 += layer1[i][j].CaM02;
	  CaM20 += layer1[i][j].CaM20;
	  CaM12 += layer1[i][j].CaM12;
	  CaM21 += layer1[i][j].CaM21;
	  CaM22 += layer1[i][j].CaM22;
	  Ng += layer1[i][j].Ng;
	  NgCaM00 += layer1[i][j].NgCaM00;
	  NgCaM01 += layer1[i][j].NgCaM01;
	  NgCaM10 += layer1[i][j].NgCaM10;
	  NgCaM11 += layer1[i][j].NgCaM11;
	  NgCaM02 += layer1[i][j].NgCaM02;
	  NgCaM20 += layer1[i][j].NgCaM20;
	  NgCaM12 += layer1[i][j].NgCaM12;
	  NgCaM21 += layer1[i][j].NgCaM21;
	  NgCaM22 += layer1[i][j].NgCaM22;
	  CaN += layer1[i][j].CaN;
	  CaNCaM00 += layer1[i][j].CaNCaM00;
	  CaNCaM01 += layer1[i][j].CaNCaM01;
	  CaNCaM10 += layer1[i][j].CaNCaM10;
	  CaNCaM11 += layer1[i][j].CaNCaM11;
	  CaNCaM02 += layer1[i][j].CaNCaM02;
	  CaNCaM20 += layer1[i][j].CaNCaM20;
	  CaNCaM12 += layer1[i][j].CaNCaM12;
	  CaNCaM21 += layer1[i][j].CaNCaM21;
	  CaNCaM22 += layer1[i][j].CaNCaM22;
	  CB0 += layer1[i][j].CB00;
          CB1 += layer1[i][j].CB10;
          CB1 += layer1[i][j].CB01;
          CB2 += layer1[i][j].CB11;
          CB2 += layer1[i][j].CB20;
          CB2 += layer1[i][j].CB02;
          CB3 += layer1[i][j].CB21;
          CB3 += layer1[i][j].CB12;
          CB4 += layer1[i][j].CB22;
	  
	  p = layer1[i][j].head;
	  while(p != NULL){
	    for(m=0; m<2; m++)
	      for(n=0; n<6; n++){
		if(p->subunit[m][n] == 0){
		  CaMKII++;
		}else if(p->subunit[m][n] == 1){
		  CaMKIICaM00++;
		}else if(p->subunit[m][n] == 2){
		  CaMKIICaM10++;
		}else if(p->subunit[m][n] == 3){
		  CaMKIICaM01++;
		}else if(p->subunit[m][n] == 4){
		  CaMKIICaM11++;
		}else if(p->subunit[m][n] == 5){
		  CaMKIICaM20++;
		}else if(p->subunit[m][n] == 6){
		  CaMKIICaM02++;
		}else if(p->subunit[m][n] == 7){
		  CaMKIICaM21++;
		}else if(p->subunit[m][n] == 8){
		  CaMKIICaM12++;
		}else if(p->subunit[m][n] == 9){
		  CaMKIICaM22++;
		}else if(p->subunit[m][n] == 10){
		  Trapped00++;
		}else if(p->subunit[m][n] == 11){
		  Trapped10++;
		}else if(p->subunit[m][n] == 12){
		  Trapped01++;
		}else if(p->subunit[m][n] == 13){
		  Trapped11++;
		}else if(p->subunit[m][n] == 14){
		  Trapped20++;
		}else if(p->subunit[m][n] == 15){
		  Trapped02++;
		}else if(p->subunit[m][n] == 16){
		  Trapped21++;
		}else if(p->subunit[m][n] == 17){
		  Trapped12++;
		}else if(p->subunit[m][n] == 18){
		  Trapped22++;
		}else if(p->subunit[m][n] == 19){
		  Auton++;
		}else if(p->subunit[m][n] == 20){
		  Capped++;
		}
	      }
	    p = p->next;
	  }
	}
        
       for(i=0; i<5; i++)
    	for(j=0; j<5; j++)
	  for(k=0; k<5; k++){
	    Ca += layer2[i][j][k].Ca;
	    Ca1 += layer2[i][j][k].Ca;
	    CaM00_1 += layer2[i][j][k].CaM00;
	    p = layer2[i][j][k].head;
	    while(p != NULL){
	      for(m=0; m<2; m++)
		for(n=0; n<6; n++){
		  if(p->subunit[m][n] == 0){
		  CaMKII++;
		}else if(p->subunit[m][n] == 1){
		  CaMKIICaM00++;
		}else if(p->subunit[m][n] == 2){
		  CaMKIICaM10++;
		}else if(p->subunit[m][n] == 3){
		  CaMKIICaM01++;
		}else if(p->subunit[m][n] == 4){
		  CaMKIICaM11++;
		}else if(p->subunit[m][n] == 5){
		  CaMKIICaM20++;
		}else if(p->subunit[m][n] == 6){
		  CaMKIICaM02++;
		}else if(p->subunit[m][n] == 7){
		  CaMKIICaM21++;
		}else if(p->subunit[m][n] == 8){
		  CaMKIICaM12++;
		}else if(p->subunit[m][n] == 9){
		  CaMKIICaM22++;
		}else if(p->subunit[m][n] == 10){
		  Trapped00++;
		}else if(p->subunit[m][n] == 11){
		  Trapped10++;
		}else if(p->subunit[m][n] == 12){
		  Trapped01++;
		}else if(p->subunit[m][n] == 13){
		  Trapped11++;
		}else if(p->subunit[m][n] == 14){
		  Trapped20++;
		}else if(p->subunit[m][n] == 15){
		  Trapped02++;
		}else if(p->subunit[m][n] == 16){
		  Trapped21++;
		}else if(p->subunit[m][n] == 17){
		  Trapped12++;
		}else if(p->subunit[m][n] == 18){
		  Trapped22++;
		}else if(p->subunit[m][n] == 19){
		  Auton++;
		}else if(p->subunit[m][n] == 20){
		  Capped++;
		}
		}
	      p = p->next;
	    }
	  }
        
       for(i=0; i<8; i++){
	 Ca1 += layer3[i].Ca;
	 p = layer3[i].head;
	 while(p != NULL){
	   for(m=0; m<2; m++)
	     for(n=0; n<6; n++){
	      if(p->subunit[m][n] == 0){
		  CaMKII++;
		}else if(p->subunit[m][n] == 1){
		  CaMKIICaM00++;
		}else if(p->subunit[m][n] == 2){
		  CaMKIICaM10++;
		}else if(p->subunit[m][n] == 3){
		  CaMKIICaM01++;
		}else if(p->subunit[m][n] == 4){
		  CaMKIICaM11++;
		}else if(p->subunit[m][n] == 5){
		  CaMKIICaM20++;
		}else if(p->subunit[m][n] == 6){
		  CaMKIICaM02++;
		}else if(p->subunit[m][n] == 7){
		  CaMKIICaM21++;
		}else if(p->subunit[m][n] == 8){
		  CaMKIICaM12++;
		}else if(p->subunit[m][n] == 9){
		  CaMKIICaM22++;
		}else if(p->subunit[m][n] == 10){
		  Trapped00++;
		}else if(p->subunit[m][n] == 11){
		  Trapped10++;
		}else if(p->subunit[m][n] == 12){
		  Trapped01++;
		}else if(p->subunit[m][n] == 13){
		  Trapped11++;
		}else if(p->subunit[m][n] == 14){
		  Trapped20++;
		}else if(p->subunit[m][n] == 15){
		  Trapped02++;
		}else if(p->subunit[m][n] == 16){
		  Trapped21++;
		}else if(p->subunit[m][n] == 17){
		  Trapped12++;
		}else if(p->subunit[m][n] == 18){
		  Trapped22++;
		}else if(p->subunit[m][n] == 19){
		  Auton++;
		}else if(p->subunit[m][n] == 20){
		  Capped++;
		}
	     }
	   p = p->next;
	  }
       }
       
       for(i=0; i<4; i++)
	for(j=0; j<2; j++)
	  for(k=0; k<2; k++){
	    Ca1 += layer4[i][j][k].Ca;
	    Ca4 += layer4[i][j][k].Ca;
	    p = layer4[i][j][k].head;
	    while(p != NULL){
	      for(m=0; m<2; m++)
		for(n=0; n<6; n++){
		  if(p->subunit[m][n] == 0){
		  CaMKII++;
		}else if(p->subunit[m][n] == 1){
		  CaMKIICaM00++;
		}else if(p->subunit[m][n] == 2){
		  CaMKIICaM10++;
		}else if(p->subunit[m][n] == 3){
		  CaMKIICaM01++;
		}else if(p->subunit[m][n] == 4){
		  CaMKIICaM11++;
		}else if(p->subunit[m][n] == 5){
		  CaMKIICaM20++;
		}else if(p->subunit[m][n] == 6){
		  CaMKIICaM02++;
		}else if(p->subunit[m][n] == 7){
		  CaMKIICaM21++;
		}else if(p->subunit[m][n] == 8){
		  CaMKIICaM12++;
		}else if(p->subunit[m][n] == 9){
		  CaMKIICaM22++;
		}else if(p->subunit[m][n] == 10){
		  Trapped00++;
		}else if(p->subunit[m][n] == 11){
		  Trapped10++;
		}else if(p->subunit[m][n] == 12){
		  Trapped01++;
		}else if(p->subunit[m][n] == 13){
		  Trapped11++;
		}else if(p->subunit[m][n] == 14){
		  Trapped20++;
		}else if(p->subunit[m][n] == 15){
		  Trapped02++;
		}else if(p->subunit[m][n] == 16){
		  Trapped21++;
		}else if(p->subunit[m][n] == 17){
		  Trapped12++;
		}else if(p->subunit[m][n] == 18){
		  Trapped22++;
		}else if(p->subunit[m][n] == 19){
		  Auton++;
		}else if(p->subunit[m][n] == 20){
		  Capped++;
		}
		}
	      p = p->next;
	    }
	  }
      
       fprintf(fp1, "%f %d %f %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", time, Ca2, Ca2/7.5, Ca, Ca1, Ca4, CaMKII, CaMKIICaM00, CaMKIICaM10, CaMKIICaM01, CaMKIICaM11, CaMKIICaM20, CaMKIICaM02, CaMKIICaM21, CaMKIICaM12, CaMKIICaM22, Trapped00, Trapped10, Trapped01, Trapped11, Trapped20, Trapped02, Trapped21, Trapped12, Trapped22, Auton, Capped, CaM00, CaM10, CaM01, CaM11, CaM20, CaM02, CaM21, CaM12, CaM22, CaM00_1, CB0, CB1, CB2, CB3, CB4, CaN, CaNCaM00, CaNCaM10, CaNCaM01, CaNCaM11, CaNCaM20, CaNCaM02, CaNCaM21, CaNCaM12, CaNCaM22, Ng, NgCaM00, NgCaM10, NgCaM01, NgCaM11, NgCaM20, NgCaM02, NgCaM21, NgCaM12, NgCaM22, total_Ca_influx, total_Ca_pumped, total_Ca_pumped1, total_Ca_pumped2, total_Ca_pumped3, total_Ca_pumped4);

       num++;

       Ca = 0;
       Ca1 = 0;
       Ca2 = 0;
       Ca4 = 0;
       CaMKII = 0;
       CaMKIICaM00 = 0;
       CaMKIICaM01 = 0;
       CaMKIICaM10 = 0;
       CaMKIICaM11 = 0;
       CaMKIICaM20 = 0;
       CaMKIICaM02 = 0;
       CaMKIICaM12 = 0;
       CaMKIICaM21 = 0;
       CaMKIICaM22 = 0;
       Trapped00 = 0;
       Trapped01 = 0;
       Trapped10 = 0;
       Trapped11 = 0;
       Trapped02 = 0;
       Trapped20 = 0;
       Trapped12 = 0;
       Trapped21 = 0;
       Trapped22 = 0;
       Auton = 0;
       Capped = 0;
       CaM00 = 0;
       CaM00_1 = 0;
       CaM01 = 0;
       CaM10 = 0;
       CaM11 = 0;
       CaM02 = 0;
       CaM20 = 0;
       CaM12 = 0;
       CaM21 = 0;
       CaM22 = 0;
       CaN = 0;
       CaNCaM00 = 0;
       CaNCaM01 = 0;
       CaNCaM10 = 0;
       CaNCaM11 = 0;
       CaNCaM02 = 0;
       CaNCaM20 = 0;
       CaNCaM12 = 0;
       CaNCaM21 = 0;
       CaNCaM22 = 0;
       Ng = 0;
       NgCaM00 = 0;
       NgCaM01 = 0;
       NgCaM10 = 0;
       NgCaM11 = 0;
       NgCaM02 = 0;
       NgCaM20 = 0;
       NgCaM12 = 0;
       NgCaM21 = 0;
       NgCaM22 = 0;
	   Trapped0 = 0;
	   Trapped1 = 0;
	   Trapped2 = 0;
	   Trapped3 = 0;
	CB0 = 0;
	CB1 = 0;
	CB2 = 0;
	CB3 = 0;
	CB4 = 0;
	bound = 0;
	trapped = 0;

       cout << "time = " << time << " " << reaction_lamda << endl;
    }
    reaction();
  }
    
  return 0;
}

void initial(){
  int i, j, k, m, inter;
  double random_number;
  struct CaMKII *head;
  /******************************************************************************/
  //initialize the Ca2+ influx 
  for(i=0; i<10; i++)
    for(j=0; j<10; j++)
      layer1[i][j].Ca_influx = 0;
  
  for(i=0; i<5; i++)
    for(j=0; j<5; j++)
      for(k=0; k<5; k++)
	layer2[i][j][k].Ca_influx = 0;
  
  for(i=0; i<8; i++)
    layer3[i].Ca_influx = 0;
  
  for(i=0; i<4; i++)
    for(j=0; j<2; j++)
      for(k=0; k<2; k++)
	layer4[i][j][k].Ca_influx = 0;

 for(i=0; i<10; i++)
    for(j=0; j<10; j++)
      layer1[i][j].pump_area = 0.05 * 0.05;

  for(i=0; i<10; i++)
    for(j=0; j<10; j++){
      if(i == 0) layer1[i][j].pump_area += 0.05 * 0.05;
      if(i == 9) layer1[i][j].pump_area += 0.05 * 0.05;
      if(j == 0) layer1[i][j].pump_area += 0.05 * 0.05;
      if(j == 9) layer1[i][j].pump_area += 0.05 * 0.05;
    }

  //****************************************************************************/
  //initialize pump area in layer2

 for(i=0; i<5; i++)
    for(j=0; j<5; j++)
      for(k=0; k<5; k++)
	layer2[i][j][k].pump_area = 0;

  for(i=0; i<5; i++)
    for(j=0; j<5; j++)
      for(k=0; k<5; k++){
	if(i == 0) layer2[i][j][k].pump_area += 0.1 * 0.1;
	if(i == 4) layer2[i][j][k].pump_area += 0.1 * 0.1;
	if(j == 0) layer2[i][j][k].pump_area += 0.1 * 0.1;
	if(j == 4) layer2[i][j][k].pump_area += 0.1 * 0.1;
	if(k == 0) layer2[i][j][k].pump_area += 0.1 * 0.1;
	if(k == 4) layer2[i][j][k].pump_area += 0.1 * 0.1;
      }

  //****************************************************************************/
  //initialize pump area in layer3
  for(i=0; i<8; i++)
    layer3[i].pump_area = 4.0 * 0.1 * 0.1;

  //****************************************************************************/
  //initialize pump area in layer4
 
 for(i=0; i<4; i++)
    for(j=0; j<2; j++)
      for(k=0; k<2; k++)
	layer4[i][j][k].pump_area = 2.0 * 0.5 * 0.5;


 for(i=0; i<4; i++)
    for(j=0; j<2; j++)
      for(k=0; k<2; k++){
	if(i == 0) layer4[i][j][k].pump_area += 0.5 * 0.5;
	if(i == 3) layer4[i][j][k].pump_area += 0.5 * 0.5;
      }

  //initialize  diffusion in layer1
  for(i=0; i<10; i++)
    for(j=0; j<10; j++){
      for(k=0; k<4; k++){
	layer1[i][j].D_CaMKII[k] = D_CaMKII_C1 * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
	layer1[i][j].D_Ca[k] = D_Ca_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
	layer1[i][j].D_CaM00[k] = D_CaM00_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
        layer1[i][j].D_CaM10[k] = D_CaM10_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
	layer1[i][j].D_CaM01[k] = D_CaM01_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
	layer1[i][j].D_CaM11[k] = D_CaM11_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
	layer1[i][j].D_CaM02[k] = D_CaM02_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
	layer1[i][j].D_CaM20[k] = D_CaM20_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
	layer1[i][j].D_CaM21[k] = D_CaM21_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
	layer1[i][j].D_CaM12[k] = D_CaM12_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
	layer1[i][j].D_CaM22[k] = D_CaM22_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
	layer1[i][j].D_CB00[k] = D_CB00 * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
	layer1[i][j].D_CB01[k] = D_CB01 * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
	layer1[i][j].D_CB10[k] = D_CB10 * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
	layer1[i][j].D_CB11[k] = D_CB11 * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
	layer1[i][j].D_CB02[k] = D_CB02 * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
	layer1[i][j].D_CB20[k] = D_CB20 * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
	layer1[i][j].D_CB21[k] = D_CB21 * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
	layer1[i][j].D_CB12[k] = D_CB12 * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
	layer1[i][j].D_CB22[k] = D_CB22 * (0.05 * 0.05) / (pow(0.05, 3) * 0.05);
      }
    }

  //connection between laryer1 and layer2	  
  for(i=0; i<10; i++)
    for(j=0; j<10; j++){
      layer1[i][j].D_CaMKII[5] = 0;
	//layer1[i][j].D_CaMKII[5] = D_CaMKII_C1 * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
      layer1[i][j].D_Ca[5] = D_Ca_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
      layer1[i][j].D_CaM00[5] = D_CaM00_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
      layer1[i][j].D_CaM10[5] = D_CaM10_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
      layer1[i][j].D_CaM01[5] = D_CaM01_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
      layer1[i][j].D_CaM11[5] = D_CaM11_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
	  layer1[i][j].D_CaM02[5] = D_CaM02_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
	  layer1[i][j].D_CaM20[5] = D_CaM20_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
	  layer1[i][j].D_CaM21[5] = D_CaM21_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
	  layer1[i][j].D_CaM12[5] = D_CaM12_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
	  layer1[i][j].D_CaM22[5] = D_CaM22_C * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
	  layer1[i][j].D_CB00[5] = D_CB00 * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
	  layer1[i][j].D_CB01[5] = D_CB01 * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
	  layer1[i][j].D_CB10[5] = D_CB10 * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
	  layer1[i][j].D_CB11[5] = D_CB11 * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
	  layer1[i][j].D_CB20[5] = D_CB20 * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
	  layer1[i][j].D_CB02[5] = D_CB02 * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
	  layer1[i][j].D_CB21[5] = D_CB21 * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
	  layer1[i][j].D_CB12[5] = D_CB12 * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
	  layer1[i][j].D_CB22[5] = D_CB22 * (0.05 * 0.05) / (pow(0.05, 3) * 0.075);
    }
  
//set the boudary
  for(j=0; j<10; j++){
    layer1[0][j].D_CaMKII[0] = 0;
    layer1[0][j].D_Ca[0] = 0;
    layer1[0][j].D_CaM00[0] = 0;
	layer1[0][j].D_CaM10[0] = 0;
	layer1[0][j].D_CaM01[0] = 0;
	layer1[0][j].D_CaM11[0] = 0;
	layer1[0][j].D_CaM02[0] = 0;
	layer1[0][j].D_CaM20[0] = 0;
	layer1[0][j].D_CaM21[0] = 0;
	layer1[0][j].D_CaM12[0] = 0;
	layer1[0][j].D_CaM22[0] = 0;
	layer1[0][j].D_CB00[0] = 0;
	layer1[0][j].D_CB01[0] = 0;
	layer1[0][j].D_CB10[0] = 0;
	layer1[0][j].D_CB11[0] = 0;	
	layer1[0][j].D_CB02[0] = 0;	
	layer1[0][j].D_CB20[0] = 0;	
	layer1[0][j].D_CB21[0] = 0;	
	layer1[0][j].D_CB12[0] = 0;
	layer1[0][j].D_CB22[0] = 0;
  }
  
  for(i=0; i<10; i++){
    layer1[i][0].D_CaMKII[1] = 0;
    layer1[i][0].D_Ca[1] = 0;
    layer1[i][0].D_CaM00[1] = 0;
	layer1[i][0].D_CaM01[1] = 0;
	layer1[i][0].D_CaM10[1] = 0;
	layer1[i][0].D_CaM11[1] = 0;
	layer1[i][0].D_CaM20[1] = 0;
	layer1[i][0].D_CaM02[1] = 0;
	layer1[i][0].D_CaM21[1] = 0;
	layer1[i][0].D_CaM12[1] = 0;
	layer1[i][0].D_CaM22[1] = 0;    
	layer1[i][0].D_CB00[1] = 0;
	layer1[i][0].D_CB10[1] = 0;
	layer1[i][0].D_CB01[1] = 0;
	layer1[i][0].D_CB11[1] = 0;
	layer1[i][0].D_CB20[1] = 0;
	layer1[i][0].D_CB02[1] = 0;
	layer1[i][0].D_CB21[1] = 0;
	layer1[i][0].D_CB12[1] = 0;
	layer1[i][0].D_CB22[1] = 0;
  }

  for(j=0; j<10; j++){
    layer1[9][j].D_CaMKII[2] = 0;
    layer1[9][j].D_Ca[2] = 0;
    layer1[9][j].D_CaM00[2] = 0;
	layer1[9][j].D_CaM10[2] = 0;
	layer1[9][j].D_CaM01[2] = 0;
	layer1[9][j].D_CaM11[2] = 0;
	layer1[9][j].D_CaM20[2] = 0;
	layer1[9][j].D_CaM02[2] = 0;
	layer1[9][j].D_CaM12[2] = 0;
	layer1[9][j].D_CaM21[2] = 0;
	layer1[9][j].D_CaM22[2] = 0;    
	layer1[9][j].D_CB00[2] = 0;
	layer1[9][j].D_CB10[2] = 0;
	layer1[9][j].D_CB01[2] = 0;
	layer1[9][j].D_CB11[2] = 0;
	layer1[9][j].D_CB20[2] = 0;
	layer1[9][j].D_CB02[2] = 0;
	layer1[9][j].D_CB21[2] = 0;
	layer1[9][j].D_CB12[2] = 0;
	layer1[9][j].D_CB22[2] = 0;
  }
  
  for(i=0; i<10; i++){
    layer1[i][9].D_CaMKII[3] = 0;
    layer1[i][9].D_Ca[3] = 0;
    layer1[i][9].D_CaM00[3] = 0;
	layer1[i][9].D_CaM01[3] = 0;
	layer1[i][9].D_CaM10[3] = 0;
	layer1[i][9].D_CaM11[3] = 0;
	layer1[i][9].D_CaM20[3] = 0;
	layer1[i][9].D_CaM02[3] = 0;
	layer1[i][9].D_CaM21[3] = 0;
	layer1[i][9].D_CaM12[3] = 0;
	layer1[i][9].D_CaM22[3] = 0;    
	layer1[i][9].D_CB00[3] = 0;
	layer1[i][9].D_CB01[3] = 0;
	layer1[i][9].D_CB10[3] = 0;
	layer1[i][9].D_CB11[3] = 0;
	layer1[i][9].D_CB20[3] = 0;
	layer1[i][9].D_CB02[3] = 0;
	layer1[i][9].D_CB21[3] = 0;
	layer1[i][9].D_CB12[3] = 0;
	layer1[i][9].D_CB22[3] = 0;
   }
  /**********************************************************************************/
  //initialize diffusion in layer2

  for(i=0; i<5; i++)
    for(j=0; j<5; j++)
      for(k=0; k<5; k++)
	for(m=0; m<6; m++){
	  layer2[i][j][k].D_CaMKII[m] = D_CaMKII_C2 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer2[i][j][k].D_Ca[m] = D_Ca_C * (0.1 * 0.1) / (pow(0.1,3) * 0.1);
	  layer2[i][j][k].D_CaM00[m] = D_CaM00_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer2[i][j][k].D_CaM10[m] = D_CaM10_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer2[i][j][k].D_CaM01[m] = D_CaM01_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer2[i][j][k].D_CaM11[m] = D_CaM11_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer2[i][j][k].D_CaM02[m] = D_CaM02_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer2[i][j][k].D_CaM20[m] = D_CaM20_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer2[i][j][k].D_CaM21[m] = D_CaM21_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer2[i][j][k].D_CaM12[m] = D_CaM12_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer2[i][j][k].D_CaM22[m] = D_CaM22_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);	  
	  layer2[i][j][k].D_CB00[m] = D_CB00 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer2[i][j][k].D_CB01[m] = D_CB01 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer2[i][j][k].D_CB10[m] = D_CB10 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer2[i][j][k].D_CB11[m] = D_CB11 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer2[i][j][k].D_CB20[m] = D_CB20 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer2[i][j][k].D_CB02[m] = D_CB02 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer2[i][j][k].D_CB21[m] = D_CB21 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer2[i][j][k].D_CB12[m] = D_CB12 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer2[i][j][k].D_CB22[m] = D_CB22 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	}				

  //connection between layer2 and layer1
  for(i=0; i<5; i++)
    for(j=0; j<5; j++){
      layer2[i][j][0].D_CaMKII[4] = D_CaMKII_C2 * (0.05 * 0.05) / (pow(0.1,3) * 0.075);
      layer2[i][j][0].D_Ca[4] = D_Ca_C * (0.05 * 0.05) / (pow(0.1,3) * 0.075);
      layer2[i][j][0].D_CaM00[4] = D_CaM00_C * (0.05 * 0.05) / (pow(0.1, 3) * 0.075);
	  layer2[i][j][0].D_CaM10[4] = D_CaM10_C * (0.05 * 0.05) / (pow(0.1, 3) * 0.075);
	  layer2[i][j][0].D_CaM01[4] = D_CaM01_C * (0.05 * 0.05) / (pow(0.1, 3) * 0.075);
	  layer2[i][j][0].D_CaM11[4] = D_CaM11_C * (0.05 * 0.05) / (pow(0.1, 3) * 0.075);
	  layer2[i][j][0].D_CaM02[4] = D_CaM02_C * (0.05 * 0.05) / (pow(0.1, 3) * 0.075);
	  layer2[i][j][0].D_CaM20[4] = D_CaM20_C * (0.05 * 0.05) / (pow(0.1, 3) * 0.075);
	  layer2[i][j][0].D_CaM12[4] = D_CaM12_C * (0.05 * 0.05) / (pow(0.1, 3) * 0.075);
	  layer2[i][j][0].D_CaM21[4] = D_CaM21_C * (0.05 * 0.05) / (pow(0.1, 3) * 0.075);
	  layer2[i][j][0].D_CaM22[4] = D_CaM22_C * (0.05 * 0.05) / (pow(0.1, 3) * 0.075);      
	  layer2[i][j][0].D_CB00[4] = D_CB00 * (0.05 * 0.05) / (pow(0.1, 3) * 0.075);
	  layer2[i][j][0].D_CB01[4] = D_CB01 * (0.05 * 0.05) / (pow(0.1, 3) * 0.075);
      layer2[i][j][0].D_CB10[4] = D_CB10 * (0.05 * 0.05) / (pow(0.1, 3) * 0.075);
      layer2[i][j][0].D_CB11[4] = D_CB11 * (0.05 * 0.05) / (pow(0.1, 3) * 0.075);
      layer2[i][j][0].D_CB20[4] = D_CB20 * (0.05 * 0.05) / (pow(0.1, 3) * 0.075);
      layer2[i][j][0].D_CB02[4] = D_CB02 * (0.05 * 0.05) / (pow(0.1, 3) * 0.075);
      layer2[i][j][0].D_CB21[4] = D_CB21 * (0.05 * 0.05) / (pow(0.1, 3) * 0.075);
      layer2[i][j][0].D_CB12[4] = D_CB12 * (0.05 * 0.05) / (pow(0.1, 3) * 0.075);
      layer2[i][j][0].D_CB22[4] = D_CB22 * (0.05 * 0.05) / (pow(0.1, 3) * 0.075);
     }

  //set boundary
  
  for(j=0; j<5; j++)
    for(k=0; k<5; k++){
      layer2[0][j][k].D_CaMKII[0] = 0;
      layer2[0][j][k].D_Ca[0] = 0;
      layer2[0][j][k].D_CaM00[0] = 0;
	  layer2[0][j][k].D_CaM10[0] = 0;
	  layer2[0][j][k].D_CaM01[0] = 0;
	  layer2[0][j][k].D_CaM11[0] = 0;
	  layer2[0][j][k].D_CaM20[0] = 0;
	  layer2[0][j][k].D_CaM02[0] = 0;
	  layer2[0][j][k].D_CaM21[0] = 0;
	  layer2[0][j][k].D_CaM12[0] = 0;
	  layer2[0][j][k].D_CaM22[0] = 0;      
	  layer2[0][j][k].D_CB00[0] = 0;
	  layer2[0][j][k].D_CB01[0] = 0;
      layer2[0][j][k].D_CB10[0] = 0;
      layer2[0][j][k].D_CB11[0] = 0;
      layer2[0][j][k].D_CB20[0] = 0;
      layer2[0][j][k].D_CB02[0] = 0;
      layer2[0][j][k].D_CB21[0] = 0;
      layer2[0][j][k].D_CB12[0] = 0;
      layer2[0][j][k].D_CB22[0] = 0;
    }
  
  for(i=0; i<5; i++)
    for(k=0; k<5; k++){
      layer2[i][0][k].D_CaMKII[1] = 0;
      layer2[i][0][k].D_Ca[1] = 0;
      layer2[i][0][k].D_CaM00[1] = 0;
      layer2[i][0][k].D_CaM01[1] = 0;
	  layer2[i][0][k].D_CaM10[1] = 0;
	  layer2[i][0][k].D_CaM11[1] = 0;
	  layer2[i][0][k].D_CaM20[1] = 0;
	  layer2[i][0][k].D_CaM02[1] = 0;
	  layer2[i][0][k].D_CaM21[1] = 0;
	  layer2[i][0][k].D_CaM12[1] = 0;
	  layer2[i][0][k].D_CaM22[1] = 0;
	  layer2[i][0][k].D_CB00[1] = 0;
      layer2[i][0][k].D_CB01[1] = 0;
      layer2[i][0][k].D_CB10[1] = 0;
      layer2[i][0][k].D_CB11[1] = 0;
      layer2[i][0][k].D_CB02[1] = 0;
      layer2[i][0][k].D_CB20[1] = 0;
      layer2[i][0][k].D_CB12[1] = 0;
      layer2[i][0][k].D_CB21[1] = 0;
      layer2[i][0][k].D_CB22[1] = 0;
    }

  for(j=0; j<5; j++)
    for(k=0; k<5; k++){
      layer2[4][j][k].D_CaMKII[2] = 0;
      layer2[4][j][k].D_Ca[2] = 0;
      layer2[4][j][k].D_CaM00[2] = 0;
	  layer2[4][j][k].D_CaM10[2] = 0;
	  layer2[4][j][k].D_CaM01[2] = 0;
	  layer2[4][j][k].D_CaM11[2] = 0;
	  layer2[4][j][k].D_CaM20[2] = 0;
	  layer2[4][j][k].D_CaM02[2] = 0;
	  layer2[4][j][k].D_CaM12[2] = 0;
	  layer2[4][j][k].D_CaM21[2] = 0;
	  layer2[4][j][k].D_CaM22[2] = 0;      
	  layer2[4][j][k].D_CB00[2] = 0;
	  layer2[4][j][k].D_CB01[2] = 0;
	  layer2[4][j][k].D_CB10[2] = 0;
	  layer2[4][j][k].D_CB11[2] = 0;
	  layer2[4][j][k].D_CB02[2] = 0;
	  layer2[4][j][k].D_CB20[2] = 0;
	  layer2[4][j][k].D_CB12[2] = 0;
	  layer2[4][j][k].D_CB21[2] = 0;
	  layer2[4][j][k].D_CB22[2] = 0;
    }
  
  for(i=0; i<5; i++)
    for(k=0; k<5; k++){
      layer2[i][4][k].D_CaMKII[3] = 0;
      layer2[i][4][k].D_Ca[3] = 0;
      layer2[i][4][k].D_CaM00[3] = 0;
	  layer2[i][4][k].D_CaM10[3] = 0;
	  layer2[i][4][k].D_CaM01[3] = 0;
	  layer2[i][4][k].D_CaM11[3] = 0;
	  layer2[i][4][k].D_CaM02[3] = 0;
	  layer2[i][4][k].D_CaM20[3] = 0;
	  layer2[i][4][k].D_CaM12[3] = 0;
	  layer2[i][4][k].D_CaM21[3] = 0;
	  layer2[i][4][k].D_CaM22[3] = 0;
	  layer2[i][4][k].D_CB00[3] = 0;
	  layer2[i][4][k].D_CB01[3] = 0;
	  layer2[i][4][k].D_CB10[3] = 0;
	  layer2[i][4][k].D_CB11[3] = 0;
	  layer2[i][4][k].D_CB20[3] = 0;
	  layer2[i][4][k].D_CB02[3] = 0;
	  layer2[i][4][k].D_CB12[3] = 0;
	  layer2[i][4][k].D_CB21[3] = 0;
	  layer2[i][4][k].D_CB22[3] = 0;
    }

  for(i=0; i<5; i++)
    for(j=0; j<5; j++){
      layer2[i][j][4].D_CaMKII[5] = 0;
      layer2[i][j][4].D_Ca[5] = 0;
      layer2[i][j][4].D_CaM00[5] = 0;
	  layer2[i][j][4].D_CaM10[5] = 0;
	  layer2[i][j][4].D_CaM01[5] = 0;
	  layer2[i][j][4].D_CaM11[5] = 0;
	  layer2[i][j][4].D_CaM02[5] = 0;
	  layer2[i][j][4].D_CaM20[5] = 0;
	  layer2[i][j][4].D_CaM21[5] = 0;
	  layer2[i][j][4].D_CaM12[5] = 0;
	  layer2[i][j][4].D_CaM22[5] = 0;      
	  layer2[i][j][4].D_CB00[5] = 0;
	  layer2[i][j][4].D_CB01[5] = 0;
	  layer2[i][j][4].D_CB10[5] = 0;
	  layer2[i][j][4].D_CB11[5] = 0;
	  layer2[i][j][4].D_CB02[5] = 0;
	  layer2[i][j][4].D_CB20[5] = 0;
	  layer2[i][j][4].D_CB21[5] = 0;
	  layer2[i][j][4].D_CB12[5] = 0;
	  layer2[i][j][4].D_CB22[5] = 0;
    }
  
  //connection between layer2 to layer3 by voxel [2,2,4]
  layer2[2][2][4].D_Ca[5] = D_Ca_C * (0.1 * 0.1) / (pow(0.1,3) * 0.1);
  layer2[2][2][4].D_CaM00[5] = D_CaM00_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
  layer2[2][2][4].D_CaM01[5] = D_CaM01_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
  layer2[2][2][4].D_CaM10[5] = D_CaM10_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
  layer2[2][2][4].D_CaM11[5] = D_CaM11_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
  layer2[2][2][4].D_CaM20[5] = D_CaM20_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
  layer2[2][2][4].D_CaM02[5] = D_CaM02_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
  layer2[2][2][4].D_CaM21[5] = D_CaM21_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
  layer2[2][2][4].D_CaM12[5] = D_CaM12_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
  layer2[2][2][4].D_CaM22[5] = D_CaM22_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1); 
  layer2[2][2][4].D_CB00[5] = D_CB00 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
  layer2[2][2][4].D_CB01[5] = D_CB01 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
  layer2[2][2][4].D_CB10[5] = D_CB10 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
  layer2[2][2][4].D_CB11[5] = D_CB11 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
  layer2[2][2][4].D_CB20[5] = D_CB20 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
  layer2[2][2][4].D_CB02[5] = D_CB02 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
  layer2[2][2][4].D_CB21[5] = D_CB21 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
  layer2[2][2][4].D_CB12[5] = D_CB12 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
  layer2[2][2][4].D_CB22[5] = D_CB22 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);

  /********************************************************************************/
  //initialize layer 3
  for(i=0; i<8; i++)
    for(j=4; j<6; j++){
      layer3[i].D_Ca[j] = D_Ca_C * (0.1 * 0.1) / (pow(0.1,3) * 0.1);
      layer3[i].D_CaM00[j] = D_CaM00_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer3[i].D_CaM01[j] = D_CaM01_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer3[i].D_CaM10[j] = D_CaM10_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer3[i].D_CaM11[j] = D_CaM11_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer3[i].D_CaM20[j] = D_CaM20_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer3[i].D_CaM02[j] = D_CaM02_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer3[i].D_CaM21[j] = D_CaM21_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer3[i].D_CaM12[j] = D_CaM12_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer3[i].D_CaM22[j] = D_CaM22_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);	  
	  layer3[i].D_CB00[j] = D_CB00 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer3[i].D_CB01[j] = D_CB01 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer3[i].D_CB10[j] = D_CB10 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer3[i].D_CB11[j] = D_CB11 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer3[i].D_CB02[j] = D_CB02 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer3[i].D_CB20[j] = D_CB20 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer3[i].D_CB21[j] = D_CB21 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer3[i].D_CB12[j] = D_CB12 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
	  layer3[i].D_CB22[j] = D_CB22 * (0.1 * 0.1) / (pow(0.1, 3) * 0.1);
    }	
  
  //connection between layer3 to layer 4
  layer3[7].D_Ca[5] = D_Ca_C * (0.1 * 0.1) / (pow(0.1,3) * 0.3);
  layer3[7].D_CaM00[5] = D_CaM00_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.3);
  layer3[7].D_CaM01[5] = D_CaM01_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.3);
  layer3[7].D_CaM10[5] = D_CaM10_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.3);
  layer3[7].D_CaM11[5] = D_CaM11_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.3);
  layer3[7].D_CaM02[5] = D_CaM02_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.3);
  layer3[7].D_CaM20[5] = D_CaM20_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.3);
  layer3[7].D_CaM21[5] = D_CaM21_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.3);
  layer3[7].D_CaM12[5] = D_CaM12_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.3);
  layer3[7].D_CaM22[5] = D_CaM22_C * (0.1 * 0.1) / (pow(0.1, 3) * 0.3);  
  layer3[7].D_CB00[5] = D_CB00 * (0.1 * 0.1) / (pow(0.1, 3) * 0.3);
  layer3[7].D_CB01[5] = D_CB01 * (0.1 * 0.1) / (pow(0.1, 3) * 0.3);
  layer3[7].D_CB10[5] = D_CB10 * (0.1 * 0.1) / (pow(0.1, 3) * 0.3);
  layer3[7].D_CB11[5] = D_CB11 * (0.1 * 0.1) / (pow(0.1, 3) * 0.3);
  layer3[7].D_CB02[5] = D_CB02 * (0.1 * 0.1) / (pow(0.1, 3) * 0.3);
  layer3[7].D_CB20[5] = D_CB20 * (0.1 * 0.1) / (pow(0.1, 3) * 0.3);
  layer3[7].D_CB21[5] = D_CB21 * (0.1 * 0.1) / (pow(0.1, 3) * 0.3);
  layer3[7].D_CB12[5] = D_CB12 * (0.1 * 0.1) / (pow(0.1, 3) * 0.3);
  layer3[7].D_CB22[5] = D_CB22 * (0.1 * 0.1) / (pow(0.1, 3) * 0.3);
 
  /*********************************************************************************/
  //initialize layer 4
  for(i=0; i<4; i++)
    for(j=0; j<2; j++)
      for(k=0; k<2; k++)
	for(m=0; m<6; m++){
	  layer4[i][j][k].D_Ca[m] = D_Ca_C * (0.5 * 0.5) / (pow(0.5,3) * 0.5);
	  layer4[i][j][k].D_CaM00[m] = D_CaM00_C * (0.5 * 0.5) / (pow(0.5, 3) * 0.5);
	  layer4[i][j][k].D_CaM01[m] = D_CaM01_C * (0.5 * 0.5) / (pow(0.5, 3) * 0.5);
	  layer4[i][j][k].D_CaM10[m] = D_CaM10_C * (0.5 * 0.5) / (pow(0.5, 3) * 0.5);
	  layer4[i][j][k].D_CaM11[m] = D_CaM11_C * (0.5 * 0.5) / (pow(0.5, 3) * 0.5);
	  layer4[i][j][k].D_CaM02[m] = D_CaM02_C * (0.5 * 0.5) / (pow(0.5, 3) * 0.5);
	  layer4[i][j][k].D_CaM20[m] = D_CaM20_C * (0.5 * 0.5) / (pow(0.5, 3) * 0.5);
	  layer4[i][j][k].D_CaM12[m] = D_CaM12_C * (0.5 * 0.5) / (pow(0.5, 3) * 0.5);
	  layer4[i][j][k].D_CaM21[m] = D_CaM21_C * (0.5 * 0.5) / (pow(0.5, 3) * 0.5);
	  layer4[i][j][k].D_CaM22[m] = D_CaM22_C * (0.5 * 0.5) / (pow(0.5, 3) * 0.5);	  
	  layer4[i][j][k].D_CB00[m] = D_CB00 * (0.5 * 0.5) / (pow(0.5, 3) * 0.5);
	  layer4[i][j][k].D_CB01[m] = D_CB01 * (0.5 * 0.5) / (pow(0.5, 3) * 0.5);
	  layer4[i][j][k].D_CB10[m] = D_CB10 * (0.5 * 0.5) / (pow(0.5, 3) * 0.5);
	  layer4[i][j][k].D_CB11[m] = D_CB11 * (0.5 * 0.5) / (pow(0.5, 3) * 0.5);
	  layer4[i][j][k].D_CB20[m] = D_CB20 * (0.5 * 0.5) / (pow(0.5, 3) * 0.5);
	  layer4[i][j][k].D_CB02[m] = D_CB02 * (0.5 * 0.5) / (pow(0.5, 3) * 0.5);
	  layer4[i][j][k].D_CB21[m] = D_CB21 * (0.5 * 0.5) / (pow(0.5, 3) * 0.5);
	  layer4[i][j][k].D_CB12[m] = D_CB12 * (0.5 * 0.5) / (pow(0.5, 3) * 0.5);
	  layer4[i][j][k].D_CB22[m] = D_CB22 * (0.5 * 0.5) / (pow(0.5, 3) * 0.5);
	}				
  
  //set boundary
  for(j=0; j<2; j++)
    for(k=0; k<2; k++){
      layer4[0][j][k].D_Ca[0] = 0;
      layer4[0][j][k].D_CaM00[0] = 0;
	  layer4[0][j][k].D_CaM01[0] = 0;
	  layer4[0][j][k].D_CaM10[0] = 0;
	  layer4[0][j][k].D_CaM11[0] = 0;
	  layer4[0][j][k].D_CaM20[0] = 0;
	  layer4[0][j][k].D_CaM02[0] = 0;
	  layer4[0][j][k].D_CaM21[0] = 0;
	  layer4[0][j][k].D_CaM12[0] = 0;
	  layer4[0][j][k].D_CaM22[0] = 0;      
	  layer4[0][j][k].D_CB00[0] = 0;
	  layer4[0][j][k].D_CB01[0] = 0;
	  layer4[0][j][k].D_CB10[0] = 0;
	  layer4[0][j][k].D_CB11[0] = 0;
	  layer4[0][j][k].D_CB20[0] = 0;
	  layer4[0][j][k].D_CB02[0] = 0;
	  layer4[0][j][k].D_CB21[0] = 0;
	  layer4[0][j][k].D_CB12[0] = 0;
	  layer4[0][j][k].D_CB22[0] = 0;
    }

  for(i=0; i<4; i++)
    for(k=0; k<2; k++){
      layer4[i][0][k].D_Ca[1] = 0;
      layer4[i][0][k].D_CaM00[1] = 0;
	  layer4[i][0][k].D_CaM01[1] = 0;
	  layer4[i][0][k].D_CaM10[1] = 0;
	  layer4[i][0][k].D_CaM11[1] = 0;
	  layer4[i][0][k].D_CaM20[1] = 0;
	  layer4[i][0][k].D_CaM02[1] = 0;
	  layer4[i][0][k].D_CaM21[1] = 0;
	  layer4[i][0][k].D_CaM12[1] = 0;
	  layer4[i][0][k].D_CaM22[1] = 0;      
	  layer4[i][0][k].D_CB00[1] = 0;
	  layer4[i][0][k].D_CB01[1] = 0;
	  layer4[i][0][k].D_CB10[1] = 0;
	  layer4[i][0][k].D_CB11[1] = 0;
	  layer4[i][0][k].D_CB02[1] = 0;
	  layer4[i][0][k].D_CB20[1] = 0;
	  layer4[i][0][k].D_CB21[1] = 0;
	  layer4[i][0][k].D_CB12[1] = 0;
	  layer4[i][0][k].D_CB22[1] = 0;
    }
  
  for(j=0; j<2; j++)
    for(k=0; k<2; k++){
      layer4[3][j][k].D_Ca[2] = 0;
      layer4[3][j][k].D_CaM00[2] = 0;
	  layer4[3][j][k].D_CaM01[2] = 0;
	  layer4[3][j][k].D_CaM10[2] = 0;
	  layer4[3][j][k].D_CaM11[2] = 0;
	  layer4[3][j][k].D_CaM20[2] = 0;
	  layer4[3][j][k].D_CaM02[2] = 0;
	  layer4[3][j][k].D_CaM21[2] = 0;
	  layer4[3][j][k].D_CaM12[2] = 0;
	  layer4[3][j][k].D_CaM22[2] = 0;   
	  layer4[3][j][k].D_CB00[2] = 0;
	  layer4[3][j][k].D_CB01[2] = 0;
	  layer4[3][j][k].D_CB10[2] = 0;
	  layer4[3][j][k].D_CB11[2] = 0;
	  layer4[3][j][k].D_CB02[2] = 0;
	  layer4[3][j][k].D_CB20[2] = 0;
	  layer4[3][j][k].D_CB12[2] = 0;
	  layer4[3][j][k].D_CB21[2] = 0;
	  layer4[3][j][k].D_CB22[2] = 0;
    }
  
  for(i=0; i<4; i++)
    for(k=0; k<2; k++){
      layer4[i][1][k].D_Ca[3] = 0;
      layer4[i][1][k].D_CaM00[3] = 0;
	  layer4[i][1][k].D_CaM01[3] = 0;
	  layer4[i][1][k].D_CaM10[3] = 0;
	  layer4[i][1][k].D_CaM11[3] = 0;
	  layer4[i][1][k].D_CaM20[3] = 0;
	  layer4[i][1][k].D_CaM02[3] = 0;
	  layer4[i][1][k].D_CaM12[3] = 0;
	  layer4[i][1][k].D_CaM21[3] = 0;
	  layer4[i][1][k].D_CaM22[3] = 0;      
	  layer4[i][1][k].D_CB00[3] = 0;
	  layer4[i][1][k].D_CB01[3] = 0;
	  layer4[i][1][k].D_CB10[3] = 0;
	  layer4[i][1][k].D_CB11[3] = 0;
	  layer4[i][1][k].D_CB02[3] = 0;
	  layer4[i][1][k].D_CB20[3] = 0;
	  layer4[i][1][k].D_CB21[3] = 0;
	  layer4[i][1][k].D_CB12[3] = 0;
	  layer4[i][1][k].D_CB22[3] = 0;
    }

  for(i=0; i<4; i++)
    for(j=0; j<2; j++){
      layer4[i][j][0].D_Ca[4] = 0;
      layer4[i][j][0].D_CaM00[4] = 0;
	  layer4[i][j][0].D_CaM01[4] = 0;
	  layer4[i][j][0].D_CaM10[4] = 0;
	  layer4[i][j][0].D_CaM11[4] = 0;
	  layer4[i][j][0].D_CaM20[4] = 0;
	  layer4[i][j][0].D_CaM02[4] = 0;
	  layer4[i][j][0].D_CaM21[4] = 0;
	  layer4[i][j][0].D_CaM12[4] = 0;
	  layer4[i][j][0].D_CaM22[4] = 0;     
	  layer4[i][j][0].D_CB00[4] = 0;
	  layer4[i][j][0].D_CB01[4] = 0;
	  layer4[i][j][0].D_CB10[4] = 0;
	  layer4[i][j][0].D_CB11[4] = 0;
	  layer4[i][j][0].D_CB02[4] = 0;
	  layer4[i][j][0].D_CB20[4] = 0;
	  layer4[i][j][0].D_CB12[4] = 0;
	  layer4[i][j][0].D_CB21[4] = 0;
	  layer4[i][j][0].D_CB22[4] = 0;
    }
  
  for(i=0; i<4; i++)
    for(j=0; j<2; j++){
      layer4[i][j][1].D_Ca[5] = 0;
      layer4[i][j][1].D_CaM00[5] = 0;
	  layer4[i][j][1].D_CaM01[5] = 0;
	  layer4[i][j][1].D_CaM10[5] = 0;
	  layer4[i][j][1].D_CaM11[5] = 0;
	  layer4[i][j][1].D_CaM02[5] = 0;
	  layer4[i][j][1].D_CaM20[5] = 0;
	  layer4[i][j][1].D_CaM21[5] = 0;
	  layer4[i][j][1].D_CaM12[5] = 0;
	  layer4[i][j][1].D_CaM22[5] = 0;      
	  layer4[i][j][1].D_CB00[5] = 0;
	  layer4[i][j][1].D_CB01[5] = 0;
	  layer4[i][j][1].D_CB10[5] = 0;
	  layer4[i][j][1].D_CB11[5] = 0;
	  layer4[i][j][1].D_CB02[5] = 0;
	  layer4[i][j][1].D_CB20[5] = 0;
	  layer4[i][j][1].D_CB12[5] = 0;
	  layer4[i][j][1].D_CB21[5] = 0;
	  layer4[i][j][1].D_CB22[5] = 0;
    }
  
  //connection between layer 4 and layer 3 by voxel[2,0,0]
  layer4[2][0][0].D_Ca[4] = D_Ca_C * (0.1 * 0.1) / (pow(0.5,3) * 0.3);
  layer4[2][0][0].D_CaM00[4] = D_CaM00_C * (0.1 * 0.1) / (pow(0.5, 3) * 0.3);
  layer4[2][0][0].D_CaM10[4] = D_CaM10_C * (0.1 * 0.1) / (pow(0.5, 3) * 0.3);
  layer4[2][0][0].D_CaM01[4] = D_CaM01_C * (0.1 * 0.1) / (pow(0.5, 3) * 0.3);
  layer4[2][0][0].D_CaM11[4] = D_CaM11_C * (0.1 * 0.1) / (pow(0.5, 3) * 0.3);
  layer4[2][0][0].D_CaM20[4] = D_CaM20_C * (0.1 * 0.1) / (pow(0.5, 3) * 0.3);
  layer4[2][0][0].D_CaM02[4] = D_CaM02_C * (0.1 * 0.1) / (pow(0.5, 3) * 0.3);
  layer4[2][0][0].D_CaM21[4] = D_CaM21_C * (0.1 * 0.1) / (pow(0.5, 3) * 0.3);
  layer4[2][0][0].D_CaM12[4] = D_CaM12_C * (0.1 * 0.1) / (pow(0.5, 3) * 0.3);
  layer4[2][0][0].D_CaM22[4] = D_CaM22_C * (0.1 * 0.1) / (pow(0.5, 3) * 0.3); 
  layer4[2][0][0].D_CB00[4] = D_CB00 * (0.1 * 0.1) / (pow(0.5, 3) * 0.3);
  layer4[2][0][0].D_CB01[4] = D_CB01 * (0.1 * 0.1) / (pow(0.5, 3) * 0.3);
  layer4[2][0][0].D_CB10[4] = D_CB10 * (0.1 * 0.1) / (pow(0.5, 3) * 0.3);
  layer4[2][0][0].D_CB11[4] = D_CB11 * (0.1 * 0.1) / (pow(0.5, 3) * 0.3);
  layer4[2][0][0].D_CB02[4] = D_CB02 * (0.1 * 0.1) / (pow(0.5, 3) * 0.3);
  layer4[2][0][0].D_CB20[4] = D_CB20 * (0.1 * 0.1) / (pow(0.5, 3) * 0.3);
  layer4[2][0][0].D_CB12[4] = D_CB12 * (0.1 * 0.1) / (pow(0.5, 3) * 0.3);
  layer4[2][0][0].D_CB21[4] = D_CB21 * (0.1 * 0.1) / (pow(0.5, 3) * 0.3);
  layer4[2][0][0].D_CB22[4] = D_CB22 * (0.1 * 0.1) / (pow(0.5, 3) * 0.3);
   
  /***************************************************************************************/
  //initialize the concentration
  for(i=0; i<10; i++)
    for(j=0; j<10; j++){
      layer1[i][j].CaMKII_num = 0;
      layer1[i][j].Ca_pump = 0;
      layer1[i][j].Ca = 0;
      layer1[i][j].CaM00 = 0;
//      layer1[i][j].CaM00 = 1;
	  layer1[i][j].CaM10 = 0;
	  layer1[i][j].CaM01 = 0;
	  layer1[i][j].CaM11 = 0;
	  layer1[i][j].CaM02 = 0;
	  layer1[i][j].CaM20 = 0;
	  layer1[i][j].CaM21 = 0;
	  layer1[i][j].CaM12 = 0;
	  layer1[i][j].CaM22 = 0;
      layer1[i][j].CaN = 2;
          layer1[i][j].CaNCaM00 = 0;
	  layer1[i][j].CaNCaM01 = 0;
	  layer1[i][j].CaNCaM10 = 0;
	  layer1[i][j].CaNCaM11 = 0;
	  layer1[i][j].CaNCaM02 = 0;
	  layer1[i][j].CaNCaM20 = 0;
	  layer1[i][j].CaNCaM21 = 0;
	  layer1[i][j].CaNCaM12 = 0;
	  layer1[i][j].CaNCaM22 = 0;      
      layer1[i][j].Ng = 1;
          layer1[i][j].NgCaM00 = 1;
	  layer1[i][j].NgCaM01 = 0;
	  layer1[i][j].NgCaM10 = 0;
	  layer1[i][j].NgCaM11 = 0;
	  layer1[i][j].NgCaM02 = 0;
	  layer1[i][j].NgCaM20 = 0;
	  layer1[i][j].NgCaM21 = 0;
	  layer1[i][j].NgCaM12 = 0;
	  layer1[i][j].NgCaM22 = 0;      
       layer1[i][j].CB00 =3;
	  layer1[i][j].CB01 =0;
	  layer1[i][j].CB10 =0;
	  layer1[i][j].CB11 =0;
	  layer1[i][j].CB20 =0;
	  layer1[i][j].CB02 =0;
	  layer1[i][j].CB12 =0;
	  layer1[i][j].CB21 =0;
	  layer1[i][j].CB22 =0;
    }
  
  for(i=0; i<5; i++)
    for(j=0; j<5; j++)
      for(k=0; k<5; k++){
	layer2[i][j][k].CaMKII_num = 0;
	layer2[i][j][k].Ca_pump = 0;
	layer2[i][j][k].Ca = 0;
	layer2[i][j][k].CaM00 = 1;
//	layer2[i][j][k].CaM00 = 8;
        layer2[i][j][k].CaM10 = 0;
	layer2[i][j][k].CaM01 = 0;
	layer2[i][j][k].CaM11 = 0;
	layer2[i][j][k].CaM02 = 0;
	layer2[i][j][k].CaM20 = 0;
	layer2[i][j][k].CaM21 = 0;
	layer2[i][j][k].CaM12 = 0;
	layer2[i][j][k].CaM22 = 0;
    layer2[i][j][k].CaN = 0;
        layer2[i][j][k].CaNCaM00 = 0;
	layer2[i][j][k].CaNCaM01 = 0;
	layer2[i][j][k].CaNCaM10 = 0;
	layer2[i][j][k].CaNCaM11 = 0;
	layer2[i][j][k].CaNCaM02 = 0;
	layer2[i][j][k].CaNCaM20 = 0;
	layer2[i][j][k].CaNCaM21 = 0;
	layer2[i][j][k].CaNCaM12 = 0;
        layer2[i][j][k].CaNCaM22 = 0;
    layer2[i][j][k].Ng = 9;
        layer2[i][j][k].NgCaM00 = 7;
	layer2[i][j][k].NgCaM01 = 0;
	layer2[i][j][k].NgCaM10 = 0;
	layer2[i][j][k].NgCaM11 = 0;
	layer2[i][j][k].NgCaM02 = 0;
	layer2[i][j][k].NgCaM20 = 0;
	layer2[i][j][k].NgCaM21 = 0;
	layer2[i][j][k].NgCaM12 = 0;
        layer2[i][j][k].NgCaM22 = 0;
     layer2[i][j][k].CB00 = 20;
	layer2[i][j][k].CB01 = 1;
	layer2[i][j][k].CB10 = 3;
	layer2[i][j][k].CB11 = 0;
	layer2[i][j][k].CB02 = 0;
	layer2[i][j][k].CB20 = 0;
	layer2[i][j][k].CB21 = 0;
	layer2[i][j][k].CB12 = 0;
	layer2[i][j][k].CB22 = 0;
      }
  
  for(i=0; i<8; i++){
    layer3[i].CaMKII_num = 0;
    layer3[i].Ca_pump = 0;
    layer3[i].Ca = 0;
    layer3[i].CaM00 = 1;
//    layer3[i].CaM00 = 8;
        layer3[i].CaM10 = 0;
	layer3[i].CaM01 = 0;
	layer3[i].CaM11 = 0;
	layer3[i].CaM02 = 0;
	layer3[i].CaM20 = 0;
	layer3[i].CaM21 = 0;
	layer3[i].CaM12 = 0;
	layer3[i].CaM22 = 0;
    layer3[i].CaN = 0;
        layer3[i].CaNCaM00 = 0;
	layer3[i].CaNCaM01 = 0;
	layer3[i].CaNCaM10 = 0;
	layer3[i].CaNCaM11 = 0;
	layer3[i].CaNCaM02 = 0;
	layer3[i].CaNCaM20 = 0;
	layer3[i].CaNCaM21 = 0;
	layer3[i].CaNCaM12 = 0;
        layer3[i].CaNCaM22 = 0;
    layer3[i].Ng = 9;
        layer3[i].NgCaM00 = 7;
	layer3[i].NgCaM01 = 0;
	layer3[i].NgCaM10 = 0;
	layer3[i].NgCaM11 = 0;
	layer3[i].NgCaM02 = 0;
	layer3[i].NgCaM20 = 0;
	layer3[i].NgCaM21 = 0;
	layer3[i].NgCaM12 = 0;
        layer3[i].NgCaM22 = 0;
     layer3[i].CB00 = 20;
	layer3[i].CB01 = 1;
	layer3[i].CB10 = 3;
	layer3[i].CB11 = 0;
	layer3[i].CB02 = 0;
	layer3[i].CB20 = 0;
	layer3[i].CB12 = 0;
	layer3[i].CB21 = 0;
	layer3[i].CB22 = 0;
  }
  
  for(i=0; i<4; i++)
    for(j=0; j<2; j++)
      for(k=0; k<2; k++){
	layer4[i][j][k].CaMKII_num = 0;
	layer4[i][j][k].Ca_pump = 0;
	layer4[i][j][k].Ca = 4;
	layer4[i][j][k].CaM00 = 126;
//	layer4[i][j][k].CaM00 = 1000;
        layer4[i][j][k].CaM10 = 0;
	layer4[i][j][k].CaM01 = 0;
	layer4[i][j][k].CaM11 = 0;
	layer4[i][j][k].CaM02 = 0;
	layer4[i][j][k].CaM20 = 0;
	layer4[i][j][k].CaM21 = 0;
	layer4[i][j][k].CaM12 = 0;
	layer4[i][j][k].CaM22 = 0;
    layer4[i][j][k].CaN = 0;
        layer4[i][j][k].CaNCaM00 = 0;
	layer4[i][j][k].CaNCaM01 = 0;
	layer4[i][j][k].CaNCaM10 = 0;
	layer4[i][j][k].CaNCaM11 = 0;
	layer4[i][j][k].CaNCaM02 = 0;
	layer4[i][j][k].CaNCaM20 = 0;
	layer4[i][j][k].CaNCaM21 = 0;
	layer4[i][j][k].CaNCaM12 = 0;
        layer4[i][j][k].CaNCaM22 = 0;
    layer4[i][j][k].Ng = 1121;
        layer4[i][j][k].NgCaM00 = 874;
	layer4[i][j][k].NgCaM01 = 0;
	layer4[i][j][k].NgCaM10 = 0;
	layer4[i][j][k].NgCaM11 = 0;
	layer4[i][j][k].NgCaM02 = 0;
	layer4[i][j][k].NgCaM20 = 0;
	layer4[i][j][k].NgCaM21 = 0;
	layer4[i][j][k].NgCaM12 = 0;
        layer4[i][j][k].NgCaM22 = 0;
     layer4[i][j][k].CB00 = 2523;
	layer4[i][j][k].CB01 = 153;
	layer4[i][j][k].CB10 = 267;
	layer4[i][j][k].CB11 = 16;
	layer4[i][j][k].CB02 = 10;
	layer4[i][j][k].CB20 = 28;
	layer4[i][j][k].CB21 = 2;
	layer4[i][j][k].CB12 = 1;
	layer4[i][j][k].CB22 = 0;
      }
  
  /*************************************************************/
  //initialize transition rate
  for(i=0; i<10; i++)
    for(j=0; j<10; j++){
		layer1[i][j].CaM00_to_CaM10 = 1.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
		layer1[i][j].CaM10_to_CaM00 = 10.0 / 1000.0;
		layer1[i][j].CaM00_to_CaM01 = 20.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
		layer1[i][j].CaM01_to_CaM00 = 1000.0 / 1000.0;
		layer1[i][j].CaM10_to_CaM11 = 20.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
		layer1[i][j].CaM11_to_CaM10 = 1000.0 / 1000.0;
		layer1[i][j].CaM10_to_CaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
		layer1[i][j].CaM20_to_CaM10 = 12.0 / 1000.0;
		layer1[i][j].CaM01_to_CaM11 = 1.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
		layer1[i][j].CaM11_to_CaM01 = 10.0 / 1000.0;
		layer1[i][j].CaM01_to_CaM02 = 100.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
		layer1[i][j].CaM02_to_CaM01 = 1000.0 / 1000.0;
		layer1[i][j].CaM11_to_CaM21 = 10.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
		layer1[i][j].CaM21_to_CaM11 = 12.0 / 1000.0;
		layer1[i][j].CaM11_to_CaM12 = 100.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
		layer1[i][j].CaM12_to_CaM11 = 1000.0 / 1000.0;
		layer1[i][j].CaM20_to_CaM21 = 20.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
		layer1[i][j].CaM21_to_CaM20 = 1000.0 / 1000.0;
		layer1[i][j].CaM02_to_CaM12 = 1.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
		layer1[i][j].CaM12_to_CaM02 = 10.0 / 1000.0;
		layer1[i][j].CaM21_to_CaM22 = 100.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
		layer1[i][j].CaM22_to_CaM21 = 1000.0 / 1000.0;
		layer1[i][j].CaM12_to_CaM22 = 10.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
		layer1[i][j].CaM22_to_CaM12 = 12.0 / 1000.0;

	layer1[i][j].CaM00_to_CaMKIICaM00 = 0.2 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM00_to_CaM00 = 90.0 / 1000.0;
        layer1[i][j].CaM10_to_CaMKIICaM10 = 2.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM10_to_CaM10 = 90.0 / 1000.0;
        layer1[i][j].CaM01_to_CaMKIICaM01 = 2.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM01_to_CaM01 = 90.0 / 1000.0;
        layer1[i][j].CaM11_to_CaMKIICaM11 = 10.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM11_to_CaM11 = 45.0 / 1000.0;
        layer1[i][j].CaM20_to_CaMKIICaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM20_to_CaM20 = 45.0 / 1000.0;
        layer1[i][j].CaM02_to_CaMKIICaM02 = 10.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM02_to_CaM02 = 45.0 / 1000.0;
        layer1[i][j].CaM12_to_CaMKIICaM12 = 20.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM12_to_CaM12 = 9.0 / 1000.0;
        layer1[i][j].CaM21_to_CaMKIICaM21 = 20.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM21_to_CaM21 = 9.0 / 1000.0;
        layer1[i][j].CaM22_to_CaMKIICaM22 = 50.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM22_to_CaM22 = 2.25 / 1000.0;

	layer1[i][j].CaMKIICaM00_to_CaMKIICaM10 = 0.48 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM10_to_CaMKIICaM00 = 0.48 / 1000.0; 
        layer1[i][j].CaMKIICaM00_to_CaMKIICaM01 = 1.6 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM01_to_CaMKIICaM00 = 8.0 / 1000.0;
        layer1[i][j].CaMKIICaM10_to_CaMKIICaM11 = 1.6 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM11_to_CaMKIICaM10 = 8.0 / 1000.0;
        layer1[i][j].CaMKIICaM10_to_CaMKIICaM20 = 66.6 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM20_to_CaMKIICaM10 = 8.0 / 1000.0;
        layer1[i][j].CaMKIICaM01_to_CaMKIICaM11 = 0.48 / (6.02 * pow(10, 5) * pow(0.05, 3));
	layer1[i][j].CaMKIICaM11_to_CaMKIICaM01 = 0.48 / 1000.0;
        layer1[i][j].CaMKIICaM01_to_CaMKIICaM02 = 50.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM02_to_CaMKIICaM01 = 50.0 / 1000.0;
        layer1[i][j].CaMKIICaM11_to_CaMKIICaM21 = 66.6 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM21_to_CaMKIICaM11 = 8.0 / 1000.0;
        layer1[i][j].CaMKIICaM11_to_CaMKIICaM12 = 50.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM12_to_CaMKIICaM11 = 50.0 / 1000.0;
        layer1[i][j].CaMKIICaM20_to_CaMKIICaM21 = 1.6 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM21_to_CaMKIICaM20 = 8.0 / 1000.0;
        layer1[i][j].CaMKIICaM02_to_CaMKIICaM12 = 0.48 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM12_to_CaMKIICaM02 = 0.48 / 1000.0;
        layer1[i][j].CaMKIICaM21_to_CaMKIICaM22 = 50.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM22_to_CaMKIICaM21 = 50.0 / 1000.0;
        layer1[i][j].CaMKIICaM12_to_CaMKIICaM22 = 66.6 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaMKIICaM22_to_CaMKIICaM12 = 8.0 / 1000.0;

        layer1[i][j].CaM00_to_CaNCaM00 = 1.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM00_to_CaM00 = 1167.0 / 1000.0;
        layer1[i][j].CaM10_to_CaNCaM10 = 1.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM10_to_CaM10 = 11.67 / 1000.0;
        layer1[i][j].CaM01_to_CaNCaM01 = 1.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM01_to_CaM01 = 7.0 / 1000.0;
        layer1[i][j].CaM11_to_CaNCaM11 = 10.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM11_to_CaM11 = 0.7 / 1000.0;
        layer1[i][j].CaM20_to_CaNCaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM20_to_CaM20 = 5.83 / 1000.0;
        layer1[i][j].CaM02_to_CaNCaM02 = 10.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM02_to_CaM02 = 0.56 / 1000.0;
        layer1[i][j].CaM12_to_CaNCaM12 = 46.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM12_to_CaM12 = 0.02576 / 1000.0;
        layer1[i][j].CaM21_to_CaNCaM21 = 46.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM21_to_CaM21 = 0.161 / 1000.0;
        layer1[i][j].CaM22_to_CaNCaM22 = 46.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM22_to_CaM22 = 0.0013 / 1000.0;

        layer1[i][j].CaNCaM00_to_CaNCaM10 = 5.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM10_to_CaNCaM00 = 0.5 / 1000.0;
        layer1[i][j].CaNCaM00_to_CaNCaM01 = 20.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM01_to_CaNCaM00 = 6.0 / 1000.0;
        layer1[i][j].CaNCaM10_to_CaNCaM11 = 20.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM11_to_CaNCaM10 = 6.0 / 1000.0;
        layer1[i][j].CaNCaM10_to_CaNCaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM20_to_CaNCaM10 = 0.6 / 1000.0;
        layer1[i][j].CaNCaM01_to_CaNCaM11 = 5.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM11_to_CaNCaM01 = 0.5 / 1000.0;
        layer1[i][j].CaNCaM01_to_CaNCaM02 = 50.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM02_to_CaNCaM01 = 4.0 / 1000.0;
        layer1[i][j].CaNCaM11_to_CaNCaM21 = 10.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM21_to_CaNCaM11 = 0.6 / 1000.0;
        layer1[i][j].CaNCaM11_to_CaNCaM12 = 50.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM12_to_CaNCaM11 = 4.0 / 1000.0;
        layer1[i][j].CaNCaM20_to_CaNCaM21 = 20.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM21_to_CaNCaM20 = 6.0 / 1000.0;
        layer1[i][j].CaNCaM02_to_CaNCaM12 = 5.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM12_to_CaNCaM02 = 0.5 / 1000.0;
        layer1[i][j].CaNCaM21_to_CaNCaM22 = 50.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM22_to_CaNCaM21 = 4.0 / 1000.0;
        layer1[i][j].CaNCaM12_to_CaNCaM22 = 10.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].CaNCaM22_to_CaNCaM12 = 0.6 / 1000.0;

        layer1[i][j].CaM00_to_NgCaM00 = 20.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM00_to_CaM00 = 30.0 / 1000.0;
        layer1[i][j].CaM10_to_NgCaM10 = 3.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM10_to_CaM10 = 30.0 / 1000.0;
        layer1[i][j].CaM01_to_NgCaM01 = 20.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM01_to_CaM01 = 30.0 / 1000.0;
        layer1[i][j].CaM11_to_NgCaM11 = 3.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM11_to_CaM11 = 30.0 / 1000.0;
        layer1[i][j].CaM20_to_NgCaM20 = 0.1 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM20_to_CaM20 = 40.0 / 1000.0;
        layer1[i][j].CaM02_to_NgCaM02 = 20.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM02_to_CaM02 = 30.0 / 1000.0;
        layer1[i][j].CaM12_to_NgCaM12 = 3.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM12_to_CaM12 = 30.0 / 1000.0;
        layer1[i][j].CaM21_to_NgCaM21 = 0.1 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM21_to_CaM21 = 40.0 / 1000.0;
        layer1[i][j].CaM22_to_NgCaM22 = 0.1 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM22_to_CaM22 = 40.0 / 1000.0;

        layer1[i][j].NgCaM00_to_NgCaM10 = 1.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM10_to_NgCaM00 = 60.0 / 1000.0;
        layer1[i][j].NgCaM00_to_NgCaM01 = 20.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM01_to_NgCaM00 = 1000.0 / 1000.0;
        layer1[i][j].NgCaM10_to_NgCaM11 = 20.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM11_to_NgCaM10 = 1000.0 / 1000.0;
        layer1[i][j].NgCaM10_to_NgCaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM20_to_NgCaM10 = 480.0 / 1000.0;
        layer1[i][j].NgCaM01_to_NgCaM11 = 1.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM11_to_NgCaM01 = 60.0 / 1000.0;
        layer1[i][j].NgCaM01_to_NgCaM02 = 100.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM02_to_NgCaM01 = 1000.0 / 1000.0;
        layer1[i][j].NgCaM11_to_NgCaM21 = 10.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM21_to_NgCaM11 = 480.0 / 1000.0;
        layer1[i][j].NgCaM11_to_NgCaM12 = 100.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM12_to_NgCaM11 = 1000.0 / 1000.0;
        layer1[i][j].NgCaM20_to_NgCaM21 = 20.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM21_to_NgCaM20 = 1000.0 / 1000.0;
        layer1[i][j].NgCaM02_to_NgCaM12 = 1.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM12_to_NgCaM02 = 60.0 / 1000.0;
        layer1[i][j].NgCaM21_to_NgCaM22 = 100.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM22_to_NgCaM21 = 1000.0 / 1000.0;
        layer1[i][j].NgCaM12_to_NgCaM22 = 10.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].NgCaM22_to_NgCaM12 = 480.0 / 1000.0;

        layer1[i][j].CaMKIICaM00_to_Trapped00 = 10.0 / 1000.0;
        layer1[i][j].Trapped00_to_CaMKIICaM00 = 0.003 / 1000.0;
        layer1[i][j].CaMKIICaM10_to_Trapped10 = 10.0 / 1000.0;
        layer1[i][j].Trapped10_to_CaMKIICaM10 = 0.003 / 1000.0;
        layer1[i][j].CaMKIICaM01_to_Trapped01 = 10.0 / 1000.0;
        layer1[i][j].Trapped01_to_CaMKIICaM01 = 0.003 / 1000.0;
        layer1[i][j].CaMKIICaM11_to_Trapped11 = 10.0 / 1000.0;
        layer1[i][j].Trapped11_to_CaMKIICaM11 = 0.003 / 1000.0;
        layer1[i][j].CaMKIICaM20_to_Trapped20 = 10.0 / 1000.0;
        layer1[i][j].Trapped20_to_CaMKIICaM20 = 0.003 / 1000.0;
        layer1[i][j].CaMKIICaM02_to_Trapped02 = 10.0 / 1000.0;
        layer1[i][j].Trapped02_to_CaMKIICaM02 = 0.003 / 1000.0;
        layer1[i][j].CaMKIICaM12_to_Trapped12 = 10.0 / 1000.0;
        layer1[i][j].Trapped12_to_CaMKIICaM12 = 0.003 / 1000.0;
        layer1[i][j].CaMKIICaM21_to_Trapped21 = 10.0 / 1000.0;
        layer1[i][j].Trapped21_to_CaMKIICaM21 = 0.003 / 1000.0;
        layer1[i][j].CaMKIICaM22_to_Trapped22 = 10.0 / 1000.0;
        layer1[i][j].Trapped22_to_CaMKIICaM22 = 0.003 / 1000.0;
  
        layer1[i][j].Trapped00_to_Auton = 0.09 / 1000.0;
        layer1[i][j].Auton_to_Trapped00 = 0.2 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped10_to_Auton = 0.09 / 1000.0;
        layer1[i][j].Auton_to_Trapped10 = 2.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped01_to_Auton = 0.09 / 1000.0;
        layer1[i][j].Auton_to_Trapped01 = 2.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped11_to_Auton = 0.045 / 1000.0;
        layer1[i][j].Auton_to_Trapped11 = 10.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped20_to_Auton = 0.045 / 1000.0;
        layer1[i][j].Auton_to_Trapped20 = 10.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped02_to_Auton = 0.045 / 1000.0;
        layer1[i][j].Auton_to_Trapped02 = 10.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped12_to_Auton = 0.009 / 1000.0;
        layer1[i][j].Auton_to_Trapped12 = 20 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped21_to_Auton = 0.009 / 1000.0;
        layer1[i][j].Auton_to_Trapped21 = 20 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped22_to_Auton = 0.00225 / 1000.0;
        layer1[i][j].Auton_to_Trapped22 = 50 / (6.02 * pow(10, 5) * pow(0.05, 3));

        layer1[i][j].Trapped00_to_Trapped10 = 0.4 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped10_to_Trapped00 = 0.4 / 1000.0;
        layer1[i][j].Trapped00_to_Trapped01 = 6.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped01_to_Trapped00 = 1.2 / 1000.0;
        layer1[i][j].Trapped10_to_Trapped11 = 6.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped11_to_Trapped10 = 1.2 / 1000.0;
        layer1[i][j].Trapped10_to_Trapped20 = 3.3 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped20_to_Trapped10 = 0.4 / 1000.0;
        layer1[i][j].Trapped01_to_Trapped11 = 0.4 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped11_to_Trapped01 = 0.4 / 1000.0;
        layer1[i][j].Trapped01_to_Trapped02 = 11.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped02_to_Trapped01 = 11.0 / 1000.0;
        layer1[i][j].Trapped11_to_Trapped21 = 3.3 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped21_to_Trapped11 = 0.4 / 1000.0;
        layer1[i][j].Trapped11_to_Trapped12 = 11.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped12_to_Trapped11 = 11.0 / 1000.0;
        layer1[i][j].Trapped20_to_Trapped21 = 6.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped21_to_Trapped20 = 1.2 / 1000.0;
        layer1[i][j].Trapped02_to_Trapped12 = 0.4 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped12_to_Trapped02 = 0.4 / 1000.0;
        layer1[i][j].Trapped21_to_Trapped22 = 11.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped22_to_Trapped21 = 11.0 / 1000.0;
        layer1[i][j].Trapped12_to_Trapped22 = 3.3  / (6.02 * pow(10, 5) * pow(0.05, 3));
        layer1[i][j].Trapped22_to_Trapped12 = 0.4 / 1000.0;

        layer1[i][j].Auton_to_CaMKII = 0.003 / 1000.0;
        layer1[i][j].Auton_to_Capped = 0.1 / 1000.0;
        layer1[i][j].Capped_to_Auton = 0.01 / 1000.0;   
		

	  layer1[i][j].CB00_to_CB10 = 5.5 / (6.02 * pow(10, 5) * pow(0.05, 3));
	  layer1[i][j].CB10_to_CB00 = 2.6 / 1000;
	  layer1[i][j].CB00_to_CB01 = 43.5 / (6.02 * pow(10, 5) * pow(0.05, 3));
	  layer1[i][j].CB01_to_CB00 = 35.8 / 1000;
	  layer1[i][j].CB10_to_CB11 = 43.5 / (6.02 * pow(10, 5) * pow(0.05, 3));
          layer1[i][j].CB11_to_CB10 = 35.8 / 1000;
	  layer1[i][j].CB10_to_CB20 = 5.5 / (6.02 * pow(10, 5) * pow(0.05, 3));
	  layer1[i][j].CB20_to_CB10 = 2.6 / 1000;
	  layer1[i][j].CB01_to_CB11 = 5.5 / (6.02 * pow(10, 5) * pow(0.05, 3));
	  layer1[i][j].CB11_to_CB01 = 2.6 / 1000;
	  layer1[i][j].CB01_to_CB02 = 43.5 / (6.02 * pow(10, 5) * pow(0.05, 3));
	  layer1[i][j].CB02_to_CB01 = 35.8 / 1000;
	  layer1[i][j].CB11_to_CB21 = 5.5 / (6.02 * pow(10, 5) * pow(0.05, 3));
	  layer1[i][j].CB21_to_CB11 = 2.6 / 1000;
	  layer1[i][j].CB11_to_CB12 = 43.5 / (6.02 * pow(10, 5) * pow(0.05, 3));
	  layer1[i][j].CB12_to_CB11 = 35.8 / 1000;
	  layer1[i][j].CB20_to_CB21 = 43.5 / (6.02 * pow(10, 5) * pow(0.05, 3));
	  layer1[i][j].CB21_to_CB20 = 35.8 / 1000;
	  layer1[i][j].CB02_to_CB12 = 5.5 / (6.02 * pow(10, 5) * pow(0.05, 3));
	  layer1[i][j].CB12_to_CB02 = 2.6 / 1000;
	  layer1[i][j].CB12_to_CB22 = 5.5 / (6.02 * pow(10, 5) * pow(0.05, 3));
	  layer1[i][j].CB22_to_CB12 = 2.6 / 1000;
	  layer1[i][j].CB21_to_CB22 = 43.5 / (6.02 * pow(10, 5) * pow(0.05, 3));
	  layer1[i][j].CB22_to_CB21 = 35.8 / 1000;
    }
  
  for(i=0; i<5; i++)
    for(j=0; j<5; j++)
      for(k=0; k<5; k++){
		layer2[i][j][k].CaM00_to_CaM10 = 1.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer2[i][j][k].CaM10_to_CaM00 = 10.0 / 1000.0;
		layer2[i][j][k].CaM00_to_CaM01 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer2[i][j][k].CaM01_to_CaM00 = 1000.0 / 1000.0;
		layer2[i][j][k].CaM10_to_CaM11 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer2[i][j][k].CaM11_to_CaM10 = 1000.0 / 1000.0;
		layer2[i][j][k].CaM10_to_CaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer2[i][j][k].CaM20_to_CaM10 = 12.0 / 1000.0;
		layer2[i][j][k].CaM01_to_CaM11 = 1.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer2[i][j][k].CaM11_to_CaM01 = 10.0 / 1000.0;
		layer2[i][j][k].CaM01_to_CaM02 = 100.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer2[i][j][k].CaM02_to_CaM01 = 1000.0 / 1000.0;
		layer2[i][j][k].CaM11_to_CaM21 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer2[i][j][k].CaM21_to_CaM11 = 12.0 / 1000.0;
		layer2[i][j][k].CaM11_to_CaM12 = 100.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer2[i][j][k].CaM12_to_CaM11 = 1000.0 / 1000.0;
		layer2[i][j][k].CaM20_to_CaM21 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer2[i][j][k].CaM21_to_CaM20 = 1000.0 / 1000.0;
		layer2[i][j][k].CaM02_to_CaM12 = 1.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer2[i][j][k].CaM12_to_CaM02 = 10.0 / 1000.0;
		layer2[i][j][k].CaM21_to_CaM22 = 100.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer2[i][j][k].CaM22_to_CaM21 = 1000.0 / 1000.0;
		layer2[i][j][k].CaM12_to_CaM22 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer2[i][j][k].CaM22_to_CaM12 = 12.0 / 1000.0;

	layer2[i][j][k].CaM00_to_CaMKIICaM00 = 0.2 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM00_to_CaM00 = 90.0 / 1000.0;
        layer2[i][j][k].CaM10_to_CaMKIICaM10 = 2.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM10_to_CaM10 = 90.0 / 1000.0;
        layer2[i][j][k].CaM01_to_CaMKIICaM01 = 2.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM01_to_CaM01 = 90.0 / 1000.0;
        layer2[i][j][k].CaM11_to_CaMKIICaM11 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM11_to_CaM11 = 45.0 / 1000.0;
        layer2[i][j][k].CaM20_to_CaMKIICaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM20_to_CaM20 = 45.0 / 1000.0;
        layer2[i][j][k].CaM02_to_CaMKIICaM02 = 10.0/ (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM02_to_CaM02 = 45.0 / 1000.0;
        layer2[i][j][k].CaM12_to_CaMKIICaM12 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM12_to_CaM12 = 9.0 / 1000.0;
        layer2[i][j][k].CaM21_to_CaMKIICaM21 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM21_to_CaM21 = 9.0 / 1000.0;
        layer2[i][j][k].CaM22_to_CaMKIICaM22 = 50.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM22_to_CaM22 = 2.25 / 1000.0;

	layer2[i][j][k].CaMKIICaM00_to_CaMKIICaM10 = 0.48 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM10_to_CaMKIICaM00 = 0.48 / 1000.0; 
        layer2[i][j][k].CaMKIICaM00_to_CaMKIICaM01 = 1.6 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM01_to_CaMKIICaM00 = 8.0 / 1000.0;
        layer2[i][j][k].CaMKIICaM10_to_CaMKIICaM11 = 1.6 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM11_to_CaMKIICaM10 = 8.0 / 1000.0;
        layer2[i][j][k].CaMKIICaM10_to_CaMKIICaM20 = 66.6 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM20_to_CaMKIICaM10 = 8.0 / 1000.0;
        layer2[i][j][k].CaMKIICaM01_to_CaMKIICaM11 = 0.48 / (6.02 * pow(10, 5) * pow(0.1, 3));
	layer2[i][j][k].CaMKIICaM11_to_CaMKIICaM01 = 0.48 / 1000.0;
        layer2[i][j][k].CaMKIICaM01_to_CaMKIICaM02 = 50.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM02_to_CaMKIICaM01 = 50.0 / 1000.0;
        layer2[i][j][k].CaMKIICaM11_to_CaMKIICaM21 = 66.6 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM21_to_CaMKIICaM11 = 8.0 / 1000.0;
        layer2[i][j][k].CaMKIICaM11_to_CaMKIICaM12 = 50.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM12_to_CaMKIICaM11 = 50.0 / 1000.0;
        layer2[i][j][k].CaMKIICaM20_to_CaMKIICaM21 = 1.6 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM21_to_CaMKIICaM20 = 8.0 / 1000.0;
        layer2[i][j][k].CaMKIICaM02_to_CaMKIICaM12 = 0.48 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM12_to_CaMKIICaM02 = 0.48 / 1000.0;
        layer2[i][j][k].CaMKIICaM21_to_CaMKIICaM22 = 50.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM22_to_CaMKIICaM21 = 50.0 / 1000.0;
        layer2[i][j][k].CaMKIICaM12_to_CaMKIICaM22 = 66.6 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaMKIICaM22_to_CaMKIICaM12 = 8.0 / 1000.0;

        layer2[i][j][k].CaM00_to_CaNCaM00 = 1.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM00_to_CaM00 = 1167.0 / 1000.0;
        layer2[i][j][k].CaM10_to_CaNCaM10 = 1.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM10_to_CaM10 = 11.670 / 1000.0;
        layer2[i][j][k].CaM01_to_CaNCaM01 = 1.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM01_to_CaM01 = 7.0 / 1000.0;
        layer2[i][j][k].CaM11_to_CaNCaM11 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM11_to_CaM11 = 0.7 / 1000.0;
        layer2[i][j][k].CaM20_to_CaNCaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM20_to_CaM20 = 5.83 / 1000.0;
        layer2[i][j][k].CaM02_to_CaNCaM02 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM02_to_CaM02 = 0.56 / 1000.0;
        layer2[i][j][k].CaM12_to_CaNCaM12 = 46.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM12_to_CaM12 = 0.02576 / 1000.0;
        layer2[i][j][k].CaM21_to_CaNCaM21 = 46.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM21_to_CaM21 = 0.161 / 1000.0;
        layer2[i][j][k].CaM22_to_CaNCaM22 = 46.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM22_to_CaM22 = 0.0013 / 1000.0;

        layer2[i][j][k].CaNCaM00_to_CaNCaM10 = 5.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM10_to_CaNCaM00 = 0.5 / 1000.0;
        layer2[i][j][k].CaNCaM00_to_CaNCaM01 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM01_to_CaNCaM00 = 6.0 / 1000.0;
        layer2[i][j][k].CaNCaM10_to_CaNCaM11 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM11_to_CaNCaM10 = 6.0 / 1000.0;
        layer2[i][j][k].CaNCaM10_to_CaNCaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM20_to_CaNCaM10 = 0.6 / 1000.0;
        layer2[i][j][k].CaNCaM01_to_CaNCaM11 = 5.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM11_to_CaNCaM01 = 0.5 / 1000.0;
        layer2[i][j][k].CaNCaM01_to_CaNCaM02 = 50.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM02_to_CaNCaM01 = 4.0 / 1000.0;
        layer2[i][j][k].CaNCaM11_to_CaNCaM21 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM21_to_CaNCaM11 = 0.6 / 1000.0;
        layer2[i][j][k].CaNCaM11_to_CaNCaM12 = 50.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM12_to_CaNCaM11 = 4.0 / 1000.0;
        layer2[i][j][k].CaNCaM20_to_CaNCaM21 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM21_to_CaNCaM20 = 6.0 / 1000.0;
        layer2[i][j][k].CaNCaM02_to_CaNCaM12 = 5.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM12_to_CaNCaM02 = 0.5 / 1000.0;
        layer2[i][j][k].CaNCaM21_to_CaNCaM22 = 50.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM22_to_CaNCaM21 = 4.0 / 1000.0;
        layer2[i][j][k].CaNCaM12_to_CaNCaM22 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].CaNCaM22_to_CaNCaM12 = 0.6 / 1000.0;

        layer2[i][j][k].CaM00_to_NgCaM00 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM00_to_CaM00 = 30.0 / 1000.0;
        layer2[i][j][k].CaM10_to_NgCaM10 = 3.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM10_to_CaM10 = 30.0 / 1000.0;
        layer2[i][j][k].CaM01_to_NgCaM01 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM01_to_CaM01 = 30.0 / 1000.0;
        layer2[i][j][k].CaM11_to_NgCaM11 = 3.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM11_to_CaM11 = 30.0 / 1000.0;
        layer2[i][j][k].CaM20_to_NgCaM20 = 0.1 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM20_to_CaM20 = 40.0 / 1000.0;
        layer2[i][j][k].CaM02_to_NgCaM02 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM02_to_CaM02 = 30.0 / 1000.0;
        layer2[i][j][k].CaM12_to_NgCaM12 = 3.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM12_to_CaM12 = 30.0 / 1000.0;
        layer2[i][j][k].CaM21_to_NgCaM21 = 0.1 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM21_to_CaM21 = 40.0 / 1000.0;
        layer2[i][j][k].CaM22_to_NgCaM22 = 0.1 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM22_to_CaM22 = 40.0 / 1000.0;

        layer2[i][j][k].NgCaM00_to_NgCaM10 = 1.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM10_to_NgCaM00 = 60.0 / 1000.0;
        layer2[i][j][k].NgCaM00_to_NgCaM01 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM01_to_NgCaM00 = 1000.0 / 1000.0;
        layer2[i][j][k].NgCaM10_to_NgCaM11 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM11_to_NgCaM10 = 1000.0 / 1000.0;
        layer2[i][j][k].NgCaM10_to_NgCaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM20_to_NgCaM10 = 480.0 / 1000.0;
        layer2[i][j][k].NgCaM01_to_NgCaM11 = 1.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM11_to_NgCaM01 = 60.0 / 1000.0;
        layer2[i][j][k].NgCaM01_to_NgCaM02 = 100.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM02_to_NgCaM01 = 1000.0 / 1000.0;
        layer2[i][j][k].NgCaM11_to_NgCaM21 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM21_to_NgCaM11 = 480.0 / 1000.0;
        layer2[i][j][k].NgCaM11_to_NgCaM12 = 100.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM12_to_NgCaM11 = 1000.0 / 1000.0;
        layer2[i][j][k].NgCaM20_to_NgCaM21 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM21_to_NgCaM20 = 1000.0 / 1000.0;
        layer2[i][j][k].NgCaM02_to_NgCaM12 = 1.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM12_to_NgCaM02 = 60.0 / 1000.0;
        layer2[i][j][k].NgCaM21_to_NgCaM22 = 100.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM22_to_NgCaM21 = 1000.0 / 1000.0;
        layer2[i][j][k].NgCaM12_to_NgCaM22 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].NgCaM22_to_NgCaM12 = 480.0 / 1000.0;

        layer2[i][j][k].CaMKIICaM00_to_Trapped00 = 10.0 / 1000.0;
        layer2[i][j][k].Trapped00_to_CaMKIICaM00 = 0.003 / 1000.0;
        layer2[i][j][k].CaMKIICaM10_to_Trapped10 = 10.0 / 1000.0;
        layer2[i][j][k].Trapped10_to_CaMKIICaM10 = 0.003 / 1000.0;
        layer2[i][j][k].CaMKIICaM01_to_Trapped01 = 10.0 / 1000.0;
        layer2[i][j][k].Trapped01_to_CaMKIICaM01 = 0.003 / 1000.0;
        layer2[i][j][k].CaMKIICaM11_to_Trapped11 = 10.0 / 1000.0;
        layer2[i][j][k].Trapped11_to_CaMKIICaM11 = 0.003 / 1000.0;
        layer2[i][j][k].CaMKIICaM20_to_Trapped20 = 10.0 / 1000.0;
        layer2[i][j][k].Trapped20_to_CaMKIICaM20 = 0.003 / 1000.0;
        layer2[i][j][k].CaMKIICaM02_to_Trapped02 = 10.0 / 1000.0;
        layer2[i][j][k].Trapped02_to_CaMKIICaM02 = 0.003 / 1000.0;
        layer2[i][j][k].CaMKIICaM12_to_Trapped12 = 10.0 / 1000.0;
        layer2[i][j][k].Trapped12_to_CaMKIICaM12 = 0.003 / 1000.0;
        layer2[i][j][k].CaMKIICaM21_to_Trapped21 = 10.0 / 1000.0;
        layer2[i][j][k].Trapped21_to_CaMKIICaM21 = 0.003 / 1000.0;
        layer2[i][j][k].CaMKIICaM22_to_Trapped22 = 10.0 / 1000.0;
        layer2[i][j][k].Trapped22_to_CaMKIICaM22 = 0.003 / 1000.0;
		  
        layer2[i][j][k].Trapped00_to_Auton = 0.09/ 1000.0;
        layer2[i][j][k].Auton_to_Trapped00 = 0.2 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].Trapped10_to_Auton = 0.09 / 1000.0;
        layer2[i][j][k].Auton_to_Trapped10 = 2.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].Trapped01_to_Auton = 0.09 / 1000.0;
        layer2[i][j][k].Auton_to_Trapped01 = 2.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].Trapped11_to_Auton = 0.045 / 1000.0;
        layer2[i][j][k].Auton_to_Trapped11 = 10 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].Trapped20_to_Auton = 0.045 / 1000.0;
        layer2[i][j][k].Auton_to_Trapped20 = 10 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].Trapped02_to_Auton = 0.045 / 1000.0;
        layer2[i][j][k].Auton_to_Trapped02 = 10 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].Trapped12_to_Auton = 0.009 / 1000.0;
        layer2[i][j][k].Auton_to_Trapped12 = 20 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].Trapped21_to_Auton = 0.009 / 1000.0;
        layer2[i][j][k].Auton_to_Trapped21 = 20 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].Trapped22_to_Auton = 0.00225 / 1000.0;
        layer2[i][j][k].Auton_to_Trapped22 = 50 / (6.02 * pow(10, 5) * pow(0.1, 3));

        layer2[i][j][k].Trapped00_to_Trapped10 = 0.4 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].Trapped10_to_Trapped00 = 0.4 / 1000.0;
        layer2[i][j][k].Trapped00_to_Trapped01 = 6.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].Trapped01_to_Trapped00 = 1.2 / 1000.0;
        layer2[i][j][k].Trapped10_to_Trapped11 = 6.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].Trapped11_to_Trapped10 = 1.2 / 1000.0;
        layer2[i][j][k].Trapped10_to_Trapped20 = 3.3 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].Trapped20_to_Trapped10 = 0.4 / 1000.0;
        layer2[i][j][k].Trapped01_to_Trapped11 = 0.4 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].Trapped11_to_Trapped01 = 0.4 / 1000.0;
        layer2[i][j][k].Trapped01_to_Trapped02 = 11.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].Trapped02_to_Trapped01 = 11.0 / 1000.0;
        layer2[i][j][k].Trapped11_to_Trapped21 = 3.3 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].Trapped21_to_Trapped11 = 0.4 / 1000.0;
        layer2[i][j][k].Trapped11_to_Trapped12 = 11.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].Trapped12_to_Trapped11 = 11.0 / 1000.0;
        layer2[i][j][k].Trapped20_to_Trapped21 = 6.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].Trapped21_to_Trapped20 = 1.2 / 1000.0;
        layer2[i][j][k].Trapped02_to_Trapped12 = 0.4 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].Trapped12_to_Trapped02 = 0.4 / 1000.0;
        layer2[i][j][k].Trapped21_to_Trapped22 = 11.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer2[i][j][k].Trapped22_to_Trapped21 = 11.0 / 1000.0;
        layer2[i][j][k].Trapped12_to_Trapped22 = 3.3 / (6.02 * pow(10, 5) * pow(0.1, 3)); 
        layer2[i][j][k].Trapped22_to_Trapped12 = 0.4 / 1000.0;

        layer2[i][j][k].Auton_to_CaMKII = 0.003 / 1000.0;
        layer2[i][j][k].Auton_to_Capped = 0.1 / 1000.0;
        layer2[i][j][k].Capped_to_Auton = 0.01 / 1000.0;        
	
	  layer2[i][j][k].CB00_to_CB10 = 5.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer2[i][j][k].CB10_to_CB00 = 2.6 / 1000;
	  layer2[i][j][k].CB00_to_CB01 = 43.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer2[i][j][k].CB01_to_CB00 = 35.8 / 1000;
	  layer2[i][j][k].CB10_to_CB11 = 43.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
          layer2[i][j][k].CB11_to_CB10 = 35.8 / 1000;
	  layer2[i][j][k].CB10_to_CB20 = 5.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer2[i][j][k].CB20_to_CB10 = 2.6 / 1000;
	  layer2[i][j][k].CB01_to_CB11 = 5.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer2[i][j][k].CB11_to_CB01 = 2.6 / 1000;
	  layer2[i][j][k].CB01_to_CB02 = 43.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer2[i][j][k].CB02_to_CB01 = 35.8 / 1000;
	  layer2[i][j][k].CB11_to_CB21 = 5.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer2[i][j][k].CB21_to_CB11 = 2.6 / 1000;
	  layer2[i][j][k].CB11_to_CB12 = 43.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer2[i][j][k].CB12_to_CB11 = 35.8 / 1000;
	  layer2[i][j][k].CB20_to_CB21 = 43.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer2[i][j][k].CB21_to_CB20 = 35.8 / 1000;
	  layer2[i][j][k].CB02_to_CB12 = 5.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer2[i][j][k].CB12_to_CB02 = 2.6 / 1000;
	  layer2[i][j][k].CB12_to_CB22 = 5.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer2[i][j][k].CB22_to_CB12 = 2.6 / 1000;
	  layer2[i][j][k].CB21_to_CB22 = 43.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer2[i][j][k].CB22_to_CB21 = 35.8 / 1000;
      }
  
  for(i=0; i<8; i++){
	        layer3[i].CaM00_to_CaM10 = 1.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer3[i].CaM10_to_CaM00 = 10.0 / 1000.0;
		layer3[i].CaM00_to_CaM01 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer3[i].CaM01_to_CaM00 = 1000.0 / 1000.0;
		layer3[i].CaM10_to_CaM11 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer3[i].CaM11_to_CaM10 = 1000.0 / 1000.0;
		layer3[i].CaM10_to_CaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer3[i].CaM20_to_CaM10 = 12.0 / 1000.0;
		layer3[i].CaM01_to_CaM11 = 1.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer3[i].CaM11_to_CaM01 = 10.0 / 1000.0;
		layer3[i].CaM01_to_CaM02 = 100.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer3[i].CaM02_to_CaM01 = 1000.0 / 1000.0;
		layer3[i].CaM11_to_CaM21 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer3[i].CaM21_to_CaM11 = 12.0 / 1000.0;
		layer3[i].CaM11_to_CaM12 = 100.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer3[i].CaM12_to_CaM11 = 1000.0 / 1000.0;
		layer3[i].CaM20_to_CaM21 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer3[i].CaM21_to_CaM20 = 1000.0 / 1000.0;
		layer3[i].CaM02_to_CaM12 = 1.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer3[i].CaM12_to_CaM02 = 10.0 / 1000.0;
		layer3[i].CaM21_to_CaM22 = 100.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer3[i].CaM22_to_CaM21 = 1000.0 / 1000.0;
		layer3[i].CaM12_to_CaM22 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
		layer3[i].CaM22_to_CaM12 = 12.0 / 1000.0;

	layer3[i].CaM00_to_CaMKIICaM00 = 0.2 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM00_to_CaM00 = 90.0 / 1000.0;
        layer3[i].CaM10_to_CaMKIICaM10 = 2.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM10_to_CaM10 = 90.0 / 1000.0;
        layer3[i].CaM01_to_CaMKIICaM01 = 2.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM01_to_CaM01 = 90.0 / 1000.0;
        layer3[i].CaM11_to_CaMKIICaM11 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM11_to_CaM11 = 45.0 / 1000.0;
        layer3[i].CaM20_to_CaMKIICaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM20_to_CaM20 = 45.0 / 1000.0;
        layer3[i].CaM02_to_CaMKIICaM02 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM02_to_CaM02 = 45.0 / 1000.0;
        layer3[i].CaM12_to_CaMKIICaM12 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM12_to_CaM12 = 9.0 / 1000.0;
        layer3[i].CaM21_to_CaMKIICaM21 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM21_to_CaM21 = 9.0 / 1000.0;
        layer3[i].CaM22_to_CaMKIICaM22 = 50.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM22_to_CaM22 = 2.25 / 1000.0;

	layer3[i].CaMKIICaM00_to_CaMKIICaM10 = 0.48 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM10_to_CaMKIICaM00 = 0.48 / 1000.0; 
        layer3[i].CaMKIICaM00_to_CaMKIICaM01 = 1.6 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM01_to_CaMKIICaM00 = 8.0 / 1000.0;
        layer3[i].CaMKIICaM10_to_CaMKIICaM11 = 1.6 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM11_to_CaMKIICaM10 = 8.0 / 1000.0;
        layer3[i].CaMKIICaM10_to_CaMKIICaM20 = 66.6 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM20_to_CaMKIICaM10 = 8.0 / 1000.0;
        layer3[i].CaMKIICaM01_to_CaMKIICaM11 = 0.48 / (6.02 * pow(10, 5) * pow(0.1, 3));
	layer3[i].CaMKIICaM11_to_CaMKIICaM01 = 0.48 / 1000.0;
        layer3[i].CaMKIICaM01_to_CaMKIICaM02 = 50.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM02_to_CaMKIICaM01 = 50.0 / 1000.0;
        layer3[i].CaMKIICaM11_to_CaMKIICaM21 = 66.6 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM21_to_CaMKIICaM11 = 8.0 / 1000.0;
        layer3[i].CaMKIICaM11_to_CaMKIICaM12 = 50.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM12_to_CaMKIICaM11 = 50.0 / 1000.0;
        layer3[i].CaMKIICaM20_to_CaMKIICaM21 = 1.6 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM21_to_CaMKIICaM20 = 8.0 / 1000.0;
        layer3[i].CaMKIICaM02_to_CaMKIICaM12 = 0.48 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM12_to_CaMKIICaM02 = 0.48 / 1000.0;
        layer3[i].CaMKIICaM21_to_CaMKIICaM22 = 50.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM22_to_CaMKIICaM21 = 50.0 / 1000.0;
        layer3[i].CaMKIICaM12_to_CaMKIICaM22 = 66.6 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaMKIICaM22_to_CaMKIICaM12 = 8.0 / 1000.0;

        layer3[i].CaM00_to_CaNCaM00 = 1.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM00_to_CaM00 = 1167.0 / 1000.0;
        layer3[i].CaM10_to_CaNCaM10 = 1.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM10_to_CaM10 = 11.67 / 1000.0;
        layer3[i].CaM01_to_CaNCaM01 = 1.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM01_to_CaM01 = 7.0 / 1000.0;
        layer3[i].CaM11_to_CaNCaM11 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM11_to_CaM11 = 0.7 / 1000.0;
        layer3[i].CaM20_to_CaNCaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM20_to_CaM20 = 5.83 / 1000.0;
        layer3[i].CaM02_to_CaNCaM02 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM02_to_CaM02 = 0.56 / 1000.0;
        layer3[i].CaM12_to_CaNCaM12 = 46.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM12_to_CaM12 = 0.02576 / 1000.0;
        layer3[i].CaM21_to_CaNCaM21 = 46.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM21_to_CaM21 = 0.161 / 1000.0;
        layer3[i].CaM22_to_CaNCaM22 = 46.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM22_to_CaM22 = 0.0013 / 1000.0;

        layer3[i].CaNCaM00_to_CaNCaM10 = 5.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM10_to_CaNCaM00 = 0.5 / 1000.0;
        layer3[i].CaNCaM00_to_CaNCaM01 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM01_to_CaNCaM00 = 6.0 / 1000.0;
        layer3[i].CaNCaM10_to_CaNCaM11 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM11_to_CaNCaM10 = 6.0 / 1000.0;
        layer3[i].CaNCaM10_to_CaNCaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM20_to_CaNCaM10 = 0.6 / 1000.0;
        layer3[i].CaNCaM01_to_CaNCaM11 = 5.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM11_to_CaNCaM01 = 0.5 / 1000.0;
        layer3[i].CaNCaM01_to_CaNCaM02 = 50.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM02_to_CaNCaM01 = 4.0 / 1000.0;
        layer3[i].CaNCaM11_to_CaNCaM21 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM21_to_CaNCaM11 = 0.6 / 1000.0;
        layer3[i].CaNCaM11_to_CaNCaM12 = 50.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM12_to_CaNCaM11 = 4.0 / 1000.0;
        layer3[i].CaNCaM20_to_CaNCaM21 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM21_to_CaNCaM20 = 6.0 / 1000.0;
        layer3[i].CaNCaM02_to_CaNCaM12 = 5.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM12_to_CaNCaM02 = 0.5 / 1000.0;
        layer3[i].CaNCaM21_to_CaNCaM22 = 50.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM22_to_CaNCaM21 = 4.0 / 1000.0;
        layer3[i].CaNCaM12_to_CaNCaM22 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].CaNCaM22_to_CaNCaM12 = 0.6 / 1000.0;

        layer3[i].CaM00_to_NgCaM00 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM00_to_CaM00 = 30.0 / 1000.0;
        layer3[i].CaM10_to_NgCaM10 = 3.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM10_to_CaM10 = 30.0 / 1000.0;
        layer3[i].CaM01_to_NgCaM01 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM01_to_CaM01 = 30.0 / 1000.0;
        layer3[i].CaM11_to_NgCaM11 = 3.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM11_to_CaM11 = 30.0 / 1000.0;
        layer3[i].CaM20_to_NgCaM20 = 0.1 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM20_to_CaM20 = 40.0 / 1000.0;
        layer3[i].CaM02_to_NgCaM02 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM02_to_CaM02 = 30.0 / 1000.0;
        layer3[i].CaM12_to_NgCaM12 = 3.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM12_to_CaM12 = 30.0 / 1000.0;
        layer3[i].CaM21_to_NgCaM21 = 0.1 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM21_to_CaM21 = 40.0 / 1000.0;
        layer3[i].CaM22_to_NgCaM22 = 0.1 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM22_to_CaM22 = 40.0 / 1000.0;

        layer3[i].NgCaM00_to_NgCaM10 = 1.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM10_to_NgCaM00 = 60.0 / 1000.0;
        layer3[i].NgCaM00_to_NgCaM01 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM01_to_NgCaM00 = 1000.0 / 1000.0;
        layer3[i].NgCaM10_to_NgCaM11 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM11_to_NgCaM10 = 1000.0 / 1000.0;
        layer3[i].NgCaM10_to_NgCaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM20_to_NgCaM10 = 480.0 / 1000.0;
        layer3[i].NgCaM01_to_NgCaM11 = 1.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM11_to_NgCaM01 = 60.0 / 1000.0;
        layer3[i].NgCaM01_to_NgCaM02 = 100.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM02_to_NgCaM01 = 1000.0 / 1000.0;
        layer3[i].NgCaM11_to_NgCaM21 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM21_to_NgCaM11 = 480.0 / 1000.0;
        layer3[i].NgCaM11_to_NgCaM12 = 100.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM12_to_NgCaM11 = 1000.0 / 1000.0;
        layer3[i].NgCaM20_to_NgCaM21 = 20.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM21_to_NgCaM20 = 1000.0 / 1000.0;
        layer3[i].NgCaM02_to_NgCaM12 = 1.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM12_to_NgCaM02 = 60.0 / 1000.0;
        layer3[i].NgCaM21_to_NgCaM22 = 100.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM22_to_NgCaM21 = 1000.0 / 1000.0;
        layer3[i].NgCaM12_to_NgCaM22 = 10.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].NgCaM22_to_NgCaM12 = 480.0 / 1000.0;

        layer3[i].CaMKIICaM00_to_Trapped00 = 10.0 / 1000.0;
        layer3[i].Trapped00_to_CaMKIICaM00 = 0.003 / 1000.0;
        layer3[i].CaMKIICaM10_to_Trapped10 = 10.0 / 1000.0;
        layer3[i].Trapped10_to_CaMKIICaM10 = 0.003 / 1000.0;
        layer3[i].CaMKIICaM01_to_Trapped01 = 10.0 / 1000.0;
        layer3[i].Trapped01_to_CaMKIICaM01 = 0.003 / 1000.0;
        layer3[i].CaMKIICaM11_to_Trapped11 = 10.0 / 1000.0;
        layer3[i].Trapped11_to_CaMKIICaM11 = 0.003 / 1000.0;
        layer3[i].CaMKIICaM20_to_Trapped20 = 10.0 / 1000.0;
        layer3[i].Trapped20_to_CaMKIICaM20 = 0.003 / 1000.0;
        layer3[i].CaMKIICaM02_to_Trapped02 = 10.0 / 1000.0;
        layer3[i].Trapped02_to_CaMKIICaM02 = 0.003 / 1000.0;
        layer3[i].CaMKIICaM12_to_Trapped12 = 10.0 / 1000.0;
        layer3[i].Trapped12_to_CaMKIICaM12 = 0.003 / 1000.0;
        layer3[i].CaMKIICaM21_to_Trapped21 = 10.0 / 1000.0;
        layer3[i].Trapped21_to_CaMKIICaM21 = 0.003 / 1000.0;
        layer3[i].CaMKIICaM22_to_Trapped22 = 10.0 / 1000.0;
        layer3[i].Trapped22_to_CaMKIICaM22 = 0.003 / 1000.0;

        layer3[i].Trapped00_to_Auton = 0.09 / 1000.0;
        layer3[i].Auton_to_Trapped00 = 0.2 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped10_to_Auton = 0.09 / 1000.0;
        layer3[i].Auton_to_Trapped10 = 2.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped01_to_Auton = 0.09 / 1000.0;
        layer3[i].Auton_to_Trapped01 = 2.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped11_to_Auton = 0.045 / 1000.0;
        layer3[i].Auton_to_Trapped11 = 10 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped20_to_Auton = 0.045 / 1000.0;
        layer3[i].Auton_to_Trapped20 = 10 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped02_to_Auton = 0.045 / 1000.0;
        layer3[i].Auton_to_Trapped02 = 10 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped12_to_Auton = 0.009 / 1000.0;
        layer3[i].Auton_to_Trapped12 = 20 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped21_to_Auton = 0.009 / 1000.0;
        layer3[i].Auton_to_Trapped21 = 20 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped22_to_Auton = 0.00225 / 1000.0;
        layer3[i].Auton_to_Trapped22 = 50 / (6.02 * pow(10, 5) * pow(0.1, 3));

        layer3[i].Trapped00_to_Trapped10 = 0.4 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped10_to_Trapped00 = 0.4 / 1000.0;
        layer3[i].Trapped00_to_Trapped01 = 6.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped01_to_Trapped00 = 1.2 / 1000.0;
        layer3[i].Trapped10_to_Trapped11 = 6.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped11_to_Trapped10 = 1.2 / 1000.0;
        layer3[i].Trapped10_to_Trapped20 = 3.3 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped20_to_Trapped10 = 0.4 / 1000.0;
        layer3[i].Trapped01_to_Trapped11 = 0.4 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped11_to_Trapped01 = 0.4 / 1000.0;
        layer3[i].Trapped01_to_Trapped02 = 11.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped02_to_Trapped01 = 11.0 / 1000.0;
        layer3[i].Trapped11_to_Trapped21 = 3.3 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped21_to_Trapped11 = 0.4 / 1000.0;
        layer3[i].Trapped11_to_Trapped12 = 11.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped12_to_Trapped11 = 11.0 / 1000.0;
        layer3[i].Trapped20_to_Trapped21 = 6.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped21_to_Trapped20 = 1.2 / 1000.0;
        layer3[i].Trapped02_to_Trapped12 = 0.4 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped12_to_Trapped02 = 0.4 / 1000.0;
        layer3[i].Trapped21_to_Trapped22 = 11.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped22_to_Trapped21 = 11.0 / 1000.0;
        layer3[i].Trapped12_to_Trapped22 = 3.3 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer3[i].Trapped22_to_Trapped12 = 0.4 / 1000.0;

        layer3[i].Auton_to_CaMKII = 0.003 / 1000.0;
        layer3[i].Auton_to_Capped = 0.1 / 1000.0;
        layer3[i].Capped_to_Auton = 0.01 / 1000.0;        

	  layer3[i].CB00_to_CB10 = 5.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer3[i].CB10_to_CB00 = 2.6 / 1000;
	  layer3[i].CB00_to_CB01 = 43.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer3[i].CB01_to_CB00 = 35.8 / 1000;
	  layer3[i].CB10_to_CB11 = 43.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
          layer3[i].CB11_to_CB10 = 35.8 / 1000;
	  layer3[i].CB10_to_CB20 = 5.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer3[i].CB20_to_CB10 = 2.6 / 1000;
	  layer3[i].CB01_to_CB11 = 5.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer3[i].CB11_to_CB01 = 2.6 / 1000;
	  layer3[i].CB01_to_CB02 = 43.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer3[i].CB02_to_CB01 = 35.8 / 1000;
	  layer3[i].CB11_to_CB21 = 5.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer3[i].CB21_to_CB11 = 2.6 / 1000;
	  layer3[i].CB11_to_CB12 = 43.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer3[i].CB12_to_CB11 = 35.8 / 1000;
	  layer3[i].CB20_to_CB21 = 43.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer3[i].CB21_to_CB20 = 35.8 / 1000;
	  layer3[i].CB02_to_CB12 = 5.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer3[i].CB12_to_CB02 = 2.6 / 1000;
	  layer3[i].CB12_to_CB22 = 5.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer3[i].CB22_to_CB12 = 2.6 / 1000;
	  layer3[i].CB21_to_CB22 = 43.5 / (6.02 * pow(10, 5) * pow(0.1, 3));
	  layer3[i].CB22_to_CB21 = 35.8 / 1000;
  }
  
  for(i=0; i<4; i++)
    for(j=0; j<2; j++)
      for(k=0; k<2; k++){
                layer4[i][j][k].CaM00_to_CaM10 = 1.0 / (6.02 * pow(10, 5) * pow(0.05, 3));
		layer4[i][j][k].CaM10_to_CaM00 = 10.0 / 1000.0;
		layer4[i][j][k].CaM00_to_CaM01 = 20.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
		layer4[i][j][k].CaM01_to_CaM00 = 1000.0 / 1000.0;
		layer4[i][j][k].CaM10_to_CaM11 = 20.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
		layer4[i][j][k].CaM11_to_CaM10 = 1000.0 / 1000.0;
		layer4[i][j][k].CaM10_to_CaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
		layer4[i][j][k].CaM20_to_CaM10 = 12.0 / 1000.0;
		layer4[i][j][k].CaM01_to_CaM11 = 1.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
		layer4[i][j][k].CaM11_to_CaM01 = 10.0 / 1000.0;
		layer4[i][j][k].CaM01_to_CaM02 = 100.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
		layer4[i][j][k].CaM02_to_CaM01 = 1000.0 / 1000.0;
		layer4[i][j][k].CaM11_to_CaM21 = 10.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
		layer4[i][j][k].CaM21_to_CaM11 = 12.0 / 1000.0;
		layer4[i][j][k].CaM11_to_CaM12 = 100.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
		layer4[i][j][k].CaM12_to_CaM11 = 1000.0 / 1000.0;
		layer4[i][j][k].CaM20_to_CaM21 = 20.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
		layer4[i][j][k].CaM21_to_CaM20 = 1000.0 / 1000.0;
		layer4[i][j][k].CaM02_to_CaM12 = 1.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
		layer4[i][j][k].CaM12_to_CaM02 = 10.0 / 1000.0;
		layer4[i][j][k].CaM21_to_CaM22 = 100.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
		layer4[i][j][k].CaM22_to_CaM21 = 1000.0 / 1000.0;
		layer4[i][j][k].CaM12_to_CaM22 = 10.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
		layer4[i][j][k].CaM22_to_CaM12 = 12.0 / 1000.0;

	layer4[i][j][k].CaM00_to_CaMKIICaM00 = 0.2 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM00_to_CaM00 = 90.0 / 1000.0;
        layer4[i][j][k].CaM10_to_CaMKIICaM10 = 2.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM10_to_CaM10 = 90.0 / 1000.0;
        layer4[i][j][k].CaM01_to_CaMKIICaM01 = 2.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM01_to_CaM01 = 90.0 / 1000.0;
        layer4[i][j][k].CaM11_to_CaMKIICaM11 = 10.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM11_to_CaM11 = 45.0 / 1000.0;
        layer4[i][j][k].CaM20_to_CaMKIICaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM20_to_CaM20 = 45.0 / 1000.0;
        layer4[i][j][k].CaM02_to_CaMKIICaM02 = 10.0/ (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM02_to_CaM02 = 45.0 / 1000.0;
        layer4[i][j][k].CaM12_to_CaMKIICaM12 = 20.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM12_to_CaM12 = 9.0 / 1000.0;
        layer4[i][j][k].CaM21_to_CaMKIICaM21 = 20.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM21_to_CaM21 = 9.0 / 1000.0;
        layer4[i][j][k].CaM22_to_CaMKIICaM22 = 50.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM22_to_CaM22 = 2.25 / 1000.0;

	layer4[i][j][k].CaMKIICaM00_to_CaMKIICaM10 = 0.48 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM10_to_CaMKIICaM00 = 0.48 / 1000.0; 
        layer4[i][j][k].CaMKIICaM00_to_CaMKIICaM01 = 1.6 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM01_to_CaMKIICaM00 = 8.0 / 1000.0;
        layer4[i][j][k].CaMKIICaM10_to_CaMKIICaM11 = 1.6 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM11_to_CaMKIICaM10 = 8.0 / 1000.0;
        layer4[i][j][k].CaMKIICaM10_to_CaMKIICaM20 = 66.6 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM20_to_CaMKIICaM10 = 8.0 / 1000.0;
        layer4[i][j][k].CaMKIICaM01_to_CaMKIICaM11 = 0.48 / (6.02 * pow(10, 5) * pow(0.5, 3));
	layer4[i][j][k].CaMKIICaM11_to_CaMKIICaM01 = 0.48 / 1000.0;
        layer4[i][j][k].CaMKIICaM01_to_CaMKIICaM02 = 50.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM02_to_CaMKIICaM01 = 50.0 / 1000.0;
        layer4[i][j][k].CaMKIICaM11_to_CaMKIICaM21 = 66.6 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM21_to_CaMKIICaM11 = 8.0 / 1000.0;
        layer4[i][j][k].CaMKIICaM11_to_CaMKIICaM12 = 50.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM12_to_CaMKIICaM11 = 50.0 / 1000.0;
        layer4[i][j][k].CaMKIICaM20_to_CaMKIICaM21 = 1.6 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM21_to_CaMKIICaM20 = 8.0 / 1000.0;
        layer4[i][j][k].CaMKIICaM02_to_CaMKIICaM12 = 0.48 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM12_to_CaMKIICaM02 = 0.48 / 1000.0;
        layer4[i][j][k].CaMKIICaM21_to_CaMKIICaM22 = 50.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM22_to_CaMKIICaM21 = 50.0 / 1000.0;
        layer4[i][j][k].CaMKIICaM12_to_CaMKIICaM22 = 66.6 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaMKIICaM22_to_CaMKIICaM12 = 8.0 / 1000.0;

        layer4[i][j][k].CaM00_to_CaNCaM00 = 1.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM00_to_CaM00 = 1167.0 / 1000.0;
        layer4[i][j][k].CaM10_to_CaNCaM10 = 1.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM10_to_CaM10 = 11.67 / 1000.0;
        layer4[i][j][k].CaM01_to_CaNCaM01 = 1.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM01_to_CaM01 = 7.0 / 1000.0;
        layer4[i][j][k].CaM11_to_CaNCaM11 = 10.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM11_to_CaM11 = 0.7 / 1000.0;
        layer4[i][j][k].CaM20_to_CaNCaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM20_to_CaM20 = 5.83 / 1000.0;
        layer4[i][j][k].CaM02_to_CaNCaM02 = 10.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM02_to_CaM02 = 0.56 / 1000.0;
        layer4[i][j][k].CaM12_to_CaNCaM12 = 46.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM12_to_CaM12 = 0.02576/ 1000.0;
        layer4[i][j][k].CaM21_to_CaNCaM21 = 46.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM21_to_CaM21 = 0.161 / 1000.0;
        layer4[i][j][k].CaM22_to_CaNCaM22 = 46.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM22_to_CaM22 = 0.0013 / 1000.0;

        layer4[i][j][k].CaNCaM00_to_CaNCaM10 = 5.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM10_to_CaNCaM00 = 0.5 / 1000.0;
        layer4[i][j][k].CaNCaM00_to_CaNCaM01 = 20.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM01_to_CaNCaM00 = 6.0 / 1000.0;
        layer4[i][j][k].CaNCaM10_to_CaNCaM11 = 20.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM11_to_CaNCaM10 = 6.0 / 1000.0;
        layer4[i][j][k].CaNCaM10_to_CaNCaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM20_to_CaNCaM10 = 0.6 / 1000.0;
        layer4[i][j][k].CaNCaM01_to_CaNCaM11 = 5.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM11_to_CaNCaM01 = 0.5 / 1000.0;
        layer4[i][j][k].CaNCaM01_to_CaNCaM02 = 50.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM02_to_CaNCaM01 = 4.0 / 1000.0;
        layer4[i][j][k].CaNCaM11_to_CaNCaM21 = 10.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM21_to_CaNCaM11 = 0.6 / 1000.0;
        layer4[i][j][k].CaNCaM11_to_CaNCaM12 = 50.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM12_to_CaNCaM11 = 4.0 / 1000.0;
        layer4[i][j][k].CaNCaM20_to_CaNCaM21 = 20.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM21_to_CaNCaM20 = 6.0 / 1000.0;
        layer4[i][j][k].CaNCaM02_to_CaNCaM12 = 5.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM12_to_CaNCaM02 = 0.5 / 1000.0;
        layer4[i][j][k].CaNCaM21_to_CaNCaM22 = 50.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM22_to_CaNCaM21 = 4.0 / 1000.0;
        layer4[i][j][k].CaNCaM12_to_CaNCaM22 = 10.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].CaNCaM22_to_CaNCaM12 = 0.6 / 1000.0;

        layer4[i][j][k].CaM00_to_NgCaM00 = 20.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM00_to_CaM00 = 30.0 / 1000.0;
        layer4[i][j][k].CaM10_to_NgCaM10 = 3.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM10_to_CaM10 = 30.0 / 1000.0;
        layer4[i][j][k].CaM01_to_NgCaM01 = 20.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM01_to_CaM01 = 30.0 / 1000.0;
        layer4[i][j][k].CaM11_to_NgCaM11 = 3.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM11_to_CaM11 = 30.0 / 1000.0;
        layer4[i][j][k].CaM20_to_NgCaM20 = 0.1 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM20_to_CaM20 = 40.0 / 1000.0;
        layer4[i][j][k].CaM02_to_NgCaM02 = 20.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM02_to_CaM02 = 30.0 / 1000.0;
        layer4[i][j][k].CaM12_to_NgCaM12 = 3.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM12_to_CaM12 = 30.0/ 1000.0;
        layer4[i][j][k].CaM21_to_NgCaM21 = 0.1 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM21_to_CaM21 = 40.0 / 1000.0;
        layer4[i][j][k].CaM22_to_NgCaM22 = 0.1 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM22_to_CaM22 = 40.0 / 1000.0;

        layer4[i][j][k].NgCaM00_to_NgCaM10 = 1.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM10_to_NgCaM00 = 60.0 / 1000.0;
        layer4[i][j][k].NgCaM00_to_NgCaM01 = 20.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM01_to_NgCaM00 = 1000.0 / 1000.0;
        layer4[i][j][k].NgCaM10_to_NgCaM11 = 20.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM11_to_NgCaM10 = 1000.0 / 1000.0;
        layer4[i][j][k].NgCaM10_to_NgCaM20 = 10.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM20_to_NgCaM10 = 480.0 / 1000.0;
        layer4[i][j][k].NgCaM01_to_NgCaM11 = 1.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM11_to_NgCaM01 = 60.0 / 1000.0;
        layer4[i][j][k].NgCaM01_to_NgCaM02 = 100.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM02_to_NgCaM01 = 1000.0 / 1000.0;
        layer4[i][j][k].NgCaM11_to_NgCaM21 = 10.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM21_to_NgCaM11 = 480.0 / 1000.0;
        layer4[i][j][k].NgCaM11_to_NgCaM12 = 100.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM12_to_NgCaM11 = 1000.0 / 1000.0;
        layer4[i][j][k].NgCaM20_to_NgCaM21 = 20.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM21_to_NgCaM20 = 1000.0 / 1000.0;
        layer4[i][j][k].NgCaM02_to_NgCaM12 = 1.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM12_to_NgCaM02 = 60.0 / 1000.0;
        layer4[i][j][k].NgCaM21_to_NgCaM22 = 100.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM22_to_NgCaM21 = 1000.0 / 1000.0;
        layer4[i][j][k].NgCaM12_to_NgCaM22 = 10.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].NgCaM22_to_NgCaM12 = 480.0 / 1000.0;

        layer4[i][j][k].CaMKIICaM00_to_Trapped00 = 10.0 / 1000.0;
        layer4[i][j][k].Trapped00_to_CaMKIICaM00 = 0.003 / 1000.0;
        layer4[i][j][k].CaMKIICaM10_to_Trapped10 = 10.0 / 1000.0;
        layer4[i][j][k].Trapped10_to_CaMKIICaM10 = 0.003 / 1000.0;
        layer4[i][j][k].CaMKIICaM01_to_Trapped01 = 10.0 / 1000.0;
        layer4[i][j][k].Trapped01_to_CaMKIICaM01 = 0.003 / 1000.0;
        layer4[i][j][k]. CaMKIICaM11_to_Trapped11 = 10.0 / 1000.0;
        layer4[i][j][k].Trapped11_to_CaMKIICaM11 = 0.003 / 1000.0;
        layer4[i][j][k].CaMKIICaM20_to_Trapped20 = 10.0 / 1000.0;
        layer4[i][j][k].Trapped20_to_CaMKIICaM20 = 0.003 / 1000.0;
        layer4[i][j][k].CaMKIICaM02_to_Trapped02 = 10.0 / 1000.0;
        layer4[i][j][k].Trapped02_to_CaMKIICaM02 = 0.003 / 1000.0;
        layer4[i][j][k].CaMKIICaM12_to_Trapped12 = 10.0 / 1000.0;
        layer4[i][j][k].Trapped12_to_CaMKIICaM12 = 0.003 / 1000.0;
        layer4[i][j][k].CaMKIICaM21_to_Trapped21 = 10.0 / 1000.0;
        layer4[i][j][k].Trapped21_to_CaMKIICaM21 = 0.003 / 1000.0;
        layer4[i][j][k].CaMKIICaM22_to_Trapped22 = 10.0 / 1000.0;
        layer4[i][j][k].Trapped22_to_CaMKIICaM22 = 0.003 / 1000.0;
		
        layer4[i][j][k].Trapped00_to_Auton = 0.09/ 1000.0;
        layer4[i][j][k].Auton_to_Trapped00 = 0.2 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer4[i][j][k].Trapped10_to_Auton = 0.09 / 1000.0;
        layer4[i][j][k].Auton_to_Trapped10 = 2.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer4[i][j][k].Trapped01_to_Auton = 0.09 / 1000.0;
        layer4[i][j][k].Auton_to_Trapped01 = 2.0 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer4[i][j][k].Trapped11_to_Auton = 0.045 / 1000.0;
        layer4[i][j][k].Auton_to_Trapped11 = 10 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer4[i][j][k].Trapped20_to_Auton = 0.045 / 1000.0;
        layer4[i][j][k].Auton_to_Trapped20 = 10 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer4[i][j][k].Trapped02_to_Auton = 0.045 / 1000.0;
        layer4[i][j][k].Auton_to_Trapped02 = 10 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer4[i][j][k].Trapped12_to_Auton = 0.009 / 1000.0;
        layer4[i][j][k].Auton_to_Trapped12 = 20 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer4[i][j][k].Trapped21_to_Auton = 0.009 / 1000.0;
        layer4[i][j][k].Auton_to_Trapped21 = 20 / (6.02 * pow(10, 5) * pow(0.1, 3));
        layer4[i][j][k].Trapped22_to_Auton = 0.00225 / 1000.0;
        layer4[i][j][k].Auton_to_Trapped22 = 50.0/ (6.02 * pow(10, 5) * pow(0.1, 3));
  
        layer4[i][j][k].Trapped00_to_Trapped10 = 0.4 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].Trapped10_to_Trapped00 = 0.4 / 1000.0;
        layer4[i][j][k].Trapped00_to_Trapped01 = 6.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].Trapped01_to_Trapped00 = 1.2 / 1000.0;
        layer4[i][j][k].Trapped10_to_Trapped11 = 6.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].Trapped11_to_Trapped10 = 1.2 / 1000.0;
        layer4[i][j][k].Trapped10_to_Trapped20 = 3.3 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].Trapped20_to_Trapped10 = 0.4 / 1000.0;
        layer4[i][j][k].Trapped01_to_Trapped11 = 0.4 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].Trapped11_to_Trapped01 = 0.4 / 1000.0;
        layer4[i][j][k].Trapped01_to_Trapped02 = 11.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].Trapped02_to_Trapped01 = 11.0 / 1000.0;
        layer4[i][j][k].Trapped11_to_Trapped21 = 3.3 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].Trapped21_to_Trapped11 = 0.4 / 1000.0;
        layer4[i][j][k].Trapped11_to_Trapped12 = 11.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].Trapped12_to_Trapped11 = 11.0 / 1000.0;
        layer4[i][j][k].Trapped20_to_Trapped21 = 6.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].Trapped21_to_Trapped20 = 1.2 / 1000.0;
        layer4[i][j][k].Trapped02_to_Trapped12 = 0.4 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].Trapped12_to_Trapped02 = 0.4 / 1000.0;
        layer4[i][j][k].Trapped21_to_Trapped22 = 11.0 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].Trapped22_to_Trapped21 = 11.0 / 1000.0;
        layer4[i][j][k].Trapped12_to_Trapped22 = 3.3 / (6.02 * pow(10, 5) * pow(0.5, 3));
        layer4[i][j][k].Trapped22_to_Trapped12 = 0.4 / 1000.0;

        layer4[i][j][k].Auton_to_CaMKII = 0.003 / 1000.0;
        layer4[i][j][k].Auton_to_Capped = 0.1 / 1000.0;
        layer4[i][j][k].Capped_to_Auton = 0.01 / 1000.0;	

	  layer4[i][j][k].CB00_to_CB10 = 5.5 / (6.02 * pow(10, 5) * pow(0.5, 3));
	  layer4[i][j][k].CB10_to_CB00 = 2.6 / 1000;
	  layer4[i][j][k].CB00_to_CB01 = 43.5 / (6.02 * pow(10, 5) * pow(0.5, 3));
	  layer4[i][j][k].CB01_to_CB00 = 35.8 / 1000;
	  layer4[i][j][k].CB10_to_CB11 = 43.5 / (6.02 * pow(10, 5) * pow(0.5, 3));
          layer4[i][j][k].CB11_to_CB10 = 35.8 / 1000;
	  layer4[i][j][k].CB10_to_CB20 = 5.5 / (6.02 * pow(10, 5) * pow(0.5, 3));
	  layer4[i][j][k].CB20_to_CB10 = 2.6 / 1000;
	  layer4[i][j][k].CB01_to_CB11 = 5.5 / (6.02 * pow(10, 5) * pow(0.5, 3));
	  layer4[i][j][k].CB11_to_CB01 = 2.6 / 1000;
	  layer4[i][j][k].CB01_to_CB02 = 43.5 / (6.02 * pow(10, 5) * pow(0.5, 3));
	  layer4[i][j][k].CB02_to_CB01 = 35.8 / 1000;
	  layer4[i][j][k].CB11_to_CB21 = 5.5 / (6.02 * pow(10, 5) * pow(0.5, 3));
	  layer4[i][j][k].CB21_to_CB11 = 2.6 / 1000;
	  layer4[i][j][k].CB11_to_CB12 = 43.5 / (6.02 * pow(10, 5) * pow(0.5, 3));
	  layer4[i][j][k].CB12_to_CB11 = 35.8 / 1000;
	  layer4[i][j][k].CB20_to_CB21 = 43.5 / (6.02 * pow(10, 5) * pow(0.5, 3));
	  layer4[i][j][k].CB21_to_CB20 = 35.8 / 1000;
	  layer4[i][j][k].CB02_to_CB12 = 5.5 / (6.02 * pow(10, 5) * pow(0.5, 3));
	  layer4[i][j][k].CB12_to_CB02 = 2.6 / 1000;
	  layer4[i][j][k].CB12_to_CB22 = 5.5 / (6.02 * pow(10, 5) * pow(0.5, 3));
	  layer4[i][j][k].CB22_to_CB12 = 2.6 / 1000;
	  layer4[i][j][k].CB21_to_CB22 = 43.5 / (6.02 * pow(10, 5) * pow(0.5, 3));
	  layer4[i][j][k].CB22_to_CB21 = 35.8 / 1000;

      }
  
  /*****************************************************************************/
  //initialize CaMKII
 
  for(i=0; i<10; i++)
    for(j=0; j<10; j++)
      layer1[i][j].head = NULL;

  for(i=0; i<5; i++)
    for(j=0; j<5; j++)
      for(k=0; k<5; k++)
	layer2[i][j][k].head = NULL;

  for(i=0; i<8; i++)
    layer3[i].head = NULL;
  
  for(i=0; i<4; i++)
    for(j=0; j<2; j++)
      for(k=0; k<2; k++)
	layer4[i][j][k].head = NULL;

//  initialize CaMKII number here, 83 or 40

   for(i=0; i<83; i++){
   random_number = 100 * ran2(idum);
   inter = (int) random_number;
    j = inter / 10;
    k = inter % 10;
    if(layer1[j][k].head == NULL){
      layer1[j][k].head = creat(i);
    }else{
      head = creat(i);
      head->next = layer1[j][k].head;
      layer1[j][k].head = head;
    }
 }

}

struct CaMKII * creat(int n){
  struct CaMKII * head;
  int j, k;
  
  head = (struct CaMKII *) malloc(sizeof(struct CaMKII));
  for(j=0; j<2; j++)
    for(k=0; k<6; k++)
      head->subunit[j][k] = 0;
  head->index = n;
  
  return(head);
}

double transition_lamda(){
  int i, j, k;
  double lamda;
  lamda = 0;

  for(i=0; i<10; i++)
    for(j=0; j<10; j++){
      layer1[i][j].diffusion_lamda = layer1_diffusion_lamda(i, j);
      layer1[i][j].chemical_lamda = voxel_chemical_lamda(layer1[i][j]);
      lamda += layer1[i][j].diffusion_lamda + layer1[i][j].chemical_lamda;
    }
 
  for(i=0; i<5; i++)
    for(j=0; j<5; j++)
      for(k=0; k<5; k++){
	layer2[i][j][k].diffusion_lamda = layer2_diffusion_lamda(i, j, k);
	layer2[i][j][k].chemical_lamda = voxel_chemical_lamda(layer2[i][j][k]);
	lamda += layer2[i][j][k].diffusion_lamda + layer2[i][j][k].chemical_lamda;
      }
  
  for(i=0; i<8; i++){
    layer3[i].diffusion_lamda = layer3_diffusion_lamda(i);
    layer3[i].chemical_lamda = voxel_chemical_lamda(layer3[i]);
    lamda += layer3[i].diffusion_lamda + layer3[i].chemical_lamda;
  }
    
  for(i=0; i<4; i++)
    for(j=0; j<2; j++)
      for(k=0; k<2; k++){
	layer4[i][j][k].diffusion_lamda = layer4_diffusion_lamda(i, j, k);
	layer4[i][j][k].chemical_lamda = voxel_chemical_lamda(layer4[i][j][k]);
	lamda += layer4[i][j][k].diffusion_lamda + layer4[i][j][k].chemical_lamda;
      }
	
  return lamda;
}

double voxel_chemical_lamda(struct voxel grid){
  double lamda;
  struct CaMKII *p1;
  int i, j;
	
  lamda = 0;
  
  lamda += grid.CaM00 * grid.Ca * grid.CaM00_to_CaM10;
  lamda += grid.CaM10 * grid.CaM10_to_CaM00;
  lamda += grid.CaM00 * grid.Ca * grid.CaM00_to_CaM01;
  lamda += grid.CaM01 * grid.CaM01_to_CaM00;
  lamda += grid.CaM10 * grid.Ca * grid.CaM10_to_CaM11;
  lamda += grid.CaM11 * grid.CaM11_to_CaM10;
  lamda += grid.CaM10 * grid.Ca * grid.CaM10_to_CaM20;
  lamda += grid.CaM20 * grid.CaM20_to_CaM10;
  lamda += grid.CaM01 * grid.Ca * grid.CaM01_to_CaM11;
  lamda += grid.CaM11 * grid.CaM11_to_CaM01;
  lamda += grid.CaM01 * grid.Ca * grid.CaM01_to_CaM02;
  lamda += grid.CaM02 * grid.CaM02_to_CaM01;
  lamda += grid.CaM11 * grid.Ca * grid.CaM11_to_CaM21;
  lamda += grid.CaM21 * grid.CaM21_to_CaM11;
  lamda += grid.CaM11 * grid.Ca * grid.CaM11_to_CaM12;
  lamda += grid.CaM12 * grid.CaM12_to_CaM11;
  lamda += grid.CaM20 * grid.Ca * grid.CaM20_to_CaM21;
  lamda += grid.CaM21 * grid.CaM21_to_CaM20;
  lamda += grid.CaM02 * grid.Ca * grid.CaM02_to_CaM12;
  lamda += grid.CaM12 * grid.CaM12_to_CaM02;
  lamda += grid.CaM21 * grid.Ca * grid.CaM21_to_CaM22;
  lamda += grid.CaM22 * grid.CaM22_to_CaM21;
  lamda += grid.CaM12 * grid.Ca * grid.CaM12_to_CaM22;
  lamda += grid.CaM22 * grid.CaM22_to_CaM12;

  lamda += grid.CaM00 * grid.CaN * grid.CaM00_to_CaNCaM00;
  lamda += grid.CaNCaM00 * grid.CaNCaM00_to_CaM00;
  lamda += grid.CaM10 * grid.CaN * grid.CaM10_to_CaNCaM10;
  lamda += grid.CaNCaM10 * grid.CaNCaM10_to_CaM10;
  lamda += grid.CaM01 * grid.CaN * grid.CaM01_to_CaNCaM01;
  lamda += grid.CaNCaM01 * grid.CaNCaM01_to_CaM01;
  lamda += grid.CaM11 * grid.CaN * grid.CaM11_to_CaNCaM11;
  lamda += grid.CaNCaM11 * grid.CaNCaM11_to_CaM11;  
  lamda += grid.CaM20 * grid.CaN * grid.CaM20_to_CaNCaM20;
  lamda += grid.CaNCaM20 * grid.CaNCaM20_to_CaM20;
  lamda += grid.CaM02 * grid.CaN * grid.CaM02_to_CaNCaM02;
  lamda += grid.CaNCaM02 * grid.CaNCaM02_to_CaM02;  
  lamda += grid.CaM12 * grid.CaN * grid.CaM12_to_CaNCaM12;
  lamda += grid.CaNCaM12 * grid.CaNCaM12_to_CaM12;
  lamda += grid.CaM21 * grid.CaN * grid.CaM21_to_CaNCaM21;
  lamda += grid.CaNCaM21 * grid.CaNCaM21_to_CaM21;
  lamda += grid.CaM22 * grid.CaN * grid.CaM22_to_CaNCaM22;
  lamda += grid.CaNCaM22 * grid.CaNCaM22_to_CaM22;

  lamda += grid.CaNCaM00 * grid.Ca * grid.CaNCaM00_to_CaNCaM10;
  lamda += grid.CaNCaM10 * grid.CaNCaM10_to_CaNCaM00;
  lamda += grid.CaNCaM00 * grid.Ca * grid.CaNCaM00_to_CaNCaM01;
  lamda += grid.CaNCaM01 * grid.CaNCaM01_to_CaNCaM00;
  lamda += grid.CaNCaM10 * grid.Ca * grid.CaNCaM10_to_CaNCaM11;
  lamda += grid.CaNCaM11 * grid.CaNCaM11_to_CaNCaM10;
  lamda += grid.CaNCaM10 * grid.Ca * grid.CaNCaM10_to_CaNCaM20;
  lamda += grid.CaNCaM20 * grid.CaNCaM20_to_CaNCaM10;
  lamda += grid.CaNCaM01 * grid.Ca * grid.CaNCaM01_to_CaNCaM11;
  lamda += grid.CaNCaM11 * grid.CaNCaM11_to_CaNCaM01;
  lamda += grid.CaNCaM01 * grid.Ca * grid.CaNCaM01_to_CaNCaM02;
  lamda += grid.CaNCaM02 * grid.CaNCaM02_to_CaNCaM01;
  lamda += grid.CaNCaM11 * grid.Ca * grid.CaNCaM11_to_CaNCaM21;
  lamda += grid.CaNCaM21 * grid.CaNCaM21_to_CaNCaM11;
  lamda += grid.CaNCaM11 * grid.Ca * grid.CaNCaM11_to_CaNCaM12;
  lamda += grid.CaNCaM21 * grid.CaNCaM21_to_CaNCaM11;
  lamda += grid.CaNCaM20 * grid.Ca * grid.CaNCaM20_to_CaNCaM21;
  lamda += grid.CaNCaM21 * grid.CaNCaM21_to_CaNCaM20;
  lamda += grid.CaNCaM02 * grid.Ca * grid.CaNCaM02_to_CaNCaM12;
  lamda += grid.CaNCaM12 * grid.CaNCaM12_to_CaNCaM02;
  lamda += grid.CaNCaM21 * grid.Ca * grid.CaNCaM21_to_CaNCaM22;
  lamda += grid.CaNCaM22 * grid.CaNCaM22_to_CaNCaM21;
  lamda += grid.CaNCaM12 * grid.Ca * grid.CaNCaM12_to_CaNCaM22;
  lamda += grid.CaNCaM22 * grid.CaNCaM22_to_CaNCaM12; 

  lamda += grid.CaM00 * grid.Ng * grid.CaM00_to_NgCaM00;
  lamda += grid.NgCaM00 * grid.NgCaM00_to_CaM00;
  lamda += grid.CaM10 * grid.Ng * grid.CaM10_to_NgCaM10;
  lamda += grid.NgCaM10 * grid.NgCaM10_to_CaM10;
  lamda += grid.CaM01 * grid.Ng * grid.CaM01_to_NgCaM01;
  lamda += grid.NgCaM01 * grid.NgCaM01_to_CaM01;
  lamda += grid.CaM11 * grid.Ng * grid.CaM11_to_NgCaM11;
  lamda += grid.NgCaM11 * grid.NgCaM11_to_CaM11;  
  lamda += grid.CaM20 * grid.Ng * grid.CaM20_to_NgCaM20;
  lamda += grid.NgCaM20 * grid.NgCaM20_to_CaM20;
  lamda += grid.CaM02 * grid.Ng * grid.CaM02_to_NgCaM02;
  lamda += grid.NgCaM02 * grid.NgCaM02_to_CaM02;  
  lamda += grid.CaM12 * grid.Ng * grid.CaM12_to_NgCaM12;
  lamda += grid.NgCaM12 * grid.NgCaM12_to_CaM12;
  lamda += grid.CaM21 * grid.Ng * grid.CaM21_to_NgCaM21;
  lamda += grid.NgCaM21 * grid.NgCaM21_to_CaM21;
  lamda += grid.CaM22 * grid.Ng * grid.CaM22_to_NgCaM22;
  lamda += grid.NgCaM22 * grid.NgCaM22_to_CaM22;

  lamda += grid.NgCaM00 * grid.Ca * grid.NgCaM00_to_NgCaM10;
  lamda += grid.NgCaM10 * grid.NgCaM10_to_NgCaM00;
  lamda += grid.NgCaM00 * grid.Ca * grid.NgCaM00_to_NgCaM01;
  lamda += grid.NgCaM01 * grid.NgCaM01_to_NgCaM00;
  lamda += grid.NgCaM10 * grid.Ca * grid.NgCaM10_to_NgCaM11;
  lamda += grid.NgCaM11 * grid.NgCaM11_to_NgCaM10;
  lamda += grid.NgCaM10 * grid.Ca * grid.NgCaM10_to_NgCaM20;
  lamda += grid.NgCaM20 * grid.NgCaM20_to_NgCaM10;
  lamda += grid.NgCaM01 * grid.Ca * grid.NgCaM01_to_NgCaM11;
  lamda += grid.NgCaM11 * grid.NgCaM11_to_NgCaM01;
  lamda += grid.NgCaM01 * grid.Ca * grid.NgCaM01_to_NgCaM02;
  lamda += grid.NgCaM02 * grid.NgCaM02_to_NgCaM01;
  lamda += grid.NgCaM11 * grid.Ca * grid.NgCaM11_to_NgCaM21;
  lamda += grid.NgCaM21 * grid.NgCaM21_to_NgCaM11;
  lamda += grid.NgCaM11 * grid.Ca * grid.NgCaM11_to_NgCaM12;
  lamda += grid.NgCaM21 * grid.NgCaM21_to_NgCaM11;
  lamda += grid.NgCaM20 * grid.Ca * grid.NgCaM20_to_NgCaM21;
  lamda += grid.NgCaM21 * grid.NgCaM21_to_NgCaM20;
  lamda += grid.NgCaM02 * grid.Ca * grid.NgCaM02_to_NgCaM12;
  lamda += grid.NgCaM12 * grid.NgCaM12_to_NgCaM02;
  lamda += grid.NgCaM21 * grid.Ca * grid.NgCaM21_to_NgCaM22;
  lamda += grid.NgCaM22 * grid.NgCaM22_to_NgCaM21;
  lamda += grid.NgCaM12 * grid.Ca * grid.NgCaM12_to_NgCaM22;
  lamda += grid.NgCaM22 * grid.NgCaM22_to_NgCaM12; 

  //chemical reaction reltated to CB
  lamda += grid.CB00 * grid.Ca * grid.CB00_to_CB10;
  lamda += grid.CB10 * grid.CB10_to_CB00;
  lamda += grid.CB00 * grid.Ca * grid.CB00_to_CB01;
  lamda += grid.CB01 * grid.CB01_to_CB00;
  lamda += grid.CB10 * grid.Ca * grid.CB10_to_CB11;
  lamda += grid.CB11 * grid.CB11_to_CB10;
  lamda += grid.CB10 * grid.Ca * grid.CB10_to_CB20;
  lamda += grid.CB20 * grid.CB20_to_CB10;
  lamda += grid.CB01 * grid.Ca * grid.CB01_to_CB11;
  lamda += grid.CB11 * grid.CB11_to_CB01;
  lamda += grid.CB01 * grid.Ca * grid.CB01_to_CB02;
  lamda += grid.CB02 * grid.CB02_to_CB01;
  lamda += grid.CB11 * grid.Ca * grid.CB11_to_CB21;
  lamda += grid.CB21 * grid.CB21_to_CB11;
  lamda += grid.CB11 * grid.Ca * grid.CB11_to_CB12;
  lamda += grid.CB12 * grid.CB12_to_CB11;
  lamda += grid.CB20 * grid.Ca * grid.CB20_to_CB21;
  lamda += grid.CB21 * grid.CB21_to_CB20;
  lamda += grid.CB02 * grid.Ca * grid.CB02_to_CB12;
  lamda += grid.CB12 * grid.CB12_to_CB02;
  lamda += grid.CB12 * grid.Ca * grid.CB12_to_CB22;
  lamda += grid.CB22 * grid.CB22_to_CB12;
  lamda += grid.CB21 * grid.Ca * grid.CB21_to_CB22;
  lamda += grid.CB22 * grid.CB22_to_CB21;

  //chemical reaction related to CaMKII
  p1 = grid.head;
  
  while(p1 != NULL){
    for(i=0; i<2; i++)
      for(j=0; j<6; j++){
	if(p1->subunit[i][j] ==0){
	  lamda += grid.CaM00 * grid.CaM00_to_CaMKIICaM00;
	  lamda += grid.CaM10 * grid.CaM10_to_CaMKIICaM10;
	  lamda += grid.CaM01 * grid.CaM01_to_CaMKIICaM01;
	  lamda += grid.CaM11 * grid.CaM11_to_CaMKIICaM11;
	  lamda += grid.CaM20 * grid.CaM20_to_CaMKIICaM20;
	  lamda += grid.CaM02 * grid.CaM02_to_CaMKIICaM02;
	  lamda += grid.CaM12 * grid.CaM12_to_CaMKIICaM12;
	  lamda += grid.CaM21 * grid.CaM21_to_CaMKIICaM21;
	  lamda += grid.CaM22 * grid.CaM22_to_CaMKIICaM22;
	}else if(p1->subunit[i][j] == 1){
	  lamda += grid.CaMKIICaM00_to_CaM00;
	  lamda += grid.Ca * grid.CaMKIICaM00_to_CaMKIICaM01;
	  lamda += grid.Ca * grid.CaMKIICaM00_to_CaMKIICaM10;
/*
	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))		 
         lamda += grid.CaMKIICaM00_to_Trapped00;
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
         (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM00_to_Trapped00;
*/
	}else if(p1->subunit[i][j] == 2){
	  lamda += grid.CaMKIICaM10_to_CaMKIICaM00;
	  lamda += grid.CaMKIICaM10_to_CaM10;
	  lamda += grid.Ca * grid.CaMKIICaM10_to_CaMKIICaM11;
          lamda += grid.Ca * grid.CaMKIICaM10_to_CaMKIICaM20;
/*
	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))		 
         lamda += grid.CaMKIICaM10_to_Trapped10;
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM10_to_Trapped10;
*/
	}else if(p1->subunit[i][j] == 3){
	  lamda += grid.CaMKIICaM01_to_CaMKIICaM00;
	  lamda += grid.CaMKIICaM01_to_CaM01;
	  lamda += grid.Ca * grid.CaMKIICaM01_to_CaMKIICaM11;
	  lamda += grid.Ca * grid.CaMKIICaM01_to_CaMKIICaM02;
/*
	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))		 
         lamda += grid.CaMKIICaM01_to_Trapped01;
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
         (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM01_to_Trapped01;
*/
	}else if(p1->subunit[i][j] == 4){
	  lamda += grid.CaMKIICaM11_to_CaMKIICaM01;
          lamda += grid.CaMKIICaM11_to_CaMKIICaM10;
	  lamda += grid.CaMKIICaM11_to_CaM11;
	  lamda += grid.Ca * grid.CaMKIICaM11_to_CaMKIICaM21;
	  lamda += grid.Ca * grid.CaMKIICaM11_to_CaMKIICaM12;
	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))		 
         lamda += grid.CaMKIICaM11_to_Trapped11;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM11_to_Trapped11;
*/
//
	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM11_to_Trapped11;
//
	}else if(p1->subunit[i][j] == 5){
	  lamda += grid.CaMKIICaM20_to_CaMKIICaM10;
	  lamda += grid.Ca * grid.CaMKIICaM20_to_CaMKIICaM21;
	  lamda += grid.CaMKIICaM20_to_CaM20;
	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))		 
         lamda += grid.CaMKIICaM20_to_Trapped20;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM20_to_Trapped20;
*/
//
	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM20_to_Trapped20;
//
	}else if(p1->subunit[i][j] == 6){
	  lamda += grid.CaMKIICaM02_to_CaMKIICaM01;
	  lamda += grid.Ca * grid.CaMKIICaM02_to_CaMKIICaM12;
	  lamda += grid.CaMKIICaM02_to_CaM02;
	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))		 
         lamda += grid.CaMKIICaM02_to_Trapped02;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM02_to_Trapped02;
*/
//
	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM02_to_Trapped02;
//
	}else if(p1->subunit[i][j] == 7){
	  lamda += grid.CaMKIICaM21_to_CaMKIICaM11;
          lamda += grid.CaMKIICaM21_to_CaMKIICaM20;
	  lamda += grid.Ca * grid.CaMKIICaM21_to_CaMKIICaM22;
	  lamda += grid.CaMKIICaM21_to_CaM21;
	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))		 
         lamda += grid.CaMKIICaM21_to_Trapped21;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM21_to_Trapped21;
*/
//
	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM21_to_Trapped21;
//	
	}else if(p1->subunit[i][j] == 8){
	  lamda += grid.CaMKIICaM12_to_CaMKIICaM11;
          lamda += grid.CaMKIICaM12_to_CaMKIICaM02;
	  lamda += grid.Ca * grid.CaMKIICaM12_to_CaMKIICaM22;
	  lamda += grid.CaMKIICaM12_to_CaM12;
	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))		 
         lamda += grid.CaMKIICaM12_to_Trapped12;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM12_to_Trapped12;
*/
//
	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM12_to_Trapped12;
//
	}else if(p1->subunit[i][j] == 9){
	  lamda += grid.CaMKIICaM22_to_CaMKIICaM21;
          lamda += grid.CaMKIICaM22_to_CaMKIICaM12;
	  lamda += grid.CaMKIICaM22_to_CaM22;
	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))		 
         lamda += grid.CaMKIICaM22_to_Trapped22;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM22_to_Trapped22;
*/
//
	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
		 lamda += 0.4 * grid.CaMKIICaM22_to_Trapped22;
//
	}else if(p1->subunit[i][j] == 10){
	  lamda += grid.Trapped00_to_CaMKIICaM00;
	  lamda += grid.Ca * grid.Trapped00_to_Trapped01;
	  lamda += grid.Ca * grid.Trapped00_to_Trapped10;
	  lamda += grid.Trapped00_to_Auton;
	 
	 }else if(p1->subunit[i][j] == 11){
	  lamda += grid.Trapped10_to_CaMKIICaM10;
	  lamda += grid.Trapped10_to_Trapped00;
          lamda += grid.Ca * grid.Trapped10_to_Trapped11;
	  lamda += grid.Ca * grid.Trapped10_to_Trapped20;
	  lamda += grid.Trapped10_to_Auton;

	}else if(p1->subunit[i][j] == 12){
	  lamda += grid.Trapped01_to_CaMKIICaM01;
	  lamda += grid.Trapped01_to_Trapped00;
          lamda += grid.Ca * grid.Trapped01_to_Trapped11;
	  lamda += grid.Ca * grid.Trapped01_to_Trapped02;
	  lamda += grid.Trapped01_to_Auton;

	}else if(p1->subunit[i][j] == 13){
	   lamda += grid.Trapped11_to_CaMKIICaM11;
	   lamda += grid.Trapped11_to_Trapped01;
	   lamda += grid.Trapped11_to_Trapped10;
	   lamda += grid.Ca * grid.Trapped11_to_Trapped21;
	   lamda += grid.Ca * grid.Trapped11_to_Trapped12;
	   lamda += grid.Trapped11_to_Auton;

	}else if(p1->subunit[i][j] == 14){
	   lamda += grid.Trapped20_to_CaMKIICaM20;
	   lamda += grid.Trapped20_to_Trapped10;
	   lamda += grid.Ca * grid.Trapped20_to_Trapped21;
	   lamda += grid.Trapped20_to_Auton;

	}else if(p1->subunit[i][j] == 15){
	   lamda += grid.Trapped02_to_CaMKIICaM02;
	   lamda += grid.Trapped02_to_Trapped01;
	   lamda += grid.Ca * grid.Trapped02_to_Trapped12;
	   lamda += grid.Trapped02_to_Auton;

	}else if(p1->subunit[i][j] == 16){
	   lamda += grid.Trapped21_to_CaMKIICaM21;
	   lamda += grid.Trapped21_to_Trapped11;
	   lamda += grid.Trapped21_to_Trapped20;
	   lamda += grid.Ca * grid.Trapped21_to_Trapped22;
	   lamda += grid.Trapped21_to_Auton;

	}else if(p1->subunit[i][j] == 17){
	   lamda += grid.Trapped12_to_CaMKIICaM12;
	   lamda += grid.Trapped12_to_Trapped11;
	   lamda += grid.Trapped12_to_Trapped02;
	   lamda += grid.Ca * grid.Trapped12_to_Trapped22;
	   lamda += grid.Trapped12_to_Auton;

	}else if(p1->subunit[i][j] == 18){
	   lamda += grid.Trapped22_to_CaMKIICaM22;
	   lamda += grid.Trapped22_to_Trapped21;
	   lamda += grid.Trapped22_to_Trapped12;
	   lamda += grid.Trapped22_to_Auton;

	}else if(p1->subunit[i][j] == 19){
	  lamda += grid.CaM00 * grid.Auton_to_Trapped00;
	  lamda += grid.CaM01 * grid.Auton_to_Trapped01;
	  lamda += grid.CaM10 * grid.Auton_to_Trapped10;
	  lamda += grid.CaM11 * grid.Auton_to_Trapped11;
	  lamda += grid.CaM20 * grid.Auton_to_Trapped20;
	  lamda += grid.CaM02 * grid.Auton_to_Trapped02;
	  lamda += grid.CaM21 * grid.Auton_to_Trapped21;
	  lamda += grid.CaM12 * grid.Auton_to_Trapped12;
	  lamda += grid.CaM22 * grid.Auton_to_Trapped22;
	  lamda += grid.Auton_to_CaMKII;

	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))		 
	     lamda += grid.Auton_to_Capped;
/*
	   if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	     lamda += 0.4 * grid.Auton_to_Capped;
*/
//
	   if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	     lamda += 0.4 * grid.Auton_to_Capped;
//
	}else{
	  lamda += grid.Capped_to_Auton;
	}
	   
      }
    p1 = p1->next;
  }
  return lamda;
}

double layer1_diffusion_lamda(int p, int q){
  double lamda;
  int k;
  int CaMKII_num;
  lamda = 0;
  
  CaMKII_num = CaMKII_number(layer1[p][q]);
  for(k=0; k<4; k++){
    lamda = lamda + CaMKII_num * layer1[p][q].D_CaMKII[k];
    lamda = lamda + layer1[p][q].Ca * layer1[p][q].D_Ca[k];

    lamda = lamda + layer1[p][q].CaM00 * layer1[p][q].D_CaM00[k];
    lamda = lamda + layer1[p][q].CaM10 * layer1[p][q].D_CaM10[k];
    lamda = lamda + layer1[p][q].CaM01 * layer1[p][q].D_CaM01[k];
    lamda = lamda + layer1[p][q].CaM11 * layer1[p][q].D_CaM11[k];
    lamda = lamda + layer1[p][q].CaM20 * layer1[p][q].D_CaM20[k];
	lamda = lamda + layer1[p][q].CaM02 * layer1[p][q].D_CaM02[k];
	lamda = lamda + layer1[p][q].CaM21 * layer1[p][q].D_CaM21[k];
	lamda = lamda + layer1[p][q].CaM12 * layer1[p][q].D_CaM12[k];
	lamda = lamda + layer1[p][q].CaM22 * layer1[p][q].D_CaM22[k];
	
	lamda = lamda + layer1[p][q].CB00 * layer1[p][q].D_CB00[k];
	lamda = lamda + layer1[p][q].CB01 * layer1[p][q].D_CB01[k];
	lamda = lamda + layer1[p][q].CB10 * layer1[p][q].D_CB10[k];
	lamda = lamda + layer1[p][q].CB11 * layer1[p][q].D_CB11[k];
	lamda = lamda + layer1[p][q].CB02 * layer1[p][q].D_CB02[k];
	lamda = lamda + layer1[p][q].CB20 * layer1[p][q].D_CB20[k];
	lamda = lamda + layer1[p][q].CB12 * layer1[p][q].D_CB12[k];
	lamda = lamda + layer1[p][q].CB21 * layer1[p][q].D_CB21[k];
	lamda = lamda + layer1[p][q].CB22 * layer1[p][q].D_CB22[k];
  }
  
  lamda = lamda + CaMKII_num * layer1[p][q].D_CaMKII[5];
  lamda = lamda + layer1[p][q].Ca * layer1[p][q].D_Ca[5];

  lamda = lamda + layer1[p][q].CaM00 * layer1[p][q].D_CaM00[5];
  lamda = lamda + layer1[p][q].CaM10 * layer1[p][q].D_CaM10[5];
  lamda = lamda + layer1[p][q].CaM01 * layer1[p][q].D_CaM01[5];
  lamda = lamda + layer1[p][q].CaM11 * layer1[p][q].D_CaM11[5];
  lamda = lamda + layer1[p][q].CaM20 * layer1[p][q].D_CaM20[5];
  lamda = lamda + layer1[p][q].CaM02 * layer1[p][q].D_CaM02[5];
  lamda = lamda + layer1[p][q].CaM21 * layer1[p][q].D_CaM21[5];
  lamda = lamda + layer1[p][q].CaM12 * layer1[p][q].D_CaM12[5];
  lamda = lamda + layer1[p][q].CaM22 * layer1[p][q].D_CaM22[5];
 
  lamda = lamda + layer1[p][q].CB00 * layer1[p][q].D_CB00[5];
  lamda = lamda + layer1[p][q].CB01 * layer1[p][q].D_CB01[5];
  lamda = lamda + layer1[p][q].CB10 * layer1[p][q].D_CB10[5];
  lamda = lamda + layer1[p][q].CB11 * layer1[p][q].D_CB11[5];
  lamda = lamda + layer1[p][q].CB02 * layer1[p][q].D_CB02[5];
  lamda = lamda + layer1[p][q].CB20 * layer1[p][q].D_CB20[5];
  lamda = lamda + layer1[p][q].CB12 * layer1[p][q].D_CB12[5];
  lamda = lamda + layer1[p][q].CB21 * layer1[p][q].D_CB21[5];
  lamda = lamda + layer1[p][q].CB22 * layer1[p][q].D_CB22[5];

  lamda = lamda + layer1[p][q].Ca_influx;
  layer1[p][q].Ca_pump = layer1[p][q].Ca * 0.0026 * layer1[p][q].pump_area / pow(0.05,3);
  lamda = lamda + layer1[p][q].Ca_pump;
  return lamda;
}

double layer2_diffusion_lamda(int p, int q, int r){
  double lamda;
  int k;
  int CaMKII_num;
  lamda = 0;
	   
  CaMKII_num = CaMKII_number(layer2[p][q][r]);

  for(k=0; k<6; k++){
    lamda = lamda + CaMKII_num * layer2[p][q][r].D_CaMKII[k];
    lamda = lamda + layer2[p][q][r].Ca * layer2[p][q][r].D_Ca[k];

    lamda = lamda + layer2[p][q][r].CaM00 * layer2[p][q][r].D_CaM00[k];
	lamda = lamda + layer2[p][q][r].CaM01 * layer2[p][q][r].D_CaM01[k];
	lamda = lamda + layer2[p][q][r].CaM10 * layer2[p][q][r].D_CaM10[k];
	lamda = lamda + layer2[p][q][r].CaM11 * layer2[p][q][r].D_CaM11[k];
	lamda = lamda + layer2[p][q][r].CaM20 * layer2[p][q][r].D_CaM20[k];
	lamda = lamda + layer2[p][q][r].CaM02 * layer2[p][q][r].D_CaM02[k];
	lamda = lamda + layer2[p][q][r].CaM21 * layer2[p][q][r].D_CaM21[k];
	lamda = lamda + layer2[p][q][r].CaM12 * layer2[p][q][r].D_CaM12[k];
	lamda = lamda + layer2[p][q][r].CaM22 * layer2[p][q][r].D_CaM22[k];
    
	lamda = lamda + layer2[p][q][r].CB00 * layer2[p][q][r].D_CB00[k];
	lamda = lamda + layer2[p][q][r].CB01 * layer2[p][q][r].D_CB01[k];
	lamda = lamda + layer2[p][q][r].CB10 * layer2[p][q][r].D_CB10[k];
	lamda = lamda + layer2[p][q][r].CB11 * layer2[p][q][r].D_CB11[k];
	lamda = lamda + layer2[p][q][r].CB02 * layer2[p][q][r].D_CB02[k];
	lamda = lamda + layer2[p][q][r].CB20 * layer2[p][q][r].D_CB20[k];
	lamda = lamda + layer2[p][q][r].CB12 * layer2[p][q][r].D_CB12[k];
	lamda = lamda + layer2[p][q][r].CB21 * layer2[p][q][r].D_CB21[k];
	lamda = lamda + layer2[p][q][r].CB22 * layer2[p][q][r].D_CB22[k];
  }
  
  if(r == 0){
    lamda = lamda + 3.0 * CaMKII_num * layer2[p][q][r].D_CaMKII[4];
    lamda = lamda + 3.0 * layer2[p][q][r].Ca * layer2[p][q][r].D_Ca[4];

	lamda = lamda + 3.0 * layer2[p][q][r].CaM00 * layer2[p][q][r].D_CaM00[4];
	lamda = lamda + 3.0 * layer2[p][q][r].CaM01 * layer2[p][q][r].D_CaM01[4];
	lamda = lamda + 3.0 * layer2[p][q][r].CaM10 * layer2[p][q][r].D_CaM10[4];
	lamda = lamda + 3.0 * layer2[p][q][r].CaM11 * layer2[p][q][r].D_CaM11[4];
	lamda = lamda + 3.0 * layer2[p][q][r].CaM20 * layer2[p][q][r].D_CaM20[4];
	lamda = lamda + 3.0 * layer2[p][q][r].CaM02 * layer2[p][q][r].D_CaM02[4];
	lamda = lamda + 3.0 * layer2[p][q][r].CaM21 * layer2[p][q][r].D_CaM21[4];
	lamda = lamda + 3.0 * layer2[p][q][r].CaM12 * layer2[p][q][r].D_CaM12[4];
	lamda = lamda + 3.0 * layer2[p][q][r].CaM22 * layer2[p][q][r].D_CaM22[4];

	lamda = lamda + 3.0 * layer2[p][q][r].CB00 * layer2[p][q][r].D_CB00[4];
	lamda = lamda + 3.0 * layer2[p][q][r].CB01 * layer2[p][q][r].D_CB01[4];
	lamda = lamda + 3.0 * layer2[p][q][r].CB10 * layer2[p][q][r].D_CB10[4];
	lamda = lamda + 3.0 * layer2[p][q][r].CB11 * layer2[p][q][r].D_CB11[4];
	lamda = lamda + 3.0 * layer2[p][q][r].CB02 * layer2[p][q][r].D_CB02[4];
	lamda = lamda + 3.0 * layer2[p][q][r].CB20 * layer2[p][q][r].D_CB20[4];
	lamda = lamda + 3.0 * layer2[p][q][r].CB12 * layer2[p][q][r].D_CB12[4];
	lamda = lamda + 3.0 * layer2[p][q][r].CB21 * layer2[p][q][r].D_CB21[4];
	lamda = lamda + 3.0 * layer2[p][q][r].CB22 * layer2[p][q][r].D_CB22[4];
  }

  layer2[p][q][r].Ca_pump = layer2[p][q][r].Ca * 0.0026 * layer2[p][q][r].pump_area / pow(0.1,3);
  lamda = lamda + layer2[p][q][r].Ca_pump;
  return lamda;
}

double layer3_diffusion_lamda(int p){
  double lamda;
  int k;
  lamda = 0;
	 
  for(k=4; k<6; k++){
    lamda = lamda + layer3[p].Ca * layer3[p].D_Ca[k];
    
	lamda = lamda + layer3[p].CaM00 * layer3[p].D_CaM00[k];
	lamda = lamda + layer3[p].CaM10 * layer3[p].D_CaM10[k];
	lamda = lamda + layer3[p].CaM01 * layer3[p].D_CaM01[k];
	lamda = lamda + layer3[p].CaM11 * layer3[p].D_CaM11[k];
	lamda = lamda + layer3[p].CaM20 * layer3[p].D_CaM20[k];
	lamda = lamda + layer3[p].CaM02 * layer3[p].D_CaM02[k];
	lamda = lamda + layer3[p].CaM21 * layer3[p].D_CaM21[k];
	lamda = lamda + layer3[p].CaM12 * layer3[p].D_CaM12[k];
	lamda = lamda + layer3[p].CaM22 * layer3[p].D_CaM22[k];

    lamda = lamda + layer3[p].CB00 * layer3[p].D_CB00[k];
	lamda = lamda + layer3[p].CB01 * layer3[p].D_CB01[k];
	lamda = lamda + layer3[p].CB10 * layer3[p].D_CB10[k];
	lamda = lamda + layer3[p].CB11 * layer3[p].D_CB11[k];
	lamda = lamda + layer3[p].CB02 * layer3[p].D_CB02[k];
	lamda = lamda + layer3[p].CB20 * layer3[p].D_CB20[k];
	lamda = lamda + layer3[p].CB12 * layer3[p].D_CB12[k];
	lamda = lamda + layer3[p].CB21 * layer3[p].D_CB21[k];
	lamda = lamda + layer3[p].CB22 * layer3[p].D_CB22[k];
  }

  layer3[p].Ca_pump = layer3[p].Ca * 0.0026 * layer3[p].pump_area / pow(0.1,3);
  lamda = lamda + layer3[p].Ca_pump;
  return lamda;
}

double layer4_diffusion_lamda(int p, int q, int r){
  double lamda;
  int k;
  lamda = 0;
      
  for(k=0; k<6; k++){
    lamda = lamda + layer4[p][q][r].Ca * layer4[p][q][r].D_Ca[k];

	lamda = lamda + layer4[p][q][r].CaM00 * layer4[p][q][r].D_CaM00[k];
	lamda = lamda + layer4[p][q][r].CaM10 * layer4[p][q][r].D_CaM10[k];
	lamda = lamda + layer4[p][q][r].CaM01 * layer4[p][q][r].D_CaM01[k];
	lamda = lamda + layer4[p][q][r].CaM11 * layer4[p][q][r].D_CaM11[k];
	lamda = lamda + layer4[p][q][r].CaM20 * layer4[p][q][r].D_CaM20[k];
	lamda = lamda + layer4[p][q][r].CaM02 * layer4[p][q][r].D_CaM02[k];
	lamda = lamda + layer4[p][q][r].CaM21 * layer4[p][q][r].D_CaM21[k];
	lamda = lamda + layer4[p][q][r].CaM12 * layer4[p][q][r].D_CaM12[k];
	lamda = lamda + layer4[p][q][r].CaM22 * layer4[p][q][r].D_CaM22[k];

	lamda = lamda + layer4[p][q][r].CB00 * layer4[p][q][r].D_CB00[k];
	lamda = lamda + layer4[p][q][r].CB01 * layer4[p][q][r].D_CB01[k];
	lamda = lamda + layer4[p][q][r].CB10 * layer4[p][q][r].D_CB10[k];
	lamda = lamda + layer4[p][q][r].CB11 * layer4[p][q][r].D_CB11[k];
	lamda = lamda + layer4[p][q][r].CB02 * layer4[p][q][r].D_CB02[k];
	lamda = lamda + layer4[p][q][r].CB20 * layer4[p][q][r].D_CB20[k];
	lamda = lamda + layer4[p][q][r].CB12 * layer4[p][q][r].D_CB12[k];
	lamda = lamda + layer4[p][q][r].CB21 * layer4[p][q][r].D_CB21[k];
	lamda = lamda + layer4[p][q][r].CB22 * layer4[p][q][r].D_CB22[k];
  }

  layer4[p][q][r].Ca_pump=0;
  if (layer4[p][q][r].Ca > 4) {
     layer4[p][q][r].Ca_pump = (layer4[p][q][r].Ca - 4) * 0.0026 * layer4[p][q][r].pump_area / pow(0.5,3);
     lamda = lamda + layer4[p][q][r].Ca_pump;
  }
  return lamda;
}

int CaMKII_number(struct voxel grid){
  struct CaMKII *p1;
  int n;
  
  n = 0;
  p1 = grid.head;
  while(p1 != NULL){
    n++;
    p1 = p1->next;
  }

  return n;
}

void reaction(){
  int i, j, k;
  double choose, choose1;

  random_num = ran2(idum) * reaction_lamda;
  choose = 0;

  for(i=0; i<10; i++){
    for(j=0; j<10; j++){
      choose1 = choose;
      choose = choose +  layer1[i][j].diffusion_lamda +  layer1[i][j].chemical_lamda;
      if((choose > random_num) && (choose1 <= random_num)){
	reaction_layer1(choose1, i, j);
	goto loop0;
      }
    }
  }

  for(i=0; i<5; i++){
    for(j=0; j<5; j++){
      for(k=0; k<5; k++){
	choose1 = choose;
	choose += layer2[i][j][k].diffusion_lamda + layer2[i][j][k].chemical_lamda;
	if((choose > random_num) && (choose1 <= random_num)){
	  reaction_layer2(choose1, i, j, k);
	  goto loop0;
	}
      }
    }
  }


  for(i=0; i<8; i++){
    choose1 = choose;
    choose += layer3[i].diffusion_lamda + layer3[i].chemical_lamda;
    if((choose > random_num) && (choose1 <= random_num)){
      reaction_layer3(choose1, i);
      goto loop0;
    }
  }

  for(i=0; i<4; i++){
    for(j=0; j<2; j++){
      for(k=0; k<2; k++){
	choose1 = choose;
	choose += layer4[i][j][k].diffusion_lamda + layer4[i][j][k].chemical_lamda;
	if((choose > random_num) && (choose1 <= random_num)){
	  reaction_layer4(choose1, i, j, k);
	  goto loop0;
	}
      }
    }
  }

 loop0: k = 0;
}

void reaction_layer1(double choose2, int p, int q){
  int k;
  struct CaMKII *p1;
  double choose, choose1, v1, v2, v3, v4;

  choose = choose2;

  choose1 = choose;
  choose += layer1[p][q].Ca_pump;
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].Ca--;
    total_Ca_pumped++;
    total_Ca_pumped1++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    reaction_lamda += layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda - v1 - v2;
    goto loop1;
  }

  choose1 = choose;
  choose += layer1[p][q].Ca_influx;
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].Ca++;
    total_Ca_influx++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    reaction_lamda += layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda - v1 - v2;
    goto loop1;
  }

  for(k=0; k<4; k++){
    p1 = layer1[p][q].head;
    while(p1 != NULL){
      choose1 = choose;
      choose += layer1[p][q].D_CaMKII[k];
      if((choose > random_num) && (choose1 <= random_num)){
	layer1[p][q].head = delete_CaMKII(layer1[p][q], p1);
	v1 = layer1[p][q].diffusion_lamda;
	v2 = layer1[p][q].chemical_lamda;
	layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
	layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);

	if(k == 0){	  
	  layer1[p-1][q].head = insert_CaMKII(layer1[p-1][q], p1);
	  v3 = layer1[p-1][q].diffusion_lamda;
	  v4 = layer1[p-1][q].chemical_lamda;
	  layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	  layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	  reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	    layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
	}else if(k == 1){
	  layer1[p][q-1].head = insert_CaMKII(layer1[p][q-1], p1);
	  v3 = layer1[p][q-1].diffusion_lamda;
	  v4 = layer1[p][q-1].chemical_lamda;
	  layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	  layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	  reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	    layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
	}else if(k == 2){
	  layer1[p+1][q].head = insert_CaMKII(layer1[p+1][q], p1);
	  v3 = layer1[p+1][q].diffusion_lamda;
	  v4 = layer1[p+1][q].chemical_lamda;
	  layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	  layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	  reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	    layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
	}else{
	  layer1[p][q+1].head = insert_CaMKII(layer1[p][q+1], p1);
	  v3 = layer1[p][q+1].diffusion_lamda;
	  v4 = layer1[p][q+1].chemical_lamda;
	  layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	  layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	  reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	    layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
	goto loop1;
      }
      p1 = p1->next;
    }

    choose1 = choose;
    choose += layer1[p][q].Ca * layer1[p][q].D_Ca[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].Ca--;
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);

      if(k == 0){
	layer1[p-1][q].Ca++;
	v3 = layer1[p-1][q].diffusion_lamda;
	v4 = layer1[p-1][q].chemical_lamda;
	layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	  layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer1[p][q-1].Ca++;
	v3 = layer1[p][q-1].diffusion_lamda;
	v4 = layer1[p][q-1].chemical_lamda;
	layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	  layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer1[p+1][q].Ca++;
	v3 = layer1[p+1][q].diffusion_lamda;
	v4 = layer1[p+1][q].chemical_lamda;
	layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	  layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else{
	layer1[p][q+1].Ca++;
	v3 = layer1[p][q+1].diffusion_lamda;
	v4 = layer1[p][q+1].chemical_lamda;
	layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	  layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
     
      goto loop1;
    }

    choose1 = choose;
    choose += layer1[p][q].CaM00 * layer1[p][q].D_CaM00[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].CaM00--;
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
      if(k == 0){
	layer1[p-1][q].CaM00++;
	v3 = layer1[p-1][q].diffusion_lamda;
	v4 = layer1[p-1][q].chemical_lamda;
	layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	  layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer1[p][q-1].CaM00++;
	v3 = layer1[p][q-1].diffusion_lamda;
	v4 = layer1[p][q-1].chemical_lamda;
	layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	  layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer1[p+1][q].CaM00++;
	v3 = layer1[p+1][q].diffusion_lamda;
	v4 = layer1[p+1][q].chemical_lamda;
	layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	  layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else{
	layer1[p][q+1].CaM00++;
	v3 = layer1[p][q+1].diffusion_lamda;
	v4 = layer1[p][q+1].chemical_lamda;
	layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	  layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += layer1[p][q].CaM10 * layer1[p][q].D_CaM10[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].CaM10--;
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
      if(k == 0){
	layer1[p-1][q].CaM10++;
	v3 = layer1[p-1][q].diffusion_lamda;
	v4 = layer1[p-1][q].chemical_lamda;
	layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	  layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer1[p][q-1].CaM10++;
	v3 = layer1[p][q-1].diffusion_lamda;
	v4 = layer1[p][q-1].chemical_lamda;
	layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	  layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer1[p+1][q].CaM10++;
	v3 = layer1[p+1][q].diffusion_lamda;
	v4 = layer1[p+1][q].chemical_lamda;
	layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	  layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else{
	layer1[p][q+1].CaM10++;
	v3 = layer1[p][q+1].diffusion_lamda;
	v4 = layer1[p][q+1].chemical_lamda;
	layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	  layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += layer1[p][q].CaM01 * layer1[p][q].D_CaM01[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].CaM01--;
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
      if(k == 0){
	layer1[p-1][q].CaM01++;
	v3 = layer1[p-1][q].diffusion_lamda;
	v4 = layer1[p-1][q].chemical_lamda;
	layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	  layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer1[p][q-1].CaM01++;
	v3 = layer1[p][q-1].diffusion_lamda;
	v4 = layer1[p][q-1].chemical_lamda;
	layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	  layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer1[p+1][q].CaM01++;
	v3 = layer1[p+1][q].diffusion_lamda;
	v4 = layer1[p+1][q].chemical_lamda;
	layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	  layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else{
	layer1[p][q+1].CaM01++;
	v3 = layer1[p][q+1].diffusion_lamda;
	v4 = layer1[p][q+1].chemical_lamda;
	layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	  layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += layer1[p][q].CaM11 * layer1[p][q].D_CaM11[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].CaM11--;
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
      if(k == 0){
	layer1[p-1][q].CaM11++;
	v3 = layer1[p-1][q].diffusion_lamda;
	v4 = layer1[p-1][q].chemical_lamda;
	layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	  layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer1[p][q-1].CaM11++;
	v3 = layer1[p][q-1].diffusion_lamda;
	v4 = layer1[p][q-1].chemical_lamda;
	layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	  layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer1[p+1][q].CaM11++;
	v3 = layer1[p+1][q].diffusion_lamda;
	v4 = layer1[p+1][q].chemical_lamda;
	layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	  layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else{
	layer1[p][q+1].CaM11++;
	v3 = layer1[p][q+1].diffusion_lamda;
	v4 = layer1[p][q+1].chemical_lamda;
	layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	  layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += layer1[p][q].CaM20 * layer1[p][q].D_CaM20[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].CaM20--;
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
      if(k == 0){
	layer1[p-1][q].CaM20++;
	v3 = layer1[p-1][q].diffusion_lamda;
	v4 = layer1[p-1][q].chemical_lamda;
	layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	  layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer1[p][q-1].CaM20++;
	v3 = layer1[p][q-1].diffusion_lamda;
	v4 = layer1[p][q-1].chemical_lamda;
	layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	  layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer1[p+1][q].CaM20++;
	v3 = layer1[p+1][q].diffusion_lamda;
	v4 = layer1[p+1][q].chemical_lamda;
	layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	  layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else{
	layer1[p][q+1].CaM20++;
	v3 = layer1[p][q+1].diffusion_lamda;
	v4 = layer1[p][q+1].chemical_lamda;
	layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	  layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += layer1[p][q].CaM02 * layer1[p][q].D_CaM02[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].CaM02--;
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
      if(k == 0){
	layer1[p-1][q].CaM02++;
	v3 = layer1[p-1][q].diffusion_lamda;
	v4 = layer1[p-1][q].chemical_lamda;
	layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	  layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer1[p][q-1].CaM02++;
	v3 = layer1[p][q-1].diffusion_lamda;
	v4 = layer1[p][q-1].chemical_lamda;
	layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	  layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer1[p+1][q].CaM02++;
	v3 = layer1[p+1][q].diffusion_lamda;
	v4 = layer1[p+1][q].chemical_lamda;
	layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	  layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else{
	layer1[p][q+1].CaM02++;
	v3 = layer1[p][q+1].diffusion_lamda;
	v4 = layer1[p][q+1].chemical_lamda;
	layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	  layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += layer1[p][q].CaM21 * layer1[p][q].D_CaM21[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].CaM21--;
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
      if(k == 0){
	layer1[p-1][q].CaM21++;
	v3 = layer1[p-1][q].diffusion_lamda;
	v4 = layer1[p-1][q].chemical_lamda;
	layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	  layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer1[p][q-1].CaM21++;
	v3 = layer1[p][q-1].diffusion_lamda;
	v4 = layer1[p][q-1].chemical_lamda;
	layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	  layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer1[p+1][q].CaM21++;
	v3 = layer1[p+1][q].diffusion_lamda;
	v4 = layer1[p+1][q].chemical_lamda;
	layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	  layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else{
	layer1[p][q+1].CaM21++;
	v3 = layer1[p][q+1].diffusion_lamda;
	v4 = layer1[p][q+1].chemical_lamda;
	layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	  layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += layer1[p][q].CaM12 * layer1[p][q].D_CaM12[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].CaM12--;
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
      if(k == 0){
	layer1[p-1][q].CaM12++;
	v3 = layer1[p-1][q].diffusion_lamda;
	v4 = layer1[p-1][q].chemical_lamda;
	layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	  layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer1[p][q-1].CaM12++;
	v3 = layer1[p][q-1].diffusion_lamda;
	v4 = layer1[p][q-1].chemical_lamda;
	layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	  layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer1[p+1][q].CaM12++;
	v3 = layer1[p+1][q].diffusion_lamda;
	v4 = layer1[p+1][q].chemical_lamda;
	layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	  layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else{
	layer1[p][q+1].CaM12++;
	v3 = layer1[p][q+1].diffusion_lamda;
	v4 = layer1[p][q+1].chemical_lamda;
	layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	  layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += layer1[p][q].CaM22 * layer1[p][q].D_CaM22[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].CaM22--;
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
      if(k == 0){
	layer1[p-1][q].CaM22++;
	v3 = layer1[p-1][q].diffusion_lamda;
	v4 = layer1[p-1][q].chemical_lamda;
	layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	  layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer1[p][q-1].CaM22++;
	v3 = layer1[p][q-1].diffusion_lamda;
	v4 = layer1[p][q-1].chemical_lamda;
	layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	  layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer1[p+1][q].CaM22++;
	v3 = layer1[p+1][q].diffusion_lamda;
	v4 = layer1[p+1][q].chemical_lamda;
	layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	  layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else{
	layer1[p][q+1].CaM22++;
	v3 = layer1[p][q+1].diffusion_lamda;
	v4 = layer1[p][q+1].chemical_lamda;
	layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	  layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }
    
	choose1 = choose;
    choose += layer1[p][q].CB00 * layer1[p][q].D_CB00[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].CB00--;
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
      
      if(k == 0){
	layer1[p-1][q].CB00++;
	v3 = layer1[p-1][q].diffusion_lamda;
	v4 = layer1[p-1][q].chemical_lamda;
	layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	  layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer1[p][q-1].CB00++;
	v3 = layer1[p][q-1].diffusion_lamda;
	v4 = layer1[p][q-1].chemical_lamda;
	layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	  layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer1[p+1][q].CB00++;
	v3 = layer1[p+1][q].diffusion_lamda;
	v4 = layer1[p+1][q].chemical_lamda;
	layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	  layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else{
	layer1[p][q+1].CB00++;
	v3 = layer1[p][q+1].diffusion_lamda;
	v4 = layer1[p][q+1].chemical_lamda;
	layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	  layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += layer1[p][q].CB01 * layer1[p][q].D_CB01[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].CB01--;
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
      
      if(k == 0){
	layer1[p-1][q].CB01++;
	v3 = layer1[p-1][q].diffusion_lamda;
	v4 = layer1[p-1][q].chemical_lamda;
	layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	  layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer1[p][q-1].CB01++;
	v3 = layer1[p][q-1].diffusion_lamda;
	v4 = layer1[p][q-1].chemical_lamda;
	layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	  layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer1[p+1][q].CB01++;
	v3 = layer1[p+1][q].diffusion_lamda;
	v4 = layer1[p+1][q].chemical_lamda;
	layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	  layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else{
	layer1[p][q+1].CB01++;
	v3 = layer1[p][q+1].diffusion_lamda;
	v4 = layer1[p][q+1].chemical_lamda;
	layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	  layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += layer1[p][q].CB10 * layer1[p][q].D_CB10[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].CB10--;
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
      
      if(k == 0){
	layer1[p-1][q].CB10++;
	v3 = layer1[p-1][q].diffusion_lamda;
	v4 = layer1[p-1][q].chemical_lamda;
	layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	  layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer1[p][q-1].CB10++;
	v3 = layer1[p][q-1].diffusion_lamda;
	v4 = layer1[p][q-1].chemical_lamda;
	layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	  layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer1[p+1][q].CB10++;
	v3 = layer1[p+1][q].diffusion_lamda;
	v4 = layer1[p+1][q].chemical_lamda;
	layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	  layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else{
	layer1[p][q+1].CB10++;
	v3 = layer1[p][q+1].diffusion_lamda;
	v4 = layer1[p][q+1].chemical_lamda;
	layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	  layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += layer1[p][q].CB11 * layer1[p][q].D_CB11[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].CB11--;
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
      
      if(k == 0){
	layer1[p-1][q].CB11++;
	v3 = layer1[p-1][q].diffusion_lamda;
	v4 = layer1[p-1][q].chemical_lamda;
	layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	  layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer1[p][q-1].CB11++;
	v3 = layer1[p][q-1].diffusion_lamda;
	v4 = layer1[p][q-1].chemical_lamda;
	layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	  layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer1[p+1][q].CB11++;
	v3 = layer1[p+1][q].diffusion_lamda;
	v4 = layer1[p+1][q].chemical_lamda;
	layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	  layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else{
	layer1[p][q+1].CB11++;
	v3 = layer1[p][q+1].diffusion_lamda;
	v4 = layer1[p][q+1].chemical_lamda;
	layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	  layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += layer1[p][q].CB20 * layer1[p][q].D_CB20[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].CB20--;
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
      
      if(k == 0){
	layer1[p-1][q].CB20++;
	v3 = layer1[p-1][q].diffusion_lamda;
	v4 = layer1[p-1][q].chemical_lamda;
	layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	  layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer1[p][q-1].CB20++;
	v3 = layer1[p][q-1].diffusion_lamda;
	v4 = layer1[p][q-1].chemical_lamda;
	layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	  layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer1[p+1][q].CB20++;
	v3 = layer1[p+1][q].diffusion_lamda;
	v4 = layer1[p+1][q].chemical_lamda;
	layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	  layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else{
	layer1[p][q+1].CB20++;
	v3 = layer1[p][q+1].diffusion_lamda;
	v4 = layer1[p][q+1].chemical_lamda;
	layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	  layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += layer1[p][q].CB02 * layer1[p][q].D_CB02[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].CB02--;
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
      
      if(k == 0){
	layer1[p-1][q].CB02++;
	v3 = layer1[p-1][q].diffusion_lamda;
	v4 = layer1[p-1][q].chemical_lamda;
	layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	  layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer1[p][q-1].CB02++;
	v3 = layer1[p][q-1].diffusion_lamda;
	v4 = layer1[p][q-1].chemical_lamda;
	layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	  layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer1[p+1][q].CB02++;
	v3 = layer1[p+1][q].diffusion_lamda;
	v4 = layer1[p+1][q].chemical_lamda;
	layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	  layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else{
	layer1[p][q+1].CB02++;
	v3 = layer1[p][q+1].diffusion_lamda;
	v4 = layer1[p][q+1].chemical_lamda;
	layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	  layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += layer1[p][q].CB21 * layer1[p][q].D_CB21[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].CB21--;
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
      
      if(k == 0){
	layer1[p-1][q].CB21++;
	v3 = layer1[p-1][q].diffusion_lamda;
	v4 = layer1[p-1][q].chemical_lamda;
	layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	  layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer1[p][q-1].CB21++;
	v3 = layer1[p][q-1].diffusion_lamda;
	v4 = layer1[p][q-1].chemical_lamda;
	layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	  layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer1[p+1][q].CB21++;
	v3 = layer1[p+1][q].diffusion_lamda;
	v4 = layer1[p+1][q].chemical_lamda;
	layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	  layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else{
	layer1[p][q+1].CB21++;
	v3 = layer1[p][q+1].diffusion_lamda;
	v4 = layer1[p][q+1].chemical_lamda;
	layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	  layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += layer1[p][q].CB12 * layer1[p][q].D_CB12[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].CB12--;
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
      
      if(k == 0){
	layer1[p-1][q].CB12++;
	v3 = layer1[p-1][q].diffusion_lamda;
	v4 = layer1[p-1][q].chemical_lamda;
	layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	  layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer1[p][q-1].CB12++;
	v3 = layer1[p][q-1].diffusion_lamda;
	v4 = layer1[p][q-1].chemical_lamda;
	layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	  layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer1[p+1][q].CB12++;
	v3 = layer1[p+1][q].diffusion_lamda;
	v4 = layer1[p+1][q].chemical_lamda;
	layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	  layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else{
	layer1[p][q+1].CB12++;
	v3 = layer1[p][q+1].diffusion_lamda;
	v4 = layer1[p][q+1].chemical_lamda;
	layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	  layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }

	choose1 = choose;
    choose += layer1[p][q].CB22 * layer1[p][q].D_CB22[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].CB22--;
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
      
      if(k == 0){
	layer1[p-1][q].CB22++;
	v3 = layer1[p-1][q].diffusion_lamda;
	v4 = layer1[p-1][q].chemical_lamda;
	layer1[p-1][q].diffusion_lamda = layer1_diffusion_lamda(p-1, q);
	layer1[p-1][q].chemical_lamda = voxel_chemical_lamda(layer1[p-1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p-1][q].diffusion_lamda +
	  layer1[p-1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer1[p][q-1].CB22++;
	v3 = layer1[p][q-1].diffusion_lamda;
	v4 = layer1[p][q-1].chemical_lamda;
	layer1[p][q-1].diffusion_lamda = layer1_diffusion_lamda(p, q-1);
	layer1[p][q-1].chemical_lamda = voxel_chemical_lamda(layer1[p][q-1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q-1].diffusion_lamda +
	  layer1[p][q-1].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer1[p+1][q].CB22++;
	v3 = layer1[p+1][q].diffusion_lamda;
	v4 = layer1[p+1][q].chemical_lamda;
	layer1[p+1][q].diffusion_lamda = layer1_diffusion_lamda(p+1, q);
	layer1[p+1][q].chemical_lamda = voxel_chemical_lamda(layer1[p+1][q]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p+1][q].diffusion_lamda +
	  layer1[p+1][q].chemical_lamda - v1 - v2 - v3 -v4;
      }else{
	layer1[p][q+1].CB22++;
	v3 = layer1[p][q+1].diffusion_lamda;
	v4 = layer1[p][q+1].chemical_lamda;
	layer1[p][q+1].diffusion_lamda = layer1_diffusion_lamda(p, q+1);
	layer1[p][q+1].chemical_lamda = voxel_chemical_lamda(layer1[p][q+1]);
	reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer1[p][q+1].diffusion_lamda +
	  layer1[p][q+1].chemical_lamda - v1 - v2 - v3 -v4;
      }
      goto loop1;
    }
  }

  p1 = layer1[p][q].head;
  while(p1 != NULL){
    choose1 = choose;
    choose += layer1[p][q].D_CaMKII[5];
    if((choose > random_num) && (choose1 <= random_num)){
      layer1[p][q].head = delete_CaMKII(layer1[p][q], p1);
      layer2[p/2][q/2][0].head = insert_CaMKII(layer2[p/2][q/2][0], p1);
      v1 = layer1[p][q].diffusion_lamda;
      v2 = layer1[p][q].chemical_lamda;
      layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
      layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
      v3 = layer2[p/2][q/2][0].diffusion_lamda;
      v4 = layer2[p/2][q/2][0].chemical_lamda;
      layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
      layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
      reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
	layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
      goto loop1;
    }
    p1 = p1->next;
  }
  

  choose1 = choose;
  choose += layer1[p][q].Ca * layer1[p][q].D_Ca[5];
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].Ca--;
    layer2[p/2][q/2][0].Ca++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    v3 = layer2[p/2][q/2][0].diffusion_lamda;
    v4 = layer2[p/2][q/2][0].chemical_lamda;
    layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
    layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
    reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
      layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
    goto loop1;
  }
  

  choose1 = choose;
  choose += layer1[p][q].CaM00 * layer1[p][q].D_CaM00[5];
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].CaM00--;
    layer2[p/2][q/2][0].CaM00++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    v3 = layer2[p/2][q/2][0].diffusion_lamda;
    v4 = layer2[p/2][q/2][0].chemical_lamda;
    layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
    layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
    reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
      layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
    goto loop1;
  }

  choose1 = choose;
  choose += layer1[p][q].CaM10 * layer1[p][q].D_CaM10[5];
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].CaM10--;
    layer2[p/2][q/2][0].CaM10++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    v3 = layer2[p/2][q/2][0].diffusion_lamda;
    v4 = layer2[p/2][q/2][0].chemical_lamda;
    layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
    layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
    reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
      layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
    goto loop1;
  }

  choose1 = choose;
  choose += layer1[p][q].CaM01 * layer1[p][q].D_CaM01[5];
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].CaM01--;
    layer2[p/2][q/2][0].CaM01++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    v3 = layer2[p/2][q/2][0].diffusion_lamda;
    v4 = layer2[p/2][q/2][0].chemical_lamda;
    layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
    layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
    reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
      layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
    goto loop1;
  }

  choose1 = choose;
  choose += layer1[p][q].CaM11 * layer1[p][q].D_CaM11[5];
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].CaM11--;
    layer2[p/2][q/2][0].CaM11++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    v3 = layer2[p/2][q/2][0].diffusion_lamda;
    v4 = layer2[p/2][q/2][0].chemical_lamda;
    layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
    layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
    reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
      layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
    goto loop1;
  }

  choose1 = choose;
  choose += layer1[p][q].CaM20 * layer1[p][q].D_CaM20[5];
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].CaM20--;
    layer2[p/2][q/2][0].CaM20++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    v3 = layer2[p/2][q/2][0].diffusion_lamda;
    v4 = layer2[p/2][q/2][0].chemical_lamda;
    layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
    layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
    reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
      layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
    goto loop1;
  }

  choose1 = choose;
  choose += layer1[p][q].CaM02 * layer1[p][q].D_CaM02[5];
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].CaM02--;
    layer2[p/2][q/2][0].CaM02++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    v3 = layer2[p/2][q/2][0].diffusion_lamda;
    v4 = layer2[p/2][q/2][0].chemical_lamda;
    layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
    layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
    reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
      layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
    goto loop1;
  }

  choose1 = choose;
  choose += layer1[p][q].CaM21 * layer1[p][q].D_CaM21[5];
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].CaM21--;
    layer2[p/2][q/2][0].CaM21++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    v3 = layer2[p/2][q/2][0].diffusion_lamda;
    v4 = layer2[p/2][q/2][0].chemical_lamda;
    layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
    layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
    reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
      layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
    goto loop1;
  }

  choose1 = choose;
  choose += layer1[p][q].CaM12 * layer1[p][q].D_CaM12[5];
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].CaM12--;
    layer2[p/2][q/2][0].CaM12++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    v3 = layer2[p/2][q/2][0].diffusion_lamda;
    v4 = layer2[p/2][q/2][0].chemical_lamda;
    layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
    layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
    reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
      layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
    goto loop1;
  }

  choose1 = choose;
  choose += layer1[p][q].CaM22 * layer1[p][q].D_CaM22[5];
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].CaM22--;
    layer2[p/2][q/2][0].CaM22++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    v3 = layer2[p/2][q/2][0].diffusion_lamda;
    v4 = layer2[p/2][q/2][0].chemical_lamda;
    layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
    layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
    reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
      layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
    goto loop1;
  }  
  
  choose1 = choose;
  choose += layer1[p][q].CB00 * layer1[p][q].D_CB00[5];
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].CB00--;
    layer2[p/2][q/2][0].CB00++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    v3 = layer2[p/2][q/2][0].diffusion_lamda;
    v4 = layer2[p/2][q/2][0].chemical_lamda;
    layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
    layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
    reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
      layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
    goto loop1;
  }

  choose1 = choose;
  choose += layer1[p][q].CB10 * layer1[p][q].D_CB10[5];
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].CB10--;
    layer2[p/2][q/2][0].CB10++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    v3 = layer2[p/2][q/2][0].diffusion_lamda;
    v4 = layer2[p/2][q/2][0].chemical_lamda;
    layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
    layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
    reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
      layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
    goto loop1;
  }

  choose1 = choose;
  choose += layer1[p][q].CB01 * layer1[p][q].D_CB01[5];
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].CB01--;
    layer2[p/2][q/2][0].CB01++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    v3 = layer2[p/2][q/2][0].diffusion_lamda;
    v4 = layer2[p/2][q/2][0].chemical_lamda;
    layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
    layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
    reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
      layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
    goto loop1;
  }

  choose1 = choose;
  choose += layer1[p][q].CB11 * layer1[p][q].D_CB11[5];
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].CB11--;
    layer2[p/2][q/2][0].CB11++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    v3 = layer2[p/2][q/2][0].diffusion_lamda;
    v4 = layer2[p/2][q/2][0].chemical_lamda;
    layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
    layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
    reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
      layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
    goto loop1;
  }

  choose1 = choose;
  choose += layer1[p][q].CB20 * layer1[p][q].D_CB20[5];
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].CB20--;
    layer2[p/2][q/2][0].CB20++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    v3 = layer2[p/2][q/2][0].diffusion_lamda;
    v4 = layer2[p/2][q/2][0].chemical_lamda;
    layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
    layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
    reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
      layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
    goto loop1;
  }

  choose1 = choose;
  choose += layer1[p][q].CB02 * layer1[p][q].D_CB02[5];
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].CB02--;
    layer2[p/2][q/2][0].CB02++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    v3 = layer2[p/2][q/2][0].diffusion_lamda;
    v4 = layer2[p/2][q/2][0].chemical_lamda;
    layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
    layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
    reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
      layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
    goto loop1;
  }

  choose1 = choose;
  choose += layer1[p][q].CB21 * layer1[p][q].D_CB21[5];
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].CB21--;
    layer2[p/2][q/2][0].CB21++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    v3 = layer2[p/2][q/2][0].diffusion_lamda;
    v4 = layer2[p/2][q/2][0].chemical_lamda;
    layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
    layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
    reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
      layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
    goto loop1;
  }

  choose1 = choose;
  choose += layer1[p][q].CB12 * layer1[p][q].D_CB12[5];
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].CB12--;
    layer2[p/2][q/2][0].CB12++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    v3 = layer2[p/2][q/2][0].diffusion_lamda;
    v4 = layer2[p/2][q/2][0].chemical_lamda;
    layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
    layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
    reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
      layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
    goto loop1;
  }

  choose1 = choose;
  choose += layer1[p][q].CB22 * layer1[p][q].D_CB22[5];
  if((choose > random_num) && (choose1 <= random_num)){
    layer1[p][q].CB22--;
    layer2[p/2][q/2][0].CB22++;
    v1 = layer1[p][q].diffusion_lamda;
    v2 = layer1[p][q].chemical_lamda;
    layer1[p][q].diffusion_lamda = layer1_diffusion_lamda(p, q);
    layer1[p][q].chemical_lamda = voxel_chemical_lamda(layer1[p][q]);
    v3 = layer2[p/2][q/2][0].diffusion_lamda;
    v4 = layer2[p/2][q/2][0].chemical_lamda;
    layer2[p/2][q/2][0].diffusion_lamda = layer2_diffusion_lamda(p/2, q/2, 0);
    layer2[p/2][q/2][0].chemical_lamda = voxel_chemical_lamda(layer2[p/2][q/2][0]);
    reaction_lamda = reaction_lamda + layer1[p][q].diffusion_lamda + layer1[p][q].chemical_lamda + layer2[p/2][q/2][0].diffusion_lamda +
      layer2[p/2][q/2][0].chemical_lamda - v1 - v2 - v3 -v4;
    goto loop1;
  }

  x = p;
  y = q;
  chemical_reaction(&layer1[p][q], choose, 1);
 loop1: k = 0;
}

void reaction_layer2(double choose2, int p, int q, int r){
  int k;
  struct CaMKII  *p1;
  double choose, choose1, v1, v2, v3, v4;
 
  choose = choose2;
  
  choose1 = choose;
  choose += layer2[p][q][r].Ca_pump;
  if((choose > random_num) && (choose1 <= random_num)){
    layer2[p][q][r].Ca--;
    total_Ca_pumped++;
    total_Ca_pumped2++;
    v1 = layer2[p][q][r].diffusion_lamda;
    v2 = layer2[p][q][r].chemical_lamda;
    layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
    layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);
    reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda - v1 - v2;
    goto loop2;
  }  

  for(k=0; k<6; k++){
    p1 = layer2[p][q][r].head;
    while(p1 != NULL){
      choose1 = choose;
      choose += layer2[p][q][r].D_CaMKII[k];
      if((choose > random_num) && (choose1 <= random_num)){
	layer2[p][q][r].head = delete_CaMKII(layer2[p][q][r], p1);
	v1 = layer2[p][q][r].diffusion_lamda;
	v2 = layer2[p][q][r].chemical_lamda;
	layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
	layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);

	if(k == 0){
	  layer2[p-1][q][r].head = insert_CaMKII(layer2[p-1][q][r], p1);
	  v3 = layer2[p-1][q][r].diffusion_lamda;
	  v4 = layer2[p-1][q][r].chemical_lamda;
	  layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	  layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	    layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
	}else if(k == 1){
	  layer2[p][q-1][r].head = insert_CaMKII(layer2[p][q-1][r], p1);
	  v3 = layer2[p][q-1][r].diffusion_lamda;
	  v4 = layer2[p][q-1][r].chemical_lamda;
	  layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	  layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	    layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
	}else if(k == 2){
	  layer2[p+1][q][r].head = insert_CaMKII(layer2[p+1][q][r], p1);
	  v3 = layer2[p+1][q][r].diffusion_lamda;
	  v4 = layer2[p+1][q][r].chemical_lamda;
	  layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	  layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	    layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
	}else if(k == 3){
	  layer2[p][q+1][r].head = insert_CaMKII(layer2[p][q+1][r], p1);
	  v3 = layer2[p][q+1][r].diffusion_lamda;
	  v4 = layer2[p][q+1][r].chemical_lamda;
	  layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	  layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	    layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
	}else if(k == 4){
	  if(r == 0){
	    layer1[2*p][2*q].head = insert_CaMKII(layer1[2*p][2*q], p1);
	    v3 = layer1[2*p][2*q].diffusion_lamda;
	    v4 = layer1[2*p][2*q].chemical_lamda;
	    layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	    layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	    reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	      layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	  }else{
	    layer2[p][q][r-1].head = insert_CaMKII(layer2[p][q][r-1], p1);
	    v3 = layer2[p][q][r-1].diffusion_lamda;
	    v4 = layer2[p][q][r-1].chemical_lamda;
	    layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	    layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	    reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	      layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	  }
	}else{
	  layer2[p][q][r+1].head = insert_CaMKII(layer2[p][q][r+1], p1);
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
	goto loop2;
      }
      p1 = p1->next;
    }
   
    choose1 = choose;
    choose += layer2[p][q][r].Ca * layer2[p][q][r].D_Ca[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][r].Ca--;
      v1 = layer2[p][q][r].diffusion_lamda;
      v2 = layer2[p][q][r].chemical_lamda;
      layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
      layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);
      if(k == 0){
	layer2[p-1][q][r].Ca++;
	v3 = layer2[p-1][q][r].diffusion_lamda;
	v4 = layer2[p-1][q][r].chemical_lamda;
	layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	  layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer2[p][q-1][r].Ca++;
	v3 = layer2[p][q-1][r].diffusion_lamda;
	v4 = layer2[p][q-1][r].chemical_lamda;
	layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	  layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer2[p+1][q][r].Ca++;
	v3 = layer2[p+1][q][r].diffusion_lamda;
	v4 = layer2[p+1][q][r].chemical_lamda;
	layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	  layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	layer2[p][q+1][r].Ca++;
	v3 = layer2[p][q+1][r].diffusion_lamda;
	v4 = layer2[p][q+1][r].chemical_lamda;
	layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	  layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	if(r == 0){
	  layer1[2*p][2*q].Ca++;
	  v3 = layer1[2*p][2*q].diffusion_lamda;
	  v4 = layer1[2*p][2*q].chemical_lamda;
	  layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	  layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	    layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r-1].Ca++;
	  v3 = layer2[p][q][r-1].diffusion_lamda;
	  v4 = layer2[p][q][r-1].chemical_lamda;
	  layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	  layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	    layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }else{
	if((p == 2) && (q == 2) && (r == 4)){
	  layer3[0].Ca++;
	  v3 = layer3[0].diffusion_lamda;
	  v4 = layer3[0].chemical_lamda;
	  layer3[0].diffusion_lamda = layer3_diffusion_lamda(0);
	  layer3[0].chemical_lamda = voxel_chemical_lamda(layer3[0]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer3[0].diffusion_lamda +
	    layer3[0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r+1].Ca++;
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][r].CaM00 * layer2[p][q][r].D_CaM00[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][r].CaM00--;
      v1 = layer2[p][q][r].diffusion_lamda;
      v2 = layer2[p][q][r].chemical_lamda;
      layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
      layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);

      if(k == 0){
	layer2[p-1][q][r].CaM00++;
	v3 = layer2[p-1][q][r].diffusion_lamda;
	v4 = layer2[p-1][q][r].chemical_lamda;
	layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	  layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer2[p][q-1][r].CaM00++;
	v3 = layer2[p][q-1][r].diffusion_lamda;
	v4 = layer2[p][q-1][r].chemical_lamda;
	layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	  layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer2[p+1][q][r].CaM00++;
	v3 = layer2[p+1][q][r].diffusion_lamda;
	v4 = layer2[p+1][q][r].chemical_lamda;
	layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	  layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	layer2[p][q+1][r].CaM00++;
	v3 = layer2[p][q+1][r].diffusion_lamda;
	v4 = layer2[p][q+1][r].chemical_lamda;
	layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	  layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	if(r == 0){
	  layer1[2*p][2*q].CaM00++;
	  v3 = layer1[2*p][2*q].diffusion_lamda;
	  v4 = layer1[2*p][2*q].chemical_lamda;
	  layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	  layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	    layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r-1].CaM00++;
	  v3 = layer2[p][q][r-1].diffusion_lamda;
	  v4 = layer2[p][q][r-1].chemical_lamda;
	  layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	  layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	    layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }else{
	if((p == 2) && (q == 2) && (r == 4)){
	  layer3[0].CaM00++;
	  v3 = layer3[0].diffusion_lamda;
	  v4 = layer3[0].chemical_lamda;
	  layer3[0].diffusion_lamda = layer3_diffusion_lamda(0);
	  layer3[0].chemical_lamda = voxel_chemical_lamda(layer3[0]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer3[0].diffusion_lamda +
	    layer3[0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r+1].CaM00++;
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }
      goto loop2;
    }

	choose1 = choose;
    choose += layer2[p][q][r].CaM10 * layer2[p][q][r].D_CaM10[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][r].CaM10--;
      v1 = layer2[p][q][r].diffusion_lamda;
      v2 = layer2[p][q][r].chemical_lamda;
      layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
      layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);

      if(k == 0){
	layer2[p-1][q][r].CaM10++;
	v3 = layer2[p-1][q][r].diffusion_lamda;
	v4 = layer2[p-1][q][r].chemical_lamda;
	layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	  layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer2[p][q-1][r].CaM10++;
	v3 = layer2[p][q-1][r].diffusion_lamda;
	v4 = layer2[p][q-1][r].chemical_lamda;
	layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	  layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer2[p+1][q][r].CaM10++;
	v3 = layer2[p+1][q][r].diffusion_lamda;
	v4 = layer2[p+1][q][r].chemical_lamda;
	layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	  layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	layer2[p][q+1][r].CaM10++;
	v3 = layer2[p][q+1][r].diffusion_lamda;
	v4 = layer2[p][q+1][r].chemical_lamda;
	layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	  layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	if(r == 0){
	  layer1[2*p][2*q].CaM10++;
	  v3 = layer1[2*p][2*q].diffusion_lamda;
	  v4 = layer1[2*p][2*q].chemical_lamda;
	  layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	  layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	    layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r-1].CaM10++;
	  v3 = layer2[p][q][r-1].diffusion_lamda;
	  v4 = layer2[p][q][r-1].chemical_lamda;
	  layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	  layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	    layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }else{
	if((p == 2) && (q == 2) && (r == 4)){
	  layer3[0].CaM10++;
	  v3 = layer3[0].diffusion_lamda;
	  v4 = layer3[0].chemical_lamda;
	  layer3[0].diffusion_lamda = layer3_diffusion_lamda(0);
	  layer3[0].chemical_lamda = voxel_chemical_lamda(layer3[0]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer3[0].diffusion_lamda +
	    layer3[0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r+1].CaM10++;
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }
      goto loop2;
    }

	choose1 = choose;
    choose += layer2[p][q][r].CaM01 * layer2[p][q][r].D_CaM01[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][r].CaM01--;
      v1 = layer2[p][q][r].diffusion_lamda;
      v2 = layer2[p][q][r].chemical_lamda;
      layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
      layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);

      if(k == 0){
	layer2[p-1][q][r].CaM01++;
	v3 = layer2[p-1][q][r].diffusion_lamda;
	v4 = layer2[p-1][q][r].chemical_lamda;
	layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	  layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer2[p][q-1][r].CaM01++;
	v3 = layer2[p][q-1][r].diffusion_lamda;
	v4 = layer2[p][q-1][r].chemical_lamda;
	layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	  layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer2[p+1][q][r].CaM01++;
	v3 = layer2[p+1][q][r].diffusion_lamda;
	v4 = layer2[p+1][q][r].chemical_lamda;
	layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	  layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	layer2[p][q+1][r].CaM01++;
	v3 = layer2[p][q+1][r].diffusion_lamda;
	v4 = layer2[p][q+1][r].chemical_lamda;
	layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	  layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	if(r == 0){
	  layer1[2*p][2*q].CaM01++;
	  v3 = layer1[2*p][2*q].diffusion_lamda;
	  v4 = layer1[2*p][2*q].chemical_lamda;
	  layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	  layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	    layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r-1].CaM01++;
	  v3 = layer2[p][q][r-1].diffusion_lamda;
	  v4 = layer2[p][q][r-1].chemical_lamda;
	  layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	  layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	    layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }else{
	if((p == 2) && (q == 2) && (r == 4)){
	  layer3[0].CaM01++;
	  v3 = layer3[0].diffusion_lamda;
	  v4 = layer3[0].chemical_lamda;
	  layer3[0].diffusion_lamda = layer3_diffusion_lamda(0);
	  layer3[0].chemical_lamda = voxel_chemical_lamda(layer3[0]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer3[0].diffusion_lamda +
	    layer3[0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r+1].CaM01++;
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }
      goto loop2;
    }

	choose1 = choose;
    choose += layer2[p][q][r].CaM11 * layer2[p][q][r].D_CaM11[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][r].CaM11--;
      v1 = layer2[p][q][r].diffusion_lamda;
      v2 = layer2[p][q][r].chemical_lamda;
      layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
      layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);

      if(k == 0){
	layer2[p-1][q][r].CaM11++;
	v3 = layer2[p-1][q][r].diffusion_lamda;
	v4 = layer2[p-1][q][r].chemical_lamda;
	layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	  layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer2[p][q-1][r].CaM11++;
	v3 = layer2[p][q-1][r].diffusion_lamda;
	v4 = layer2[p][q-1][r].chemical_lamda;
	layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	  layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer2[p+1][q][r].CaM11++;
	v3 = layer2[p+1][q][r].diffusion_lamda;
	v4 = layer2[p+1][q][r].chemical_lamda;
	layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	  layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	layer2[p][q+1][r].CaM11++;
	v3 = layer2[p][q+1][r].diffusion_lamda;
	v4 = layer2[p][q+1][r].chemical_lamda;
	layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	  layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	if(r == 0){
	  layer1[2*p][2*q].CaM11++;
	  v3 = layer1[2*p][2*q].diffusion_lamda;
	  v4 = layer1[2*p][2*q].chemical_lamda;
	  layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	  layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	    layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r-1].CaM11++;
	  v3 = layer2[p][q][r-1].diffusion_lamda;
	  v4 = layer2[p][q][r-1].chemical_lamda;
	  layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	  layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	    layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }else{
	if((p == 2) && (q == 2) && (r == 4)){
	  layer3[0].CaM11++;
	  v3 = layer3[0].diffusion_lamda;
	  v4 = layer3[0].chemical_lamda;
	  layer3[0].diffusion_lamda = layer3_diffusion_lamda(0);
	  layer3[0].chemical_lamda = voxel_chemical_lamda(layer3[0]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer3[0].diffusion_lamda +
	    layer3[0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r+1].CaM11++;
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }
      goto loop2;
    }

	choose1 = choose;
    choose += layer2[p][q][r].CaM20 * layer2[p][q][r].D_CaM20[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][r].CaM20--;
      v1 = layer2[p][q][r].diffusion_lamda;
      v2 = layer2[p][q][r].chemical_lamda;
      layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
      layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);

      if(k == 0){
	layer2[p-1][q][r].CaM20++;
	v3 = layer2[p-1][q][r].diffusion_lamda;
	v4 = layer2[p-1][q][r].chemical_lamda;
	layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	  layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer2[p][q-1][r].CaM20++;
	v3 = layer2[p][q-1][r].diffusion_lamda;
	v4 = layer2[p][q-1][r].chemical_lamda;
	layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	  layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer2[p+1][q][r].CaM20++;
	v3 = layer2[p+1][q][r].diffusion_lamda;
	v4 = layer2[p+1][q][r].chemical_lamda;
	layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	  layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	layer2[p][q+1][r].CaM20++;
	v3 = layer2[p][q+1][r].diffusion_lamda;
	v4 = layer2[p][q+1][r].chemical_lamda;
	layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	  layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	if(r == 0){
	  layer1[2*p][2*q].CaM20++;
	  v3 = layer1[2*p][2*q].diffusion_lamda;
	  v4 = layer1[2*p][2*q].chemical_lamda;
	  layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	  layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	    layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r-1].CaM20++;
	  v3 = layer2[p][q][r-1].diffusion_lamda;
	  v4 = layer2[p][q][r-1].chemical_lamda;
	  layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	  layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	    layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }else{
	if((p == 2) && (q == 2) && (r == 4)){
	  layer3[0].CaM20++;
	  v3 = layer3[0].diffusion_lamda;
	  v4 = layer3[0].chemical_lamda;
	  layer3[0].diffusion_lamda = layer3_diffusion_lamda(0);
	  layer3[0].chemical_lamda = voxel_chemical_lamda(layer3[0]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer3[0].diffusion_lamda +
	    layer3[0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r+1].CaM20++;
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }
      goto loop2;
    }

	choose1 = choose;
    choose += layer2[p][q][r].CaM02 * layer2[p][q][r].D_CaM02[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][r].CaM02--;
      v1 = layer2[p][q][r].diffusion_lamda;
      v2 = layer2[p][q][r].chemical_lamda;
      layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
      layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);

      if(k == 0){
	layer2[p-1][q][r].CaM02++;
	v3 = layer2[p-1][q][r].diffusion_lamda;
	v4 = layer2[p-1][q][r].chemical_lamda;
	layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	  layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer2[p][q-1][r].CaM02++;
	v3 = layer2[p][q-1][r].diffusion_lamda;
	v4 = layer2[p][q-1][r].chemical_lamda;
	layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	  layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer2[p+1][q][r].CaM02++;
	v3 = layer2[p+1][q][r].diffusion_lamda;
	v4 = layer2[p+1][q][r].chemical_lamda;
	layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	  layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	layer2[p][q+1][r].CaM02++;
	v3 = layer2[p][q+1][r].diffusion_lamda;
	v4 = layer2[p][q+1][r].chemical_lamda;
	layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	  layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	if(r == 0){
	  layer1[2*p][2*q].CaM02++;
	  v3 = layer1[2*p][2*q].diffusion_lamda;
	  v4 = layer1[2*p][2*q].chemical_lamda;
	  layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	  layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	    layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r-1].CaM02++;
	  v3 = layer2[p][q][r-1].diffusion_lamda;
	  v4 = layer2[p][q][r-1].chemical_lamda;
	  layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	  layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	    layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }else{
	if((p == 2) && (q == 2) && (r == 4)){
	  layer3[0].CaM02++;
	  v3 = layer3[0].diffusion_lamda;
	  v4 = layer3[0].chemical_lamda;
	  layer3[0].diffusion_lamda = layer3_diffusion_lamda(0);
	  layer3[0].chemical_lamda = voxel_chemical_lamda(layer3[0]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer3[0].diffusion_lamda +
	    layer3[0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r+1].CaM02++;
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }
      goto loop2;
    }

	choose1 = choose;
    choose += layer2[p][q][r].CaM21 * layer2[p][q][r].D_CaM21[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][r].CaM21--;
      v1 = layer2[p][q][r].diffusion_lamda;
      v2 = layer2[p][q][r].chemical_lamda;
      layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
      layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);

      if(k == 0){
	layer2[p-1][q][r].CaM21++;
	v3 = layer2[p-1][q][r].diffusion_lamda;
	v4 = layer2[p-1][q][r].chemical_lamda;
	layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	  layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer2[p][q-1][r].CaM21++;
	v3 = layer2[p][q-1][r].diffusion_lamda;
	v4 = layer2[p][q-1][r].chemical_lamda;
	layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	  layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer2[p+1][q][r].CaM21++;
	v3 = layer2[p+1][q][r].diffusion_lamda;
	v4 = layer2[p+1][q][r].chemical_lamda;
	layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	  layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	layer2[p][q+1][r].CaM21++;
	v3 = layer2[p][q+1][r].diffusion_lamda;
	v4 = layer2[p][q+1][r].chemical_lamda;
	layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	  layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	if(r == 0){
	  layer1[2*p][2*q].CaM21++;
	  v3 = layer1[2*p][2*q].diffusion_lamda;
	  v4 = layer1[2*p][2*q].chemical_lamda;
	  layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	  layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	    layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r-1].CaM21++;
	  v3 = layer2[p][q][r-1].diffusion_lamda;
	  v4 = layer2[p][q][r-1].chemical_lamda;
	  layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	  layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	    layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }else{
	if((p == 2) && (q == 2) && (r == 4)){
	  layer3[0].CaM21++;
	  v3 = layer3[0].diffusion_lamda;
	  v4 = layer3[0].chemical_lamda;
	  layer3[0].diffusion_lamda = layer3_diffusion_lamda(0);
	  layer3[0].chemical_lamda = voxel_chemical_lamda(layer3[0]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer3[0].diffusion_lamda +
	    layer3[0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r+1].CaM21++;
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }
      goto loop2;
    }

	choose1 = choose;
    choose += layer2[p][q][r].CaM12 * layer2[p][q][r].D_CaM12[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][r].CaM12--;
      v1 = layer2[p][q][r].diffusion_lamda;
      v2 = layer2[p][q][r].chemical_lamda;
      layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
      layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);

      if(k == 0){
	layer2[p-1][q][r].CaM12++;
	v3 = layer2[p-1][q][r].diffusion_lamda;
	v4 = layer2[p-1][q][r].chemical_lamda;
	layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	  layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer2[p][q-1][r].CaM12++;
	v3 = layer2[p][q-1][r].diffusion_lamda;
	v4 = layer2[p][q-1][r].chemical_lamda;
	layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	  layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer2[p+1][q][r].CaM12++;
	v3 = layer2[p+1][q][r].diffusion_lamda;
	v4 = layer2[p+1][q][r].chemical_lamda;
	layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	  layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	layer2[p][q+1][r].CaM12++;
	v3 = layer2[p][q+1][r].diffusion_lamda;
	v4 = layer2[p][q+1][r].chemical_lamda;
	layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	  layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	if(r == 0){
	  layer1[2*p][2*q].CaM12++;
	  v3 = layer1[2*p][2*q].diffusion_lamda;
	  v4 = layer1[2*p][2*q].chemical_lamda;
	  layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	  layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	    layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r-1].CaM12++;
	  v3 = layer2[p][q][r-1].diffusion_lamda;
	  v4 = layer2[p][q][r-1].chemical_lamda;
	  layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	  layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	    layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }else{
	if((p == 2) && (q == 2) && (r == 4)){
	  layer3[0].CaM12++;
	  v3 = layer3[0].diffusion_lamda;
	  v4 = layer3[0].chemical_lamda;
	  layer3[0].diffusion_lamda = layer3_diffusion_lamda(0);
	  layer3[0].chemical_lamda = voxel_chemical_lamda(layer3[0]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer3[0].diffusion_lamda +
	    layer3[0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r+1].CaM12++;
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }
      goto loop2;
    }

	choose1 = choose;
    choose += layer2[p][q][r].CaM22 * layer2[p][q][r].D_CaM22[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][r].CaM22--;
      v1 = layer2[p][q][r].diffusion_lamda;
      v2 = layer2[p][q][r].chemical_lamda;
      layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
      layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);

      if(k == 0){
	layer2[p-1][q][r].CaM22++;
	v3 = layer2[p-1][q][r].diffusion_lamda;
	v4 = layer2[p-1][q][r].chemical_lamda;
	layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	  layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer2[p][q-1][r].CaM22++;
	v3 = layer2[p][q-1][r].diffusion_lamda;
	v4 = layer2[p][q-1][r].chemical_lamda;
	layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	  layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer2[p+1][q][r].CaM22++;
	v3 = layer2[p+1][q][r].diffusion_lamda;
	v4 = layer2[p+1][q][r].chemical_lamda;
	layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	  layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	layer2[p][q+1][r].CaM22++;
	v3 = layer2[p][q+1][r].diffusion_lamda;
	v4 = layer2[p][q+1][r].chemical_lamda;
	layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	  layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	if(r == 0){
	  layer1[2*p][2*q].CaM22++;
	  v3 = layer1[2*p][2*q].diffusion_lamda;
	  v4 = layer1[2*p][2*q].chemical_lamda;
	  layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	  layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	    layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r-1].CaM22++;
	  v3 = layer2[p][q][r-1].diffusion_lamda;
	  v4 = layer2[p][q][r-1].chemical_lamda;
	  layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	  layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	    layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }else{
	if((p == 2) && (q == 2) && (r == 4)){
	  layer3[0].CaM22++;
	  v3 = layer3[0].diffusion_lamda;
	  v4 = layer3[0].chemical_lamda;
	  layer3[0].diffusion_lamda = layer3_diffusion_lamda(0);
	  layer3[0].chemical_lamda = voxel_chemical_lamda(layer3[0]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer3[0].diffusion_lamda +
	    layer3[0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r+1].CaM22++;
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }
      goto loop2;
    }
    
	choose1 = choose;
    choose += layer2[p][q][r].CB00 * layer2[p][q][r].D_CB00[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][r].CB00--;
      v1 = layer2[p][q][r].diffusion_lamda;
      v2 = layer2[p][q][r].chemical_lamda;
      layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
      layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);

      if(k == 0){
	layer2[p-1][q][r].CB00++;
	v3 = layer2[p-1][q][r].diffusion_lamda;
	v4 = layer2[p-1][q][r].chemical_lamda;
	layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	  layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer2[p][q-1][r].CB00++;
	v3 = layer2[p][q-1][r].diffusion_lamda;
	v4 = layer2[p][q-1][r].chemical_lamda;
	layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	  layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer2[p+1][q][r].CB00++;
	v3 = layer2[p+1][q][r].diffusion_lamda;
	v4 = layer2[p+1][q][r].chemical_lamda;
	layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	  layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	layer2[p][q+1][r].CB00++;
	v3 = layer2[p][q+1][r].diffusion_lamda;
	v4 = layer2[p][q+1][r].chemical_lamda;
	layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	  layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	if(r == 0){
	  layer1[2*p][2*q].CB00++;
	  v3 = layer1[2*p][2*q].diffusion_lamda;
	  v4 = layer1[2*p][2*q].chemical_lamda;
	  layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	  layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	    layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r-1].CB00++;
	  v3 = layer2[p][q][r-1].diffusion_lamda;
	  v4 = layer2[p][q][r-1].chemical_lamda;
	  layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	  layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	    layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }else{
	if((p == 2) && (q == 2) && (r == 4)){
	  layer3[0].CB00++;
	  v3 = layer3[0].diffusion_lamda;
	  v4 = layer3[0].chemical_lamda;
	  layer3[0].diffusion_lamda = layer3_diffusion_lamda(0);
	  layer3[0].chemical_lamda = voxel_chemical_lamda(layer3[0]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer3[0].diffusion_lamda +
	    layer3[0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r+1].CB00++;
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }
      goto loop2;
    }

	choose1 = choose;
    choose += layer2[p][q][r].CB10 * layer2[p][q][r].D_CB10[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][r].CB10--;
      v1 = layer2[p][q][r].diffusion_lamda;
      v2 = layer2[p][q][r].chemical_lamda;
      layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
      layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);

      if(k == 0){
	layer2[p-1][q][r].CB10++;
	v3 = layer2[p-1][q][r].diffusion_lamda;
	v4 = layer2[p-1][q][r].chemical_lamda;
	layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	  layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer2[p][q-1][r].CB10++;
	v3 = layer2[p][q-1][r].diffusion_lamda;
	v4 = layer2[p][q-1][r].chemical_lamda;
	layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	  layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer2[p+1][q][r].CB10++;
	v3 = layer2[p+1][q][r].diffusion_lamda;
	v4 = layer2[p+1][q][r].chemical_lamda;
	layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	  layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	layer2[p][q+1][r].CB10++;
	v3 = layer2[p][q+1][r].diffusion_lamda;
	v4 = layer2[p][q+1][r].chemical_lamda;
	layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	  layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	if(r == 0){
	  layer1[2*p][2*q].CB10++;
	  v3 = layer1[2*p][2*q].diffusion_lamda;
	  v4 = layer1[2*p][2*q].chemical_lamda;
	  layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	  layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	    layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r-1].CB10++;
	  v3 = layer2[p][q][r-1].diffusion_lamda;
	  v4 = layer2[p][q][r-1].chemical_lamda;
	  layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	  layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	    layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }else{
	if((p == 2) && (q == 2) && (r == 4)){
	  layer3[0].CB10++;
	  v3 = layer3[0].diffusion_lamda;
	  v4 = layer3[0].chemical_lamda;
	  layer3[0].diffusion_lamda = layer3_diffusion_lamda(0);
	  layer3[0].chemical_lamda = voxel_chemical_lamda(layer3[0]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer3[0].diffusion_lamda +
	    layer3[0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r+1].CB10++;
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }
      goto loop2;
    }

	choose1 = choose;
    choose += layer2[p][q][r].CB01 * layer2[p][q][r].D_CB01[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][r].CB01--;
      v1 = layer2[p][q][r].diffusion_lamda;
      v2 = layer2[p][q][r].chemical_lamda;
      layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
      layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);

      if(k == 0){
	layer2[p-1][q][r].CB01++;
	v3 = layer2[p-1][q][r].diffusion_lamda;
	v4 = layer2[p-1][q][r].chemical_lamda;
	layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	  layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer2[p][q-1][r].CB01++;
	v3 = layer2[p][q-1][r].diffusion_lamda;
	v4 = layer2[p][q-1][r].chemical_lamda;
	layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	  layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer2[p+1][q][r].CB01++;
	v3 = layer2[p+1][q][r].diffusion_lamda;
	v4 = layer2[p+1][q][r].chemical_lamda;
	layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	  layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	layer2[p][q+1][r].CB01++;
	v3 = layer2[p][q+1][r].diffusion_lamda;
	v4 = layer2[p][q+1][r].chemical_lamda;
	layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	  layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	if(r == 0){
	  layer1[2*p][2*q].CB01++;
	  v3 = layer1[2*p][2*q].diffusion_lamda;
	  v4 = layer1[2*p][2*q].chemical_lamda;
	  layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	  layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	    layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r-1].CB01++;
	  v3 = layer2[p][q][r-1].diffusion_lamda;
	  v4 = layer2[p][q][r-1].chemical_lamda;
	  layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	  layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	    layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }else{
	if((p == 2) && (q == 2) && (r == 4)){
	  layer3[0].CB01++;
	  v3 = layer3[0].diffusion_lamda;
	  v4 = layer3[0].chemical_lamda;
	  layer3[0].diffusion_lamda = layer3_diffusion_lamda(0);
	  layer3[0].chemical_lamda = voxel_chemical_lamda(layer3[0]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer3[0].diffusion_lamda +
	    layer3[0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r+1].CB01++;
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }
      goto loop2;
    }

	choose1 = choose;
    choose += layer2[p][q][r].CB11 * layer2[p][q][r].D_CB11[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][r].CB11--;
      v1 = layer2[p][q][r].diffusion_lamda;
      v2 = layer2[p][q][r].chemical_lamda;
      layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
      layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);

      if(k == 0){
	layer2[p-1][q][r].CB11++;
	v3 = layer2[p-1][q][r].diffusion_lamda;
	v4 = layer2[p-1][q][r].chemical_lamda;
	layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	  layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer2[p][q-1][r].CB11++;
	v3 = layer2[p][q-1][r].diffusion_lamda;
	v4 = layer2[p][q-1][r].chemical_lamda;
	layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	  layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer2[p+1][q][r].CB11++;
	v3 = layer2[p+1][q][r].diffusion_lamda;
	v4 = layer2[p+1][q][r].chemical_lamda;
	layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	  layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	layer2[p][q+1][r].CB11++;
	v3 = layer2[p][q+1][r].diffusion_lamda;
	v4 = layer2[p][q+1][r].chemical_lamda;
	layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	  layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	if(r == 0){
	  layer1[2*p][2*q].CB11++;
	  v3 = layer1[2*p][2*q].diffusion_lamda;
	  v4 = layer1[2*p][2*q].chemical_lamda;
	  layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	  layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	    layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r-1].CB11++;
	  v3 = layer2[p][q][r-1].diffusion_lamda;
	  v4 = layer2[p][q][r-1].chemical_lamda;
	  layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	  layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	    layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }else{
	if((p == 2) && (q == 2) && (r == 4)){
	  layer3[0].CB11++;
	  v3 = layer3[0].diffusion_lamda;
	  v4 = layer3[0].chemical_lamda;
	  layer3[0].diffusion_lamda = layer3_diffusion_lamda(0);
	  layer3[0].chemical_lamda = voxel_chemical_lamda(layer3[0]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer3[0].diffusion_lamda +
	    layer3[0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r+1].CB11++;
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }
      goto loop2;
    }

	choose1 = choose;
    choose += layer2[p][q][r].CB20 * layer2[p][q][r].D_CB20[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][r].CB20--;
      v1 = layer2[p][q][r].diffusion_lamda;
      v2 = layer2[p][q][r].chemical_lamda;
      layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
      layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);

      if(k == 0){
	layer2[p-1][q][r].CB20++;
	v3 = layer2[p-1][q][r].diffusion_lamda;
	v4 = layer2[p-1][q][r].chemical_lamda;
	layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	  layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer2[p][q-1][r].CB20++;
	v3 = layer2[p][q-1][r].diffusion_lamda;
	v4 = layer2[p][q-1][r].chemical_lamda;
	layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	  layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer2[p+1][q][r].CB20++;
	v3 = layer2[p+1][q][r].diffusion_lamda;
	v4 = layer2[p+1][q][r].chemical_lamda;
	layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	  layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	layer2[p][q+1][r].CB20++;
	v3 = layer2[p][q+1][r].diffusion_lamda;
	v4 = layer2[p][q+1][r].chemical_lamda;
	layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	  layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	if(r == 0){
	  layer1[2*p][2*q].CB20++;
	  v3 = layer1[2*p][2*q].diffusion_lamda;
	  v4 = layer1[2*p][2*q].chemical_lamda;
	  layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	  layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	    layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r-1].CB20++;
	  v3 = layer2[p][q][r-1].diffusion_lamda;
	  v4 = layer2[p][q][r-1].chemical_lamda;
	  layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	  layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	    layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }else{
	if((p == 2) && (q == 2) && (r == 4)){
	  layer3[0].CB20++;
	  v3 = layer3[0].diffusion_lamda;
	  v4 = layer3[0].chemical_lamda;
	  layer3[0].diffusion_lamda = layer3_diffusion_lamda(0);
	  layer3[0].chemical_lamda = voxel_chemical_lamda(layer3[0]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer3[0].diffusion_lamda +
	    layer3[0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r+1].CB20++;
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }
      goto loop2;
    }

	choose1 = choose;
    choose += layer2[p][q][r].CB02 * layer2[p][q][r].D_CB02[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][r].CB02--;
      v1 = layer2[p][q][r].diffusion_lamda;
      v2 = layer2[p][q][r].chemical_lamda;
      layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
      layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);

      if(k == 0){
	layer2[p-1][q][r].CB02++;
	v3 = layer2[p-1][q][r].diffusion_lamda;
	v4 = layer2[p-1][q][r].chemical_lamda;
	layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	  layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer2[p][q-1][r].CB02++;
	v3 = layer2[p][q-1][r].diffusion_lamda;
	v4 = layer2[p][q-1][r].chemical_lamda;
	layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	  layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer2[p+1][q][r].CB02++;
	v3 = layer2[p+1][q][r].diffusion_lamda;
	v4 = layer2[p+1][q][r].chemical_lamda;
	layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	  layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	layer2[p][q+1][r].CB02++;
	v3 = layer2[p][q+1][r].diffusion_lamda;
	v4 = layer2[p][q+1][r].chemical_lamda;
	layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	  layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	if(r == 0){
	  layer1[2*p][2*q].CB02++;
	  v3 = layer1[2*p][2*q].diffusion_lamda;
	  v4 = layer1[2*p][2*q].chemical_lamda;
	  layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	  layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	    layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r-1].CB02++;
	  v3 = layer2[p][q][r-1].diffusion_lamda;
	  v4 = layer2[p][q][r-1].chemical_lamda;
	  layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	  layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	    layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }else{
	if((p == 2) && (q == 2) && (r == 4)){
	  layer3[0].CB02++;
	  v3 = layer3[0].diffusion_lamda;
	  v4 = layer3[0].chemical_lamda;
	  layer3[0].diffusion_lamda = layer3_diffusion_lamda(0);
	  layer3[0].chemical_lamda = voxel_chemical_lamda(layer3[0]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer3[0].diffusion_lamda +
	    layer3[0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r+1].CB02++;
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }
      goto loop2;
    }

	choose1 = choose;
    choose += layer2[p][q][r].CB21 * layer2[p][q][r].D_CB21[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][r].CB21--;
      v1 = layer2[p][q][r].diffusion_lamda;
      v2 = layer2[p][q][r].chemical_lamda;
      layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
      layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);

      if(k == 0){
	layer2[p-1][q][r].CB21++;
	v3 = layer2[p-1][q][r].diffusion_lamda;
	v4 = layer2[p-1][q][r].chemical_lamda;
	layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	  layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer2[p][q-1][r].CB21++;
	v3 = layer2[p][q-1][r].diffusion_lamda;
	v4 = layer2[p][q-1][r].chemical_lamda;
	layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	  layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer2[p+1][q][r].CB21++;
	v3 = layer2[p+1][q][r].diffusion_lamda;
	v4 = layer2[p+1][q][r].chemical_lamda;
	layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	  layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	layer2[p][q+1][r].CB21++;
	v3 = layer2[p][q+1][r].diffusion_lamda;
	v4 = layer2[p][q+1][r].chemical_lamda;
	layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	  layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	if(r == 0){
	  layer1[2*p][2*q].CB21++;
	  v3 = layer1[2*p][2*q].diffusion_lamda;
	  v4 = layer1[2*p][2*q].chemical_lamda;
	  layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	  layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	    layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r-1].CB21++;
	  v3 = layer2[p][q][r-1].diffusion_lamda;
	  v4 = layer2[p][q][r-1].chemical_lamda;
	  layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	  layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	    layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }else{
	if((p == 2) && (q == 2) && (r == 4)){
	  layer3[0].CB21++;
	  v3 = layer3[0].diffusion_lamda;
	  v4 = layer3[0].chemical_lamda;
	  layer3[0].diffusion_lamda = layer3_diffusion_lamda(0);
	  layer3[0].chemical_lamda = voxel_chemical_lamda(layer3[0]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer3[0].diffusion_lamda +
	    layer3[0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r+1].CB21++;
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }
      goto loop2;
    }

	choose1 = choose;
    choose += layer2[p][q][r].CB12 * layer2[p][q][r].D_CB12[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][r].CB12--;
      v1 = layer2[p][q][r].diffusion_lamda;
      v2 = layer2[p][q][r].chemical_lamda;
      layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
      layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);

      if(k == 0){
	layer2[p-1][q][r].CB12++;
	v3 = layer2[p-1][q][r].diffusion_lamda;
	v4 = layer2[p-1][q][r].chemical_lamda;
	layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	  layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer2[p][q-1][r].CB12++;
	v3 = layer2[p][q-1][r].diffusion_lamda;
	v4 = layer2[p][q-1][r].chemical_lamda;
	layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	  layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer2[p+1][q][r].CB12++;
	v3 = layer2[p+1][q][r].diffusion_lamda;
	v4 = layer2[p+1][q][r].chemical_lamda;
	layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	  layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	layer2[p][q+1][r].CB12++;
	v3 = layer2[p][q+1][r].diffusion_lamda;
	v4 = layer2[p][q+1][r].chemical_lamda;
	layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	  layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	if(r == 0){
	  layer1[2*p][2*q].CB12++;
	  v3 = layer1[2*p][2*q].diffusion_lamda;
	  v4 = layer1[2*p][2*q].chemical_lamda;
	  layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	  layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	    layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r-1].CB12++;
	  v3 = layer2[p][q][r-1].diffusion_lamda;
	  v4 = layer2[p][q][r-1].chemical_lamda;
	  layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	  layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	    layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }else{
	if((p == 2) && (q == 2) && (r == 4)){
	  layer3[0].CB12++;
	  v3 = layer3[0].diffusion_lamda;
	  v4 = layer3[0].chemical_lamda;
	  layer3[0].diffusion_lamda = layer3_diffusion_lamda(0);
	  layer3[0].chemical_lamda = voxel_chemical_lamda(layer3[0]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer3[0].diffusion_lamda +
	    layer3[0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r+1].CB12++;
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }
      goto loop2;
    }

	choose1 = choose;
    choose += layer2[p][q][r].CB22 * layer2[p][q][r].D_CB22[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][r].CB22--;
      v1 = layer2[p][q][r].diffusion_lamda;
      v2 = layer2[p][q][r].chemical_lamda;
      layer2[p][q][r].diffusion_lamda = layer2_diffusion_lamda(p, q, r);
      layer2[p][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r]);

      if(k == 0){
	layer2[p-1][q][r].CB22++;
	v3 = layer2[p-1][q][r].diffusion_lamda;
	v4 = layer2[p-1][q][r].chemical_lamda;
	layer2[p-1][q][r].diffusion_lamda = layer2_diffusion_lamda(p-1, q, r);
	layer2[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p-1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p-1][q][r].diffusion_lamda + 
	  layer2[p-1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 1){
	layer2[p][q-1][r].CB22++;
	v3 = layer2[p][q-1][r].diffusion_lamda;
	v4 = layer2[p][q-1][r].chemical_lamda;
	layer2[p][q-1][r].diffusion_lamda = layer2_diffusion_lamda(p, q-1, r);
	layer2[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q-1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q-1][r].diffusion_lamda + 
	  layer2[p][q-1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 2){
	layer2[p+1][q][r].CB22++;
	v3 = layer2[p+1][q][r].diffusion_lamda;
	v4 = layer2[p+1][q][r].chemical_lamda;
	layer2[p+1][q][r].diffusion_lamda = layer2_diffusion_lamda(p+1, q, r);
	layer2[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer2[p+1][q][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p+1][q][r].diffusion_lamda + 
	  layer2[p+1][q][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 3){
	layer2[p][q+1][r].CB22++;
	v3 = layer2[p][q+1][r].diffusion_lamda;
	v4 = layer2[p][q+1][r].chemical_lamda;
	layer2[p][q+1][r].diffusion_lamda = layer2_diffusion_lamda(p, q+1, r);
	layer2[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer2[p][q+1][r]);
	reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q+1][r].diffusion_lamda + 
	  layer2[p][q+1][r].chemical_lamda - v1 - v2 - v3 -v4;
      }else if(k == 4){
	if(r == 0){
	  layer1[2*p][2*q].CB22++;
	  v3 = layer1[2*p][2*q].diffusion_lamda;
	  v4 = layer1[2*p][2*q].chemical_lamda;
	  layer1[2*p][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q);
	  layer1[2*p][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer1[2*p][2*q].diffusion_lamda +
	    layer1[2*p][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r-1].CB22++;
	  v3 = layer2[p][q][r-1].diffusion_lamda;
	  v4 = layer2[p][q][r-1].chemical_lamda;
	  layer2[p][q][r-1].diffusion_lamda = layer2_diffusion_lamda(p, q, r-1);
	  layer2[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r-1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r-1].diffusion_lamda + 
	    layer2[p][q][r-1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }else{
	if((p == 2) && (q == 2) && (r == 4)){
	  layer3[0].CB22++;
	  v3 = layer3[0].diffusion_lamda;
	  v4 = layer3[0].chemical_lamda;
	  layer3[0].diffusion_lamda = layer3_diffusion_lamda(0);
	  layer3[0].chemical_lamda = voxel_chemical_lamda(layer3[0]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer3[0].diffusion_lamda +
	    layer3[0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer2[p][q][r+1].CB22++;
	  v3 = layer2[p][q][r+1].diffusion_lamda;
	  v4 = layer2[p][q][r+1].chemical_lamda;
	  layer2[p][q][r+1].diffusion_lamda = layer2_diffusion_lamda(p, q, r+1);
	  layer2[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer2[p][q][r+1]);
	  reaction_lamda += layer2[p][q][r].diffusion_lamda + layer2[p][q][r].chemical_lamda + layer2[p][q][r+1].diffusion_lamda + 
	    layer2[p][q][r+1].chemical_lamda - v1 - v2 - v3 -v4;
	}
      }
      goto loop2;
    }

  }

  if(r == 0){
    p1 = layer2[p][q][0].head;
    while(p1 != NULL){
      choose1 = choose;
      choose += layer2[p][q][0].D_CaMKII[4];
      if((choose > random_num) && (choose1 <= random_num)){
	layer2[p][q][0].head = delete_CaMKII(layer2[p][q][0], p1);
	layer1[2*p][2*q+1].head = insert_CaMKII(layer1[2*p][2*q+1], p1);
	v1 = layer2[p][q][0].diffusion_lamda;
	v2 = layer2[p][q][0].chemical_lamda;
	layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
	layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
	v3 = layer1[2*p][2*q+1].diffusion_lamda;
	v4 = layer1[2*p][2*q+1].chemical_lamda;
	layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
	layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
	reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	  layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
	goto loop2;
      }

      choose1 = choose;
      choose += layer2[p][q][0].D_CaMKII[4];
      if((choose > random_num) && (choose1 <= random_num)){
	layer2[p][q][0].head = delete_CaMKII(layer2[p][q][0], p1);
	layer1[2*p+1][2*q].head = insert_CaMKII(layer1[2*p+1][2*q], p1);
	v1 = layer2[p][q][0].diffusion_lamda;
	v2 = layer2[p][q][0].chemical_lamda;
	layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
	layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
	v3 = layer1[2*p+1][2*q].diffusion_lamda;
	v4 = layer1[2*p+1][2*q].chemical_lamda;
	layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
	layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
	reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	  layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
	goto loop2;
      }
      
      choose1 = choose;
      choose += layer2[p][q][0].D_CaMKII[4];
      if((choose > random_num) && (choose1 <= random_num)){
	layer2[p][q][0].head = delete_CaMKII(layer2[p][q][0], p1);
	layer1[2*p+1][2*q+1].head = insert_CaMKII(layer1[2*p+1][2*q+1], p1);
	v1 = layer2[p][q][0].diffusion_lamda;
	v2 = layer2[p][q][0].chemical_lamda;
	layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
	layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
	v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
	v4 = layer1[2*p+1][2*q+1].chemical_lamda;
	layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
	layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
	reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	  layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
	goto loop2;
      }
      
      p1 = p1->next;
    }
    /***************************************************/  

  
    choose1 = choose;
    choose += layer2[p][q][0].Ca * layer2[p][q][0].D_Ca[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].Ca--;
      layer1[2*p][2*q+1].Ca++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p][2*q+1].diffusion_lamda;
      v4 = layer1[2*p][2*q+1].chemical_lamda;
      layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
      layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].Ca * layer2[p][q][0].D_Ca[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].Ca--;
      layer1[2*p+1][2*q].Ca++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q].diffusion_lamda;
      v4 = layer1[2*p+1][2*q].chemical_lamda;
      layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
      layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }
  
    choose1 = choose;
    choose += layer2[p][q][0].Ca * layer2[p][q][0].D_Ca[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].Ca--;
      layer1[2*p+1][2*q+1].Ca++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
      v4 = layer1[2*p+1][2*q+1].chemical_lamda;
      layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
      layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }
    /*******************************************************/
    
    choose1 = choose;
    choose += layer2[p][q][0].CaM00 * layer2[p][q][0].D_CaM00[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM00--;
      layer1[2*p][2*q+1].CaM00++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p][2*q+1].diffusion_lamda;
      v4 = layer1[2*p][2*q+1].chemical_lamda;
      layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
      layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CaM00 * layer2[p][q][0].D_CaM00[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM00--;
      layer1[2*p+1][2*q].CaM00++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q].diffusion_lamda;
      v4 = layer1[2*p+1][2*q].chemical_lamda;
      layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
      layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CaM00 * layer2[p][q][0].D_CaM00[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM00--;
      layer1[2*p+1][2*q+1].CaM00++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
      v4 = layer1[2*p+1][2*q+1].chemical_lamda;
      layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
      layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    /*******************************************************/
    
    choose1 = choose;
    choose += layer2[p][q][0].CaM10 * layer2[p][q][0].D_CaM10[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM10--;
      layer1[2*p][2*q+1].CaM10++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p][2*q+1].diffusion_lamda;
      v4 = layer1[2*p][2*q+1].chemical_lamda;
      layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
      layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CaM10 * layer2[p][q][0].D_CaM10[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM10--;
      layer1[2*p+1][2*q].CaM10++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q].diffusion_lamda;
      v4 = layer1[2*p+1][2*q].chemical_lamda;
      layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
      layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CaM10 * layer2[p][q][0].D_CaM10[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM10--;
      layer1[2*p+1][2*q+1].CaM10++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
      v4 = layer1[2*p+1][2*q+1].chemical_lamda;
      layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
      layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }
    
	/*******************************************************/
    
    choose1 = choose;
    choose += layer2[p][q][0].CaM01 * layer2[p][q][0].D_CaM01[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM01--;
      layer1[2*p][2*q+1].CaM01++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p][2*q+1].diffusion_lamda;
      v4 = layer1[2*p][2*q+1].chemical_lamda;
      layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
      layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CaM01 * layer2[p][q][0].D_CaM01[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM01--;
      layer1[2*p+1][2*q].CaM01++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q].diffusion_lamda;
      v4 = layer1[2*p+1][2*q].chemical_lamda;
      layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
      layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CaM01 * layer2[p][q][0].D_CaM01[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM01--;
      layer1[2*p+1][2*q+1].CaM01++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
      v4 = layer1[2*p+1][2*q+1].chemical_lamda;
      layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
      layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }
    
	/*******************************************************/
    
    choose1 = choose;
    choose += layer2[p][q][0].CaM11 * layer2[p][q][0].D_CaM11[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM11--;
      layer1[2*p][2*q+1].CaM11++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p][2*q+1].diffusion_lamda;
      v4 = layer1[2*p][2*q+1].chemical_lamda;
      layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
      layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CaM11 * layer2[p][q][0].D_CaM11[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM11--;
      layer1[2*p+1][2*q].CaM11++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q].diffusion_lamda;
      v4 = layer1[2*p+1][2*q].chemical_lamda;
      layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
      layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CaM11 * layer2[p][q][0].D_CaM11[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM11--;
      layer1[2*p+1][2*q+1].CaM11++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
      v4 = layer1[2*p+1][2*q+1].chemical_lamda;
      layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
      layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }
    
	/*******************************************************/
    
    choose1 = choose;
    choose += layer2[p][q][0].CaM20 * layer2[p][q][0].D_CaM20[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM20--;
      layer1[2*p][2*q+1].CaM20++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p][2*q+1].diffusion_lamda;
      v4 = layer1[2*p][2*q+1].chemical_lamda;
      layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
      layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CaM20 * layer2[p][q][0].D_CaM20[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM20--;
      layer1[2*p+1][2*q].CaM20++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q].diffusion_lamda;
      v4 = layer1[2*p+1][2*q].chemical_lamda;
      layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
      layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CaM20 * layer2[p][q][0].D_CaM20[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM20--;
      layer1[2*p+1][2*q+1].CaM20++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
      v4 = layer1[2*p+1][2*q+1].chemical_lamda;
      layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
      layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }
    
	/*******************************************************/
    
    choose1 = choose;
    choose += layer2[p][q][0].CaM02 * layer2[p][q][0].D_CaM02[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM02--;
      layer1[2*p][2*q+1].CaM02++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p][2*q+1].diffusion_lamda;
      v4 = layer1[2*p][2*q+1].chemical_lamda;
      layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
      layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CaM02 * layer2[p][q][0].D_CaM02[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM02--;
      layer1[2*p+1][2*q].CaM02++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q].diffusion_lamda;
      v4 = layer1[2*p+1][2*q].chemical_lamda;
      layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
      layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CaM02 * layer2[p][q][0].D_CaM02[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM02--;
      layer1[2*p+1][2*q+1].CaM02++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
      v4 = layer1[2*p+1][2*q+1].chemical_lamda;
      layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
      layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }
    
	/*******************************************************/
    
    choose1 = choose;
    choose += layer2[p][q][0].CaM21 * layer2[p][q][0].D_CaM21[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM21--;
      layer1[2*p][2*q+1].CaM21++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p][2*q+1].diffusion_lamda;
      v4 = layer1[2*p][2*q+1].chemical_lamda;
      layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
      layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CaM21 * layer2[p][q][0].D_CaM21[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM21--;
      layer1[2*p+1][2*q].CaM21++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q].diffusion_lamda;
      v4 = layer1[2*p+1][2*q].chemical_lamda;
      layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
      layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CaM21 * layer2[p][q][0].D_CaM21[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM21--;
      layer1[2*p+1][2*q+1].CaM21++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
      v4 = layer1[2*p+1][2*q+1].chemical_lamda;
      layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
      layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }
    
	/*******************************************************/
    
    choose1 = choose;
    choose += layer2[p][q][0].CaM12 * layer2[p][q][0].D_CaM12[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM12--;
      layer1[2*p][2*q+1].CaM12++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p][2*q+1].diffusion_lamda;
      v4 = layer1[2*p][2*q+1].chemical_lamda;
      layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
      layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CaM12 * layer2[p][q][0].D_CaM12[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM12--;
      layer1[2*p+1][2*q].CaM12++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q].diffusion_lamda;
      v4 = layer1[2*p+1][2*q].chemical_lamda;
      layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
      layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CaM12 * layer2[p][q][0].D_CaM12[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM12--;
      layer1[2*p+1][2*q+1].CaM12++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
      v4 = layer1[2*p+1][2*q+1].chemical_lamda;
      layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
      layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }
    
	/*******************************************************/
    
    choose1 = choose;
    choose += layer2[p][q][0].CaM22 * layer2[p][q][0].D_CaM22[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM22--;
      layer1[2*p][2*q+1].CaM22++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p][2*q+1].diffusion_lamda;
      v4 = layer1[2*p][2*q+1].chemical_lamda;
      layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
      layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CaM22 * layer2[p][q][0].D_CaM22[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM22--;
      layer1[2*p+1][2*q].CaM22++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q].diffusion_lamda;
      v4 = layer1[2*p+1][2*q].chemical_lamda;
      layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
      layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CaM22 * layer2[p][q][0].D_CaM22[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CaM22--;
      layer1[2*p+1][2*q+1].CaM22++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
      v4 = layer1[2*p+1][2*q+1].chemical_lamda;
      layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
      layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }
    
    
	/*******************************************************/
    
    choose1 = choose;
    choose += layer2[p][q][0].CB00 * layer2[p][q][0].D_CB00[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB00--;
      layer1[2*p][2*q+1].CB00++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p][2*q+1].diffusion_lamda;
      v4 = layer1[2*p][2*q+1].chemical_lamda;
      layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
      layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CB00 * layer2[p][q][0].D_CB00[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB00--;
      layer1[2*p+1][2*q].CB00++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q].diffusion_lamda;
      v4 = layer1[2*p+1][2*q].chemical_lamda;
      layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
      layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CB00 * layer2[p][q][0].D_CB00[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB00--;
      layer1[2*p+1][2*q+1].CB00++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
      v4 = layer1[2*p+1][2*q+1].chemical_lamda;
      layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
      layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

	/*******************************************************/
    
    choose1 = choose;
    choose += layer2[p][q][0].CB10 * layer2[p][q][0].D_CB10[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB10--;
      layer1[2*p][2*q+1].CB10++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p][2*q+1].diffusion_lamda;
      v4 = layer1[2*p][2*q+1].chemical_lamda;
      layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
      layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CB10 * layer2[p][q][0].D_CB10[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB10--;
      layer1[2*p+1][2*q].CB10++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q].diffusion_lamda;
      v4 = layer1[2*p+1][2*q].chemical_lamda;
      layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
      layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CB10 * layer2[p][q][0].D_CB10[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB10--;
      layer1[2*p+1][2*q+1].CB10++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
      v4 = layer1[2*p+1][2*q+1].chemical_lamda;
      layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
      layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

	/*******************************************************/
    
    choose1 = choose;
    choose += layer2[p][q][0].CB01 * layer2[p][q][0].D_CB01[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB01--;
      layer1[2*p][2*q+1].CB01++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p][2*q+1].diffusion_lamda;
      v4 = layer1[2*p][2*q+1].chemical_lamda;
      layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
      layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CB01 * layer2[p][q][0].D_CB01[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB01--;
      layer1[2*p+1][2*q].CB01++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q].diffusion_lamda;
      v4 = layer1[2*p+1][2*q].chemical_lamda;
      layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
      layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CB01 * layer2[p][q][0].D_CB01[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB01--;
      layer1[2*p+1][2*q+1].CB01++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
      v4 = layer1[2*p+1][2*q+1].chemical_lamda;
      layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
      layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }
 
	/*******************************************************/
    
    choose1 = choose;
    choose += layer2[p][q][0].CB11 * layer2[p][q][0].D_CB11[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB11--;
      layer1[2*p][2*q+1].CB11++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p][2*q+1].diffusion_lamda;
      v4 = layer1[2*p][2*q+1].chemical_lamda;
      layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
      layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CB11 * layer2[p][q][0].D_CB11[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB11--;
      layer1[2*p+1][2*q].CB11++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q].diffusion_lamda;
      v4 = layer1[2*p+1][2*q].chemical_lamda;
      layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
      layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CB11 * layer2[p][q][0].D_CB11[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB11--;
      layer1[2*p+1][2*q+1].CB11++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
      v4 = layer1[2*p+1][2*q+1].chemical_lamda;
      layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
      layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }


	/*******************************************************/
    
    choose1 = choose;
    choose += layer2[p][q][0].CB20 * layer2[p][q][0].D_CB20[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB20--;
      layer1[2*p][2*q+1].CB20++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p][2*q+1].diffusion_lamda;
      v4 = layer1[2*p][2*q+1].chemical_lamda;
      layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
      layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CB20 * layer2[p][q][0].D_CB20[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB20--;
      layer1[2*p+1][2*q].CB20++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q].diffusion_lamda;
      v4 = layer1[2*p+1][2*q].chemical_lamda;
      layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
      layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CB20 * layer2[p][q][0].D_CB20[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB20--;
      layer1[2*p+1][2*q+1].CB20++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
      v4 = layer1[2*p+1][2*q+1].chemical_lamda;
      layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
      layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

	/*******************************************************/
    
    choose1 = choose;
    choose += layer2[p][q][0].CB02 * layer2[p][q][0].D_CB02[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB02--;
      layer1[2*p][2*q+1].CB02++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p][2*q+1].diffusion_lamda;
      v4 = layer1[2*p][2*q+1].chemical_lamda;
      layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
      layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CB02 * layer2[p][q][0].D_CB02[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB02--;
      layer1[2*p+1][2*q].CB02++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q].diffusion_lamda;
      v4 = layer1[2*p+1][2*q].chemical_lamda;
      layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
      layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CB02 * layer2[p][q][0].D_CB02[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB02--;
      layer1[2*p+1][2*q+1].CB02++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
      v4 = layer1[2*p+1][2*q+1].chemical_lamda;
      layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
      layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

	/*******************************************************/
    
    choose1 = choose;
    choose += layer2[p][q][0].CB21 * layer2[p][q][0].D_CB21[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB21--;
      layer1[2*p][2*q+1].CB21++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p][2*q+1].diffusion_lamda;
      v4 = layer1[2*p][2*q+1].chemical_lamda;
      layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
      layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CB21 * layer2[p][q][0].D_CB21[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB21--;
      layer1[2*p+1][2*q].CB21++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q].diffusion_lamda;
      v4 = layer1[2*p+1][2*q].chemical_lamda;
      layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
      layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CB21 * layer2[p][q][0].D_CB21[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB21--;
      layer1[2*p+1][2*q+1].CB21++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
      v4 = layer1[2*p+1][2*q+1].chemical_lamda;
      layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
      layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

	/*******************************************************/
    
    choose1 = choose;
    choose += layer2[p][q][0].CB12 * layer2[p][q][0].D_CB12[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB12--;
      layer1[2*p][2*q+1].CB12++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p][2*q+1].diffusion_lamda;
      v4 = layer1[2*p][2*q+1].chemical_lamda;
      layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
      layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CB12 * layer2[p][q][0].D_CB12[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB12--;
      layer1[2*p+1][2*q].CB12++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q].diffusion_lamda;
      v4 = layer1[2*p+1][2*q].chemical_lamda;
      layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
      layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CB12 * layer2[p][q][0].D_CB12[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB12--;
      layer1[2*p+1][2*q+1].CB12++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
      v4 = layer1[2*p+1][2*q+1].chemical_lamda;
      layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
      layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

	/*******************************************************/
    
    choose1 = choose;
    choose += layer2[p][q][0].CB22 * layer2[p][q][0].D_CB22[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB22--;
      layer1[2*p][2*q+1].CB22++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p][2*q+1].diffusion_lamda;
      v4 = layer1[2*p][2*q+1].chemical_lamda;
      layer1[2*p][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p, 2*q+1);
      layer1[2*p][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p][2*q+1].diffusion_lamda +
	layer1[2*p][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CB22 * layer2[p][q][0].D_CB22[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB22--;
      layer1[2*p+1][2*q].CB22++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q].diffusion_lamda;
      v4 = layer1[2*p+1][2*q].chemical_lamda;
      layer1[2*p+1][2*q].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q);
      layer1[2*p+1][2*q].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q].diffusion_lamda +
	layer1[2*p+1][2*q].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }

    choose1 = choose;
    choose += layer2[p][q][0].CB22 * layer2[p][q][0].D_CB22[4];
    if((choose > random_num) && (choose1 <= random_num)){
      layer2[p][q][0].CB22--;
      layer1[2*p+1][2*q+1].CB22++;
      v1 = layer2[p][q][0].diffusion_lamda;
      v2 = layer2[p][q][0].chemical_lamda;
      layer2[p][q][0].diffusion_lamda = layer2_diffusion_lamda(p, q, 0);
      layer2[p][q][0].chemical_lamda = voxel_chemical_lamda(layer2[p][q][0]);
      v3 = layer1[2*p+1][2*q+1].diffusion_lamda;
      v4 = layer1[2*p+1][2*q+1].chemical_lamda;
      layer1[2*p+1][2*q+1].diffusion_lamda = layer1_diffusion_lamda(2*p+1, 2*q+1);
      layer1[2*p+1][2*q+1].chemical_lamda = voxel_chemical_lamda(layer1[2*p+1][2*q+1]);
      reaction_lamda += layer2[p][q][0].diffusion_lamda + layer2[p][q][0].chemical_lamda + layer1[2*p+1][2*q+1].diffusion_lamda +
	layer1[2*p+1][2*q+1].chemical_lamda - v1 - v2 - v3 - v4;
      goto loop2;
    }
  }

  x = p;
  y = q;
  z = r;
  chemical_reaction(&layer2[p][q][r], choose, 2);
  
 loop2: k = 0;
}

void reaction_layer3(double choose2, int p){
  int k;
  double choose, choose1, v1, v2, v3, v4;

  choose = choose2;

  choose1 = choose;
  choose += layer3[p].Ca_pump;
  if((choose > random_num) && (choose1 <= random_num)){
    layer3[p].Ca--;
    total_Ca_pumped++;
    total_Ca_pumped3++;
    v1 = layer3[p].diffusion_lamda;
    v2 = layer3[p].chemical_lamda;
    layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
    layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
    reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda - v1 - v2;
  }

  for(k=4; k<6; k++){

    choose1 = choose;
    choose += layer3[p].Ca * layer3[p].D_Ca[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer3[p].Ca--;
      v1 = layer3[p].diffusion_lamda;
      v2 = layer3[p].chemical_lamda;
      layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
      layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
      
      if(k == 4){
	if(p == 0){
	  layer2[2][2][4].Ca++;
	  v3 = layer2[2][2][4].diffusion_lamda;
	  v4 = layer2[2][2][4].chemical_lamda;
	  layer2[2][2][4].diffusion_lamda = layer2_diffusion_lamda(2, 2, 4);
	  layer2[2][2][4].chemical_lamda = voxel_chemical_lamda(layer2[2][2][4]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer2[2][2][4].diffusion_lamda +
	    layer2[2][2][4].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p-1].Ca++;
	  v3 = layer3[p-1].diffusion_lamda;
	  v4 = layer3[p-1].chemical_lamda;
	  layer3[p-1].diffusion_lamda = layer3_diffusion_lamda(p-1);
	  layer3[p-1].chemical_lamda = voxel_chemical_lamda(layer3[p-1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p-1].diffusion_lamda +
	    layer3[p-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	if(p == 7){
	  layer4[2][0][0].Ca++;
	  v3 = layer4[2][0][0].diffusion_lamda;
	  v4 = layer4[2][0][0].chemical_lamda;
	  layer4[2][0][0].diffusion_lamda = layer4_diffusion_lamda(2, 0, 0);
	  layer4[2][0][0].chemical_lamda = voxel_chemical_lamda(layer4[2][0][0]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer4[2][0][0].diffusion_lamda +
	    layer4[2][0][0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p+1].Ca++;
	  v3 = layer3[p+1].diffusion_lamda;
	  v4 = layer3[p+1].chemical_lamda;
	  layer3[p+1].diffusion_lamda = layer3_diffusion_lamda(p+1);
	  layer3[p+1].chemical_lamda = voxel_chemical_lamda(layer3[p+1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p+1].diffusion_lamda +
	    layer3[p+1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }
      goto loop3;
    }
    
    choose1 = choose;
    choose += layer3[p].CaM00 * layer3[p].D_CaM00[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer3[p].CaM00--;
      v1 = layer3[p].diffusion_lamda;
      v2 = layer3[p].chemical_lamda;
      layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
      layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
      
      if(k == 4){
	if(p == 0){
	  layer2[2][2][4].CaM00++;
	  v3 = layer2[2][2][4].diffusion_lamda;
	  v4 = layer2[2][2][4].chemical_lamda;
	  layer2[2][2][4].diffusion_lamda = layer2_diffusion_lamda(2, 2, 4);
	  layer2[2][2][4].chemical_lamda = voxel_chemical_lamda(layer2[2][2][4]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer2[2][2][4].diffusion_lamda +
	    layer2[2][2][4].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer3[p-1].CaM00++;
	  v3 = layer3[p-1].diffusion_lamda;
	  v4 = layer3[p-1].chemical_lamda;
	  layer3[p-1].diffusion_lamda = layer3_diffusion_lamda(p-1);
	  layer3[p-1].chemical_lamda = voxel_chemical_lamda(layer3[p-1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p-1].diffusion_lamda +
	    layer3[p-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	if(p == 7){
	  layer4[2][0][0].CaM00++;
	  v3 = layer4[2][0][0].diffusion_lamda;
	  v4 = layer4[2][0][0].chemical_lamda;
	  layer4[2][0][0].diffusion_lamda = layer4_diffusion_lamda(2, 0, 0);
	  layer4[2][0][0].chemical_lamda = voxel_chemical_lamda(layer4[2][0][0]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer4[2][0][0].diffusion_lamda +
	    layer4[2][0][0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p+1].CaM00++;
	  v3 = layer3[p+1].diffusion_lamda;
	  v4 = layer3[p+1].chemical_lamda;
	  layer3[p+1].diffusion_lamda = layer3_diffusion_lamda(p+1);
	  layer3[p+1].chemical_lamda = voxel_chemical_lamda(layer3[p+1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p+1].diffusion_lamda +
	    layer3[p+1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }
      goto loop3;
    }

	choose1 = choose;
    choose += layer3[p].CaM10 * layer3[p].D_CaM10[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer3[p].CaM10--;
      v1 = layer3[p].diffusion_lamda;
      v2 = layer3[p].chemical_lamda;
      layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
      layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
      
      if(k == 4){
	if(p == 0){
	  layer2[2][2][4].CaM10++;
	  v3 = layer2[2][2][4].diffusion_lamda;
	  v4 = layer2[2][2][4].chemical_lamda;
	  layer2[2][2][4].diffusion_lamda = layer2_diffusion_lamda(2, 2, 4);
	  layer2[2][2][4].chemical_lamda = voxel_chemical_lamda(layer2[2][2][4]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer2[2][2][4].diffusion_lamda +
	    layer2[2][2][4].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer3[p-1].CaM10++;
	  v3 = layer3[p-1].diffusion_lamda;
	  v4 = layer3[p-1].chemical_lamda;
	  layer3[p-1].diffusion_lamda = layer3_diffusion_lamda(p-1);
	  layer3[p-1].chemical_lamda = voxel_chemical_lamda(layer3[p-1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p-1].diffusion_lamda +
	    layer3[p-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	if(p == 7){
	  layer4[2][0][0].CaM10++;
	  v3 = layer4[2][0][0].diffusion_lamda;
	  v4 = layer4[2][0][0].chemical_lamda;
	  layer4[2][0][0].diffusion_lamda = layer4_diffusion_lamda(2, 0, 0);
	  layer4[2][0][0].chemical_lamda = voxel_chemical_lamda(layer4[2][0][0]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer4[2][0][0].diffusion_lamda +
	    layer4[2][0][0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p+1].CaM10++;
	  v3 = layer3[p+1].diffusion_lamda;
	  v4 = layer3[p+1].chemical_lamda;
	  layer3[p+1].diffusion_lamda = layer3_diffusion_lamda(p+1);
	  layer3[p+1].chemical_lamda = voxel_chemical_lamda(layer3[p+1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p+1].diffusion_lamda +
	    layer3[p+1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }
      goto loop3;
    }

	choose1 = choose;
    choose += layer3[p].CaM01 * layer3[p].D_CaM01[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer3[p].CaM01--;
      v1 = layer3[p].diffusion_lamda;
      v2 = layer3[p].chemical_lamda;
      layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
      layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
      
      if(k == 4){
	if(p == 0){
	  layer2[2][2][4].CaM01++;
	  v3 = layer2[2][2][4].diffusion_lamda;
	  v4 = layer2[2][2][4].chemical_lamda;
	  layer2[2][2][4].diffusion_lamda = layer2_diffusion_lamda(2, 2, 4);
	  layer2[2][2][4].chemical_lamda = voxel_chemical_lamda(layer2[2][2][4]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer2[2][2][4].diffusion_lamda +
	    layer2[2][2][4].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer3[p-1].CaM01++;
	  v3 = layer3[p-1].diffusion_lamda;
	  v4 = layer3[p-1].chemical_lamda;
	  layer3[p-1].diffusion_lamda = layer3_diffusion_lamda(p-1);
	  layer3[p-1].chemical_lamda = voxel_chemical_lamda(layer3[p-1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p-1].diffusion_lamda +
	    layer3[p-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	if(p == 7){
	  layer4[2][0][0].CaM01++;
	  v3 = layer4[2][0][0].diffusion_lamda;
	  v4 = layer4[2][0][0].chemical_lamda;
	  layer4[2][0][0].diffusion_lamda = layer4_diffusion_lamda(2, 0, 0);
	  layer4[2][0][0].chemical_lamda = voxel_chemical_lamda(layer4[2][0][0]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer4[2][0][0].diffusion_lamda +
	    layer4[2][0][0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p+1].CaM01++;
	  v3 = layer3[p+1].diffusion_lamda;
	  v4 = layer3[p+1].chemical_lamda;
	  layer3[p+1].diffusion_lamda = layer3_diffusion_lamda(p+1);
	  layer3[p+1].chemical_lamda = voxel_chemical_lamda(layer3[p+1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p+1].diffusion_lamda +
	    layer3[p+1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }
      goto loop3;
    }

	choose1 = choose;
    choose += layer3[p].CaM11 * layer3[p].D_CaM11[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer3[p].CaM11--;
      v1 = layer3[p].diffusion_lamda;
      v2 = layer3[p].chemical_lamda;
      layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
      layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
      
      if(k == 4){
	if(p == 0){
	  layer2[2][2][4].CaM11++;
	  v3 = layer2[2][2][4].diffusion_lamda;
	  v4 = layer2[2][2][4].chemical_lamda;
	  layer2[2][2][4].diffusion_lamda = layer2_diffusion_lamda(2, 2, 4);
	  layer2[2][2][4].chemical_lamda = voxel_chemical_lamda(layer2[2][2][4]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer2[2][2][4].diffusion_lamda +
	    layer2[2][2][4].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer3[p-1].CaM11++;
	  v3 = layer3[p-1].diffusion_lamda;
	  v4 = layer3[p-1].chemical_lamda;
	  layer3[p-1].diffusion_lamda = layer3_diffusion_lamda(p-1);
	  layer3[p-1].chemical_lamda = voxel_chemical_lamda(layer3[p-1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p-1].diffusion_lamda +
	    layer3[p-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	if(p == 7){
	  layer4[2][0][0].CaM11++;
	  v3 = layer4[2][0][0].diffusion_lamda;
	  v4 = layer4[2][0][0].chemical_lamda;
	  layer4[2][0][0].diffusion_lamda = layer4_diffusion_lamda(2, 0, 0);
	  layer4[2][0][0].chemical_lamda = voxel_chemical_lamda(layer4[2][0][0]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer4[2][0][0].diffusion_lamda +
	    layer4[2][0][0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p+1].CaM11++;
	  v3 = layer3[p+1].diffusion_lamda;
	  v4 = layer3[p+1].chemical_lamda;
	  layer3[p+1].diffusion_lamda = layer3_diffusion_lamda(p+1);
	  layer3[p+1].chemical_lamda = voxel_chemical_lamda(layer3[p+1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p+1].diffusion_lamda +
	    layer3[p+1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }
      goto loop3;
    }

	choose1 = choose;
    choose += layer3[p].CaM20 * layer3[p].D_CaM20[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer3[p].CaM20--;
      v1 = layer3[p].diffusion_lamda;
      v2 = layer3[p].chemical_lamda;
      layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
      layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
      
      if(k == 4){
	if(p == 0){
	  layer2[2][2][4].CaM20++;
	  v3 = layer2[2][2][4].diffusion_lamda;
	  v4 = layer2[2][2][4].chemical_lamda;
	  layer2[2][2][4].diffusion_lamda = layer2_diffusion_lamda(2, 2, 4);
	  layer2[2][2][4].chemical_lamda = voxel_chemical_lamda(layer2[2][2][4]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer2[2][2][4].diffusion_lamda +
	    layer2[2][2][4].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer3[p-1].CaM20++;
	  v3 = layer3[p-1].diffusion_lamda;
	  v4 = layer3[p-1].chemical_lamda;
	  layer3[p-1].diffusion_lamda = layer3_diffusion_lamda(p-1);
	  layer3[p-1].chemical_lamda = voxel_chemical_lamda(layer3[p-1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p-1].diffusion_lamda +
	    layer3[p-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	if(p == 7){
	  layer4[2][0][0].CaM20++;
	  v3 = layer4[2][0][0].diffusion_lamda;
	  v4 = layer4[2][0][0].chemical_lamda;
	  layer4[2][0][0].diffusion_lamda = layer4_diffusion_lamda(2, 0, 0);
	  layer4[2][0][0].chemical_lamda = voxel_chemical_lamda(layer4[2][0][0]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer4[2][0][0].diffusion_lamda +
	    layer4[2][0][0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p+1].CaM20++;
	  v3 = layer3[p+1].diffusion_lamda;
	  v4 = layer3[p+1].chemical_lamda;
	  layer3[p+1].diffusion_lamda = layer3_diffusion_lamda(p+1);
	  layer3[p+1].chemical_lamda = voxel_chemical_lamda(layer3[p+1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p+1].diffusion_lamda +
	    layer3[p+1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }
      goto loop3;
    }

	choose1 = choose;
    choose += layer3[p].CaM02 * layer3[p].D_CaM02[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer3[p].CaM02--;
      v1 = layer3[p].diffusion_lamda;
      v2 = layer3[p].chemical_lamda;
      layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
      layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
      
      if(k == 4){
	if(p == 0){
	  layer2[2][2][4].CaM02++;
	  v3 = layer2[2][2][4].diffusion_lamda;
	  v4 = layer2[2][2][4].chemical_lamda;
	  layer2[2][2][4].diffusion_lamda = layer2_diffusion_lamda(2, 2, 4);
	  layer2[2][2][4].chemical_lamda = voxel_chemical_lamda(layer2[2][2][4]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer2[2][2][4].diffusion_lamda +
	    layer2[2][2][4].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer3[p-1].CaM02++;
	  v3 = layer3[p-1].diffusion_lamda;
	  v4 = layer3[p-1].chemical_lamda;
	  layer3[p-1].diffusion_lamda = layer3_diffusion_lamda(p-1);
	  layer3[p-1].chemical_lamda = voxel_chemical_lamda(layer3[p-1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p-1].diffusion_lamda +
	    layer3[p-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	if(p == 7){
	  layer4[2][0][0].CaM02++;
	  v3 = layer4[2][0][0].diffusion_lamda;
	  v4 = layer4[2][0][0].chemical_lamda;
	  layer4[2][0][0].diffusion_lamda = layer4_diffusion_lamda(2, 0, 0);
	  layer4[2][0][0].chemical_lamda = voxel_chemical_lamda(layer4[2][0][0]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer4[2][0][0].diffusion_lamda +
	    layer4[2][0][0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p+1].CaM02++;
	  v3 = layer3[p+1].diffusion_lamda;
	  v4 = layer3[p+1].chemical_lamda;
	  layer3[p+1].diffusion_lamda = layer3_diffusion_lamda(p+1);
	  layer3[p+1].chemical_lamda = voxel_chemical_lamda(layer3[p+1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p+1].diffusion_lamda +
	    layer3[p+1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }
      goto loop3;
    }

	choose1 = choose;
    choose += layer3[p].CaM21 * layer3[p].D_CaM21[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer3[p].CaM21--;
      v1 = layer3[p].diffusion_lamda;
      v2 = layer3[p].chemical_lamda;
      layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
      layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
      
      if(k == 4){
	if(p == 0){
	  layer2[2][2][4].CaM21++;
	  v3 = layer2[2][2][4].diffusion_lamda;
	  v4 = layer2[2][2][4].chemical_lamda;
	  layer2[2][2][4].diffusion_lamda = layer2_diffusion_lamda(2, 2, 4);
	  layer2[2][2][4].chemical_lamda = voxel_chemical_lamda(layer2[2][2][4]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer2[2][2][4].diffusion_lamda +
	    layer2[2][2][4].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer3[p-1].CaM21++;
	  v3 = layer3[p-1].diffusion_lamda;
	  v4 = layer3[p-1].chemical_lamda;
	  layer3[p-1].diffusion_lamda = layer3_diffusion_lamda(p-1);
	  layer3[p-1].chemical_lamda = voxel_chemical_lamda(layer3[p-1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p-1].diffusion_lamda +
	    layer3[p-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	if(p == 7){
	  layer4[2][0][0].CaM21++;
	  v3 = layer4[2][0][0].diffusion_lamda;
	  v4 = layer4[2][0][0].chemical_lamda;
	  layer4[2][0][0].diffusion_lamda = layer4_diffusion_lamda(2, 0, 0);
	  layer4[2][0][0].chemical_lamda = voxel_chemical_lamda(layer4[2][0][0]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer4[2][0][0].diffusion_lamda +
	    layer4[2][0][0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p+1].CaM21++;
	  v3 = layer3[p+1].diffusion_lamda;
	  v4 = layer3[p+1].chemical_lamda;
	  layer3[p+1].diffusion_lamda = layer3_diffusion_lamda(p+1);
	  layer3[p+1].chemical_lamda = voxel_chemical_lamda(layer3[p+1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p+1].diffusion_lamda +
	    layer3[p+1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }
      goto loop3;
    }

	choose1 = choose;
    choose += layer3[p].CaM12 * layer3[p].D_CaM12[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer3[p].CaM12--;
      v1 = layer3[p].diffusion_lamda;
      v2 = layer3[p].chemical_lamda;
      layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
      layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
      
      if(k == 4){
	if(p == 0){
	  layer2[2][2][4].CaM12++;
	  v3 = layer2[2][2][4].diffusion_lamda;
	  v4 = layer2[2][2][4].chemical_lamda;
	  layer2[2][2][4].diffusion_lamda = layer2_diffusion_lamda(2, 2, 4);
	  layer2[2][2][4].chemical_lamda = voxel_chemical_lamda(layer2[2][2][4]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer2[2][2][4].diffusion_lamda +
	    layer2[2][2][4].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer3[p-1].CaM12++;
	  v3 = layer3[p-1].diffusion_lamda;
	  v4 = layer3[p-1].chemical_lamda;
	  layer3[p-1].diffusion_lamda = layer3_diffusion_lamda(p-1);
	  layer3[p-1].chemical_lamda = voxel_chemical_lamda(layer3[p-1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p-1].diffusion_lamda +
	    layer3[p-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	if(p == 7){
	  layer4[2][0][0].CaM12++;
	  v3 = layer4[2][0][0].diffusion_lamda;
	  v4 = layer4[2][0][0].chemical_lamda;
	  layer4[2][0][0].diffusion_lamda = layer4_diffusion_lamda(2, 0, 0);
	  layer4[2][0][0].chemical_lamda = voxel_chemical_lamda(layer4[2][0][0]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer4[2][0][0].diffusion_lamda +
	    layer4[2][0][0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p+1].CaM12++;
	  v3 = layer3[p+1].diffusion_lamda;
	  v4 = layer3[p+1].chemical_lamda;
	  layer3[p+1].diffusion_lamda = layer3_diffusion_lamda(p+1);
	  layer3[p+1].chemical_lamda = voxel_chemical_lamda(layer3[p+1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p+1].diffusion_lamda +
	    layer3[p+1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }
      goto loop3;
    }

	choose1 = choose;
    choose += layer3[p].CaM22 * layer3[p].D_CaM22[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer3[p].CaM22--;
      v1 = layer3[p].diffusion_lamda;
      v2 = layer3[p].chemical_lamda;
      layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
      layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
      
      if(k == 4){
	if(p == 0){
	  layer2[2][2][4].CaM22++;
	  v3 = layer2[2][2][4].diffusion_lamda;
	  v4 = layer2[2][2][4].chemical_lamda;
	  layer2[2][2][4].diffusion_lamda = layer2_diffusion_lamda(2, 2, 4);
	  layer2[2][2][4].chemical_lamda = voxel_chemical_lamda(layer2[2][2][4]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer2[2][2][4].diffusion_lamda +
	    layer2[2][2][4].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer3[p-1].CaM22++;
	  v3 = layer3[p-1].diffusion_lamda;
	  v4 = layer3[p-1].chemical_lamda;
	  layer3[p-1].diffusion_lamda = layer3_diffusion_lamda(p-1);
	  layer3[p-1].chemical_lamda = voxel_chemical_lamda(layer3[p-1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p-1].diffusion_lamda +
	    layer3[p-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	if(p == 7){
	  layer4[2][0][0].CaM22++;
	  v3 = layer4[2][0][0].diffusion_lamda;
	  v4 = layer4[2][0][0].chemical_lamda;
	  layer4[2][0][0].diffusion_lamda = layer4_diffusion_lamda(2, 0, 0);
	  layer4[2][0][0].chemical_lamda = voxel_chemical_lamda(layer4[2][0][0]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer4[2][0][0].diffusion_lamda +
	    layer4[2][0][0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p+1].CaM22++;
	  v3 = layer3[p+1].diffusion_lamda;
	  v4 = layer3[p+1].chemical_lamda;
	  layer3[p+1].diffusion_lamda = layer3_diffusion_lamda(p+1);
	  layer3[p+1].chemical_lamda = voxel_chemical_lamda(layer3[p+1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p+1].diffusion_lamda +
	    layer3[p+1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }
      goto loop3;
    }

	choose1 = choose;
    choose += layer3[p].CB00 * layer3[p].D_CB00[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer3[p].CB00--;
      v1 = layer3[p].diffusion_lamda;
      v2 = layer3[p].chemical_lamda;
      layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
      layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
      if(k == 4){
	if(p == 0){
	  layer2[2][2][4].CB00++;
	  v3 = layer2[2][2][4].diffusion_lamda;
	  v4 = layer2[2][2][4].chemical_lamda;
	  layer2[2][2][4].diffusion_lamda = layer2_diffusion_lamda(2, 2, 4);
	  layer2[2][2][4].chemical_lamda = voxel_chemical_lamda(layer2[2][2][4]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer2[2][2][4].diffusion_lamda +
	    layer2[2][2][4].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p-1].CB00++;
	  v3 = layer3[p-1].diffusion_lamda;
	  v4 = layer3[p-1].chemical_lamda;
	  layer3[p-1].diffusion_lamda = layer3_diffusion_lamda(p-1);
	  layer3[p-1].chemical_lamda = voxel_chemical_lamda(layer3[p-1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p-1].diffusion_lamda +
	    layer3[p-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	if(p == 7){
	  layer4[2][0][0].CB00++;
	  v3 = layer4[2][0][0].diffusion_lamda;
	  v4 = layer4[2][0][0].chemical_lamda;
	  layer4[2][0][0].diffusion_lamda = layer4_diffusion_lamda(2, 0, 0);
	  layer4[2][0][0].chemical_lamda = voxel_chemical_lamda(layer4[2][0][0]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer4[2][0][0].diffusion_lamda +
	    layer4[2][0][0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p+1].CB00++;
	  v3 = layer3[p+1].diffusion_lamda;
	  v4 = layer3[p+1].chemical_lamda;
	  layer3[p+1].diffusion_lamda = layer3_diffusion_lamda(p+1);
	  layer3[p+1].chemical_lamda = voxel_chemical_lamda(layer3[p+1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p+1].diffusion_lamda +
	    layer3[p+1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }
      goto loop3; 
    }

	choose1 = choose;
    choose += layer3[p].CB10 * layer3[p].D_CB10[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer3[p].CB10--;
      v1 = layer3[p].diffusion_lamda;
      v2 = layer3[p].chemical_lamda;
      layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
      layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
      if(k == 4){
	if(p == 0){
	  layer2[2][2][4].CB10++;
	  v3 = layer2[2][2][4].diffusion_lamda;
	  v4 = layer2[2][2][4].chemical_lamda;
	  layer2[2][2][4].diffusion_lamda = layer2_diffusion_lamda(2, 2, 4);
	  layer2[2][2][4].chemical_lamda = voxel_chemical_lamda(layer2[2][2][4]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer2[2][2][4].diffusion_lamda +
	    layer2[2][2][4].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p-1].CB10++;
	  v3 = layer3[p-1].diffusion_lamda;
	  v4 = layer3[p-1].chemical_lamda;
	  layer3[p-1].diffusion_lamda = layer3_diffusion_lamda(p-1);
	  layer3[p-1].chemical_lamda = voxel_chemical_lamda(layer3[p-1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p-1].diffusion_lamda +
	    layer3[p-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	if(p == 7){
	  layer4[2][0][0].CB10++;
	  v3 = layer4[2][0][0].diffusion_lamda;
	  v4 = layer4[2][0][0].chemical_lamda;
	  layer4[2][0][0].diffusion_lamda = layer4_diffusion_lamda(2, 0, 0);
	  layer4[2][0][0].chemical_lamda = voxel_chemical_lamda(layer4[2][0][0]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer4[2][0][0].diffusion_lamda +
	    layer4[2][0][0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p+1].CB10++;
	  v3 = layer3[p+1].diffusion_lamda;
	  v4 = layer3[p+1].chemical_lamda;
	  layer3[p+1].diffusion_lamda = layer3_diffusion_lamda(p+1);
	  layer3[p+1].chemical_lamda = voxel_chemical_lamda(layer3[p+1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p+1].diffusion_lamda +
	    layer3[p+1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }
      goto loop3; 
    }

	choose1 = choose;
    choose += layer3[p].CB01 * layer3[p].D_CB01[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer3[p].CB01--;
      v1 = layer3[p].diffusion_lamda;
      v2 = layer3[p].chemical_lamda;
      layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
      layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
      if(k == 4){
	if(p == 0){
	  layer2[2][2][4].CB01++;
	  v3 = layer2[2][2][4].diffusion_lamda;
	  v4 = layer2[2][2][4].chemical_lamda;
	  layer2[2][2][4].diffusion_lamda = layer2_diffusion_lamda(2, 2, 4);
	  layer2[2][2][4].chemical_lamda = voxel_chemical_lamda(layer2[2][2][4]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer2[2][2][4].diffusion_lamda +
	    layer2[2][2][4].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p-1].CB01++;
	  v3 = layer3[p-1].diffusion_lamda;
	  v4 = layer3[p-1].chemical_lamda;
	  layer3[p-1].diffusion_lamda = layer3_diffusion_lamda(p-1);
	  layer3[p-1].chemical_lamda = voxel_chemical_lamda(layer3[p-1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p-1].diffusion_lamda +
	    layer3[p-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	if(p == 7){
	  layer4[2][0][0].CB01++;
	  v3 = layer4[2][0][0].diffusion_lamda;
	  v4 = layer4[2][0][0].chemical_lamda;
	  layer4[2][0][0].diffusion_lamda = layer4_diffusion_lamda(2, 0, 0);
	  layer4[2][0][0].chemical_lamda = voxel_chemical_lamda(layer4[2][0][0]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer4[2][0][0].diffusion_lamda +
	    layer4[2][0][0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p+1].CB01++;
	  v3 = layer3[p+1].diffusion_lamda;
	  v4 = layer3[p+1].chemical_lamda;
	  layer3[p+1].diffusion_lamda = layer3_diffusion_lamda(p+1);
	  layer3[p+1].chemical_lamda = voxel_chemical_lamda(layer3[p+1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p+1].diffusion_lamda +
	    layer3[p+1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }
      goto loop3; 
    }

	choose1 = choose;
    choose += layer3[p].CB11 * layer3[p].D_CB11[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer3[p].CB11--;
      v1 = layer3[p].diffusion_lamda;
      v2 = layer3[p].chemical_lamda;
      layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
      layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
      if(k == 4){
	if(p == 0){
	  layer2[2][2][4].CB11++;
	  v3 = layer2[2][2][4].diffusion_lamda;
	  v4 = layer2[2][2][4].chemical_lamda;
	  layer2[2][2][4].diffusion_lamda = layer2_diffusion_lamda(2, 2, 4);
	  layer2[2][2][4].chemical_lamda = voxel_chemical_lamda(layer2[2][2][4]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer2[2][2][4].diffusion_lamda +
	    layer2[2][2][4].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p-1].CB11++;
	  v3 = layer3[p-1].diffusion_lamda;
	  v4 = layer3[p-1].chemical_lamda;
	  layer3[p-1].diffusion_lamda = layer3_diffusion_lamda(p-1);
	  layer3[p-1].chemical_lamda = voxel_chemical_lamda(layer3[p-1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p-1].diffusion_lamda +
	    layer3[p-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	if(p == 7){
	  layer4[2][0][0].CB11++;
	  v3 = layer4[2][0][0].diffusion_lamda;
	  v4 = layer4[2][0][0].chemical_lamda;
	  layer4[2][0][0].diffusion_lamda = layer4_diffusion_lamda(2, 0, 0);
	  layer4[2][0][0].chemical_lamda = voxel_chemical_lamda(layer4[2][0][0]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer4[2][0][0].diffusion_lamda +
	    layer4[2][0][0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p+1].CB11++;
	  v3 = layer3[p+1].diffusion_lamda;
	  v4 = layer3[p+1].chemical_lamda;
	  layer3[p+1].diffusion_lamda = layer3_diffusion_lamda(p+1);
	  layer3[p+1].chemical_lamda = voxel_chemical_lamda(layer3[p+1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p+1].diffusion_lamda +
	    layer3[p+1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }
      goto loop3; 
    }

	choose1 = choose;
    choose += layer3[p].CB20 * layer3[p].D_CB20[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer3[p].CB20--;
      v1 = layer3[p].diffusion_lamda;
      v2 = layer3[p].chemical_lamda;
      layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
      layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
      if(k == 4){
	if(p == 0){
	  layer2[2][2][4].CB20++;
	  v3 = layer2[2][2][4].diffusion_lamda;
	  v4 = layer2[2][2][4].chemical_lamda;
	  layer2[2][2][4].diffusion_lamda = layer2_diffusion_lamda(2, 2, 4);
	  layer2[2][2][4].chemical_lamda = voxel_chemical_lamda(layer2[2][2][4]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer2[2][2][4].diffusion_lamda +
	    layer2[2][2][4].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p-1].CB20++;
	  v3 = layer3[p-1].diffusion_lamda;
	  v4 = layer3[p-1].chemical_lamda;
	  layer3[p-1].diffusion_lamda = layer3_diffusion_lamda(p-1);
	  layer3[p-1].chemical_lamda = voxel_chemical_lamda(layer3[p-1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p-1].diffusion_lamda +
	    layer3[p-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	if(p == 7){
	  layer4[2][0][0].CB20++;
	  v3 = layer4[2][0][0].diffusion_lamda;
	  v4 = layer4[2][0][0].chemical_lamda;
	  layer4[2][0][0].diffusion_lamda = layer4_diffusion_lamda(2, 0, 0);
	  layer4[2][0][0].chemical_lamda = voxel_chemical_lamda(layer4[2][0][0]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer4[2][0][0].diffusion_lamda +
	    layer4[2][0][0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p+1].CB20++;
	  v3 = layer3[p+1].diffusion_lamda;
	  v4 = layer3[p+1].chemical_lamda;
	  layer3[p+1].diffusion_lamda = layer3_diffusion_lamda(p+1);
	  layer3[p+1].chemical_lamda = voxel_chemical_lamda(layer3[p+1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p+1].diffusion_lamda +
	    layer3[p+1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }
      goto loop3; 
    }

	choose1 = choose;
    choose += layer3[p].CB02 * layer3[p].D_CB02[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer3[p].CB02--;
      v1 = layer3[p].diffusion_lamda;
      v2 = layer3[p].chemical_lamda;
      layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
      layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
      if(k == 4){
	if(p == 0){
	  layer2[2][2][4].CB02++;
	  v3 = layer2[2][2][4].diffusion_lamda;
	  v4 = layer2[2][2][4].chemical_lamda;
	  layer2[2][2][4].diffusion_lamda = layer2_diffusion_lamda(2, 2, 4);
	  layer2[2][2][4].chemical_lamda = voxel_chemical_lamda(layer2[2][2][4]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer2[2][2][4].diffusion_lamda +
	    layer2[2][2][4].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p-1].CB02++;
	  v3 = layer3[p-1].diffusion_lamda;
	  v4 = layer3[p-1].chemical_lamda;
	  layer3[p-1].diffusion_lamda = layer3_diffusion_lamda(p-1);
	  layer3[p-1].chemical_lamda = voxel_chemical_lamda(layer3[p-1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p-1].diffusion_lamda +
	    layer3[p-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	if(p == 7){
	  layer4[2][0][0].CB02++;
	  v3 = layer4[2][0][0].diffusion_lamda;
	  v4 = layer4[2][0][0].chemical_lamda;
	  layer4[2][0][0].diffusion_lamda = layer4_diffusion_lamda(2, 0, 0);
	  layer4[2][0][0].chemical_lamda = voxel_chemical_lamda(layer4[2][0][0]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer4[2][0][0].diffusion_lamda +
	    layer4[2][0][0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p+1].CB02++;
	  v3 = layer3[p+1].diffusion_lamda;
	  v4 = layer3[p+1].chemical_lamda;
	  layer3[p+1].diffusion_lamda = layer3_diffusion_lamda(p+1);
	  layer3[p+1].chemical_lamda = voxel_chemical_lamda(layer3[p+1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p+1].diffusion_lamda +
	    layer3[p+1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }
      goto loop3; 
    }

	choose1 = choose;
    choose += layer3[p].CB21 * layer3[p].D_CB21[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer3[p].CB21--;
      v1 = layer3[p].diffusion_lamda;
      v2 = layer3[p].chemical_lamda;
      layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
      layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
      if(k == 4){
	if(p == 0){
	  layer2[2][2][4].CB21++;
	  v3 = layer2[2][2][4].diffusion_lamda;
	  v4 = layer2[2][2][4].chemical_lamda;
	  layer2[2][2][4].diffusion_lamda = layer2_diffusion_lamda(2, 2, 4);
	  layer2[2][2][4].chemical_lamda = voxel_chemical_lamda(layer2[2][2][4]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer2[2][2][4].diffusion_lamda +
	    layer2[2][2][4].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p-1].CB21++;
	  v3 = layer3[p-1].diffusion_lamda;
	  v4 = layer3[p-1].chemical_lamda;
	  layer3[p-1].diffusion_lamda = layer3_diffusion_lamda(p-1);
	  layer3[p-1].chemical_lamda = voxel_chemical_lamda(layer3[p-1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p-1].diffusion_lamda +
	    layer3[p-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	if(p == 7){
	  layer4[2][0][0].CB21++;
	  v3 = layer4[2][0][0].diffusion_lamda;
	  v4 = layer4[2][0][0].chemical_lamda;
	  layer4[2][0][0].diffusion_lamda = layer4_diffusion_lamda(2, 0, 0);
	  layer4[2][0][0].chemical_lamda = voxel_chemical_lamda(layer4[2][0][0]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer4[2][0][0].diffusion_lamda +
	    layer4[2][0][0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p+1].CB21++;
	  v3 = layer3[p+1].diffusion_lamda;
	  v4 = layer3[p+1].chemical_lamda;
	  layer3[p+1].diffusion_lamda = layer3_diffusion_lamda(p+1);
	  layer3[p+1].chemical_lamda = voxel_chemical_lamda(layer3[p+1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p+1].diffusion_lamda +
	    layer3[p+1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }
      goto loop3; 
    }

	choose1 = choose;
    choose += layer3[p].CB12 * layer3[p].D_CB12[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer3[p].CB12--;
      v1 = layer3[p].diffusion_lamda;
      v2 = layer3[p].chemical_lamda;
      layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
      layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
      if(k == 4){
	if(p == 0){
	  layer2[2][2][4].CB12++;
	  v3 = layer2[2][2][4].diffusion_lamda;
	  v4 = layer2[2][2][4].chemical_lamda;
	  layer2[2][2][4].diffusion_lamda = layer2_diffusion_lamda(2, 2, 4);
	  layer2[2][2][4].chemical_lamda = voxel_chemical_lamda(layer2[2][2][4]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer2[2][2][4].diffusion_lamda +
	    layer2[2][2][4].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p-1].CB12++;
	  v3 = layer3[p-1].diffusion_lamda;
	  v4 = layer3[p-1].chemical_lamda;
	  layer3[p-1].diffusion_lamda = layer3_diffusion_lamda(p-1);
	  layer3[p-1].chemical_lamda = voxel_chemical_lamda(layer3[p-1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p-1].diffusion_lamda +
	    layer3[p-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	if(p == 7){
	  layer4[2][0][0].CB12++;
	  v3 = layer4[2][0][0].diffusion_lamda;
	  v4 = layer4[2][0][0].chemical_lamda;
	  layer4[2][0][0].diffusion_lamda = layer4_diffusion_lamda(2, 0, 0);
	  layer4[2][0][0].chemical_lamda = voxel_chemical_lamda(layer4[2][0][0]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer4[2][0][0].diffusion_lamda +
	    layer4[2][0][0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p+1].CB12++;
	  v3 = layer3[p+1].diffusion_lamda;
	  v4 = layer3[p+1].chemical_lamda;
	  layer3[p+1].diffusion_lamda = layer3_diffusion_lamda(p+1);
	  layer3[p+1].chemical_lamda = voxel_chemical_lamda(layer3[p+1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p+1].diffusion_lamda +
	    layer3[p+1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }
      goto loop3; 
    }

	choose1 = choose;
    choose += layer3[p].CB22 * layer3[p].D_CB22[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer3[p].CB22--;
      v1 = layer3[p].diffusion_lamda;
      v2 = layer3[p].chemical_lamda;
      layer3[p].diffusion_lamda = layer3_diffusion_lamda(p);
      layer3[p].chemical_lamda = voxel_chemical_lamda(layer3[p]);
      if(k == 4){
	if(p == 0){
	  layer2[2][2][4].CB22++;
	  v3 = layer2[2][2][4].diffusion_lamda;
	  v4 = layer2[2][2][4].chemical_lamda;
	  layer2[2][2][4].diffusion_lamda = layer2_diffusion_lamda(2, 2, 4);
	  layer2[2][2][4].chemical_lamda = voxel_chemical_lamda(layer2[2][2][4]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer2[2][2][4].diffusion_lamda +
	    layer2[2][2][4].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p-1].CB22++;
	  v3 = layer3[p-1].diffusion_lamda;
	  v4 = layer3[p-1].chemical_lamda;
	  layer3[p-1].diffusion_lamda = layer3_diffusion_lamda(p-1);
	  layer3[p-1].chemical_lamda = voxel_chemical_lamda(layer3[p-1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p-1].diffusion_lamda +
	    layer3[p-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	if(p == 7){
	  layer4[2][0][0].CB22++;
	  v3 = layer4[2][0][0].diffusion_lamda;
	  v4 = layer4[2][0][0].chemical_lamda;
	  layer4[2][0][0].diffusion_lamda = layer4_diffusion_lamda(2, 0, 0);
	  layer4[2][0][0].chemical_lamda = voxel_chemical_lamda(layer4[2][0][0]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer4[2][0][0].diffusion_lamda +
	    layer4[2][0][0].chemical_lamda - v1 - v2 - v3 - v4;
	}else{
	  layer3[p+1].CB22++;
	  v3 = layer3[p+1].diffusion_lamda;
	  v4 = layer3[p+1].chemical_lamda;
	  layer3[p+1].diffusion_lamda = layer3_diffusion_lamda(p+1);
	  layer3[p+1].chemical_lamda = voxel_chemical_lamda(layer3[p+1]);
	  reaction_lamda += layer3[p].diffusion_lamda + layer3[p].chemical_lamda + layer3[p+1].diffusion_lamda +
	    layer3[p+1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }
      goto loop3; 
    }
  }

  x = p;
  chemical_reaction(&layer3[p], choose, 3);
  
 loop3: k = 0;
}

void reaction_layer4(double choose2, int p, int q, int r){
  int k;
  double choose, choose1, v1, v2, v3, v4;
  
  choose = choose2;

  choose1 = choose;
  choose += layer4[p][q][r].Ca_pump;
  if((choose > random_num) && (choose1 <= random_num)){
    layer4[p][q][r].Ca--;
    total_Ca_pumped++;
    total_Ca_pumped4++;
    v1 = layer4[p][q][r].diffusion_lamda;
    v2 = layer4[p][q][r].chemical_lamda;
    layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
    layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);
    reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda - v1 - v2;
  }

  for(k=0; k<6; k++){
    
    choose1 = choose;
    choose += layer4[p][q][r].Ca * layer4[p][q][r].D_Ca[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer4[p][q][r].Ca--;
      v1 = layer4[p][q][r].diffusion_lamda;
      v2 = layer4[p][q][r].chemical_lamda;
      layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
      layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);

      if(k == 0){
	layer4[p-1][q][r].Ca++;
	v3 = layer4[p-1][q][r].diffusion_lamda;
	v4 = layer4[p-1][q][r].chemical_lamda;
	layer4[p-1][q][r].diffusion_lamda = layer4_diffusion_lamda(p-1, q, r);
	layer4[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p-1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p-1][q][r].diffusion_lamda +
	  layer4[p-1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 1){
	layer4[p][q-1][r].Ca++;
	v3 = layer4[p][q-1][r].diffusion_lamda;
	v4 = layer4[p][q-1][r].chemical_lamda;
	layer4[p][q-1][r].diffusion_lamda = layer4_diffusion_lamda(p, q-1, r);
	layer4[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q-1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q-1][r].diffusion_lamda +
	  layer4[p][q-1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 2){
	layer4[p+1][q][r].Ca++;
	v3 = layer4[p+1][q][r].diffusion_lamda;
	v4 = layer4[p+1][q][r].chemical_lamda;
	layer4[p+1][q][r].diffusion_lamda = layer4_diffusion_lamda(p+1, q, r);
	layer4[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p+1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p+1][q][r].diffusion_lamda +
	  layer4[p+1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 3){
	layer4[p][q+1][r].Ca++;
	v3 = layer4[p][q+1][r].diffusion_lamda;
	v4 = layer4[p][q+1][r].chemical_lamda;
	layer4[p][q+1][r].diffusion_lamda = layer4_diffusion_lamda(p, q+1, r);
	layer4[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q+1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q+1][r].diffusion_lamda +
	  layer4[p][q+1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 4){
	if((p == 2) && (q == 0) && (r == 0)){
	  layer3[7].Ca++;
	  v3 = layer3[7].diffusion_lamda;
	  v4 = layer3[7].chemical_lamda;
	  layer3[7].diffusion_lamda = layer3_diffusion_lamda(7);
	  layer3[7].chemical_lamda = voxel_chemical_lamda(layer3[7]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer3[7].diffusion_lamda +
	    layer3[7].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer4[p][q][r-1].Ca++;
	  v3 = layer4[p][q][r-1].diffusion_lamda;
	  v4 = layer4[p][q][r-1].chemical_lamda;
	  layer4[p][q][r-1].diffusion_lamda = layer4_diffusion_lamda(p, q, r-1);
	  layer4[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r-1]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r-1].diffusion_lamda +
	    layer4[p][q][r-1].chemical_lamda - v1 - v2 - v3 - v4;

	}
      }else{
	layer4[p][q][r+1].Ca++;
	v3 = layer4[p][q][r+1].diffusion_lamda;
	v4 = layer4[p][q][r+1].chemical_lamda;
	layer4[p][q][r+1].diffusion_lamda = layer4_diffusion_lamda(p, q, r+1);
	layer4[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r+1]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r+1].diffusion_lamda +
	  layer4[p][q][r+1].chemical_lamda - v1 - v2 - v3 - v4;
      }
      goto loop4;
    }

    choose1 = choose;
    choose += layer4[p][q][r].CaM00 * layer4[p][q][r].D_CaM00[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer4[p][q][r].CaM00--;
      v1 = layer4[p][q][r].diffusion_lamda;
      v2 = layer4[p][q][r].chemical_lamda;
      layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
      layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);

      if(k == 0){
	layer4[p-1][q][r].CaM00++;
	v3 = layer4[p-1][q][r].diffusion_lamda;
	v4 = layer4[p-1][q][r].chemical_lamda;
	layer4[p-1][q][r].diffusion_lamda = layer4_diffusion_lamda(p-1, q, r);
	layer4[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p-1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p-1][q][r].diffusion_lamda +
	  layer4[p-1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 1){
	layer4[p][q-1][r].CaM00++;
	v3 = layer4[p][q-1][r].diffusion_lamda;
	v4 = layer4[p][q-1][r].chemical_lamda;
	layer4[p][q-1][r].diffusion_lamda = layer4_diffusion_lamda(p, q-1, r);
	layer4[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q-1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q-1][r].diffusion_lamda +
	  layer4[p][q-1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 2){
	layer4[p+1][q][r].CaM00++;
	v3 = layer4[p+1][q][r].diffusion_lamda;
	v4 = layer4[p+1][q][r].chemical_lamda;
	layer4[p+1][q][r].diffusion_lamda = layer4_diffusion_lamda(p+1, q, r);
	layer4[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p+1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p+1][q][r].diffusion_lamda +
	  layer4[p+1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 3){
	layer4[p][q+1][r].CaM00++;
	v3 = layer4[p][q+1][r].diffusion_lamda;
	v4 = layer4[p][q+1][r].chemical_lamda;
	layer4[p][q+1][r].diffusion_lamda = layer4_diffusion_lamda(p, q+1, r);
	layer4[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q+1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q+1][r].diffusion_lamda +
	  layer4[p][q+1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 4){
	if((p == 2) && (q == 0) && (r == 0)){
	  layer3[7].CaM00++;
	  v3 = layer3[7].diffusion_lamda;
	  v4 = layer3[7].chemical_lamda;
	  layer3[7].diffusion_lamda = layer3_diffusion_lamda(7);
	  layer3[7].chemical_lamda = voxel_chemical_lamda(layer3[7]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer3[7].diffusion_lamda +
	    layer3[7].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer4[p][q][r-1].CaM00++;
	  v3 = layer4[p][q][r-1].diffusion_lamda;
	  v4 = layer4[p][q][r-1].chemical_lamda;
	  layer4[p][q][r-1].diffusion_lamda = layer4_diffusion_lamda(p, q, r-1);
	  layer4[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r-1]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r-1].diffusion_lamda +
	    layer4[p][q][r-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	layer4[p][q][r+1].CaM00++;
	v3 = layer4[p][q][r+1].diffusion_lamda;
	v4 = layer4[p][q][r+1].chemical_lamda;
	layer4[p][q][r+1].diffusion_lamda = layer4_diffusion_lamda(p, q, r+1);
	layer4[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r+1]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r+1].diffusion_lamda +
	  layer4[p][q][r+1].chemical_lamda - v1 - v2 - v3 - v4;
      }
      goto loop4;
    }
   
	choose1 = choose;
    choose += layer4[p][q][r].CaM10 * layer4[p][q][r].D_CaM10[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer4[p][q][r].CaM10--;
      v1 = layer4[p][q][r].diffusion_lamda;
      v2 = layer4[p][q][r].chemical_lamda;
      layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
      layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);

      if(k == 0){
	layer4[p-1][q][r].CaM10++;
	v3 = layer4[p-1][q][r].diffusion_lamda;
	v4 = layer4[p-1][q][r].chemical_lamda;
	layer4[p-1][q][r].diffusion_lamda = layer4_diffusion_lamda(p-1, q, r);
	layer4[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p-1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p-1][q][r].diffusion_lamda +
	  layer4[p-1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 1){
	layer4[p][q-1][r].CaM10++;
	v3 = layer4[p][q-1][r].diffusion_lamda;
	v4 = layer4[p][q-1][r].chemical_lamda;
	layer4[p][q-1][r].diffusion_lamda = layer4_diffusion_lamda(p, q-1, r);
	layer4[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q-1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q-1][r].diffusion_lamda +
	  layer4[p][q-1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 2){
	layer4[p+1][q][r].CaM10++;
	v3 = layer4[p+1][q][r].diffusion_lamda;
	v4 = layer4[p+1][q][r].chemical_lamda;
	layer4[p+1][q][r].diffusion_lamda = layer4_diffusion_lamda(p+1, q, r);
	layer4[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p+1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p+1][q][r].diffusion_lamda +
	  layer4[p+1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 3){
	layer4[p][q+1][r].CaM10++;
	v3 = layer4[p][q+1][r].diffusion_lamda;
	v4 = layer4[p][q+1][r].chemical_lamda;
	layer4[p][q+1][r].diffusion_lamda = layer4_diffusion_lamda(p, q+1, r);
	layer4[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q+1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q+1][r].diffusion_lamda +
	  layer4[p][q+1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 4){
	if((p == 2) && (q == 0) && (r == 0)){
	  layer3[7].CaM10++;
	  v3 = layer3[7].diffusion_lamda;
	  v4 = layer3[7].chemical_lamda;
	  layer3[7].diffusion_lamda = layer3_diffusion_lamda(7);
	  layer3[7].chemical_lamda = voxel_chemical_lamda(layer3[7]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer3[7].diffusion_lamda +
	    layer3[7].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer4[p][q][r-1].CaM10++;
	  v3 = layer4[p][q][r-1].diffusion_lamda;
	  v4 = layer4[p][q][r-1].chemical_lamda;
	  layer4[p][q][r-1].diffusion_lamda = layer4_diffusion_lamda(p, q, r-1);
	  layer4[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r-1]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r-1].diffusion_lamda +
	    layer4[p][q][r-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	layer4[p][q][r+1].CaM10++;
	v3 = layer4[p][q][r+1].diffusion_lamda;
	v4 = layer4[p][q][r+1].chemical_lamda;
	layer4[p][q][r+1].diffusion_lamda = layer4_diffusion_lamda(p, q, r+1);
	layer4[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r+1]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r+1].diffusion_lamda +
	  layer4[p][q][r+1].chemical_lamda - v1 - v2 - v3 - v4;
      }
      goto loop4;
    }

	choose1 = choose;
    choose += layer4[p][q][r].CaM01 * layer4[p][q][r].D_CaM01[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer4[p][q][r].CaM01--;
      v1 = layer4[p][q][r].diffusion_lamda;
      v2 = layer4[p][q][r].chemical_lamda;
      layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
      layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);

      if(k == 0){
	layer4[p-1][q][r].CaM01++;
	v3 = layer4[p-1][q][r].diffusion_lamda;
	v4 = layer4[p-1][q][r].chemical_lamda;
	layer4[p-1][q][r].diffusion_lamda = layer4_diffusion_lamda(p-1, q, r);
	layer4[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p-1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p-1][q][r].diffusion_lamda +
	  layer4[p-1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 1){
	layer4[p][q-1][r].CaM01++;
	v3 = layer4[p][q-1][r].diffusion_lamda;
	v4 = layer4[p][q-1][r].chemical_lamda;
	layer4[p][q-1][r].diffusion_lamda = layer4_diffusion_lamda(p, q-1, r);
	layer4[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q-1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q-1][r].diffusion_lamda +
	  layer4[p][q-1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 2){
	layer4[p+1][q][r].CaM01++;
	v3 = layer4[p+1][q][r].diffusion_lamda;
	v4 = layer4[p+1][q][r].chemical_lamda;
	layer4[p+1][q][r].diffusion_lamda = layer4_diffusion_lamda(p+1, q, r);
	layer4[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p+1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p+1][q][r].diffusion_lamda +
	  layer4[p+1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 3){
	layer4[p][q+1][r].CaM01++;
	v3 = layer4[p][q+1][r].diffusion_lamda;
	v4 = layer4[p][q+1][r].chemical_lamda;
	layer4[p][q+1][r].diffusion_lamda = layer4_diffusion_lamda(p, q+1, r);
	layer4[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q+1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q+1][r].diffusion_lamda +
	  layer4[p][q+1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 4){
	if((p == 2) && (q == 0) && (r == 0)){
	  layer3[7].CaM01++;
	  v3 = layer3[7].diffusion_lamda;
	  v4 = layer3[7].chemical_lamda;
	  layer3[7].diffusion_lamda = layer3_diffusion_lamda(7);
	  layer3[7].chemical_lamda = voxel_chemical_lamda(layer3[7]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer3[7].diffusion_lamda +
	    layer3[7].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer4[p][q][r-1].CaM01++;
	  v3 = layer4[p][q][r-1].diffusion_lamda;
	  v4 = layer4[p][q][r-1].chemical_lamda;
	  layer4[p][q][r-1].diffusion_lamda = layer4_diffusion_lamda(p, q, r-1);
	  layer4[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r-1]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r-1].diffusion_lamda +
	    layer4[p][q][r-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	layer4[p][q][r+1].CaM01++;
	v3 = layer4[p][q][r+1].diffusion_lamda;
	v4 = layer4[p][q][r+1].chemical_lamda;
	layer4[p][q][r+1].diffusion_lamda = layer4_diffusion_lamda(p, q, r+1);
	layer4[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r+1]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r+1].diffusion_lamda +
	  layer4[p][q][r+1].chemical_lamda - v1 - v2 - v3 - v4;
      }
      goto loop4;
    }

	choose1 = choose;
    choose += layer4[p][q][r].CaM11 * layer4[p][q][r].D_CaM11[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer4[p][q][r].CaM11--;
      v1 = layer4[p][q][r].diffusion_lamda;
      v2 = layer4[p][q][r].chemical_lamda;
      layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
      layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);

      if(k == 0){
	layer4[p-1][q][r].CaM11++;
	v3 = layer4[p-1][q][r].diffusion_lamda;
	v4 = layer4[p-1][q][r].chemical_lamda;
	layer4[p-1][q][r].diffusion_lamda = layer4_diffusion_lamda(p-1, q, r);
	layer4[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p-1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p-1][q][r].diffusion_lamda +
	  layer4[p-1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 1){
	layer4[p][q-1][r].CaM11++;
	v3 = layer4[p][q-1][r].diffusion_lamda;
	v4 = layer4[p][q-1][r].chemical_lamda;
	layer4[p][q-1][r].diffusion_lamda = layer4_diffusion_lamda(p, q-1, r);
	layer4[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q-1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q-1][r].diffusion_lamda +
	  layer4[p][q-1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 2){
	layer4[p+1][q][r].CaM11++;
	v3 = layer4[p+1][q][r].diffusion_lamda;
	v4 = layer4[p+1][q][r].chemical_lamda;
	layer4[p+1][q][r].diffusion_lamda = layer4_diffusion_lamda(p+1, q, r);
	layer4[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p+1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p+1][q][r].diffusion_lamda +
	  layer4[p+1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 3){
	layer4[p][q+1][r].CaM11++;
	v3 = layer4[p][q+1][r].diffusion_lamda;
	v4 = layer4[p][q+1][r].chemical_lamda;
	layer4[p][q+1][r].diffusion_lamda = layer4_diffusion_lamda(p, q+1, r);
	layer4[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q+1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q+1][r].diffusion_lamda +
	  layer4[p][q+1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 4){
	if((p == 2) && (q == 0) && (r == 0)){
	  layer3[7].CaM11++;
	  v3 = layer3[7].diffusion_lamda;
	  v4 = layer3[7].chemical_lamda;
	  layer3[7].diffusion_lamda = layer3_diffusion_lamda(7);
	  layer3[7].chemical_lamda = voxel_chemical_lamda(layer3[7]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer3[7].diffusion_lamda +
	    layer3[7].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer4[p][q][r-1].CaM11++;
	  v3 = layer4[p][q][r-1].diffusion_lamda;
	  v4 = layer4[p][q][r-1].chemical_lamda;
	  layer4[p][q][r-1].diffusion_lamda = layer4_diffusion_lamda(p, q, r-1);
	  layer4[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r-1]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r-1].diffusion_lamda +
	    layer4[p][q][r-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	layer4[p][q][r+1].CaM11++;
	v3 = layer4[p][q][r+1].diffusion_lamda;
	v4 = layer4[p][q][r+1].chemical_lamda;
	layer4[p][q][r+1].diffusion_lamda = layer4_diffusion_lamda(p, q, r+1);
	layer4[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r+1]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r+1].diffusion_lamda +
	  layer4[p][q][r+1].chemical_lamda - v1 - v2 - v3 - v4;
      }
      goto loop4;
    }

	choose1 = choose;
    choose += layer4[p][q][r].CaM20 * layer4[p][q][r].D_CaM20[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer4[p][q][r].CaM20--;
      v1 = layer4[p][q][r].diffusion_lamda;
      v2 = layer4[p][q][r].chemical_lamda;
      layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
      layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);

      if(k == 0){
	layer4[p-1][q][r].CaM20++;
	v3 = layer4[p-1][q][r].diffusion_lamda;
	v4 = layer4[p-1][q][r].chemical_lamda;
	layer4[p-1][q][r].diffusion_lamda = layer4_diffusion_lamda(p-1, q, r);
	layer4[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p-1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p-1][q][r].diffusion_lamda +
	  layer4[p-1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 1){
	layer4[p][q-1][r].CaM20++;
	v3 = layer4[p][q-1][r].diffusion_lamda;
	v4 = layer4[p][q-1][r].chemical_lamda;
	layer4[p][q-1][r].diffusion_lamda = layer4_diffusion_lamda(p, q-1, r);
	layer4[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q-1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q-1][r].diffusion_lamda +
	  layer4[p][q-1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 2){
	layer4[p+1][q][r].CaM20++;
	v3 = layer4[p+1][q][r].diffusion_lamda;
	v4 = layer4[p+1][q][r].chemical_lamda;
	layer4[p+1][q][r].diffusion_lamda = layer4_diffusion_lamda(p+1, q, r);
	layer4[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p+1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p+1][q][r].diffusion_lamda +
	  layer4[p+1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 3){
	layer4[p][q+1][r].CaM20++;
	v3 = layer4[p][q+1][r].diffusion_lamda;
	v4 = layer4[p][q+1][r].chemical_lamda;
	layer4[p][q+1][r].diffusion_lamda = layer4_diffusion_lamda(p, q+1, r);
	layer4[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q+1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q+1][r].diffusion_lamda +
	  layer4[p][q+1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 4){
	if((p == 2) && (q == 0) && (r == 0)){
	  layer3[7].CaM20++;
	  v3 = layer3[7].diffusion_lamda;
	  v4 = layer3[7].chemical_lamda;
	  layer3[7].diffusion_lamda = layer3_diffusion_lamda(7);
	  layer3[7].chemical_lamda = voxel_chemical_lamda(layer3[7]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer3[7].diffusion_lamda +
	    layer3[7].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer4[p][q][r-1].CaM20++;
	  v3 = layer4[p][q][r-1].diffusion_lamda;
	  v4 = layer4[p][q][r-1].chemical_lamda;
	  layer4[p][q][r-1].diffusion_lamda = layer4_diffusion_lamda(p, q, r-1);
	  layer4[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r-1]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r-1].diffusion_lamda +
	    layer4[p][q][r-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	layer4[p][q][r+1].CaM20++;
	v3 = layer4[p][q][r+1].diffusion_lamda;
	v4 = layer4[p][q][r+1].chemical_lamda;
	layer4[p][q][r+1].diffusion_lamda = layer4_diffusion_lamda(p, q, r+1);
	layer4[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r+1]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r+1].diffusion_lamda +
	  layer4[p][q][r+1].chemical_lamda - v1 - v2 - v3 - v4;
      }
      goto loop4;
    }

	choose1 = choose;
    choose += layer4[p][q][r].CaM02 * layer4[p][q][r].D_CaM02[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer4[p][q][r].CaM02--;
      v1 = layer4[p][q][r].diffusion_lamda;
      v2 = layer4[p][q][r].chemical_lamda;
      layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
      layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);

      if(k == 0){
	layer4[p-1][q][r].CaM02++;
	v3 = layer4[p-1][q][r].diffusion_lamda;
	v4 = layer4[p-1][q][r].chemical_lamda;
	layer4[p-1][q][r].diffusion_lamda = layer4_diffusion_lamda(p-1, q, r);
	layer4[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p-1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p-1][q][r].diffusion_lamda +
	  layer4[p-1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 1){
	layer4[p][q-1][r].CaM02++;
	v3 = layer4[p][q-1][r].diffusion_lamda;
	v4 = layer4[p][q-1][r].chemical_lamda;
	layer4[p][q-1][r].diffusion_lamda = layer4_diffusion_lamda(p, q-1, r);
	layer4[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q-1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q-1][r].diffusion_lamda +
	  layer4[p][q-1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 2){
	layer4[p+1][q][r].CaM02++;
	v3 = layer4[p+1][q][r].diffusion_lamda;
	v4 = layer4[p+1][q][r].chemical_lamda;
	layer4[p+1][q][r].diffusion_lamda = layer4_diffusion_lamda(p+1, q, r);
	layer4[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p+1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p+1][q][r].diffusion_lamda +
	  layer4[p+1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 3){
	layer4[p][q+1][r].CaM02++;
	v3 = layer4[p][q+1][r].diffusion_lamda;
	v4 = layer4[p][q+1][r].chemical_lamda;
	layer4[p][q+1][r].diffusion_lamda = layer4_diffusion_lamda(p, q+1, r);
	layer4[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q+1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q+1][r].diffusion_lamda +
	  layer4[p][q+1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 4){
	if((p == 2) && (q == 0) && (r == 0)){
	  layer3[7].CaM02++;
	  v3 = layer3[7].diffusion_lamda;
	  v4 = layer3[7].chemical_lamda;
	  layer3[7].diffusion_lamda = layer3_diffusion_lamda(7);
	  layer3[7].chemical_lamda = voxel_chemical_lamda(layer3[7]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer3[7].diffusion_lamda +
	    layer3[7].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer4[p][q][r-1].CaM02++;
	  v3 = layer4[p][q][r-1].diffusion_lamda;
	  v4 = layer4[p][q][r-1].chemical_lamda;
	  layer4[p][q][r-1].diffusion_lamda = layer4_diffusion_lamda(p, q, r-1);
	  layer4[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r-1]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r-1].diffusion_lamda +
	    layer4[p][q][r-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	layer4[p][q][r+1].CaM02++;
	v3 = layer4[p][q][r+1].diffusion_lamda;
	v4 = layer4[p][q][r+1].chemical_lamda;
	layer4[p][q][r+1].diffusion_lamda = layer4_diffusion_lamda(p, q, r+1);
	layer4[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r+1]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r+1].diffusion_lamda +
	  layer4[p][q][r+1].chemical_lamda - v1 - v2 - v3 - v4;
      }
      goto loop4;
    }

	choose1 = choose;
    choose += layer4[p][q][r].CaM21 * layer4[p][q][r].D_CaM21[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer4[p][q][r].CaM21--;
      v1 = layer4[p][q][r].diffusion_lamda;
      v2 = layer4[p][q][r].chemical_lamda;
      layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
      layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);

      if(k == 0){
	layer4[p-1][q][r].CaM21++;
	v3 = layer4[p-1][q][r].diffusion_lamda;
	v4 = layer4[p-1][q][r].chemical_lamda;
	layer4[p-1][q][r].diffusion_lamda = layer4_diffusion_lamda(p-1, q, r);
	layer4[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p-1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p-1][q][r].diffusion_lamda +
	  layer4[p-1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 1){
	layer4[p][q-1][r].CaM21++;
	v3 = layer4[p][q-1][r].diffusion_lamda;
	v4 = layer4[p][q-1][r].chemical_lamda;
	layer4[p][q-1][r].diffusion_lamda = layer4_diffusion_lamda(p, q-1, r);
	layer4[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q-1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q-1][r].diffusion_lamda +
	  layer4[p][q-1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 2){
	layer4[p+1][q][r].CaM21++;
	v3 = layer4[p+1][q][r].diffusion_lamda;
	v4 = layer4[p+1][q][r].chemical_lamda;
	layer4[p+1][q][r].diffusion_lamda = layer4_diffusion_lamda(p+1, q, r);
	layer4[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p+1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p+1][q][r].diffusion_lamda +
	  layer4[p+1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 3){
	layer4[p][q+1][r].CaM21++;
	v3 = layer4[p][q+1][r].diffusion_lamda;
	v4 = layer4[p][q+1][r].chemical_lamda;
	layer4[p][q+1][r].diffusion_lamda = layer4_diffusion_lamda(p, q+1, r);
	layer4[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q+1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q+1][r].diffusion_lamda +
	  layer4[p][q+1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 4){
	if((p == 2) && (q == 0) && (r == 0)){
	  layer3[7].CaM21++;
	  v3 = layer3[7].diffusion_lamda;
	  v4 = layer3[7].chemical_lamda;
	  layer3[7].diffusion_lamda = layer3_diffusion_lamda(7);
	  layer3[7].chemical_lamda = voxel_chemical_lamda(layer3[7]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer3[7].diffusion_lamda +
	    layer3[7].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer4[p][q][r-1].CaM21++;
	  v3 = layer4[p][q][r-1].diffusion_lamda;
	  v4 = layer4[p][q][r-1].chemical_lamda;
	  layer4[p][q][r-1].diffusion_lamda = layer4_diffusion_lamda(p, q, r-1);
	  layer4[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r-1]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r-1].diffusion_lamda +
	    layer4[p][q][r-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	layer4[p][q][r+1].CaM21++;
	v3 = layer4[p][q][r+1].diffusion_lamda;
	v4 = layer4[p][q][r+1].chemical_lamda;
	layer4[p][q][r+1].diffusion_lamda = layer4_diffusion_lamda(p, q, r+1);
	layer4[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r+1]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r+1].diffusion_lamda +
	  layer4[p][q][r+1].chemical_lamda - v1 - v2 - v3 - v4;
      }
      goto loop4;
    }

	choose1 = choose;
    choose += layer4[p][q][r].CaM12 * layer4[p][q][r].D_CaM12[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer4[p][q][r].CaM12--;
      v1 = layer4[p][q][r].diffusion_lamda;
      v2 = layer4[p][q][r].chemical_lamda;
      layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
      layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);

      if(k == 0){
	layer4[p-1][q][r].CaM12++;
	v3 = layer4[p-1][q][r].diffusion_lamda;
	v4 = layer4[p-1][q][r].chemical_lamda;
	layer4[p-1][q][r].diffusion_lamda = layer4_diffusion_lamda(p-1, q, r);
	layer4[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p-1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p-1][q][r].diffusion_lamda +
	  layer4[p-1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 1){
	layer4[p][q-1][r].CaM12++;
	v3 = layer4[p][q-1][r].diffusion_lamda;
	v4 = layer4[p][q-1][r].chemical_lamda;
	layer4[p][q-1][r].diffusion_lamda = layer4_diffusion_lamda(p, q-1, r);
	layer4[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q-1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q-1][r].diffusion_lamda +
	  layer4[p][q-1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 2){
	layer4[p+1][q][r].CaM12++;
	v3 = layer4[p+1][q][r].diffusion_lamda;
	v4 = layer4[p+1][q][r].chemical_lamda;
	layer4[p+1][q][r].diffusion_lamda = layer4_diffusion_lamda(p+1, q, r);
	layer4[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p+1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p+1][q][r].diffusion_lamda +
	  layer4[p+1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 3){
	layer4[p][q+1][r].CaM12++;
	v3 = layer4[p][q+1][r].diffusion_lamda;
	v4 = layer4[p][q+1][r].chemical_lamda;
	layer4[p][q+1][r].diffusion_lamda = layer4_diffusion_lamda(p, q+1, r);
	layer4[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q+1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q+1][r].diffusion_lamda +
	  layer4[p][q+1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 4){
	if((p == 2) && (q == 0) && (r == 0)){
	  layer3[7].CaM12++;
	  v3 = layer3[7].diffusion_lamda;
	  v4 = layer3[7].chemical_lamda;
	  layer3[7].diffusion_lamda = layer3_diffusion_lamda(7);
	  layer3[7].chemical_lamda = voxel_chemical_lamda(layer3[7]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer3[7].diffusion_lamda +
	    layer3[7].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer4[p][q][r-1].CaM12++;
	  v3 = layer4[p][q][r-1].diffusion_lamda;
	  v4 = layer4[p][q][r-1].chemical_lamda;
	  layer4[p][q][r-1].diffusion_lamda = layer4_diffusion_lamda(p, q, r-1);
	  layer4[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r-1]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r-1].diffusion_lamda +
	    layer4[p][q][r-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	layer4[p][q][r+1].CaM12++;
	v3 = layer4[p][q][r+1].diffusion_lamda;
	v4 = layer4[p][q][r+1].chemical_lamda;
	layer4[p][q][r+1].diffusion_lamda = layer4_diffusion_lamda(p, q, r+1);
	layer4[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r+1]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r+1].diffusion_lamda +
	  layer4[p][q][r+1].chemical_lamda - v1 - v2 - v3 - v4;
      }
      goto loop4;
    }

	choose1 = choose;
    choose += layer4[p][q][r].CaM22 * layer4[p][q][r].D_CaM22[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer4[p][q][r].CaM22--;
      v1 = layer4[p][q][r].diffusion_lamda;
      v2 = layer4[p][q][r].chemical_lamda;
      layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
      layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);

      if(k == 0){
	layer4[p-1][q][r].CaM22++;
	v3 = layer4[p-1][q][r].diffusion_lamda;
	v4 = layer4[p-1][q][r].chemical_lamda;
	layer4[p-1][q][r].diffusion_lamda = layer4_diffusion_lamda(p-1, q, r);
	layer4[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p-1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p-1][q][r].diffusion_lamda +
	  layer4[p-1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 1){
	layer4[p][q-1][r].CaM22++;
	v3 = layer4[p][q-1][r].diffusion_lamda;
	v4 = layer4[p][q-1][r].chemical_lamda;
	layer4[p][q-1][r].diffusion_lamda = layer4_diffusion_lamda(p, q-1, r);
	layer4[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q-1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q-1][r].diffusion_lamda +
	  layer4[p][q-1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 2){
	layer4[p+1][q][r].CaM22++;
	v3 = layer4[p+1][q][r].diffusion_lamda;
	v4 = layer4[p+1][q][r].chemical_lamda;
	layer4[p+1][q][r].diffusion_lamda = layer4_diffusion_lamda(p+1, q, r);
	layer4[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p+1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p+1][q][r].diffusion_lamda +
	  layer4[p+1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 3){
	layer4[p][q+1][r].CaM22++;
	v3 = layer4[p][q+1][r].diffusion_lamda;
	v4 = layer4[p][q+1][r].chemical_lamda;
	layer4[p][q+1][r].diffusion_lamda = layer4_diffusion_lamda(p, q+1, r);
	layer4[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q+1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q+1][r].diffusion_lamda +
	  layer4[p][q+1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 4){
	if((p == 2) && (q == 0) && (r == 0)){
	  layer3[7].CaM22++;
	  v3 = layer3[7].diffusion_lamda;
	  v4 = layer3[7].chemical_lamda;
	  layer3[7].diffusion_lamda = layer3_diffusion_lamda(7);
	  layer3[7].chemical_lamda = voxel_chemical_lamda(layer3[7]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer3[7].diffusion_lamda +
	    layer3[7].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer4[p][q][r-1].CaM22++;
	  v3 = layer4[p][q][r-1].diffusion_lamda;
	  v4 = layer4[p][q][r-1].chemical_lamda;
	  layer4[p][q][r-1].diffusion_lamda = layer4_diffusion_lamda(p, q, r-1);
	  layer4[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r-1]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r-1].diffusion_lamda +
	    layer4[p][q][r-1].chemical_lamda - v1 - v2 - v3 - v4;
	}
      }else{
	layer4[p][q][r+1].CaM22++;
	v3 = layer4[p][q][r+1].diffusion_lamda;
	v4 = layer4[p][q][r+1].chemical_lamda;
	layer4[p][q][r+1].diffusion_lamda = layer4_diffusion_lamda(p, q, r+1);
	layer4[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r+1]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r+1].diffusion_lamda +
	  layer4[p][q][r+1].chemical_lamda - v1 - v2 - v3 - v4;
      }
      goto loop4;
    }
    
	choose1 = choose;
    choose += layer4[p][q][r].CB00 * layer4[p][q][r].D_CB00[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer4[p][q][r].CB00--;
      v1 = layer4[p][q][r].diffusion_lamda;
      v2 = layer4[p][q][r].chemical_lamda;
      layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
      layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);

      if(k == 0){
	layer4[p-1][q][r].CB00++;
	v3 = layer4[p-1][q][r].diffusion_lamda;
	v4 = layer4[p-1][q][r].chemical_lamda;
	layer4[p-1][q][r].diffusion_lamda = layer4_diffusion_lamda(p-1, q, r);
	layer4[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p-1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p-1][q][r].diffusion_lamda +
	  layer4[p-1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 1){
	layer4[p][q-1][r].CB00++;
	v3 = layer4[p][q-1][r].diffusion_lamda;
	v4 = layer4[p][q-1][r].chemical_lamda;
	layer4[p][q-1][r].diffusion_lamda = layer4_diffusion_lamda(p, q-1, r);
	layer4[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q-1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q-1][r].diffusion_lamda +
	  layer4[p][q-1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 2){
	layer4[p+1][q][r].CB00++;
	v3 = layer4[p+1][q][r].diffusion_lamda;
	v4 = layer4[p+1][q][r].chemical_lamda;
	layer4[p+1][q][r].diffusion_lamda = layer4_diffusion_lamda(p+1, q, r);
	layer4[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p+1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p+1][q][r].diffusion_lamda +
	  layer4[p+1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 3){
	layer4[p][q+1][r].CB00++;
	v3 = layer4[p][q+1][r].diffusion_lamda;
	v4 = layer4[p][q+1][r].chemical_lamda;
	layer4[p][q+1][r].diffusion_lamda = layer4_diffusion_lamda(p, q+1, r);
	layer4[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q+1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q+1][r].diffusion_lamda +
	  layer4[p][q+1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 4){
	if((p == 2) && (q == 0) && (r == 0)){
	  layer3[7].CB00++;
	  v3 = layer3[7].diffusion_lamda;
	  v4 = layer3[7].chemical_lamda;
	  layer3[7].diffusion_lamda = layer3_diffusion_lamda(7);
	  layer3[7].chemical_lamda = voxel_chemical_lamda(layer3[7]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer3[7].diffusion_lamda +
	    layer3[7].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer4[p][q][r-1].CB00++;
	  v3 = layer4[p][q][r-1].diffusion_lamda;
	  v4 = layer4[p][q][r-1].chemical_lamda;
	  layer4[p][q][r-1].diffusion_lamda = layer4_diffusion_lamda(p, q, r-1);
	  layer4[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r-1]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r-1].diffusion_lamda +
	    layer4[p][q][r-1].chemical_lamda - v1 - v2 - v3 - v4;

	}
      }else{
	layer4[p][q][r+1].CB00++;
	v3 = layer4[p][q][r+1].diffusion_lamda;
	v4 = layer4[p][q][r+1].chemical_lamda;
	layer4[p][q][r+1].diffusion_lamda = layer4_diffusion_lamda(p, q, r+1);
	layer4[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r+1]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r+1].diffusion_lamda +
	  layer4[p][q][r+1].chemical_lamda - v1 - v2 - v3 - v4;
      }
      goto loop4;
    }

	choose1 = choose;
    choose += layer4[p][q][r].CB10 * layer4[p][q][r].D_CB10[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer4[p][q][r].CB10--;
      v1 = layer4[p][q][r].diffusion_lamda;
      v2 = layer4[p][q][r].chemical_lamda;
      layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
      layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);

      if(k == 0){
	layer4[p-1][q][r].CB10++;
	v3 = layer4[p-1][q][r].diffusion_lamda;
	v4 = layer4[p-1][q][r].chemical_lamda;
	layer4[p-1][q][r].diffusion_lamda = layer4_diffusion_lamda(p-1, q, r);
	layer4[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p-1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p-1][q][r].diffusion_lamda +
	  layer4[p-1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 1){
	layer4[p][q-1][r].CB10++;
	v3 = layer4[p][q-1][r].diffusion_lamda;
	v4 = layer4[p][q-1][r].chemical_lamda;
	layer4[p][q-1][r].diffusion_lamda = layer4_diffusion_lamda(p, q-1, r);
	layer4[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q-1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q-1][r].diffusion_lamda +
	  layer4[p][q-1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 2){
	layer4[p+1][q][r].CB10++;
	v3 = layer4[p+1][q][r].diffusion_lamda;
	v4 = layer4[p+1][q][r].chemical_lamda;
	layer4[p+1][q][r].diffusion_lamda = layer4_diffusion_lamda(p+1, q, r);
	layer4[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p+1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p+1][q][r].diffusion_lamda +
	  layer4[p+1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 3){
	layer4[p][q+1][r].CB10++;
	v3 = layer4[p][q+1][r].diffusion_lamda;
	v4 = layer4[p][q+1][r].chemical_lamda;
	layer4[p][q+1][r].diffusion_lamda = layer4_diffusion_lamda(p, q+1, r);
	layer4[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q+1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q+1][r].diffusion_lamda +
	  layer4[p][q+1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 4){
	if((p == 2) && (q == 0) && (r == 0)){
	  layer3[7].CB10++;
	  v3 = layer3[7].diffusion_lamda;
	  v4 = layer3[7].chemical_lamda;
	  layer3[7].diffusion_lamda = layer3_diffusion_lamda(7);
	  layer3[7].chemical_lamda = voxel_chemical_lamda(layer3[7]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer3[7].diffusion_lamda +
	    layer3[7].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer4[p][q][r-1].CB10++;
	  v3 = layer4[p][q][r-1].diffusion_lamda;
	  v4 = layer4[p][q][r-1].chemical_lamda;
	  layer4[p][q][r-1].diffusion_lamda = layer4_diffusion_lamda(p, q, r-1);
	  layer4[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r-1]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r-1].diffusion_lamda +
	    layer4[p][q][r-1].chemical_lamda - v1 - v2 - v3 - v4;

	}
      }else{
	layer4[p][q][r+1].CB10++;
	v3 = layer4[p][q][r+1].diffusion_lamda;
	v4 = layer4[p][q][r+1].chemical_lamda;
	layer4[p][q][r+1].diffusion_lamda = layer4_diffusion_lamda(p, q, r+1);
	layer4[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r+1]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r+1].diffusion_lamda +
	  layer4[p][q][r+1].chemical_lamda - v1 - v2 - v3 - v4;
      }
      goto loop4;
    }

	choose1 = choose;
    choose += layer4[p][q][r].CB01 * layer4[p][q][r].D_CB01[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer4[p][q][r].CB01--;
      v1 = layer4[p][q][r].diffusion_lamda;
      v2 = layer4[p][q][r].chemical_lamda;
      layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
      layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);

      if(k == 0){
	layer4[p-1][q][r].CB01++;
	v3 = layer4[p-1][q][r].diffusion_lamda;
	v4 = layer4[p-1][q][r].chemical_lamda;
	layer4[p-1][q][r].diffusion_lamda = layer4_diffusion_lamda(p-1, q, r);
	layer4[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p-1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p-1][q][r].diffusion_lamda +
	  layer4[p-1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 1){
	layer4[p][q-1][r].CB01++;
	v3 = layer4[p][q-1][r].diffusion_lamda;
	v4 = layer4[p][q-1][r].chemical_lamda;
	layer4[p][q-1][r].diffusion_lamda = layer4_diffusion_lamda(p, q-1, r);
	layer4[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q-1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q-1][r].diffusion_lamda +
	  layer4[p][q-1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 2){
	layer4[p+1][q][r].CB01++;
	v3 = layer4[p+1][q][r].diffusion_lamda;
	v4 = layer4[p+1][q][r].chemical_lamda;
	layer4[p+1][q][r].diffusion_lamda = layer4_diffusion_lamda(p+1, q, r);
	layer4[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p+1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p+1][q][r].diffusion_lamda +
	  layer4[p+1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 3){
	layer4[p][q+1][r].CB01++;
	v3 = layer4[p][q+1][r].diffusion_lamda;
	v4 = layer4[p][q+1][r].chemical_lamda;
	layer4[p][q+1][r].diffusion_lamda = layer4_diffusion_lamda(p, q+1, r);
	layer4[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q+1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q+1][r].diffusion_lamda +
	  layer4[p][q+1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 4){
	if((p == 2) && (q == 0) && (r == 0)){
	  layer3[7].CB01++;
	  v3 = layer3[7].diffusion_lamda;
	  v4 = layer3[7].chemical_lamda;
	  layer3[7].diffusion_lamda = layer3_diffusion_lamda(7);
	  layer3[7].chemical_lamda = voxel_chemical_lamda(layer3[7]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer3[7].diffusion_lamda +
	    layer3[7].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer4[p][q][r-1].CB01++;
	  v3 = layer4[p][q][r-1].diffusion_lamda;
	  v4 = layer4[p][q][r-1].chemical_lamda;
	  layer4[p][q][r-1].diffusion_lamda = layer4_diffusion_lamda(p, q, r-1);
	  layer4[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r-1]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r-1].diffusion_lamda +
	    layer4[p][q][r-1].chemical_lamda - v1 - v2 - v3 - v4;

	}
      }else{
	layer4[p][q][r+1].CB01++;
	v3 = layer4[p][q][r+1].diffusion_lamda;
	v4 = layer4[p][q][r+1].chemical_lamda;
	layer4[p][q][r+1].diffusion_lamda = layer4_diffusion_lamda(p, q, r+1);
	layer4[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r+1]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r+1].diffusion_lamda +
	  layer4[p][q][r+1].chemical_lamda - v1 - v2 - v3 - v4;
      }
      goto loop4;
    }

	choose1 = choose;
    choose += layer4[p][q][r].CB11 * layer4[p][q][r].D_CB11[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer4[p][q][r].CB11--;
      v1 = layer4[p][q][r].diffusion_lamda;
      v2 = layer4[p][q][r].chemical_lamda;
      layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
      layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);

      if(k == 0){
	layer4[p-1][q][r].CB11++;
	v3 = layer4[p-1][q][r].diffusion_lamda;
	v4 = layer4[p-1][q][r].chemical_lamda;
	layer4[p-1][q][r].diffusion_lamda = layer4_diffusion_lamda(p-1, q, r);
	layer4[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p-1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p-1][q][r].diffusion_lamda +
	  layer4[p-1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 1){
	layer4[p][q-1][r].CB11++;
	v3 = layer4[p][q-1][r].diffusion_lamda;
	v4 = layer4[p][q-1][r].chemical_lamda;
	layer4[p][q-1][r].diffusion_lamda = layer4_diffusion_lamda(p, q-1, r);
	layer4[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q-1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q-1][r].diffusion_lamda +
	  layer4[p][q-1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 2){
	layer4[p+1][q][r].CB11++;
	v3 = layer4[p+1][q][r].diffusion_lamda;
	v4 = layer4[p+1][q][r].chemical_lamda;
	layer4[p+1][q][r].diffusion_lamda = layer4_diffusion_lamda(p+1, q, r);
	layer4[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p+1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p+1][q][r].diffusion_lamda +
	  layer4[p+1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 3){
	layer4[p][q+1][r].CB11++;
	v3 = layer4[p][q+1][r].diffusion_lamda;
	v4 = layer4[p][q+1][r].chemical_lamda;
	layer4[p][q+1][r].diffusion_lamda = layer4_diffusion_lamda(p, q+1, r);
	layer4[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q+1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q+1][r].diffusion_lamda +
	  layer4[p][q+1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 4){
	if((p == 2) && (q == 0) && (r == 0)){
	  layer3[7].CB11++;
	  v3 = layer3[7].diffusion_lamda;
	  v4 = layer3[7].chemical_lamda;
	  layer3[7].diffusion_lamda = layer3_diffusion_lamda(7);
	  layer3[7].chemical_lamda = voxel_chemical_lamda(layer3[7]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer3[7].diffusion_lamda +
	    layer3[7].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer4[p][q][r-1].CB11++;
	  v3 = layer4[p][q][r-1].diffusion_lamda;
	  v4 = layer4[p][q][r-1].chemical_lamda;
	  layer4[p][q][r-1].diffusion_lamda = layer4_diffusion_lamda(p, q, r-1);
	  layer4[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r-1]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r-1].diffusion_lamda +
	    layer4[p][q][r-1].chemical_lamda - v1 - v2 - v3 - v4;

	}
      }else{
	layer4[p][q][r+1].CB11++;
	v3 = layer4[p][q][r+1].diffusion_lamda;
	v4 = layer4[p][q][r+1].chemical_lamda;
	layer4[p][q][r+1].diffusion_lamda = layer4_diffusion_lamda(p, q, r+1);
	layer4[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r+1]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r+1].diffusion_lamda +
	  layer4[p][q][r+1].chemical_lamda - v1 - v2 - v3 - v4;
      }
      goto loop4;
    }

	choose1 = choose;
    choose += layer4[p][q][r].CB20 * layer4[p][q][r].D_CB20[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer4[p][q][r].CB20--;
      v1 = layer4[p][q][r].diffusion_lamda;
      v2 = layer4[p][q][r].chemical_lamda;
      layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
      layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);

      if(k == 0){
	layer4[p-1][q][r].CB20++;
	v3 = layer4[p-1][q][r].diffusion_lamda;
	v4 = layer4[p-1][q][r].chemical_lamda;
	layer4[p-1][q][r].diffusion_lamda = layer4_diffusion_lamda(p-1, q, r);
	layer4[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p-1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p-1][q][r].diffusion_lamda +
	  layer4[p-1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 1){
	layer4[p][q-1][r].CB20++;
	v3 = layer4[p][q-1][r].diffusion_lamda;
	v4 = layer4[p][q-1][r].chemical_lamda;
	layer4[p][q-1][r].diffusion_lamda = layer4_diffusion_lamda(p, q-1, r);
	layer4[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q-1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q-1][r].diffusion_lamda +
	  layer4[p][q-1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 2){
	layer4[p+1][q][r].CB20++;
	v3 = layer4[p+1][q][r].diffusion_lamda;
	v4 = layer4[p+1][q][r].chemical_lamda;
	layer4[p+1][q][r].diffusion_lamda = layer4_diffusion_lamda(p+1, q, r);
	layer4[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p+1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p+1][q][r].diffusion_lamda +
	  layer4[p+1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 3){
	layer4[p][q+1][r].CB20++;
	v3 = layer4[p][q+1][r].diffusion_lamda;
	v4 = layer4[p][q+1][r].chemical_lamda;
	layer4[p][q+1][r].diffusion_lamda = layer4_diffusion_lamda(p, q+1, r);
	layer4[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q+1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q+1][r].diffusion_lamda +
	  layer4[p][q+1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 4){
	if((p == 2) && (q == 0) && (r == 0)){
	  layer3[7].CB20++;
	  v3 = layer3[7].diffusion_lamda;
	  v4 = layer3[7].chemical_lamda;
	  layer3[7].diffusion_lamda = layer3_diffusion_lamda(7);
	  layer3[7].chemical_lamda = voxel_chemical_lamda(layer3[7]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer3[7].diffusion_lamda +
	    layer3[7].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer4[p][q][r-1].CB20++;
	  v3 = layer4[p][q][r-1].diffusion_lamda;
	  v4 = layer4[p][q][r-1].chemical_lamda;
	  layer4[p][q][r-1].diffusion_lamda = layer4_diffusion_lamda(p, q, r-1);
	  layer4[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r-1]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r-1].diffusion_lamda +
	    layer4[p][q][r-1].chemical_lamda - v1 - v2 - v3 - v4;

	}
      }else{
	layer4[p][q][r+1].CB20++;
	v3 = layer4[p][q][r+1].diffusion_lamda;
	v4 = layer4[p][q][r+1].chemical_lamda;
	layer4[p][q][r+1].diffusion_lamda = layer4_diffusion_lamda(p, q, r+1);
	layer4[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r+1]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r+1].diffusion_lamda +
	  layer4[p][q][r+1].chemical_lamda - v1 - v2 - v3 - v4;
      }
      goto loop4;
    }

	choose1 = choose;
    choose += layer4[p][q][r].CB02 * layer4[p][q][r].D_CB02[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer4[p][q][r].CB02--;
      v1 = layer4[p][q][r].diffusion_lamda;
      v2 = layer4[p][q][r].chemical_lamda;
      layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
      layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);

      if(k == 0){
	layer4[p-1][q][r].CB02++;
	v3 = layer4[p-1][q][r].diffusion_lamda;
	v4 = layer4[p-1][q][r].chemical_lamda;
	layer4[p-1][q][r].diffusion_lamda = layer4_diffusion_lamda(p-1, q, r);
	layer4[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p-1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p-1][q][r].diffusion_lamda +
	  layer4[p-1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 1){
	layer4[p][q-1][r].CB02++;
	v3 = layer4[p][q-1][r].diffusion_lamda;
	v4 = layer4[p][q-1][r].chemical_lamda;
	layer4[p][q-1][r].diffusion_lamda = layer4_diffusion_lamda(p, q-1, r);
	layer4[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q-1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q-1][r].diffusion_lamda +
	  layer4[p][q-1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 2){
	layer4[p+1][q][r].CB02++;
	v3 = layer4[p+1][q][r].diffusion_lamda;
	v4 = layer4[p+1][q][r].chemical_lamda;
	layer4[p+1][q][r].diffusion_lamda = layer4_diffusion_lamda(p+1, q, r);
	layer4[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p+1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p+1][q][r].diffusion_lamda +
	  layer4[p+1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 3){
	layer4[p][q+1][r].CB02++;
	v3 = layer4[p][q+1][r].diffusion_lamda;
	v4 = layer4[p][q+1][r].chemical_lamda;
	layer4[p][q+1][r].diffusion_lamda = layer4_diffusion_lamda(p, q+1, r);
	layer4[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q+1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q+1][r].diffusion_lamda +
	  layer4[p][q+1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 4){
	if((p == 2) && (q == 0) && (r == 0)){
	  layer3[7].CB02++;
	  v3 = layer3[7].diffusion_lamda;
	  v4 = layer3[7].chemical_lamda;
	  layer3[7].diffusion_lamda = layer3_diffusion_lamda(7);
	  layer3[7].chemical_lamda = voxel_chemical_lamda(layer3[7]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer3[7].diffusion_lamda +
	    layer3[7].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer4[p][q][r-1].CB02++;
	  v3 = layer4[p][q][r-1].diffusion_lamda;
	  v4 = layer4[p][q][r-1].chemical_lamda;
	  layer4[p][q][r-1].diffusion_lamda = layer4_diffusion_lamda(p, q, r-1);
	  layer4[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r-1]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r-1].diffusion_lamda +
	    layer4[p][q][r-1].chemical_lamda - v1 - v2 - v3 - v4;

	}
      }else{
	layer4[p][q][r+1].CB02++;
	v3 = layer4[p][q][r+1].diffusion_lamda;
	v4 = layer4[p][q][r+1].chemical_lamda;
	layer4[p][q][r+1].diffusion_lamda = layer4_diffusion_lamda(p, q, r+1);
	layer4[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r+1]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r+1].diffusion_lamda +
	  layer4[p][q][r+1].chemical_lamda - v1 - v2 - v3 - v4;
      }
      goto loop4;
    }

	choose1 = choose;
    choose += layer4[p][q][r].CB21 * layer4[p][q][r].D_CB21[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer4[p][q][r].CB21--;
      v1 = layer4[p][q][r].diffusion_lamda;
      v2 = layer4[p][q][r].chemical_lamda;
      layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
      layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);

      if(k == 0){
	layer4[p-1][q][r].CB21++;
	v3 = layer4[p-1][q][r].diffusion_lamda;
	v4 = layer4[p-1][q][r].chemical_lamda;
	layer4[p-1][q][r].diffusion_lamda = layer4_diffusion_lamda(p-1, q, r);
	layer4[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p-1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p-1][q][r].diffusion_lamda +
	  layer4[p-1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 1){
	layer4[p][q-1][r].CB21++;
	v3 = layer4[p][q-1][r].diffusion_lamda;
	v4 = layer4[p][q-1][r].chemical_lamda;
	layer4[p][q-1][r].diffusion_lamda = layer4_diffusion_lamda(p, q-1, r);
	layer4[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q-1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q-1][r].diffusion_lamda +
	  layer4[p][q-1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 2){
	layer4[p+1][q][r].CB21++;
	v3 = layer4[p+1][q][r].diffusion_lamda;
	v4 = layer4[p+1][q][r].chemical_lamda;
	layer4[p+1][q][r].diffusion_lamda = layer4_diffusion_lamda(p+1, q, r);
	layer4[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p+1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p+1][q][r].diffusion_lamda +
	  layer4[p+1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 3){
	layer4[p][q+1][r].CB21++;
	v3 = layer4[p][q+1][r].diffusion_lamda;
	v4 = layer4[p][q+1][r].chemical_lamda;
	layer4[p][q+1][r].diffusion_lamda = layer4_diffusion_lamda(p, q+1, r);
	layer4[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q+1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q+1][r].diffusion_lamda +
	  layer4[p][q+1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 4){
	if((p == 2) && (q == 0) && (r == 0)){
	  layer3[7].CB21++;
	  v3 = layer3[7].diffusion_lamda;
	  v4 = layer3[7].chemical_lamda;
	  layer3[7].diffusion_lamda = layer3_diffusion_lamda(7);
	  layer3[7].chemical_lamda = voxel_chemical_lamda(layer3[7]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer3[7].diffusion_lamda +
	    layer3[7].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer4[p][q][r-1].CB21++;
	  v3 = layer4[p][q][r-1].diffusion_lamda;
	  v4 = layer4[p][q][r-1].chemical_lamda;
	  layer4[p][q][r-1].diffusion_lamda = layer4_diffusion_lamda(p, q, r-1);
	  layer4[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r-1]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r-1].diffusion_lamda +
	    layer4[p][q][r-1].chemical_lamda - v1 - v2 - v3 - v4;

	}
      }else{
	layer4[p][q][r+1].CB21++;
	v3 = layer4[p][q][r+1].diffusion_lamda;
	v4 = layer4[p][q][r+1].chemical_lamda;
	layer4[p][q][r+1].diffusion_lamda = layer4_diffusion_lamda(p, q, r+1);
	layer4[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r+1]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r+1].diffusion_lamda +
	  layer4[p][q][r+1].chemical_lamda - v1 - v2 - v3 - v4;
      }
      goto loop4;
    }

	choose1 = choose;
    choose += layer4[p][q][r].CB12 * layer4[p][q][r].D_CB12[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer4[p][q][r].CB12--;
      v1 = layer4[p][q][r].diffusion_lamda;
      v2 = layer4[p][q][r].chemical_lamda;
      layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
      layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);

      if(k == 0){
	layer4[p-1][q][r].CB12++;
	v3 = layer4[p-1][q][r].diffusion_lamda;
	v4 = layer4[p-1][q][r].chemical_lamda;
	layer4[p-1][q][r].diffusion_lamda = layer4_diffusion_lamda(p-1, q, r);
	layer4[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p-1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p-1][q][r].diffusion_lamda +
	  layer4[p-1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 1){
	layer4[p][q-1][r].CB12++;
	v3 = layer4[p][q-1][r].diffusion_lamda;
	v4 = layer4[p][q-1][r].chemical_lamda;
	layer4[p][q-1][r].diffusion_lamda = layer4_diffusion_lamda(p, q-1, r);
	layer4[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q-1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q-1][r].diffusion_lamda +
	  layer4[p][q-1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 2){
	layer4[p+1][q][r].CB12++;
	v3 = layer4[p+1][q][r].diffusion_lamda;
	v4 = layer4[p+1][q][r].chemical_lamda;
	layer4[p+1][q][r].diffusion_lamda = layer4_diffusion_lamda(p+1, q, r);
	layer4[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p+1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p+1][q][r].diffusion_lamda +
	  layer4[p+1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 3){
	layer4[p][q+1][r].CB12++;
	v3 = layer4[p][q+1][r].diffusion_lamda;
	v4 = layer4[p][q+1][r].chemical_lamda;
	layer4[p][q+1][r].diffusion_lamda = layer4_diffusion_lamda(p, q+1, r);
	layer4[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q+1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q+1][r].diffusion_lamda +
	  layer4[p][q+1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 4){
	if((p == 2) && (q == 0) && (r == 0)){
	  layer3[7].CB12++;
	  v3 = layer3[7].diffusion_lamda;
	  v4 = layer3[7].chemical_lamda;
	  layer3[7].diffusion_lamda = layer3_diffusion_lamda(7);
	  layer3[7].chemical_lamda = voxel_chemical_lamda(layer3[7]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer3[7].diffusion_lamda +
	    layer3[7].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer4[p][q][r-1].CB12++;
	  v3 = layer4[p][q][r-1].diffusion_lamda;
	  v4 = layer4[p][q][r-1].chemical_lamda;
	  layer4[p][q][r-1].diffusion_lamda = layer4_diffusion_lamda(p, q, r-1);
	  layer4[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r-1]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r-1].diffusion_lamda +
	    layer4[p][q][r-1].chemical_lamda - v1 - v2 - v3 - v4;

	}
      }else{
	layer4[p][q][r+1].CB12++;
	v3 = layer4[p][q][r+1].diffusion_lamda;
	v4 = layer4[p][q][r+1].chemical_lamda;
	layer4[p][q][r+1].diffusion_lamda = layer4_diffusion_lamda(p, q, r+1);
	layer4[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r+1]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r+1].diffusion_lamda +
	  layer4[p][q][r+1].chemical_lamda - v1 - v2 - v3 - v4;
      }
      goto loop4;
    }

	choose1 = choose;
    choose += layer4[p][q][r].CB22 * layer4[p][q][r].D_CB22[k];
    if((choose > random_num) && (choose1 <= random_num)){
      layer4[p][q][r].CB22--;
      v1 = layer4[p][q][r].diffusion_lamda;
      v2 = layer4[p][q][r].chemical_lamda;
      layer4[p][q][r].diffusion_lamda = layer4_diffusion_lamda(p, q, r);
      layer4[p][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r]);

      if(k == 0){
	layer4[p-1][q][r].CB22++;
	v3 = layer4[p-1][q][r].diffusion_lamda;
	v4 = layer4[p-1][q][r].chemical_lamda;
	layer4[p-1][q][r].diffusion_lamda = layer4_diffusion_lamda(p-1, q, r);
	layer4[p-1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p-1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p-1][q][r].diffusion_lamda +
	  layer4[p-1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 1){
	layer4[p][q-1][r].CB22++;
	v3 = layer4[p][q-1][r].diffusion_lamda;
	v4 = layer4[p][q-1][r].chemical_lamda;
	layer4[p][q-1][r].diffusion_lamda = layer4_diffusion_lamda(p, q-1, r);
	layer4[p][q-1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q-1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q-1][r].diffusion_lamda +
	  layer4[p][q-1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 2){
	layer4[p+1][q][r].CB22++;
	v3 = layer4[p+1][q][r].diffusion_lamda;
	v4 = layer4[p+1][q][r].chemical_lamda;
	layer4[p+1][q][r].diffusion_lamda = layer4_diffusion_lamda(p+1, q, r);
	layer4[p+1][q][r].chemical_lamda = voxel_chemical_lamda(layer4[p+1][q][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p+1][q][r].diffusion_lamda +
	  layer4[p+1][q][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 3){
	layer4[p][q+1][r].CB22++;
	v3 = layer4[p][q+1][r].diffusion_lamda;
	v4 = layer4[p][q+1][r].chemical_lamda;
	layer4[p][q+1][r].diffusion_lamda = layer4_diffusion_lamda(p, q+1, r);
	layer4[p][q+1][r].chemical_lamda = voxel_chemical_lamda(layer4[p][q+1][r]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q+1][r].diffusion_lamda +
	  layer4[p][q+1][r].chemical_lamda - v1 - v2 - v3 - v4;

      }else if(k == 4){
	if((p == 2) && (q == 0) && (r == 0)){
	  layer3[7].CB22++;
	  v3 = layer3[7].diffusion_lamda;
	  v4 = layer3[7].chemical_lamda;
	  layer3[7].diffusion_lamda = layer3_diffusion_lamda(7);
	  layer3[7].chemical_lamda = voxel_chemical_lamda(layer3[7]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer3[7].diffusion_lamda +
	    layer3[7].chemical_lamda - v1 - v2 - v3 - v4;

	}else{
	  layer4[p][q][r-1].CB22++;
	  v3 = layer4[p][q][r-1].diffusion_lamda;
	  v4 = layer4[p][q][r-1].chemical_lamda;
	  layer4[p][q][r-1].diffusion_lamda = layer4_diffusion_lamda(p, q, r-1);
	  layer4[p][q][r-1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r-1]);
	  reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r-1].diffusion_lamda +
	    layer4[p][q][r-1].chemical_lamda - v1 - v2 - v3 - v4;

	}
      }else{
	layer4[p][q][r+1].CB22++;
	v3 = layer4[p][q][r+1].diffusion_lamda;
	v4 = layer4[p][q][r+1].chemical_lamda;
	layer4[p][q][r+1].diffusion_lamda = layer4_diffusion_lamda(p, q, r+1);
	layer4[p][q][r+1].chemical_lamda = voxel_chemical_lamda(layer4[p][q][r+1]);
	reaction_lamda += layer4[p][q][r].diffusion_lamda + layer4[p][q][r].chemical_lamda + layer4[p][q][r+1].diffusion_lamda +
	  layer4[p][q][r+1].chemical_lamda - v1 - v2 - v3 - v4;
      }
      goto loop4;
    }
  }

  x = p;
  y = q;
  z = r;
  chemical_reaction(&layer4[p][q][r], choose, 4);
  
 loop4: k = 0;
}

void chemical_reaction(struct voxel *point, double choose2, int index){
  int i, j;
  double choose, choose1, v1, v2;
  struct CaMKII *p1;

  choose = choose2;

  choose1 = choose;
  choose += point->CaM00 * point->Ca * point->CaM00_to_CaM10;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM00--;
    point->Ca--;
    point->CaM10++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM10 * point->CaM10_to_CaM00;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM00++;
    point->Ca++;
    point->CaM10--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM00 * point->Ca * point->CaM00_to_CaM01;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM00--;
    point->Ca--;
    point->CaM01++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM01 * point->CaM01_to_CaM00;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM00++;
    point->Ca++;
    point->CaM01--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM10 * point->Ca * point->CaM10_to_CaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM10--;
    point->Ca--;
    point->CaM11++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM11 * point->CaM11_to_CaM10;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM10++;
    point->Ca++;
    point->CaM11--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM10 * point->Ca * point->CaM10_to_CaM20;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM10--;
    point->Ca--;
    point->CaM20++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM20 * point->CaM20_to_CaM10;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM10++;
    point->Ca++;
    point->CaM20--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM01 * point->Ca * point->CaM01_to_CaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM01--;
    point->Ca--;
    point->CaM11++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM11 * point->CaM11_to_CaM01;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM01++;
    point->Ca++;
    point->CaM11--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM01 * point->Ca * point->CaM01_to_CaM02;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM01--;
    point->Ca--;
    point->CaM02++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM02 * point->CaM02_to_CaM01;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM01++;
    point->Ca++;
    point->CaM02--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM11 * point->Ca * point->CaM11_to_CaM21;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM11--;
    point->Ca--;
    point->CaM21++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM21 * point->CaM21_to_CaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM11++;
    point->Ca++;
    point->CaM21--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM11 * point->Ca * point->CaM11_to_CaM12;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM11--;
    point->Ca--;
    point->CaM12++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM12 * point->CaM12_to_CaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM11++;
    point->Ca++;
    point->CaM12--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM20 * point->Ca * point->CaM20_to_CaM21;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM20--;
    point->Ca--;
    point->CaM21++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM21 * point->CaM21_to_CaM20;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM20++;
    point->Ca++;
    point->CaM21--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM02 * point->Ca * point->CaM02_to_CaM12;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM02--;
    point->Ca--;
    point->CaM12++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM12 * point->CaM12_to_CaM02;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM02++;
    point->Ca++;
    point->CaM12--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM21 * point->Ca * point->CaM21_to_CaM22;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM21--;
    point->Ca--;
    point->CaM22++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM22 * point->CaM22_to_CaM21;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM21++;
    point->Ca++;
    point->CaM22--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM12 * point->Ca * point->CaM12_to_CaM22;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM12--;
    point->Ca--;
    point->CaM22++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaM22 * point->CaM22_to_CaM12;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaM12++;
    point->Ca++;
    point->CaM22--;
    goto loop5;
  }
  
/*******************************************************************/ 

  choose1 = choose;
  choose += point->CaN * point->CaM00 * point->CaM00_to_CaNCaM00;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaN--;
    point->CaM00--;
    point->CaNCaM00++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM00 * point->CaNCaM00_to_CaM00;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaN++;
    point->CaM00++;
    point->CaNCaM00--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaN * point->CaM10 * point->CaM10_to_CaNCaM10;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaN--;
    point->CaM10--;
    point->CaNCaM10++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM10 * point->CaNCaM10_to_CaM10;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaN++;
    point->CaM10++;
    point->CaNCaM10--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaN * point->CaM01 * point->CaM01_to_CaNCaM01;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaN--;
    point->CaM01--;
    point->CaNCaM01++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM01 * point->CaNCaM01_to_CaM01;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaN++;
    point->CaM01++;
    point->CaNCaM01--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaN * point->CaM11 * point->CaM11_to_CaNCaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaN--;
    point->CaM11--;
    point->CaNCaM11++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM11 * point->CaNCaM11_to_CaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaN++;
    point->CaM11++;
    point->CaNCaM11--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaN * point->CaM20 * point->CaM20_to_CaNCaM20;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaN--;
    point->CaM20--;
    point->CaNCaM20++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM20 * point->CaNCaM20_to_CaM20;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaN++;
    point->CaM20++;
    point->CaNCaM20--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaN * point->CaM02 * point->CaM02_to_CaNCaM02;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaN--;
    point->CaM02--;
    point->CaNCaM02++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM02 * point->CaNCaM02_to_CaM02;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaN++;
    point->CaM02++;
    point->CaNCaM02--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaN * point->CaM12 * point->CaM12_to_CaNCaM12;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaN--;
    point->CaM12--;
    point->CaNCaM12++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM12 * point->CaNCaM12_to_CaM12;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaN++;
    point->CaM12++;
    point->CaNCaM12--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaN * point->CaM21 * point->CaM21_to_CaNCaM21;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaN--;
    point->CaM21--;
    point->CaNCaM21++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM21 * point->CaNCaM21_to_CaM21;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaN++;
    point->CaM21++;
    point->CaNCaM21--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaN * point->CaM22 * point->CaM22_to_CaNCaM22;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaN--;
    point->CaM22--;
    point->CaNCaM22++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM22 * point->CaNCaM22_to_CaM22;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaN++;
    point->CaM22++;
    point->CaNCaM22--;
    goto loop5;
  }

/*********************************************************************/

  choose1 = choose;
  choose += point->CaNCaM00 * point->Ca * point->CaNCaM00_to_CaNCaM10;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaNCaM00--;
    point->Ca--;
    point->CaNCaM10++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM10 * point->CaNCaM10_to_CaNCaM00;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->CaNCaM00++;
    point->CaNCaM10--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM00 * point->Ca * point->CaNCaM00_to_CaNCaM01;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaNCaM00--;
    point->Ca--;
    point->CaNCaM01++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM01 * point->CaNCaM01_to_CaNCaM00;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->CaNCaM00++;
    point->CaNCaM01--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM10 * point->Ca * point->CaNCaM10_to_CaNCaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaNCaM10--;
    point->Ca--;
    point->CaNCaM11++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM11 * point->CaNCaM11_to_CaNCaM10;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->CaNCaM10++;
    point->CaNCaM11--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM10 * point->Ca * point->CaNCaM10_to_CaNCaM20;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaNCaM10--;
    point->Ca--;
    point->CaNCaM20++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM20 * point->CaNCaM20_to_CaNCaM10;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->CaNCaM10++;
    point->CaNCaM20--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM01 * point->Ca * point->CaNCaM01_to_CaNCaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaNCaM01--;
    point->Ca--;
    point->CaNCaM11++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM11 * point->CaNCaM11_to_CaNCaM01;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->CaNCaM01++;
    point->CaNCaM11--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM01 * point->Ca * point->CaNCaM01_to_CaNCaM02;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaNCaM01--;
    point->Ca--;
    point->CaNCaM02++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM02 * point->CaNCaM02_to_CaNCaM01;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->CaNCaM01++;
    point->CaNCaM02--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM11 * point->Ca * point->CaNCaM11_to_CaNCaM21;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaNCaM11--;
    point->Ca--;
    point->CaNCaM21++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM21 * point->CaNCaM21_to_CaNCaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->CaNCaM11++;
    point->CaNCaM21--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM11 * point->Ca * point->CaNCaM11_to_CaNCaM12;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaNCaM11--;
    point->Ca--;
    point->CaNCaM12++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM12 * point->CaNCaM12_to_CaNCaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->CaNCaM11++;
    point->CaNCaM12--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM20 * point->Ca * point->CaNCaM20_to_CaNCaM21;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaNCaM20--;
    point->Ca--;
    point->CaNCaM21++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM21 * point->CaNCaM21_to_CaNCaM20;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->CaNCaM20++;
    point->CaNCaM21--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM02 * point->Ca * point->CaNCaM02_to_CaNCaM12;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaNCaM02--;
    point->Ca--;
    point->CaNCaM12++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM12 * point->CaNCaM12_to_CaNCaM02;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->CaNCaM02++;
    point->CaNCaM12--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM21 * point->Ca * point->CaNCaM21_to_CaNCaM22;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaNCaM21--;
    point->Ca--;
    point->CaNCaM22++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM22 * point->CaNCaM22_to_CaNCaM21;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->CaNCaM21++;
    point->CaNCaM22--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM12 * point->Ca * point->CaNCaM12_to_CaNCaM22;
  if((choose > random_num) && (choose1 <= random_num)){
    point->CaNCaM12--;
    point->Ca--;
    point->CaNCaM22++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->CaNCaM22 * point->CaNCaM22_to_CaNCaM12;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->CaNCaM12++;
    point->CaNCaM22--;
    goto loop5;
  }

/*******************************************************************/ 

  choose1 = choose;
  choose += point->Ng * point->CaM00 * point->CaM00_to_NgCaM00;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ng--;
    point->CaM00--;
    point->NgCaM00++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM00 * point->NgCaM00_to_CaM00;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ng++;
    point->CaM00++;
    point->NgCaM00--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->Ng * point->CaM10 * point->CaM10_to_NgCaM10;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ng--;
    point->CaM10--;
    point->NgCaM10++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM10 * point->NgCaM10_to_CaM10;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ng++;
    point->CaM10++;
    point->NgCaM10--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->Ng * point->CaM01 * point->CaM01_to_NgCaM01;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ng--;
    point->CaM01--;
    point->NgCaM01++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM01 * point->NgCaM01_to_CaM01;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ng++;
    point->CaM01++;
    point->NgCaM01--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->Ng * point->CaM11 * point->CaM11_to_NgCaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ng--;
    point->CaM11--;
    point->NgCaM11++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM11 * point->NgCaM11_to_CaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ng++;
    point->CaM11++;
    point->NgCaM11--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->Ng * point->CaM20 * point->CaM20_to_NgCaM20;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ng--;
    point->CaM20--;
    point->NgCaM20++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM20 * point->NgCaM20_to_CaM20;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ng++;
    point->CaM20++;
    point->NgCaM20--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->Ng * point->CaM02 * point->CaM02_to_NgCaM02;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ng--;
    point->CaM02--;
    point->NgCaM02++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM02 * point->NgCaM02_to_CaM02;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ng++;
    point->CaM02++;
    point->NgCaM02--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->Ng * point->CaM12 * point->CaM12_to_NgCaM12;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ng--;
    point->CaM12--;
    point->NgCaM12++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM12 * point->NgCaM12_to_CaM12;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ng++;
    point->CaM12++;
    point->NgCaM12--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->Ng * point->CaM21 * point->CaM21_to_NgCaM21;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ng--;
    point->CaM21--;
    point->NgCaM21++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM21 * point->NgCaM21_to_CaM21;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ng++;
    point->CaM21++;
    point->NgCaM21--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->Ng * point->CaM22 * point->CaM22_to_NgCaM22;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ng--;
    point->CaM22--;
    point->NgCaM22++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM22 * point->NgCaM22_to_CaM22;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ng++;
    point->CaM22++;
    point->NgCaM22--;
    goto loop5;
  }

/*********************************************************************/

  choose1 = choose;
  choose += point->NgCaM00 * point->Ca * point->NgCaM00_to_NgCaM10;
  if((choose > random_num) && (choose1 <= random_num)){
    point->NgCaM00--;
    point->Ca--;
    point->NgCaM10++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM10 * point->NgCaM10_to_NgCaM00;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->NgCaM00++;
    point->NgCaM10--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM00 * point->Ca * point->NgCaM00_to_NgCaM01;
  if((choose > random_num) && (choose1 <= random_num)){
    point->NgCaM00--;
    point->Ca--;
    point->NgCaM01++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM01 * point->NgCaM01_to_NgCaM00;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->NgCaM00++;
    point->NgCaM01--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM10 * point->Ca * point->NgCaM10_to_NgCaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->NgCaM10--;
    point->Ca--;
    point->NgCaM11++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM11 * point->NgCaM11_to_NgCaM10;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->NgCaM10++;
    point->NgCaM11--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM10 * point->Ca * point->NgCaM10_to_NgCaM20;
  if((choose > random_num) && (choose1 <= random_num)){
    point->NgCaM10--;
    point->Ca--;
    point->NgCaM20++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM20 * point->NgCaM20_to_NgCaM10;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->NgCaM10++;
    point->NgCaM20--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM01 * point->Ca * point->NgCaM01_to_NgCaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->NgCaM01--;
    point->Ca--;
    point->NgCaM11++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM11 * point->NgCaM11_to_NgCaM01;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->NgCaM01++;
    point->NgCaM11--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM01 * point->Ca * point->NgCaM01_to_NgCaM02;
  if((choose > random_num) && (choose1 <= random_num)){
    point->NgCaM01--;
    point->Ca--;
    point->NgCaM02++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM02 * point->NgCaM02_to_NgCaM01;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->NgCaM01++;
    point->NgCaM02--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM11 * point->Ca * point->NgCaM11_to_NgCaM21;
  if((choose > random_num) && (choose1 <= random_num)){
    point->NgCaM11--;
    point->Ca--;
    point->NgCaM21++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM21 * point->NgCaM21_to_NgCaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->NgCaM11++;
    point->NgCaM21--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM11 * point->Ca * point->NgCaM11_to_NgCaM12;
  if((choose > random_num) && (choose1 <= random_num)){
    point->NgCaM11--;
    point->Ca--;
    point->NgCaM12++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM12 * point->NgCaM12_to_NgCaM11;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->NgCaM11++;
    point->NgCaM12--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM20 * point->Ca * point->NgCaM20_to_NgCaM21;
  if((choose > random_num) && (choose1 <= random_num)){
    point->NgCaM20--;
    point->Ca--;
    point->NgCaM21++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM21 * point->NgCaM21_to_NgCaM20;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->NgCaM20++;
    point->NgCaM21--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM02 * point->Ca * point->NgCaM02_to_NgCaM12;
  if((choose > random_num) && (choose1 <= random_num)){
    point->NgCaM02--;
    point->Ca--;
    point->NgCaM12++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM12 * point->NgCaM12_to_NgCaM02;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->NgCaM02++;
    point->NgCaM12--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM21 * point->Ca * point->NgCaM21_to_NgCaM22;
  if((choose > random_num) && (choose1 <= random_num)){
    point->NgCaM21--;
    point->Ca--;
    point->NgCaM22++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM22 * point->NgCaM22_to_NgCaM21;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->NgCaM21++;
    point->NgCaM22--;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM12 * point->Ca * point->NgCaM12_to_NgCaM22;
  if((choose > random_num) && (choose1 <= random_num)){
    point->NgCaM12--;
    point->Ca--;
    point->NgCaM22++;
    goto loop5;
  }

  choose1 = choose;
  choose += point->NgCaM22 * point->NgCaM22_to_NgCaM12;
  if((choose > random_num) && (choose1 <= random_num)){
    point->Ca++;
    point->NgCaM12++;
    point->NgCaM22--;
    goto loop5;
  }

/*****************************************************************/


  choose1 = choose;
  choose += point->CB00 * point->Ca * point->CB00_to_CB10;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB00--;
	  point->Ca--;
	  point->CB10++;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB10 * point->CB10_to_CB00;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB00++;
	  point->Ca++;
	  point->CB10--;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB00 * point->Ca * point->CB00_to_CB01;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB00--;
	  point->Ca--;
	  point->CB01++;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB01 * point->CB01_to_CB00;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB00++;
	  point->Ca++;
	  point->CB01--;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB10 * point->Ca * point->CB10_to_CB11;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB10--;
	  point->Ca--;
	  point->CB11++;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB11 * point->CB11_to_CB10;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB10++;
	  point->Ca++;
	  point->CB11--;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB10 * point->Ca * point->CB10_to_CB20;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB10--;
	  point->Ca--;
	  point->CB20++;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB20 * point->CB20_to_CB10;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB10++;
	  point->Ca++;
	  point->CB20--;
	  goto loop5;
  }
  
  choose1 = choose;
  choose += point->CB01 * point->Ca * point->CB01_to_CB11;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB01--;
	  point->Ca--;
	  point->CB11++;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB11 * point->CB11_to_CB01;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB01++;
	  point->Ca++;
	  point->CB11--;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB01 * point->Ca * point->CB01_to_CB02;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB01--;
	  point->Ca--;
	  point->CB02++;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB02 * point->CB02_to_CB01;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB01++;
	  point->Ca++;
	  point->CB02--;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB11 * point->Ca * point->CB11_to_CB21;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB11--;
	  point->Ca--;
	  point->CB21++;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB21 * point->CB21_to_CB11;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB11++;
	  point->Ca++;
	  point->CB21--;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB11 * point->Ca * point->CB11_to_CB12;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB11--;
	  point->Ca--;
	  point->CB12++;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB12 * point->CB12_to_CB11;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB11++;
	  point->Ca++;
	  point->CB12--;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB20 * point->Ca * point->CB20_to_CB21;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB20--;
	  point->Ca--;
	  point->CB21++;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB21 * point->CB21_to_CB20;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB20++;
	  point->Ca++;
	  point->CB21--;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB02 * point->Ca * point->CB02_to_CB12;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB02--;
	  point->Ca--;
	  point->CB12++;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB12 * point->CB12_to_CB02;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB02++;
	  point->Ca++;
	  point->CB12--;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB12 * point->Ca * point->CB12_to_CB22;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB12--;
	  point->Ca--;
	  point->CB22++;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB22 * point->CB22_to_CB12;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB12++;
	  point->Ca++;
	  point->CB22--;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB21 * point->Ca * point->CB21_to_CB22;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB21--;
	  point->Ca--;
	  point->CB22++;
	  goto loop5;
  }

  choose1 = choose;
  choose += point->CB22 * point->CB22_to_CB21;
  if((choose > random_num) && (choose1 <= random_num)){
	  point->CB21++;
	  point->Ca++;
	  point->CB22--;
	  goto loop5;
  }

/*****************************************************************/

  p1 = point->head;
  while(p1 != NULL){
    for(i=0; i<2; i++){
      for(j=0; j<6; j++){

	if(p1->subunit[i][j] ==0){
		
	  choose1 = choose;
	  choose += point->CaM00 * point->CaM00_to_CaMKIICaM00;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM00--;
	    p1->subunit[i][j] = 1;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM10 * point->CaM10_to_CaMKIICaM10;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM10--;
	    p1->subunit[i][j] = 2;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM01 * point->CaM01_to_CaMKIICaM01;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM01--;
	    p1->subunit[i][j] = 3;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM11 * point->CaM11_to_CaMKIICaM11;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM11--;
	    p1->subunit[i][j] = 4;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM20 * point->CaM20_to_CaMKIICaM20;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM20--;
	    p1->subunit[i][j] = 5;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM02 * point->CaM02_to_CaMKIICaM02;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM02--;
	    p1->subunit[i][j] = 6;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM21 * point->CaM21_to_CaMKIICaM21;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM21--;
	    p1->subunit[i][j] = 7;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM12 * point->CaM12_to_CaMKIICaM12;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM12--;
	    p1->subunit[i][j] = 8;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM22 * point->CaM22_to_CaMKIICaM22;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM22--;
	    p1->subunit[i][j] = 9;
	    goto loop5;
	  }
 
	}else if(p1->subunit[i][j] == 1){

	  choose1 = choose;
	  choose += point->CaMKIICaM00_to_CaM00;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM00++;
	    p1->subunit[i][j] = 0;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM00_to_CaMKIICaM10;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 2;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM00_to_CaMKIICaM01;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 3;
	    goto loop5;
	  }
/*
	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))
     	 choose += point->CaMKIICaM00_to_Trapped00;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 10;
	    goto loop5;
	  }

	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
         (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
         choose += 0.4 * point->CaMKIICaM00_to_Trapped00;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 10;
	    goto loop5;
	  }
*/
	}else if(p1->subunit[i][j] == 2){

	  choose1 = choose;
	  choose += point->CaMKIICaM10_to_CaMKIICaM00;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 1;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaMKIICaM10_to_CaM10;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM10++;
	    p1->subunit[i][j] = 0;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM10_to_CaMKIICaM11;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 4;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM10_to_CaMKIICaM20;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 5;
	    goto loop5;
	  }
/*
	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))
     	 choose += point->CaMKIICaM10_to_Trapped10;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 11;
	    goto loop5;
	  }

	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
         (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	     choose += 0.4 * point->CaMKIICaM10_to_Trapped10;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 11;
	    goto loop5;
	  }
*/
	}else if(p1->subunit[i][j] == 3){

	  choose1 = choose;
	  choose += point->CaMKIICaM01_to_CaMKIICaM00;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 1;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaMKIICaM01_to_CaM01;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM01++;
	    p1->subunit[i][j] = 0;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM01_to_CaMKIICaM11;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 4;
	    goto loop5;
	  }

      choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM01_to_CaMKIICaM02;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 6;
	    goto loop5;
	  }
/*
	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))	
     	 choose += point->CaMKIICaM01_to_Trapped01;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 12;
	    goto loop5;
	  }

	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
         (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	     choose += 0.4 * point->CaMKIICaM01_to_Trapped01;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 12;
	    goto loop5;
	  }
*/
	}else if(p1->subunit[i][j] == 4){
	  
	  choose1 = choose;
	  choose += point->CaMKIICaM11_to_CaMKIICaM01;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 3;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaMKIICaM11_to_CaMKIICaM10;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 2;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose1 = point->CaMKIICaM11_to_CaM11;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM11++;
	    p1->subunit[i][j] = 0;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose1 = point->Ca * point->CaMKIICaM11_to_CaMKIICaM21;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 7;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose1 = point->Ca * point->CaMKIICaM11_to_CaMKIICaM12;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 8;
	    goto loop5;
	  }

	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))		 
         choose += point->CaMKIICaM11_to_Trapped11;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 13;
	    goto loop5;
	  }

	  choose1 = choose;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	     choose += 0.4 * point->CaMKIICaM11_to_Trapped11;
*/
//
	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	     choose += 0.4 * point->CaMKIICaM11_to_Trapped11;
//
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 13;
	    goto loop5;
	  }
	  
	}else if(p1->subunit[i][j] == 5){

	  choose1 = choose;
	  choose += point->CaMKIICaM20_to_CaMKIICaM10;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 2;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM20_to_CaMKIICaM21;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 7;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaMKIICaM20_to_CaM20;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM20++;
	    p1->subunit[i][j] = 0;
	    goto loop5;
	  }

	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))		 
        choose += point->CaMKIICaM20_to_Trapped20;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 14;
	    goto loop5;
	  }

	  choose1 = choose;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM20_to_Trapped20;
*/
//
	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM20_to_Trapped20;
//
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 14;
	    goto loop5;
	  }

	}else if(p1->subunit[i][j] == 6){
	
	  choose1 = choose;
	  choose += point->CaMKIICaM02_to_CaMKIICaM01;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 3;
	  point->Ca++;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM02_to_CaMKIICaM12;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 8;
	  point->Ca--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaMKIICaM02_to_CaM02;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 0;
	  point->CaM02++;
	  goto loop5;
	  }

	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))		 
        choose += point->CaMKIICaM02_to_Trapped02;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 15;
	    goto loop5;
	  }

	  choose1 = choose;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM02_to_Trapped02;
*/
//
	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM02_to_Trapped02;
//
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 15;
	    goto loop5;
	  }
	  
	}else if(p1->subunit[i][j] == 7){
		
	  choose1 = choose;
	  choose += point->CaMKIICaM21_to_CaMKIICaM11;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 4;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaMKIICaM21_to_CaMKIICaM20;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 5;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM21_to_CaMKIICaM22;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 9;
	    goto loop5;
	  }

      choose1 = choose;
	  choose += point->CaMKIICaM21_to_CaM21;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM21++;
	    p1->subunit[i][j] = 0;
	    goto loop5;
	  }

	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))		 
        choose += point->CaMKIICaM21_to_Trapped21;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 16;
	    goto loop5;
	  }

	  choose1 = choose;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM21_to_Trapped21;
*/
//
	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM21_to_Trapped21;
//
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 16;
	    goto loop5;
	  }
	}else if(p1->subunit[i][j] == 8){

	  choose1 = choose;
	  choose += point->CaMKIICaM12_to_CaMKIICaM11;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 4;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaMKIICaM12_to_CaMKIICaM02;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 6;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->CaMKIICaM12_to_CaMKIICaM22;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca--;
	    p1->subunit[i][j] = 9;
	    goto loop5;
	  }

      choose1 = choose;
	  choose += point->CaMKIICaM12_to_CaM12;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM12++;
	    p1->subunit[i][j] = 0;
	    goto loop5;
	  }

	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))		 
        choose += point->CaMKIICaM12_to_Trapped12;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 17;
	    goto loop5;
	  }

	  choose1 = choose;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM12_to_Trapped12;
*/
//
	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM12_to_Trapped12;
//
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 17;
	    goto loop5;
	  }
	}else if(p1->subunit[i][j] == 9){
      
	  choose1 = choose;
	  choose += point->CaMKIICaM22_to_CaMKIICaM21;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 7;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaMKIICaM22_to_CaMKIICaM12;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->Ca++;
	    p1->subunit[i][j] = 8;
	    goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaMKIICaM22_to_CaM22;
	  if((choose > random_num) && (choose1 <= random_num)){
	    point->CaM22++;
	    p1->subunit[i][j] = 0;
	    goto loop5;
	  }

	  choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))		 
        choose += point->CaMKIICaM22_to_Trapped22;
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 18;
	    goto loop5;
	  }

	  choose1 = choose;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM22_to_Trapped22;
*/
//
	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->CaMKIICaM22_to_Trapped22;
//
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 18;
	    goto loop5;
	  }
	}else if(p1->subunit[i][j] == 10){

      choose1 = choose;
	  choose += point->Trapped00_to_CaMKIICaM00;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 1;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->Trapped00_to_Trapped01;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 12;
		 point->Ca--;
	     goto loop5;
	   }

	  choose1 = choose;
	  choose += point->Ca * point->Trapped00_to_Trapped10;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 11;
		 point->Ca--;
	     goto loop5;
	   }

	  choose1 = choose;
	  choose += point->Trapped00_to_Auton;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 19;
		 point->CaM00++;
	     goto loop5;
	   }

      }else if(p1->subunit[i][j] == 11){
      choose1 = choose;
	  choose += point->Trapped10_to_CaMKIICaM10;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 2;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped10_to_Trapped00;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 10;
	  goto loop5;
	  }

      choose1 = choose;
	  choose += point->Ca * point->Trapped10_to_Trapped11;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 13;
		 point->Ca--;
	     goto loop5;
	   }

	  choose1 = choose;
	  choose += point->Ca * point->Trapped10_to_Trapped20;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 14;
		 point->Ca--;
	     goto loop5;
	   }

      choose1 = choose;
	  choose += point->Trapped10_to_Auton;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 19;
		 point->CaM10++;
	     goto loop5;
	   }

	  }else if(p1->subunit[i][j] == 12){
      choose1 = choose;
	  choose += point->Trapped01_to_CaMKIICaM01;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 3;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped01_to_Trapped00;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 10;
	  goto loop5;
	  }

      choose1 = choose;
	  choose += point->Ca * point->Trapped01_to_Trapped11;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 13;
		 point->Ca--;
	     goto loop5;
	   }

	  choose1 = choose;
	  choose += point->Ca * point->Trapped01_to_Trapped02;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 15;
		 point->Ca--;
	     goto loop5;
	   }

      choose1 = choose;
	  choose += point->Trapped01_to_Auton;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 19;
		 point->CaM01++;
	     goto loop5;
	   }

	  }else if(p1->subunit[i][j] == 13){
      choose1 = choose;
	  choose += point->Trapped11_to_CaMKIICaM11;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 4;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped11_to_Trapped10;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 11;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped11_to_Trapped01;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 12;
	  goto loop5;
	  }

      choose1 = choose;
	  choose += point->Ca * point->Trapped11_to_Trapped21;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 16;
		 point->Ca--;
	     goto loop5;
	   }

	  choose1 = choose;
	  choose += point->Ca * point->Trapped11_to_Trapped12;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 17;
		 point->Ca--;
	     goto loop5;
	   }

	  choose1 = choose;
	  choose += point->Trapped11_to_Auton;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 19;
		 point->CaM11++;
	     goto loop5;
	   }

	  }else if(p1->subunit[i][j] == 14){
      choose1 = choose;
	  choose += point->Trapped20_to_CaMKIICaM20;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 5;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped20_to_Trapped10;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 11;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->Trapped20_to_Trapped21;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 16;
		 point->Ca--;
	     goto loop5;
	   }
	  
	  choose1 = choose;
	  choose += point->Trapped20_to_Auton;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 19;
		 point->CaM20++;
	     goto loop5;
	   }

	  }else if(p1->subunit[i][j] == 15){
      choose1 = choose;
	  choose += point->Trapped02_to_CaMKIICaM02;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 6;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped02_to_Trapped01;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 12;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->Trapped02_to_Trapped12;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 17;
		 point->Ca--;
	     goto loop5;
	   }
	  
	  choose1 = choose;
	  choose += point->Trapped02_to_Auton;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 19;
		 point->CaM02++;
	     goto loop5;
	   }

	  }else if(p1->subunit[i][j] == 16){
      choose1 = choose;
	  choose += point->Trapped21_to_CaMKIICaM21;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 7;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped21_to_Trapped20;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 14;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped21_to_Trapped11;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 13;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->Trapped21_to_Trapped22;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 18;
		 point->Ca--;
	     goto loop5;
	   }
	  
	  choose1 = choose;
	  choose += point->Trapped21_to_Auton;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 19;
		 point->CaM21++;
	     goto loop5;
	   }

	  }else if(p1->subunit[i][j] == 17){
      choose1 = choose;
	  choose += point->Trapped12_to_CaMKIICaM12;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 8;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped12_to_Trapped02;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 15;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped12_to_Trapped11;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 13;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Ca * point->Trapped12_to_Trapped22;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 18;
		 point->Ca--;
	     goto loop5;
	   }

	  choose1 = choose;
	  choose += point->Trapped21_to_Auton;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 19;
		 point->CaM21++;
	     goto loop5;
	   }

	  }else if(p1->subunit[i][j] == 18){
      choose1 = choose;
	  choose += point->Trapped22_to_CaMKIICaM22;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 9;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped22_to_Trapped12;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 17;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Trapped22_to_Trapped21;
	  if((choose > random_num) && (choose1 <= random_num)){
	  point->Ca++;
	  p1->subunit[i][j] = 16;
	  goto loop5;
	  }
	  
	  choose1 = choose;
	  choose += point->Trapped22_to_Auton;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 19;
		 point->CaM22++;
	     goto loop5;
	   }

	  }else if(p1->subunit[i][j] == 19){
          choose1 = choose;
	  choose += point->CaM00 * point->Auton_to_Trapped00;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 10;
	  point->CaM00--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM10 * point->Auton_to_Trapped10;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 11;
	  point->CaM10--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM01 * point->Auton_to_Trapped01;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 12;
	  point->CaM01--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM11 * point->Auton_to_Trapped11;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 13;
	  point->CaM11--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM20 * point->Auton_to_Trapped20;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 14;
	  point->CaM20--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM02 * point->Auton_to_Trapped02;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 15;
	  point->CaM02--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM21 * point->Auton_to_Trapped21;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 16;
	  point->CaM21--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM12 * point->Auton_to_Trapped12;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 17;
	  point->CaM12--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->CaM22 * point->Auton_to_Trapped22;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 18;
	  point->CaM22--;
	  goto loop5;
	  }

	  choose1 = choose;
	  choose += point->Auton_to_CaMKII;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 0;
	  goto loop5;
	  }

          choose1 = choose;
	  if((p1->subunit[i][(j+5)%6] == 4) || (p1->subunit[i][(j+5)%6] == 5) ||
		 (p1->subunit[i][(j+5)%6] == 6) || (p1->subunit[i][(j+5)%6] == 7) ||
		 (p1->subunit[i][(j+5)%6] == 8) || (p1->subunit[i][(j+5)%6] == 9) ||
		 (p1->subunit[i][(j+5)%6] == 10) || (p1->subunit[i][(j+5)%6] == 11) ||
		 (p1->subunit[i][(j+5)%6] == 12) || (p1->subunit[i][(j+5)%6] == 13) ||
		 (p1->subunit[i][(j+5)%6] == 14) || (p1->subunit[i][(j+5)%6] == 15) ||
		 (p1->subunit[i][(j+5)%6] == 16) || (p1->subunit[i][(j+5)%6] == 17) ||
		 (p1->subunit[i][(j+5)%6] == 18) ||(p1->subunit[i][(j+7)%6] == 4) ||
		 (p1->subunit[i][(j+7)%6] == 5) || (p1->subunit[i][(j+7)%6] == 6) || 
		 (p1->subunit[i][(j+7)%6] == 7) || (p1->subunit[i][(j+7)%6] == 8) ||
		 (p1->subunit[i][(j+7)%6] == 9) || (p1->subunit[i][(j+7)%6] == 10) || 
		 (p1->subunit[i][(j+7)%6] == 11) || (p1->subunit[i][(j+7)%6] == 12) || 
		 (p1->subunit[i][(j+7)%6] == 13) || (p1->subunit[i][(j+7)%6] == 14) || 
		 (p1->subunit[i][(j+7)%6] == 15) || (p1->subunit[i][(j+7)%6] == 16) || 
		 (p1->subunit[i][(j+7)%6] == 17) || (p1->subunit[i][(j+7)%6] == 18))		 
   	    choose += point->Auton_to_Capped;
	   if((choose > random_num) && (choose1 <= random_num)){
	     p1->subunit[i][j] = 20;
	     goto loop5;
	   }

	  choose1 = choose;
/*
	  if((p1->subunit[i][(j+5)%6] == 1) || (p1->subunit[i][(j+5)%6] == 2) ||
		 (p1->subunit[i][(j+5)%6] == 3) || (p1->subunit[i][(j+5)%6] == 19) ||
		 (p1->subunit[i][(j+5)%6] == 20) || (p1->subunit[i][(j+7)%6] == 1) ||
		 (p1->subunit[i][(j+7)%6] == 2) || (p1->subunit[i][(j+7)%6] == 3) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->Auton_to_Capped;
*/
//
	  if((p1->subunit[i][(j+5)%6] == 19) || (p1->subunit[i][(j+5)%6] == 20) ||
                 (p1->subunit[i][(j+7)%6] == 19) || (p1->subunit[i][(j+7)%6] == 20))
	    choose += 0.4 * point->Auton_to_Capped;
//
	  if((choose > random_num) && (choose1 <= random_num)){
	    p1->subunit[i][j] = 20;
	    goto loop5;
	  }

	  }else{
          choose1 = choose;
	  choose += point->Capped_to_Auton;
	  if((choose > random_num) && (choose1 <= random_num)){
	  p1->subunit[i][j] = 19;
	  goto loop5;
	  }
	  }
	  }
    }
    p1 = p1->next;
  }    

 loop5: v1 = point->diffusion_lamda;
  v2 = point->chemical_lamda;
  if(index == 1){
    point->diffusion_lamda = layer1_diffusion_lamda(x, y);
    point->chemical_lamda = voxel_chemical_lamda(layer1[x][y]);
  }else if(index == 2){
    point->diffusion_lamda = layer2_diffusion_lamda(x, y, z);
    point->chemical_lamda = voxel_chemical_lamda(layer2[x][y][z]);
  }else if(index == 3){
    point->diffusion_lamda = layer3_diffusion_lamda(x);
    point->chemical_lamda = voxel_chemical_lamda(layer3[x]);
  }else{
    point->diffusion_lamda = layer4_diffusion_lamda(x, y, z);
    point->chemical_lamda = voxel_chemical_lamda(layer4[x][y][z]);
  }
  reaction_lamda += point->diffusion_lamda + point->chemical_lamda - v1 - v2;
}   


struct CaMKII * delete_CaMKII(struct voxel grid, struct CaMKII *point){
  struct CaMKII *p1, *p2, *head;

  head = grid.head;
  if(head == NULL) cout << "empty list" << endl;
  p1 = head;
  while(point->index != p1->index && p1->next != NULL){
    p2 = p1;
    p1 = p1->next;
  }
  if(point->index == p1->index){
    if(p1 == head) head = p1->next;
    else p2->next = p1->next;
  }

  return head;
}

struct CaMKII * insert_CaMKII(struct voxel grid, struct CaMKII *point){
 
  point->next = grid.head;
  return point;
}
	
double ran2(long *idum){
  long const IM1 = 2147483563;
  long const IM2 = 2147483399;
  double const AM = 1.0/IM1;
  long const IMM1 = IM1-1;
  int const IA1 = 40014;
  int const IA2 = 40692;
  int const IQ1 = 53668;
  int const IQ2 = 52774;
  int const IR1 = 12211;
  int const IR2 = 3791;
  int const NTAB = 32;
  long const NDIV = (1 + IMM1/NTAB);
  double const EPS = 1.2e-7;
  double const RNMX = 1.0 - EPS;

  int j;
  long k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  double temp;
  
  if (*idum <= 0){
    if ( -(*idum) < 1) *idum = 1;
    else *idum = -(*idum);
    idum2 = (*idum);
    for ( j = NTAB + 7; j >= 0; j--){
      k = (*idum)/ IQ1;
      *idum = IA1 * (*idum - k*IQ1) - k * IR1;
      if ( *idum < 0) *idum += IM1;
      if ( j < NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }
  
  k = (*idum) / IQ1;
  *idum = IA1 * (*idum - k*IQ1) - k * IR1;
  if ( *idum < 0) *idum += IM1;
  k = idum2 / IQ2;
  idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
  if (idum2 < 0) idum2 += IM2;
  j = iy/ NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  
  if ( iy < 1) iy += IMM1;
  if (( temp = AM*iy) > RNMX){
    return RNMX;
  }
  else{
    return temp;
  }
}
