#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdint.h>
//#include"header.h"
static inline float faster_exp (float x) {
  float p;
  p=1.442695040f *x;
  union { uint32_t i; float f; } v = { (1 << 23) * (p + 126.94269504f) };
  return v.f;
}
void parallel(double r,double *t,int l,int *replicaenergy,unsigned char **config,unsigned char **clauseenergy,int **clauseenergyone,int **clauseenergyzero){
  extern int k,m,n;
  int i=0,j=0;
  double betadiff=0, energydiff=0;
  int w=0;
  unsigned char u=0;
  unsigned char *point;
  unsigned char  *clauseenergyflip;
  int *clauseflip;
  
  /*switching*/
  energydiff=(double)(replicaenergy[l]-replicaenergy[l+1]);
  betadiff=1/t[l]-1/t[l+1];
  if (faster_exp(betadiff*energydiff)<1){
    if ((faster_exp((betadiff*energydiff)))>r){
      
      //  raterecord[l]++;          
      /*update config*/
      point=config[l];
      config[l]=config[l+1];
      config[l+1]=point;
      /*update replicaenergy*/
      w=*(replicaenergy+l);
      *(replicaenergy+l)=*(replicaenergy+l+1);
      *(replicaenergy+l+1)=w;
      /*update clauseenergy*/
      clauseenergyflip=*(clauseenergy+l);
      *(clauseenergy+l)=*(clauseenergy+l+1);
      *(clauseenergy+l+1)=clauseenergyflip;
      /*update clauseenergyone and clauseenergyzero*/
      clauseflip=*(clauseenergyone+l);
      *(clauseenergyone+l)=*(clauseenergyone+l+1);
      *(clauseenergyone+l+1)=clauseflip;
      clauseflip=*(clauseenergyzero+l);
      *(clauseenergyzero+l)=*(clauseenergyzero+l+1);
      *(clauseenergyzero+l+1)=clauseflip;
    }
  }
  if (faster_exp(betadiff*energydiff)>=1){
    
    //  raterecord[l]++;      
    /*update config*/
    point=*(config+l);
    *(config+l)=*(config+l+1);
    *(config+l+1)=point;  
    /*update replicaenergy*/
    w=*(replicaenergy+l);
    *(replicaenergy+l)=*(replicaenergy+l+1);
    *(replicaenergy+l+1)=w;
    /*update clauseenergy*/
    clauseenergyflip=*(clauseenergy+l);
    *(clauseenergy+l)=*(clauseenergy+l+1);
    *(clauseenergy+l+1)=clauseenergyflip;
    /*update clauseenergyone and clauseenergyzero*/
    clauseflip=*(clauseenergyone+l);
    *(clauseenergyone+l)=*(clauseenergyone+l+1);
    *(clauseenergyone+l+1)=clauseflip;
    clauseflip=*(clauseenergyzero+l);
    *(clauseenergyzero+l)=*(clauseenergyzero+l+1);
    *(clauseenergyzero+l+1)=clauseflip;
  }
}
