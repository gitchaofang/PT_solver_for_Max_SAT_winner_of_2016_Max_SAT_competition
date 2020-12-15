#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdint.h>
//#include"header.h"
int energy(int *clauseenergyzero,int *clauseenergyone,unsigned char *config,int *ksat,int *absksat,unsigned char *ksatclause,unsigned char *clauseenergy){
  extern int k,m,n;
  int energy=0;
  int i=0,j=0;
  for (i=0;i<=m-1;++i){
    for (j=0;j<=k-1;++j){
      *(clauseenergy+i)+=((*(config+absksat[i*k+j]-1))^ksatclause[i*k+j]); 
    }
    if (*(clauseenergy+i)==k){
      for (j=0;j<=k-1;++j){  
        (*(clauseenergyzero+(absksat[i*k+j]-1)))++;
      }
    }
    if (*(clauseenergy+i)==(k-1)){
      for (j=0;j<=k-1;++j){      
        if (((*(config+absksat[i*k+j]-1))^ksatclause[i*k+j])==0){         
          (*(clauseenergyone+absksat[i*k+j]-1))++;
        }
      }
    }
    energy+=(*(clauseenergy+i)/k);
  }
  return energy;
}
