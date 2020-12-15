#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdint.h>
static inline float faster_exp (float x) {
  float p;
  p=1.442695040f *x;
  union { uint32_t i; float f; } v = { (1 << 23) * (p + 126.94269504f) };
  return v.f;
}
void metro (double rf,int *replicaenergy,int l,unsigned char *config, int *ksat,int *absksat,int **ksatm,int **ksatk,int **ksatmpos,int **ksatkpos,int **ksatmneg,int **ksatkneg,int *repconfig,int *repconfigpos,int *repconfigneg,double t,int rn,unsigned char  *clauseenergy,int *clauseenergyone,int *clauseenergyzero,unsigned char *ksatclause){
  extern int k,m,n;
  int i=0,j=0,tempclause=0,y=0,temp=0,w=0,wpos=0,wneg=0;
  unsigned char tempconfig=0,tempclauseenergy=0,crit=0,flip=0,base=1;
  int **literalmpos,**literalmneg,**literalkpos,**literalkneg; 
  double energydiff=0;
  w=repconfig[rn];
  tempconfig=(*(config+rn));
  temp=*(clauseenergyone+rn);
  energydiff=(double)(*(clauseenergyone+rn)-*(clauseenergyzero+rn));
  if (*(config+rn)==1){
    wpos=repconfigpos[rn];
    wneg=repconfigneg[rn];
    literalmpos=ksatmpos;
    literalmneg=ksatmneg;
    literalkpos=ksatkpos;
    literalkneg=ksatkneg;
  }
  else{
    wpos=repconfigneg[rn];
    wneg=repconfigpos[rn];
    literalmpos=ksatmneg;
    literalmneg=ksatmpos;
    literalkpos=ksatkneg;
    literalkneg=ksatkpos;
  }
   
  if (energydiff<=0){
    replicaenergy[l]+=energydiff;  /*update replicaenergy*/
    /*update clause related quantities*/
    for (i=0;i<=wpos-1;++i){
      tempclause=*(literalmpos[rn]+i);    
      if (*(clauseenergy+tempclause)==(k-1)){
        for (j=0;j<=k-1;++j){
          (*(clauseenergyzero+absksat[tempclause*k+j]-1))++;
        }
        (*(clauseenergyone+rn))--;
      }
      else if (*(clauseenergy+tempclause)==(k-2)){
        j=0; 
        while ((((*(config+absksat[tempclause*k+j]-1)) ^ (ksatclause[tempclause*k+j]))!=0) || (j==*(literalkpos[rn]+i))){
          j++;  /*ksatclause[abs(ksat[tempclause*k+j])]*/
        }
        (*(clauseenergyone+absksat[tempclause*k+j]-1))++;
      }
      (*(clauseenergy+tempclause))++;    /*update clauseenergy*/    
    }    
    for (i=0;i<=wneg-1;++i){
      tempclause=*(literalmneg[rn]+i); 
      if (*(clauseenergy+tempclause)==k){
        for (j=0;j<=k-1;j++){
          (*(clauseenergyzero+absksat[tempclause*k+j]-1))--;  
        }
        (*(clauseenergyone+rn))++;
      }
      else if ((*(clauseenergy+tempclause))==(k-1)){
        j=0;
        while ((((*(config+absksat[tempclause*k+j]-1)) ^ (ksatclause[tempclause*k+j]))!=0)){
          j++;
        }
        (*(clauseenergyone+absksat[tempclause*k+j]-1))--;
      }
      (*(clauseenergy+tempclause))--;    /*update clauseenergy*/
    }
    /*update config*/
    *(config+rn)^=base;
  }

  else if (energydiff>0){      
    if (rf<((faster_exp(-energydiff/t)))){ 
      //-----------check------------------
    //  printf("haha\n`");
      //---------------------------------
      replicaenergy[l]+=energydiff;//(*(clauseenergyone+rn)-*(clauseenergyzero+rn));  /*update replicaenergy*/
      //printf("%d%s%d\n",l," ",replicaenergy[l]);
      /*update clause related quantities*/
      for (i=0;i<=wpos-1;++i){
        tempclause=*(literalmpos[rn]+i);
        if (*(clauseenergy+tempclause)==(k-1)){
          for (j=0;j<=k-1;++j){
            (*(clauseenergyzero+absksat[tempclause*k+j]-1))++;
          }
          (*(clauseenergyone+rn))--;
        }
        else if (*(clauseenergy+tempclause)==(k-2)){
          j=0;
          while ((((*(config+absksat[tempclause*k+j]-1)) ^ (ksatclause[tempclause*k+j]))!=0) || (j==*(literalkpos[rn]+i))){
            j++;  /*ksatclause[abs(ksat[tempclause*k+j])]*/
          }
          (*(clauseenergyone+absksat[tempclause*k+j]-1))++;
        }
        (*(clauseenergy+tempclause))++;    /*update clauseenergy*/
      }
      for (i=0;i<=wneg-1;++i){
        tempclause=*(literalmneg[rn]+i);
        if (*(clauseenergy+tempclause)==k){
          for (j=0;j<=k-1;j++){
            (*(clauseenergyzero+absksat[tempclause*k+j]-1))--;
          }
          (*(clauseenergyone+rn))++;
        }
        else if ((*(clauseenergy+tempclause))==(k-1)){
          j=0;
          while ((((*(config+absksat[tempclause*k+j]-1)) ^ (ksatclause[tempclause*k+j]))!=0)){
            j++;
          }
          (*(clauseenergyone+absksat[tempclause*k+j]-1))--;
        }
        (*(clauseenergy+tempclause))--;    /*update clauseenergy*/
      }
      /*update config*/
      *(config+rn)^=base;
    }
  }
}
