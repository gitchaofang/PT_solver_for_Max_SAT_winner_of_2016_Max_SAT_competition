#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include<stdint.h>
//#include"header.h"
static inline float faster_exp (float x) {
  float p;
  p=1.442695040f *x;
  union { uint32_t i2; float f; } v = { (1 << 23) * (p + 126.94269504f) };
  return v.f;
}
int k=0,m=0,n=0,actualn=0,replica=1,key=350,norreplica=1,norkey=350,tnumb=50,nortnumb=25;


void metro (double rf,int *replicaenergy,int l,unsigned char *config,int *ksat,int *absksat,int **ksatm,int **ksatk,int **ksatmpos,int **ksatkpos,int **ksatmneg,int **ksatkneg,int *repconfig,int *repconfigpos,int *repconfigneg,double t,int rn,unsigned char *clauseenergy,int *clauseenergyone,int *clauseenergyzero,unsigned char *ksatclause);
int  energy(int *clauseenergyzero,int *clauseenergyone,unsigned char *config,int *ksat,int *absksat,unsigned char *ksatclause,unsigned char *clauseenergy);
void parallel(double r,double *t,int l,int *replicaenergy,unsigned char **config,unsigned char  **clauseenergy,int **clauseenergyone,int **clauseenergyzero);
long seedgen(void);

int main(int argc,char  *argv[]){
  clock_t tic=clock();
  long seed;
  seed=seedgen();
  srand(seed);
  /*set up the variables and arrays*/
  int i=0,j=0,l=0,d=0,f=0,u=0,v=0,digit=1,countkey=0,best=0,temptype=0;
  unsigned char y=0,base=1,crit=1;
  long energydiff=0;
  double *t,rf=0,time=0;
  FILE *fp,*fp1;
  char str[30],type[10],*cnf;
/*-------------------cnf reading 1------------------------*/
  fp=fopen(argv[1],"r");
/*------------------identify the type---------------------*/
char line[50];
unsigned char pan=0;
i=0;
while (*(argv[1]+i)!='H'){
  i++;
}
if ((*(argv[1]+i)=='H')&&(*(argv[1]+i+1)=='G')){
  pan++;;
}
line[0]='q';
while (line[0]!='p'){
  fscanf(fp,"%s",line);
  if ((line[0]=='g') && (line[1]=='i') && (line[2]=='r') && (line[3]=='t') && (line[4]=='h')){
    pan++;
    line[0]='p';
  } 
}
if (pan==0){
  key=norkey;
  replica=norreplica;
  temptype=1;
  tnumb=nortnumb;
}
/*--------------------------------------------------------*/
  rewind(fp);
  while (str[0]!='p'){
    fscanf(fp,"%s",str);
  }
  fscanf(fp,"%s",type);
  /*set up n m and k*/
  fscanf(fp,"%d",&actualn);
  fscanf(fp,"%d",&m);
  /*---------set up for k----------- */
  for (i=0;i<=m-1;i++){
    j=0;
    fscanf(fp,"%d",&d); 
    while (d!=0){
      j++;  
      fscanf(fp,"%d",&d);
    }
    k=((j>k)*(j-k)+k);
  }
  rewind(fp);
  fscanf(fp,"%s",str);
  while(str[0]!='p'){
    fscanf(fp,"%s",str);
  }
  fscanf(fp,"%s",str);
  fscanf(fp,"%s",str);
  fscanf(fp,"%s",str);   
  n=actualn+k;
  u=n;
  while(u>=1){
    digit++;
    u=u/10;
  }
  cnf=(char*)calloc(m*((digit+1)*k+k+2),sizeof(char));
  fread(cnf,m*((digit+1)*k+k+2),1,fp);
  fclose(fp);
/*---------------------------------------------------*/
  t=(double*)calloc(tnumb,sizeof(double));  

  int *ksat;
  ksat=(int*)calloc(m*k,sizeof(int));
  

  int *absksat;
  absksat=(int*)calloc(m*k,sizeof(int));

  unsigned char ***config,*solution;
  config=(unsigned char***)calloc(replica,sizeof(unsigned char**));
  for (i=0;i<=replica-1;++i){
    config[i]=(unsigned char**)calloc(tnumb,sizeof(unsigned char*));
    for (j=0;j<=tnumb-1;++j){
      *(config[i]+j)=(unsigned char*)calloc(n,sizeof(unsigned char));
    }
  }
  solution=(unsigned char*)calloc(actualn,sizeof(unsigned char));
  
  int *repconfig,*repconfigpos,*repconfigneg;
  repconfig=(int*)calloc(n,sizeof(int));
  repconfigpos=(int*)calloc(n,sizeof(int));
  repconfigneg=(int*)calloc(n,sizeof(int));

  int **ksatm,**ksatmpos,**ksatmneg;  
  ksatm=(int**)calloc(n,sizeof(int*));
  ksatmpos=(int**)calloc(n,sizeof(int*));
  ksatmneg=(int**)calloc(n,sizeof(int*));
  
  int **ksatk,**ksatkpos,**ksatkneg;
  ksatk=(int**)calloc(n,sizeof(int*));
  ksatkpos=(int**)calloc(n,sizeof(int*));
  ksatkneg=(int**)calloc(n,sizeof(int*));

  unsigned char *ksatclause;
  ksatclause=(unsigned char*)calloc(m*k,sizeof(unsigned char));
  
  unsigned char  ***clauseenergy;
  clauseenergy=(unsigned char***)calloc(replica,sizeof(unsigned char**));
  int ***clauseenergyzero;
  clauseenergyzero=(int***)calloc(replica,sizeof(int**));
  int ***clauseenergyone;
  clauseenergyone=(int***)calloc(replica,sizeof(int**));
  for (i=0;i<=replica-1;++i){
    clauseenergy[i]=(unsigned char**)calloc(tnumb,sizeof(unsigned char*));
    clauseenergyzero[i]=(int**)calloc(tnumb,sizeof(int*));
    clauseenergyone[i]=(int**)calloc(tnumb,sizeof(int*));
    for (j=0;j<=tnumb-1;++j){
      *(clauseenergy[i]+j)=(unsigned char*)calloc(m,sizeof(unsigned char));
      *(clauseenergyzero[i]+j)=(int*)calloc(n,sizeof(int));
      *(clauseenergyone[i]+j)=(int*)calloc(n,sizeof(int));
    }
  }

  int **replicaenergy;
  replicaenergy=(int**)calloc(replica,sizeof(int*));
  for (i=0;i<=replica-1;++i){
    replicaenergy[i]=(int*)calloc(tnumb,sizeof(int));
  }

  int *count,*countpos,*countneg;
  count=(int*)calloc(n,sizeof(int));
  countpos=(int*)calloc(n,sizeof(int));
  countneg=(int*)calloc(n,sizeof(int)); 
 
  int *sign;
  sign=(int*)calloc(n,sizeof(int));

  /*----------------------------partial--------------------------*/
  int *tempreplicaenergy;
  tempreplicaenergy=(int*)calloc(replica,sizeof(int));
/*---------------------------------cnf reading 2------------------------------------*/  
  if ((type[0]=='c') && (type[1]=='n') && (type[2]=='f')){
  //  printf("cnf\n");
    int *literalcount;
    literalcount=(int*)calloc(n,sizeof(int));
    v=0;
    while (cnf[v]!='\n'){
      v++;
    }
    v++;
    for (i=0;i<=m-1;++i){
      j=0;
      /*------read a string------*/
      u=0;
      str[u]=cnf[v];
      while(str[u]!=' '){
        u++;
        v++;
        str[u]=cnf[v];
      }
      v++;
      str[u]='\0';
      d=atoi(str);
      f=1;
      while (d!=0){
        if(literalcount[abs(d)-1]!=i+1){
          ksat[i*k+j]=d;
          absksat[i*k+j]=abs(ksat[i*k+j]);
          literalcount[abs(d)-1]=i+1;
          ++repconfig[abs(ksat[i*k+j])-1];
          repconfigpos[abs(ksat[i*k+j])-1]+=(ksat[i*k+j]>0);
          repconfigneg[abs(ksat[i*k+j])-1]+=(ksat[i*k+j]<0);
          sign[abs(d)-1]=d/abs(d);
          j++;
          /*------read a string------*/
          u=0;
          str[u]=cnf[v];
          while((str[u]!=' ')&&(str[u]!='\n')){
            u++;
            v++;
            str[u]=cnf[v];
          } 
          v++;
          str[u]='\0';
          d=atoi(str);
        }
        else if(sign[abs(d)-1]*d>0){
          ksat[i*k+j]=absksat[i*k+j]=(actualn+f);  
          literalcount[absksat[i*k+j]-1]=i+1;
          ++repconfig[abs(ksat[i*k+j])-1];
          repconfigpos[abs(ksat[i*k+j])-1]+=(ksat[i*k+j]>0);
          repconfigneg[abs(ksat[i*k+j])-1]+=(ksat[i*k+j]<0);
          j++;
          f++;
          /*------read a string------*/
          u=0;
          str[u]=cnf[v];
          while((str[u]!=' ')&&(str[u]!='\n')){
            u++;
            v++;
            str[u]=cnf[v];
          }
          v++;
          str[u]='\0';
          d=atoi(str);
        }
        else {
          for (l=0;l<=j-1;++l){
            --repconfig[abs(ksat[i*k+l])-1];
            repconfigpos[abs(ksat[i*k+l])-1]-=(ksat[i*k+l]>0);
            repconfigneg[abs(ksat[i*k+l])-1]-=(ksat[i*k+l]<0);
            literalcount[abs(ksat[i*k+l])-1]=0;
          }
          while(cnf[v]!='\n'){
            v++;
          }
          --m;
          --i;
          v++;
          d=0;
          j=k;
        }  
      }
      while (j<k){
        ksat[i*k+j]=absksat[i*k+j]=(actualn+f);
        ++repconfig[actualn+f-1];
        repconfigpos[actualn+f-1]+=(ksat[i*k+j]>0);
        repconfigneg[actualn+f-1]+=(ksat[i*k+j]<0);
        j++;
        f++;
      }     
    }
  }
  free(cnf);
  cnf=NULL;
/*----------------------------------------------------------------------------------------------------------*/  
  /*read optimized tempratures*/
  if (temptype==0){
     t[0]=0.1000;
     t[1]=0.1042;
     t[2]=0.1087;
     t[3]=0.1138;
     t[4]=0.1197;
     t[5]=0.1262;
     t[6]=0.1331;
     t[7]=0.1399;
     t[8]=0.1468;
     t[9]=0.1537;
     t[10]=0.1604;
     t[11]=0.1668;
     t[12]=0.1731;
     t[13]=0.1794;
     t[14]=0.1855;
     t[15]=0.1916;
     t[16]=0.1979;
     t[17]=0.2040;
     t[18]=0.2100;
     t[19]=0.2158;
     t[20]=0.2218;
     t[21]=0.2277;
     t[22]=0.2339;
     t[23]=0.2399;
     t[24]=0.2458;
     t[25]=0.2515;
     t[26]=0.2574;
     t[27]=0.2634;
     t[28]=0.2694;
     t[29]=0.2755;
     t[30]=0.2816;
     t[31]=0.2877;
     t[32]=0.2937;
     t[33]=0.2996;
     t[34]=0.3054;
     t[35]=0.3114;
     t[36]=0.3174;
     t[37]=0.3235;
     t[38]=0.3297;
     t[39]=0.3361;
     t[40]=0.3427;
     t[41]=0.3492;
     t[42]=0.3558;
     t[43]=0.3622;
     t[44]=0.3684;
     t[45]=0.3745;
     t[46]=0.3806;
     t[47]=0.3868;
     t[48]=0.3933;
     t[49]=0.4000;
  }
  else{
    t[0]=0.0500;
    t[1]=0.0639;
    t[2]=0.0815;
    t[3]=0.1035;
    t[4]=0.1305;
    t[5]=0.1624;
    t[6]=0.1960;
    t[7]=0.2273;
    t[8]=0.2560;
    t[9]=0.2834;
    t[10]=0.3112;
    t[11]=0.3405;
    t[12]=0.3729;
    t[13]=0.4090;
    t[14]=0.4485;
    t[15]=0.4910;
    t[16]=0.5368;
    t[17]=0.5871;
    t[18]=0.6437;
    t[19]=0.7104;
    t[20]=0.7848;
    t[21]=0.8680;
    t[22]=0.9709;
    t[23]=1.0909;
    t[24]=1.2345;
  }

  
  /*allocate storae for 'ksatm', 'ksatk'*/
  for (l=0;l<=n-1;++l){
    ksatm[l]=(int*)calloc(repconfig[l],sizeof(int));
    ksatk[l]=(int*)calloc(repconfig[l],sizeof(int));
    ksatmpos[l]=(int*)calloc(repconfigpos[l],sizeof(int));
    ksatkpos[l]=(int*)calloc(repconfigpos[l],sizeof(int));
    ksatmneg[l]=(int*)calloc(repconfigneg[l],sizeof(int));
    ksatkneg[l]=(int*)calloc(repconfigneg[l],sizeof(int));
  }
  
  /*set up for ksatclause*/
  for (i=0;i<=m-1;++i){
    for (j=0;j<=k-1;++j){
      l=(absksat[i*k+j]-1);
      *(ksatm[l]+count[l])=i;
      *(ksatk[l]+count[l])=j;
      ++count[l];
      if (ksat[i*k+j]>0){ 
        ksatclause[i*k+j]=1;
        *(ksatmpos[l]+countpos[l])=i;
        *(ksatkpos[l]+countpos[l])=j;
        ++countpos[l];
      }
      else if(ksat[i*k+j]<0){
        ksatclause[i*k+j]=0;
        *(ksatmneg[l]+countneg[l])=i;
        *(ksatkneg[l]+countneg[l])=j;
        ++countneg[l];
      }
    }
  }
  free(count);
  free(countpos);
  free(countneg);
  countneg=NULL; 
  countpos=NULL; 
  count=NULL;               

  /*read initial into configs*/
/*  for (u=0;u<=replica-1;++u){
    for (j=0;j<=tnumb-1;++j){
      for (i=0;i<=actualn-1;++i){
        *(*(config[u]+j)+1)=rand()%2;
      }
    }
  }*/

  /*initialize the energies for all replicas*/
  for (i=0;i<=replica-1;i++){
    for (j=0;j<=tnumb-1;j++){
      *(replicaenergy[i]+j)=energy(*(clauseenergyzero[i]+j),*(clauseenergyone[i]+j),*(config[i]+j),ksat,absksat,ksatclause,*(clauseenergy[i]+j));
    }
    tempreplicaenergy[i]=*(replicaenergy[i]+0);
  }
  best=tempreplicaenergy[0];

  /*MCS*/
   clock_t toc=clock();
   time=((double)(toc-tic)/CLOCKS_PER_SEC);
   while(time<=280.00){
     crit=1;
   /*-----------------------------------------------------the first set of replica-----------------------------------------------*/  
   /*metro*/
    for (l=0;l<=actualn-1;++l){
    //  rf=(double)rand()/RAND_MAX;
      for (j=0;j<=(tnumb-1);++j) {
        rf=(double)rand()/RAND_MAX;
        metro(rf,replicaenergy[0],j,*(config[0]+j),ksat,absksat,ksatm,ksatk,ksatmpos,ksatkpos,ksatmneg,ksatkneg,repconfig,repconfigpos,repconfigneg,t[j],l,*(clauseenergy[0]+j),*(clauseenergyone[0]+j),*(clauseenergyzero[0]+j),ksatclause); 
        if (*(replicaenergy[0]+j)<best){
          best=(*(replicaenergy[0]+j));
          printf ("o %d\n",best);
          /*-----------------record solution---------------*/
          for (u=0;u<=actualn-1;++u){
            solution[u]=*(*(config[0]+j)+u);
          }
          /*-----------------------------------------------*/
        }
      }
    }
    /*parallel tempering*/
    for (j=0;j<=(tnumb-2);++j){
      rf=(double)rand()/RAND_MAX;
      parallel(rf,t,j,replicaenergy[0],config[0],clauseenergy[0],clauseenergyone[0],clauseenergyzero[0]);
    }
    toc=clock();
    time=((double)(toc-tic)/CLOCKS_PER_SEC);
 //   printf("%d   %d\n",i,tempreplicaenergy[0]);
  }
  
  /*---------------------------result----------------------*/
   printf("v");
   for (j=0;j<=actualn-1;++j){
     printf(" %d",(solution[j]*2-1)*(j+1));
   }
 /*free arrays */
  free(clauseenergy);
  clauseenergy=NULL;

  free(clauseenergyone);
  clauseenergyone=NULL; 
  
  free(clauseenergyzero);
  clauseenergyzero=NULL;
  
  free(absksat);
  absksat=NULL;
  
  free(ksat);
  ksat=NULL;
  
  free(repconfig);
  repconfig=NULL; 
  
  
  for (i=0;i<=replica-1;++i){
    free(config[i]);
    config[i]=NULL;
  }
  free(config);
  config=NULL;
  
  for (l=0;l<=n-1;++l){
    free(ksatm[l]);
    ksatm[l]=NULL;
    free(ksatk[l]);
    ksatk[l]=NULL;
  }
  free(ksatm);
  ksatm=NULL;
  free(ksatk);
  ksatk=NULL;
  
  free(replicaenergy);
  replicaenergy=NULL;
  
  free(t);
  t=NULL;

  free(ksatclause);
  ksatclause=NULL;
  
  free(sign);
  sign=NULL;
  
  return 0;
}

