#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NF 10000 // número de archivos
#define L 24 // tamaño lineal de la red
#define N L*L
#define K 0.84 // kappa correspondiente a la configuración
#define F 10000 //número de archivos termalizados
#define NP 50 //numero de puntos extrapolaciones
#define A 1 // factor sigma de se
#define NJK 20 //tamaño bloque jacknife 
#define BJK F/NJK // número de bloques jacknife

double sc[F],se[F],sc2[F],Sc[F],Sc2[F],rg2[F],rg2se[F];//observables de lectura
double se_media[NP],rg2_media[NP],rg2se_media[NP],drg2_pro[NP],k[NP];
double drg2_critico,k_critico;

double se_estimado(double);
double rg2_estimado(double);
double rg2se_estimado(double);

int index_max(double*, int);

double se_estimado_JK(double ,int );
double rg2_estimado_JK(double ,int );
double rg2se_estimado_JK(double ,int );

int main(void )
{
  FILE *input_drg2,*output_drg2pro;
  char drg2_data[256],drg2pro_data[256];
  
  double l;
  double sigma_se;
  int i;
  int b;
  int imax;
  double trash;

  double se_JK[NP],rg2_JK[NP],rg2se_JK[NP],drg2_JK[NP];
  double error_drg2,error_k;

  //Lee el valor de se para cada configuración

  sprintf(drg2_data,"./MEDIDAS_Drg2/L%d/K%.2f/Drg2_L%d_K%.2f.dat",L,K,L,K);
  if((input_drg2=fopen(drg2_data,"r"))==NULL)
    {
      printf("Error existencial: El fichero %s no existe\n",drg2_data);
      return 1;
    }
  else
    {
      i=0;
      while((fscanf(input_drg2,"%lf %lf %lf %lf %lf %lf %lf",&rg2[0],&trash,&se[0],&trash,&rg2se[0],&trash,&trash)!=EOF)&&(i<(NF-F)))
	i++;
	
      
      sc[0]=se[0]*((double) N);
      sc2[0]=sc[0]*sc[0];

      Sc[F-1]+=sc[0];
      Sc[0]=Sc[F-1];

      Sc2[F-1]+=sc2[0];
      Sc2[0]=Sc2[F-1];

      i=1;
      while(fscanf(input_drg2,"%lf %lf %lf %lf %lf %lf %lf",&rg2[i],&trash,&se[i],&trash,&rg2se[i],&trash,&trash)!=EOF)
	{
	  
	  sc[i]=se[i]*((double) N);
	  sc2[i]=sc[i]*sc[i];

	  Sc[F-1]+=sc[i];
	  Sc[i]=Sc[F-1];

	  Sc2[F-1]+=sc2[i];
	  Sc2[i]=Sc2[F-1];;
	  i++;
	}
      fclose(input_drg2);
      if(i!=F)
	{
	  printf("Error lineal: El fichero %s no contiene %d líneas, tiene %d\n",drg2_data,F,i);
	  return 1;
      	}
    }
  
  
  sigma_se=sqrt(Sc2[F-1]/(double)F-Sc[F-1]*Sc[F-1]/(double)F/(double)F);
  

  for(i=0,l=-1; l<=1; i++,l+=2.0F/(double) NP)
    {
      k[i]=K+l*((double)A/sigma_se);// kappas en dode vamos a extrapolar
      printf("%d %lf \n",i,k[i]);
    }
  
  
  //Cálculo de los nuevos observables

  sprintf(drg2pro_data,"./MEDIDAS_Drg2/L%d/extrapolacion_Drg2_L%d_K%.2f.dat",L,K,L,K);
  output_drg2pro=fopen(drg2pro_data,"w");

  for(i=0; i<NP; i++)
    {
      se_media[i]=se_estimado(k[i]);
      rg2_media[i]=rg2_estimado(k[i]);
      rg2se_media[i]=rg2se_estimado(k[i]);
      drg2_pro[i]=rg2se_media[i]-se_media[i]*rg2_media[i];
      fprintf(output_drg2pro,"%lf %lf %lf %lf %lf  \n",k[i],drg2_pro[i],se_media[i],rg2_media[i],rg2se_media[i]);
    }
  close(output_drg2pro);

  //Observables críticos

  imax=index_max(drg2_pro,NP);
  drg2_critico=drg2_pro[imax];
  k_critico=k[imax];  
  
  //Error Jacknife

  error_drg2=0.0F;
  error_k=0.0F;

  for(b=1;b<=BJK;b++)
    {
      
      for(i=0; i<NP; i++)
	{
	  se_JK[i]=se_estimado_JK(k[i],b);
	  rg2_JK[i]=rg2_estimado_JK(k[i],b);
	  rg2se_JK[i]=rg2se_estimado_JK(k[i],b);
	  drg2_JK[i]=rg2se_JK[i]-se_JK[i]*rg2_JK[i];
	}

      imax=index_max(drg2_JK,NP);
      
      error_drg2+=pow(drg2_critico-drg2_JK[imax],2.0);
      error_k+=pow(k_critico-k[imax],2.0);
    }
  
  error_drg2=(double)(BJK-1.0F)/((double) BJK) * error_drg2;
  error_drg2=sqrt(error_drg2);
  error_k=(double)(BJK-1.0F)/((double) BJK) * error_k;
  error_k=sqrt(error_k);

  //resultados
  printf("Kappa critico = %lf +- %lf\n",k_critico,error_k);
  printf("Drg2 critico = %lf+-%lf\n",drg2_critico,error_drg2);
  
  return 0;
}

double se_estimado(double k)
{
  int i;  
  double num,den;

  num=den=0.0F;
  for(i=0; i<F; i++)
    {
      num+=se[i]*exp(-(K-k)*sc[i]);
      den+=exp(-(K-k)*sc[i]);
    }
  return num/den;
}
double rg2_estimado(double k)
{
  int i;  
  double num,den;

  num=den=0.0F;
  for(i=0; i<F; i++)
    {
      num+=rg2[i]*exp(-(K-k)*sc[i]);
      den+=exp(-(K-k)*sc[i]);
    }
  return num/den;
}
double rg2se_estimado(double k)
{
  int i;  
  double num,den;

  num=den=0.0F;
  for(i=0; i<F; i++)
    {
      num+=rg2se[i]*exp(-(K-k)*sc[i]);
      den+=exp(-(K-k)*sc[i]);
    }
  return num/den;
}
int index_max(double* valores, int num) 
{ 
  int i; 
  double max; 
  int imax;

  max = valores[0]; 
  imax=0;
  for (i = 1; i < num; i++) 
    if (valores[i] > max) 
      {
	  max = valores[i];
	  imax=i;
      }
  return imax; 
} 
double se_estimado_JK(double k,int b)
{
  //b es el ordinal del bloque = [0,BJK]
  int i;  
  double num,den;

  num=den=0.0F;
  for(i=0; i<NJK*(b-1); i++)
    {
      num+=se[i]*exp(-(K-k)*sc[i]);
      den+=exp(-(K-k)*sc[i]);
    }
  for(i=b*NJK; i<F; i++)
    {
      num+=se[i]*exp(-(K-k)*sc[i]);
      den+=exp(-(K-k)*sc[i]);
    }
  return num/den;
}
double rg2_estimado_JK(double k,int b)
{
  //b es el ordinal del bloque = [0,BJK]
  int i;  
  double num,den;

  num=den=0.0F;
  for(i=0; i<NJK*(b-1); i++)
    {
      num+=rg2[i]*exp(-(K-k)*sc[i]);
      den+=exp(-(K-k)*sc[i]);
    }
  for(i=b*NJK; i<F; i++)
    {
      num+=rg2[i]*exp(-(K-k)*sc[i]);
      den+=exp(-(K-k)*sc[i]);
    }
  return num/den;
}
double rg2se_estimado_JK(double k,int b)
{
  //b es el ordinal del bloque = [0,BJK]
  int i;  
  double num,den;

  num=den=0.0F;
  for(i=0; i<NJK*(b-1); i++)
    {
      num+=rg2se[i]*exp(-(K-k)*sc[i]);
      den+=exp(-(K-k)*sc[i]);
    }
  for(i=b*NJK; i<F; i++)
    {
      num+=rg2se[i]*exp(-(K-k)*sc[i]);
      den+=exp(-(K-k)*sc[i]);
    }
  return num/den;
}
