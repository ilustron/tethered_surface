#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NF 10000 // número de archivos
#define L 32 // tamaño lineal de la red
#define N L*L
#define K 0.82 // kappa correspondiente a la configuración
#define F 6000 //número de archivos termalizados
#define NP 50 //numero de puntos extrapolaciones
#define A 2 // factor sigma de se
#define NJK 10 //tamaño bloque jacknife 
#define BJK F/NJK // número de bloques jacknife

double se[F],se2[F],Se[F],Se2[F];//observables de lectura
double se_media[NP],se2_media[NP],cv_pro[NP],k[NP];
double cv_critico,k_critico;

double se_estimado(double);
double se2_estimado(double);

int index_max(double*, int);

double se_estimado_JK(double ,int );
double se2_estimado_JK(double ,int );

int main(void )
{
  FILE *input_se,*output_cvpro;
  char se_data[256],cvpro_data[256];
  
  double l;
  double sigma_se;
  int i;
  int b;
  int imax;
  double trash;

  double se_JK[NP],se2_JK[NP],cv_JK[NP];
  double error_cv,error_k;

  //Lee el valor de se para cada configuración

  sprintf(se_data,"./MEDIDAS_Cv/L%d/K%.2f/Cv_L%d_K%.2f.dat",L,K,L,K);
  if((input_se=fopen(se_data,"r"))==NULL)
    {
      printf("Error existencial: El fichero %s no existe\n",se_data);
      return 1;
    }
  else
    {
      i=0;
      while((fscanf(input_se,"%lf %lf %lf %lf",&se[0],&trash,&trash,&trash)!=EOF)&&(i<(NF-F)))
	i++;
		
      se2[0]=se[0]*se[0];

      Se[F-1]+=se[0];
      Se[0]=Se[F-1];

      Se2[F-1]+=se2[0];
      Se2[0]=Se2[F-1];

      i=1;
      while(fscanf(input_se,"%lf %lf %lf %lf",&se[i],&trash,&trash,&trash)!=EOF)
	{
	  se2[i]=se[i]*se[i];

	  Se[F-1]+=se[i];
	  Se[i]=Se[F-1];

	  Se2[F-1]+=se2[i];
	  Se2[i]=Se2[F-1];;
	  i++;
	}
      fclose(input_se);
      if(i!=F)
	{
	  printf("Error lineal: El fichero %s no contiene %d líneas, tiene %d\n",se_data,F,i);
	  return 1;
      	}
    }
  
  sigma_se=sqrt(Se2[F-1]/(double)F-Se[F-1]*Se[F-1]/(double)F/(double)F);

  for(i=0,l=-1; l<=1; i++,l+=2.0F/(double) NP)
    {
      k[i]=K+l*((double)A/sigma_se);// kappas en dode vamos a extrapolar
    }
  
  
  //Cálculo de los nuevos observables

  sprintf(cvpro_data,"./MEDIDAS_Cv/L%d/extrapolacion_Cv_L%d_K%.2f.dat",L,K,L,K);
  output_cvpro=fopen(cvpro_data,"w");

  for(i=0; i<NP; i++)
    {
      se_media[i]=se_estimado(k[i]);
      se2_media[i]=se2_estimado(k[i]);
      cv_pro[i]=se2_media[i]-se_media[i]*se_media[i];
      cv_pro[i]*=(k[i]*k[i])/((double) N);
      fprintf(output_cvpro,"%lf %lf \n",k[i],cv_pro[i]);
    }
  close(output_cvpro);

  //Observables críticos

  imax=index_max(cv_pro,NP);
  cv_critico=cv_pro[imax];
  k_critico=k[imax];  
  
  //Error Jacknife

  error_cv=0.0F;
  error_k=0.0F;

  for(b=1;b<=BJK;b++)
    {
      
      for(i=0; i<NP; i++)
	{
	  se_JK[i]=se_estimado_JK(k[i],b);
	  se2_JK[i]=se2_estimado_JK(k[i],b);
	  cv_JK[i]=se2_JK[i]-se_JK[i]*se_JK[i];
	  cv_JK[i]*=(k[i]*k[i])/((double) N);
	}

      imax=index_max(cv_JK,NP);
      
      error_cv+=pow(cv_critico-cv_JK[imax],2.0);
      error_k+=pow(k_critico-k[imax],2.0);
    }
  
  error_cv=(double)(BJK-1.0F)/((double) BJK) * error_cv;
  error_cv=sqrt(error_cv);
  error_k=(double)(BJK-1.0F)/((double) BJK) * error_k;
  error_k=sqrt(error_k);

  //resultados
  printf("Kappa critico = %lf +- %lf\n",k_critico,error_k);
  printf("Cv critico = %lf+-%lf\n",cv_critico,error_cv);
  
  return 0;
}

double se_estimado(double k)
{
  int i;  
  double num,den;

  num=den=0.0F;
  for(i=0; i<F; i++)
    {
      num+=se[i]*exp(-(K-k)*se[i]);
      den+=exp(-(K-k)*se[i]);
    }
  return num/den;
}
double se2_estimado(double k)
{
  int i;  
  double num,den;

  num=den=0.0F;
  for(i=0; i<F; i++)
    {
      num+=se2[i]*exp(-(K-k)*se[i]);
      den+=exp(-(K-k)*se[i]);
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
      num+=se[i]*exp(-(K-k)*se[i]);
      den+=exp(-(K-k)*se[i]);
    }
  for(i=b*NJK; i<F; i++)
    {
      num+=se[i]*exp(-(K-k)*se[i]);
      den+=exp(-(K-k)*se[i]);
    }
  return num/den;
}
double se2_estimado_JK(double k,int b)
{
  //b es el ordinal del bloque = [0,BJK]
  int i;  
  double num,den;

  num=den=0.0F;
  for(i=0; i<NJK*(b-1); i++)
    {
      num+=se2[i]*exp(-(K-k)*se[i]);
      den+=exp(-(K-k)*se[i]);
    }
  for(i=b*NJK; i<F; i++)
    {
      num+=se2[i]*exp(-(K-k)*se[i]);
      den+=exp(-(K-k)*se[i]);
    }
  return num/den;
}
