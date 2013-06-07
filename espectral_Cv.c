#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NF 1000 // número de archivos
#define L 16 // tamaño lineal de la red
#define N L*L
#define K 0.8 // kappa correspondiente a la configuración
#define F 1000 //número de archivos termalizados
#define NP 50 //numero de extrapolaciones
#define A 3 // factor sigma de se
#define NJK 10 //tamaño bloque jacknife 
double se[F],se2[F],Se[F],Se2[F],cv;//observables de lectura
double se_media[NP],se2_media[NP];

double se_estimado(double);
double se2_estimado(double);

int main(void )
{
  FILE *input_se,*output_cvpro;
  char se_data[256],cvpro_data[256];
  
  double se_new,se2_new,cv_new;
  double k,l;
  double sigma_se;
  int i;

  double trash;

  //Lee el valor de se para cada configuración

  sprintf(se_data,"./MEDIDAS_Cv/L%d/K%.1f/Cv_L%d_K%.1f.dat",L,K,L,K);
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

  //Cálculo de los nuevos observables

  sprintf(cvpro_data,"./MEDIDAS_Cv/L%d/extrapolacion_Cv_L%d_K%.1f.dat",L,K,L,K);
  output_cvpro=fopen(cvpro_data,"w");
  for(l=-1; l<=1; l+=(1.0F/(double) NP))
    {
      k=K+l*((double)A/sigma_se);     
      se_new=se_estimado(k);
      se2_new=se2_estimado(k);
      cv_new=se2_new-se_new*se_new;
      cv_new*=(k*k)/((double) N);
      fprintf(output_cvpro,"%lf %lf \n",k,cv_new);
    }
  close(output_cvpro);
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
