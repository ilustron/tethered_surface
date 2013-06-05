#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NF 1000 // número de archivos
#define L 16 // tamaño lineal de la red
#define N L*L
#define K 0.8 // kappa correspondiente a la configuración
#define ERROR 0.0649 // error de Cv para estimar el rango de extrapolación

double se[NF],se2[NF],Se[NF],Se2[NF],cv;

double se_estimado(double);
double se2_estimado(double);

int main(void )
{
  FILE *inputse,*output;
  char sedata[256],cvdata[256];
  
  double se_new,se2_new;
  double cv;
  double k,l;
  int i;

  //Lee el valor de se para cada confguración

  sprintf(sedata,"./MEDIDAS_Cv/L%d/K%.1f/Cv_L%d_K%.1f.dat",L,K,L,K);
  if((inputse=fopen(sedata,"r"))==NULL)
    {
      printf("Error existencial: El fichero %s no existe\n",sedata);
      return 1;
    }
  else
    {
      i=0;      
      while(fscanf(inputse,"%lf %lf %lf %lf",&se[i],&Se[i],&Se2[i],&cv)!=EOF)
	{
	  se2[i]=se[i]*se[i];
	  i++;
	}
      fclose(inputse);     
      if(i!=NF)
	{
	  printf("Error lineal: El fichero %s no contiene %d líneas\n",sedata,NF);
	  return 1;
      	}
    }
  
  printf("%lf\n",cv);
  //Cálculo de los nuevos observables

  sprintf(cvdata,"./MEDIDAS_Cv/L%d/extrapolacion_Cv_L%d_K%.1f.dat",L,K,L,K);
  output=fopen(cvdata,"w");
  for(l=-0.1; l<=0.1; l+=0.001)
    {
      k=K+l;     
      se_new=se_estimado(k);
      se2_new=se2_estimado(k);
      cv=se2_new-se_new*se_new;
      cv*=(k*k)/((double) N);
      fprintf(output,"%lf %lf \n",k,cv);
    }
  close(output);
  return 0;

  
}

double se_estimado(double k)
{
  int i;  
  double num,den;

  num=den=0.0F;
  for(i=0; i<NF; i++)
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
  for(i=0; i<NF; i++)
    {
      num+=se2[i]*exp(-(K-k)*se[i]);
      den+=exp(-(K-k)*se[i]);
    }
  return num/den;
}
