#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double radio_giraton(void);
int radio_giraton_inicial(void);
int Rg2_file(void);

void plot_historico(void);
void plot_termalizacion(void);
void error_Jacknife(int );

#define N L*L 
#define M 2*(L-1)*(L-1)


#define _escala(u,x) \
                         u.a *= (x);\
                         u.b *= (x);\
                         u.c *= (x);

#define _norma2(u)  (u.a*u.a + u.b*u.b + u.c*u.c)


#define _suma(u,v,w) \
			u.a = v.a + w.a; \
                        u.b = v.b + w.b; \
                        u.c = v.c + w.c;
#define _resta(u,v,w) \
                         u.a = v.a - w.a; \
                         u.b = v.b - w.b; \
                         u.c = v.c - w.c; 
#define _prodesc(w,u) (w.a*u.a+w.b*u.b+w.c*u.c)
#define _inverso(u,v) \
		        u.a = -v.a; \
                        u.b = -v.b; \
                        u.c = -v.c;

typedef struct{double a,b,c;} vector;

vector x[N]; // Vectores de posición de los nodos

double rg2flat;//Observable fase completamente plana
double rg2[NF];//observable 
double Rg2[NF];//valor acumulado  
double media_rg2;// media del radio del giraton al cuadrado
double error_rg2;

int main(void )
{
  FILE *output;
  FILE *input;

  char rg2data[255];
  char namein[255];
  char nameout[255];

  //Indices y contadores enteros
  int i,k,f,n;
  int indtermal;
  int nbloq;

  
  //CONFIGURACION PLANA
  if(radio_giraton_inicial()==1)// Rgiraton² para la conf. incial plana
    return 1;   

  

  // MEDIDA/LECTURA del Radio del giratón al para cada configuración:
  sprintf(rg2data,"./MEDIDAS_Rg2/L%d/K%.1f/Rg2_L%d_K%.1f.dat",L,K,L,K);
  if((input=fopen(rg2data,"r"))==NULL)
    {
      if(Rg2_file()==1)// Cálcula el rg2 para cada archivo de posicines
	return 1;
    }
  else
    {
      i=0;      
      while(fscanf(input,"%lf %lf",&rg2[i],&Rg2[i])!=EOF)
	{
	  i++;
	}
      fclose(input);     
      if(i!=NF)
	{
	  printf("Error lineal: El fichero %s no contiene %d líneas\n",rg2data,NF);
	  return 1;
      	}
    }

  //GRÁFICO Histórico
  plot_historico();

  //TERMALIZACION: 
  plot_termalizacion();
  
  //LECTURA del índice del archivo en donde comienza la termalización
  printf("\n Valor Rg2 configuración plana = %lf\n",rg2flat);
  printf("\n Escribe el valor del índice del archivo correspondiente a la termalización=");
  while(scanf("%d",&indtermal)==0 || indtermal<0)
    {
      while (getchar()!= '\n');// para leer un único dato por línea
      printf("\n Debe ser un número entero positivo="); 
    }
      
  // ERROR JACKNIFE
  error_Jacknife(indtermal);

  // LECTURA del tamaño del bloque jacknife donde comienza el plateau:
  printf("\n Escribe el tamaño del bloque Jacknife en donde se estabiliza el error=");  
  while(scanf("%d",&nbloq)==0 || nbloq<0)
    {
      while (getchar()!= '\n');// para leer un único dato por línea
      printf("\n Debe ser un número entero positivo="); 
    } 
  sprintf(namein,"./MEDIDAS_Rg2/L%d/K%.1f/error_Rg2_L%d_K%.1f.dat",L,K,L,K); //Abrimos el archivo de errores para lectura

  input=fopen(namein,"r");
  n=0;
  while(nbloq!=n)
    {
      if(fscanf(input,"%d %lf",&n,&error_rg2)==EOF)
	{
	  close(input);
	  printf("\n Debe ser un número válido de elementos del bloque Jacknife=");  
	  while(scanf("%d",&nbloq)==0 || nbloq<0)
	    {
	      while (getchar()!= '\n');// para leer un único dato por línea
	      printf("\n Escribe un número entero positivo="); 
	    }
	  input=fopen(namein,"r");
	}
    }
  close(input);
  printf("\n -> Rg2=%lf +- %lf\n", media_rg2, error_rg2);

  //ESCRITURA RESULTADOS: 
  sprintf(nameout,"./MEDIDAS_Rg2/L%d/Medidas_Rg2_L%d.dat",L,L);
  output=fopen(nameout,"a");
  fprintf(output,"%lf %lf %lf %lf %d %d\n", K, media_rg2, error_rg2,rg2flat,NF-indtermal,nbloq);
  close(output);
  return 0;
}

double radio_giraton(void )
{
  int i;
  double sum_norma2,norma2,norma2sumx,norma2xcm;
  double densidad,momento2,r2giraton;
  vector sum_x,xcm;
  

  sum_norma2=0.0F; 

  sum_x.a=0.0F;
  sum_x.b=0.0F;  
  sum_x.c=0.0F;

  norma2=0.0F;

  for(i=0; i<N; i++)
    {
      norma2=_norma2(x[i]);
      _suma(sum_x,sum_x,x[i]);
      sum_norma2+=norma2;      
    }
  
  densidad=1.0F/ ((double) N);
  
  xcm=sum_x;
  _escala(xcm,densidad);   
  norma2xcm=_norma2(xcm);
  momento2=sum_norma2/ ((double) N);

  return r2giraton=(momento2-norma2xcm)/3.0F;
}

int radio_giraton_inicial(void )
{
  int i;
  FILE *input;
  char namein[255];

  sprintf(namein,"./RUNS/L%d/K%.1f/xpos%d_inicial.dat",L,K,L);
  if((input=fopen(namein,"r"))==NULL)
    {
      printf("Error existencial: El archivo %s no existe\n",namein); 
      return 1;
    }      
  i=0;      
  while(fscanf(input,"%lf %lf %lf",&x[i].a,&x[i].b,&x[i].c)!=EOF)
    i++;
  
  fclose(input);     
  if(i!=N)
    {
      printf("Error lineal: El fichero %s no contiene %d líneas\n",namein,N);
      return 1;
    }
  rg2flat=radio_giraton();
  return 0;
}
int Rg2_file(void )
{
  int f,i;
  FILE *input;
  FILE *output;
  char namein[255];
  char nameout[255];

  f=0;
  Rg2[NF-1]=0.0F;
  
  sprintf(nameout,"./MEDIDAS_Rg2/L%d/K%.1f/Rg2_L%d_K%.1f.dat",L,K,L,K);
  output=fopen(nameout,"w");
  sprintf(namein,"./RUNS/L%d/K%.1f/xpos_L%d_K%.1f-%d.dat",L,K,L,K,f);
  while((input=fopen(namein,"r"))!=NULL)
    {
      i=0;      
      while(fscanf(input,"%lf %lf %lf",&x[i].a,&x[i].b,&x[i].c)!=EOF)
	{
	  i++;
	}
      fclose(input);     
      if(i!=N)
	{
	  printf("Error lineal: El fichero %s no contiene %d líneas\n",namein,N);
	  return 1;
	}
      rg2[f]=radio_giraton();
      Rg2[NF-1]+=rg2[f];
      Rg2[f]=Rg2[NF-1];
      // Escribimos los resultados en el archivo
      fprintf(output,"%lf %lf\n",rg2[f],Rg2[f]);
      f++;
      sprintf(namein,"./RUNS/L%d/K%.1f/xpos_L%d_K%.1f-%d.dat",L,K,L,K,f);
    }
  fclose(output);
  if(f!=NF)
    {
      printf("Error lectura: El número de archivos leidos en ./RUNS/L%d/K%.1f/ no es %d\n",L,K,NF); 
      return 1;
    }
  return 0;
}
void plot_historico(void )
{
  FILE *pipegp = popen("gnuplot -persist","w");// Tubería a gnuplot (gráficas)  
  char nameout[255];

  sprintf(nameout,"./MEDIDAS_Rg2/L%d/K%.1f/Rg2_L%d_K%.1f.dat",L,K,L,K);
  fprintf(pipegp,"set terminal qt 1\n");//Guarda la terminal gp por defecto
  fprintf(pipegp,"set terminal push\n");//Guarda la terminal gp por defecto

  //Guarda el gráfico en formato latex en el disco:
  fprintf(pipegp,"set terminal epslatex color colortext\n");
  fprintf(pipegp,"set output \"./MEDIDAS_Rg2/L%d/K%.1f/Phistorico_Rg2-L%d-K%d.tex\" \n",L,K,L,(int)(10*K));
  fprintf(pipegp, "set title \' Histórico medidas ($L=%d$ $K=%.1f$) \' \n",L,K);
  fprintf(pipegp, "set xlabel \'Sweeps\' \n");
  fprintf(pipegp, "set ylabel \' $\\langle R^2_g \\rangle$ \'\n");
  fprintf(pipegp, "plot \"%s\" u 1 title \'$R_g^2$\' w lp\n",nameout);
  
  //vuelve a la anterior terminal gnuplot:
  fprintf(pipegp,"set output\n");
  fprintf(pipegp,"set terminal pop\n");

  //gráfico en pantalla:
  fprintf(pipegp, "set title \" Histórico medidas L=%d K=%.1f\" \n",L,K);
  fprintf(pipegp, "set xlabel\" Sweeps\"\n");
  fprintf(pipegp, "set ylabel\" Rg2 \"\n");
  fprintf(pipegp, "plot \"%s\"  u 1 title \"Rg^2\" w lp\n",nameout);
    
  fflush(pipegp);//vacía el buffer de la tubería gnuplot:
  close(pipegp); 
}
void plot_termalizacion(void)
{
  int f;
  double media;
  FILE *pipegp = popen("gnuplot -persist","w");// Tubería a gnuplot (gráficas)  
  FILE *output;
  
  char nameout[255];

  sprintf(nameout,"./MEDIDAS_Rg2/L%d/K%.1f/logtermal_Rg2_L%d_K%.1f.dat",L,K,L,K);
  output=fopen(nameout,"w");
  printf("\n Plot: Termalización L=%d K=%.1f \n",L,K);
  printf(" X=1/(nº archivos conservados)\n");
  printf(" Y=promedio del radio del giratón al cuadrado\n");
  printf(" Ind.=índice del primer archivo conservado \n");
  printf("\n Ind. X  Y\n");
  for(f=NF/20; f<NF; f++)
    {      
      if((NF%f)==0)
      {
	      media=(Rg2[NF-1]-Rg2[NF-f-1])/(double) f;
	      fprintf(output,"%lf %lf\n",1.0F/(double) f,media);
	      printf(" %d %lf %lf\n",NF-f,1.0F/(double) f,media);
      }
    }
  media=Rg2[NF-1]/(double)NF;
  fprintf(output,"%lf %lf\n",1.0F/(double)NF,media);
  printf(" %d %lf %lf\n",0,1.0F/(double) f,media);
  fclose(output);

  //GRÁFICA GNUPLOT TERMALIZACION:

  fprintf(pipegp,"set terminal qt 2\n");
  fprintf(pipegp,"set terminal push\n");//Guarda la terminal gp por defecto

  //Guarda el gráfico en formato latex en el disco:
  fprintf(pipegp,"set terminal epslatex color colortext\n");
  fprintf(pipegp,"set output \"./MEDIDAS_Rg2/L%d/K%.1f/Plogtermal_Rg2-L%d-K%d.tex\" \n",L,K,L,(int)(10*K));
  fprintf(pipegp, "set title \' Gráfico Termalización ($L=%d$ $K=%.1f$) \' \n",L,K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel \'1/ (n\\textdegree de archivos conservados)\' \n");
  fprintf(pipegp, "set ylabel \' $\\langle R^2_g \\rangle$ \'\n");
  fprintf(pipegp, "plot \"%s\" title \'$R_g^2$\' w lp\n",nameout);
  
  //vuelve a la anterior terminal gnuplot:
  fprintf(pipegp,"set output\n");
  fprintf(pipegp,"set terminal pop\n");

  //gráfico en pantalla:
  fprintf(pipegp, "set title \" Termalización L=%d K=%.1f\" \n",L,K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel\" nº de archivos conservados\"\n");
  fprintf(pipegp, "set ylabel\" promedio rg2 \"\n");
  fprintf(pipegp, "plot \"%s\" title \"Rg^2\" w lp\n",nameout);
    
  fflush(pipegp);//vacía el buffer de la tubería gnuplot:
  close(pipegp);
}
void error_Jacknife(int index)
{
  int b,n,f,k;
  int fmax;

  double rg2JK[NF];
  double sumJK_rg2;

  double Rg2new[NF];
  double Rg2ref;

  FILE *output;// Archivo de escritura = Error vs. Tamaño del bloque Jacknife 
  char nameout[255];

  FILE *pipegp = popen("gnuplot -persist","w");// Tubería a gnuplot (gráficas)

  fmax=NF-index;//fmax es ahora el nº total de archivos para los cálculos
 
  Rg2ref=0.0F;
  if(fmax<NF)
  {
    Rg2ref=Rg2[index-1];
  }

  for(f=0; f<fmax; f++)
    {
      Rg2new[f]=Rg2[f+index]-Rg2ref;
    }   
  
  media_rg2=Rg2new[fmax-1]/(double)fmax;
  
  sprintf(nameout,"./MEDIDAS_Rg2/L%d/K%.1f/error_Rg2_L%d_K%.1f.dat",L,K,L,K);
  output=fopen(nameout,"w");
  printf("\n Plot: Error vs. tamaño bloque Jacknife L=%d K=%.1f \n",L,K);
  printf(" X=tamaño del bloque Jacknife\n");
  printf(" Y=error \n");
  
  printf("\n X  Y\n"); 
  for(b=2; b<=fmax; b++)
    {
      if((fmax%b)==0)
	{
	  n=fmax/b;
	  error_rg2=0.0F;

	  rg2JK[0]=Rg2new[n-1];

	  sumJK_rg2=(Rg2new[fmax-1]-rg2JK[0])/(double)(fmax-n);
	  error_rg2+=pow(sumJK_rg2-media_rg2,2.0);
	  
	  for (k=1; k<b; k++)
	    {
	      rg2JK[k]=Rg2new[(k+1)*n-1]-Rg2new[k*n-1];
	      sumJK_rg2=(Rg2new[fmax-1]-rg2JK[k])/(double)(fmax-n);
	      error_rg2+=pow(sumJK_rg2-media_rg2,2.0);
	    }

	  error_rg2=(double)(b-1)/((double) b) * error_rg2;
	  error_rg2=sqrt(error_rg2);
	  fprintf(output,"%d %lf\n",n,error_rg2);  
	  printf(" %d %lf\n",n,error_rg2);  
	}
    }
  fclose(output);
      
  //GRAFICA GNUPLOT ERROR
  fprintf(pipegp,"set terminal qt 3\n");
  fprintf(pipegp,"set terminal push\n");//Almacena el anterior terminal gp

  //Guardamos el gráfico en formato latex en el disco:
  fprintf(pipegp,"set terminal epslatex color colortext\n");
  fprintf(pipegp,"set output \"./MEDIDAS_Rg2/L%d/K%.1f/Perror_Rg2-L%d-K%d.tex\" \n",L,K,L,(int)(10*K));
  fprintf(pipegp, "set title \' Error vs Tamaño Jacknife ($L=%d$ $K=%.1f$) \' \n",L,K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel \'Tamaño del bloque Jacknife (sweeps/$\\tau_0$)\' \n");
  fprintf(pipegp, "set ylabel \' Error $R^2_g$ \'\n");
  fprintf(pipegp, "plot \"%s\" title \'Error $R_g^2$\' w lp\n",nameout);

  //Vuelve a la anterior terminal gnuplot:
  fprintf(pipegp,"set output\n");
  fprintf(pipegp,"set terminal pop\n");

  //Gráfico en pantalla:  
  fprintf(pipegp, "set title \" Error vs. tamaño del bloque Jacknife L=%d K=%.1f\" \n",L, K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel\" tamaño del bloque Jacknife (sweeps/tau)\"\n");
  fprintf(pipegp, "set ylabel\"Error rg2 \"\n");  
  fprintf(pipegp, "plot \"%s\" title \"Rg2\" w lp\n",nameout);


  //Vacía el buffer de la tubería gnuplot y se cierra:
  fflush(pipegp);
  close(pipegp);
}
 
