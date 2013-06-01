#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double radio_giraton(void);
int radio_giraton_inicial(void);
int Rg2_files(void);

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

double radio2gflat;//Observable fase completamente plana
double radio2g[NF];//observable 
double Radio2g[NF];//valor acumulado  

int main(void )
{
  FILE *frg2_file;//Archivo de escritura/lectura = Rg2 para cada archivo 
  FILE *flogtermal;// Archivo de escritura = Termalización
  FILE *ferror;// Archivo de escritura = Error vs. Tamaño del bloque Jacknife 
  FILE *frg2_LK;// Archivo de escritura = Valores Rg2 con errores
  FILE *pipegp = popen("gnuplot -persist","w");// Tubería a gnuplot (gráficas)  
  FILE *output;
  
  char nameout[255];
  char rg2_file[255];
  char logtermal[255];
  char error[255];
  char rg2_LK[255];

  //OBSERVABLES:
  double media_radio2g;// media del radio del giraton al cuadrado
  double Radio2gref;// Valor de referencia necesario para redefinir los observables al definir el momento de termalización

  //Valores bloque Jacknife
  double radio2gJK[NF];
  double sumJK_radio2g;

  //Errores
  double error_radio2g;

  //Indices y contadores enteros
  int i,k,f,b,n;
  int indtermal;
  int fmax;
  int nbloq;



  if(radio_giraton_inicial()==1)// Rgiraton² para la conf. incial plana
    return 1;   


  // Radio del giratón al para cada configuración:

  sprintf(rg2_file,"./MEDIDAS_Rg2/L%d/K%.1f/Rg2_L%d_K%.1f.dat",L,K,L,K,f);
  if((frg2_file=fopen(rg2_file,"r"))==NULL)
    {
      if(Rg2_files()==1)// Cálcula el rg2 para cada archivo de posicines
	return 1;
    }
  else
    {
      i=0;      
      while(fscanf(frg2_file,"%lf %lf",&radio2g[i],&Radio2g[i])!=EOF)
	{
	  i++;
	}
      fclose(frg2_file);     
      if(i!=NF)
	{
	  printf("Error lineal: El fichero %s no contiene %d líneas\n",rg2_file,NF);
	  return 1;
      	}
    }
  //GRÁFICO Histórico

  sprintf(nameout,"./MEDIDAS_Rg2/L%d/K%.1f/Rg2_L%d_K%.1f.dat",L,K,L,K,f);
  fprintf(pipegp,"set terminal qt 1\n");//Guarda la terminal gp por defecto
  fprintf(pipegp,"set terminal push\n");//Guarda la terminal gp por defecto

  //Guarda el gráfico en formato latex en el disco:
  fprintf(pipegp,"set terminal epslatex color colortext\n");
  fprintf(pipegp,"set output \"./MEDIDAS_Rg2/L%d/K%.1f/Phistorico_Rg2-L%d-K%d.tex\" \n",L,K,L,(int)(10*K));
  fprintf(pipegp, "set title \' Termalización($L=%d$ $K=%.1f$) \' \n",L,K);
  fprintf(pipegp, "set xlabel \'Sweeps\' \n");
  fprintf(pipegp, "set ylabel \' $\\langle R^2_g \\rangle$ \'\n");
  fprintf(pipegp, "plot \"%s\" title \'$R_g^2$\' w lp\n",nameout);
  
  //vuelve a la anterior terminal gnuplot:
  fprintf(pipegp,"set output\n");
  fprintf(pipegp,"set terminal pop\n");

  //gráfico en pantalla:
  fprintf(pipegp, "set title \" Termalización L=%d K=%.1f\" \n",L,K);
  fprintf(pipegp, "set xlabel\" Sweeps\"\n");
  fprintf(pipegp, "set ylabel\" Radio2g \"\n");
  fprintf(pipegp, "plot \"%s\"  u 1 title \"Rg^2\"\n",nameout);
    
  fflush(pipegp);//vacía el buffer de la tubería gnuplot:
  close(pipegp);
  

  //TERMALIZACION: Escala logarítmica

  sprintf(logtermal,"./MEDIDAS_Rg2/L%d/K%.1f/logtermal_Rg2_L%d_K%.1f.dat",L,K,L,K);
  flogtermal=fopen(logtermal,"w");
  printf("\n Plot: Termalización L=%d K=%.1f \n",L,K);
  printf(" X=1/(nº archivos conservados)\n");
  printf(" Y=promedio del radio del giratón al cuadrado\n");
  printf(" Ind.=índice del primer archivo conservado \n");
  printf("\n Ind. X  Y\n");
  for(f=NF/20; f<NF; f++)
    {      
      if((NF%f)==0)
      {
	      media_radio2g=(Radio2g[NF-1]-Radio2g[NF-f-1])/(double) f;
	      fprintf(flogtermal,"%lf %lf\n",1.0F/(double) f,media_radio2g);
	      printf(" %d %lf %lf\n",NF-f,1.0F/(double) f,media_radio2g);
      }
    }
  media_radio2g=Radio2g[NF-1]/(double)NF;
  fprintf(flogtermal,"%lf %lf\n",1.0F/(double)NF,media_radio2g);
  printf(" %d %lf %lf\n",0,1.0F/(double) f,media_radio2g);
  fclose(flogtermal);

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
  fprintf(pipegp, "plot \"%s\" title \'$R_g^2$\' w lp\n",logtermal);
  
  //vuelve a la anterior terminal gnuplot:
  fprintf(pipegp,"set output\n");
  fprintf(pipegp,"set terminal pop\n");

  //gráfico en pantalla:
  fprintf(pipegp, "set title \" Termalización L=%d K=%.1f\" \n",L,K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel\" nº de archivos conservados\"\n");
  fprintf(pipegp, "set ylabel\" promedio radio2g \"\n");
  fprintf(pipegp, "plot \"%s\" title \"Rg^2\" w lp\n",logtermal);
    
  fflush(pipegp);//vacía el buffer de la tubería gnuplot:
  close(pipegp);


  
  //Entrada del valor del índice del archivo en el que se produce  printf("\n Valor Rg2 configuración plana = %lf\n",radio2gflat);
  printf("\n Escribe el valor del índice del archivo correspondiente a la termalización=");
 
  while(scanf("%d",&indtermal)==0 || indtermal<0)
    {
      while (getchar()!= '\n');// para leer un único dato por línea
      printf("\n Debe ser un número entero positivo="); 
    }

  fmax=NF-indtermal;//fmax es ahora el nº total de archivos para los cálculos
  if(fmax!=NF)
    {
      Radio2gref=Radio2g[indtermal-1];
      for(f=0; f<fmax; f++)//redefinimos los observables
	{
	  Radio2g[f]=Radio2g[f+indtermal]-Radio2gref;
	}   
      media_radio2g=Radio2g[fmax-1]/(double)fmax;
    }

  printf("\n -> Media rg2= %lf\n",media_radio2g);
    
 

  
  
  // Error en función del nº de bloques

  sprintf(error,"./MEDIDAS_Rg2/L%d/K%.1f/error_Rg2_L%d_K%.1.dat",L,K,L,K);
  ferror=fopen(error,"w");
  printf("\n Plot: Error vs. tamaño bloque Jacknife L=%d K=%.1f \n",L,K);
  printf(" X=tamaño del bloque Jacknife\n");
  printf(" Y=error \n");
  
  printf("\n X  Y\n"); 
  for(b=2; b<=fmax; b++)
    {
      if((fmax%b)==0)
	{
	  n=fmax/b;
	  error_radio2g=0.0F;

	  radio2gJK[0]=Radio2g[n-1];

	  sumJK_radio2g=(Radio2g[fmax-1]-radio2gJK[0])/(double)(fmax-n);
	  error_radio2g+=pow(sumJK_radio2g-media_radio2g,2.0);
	  
	  for (k=1; k<b; k++)
	    {
	      radio2gJK[k]=Radio2g[(k+1)*n-1]-Radio2g[k*n-1];
	      sumJK_radio2g=(Radio2g[fmax-1]-radio2gJK[k])/(double)(fmax-n);
	      error_radio2g+=pow(sumJK_radio2g-media_radio2g,2.0);
	    }

	  error_radio2g=(double)(b-1)/((double) b) * error_radio2g;
	  error_radio2g=sqrt(error_radio2g);
	  fprintf(ferror,"%d %lf\n",n,error_radio2g);  
	  printf(" %d %lf\n",n,error_radio2g);  
	}
    }
  fclose(ferror);
      
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
  fprintf(pipegp, "plot \"%s\" title \'Error $R_g^2$\' w lp\n",error);

  //Vuelve a la anterior terminal gnuplot:
  fprintf(pipegp,"set output\n");
  fprintf(pipegp,"set terminal pop\n");

  //Gráfico en pantalla:  
  fprintf(pipegp, "set title \" Error vs. tamaño del bloque Jacknife L=%d K=%.1f\" \n",L, K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel\" tamaño del bloque Jacknife (sweeps/tau)\"\n");
  fprintf(pipegp, "set ylabel\"Error radio2g \"\n");  
  fprintf(pipegp, "plot \"%s\" title \"Radio2g\" w lp\n",error);


  //Vacía el buffer de la tubería gnuplot y se cierra:
  fflush(pipegp);
  close(pipegp);

  
  // Radio del giratón vs. Valor del error:

  printf("\n Escribe el tamaño del bloque Jacknife en donde se estabiliza el error=");  
  while(scanf("%d",&nbloq)==0 || nbloq<0)
    {
      while (getchar()!= '\n');// para leer un único dato por línea
      printf("\n Debe ser un número entero positivo="); 
    }
  
  ferror=fopen(error,"r");//Abrimos el archivo de errores para lectura
  while(nbloq!=n)
    {
      if(fscanf(ferror,"%d %lf",&n,&error_radio2g)==EOF)
	{
	  close(ferror);
	  printf("\n Debe ser un número válido de elementos del bloque Jacknife=");  
	  while(scanf("%d",&nbloq)==0 || nbloq<0)
	    {
	      while (getchar()!= '\n');// para leer un único dato por línea
	      printf("\n Escribe un número entero positivo="); 
	    }
	  ferror=fopen(error,"r");
	}
    }
  close(ferror);

  printf("\n -> Rg2=%lf +- %lf\n", media_radio2g, error_radio2g);

  //Escribe el valor de la medida con su error en el archivo Medidas_Rg2
  sprintf(rg2_LK,"./MEDIDAS_Rg2/L%d/Medidas_Rg2_L%d.dat",L,L);
  frg2_LK=fopen(rg2_LK,"a");
  fprintf(frg2_LK,"%lf %lf %lf %lf %d %d\n", K, media_radio2g, error_radio2g,radio2gflat,indtermal,nbloq);
  close(frg2_LK);
  close(pipegp);
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
  radio2gflat=radio_giraton();
  return 0;
}
int Rg2_files(void )
{
  int f,i;
  FILE *input;
  FILE *output;
  char namein[255];
  char nameout[255];

  f=0;
  Radio2g[NF-1]=0.0F;
  
  sprintf(nameout,"./MEDIDAS_Rg2/L%d/K%.1f/Rg2_L%d_K%.1f.dat",L,K,L,K,f);
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
      radio2g[f]=radio_giraton();
      Radio2g[NF-1]+=radio2g[f];
      Radio2g[f]=Radio2g[NF-1];
      // Escribimos los resultados en un archivo
      fprintf(output,"%lf %lf\n",radio2g[f],Radio2g[f]);
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
