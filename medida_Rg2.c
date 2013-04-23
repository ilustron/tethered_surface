//COMPILAR CON:
//gcc -O2 -DK=2.0 -DL=16 -DNF=1000 -DTAU=16000 -DTERMAL=8000000 medida_Rg2.c -lm -o medida_Rg2.out 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double radio_giraton(void);

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

typedef struct{double a,b,c;} vector;

vector x[N];

int main(void )
{
  FILE *input;// -> Archivo de lectura = Posiciones 3d de los nodos  
  FILE *ftermal;// -> Archivo de escritura = Termalización
  FILE *ferror;// -> Archivo de escritura = Error vs. Tamaño del bloque Jacknife 
  FILE *filerg2;// -> Archivo de escritura = Valores Rg2 con errores
  FILE *pipegp = popen("gnuplot -persist","w");//Tubería a gnuplot (gráficas)  
  
  char namein[255];
  char nametermal[255];
  char nameerror[255];
  char namerg2[255];

  //OBSERVABLES
  double radio2gflat;//Observable fase completamente plana
  double radio2g[NF];//observable 
  double Radio2g[NF];//valor acumulado  
  double media_radio2g;// media del radio del giraton al cuadrado

  //Valores -bloque Jacknife

  double radio2gJK[NF];
  double sumJK_radio2g;

  //Errores
  double error_radio2g;

  int i,k,f,b,n;
  int indtermal;

  int fmax;
  int nbloq;


  // LECTURA ARCHIVOS FUENTE (Correspondientes a las posiciones de los nodos de la membrana)

  //Archivo de configuración inicial = fase completamente plana
  
  sprintf(namein,"./RUNS/L%d/K%.1f/xpos%d_inicial.dat",L,K,L);
  if((input=fopen(namein,"r"))==NULL)
    {
      printf("Error existencial: El archivo $s no existe\n",namein); 
      return 0;
    }      
  i=0;      
  while(fscanf(input,"%lf %lf %lf",&x[i].a,&x[i].b,&x[i].c)!=EOF)
    i++;
  
  fclose(input);     
  if(i!=N)
    {
      printf("Error lineal: El fichero %s no contiene %d líneas\n",namein,N);
      return 0;
    }
  radio2gflat=radio_giraton();
    
  //Lectura de lo archivos de posiciones resultantes de la simulación
  f=0;
  Radio2g[NF-1]=0.0F;
  
  sprintf(namein,"./RUNS/L%d/K%.1f/xpos_L%d_K%.1f-%d.dat",L,K,L,K,f);
  while((input=fopen(namein,"r"))!=NULL)
    {
      i=0;      
      while(fscanf(input,"%lf %lf %lf",&x[i].a,&x[i].b,&x[i].c)!=EOF)
	i++;

      fclose(input);     
      if(i!=N)
	{
	  printf("Error lineal: El fichero %s no contiene %d líneas\n",namein,N);
	  return 0;
	}
      radio2g[f]=radio_giraton();
      Radio2g[NF-1]+=radio2g[f];
      Radio2g[f]=Radio2g[NF-1];
      f++;
      sprintf(namein,"./RUNS/L%d/K%.1f/xpos_L%d_K%.1f-%d.dat",L,K,L,K,f);
    }
  
  if(f!=NF)
    {
      printf("Error lectura: El número de archivos leidos en ./RUNS/L%d/K%.1f/ no es %d\n",L,K,NF); 
      return 0;
    }

  //TERMALIZACION:

  sprintf(nametermal,"./MEDIDAS_Rg2/L%d/K%.1f/termalizacion_Rg2_L%d_K%.1f.dat",L,K,L,K);
  ftermal=fopen(nametermal,"w");
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
	      fprintf(ftermal,"%lf %lf\n",1.0F/(double) f,media_radio2g);
	      printf(" %d %lf %lf\n",NF-f,1.0F/(double) f,media_radio2g);
      }
    }
  media_radio2g=Radio2g[NF-1]/(double)NF;
  fprintf(ftermal,"%lf %lf\n",1.0F/(double)NF,media_radio2g);
  printf(" %d %lf %lf\n",0,1.0F/(double) f,media_radio2g);
  fclose(ftermal);

  //GRÁFICA GNUPLOT TERMALIZACION:
 
  fprintf(pipegp,"set terminal push\n");//Guarda la terminal gp por defecto

  //Guarda el gráfico en formato latex en el disco:
  fprintf(pipegp,"set terminal epslatex color colortext\n");
  fprintf(pipegp,"set output \"./MEDIDAS_Rg2/L%d/K%.1f/Plottermal_Rg2-L%d-K%d.tex\" \n",L,K,L,(int)(10*K));
  fprintf(pipegp, "set title \' Gráfico Termalización ($L=%d$ $K=%.1f$) \' \n",L,K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel \'1/ (n\\textdegree de archivos conservados)\' \n");
  fprintf(pipegp, "set ylabel \' $\\langle R^2_g \\rangle$ \'\n");
  fprintf(pipegp, "plot \"%s\" title \'$R_g^2$\' w lp\n",nametermal);
  
  //vuelve a la anterior terminal gnuplot:
  fprintf(pipegp,"set output\n");
  fprintf(pipegp,"set terminal pop\n");

  //gráfico en pantalla:
  fprintf(pipegp, "set title \" Termalización L=%d K=%.1f\" \n",L,K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel\" nº de archivos conservados\"\n",TAU);
  fprintf(pipegp, "set ylabel\" promedio radio2g acumulado\"\n");
  fprintf(pipegp, "plot \"%s\" title \"Rg^2\" w lp\n",nametermal);
    
  fflush(pipegp);//vacía el buffer de la tubería gnuplot:
  //close(pipegp);
  
  //Entrada del valor del índice del archivo en el que se produce la termalización

  printf("\n Valor Rg2 configuración plana = %lf\n",radio2gflat);
  
  printf("\n Escribe el valor del índice del archivo correspondiente a la termalización=");
 
  while(scanf("%d",&indtermal)==0 || indtermal<0)
    {
      while (getchar()!= '\n');// para leer un único dato por línea
      printf("\n Debe ser un número entero positivo="); 
    }

  fmax=NF-indtermal;//fmax es ahora el nº total de archivos para los cálculos
  for(f=0; f<fmax; f++)//redefinimos los observables
    {
      Radio2g[f]=Radio2g[f+indtermal]-Radio2g[indtermal-1];
    }   
  media_radio2g=Radio2g[fmax-1]/(double)fmax;
  printf("\n -> Media rg2= %lf\n",media_radio2g);
      
  // Error en función del nº de bloques

  sprintf(nameerror,"./MEDIDAS_Rg2/L%d/K%.1f/error_Rg2_L%d_K%d.dat",L,K,L,K);
  ferror=fopen(nameerror,"w");
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

  fprintf(pipegp,"set terminal push\n");//Almacena el anterior terminal gp

  //Guardamos el gráfico en formato latex en el disco:
  fprintf(pipegp,"set terminal epslatex color colortext\n");
  fprintf(pipegp,"set output \"./MEDIDAS_Rg2/L%d/K%.1f/Ploterror_Rg2-L%d-K%d.tex\" \n",L,K,L,(int)(10*K));
  fprintf(pipegp, "set title \' Error vs Tamaño Jacknife ($L=%d$ $K=%.1f$) \' \n",L,K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel \'Tamaño del bloque Jacknife (sweeps/$\\tau_0$)\' \n");
  fprintf(pipegp, "set ylabel \' Error $R^2_g$ \'\n");
  fprintf(pipegp, "plot \"%s\" title \'Error $R_g^2$\' w lp\n",nameerror);

  //Vuelve a la anterior terminal gnuplot:
  fprintf(pipegp,"set output\n");
  fprintf(pipegp,"set terminal pop\n");

  //Gráfico en pantalla:  
  fprintf(pipegp, "set title \" Error vs. tamaño del bloque Jacknife L=%d K=%.1f\" \n",L, K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel\" tamaño del bloque Jacknife (sweeps/tau)\"\n");
  fprintf(pipegp, "set ylabel\"Error radio2g \"\n");  
  fprintf(pipegp, "plot \"%s\" title \"Radio2 g\" w lp\n",nameerror);


  //Vacía el buffer de la tubería gnuplot y se cierra:
  fflush(pipegp);
  //close(pipegp);

  
  // Radio del giratón vs. Valor del error:

  printf("\n Escribe el tamaño del bloque Jacknife en donde se estabiliza el error=");  
  while(scanf("%d",&nbloq)==0 || nbloq<0)
    {
      while (getchar()!= '\n');// para leer un único dato por línea
      printf("\n Debe ser un número entero positivo="); 
    }
  
  ferror=fopen(nameerror,"r");//Abrimos el archivo de errores para lectura
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
	  ferror=fopen(nameerror,"r");
	}
    }
  close(ferror);

  printf("\n -> Rg2=%lf +- %lf\n", media_radio2g, error_radio2g);

  //Escribe el valor de la medida con su error en el archivo valores_Rg2
  sprintf(namerg2,"./MEDIDAS_Rg2/L%d/Medidas_Rg2_L%d.dat",L,L);
  filerg2=fopen(namerg2,"a");
  fprintf(filerg2,"%lf %lf %lf %d %d\n", K, media_radio2g, error_radio2g,indtermal,nbloq);
  close(filerg2);
  close(pipegp);
  return 1;
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
