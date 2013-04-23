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
  FILE *input;
  FILE *output;
  FILE *filerg2;
  FILE *pipe;
  
  char namein[255];
  char nameout[255];
  char namerg2[255];

  //OBSERVABLES
  double radio2g[NF];//observable 
  double Radio2g[NF];//valor acumulado  
  double media_radio2g;// media del radio del giraton al cuadrado

  //Valores -bloque Jacknife

  double radio2gJK[NF];
  double sumJK_radio2g;

  //Errores
  double error_radio2g;

  int i,k,f,b,n;
  int ftermal,ntermal;

  int fmax;
  int nbloq;


  // LECTURA ARCHIVOS FUENTE (Correspondientes a las posiciones de los nodos de la membrana)
 
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
	  printf("Error: El fichero %s no contiene %d líneas\n",namein,N);
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
      printf("Error:El número de archivos leidos en ./RUNS/L%d/K%.1f/ no es %d\n",L,K,NF); 
      return 0;
    }

  //TERMALIZACION:

  sprintf(nameout,"./MEDIDAS_Rg2/L%d/K%.1f/termalizacionRg2_L%d_K%.1f.dat",L,K,L,K);
  output=fopen(nameout,"w");
  printf("\n Plot Termalización L=%d K=%.1f: \n",L,K);
  printf(" X=1/(nº archivos conservados)\n");
  printf(" Y=promedio del radio del giratón al cuadrado\n");
  printf(" Ind.=índice del primer archivo conservado \n");
  printf("\n Ind. X  Y(promedio)\n");
  for(f=NF/20; f<NF; f++)
    {      
      if((NF%f)==0)
      {
	      media_radio2g=(Radio2g[NF-1]-Radio2g[NF-f-1])/(double) f;
	      fprintf(output,"%lf %lf\n",1.0F/(double) f,media_radio2g);
	      printf(" %d %lf %lf\n",NF-f,1.0F/(double) f,media_radio2g);
      }
    }
  media_radio2g=Radio2g[NF-1]/(double)NF;
  fprintf(output,"%lf %lf\n",1.0F/(double)NF,media_radio2g);
  printf(" %d %lf %lf\n",0,1.0F/(double) f,media_radio2g);
  fclose(output);

  //GRÁFICA GNUPLOT TERMALIZACION:
  pipe = popen("gnuplot -persist","w");//tubería a gnuplot

  //gráfico en pantalla:
  fprintf(pipe, "set title \" Termalización L=%d K=%.1f\" \n",L,K);
  fprintf(pipe, "set logscale x\n");
  fprintf(pipe, "set xlabel\" número de sweeps/(tau=%d)\"\n",TAU);
  fprintf(pipe, "set ylabel\" promedio radio2g acumulado\"\n");
  fprintf(pipe, "plot \"%s\" title \"Rg^2\" w lp\n",nameout);
  fprintf(pipe,"set terminal push\n");

  //guardamos el gráfico en formato latex en el disco:
  fprintf(pipe,"set terminal epslatex color colortext\n");
  fprintf(pipe,"set output \"./MEDIDAS_Rg2/L%d/K%.1f/Plottermal_Rg2_L%d.tex\" \n",L,K,L,K);
  fprintf(pipe, "set title \' Termalización $L=%d$ $K=%.1f$ \' \n",L,K);
  fprintf(pipe, "set logscale x\n");
  fprintf(pipe, "set xlabel \'número de sweeps/($\\tau$=%d)\' \n",TAU);
  fprintf(pipe, "set ylabel \' promedio $R^2_g$ acumulado\'\n");
  fprintf(pipe, "plot \"%s\" title \'$R_g^2$\' w lp\n",nameout);
  
  //vuelve a la anterior terminal gnuplot:
  fprintf(pipe,"set output\n");
  fprintf(pipe,"set terminal pop\n");

  //vacía el buffer de la tubería gnuplot y se cierra:
  fflush(pipe);
  //close(pipe);
  //Entrada del valor del índice del archivo en el que se produce la termalización

  ftermal=0;

  printf("Escribe el valor de sweep/TAU correspondiente a la termalización:\n");
 
  while(scanf("%d",&ftermal)==0 || ftermal<0)
    {
      while (getchar()!= '\n');// para leer un único dato por línea
      printf("Debe ser un número entero positivo:\n"); 
    }

  fmax=NF-ftermal;

  for(f=0; f<fmax; f++)//redefinimos los observables
    {
      Radio2g[f]=Radio2g[f+ftermal]-Radio2g[ftermal-1];
    }
      
  media_radio2g=Radio2g[fmax-1]/(double)fmax;
  printf(" %d %lf %lf\n",0,1.0F/(double) f,media_radio2g);
      
  
  // Error en función del nº de bloques

  sprintf(nameout,"./MEDIDAS_Rg2/L%d/K%.1f/error_Rg2_L%d_K%.1f.dat",L,K,L,K);
  output=fopen(nameout,"w");
  printf("\n Plot Error L=%d K=%.1f: \n",L,K);
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
	  fprintf(output,"%d %lf\n",n,error_radio2g);  
	  printf(" %d %lf\n",n,error_radio2g);  
	}
    }
  fclose(output);
      
  //GRAFICA GNUPLOT ERROR
  pipe = popen("gnuplot -persist","w");//tubería a gnuplot

  //grafico en pantalla:  
  fprintf(pipe, "set title \" Error vs. tamaño del bloque Jacknife L=%d K=%.1f\" \n",L, K);
  fprintf(pipe, "set logscale x\n");
  fprintf(pipe, "set xlabel\" tamaño del bloque jacknife  (sweeps/tau)\"\n",TAU);
  fprintf(pipe, "set ylabel\"error radio2g \"\n");  
  fprintf(pipe, "plot \"%s\" title \"Radio2 g\" w lp\n",nameout);
  fprintf(pipe,"set terminal push\n");

  //guardamos el gráfico en formato latex en el disco:
  fprintf(pipe,"set terminal epslatex size 10cm,20cm color colortext\n");
  fprintf(pipe,"set output \"./MEDIDAS_Rg2/L%d/K%.1f/errorRg2_L%d_K%.1f.tex\" \n",L,K,L,K);
  fprintf(pipe, "set title \' Error $L=%d$ $K=%.1f$ \' \n",L,K);
  fprintf(pipe, "set logscale x\n");
  fprintf(pipe, "set xlabel \'tamaño del bloque Jacknife (sweeps/($\\tau_{art}$=%d))\' \n",TAU);
  fprintf(pipe, "set ylabel \' Error $R^2_g$ \'\n");
  fprintf(pipe, "plot \"%s\" title \'Error $R_g^2$\' w lp\n",nameout);

  //Vuelve a la anterior terminal gnuplot:
  fprintf(pipe,"set output\n");
  fprintf(pipe,"set terminal pop\n");

  //Vacía el buffer de la tubería gnuplot y se cierra:
  fflush(pipe);
  //close(pipe);

  
  // Radio del giratón Valor del error
  printf("Escribe el tamaño del bloque Jacknife en donde se estabiliza el error:\n");  
  while(scanf("%d",&nbloq)==0 || nbloq<0)
    {
      while (getchar()!= '\n');// para leer un único dato por línea
      printf("Debe ser un número entero positivo:\n"); 
    }
  
  input=fopen(nameout,"r");
  while(nbloq!=n)
    {
      if(fscanf(input,"%d %lf",&n,&error_radio2g)==EOF)
	{
	  close(nameout);
	  printf("Error: Escribe el número de elementos del bloque Jacknife:\n");  
	  while(scanf("%d",&nbloq)==0 || nbloq<0)
	    {
	      while (getchar()!= '\n');// para leer un único dato por línea
	      printf("Escribe un número entero positivo:\n"); 
	    }
	  input=fopen(nameout,"r");
	}
    }
  close(input);

  printf("media Rg2=%lf +- %lf\n", media_radio2g, error_radio2g);//El último valor del error corresponde a un tamaño de bloque 1

  //Escribe el valor de la medida con su error en el archivo valores_Rg2
  sprintf(namerg2,"./MEDIDAS_Rg2/L%d/valores_Rg2_L%d.dat",L,L);
  filerg2=fopen(namerg2,"a");
  fprintf(filerg2,"%lf %lf %lf\n", K, media_radio2g, error_radio2g);
  close(filerg2);
  close(pipe);
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
