#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double energia_gauss(void );

void indice_vecnos_prox();

int sigma_min(int );
int sigma_max(int );

//#define K=2.0 // O compilar con la opcion -DK=2.0 
//#define L=16 // O compilar con la opción -DL=16 
#define N L*L 
#define M 2*(L-1)*(L-1)
//#define NF=1000 // O compilar con la opción -DNF=1000 

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
#define SI 1
#define NO 0

typedef struct{int s1,s2;}vector2D; 

vector2D basev1[6],basev2[6];

typedef struct{vector2D min,max;}rect;

rect rombov1[6],rombov2[6];

int v1[N][6];
int dirv1_ini[N],dirv1_fin[N],zv1[N];

typedef struct{double a,b,c;} vector;

vector x[N]; // Vectores de posición de los nodos
int v1[N][6]; // Indices de los nodos vecinos

int main(void )
{
  FILE *input;// -> Archivo de lectura = Posiciones 3d de los nodos  
  FILE *ftermal;// -> Archivo de escritura = Termalización
  FILE *ferror;// -> Archivo de escritura = Error vs. Tamaño del bloque Jacknife 
  FILE *filegauss;// -> Archivo de escritura = Valores Rg2 con errores
  FILE *pipegp = popen("gnuplot -persist","w");//Tubería a gnuplot (gráficas)  
  
  char namein[255];
  char nametermal[255];
  char nameerror[255];
  char namegauss[255];

  //OBSERVABLES

  double egauss[NF];//observable 
  double Egauss[NF];//valor acumulado  
  double media_egauss;// media del radio del giraton al cuadrado
  double Egaussref;//valor de referencia necesario para redefinir los observables

  //Valores bloque Jacknife

  double egaussJK[NF];
  double sumJK_egauss;

  //Errores
  double error_egauss;

  int i,k,f,b,n;
  int indtermal;

  int fmax;
  int nbloq;

  indice_vecnos_prox(); // Cargamos los indices de v1[N][6]

  // LECTURA ARCHIVOS FUENTE (Correspondientes a las posiciones de los nodos de la membrana)

  //Lectura de lo archivos de posiciones resultantes de la simulación
  f=0;
  Egauss[NF-1]=0.0F;
  
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
	  return 0;
	}
      egauss[f]=energia_gauss();
      Egauss[NF-1]+=egauss[f];
      Egauss[f]=Egauss[NF-1];
      f++;
      sprintf(namein,"./RUNS/L%d/K%.1f/xpos_L%d_K%.1f-%d.dat",L,K,L,K,f);
    }
  
  if(f!=NF)
    {
      printf("Error lectura: El número de archivos leidos en ./RUNS/L%d/K%.1f/ no es %d\n",L,K,NF); 
      return 0;
    }

  //TERMALIZACION:

  sprintf(nametermal,"./MEDIDAS_Egauss/L%d/K%.1f/termalizacion_Egauss_L%d_K%.1f.dat",L,K,L,K);
  ftermal=fopen(nametermal,"w");
  printf("\n Plot: Termalización L=%d K=%.1f \n",L,K);
  printf(" X=1/(nº archivos conservados)\n");
  printf(" Y=promedio de la energía elástica\n");
  printf(" Ind.=índice del primer archivo conservado \n");
  printf("\n Ind. X  Y\n");
  for(f=1; f<NF; f++)
    {      
      if((f%10)==0)
      {
	      media_egauss=(Egauss[NF-1]-Egauss[NF-f-1])/(double) f;
	      fprintf(ftermal,"%lf %lf\n",1.0F/(double) f,media_egauss);
	      printf(" %d %lf %lf\n",NF-f,1.0F/(double) f,media_egauss);
      }
    }
  media_egauss=Egauss[NF-1]/(double)NF;
  fprintf(ftermal,"%lf %lf\n",1.0F/(double)NF,media_egauss);
  printf(" %d %lf %lf\n",0,1.0F/(double) f,media_egauss);
  fclose(ftermal);

  //GRÁFICA GNUPLOT TERMALIZACION:
 
  fprintf(pipegp,"set terminal push\n");//Guarda la terminal gp por defecto

  //Guarda el gráfico en formato latex en el disco:
  fprintf(pipegp,"set terminal epslatex color colortext\n");
  fprintf(pipegp,"set output \"./MEDIDAS_Egauss/L%d/K%.1f/Plottermal_Egauss-L%d-K%d.tex\" \n",L,K,L,(int)(10*K));
  fprintf(pipegp, "set title \' Gráfico Termalización ($L=%d$ $K=%.1f$) \' \n",L,K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel \'1/ (n\\textdegree de archivos conservados)\' \n");
  fprintf(pipegp, "set ylabel \' $\\langle R^2_g \\rangle$ \'\n");
  fprintf(pipegp, "plot [][0:] \"%s\" title \'$R_g^2$\' w lp\n",nametermal);
  
  //vuelve a la anterior terminal gnuplot:
  fprintf(pipegp,"set output\n");
  fprintf(pipegp,"set terminal pop\n");

  //gráfico en pantalla:
  fprintf(pipegp, "set title \" Termalización L=%d K=%.1f\" \n",L,K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel\" nº de archivos conservados\"\n");
  fprintf(pipegp, "set ylabel\" promedio radio2g \"\n");
  fprintf(pipegp, "plot [][0:] \"%s\" title \"Rg^2\" w lp\n",nametermal);
    
  fflush(pipegp);//vacía el buffer de la tubería gnuplot:
  //close(pipegp);
  
  //Entrada del valor del índice del archivo en el que se produce la termalización

  printf("\n Escribe el valor del índice del archivo correspondiente a la termalización=");
 
  while(scanf("%d",&indtermal)==0 || indtermal<0)
    {
      while (getchar()!= '\n');// para leer un único dato por línea
      printf("\n Debe ser un número entero positivo="); 
    }

  fmax=NF-indtermal;//fmax es ahora el nº total de archivos para los cálculos
  if(fmax!=NF)
    {
      Egaussref=Egauss[indtermal-1];
      for(f=0; f<fmax; f++)//redefinimos los observables
	{
	  Egauss[f]=Egauss[f+indtermal]-Egaussref;
	}   
      media_egauss=Egauss[fmax-1]/(double)fmax;
    }

  printf("\n -> Media rg2= %lf\n",media_egauss);
      
  // Error en función del nº de bloques

  sprintf(nameerror,"./MEDIDAS_Egauss/L%d/K%.1f/error_Egauss_L%d_K%d.dat",L,K,L,K);
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
	  error_egauss=0.0F;

	  egaussJK[0]=Egauss[n-1];

	  sumJK_egauss=(Egauss[fmax-1]-egaussJK[0])/(double)(fmax-n);
	  error_egauss+=pow(sumJK_egauss-media_egauss,2.0);
	  
	  for (k=1; k<b; k++)
	    {
	      egaussJK[k]=Egauss[(k+1)*n-1]-Egauss[k*n-1];
	      sumJK_egauss=(Egauss[fmax-1]-egaussJK[k])/(double)(fmax-n);
	      error_egauss+=pow(sumJK_egauss-media_egauss,2.0);
	    }

	  error_egauss=(double)(b-1)/((double) b) * error_egauss;
	  error_egauss=sqrt(error_egauss);
	  fprintf(ferror,"%d %lf\n",n,error_egauss);  
	  printf(" %d %lf\n",n,error_egauss);  
	}
    }
  fclose(ferror);
      
  //GRAFICA GNUPLOT ERROR

  fprintf(pipegp,"set terminal push\n");//Almacena el anterior terminal gp

  //Guardamos el gráfico en formato latex en el disco:
  fprintf(pipegp,"set terminal epslatex color colortext\n");
  fprintf(pipegp,"set output \"./MEDIDAS_Egauss/L%d/K%.1f/Ploterror_Egauss-L%d-K%d.tex\" \n",L,K,L,(int)(10*K));
  fprintf(pipegp, "set title \' Error vs Tamaño Jacknife ($L=%d$ $K=%.1f$) \' \n",L,K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel \'Tamaño del bloque Jacknife (sweeps/$\\tau_0$)\' \n");
  fprintf(pipegp, "set ylabel \' Error $E_g$ \'\n");
  fprintf(pipegp, "plot \"%s\" title \'Error $E_g$\' w lp\n",nameerror);

  //Vuelve a la anterior terminal gnuplot:
  fprintf(pipegp,"set output\n");
  fprintf(pipegp,"set terminal pop\n");

  //Gráfico en pantalla:  
  fprintf(pipegp, "set title \" Error vs. tamaño del bloque Jacknife L=%d K=%.1f\" \n",L, K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel\" tamaño del bloque Jacknife (sweeps/tau)\"\n");
  fprintf(pipegp, "set ylabel\"Error energía elástica \"\n");  
  fprintf(pipegp, "plot \"%s\" title \"E. elástica\" w lp\n",nameerror);


  //Vacía el buffer de la tubería gnuplot y se cierra:
  fflush(pipegp);
  //close(pipegp);

  
  // Radio de la energía elástica vs. Valor del error:

  printf("\n Escribe el tamaño del bloque Jacknife en donde se estabiliza el error=");  
  while(scanf("%d",&nbloq)==0 || nbloq<0)
    {
      while (getchar()!= '\n');// para leer un único dato por línea
      printf("\n Debe ser un número entero positivo="); 
    }
  
  ferror=fopen(nameerror,"r");//Abrimos el archivo de errores para lectura
  while(nbloq!=n)
    {
      if(fscanf(ferror,"%d %lf",&n,&error_egauss)==EOF)
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

  printf("\n -> Rg2=%lf +- %lf\n", media_egauss, error_egauss);

  //Escribe el valor de la medida con su error en el archivo Medidas_Rg2
  sprintf(namegauss,"./MEDIDAS_Egauss/L%d/Medidas_Egauss_L%d.dat",L,L);
  filegauss=fopen(namegauss,"a");
  fprintf(filegauss,"%lf %lf %lf %d %d\n", K, media_egauss, error_egauss,indtermal,nbloq);
  close(filegauss);
  close(pipegp);
  return 1;
}

double energia_gauss(void )
{

  int dir,s1,s2,i;
  double norma2,sum_norma2,egauss;
  vector r;

  sum_norma2=0.0F;
  for(dir=0; dir<3; dir++)
    {
      for(s1=rombov1[dir].min.s1; s1<rombov1[dir].max.s1; s1++)
	{
	  for(s2=rombov1[dir].min.s2; s2<rombov1[dir].max.s2; s2++)
	    {
	      i=s1+L*s2;
	      _resta(r,x[v1[i][dir]],x[i]); 
	      norma2=_norma2(r);
	      sum_norma2+=norma2;
	    }
	}
    }
  return egauss=sum_norma2;
}

int sigma_min(int si)
{
  int simin;

  if(si<0)
    simin=-si;
  else
    simin=0;
  return simin;
}

int sigma_max(int si)
{
  int simax;

  if(si>0)
    simax=L-si;
  else
    simax=L;
  return simax;
}

void indice_vecnos_prox()
{
  int s1,s2;
  int s1min,s1max,s2min,s2max;
  int dir;
  int s1_dir,s2_dir,s1_old,s2_old;
  int i;
  int vec_new,vec_old;
  int ord;
 
  for(s1_dir=1,s2_dir=0,dir=0; dir<6; dir++)
    {
      basev1[dir].s1=s1_dir;
      basev1[dir].s2=s2_dir;
      
      s1_old=s1_dir;
      s2_old=s2_dir;
      
      s1_dir=-s2_old;
      s2_dir= s1_old + s2_old;
    }
  
  for(dir=0; dir<6; dir++)
    {
      rombov1[dir].min.s1=sigma_min(basev1[dir].s1);
      rombov1[dir].max.s1=sigma_max(basev1[dir].s1);
      rombov1[dir].min.s2=sigma_min(basev1[dir].s2);
      rombov1[dir].max.s2=sigma_max(basev1[dir].s2);
    }

  // Cálculo de v1[N][6],zv1[N]

  for(s2=0; s2<L; s2++)
    {
      for(s1=0; s1<L; s1++)
	{
	  i=s1+L*s2;
	  vec_old=NO;
	  zv1[i]=0;
	  dirv1_ini[i]=0;
	  dirv1_fin[i]=0;
	  for(dir=0; dir<6; dir++)
	    {
	      if(s1>=rombov1[dir].min.s1 && s1<rombov1[dir].max.s1 && s2>=rombov1[dir].min.s2 && s2<rombov1[dir].max.s2)
		{
		  v1[i][dir] = i + basev1[dir].s1 + L * basev1[dir].s2;
		  zv1[i]+=1;
		  vec_new=SI;
		}
	      else
		{
		  v1[i][dir]=-1;
		  vec_new=NO;
		}
	      if(vec_old-vec_new<0)
		dirv1_ini[i]=dir;
	      if(vec_old-vec_new>0)
		dirv1_fin[i]=dir-1;
	      vec_old=vec_new; 
	    }
	  if(vec_old==SI && v1[i][0]==-1)//Comparamos dir 0 con dir 5
		dirv1_fin[i]=5;
	}
    }
}
