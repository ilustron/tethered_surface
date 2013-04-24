//COMPILAR CON:
//gcc -O2 -DK=0.9 -DL=16 -DNF=1000 -DTAU=16000 -DTERMAL=8000000 medida_Cv.c -lm -o medida_Cv.out 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double energia_curvatura(void);

void indice_vecnos_prox();

int sigma_min(int );
int sigma_max(int );


int index_plqta_dir0(int i);
int index_plqta_dir1(int i);
int index_plqta_prox(int i, int dir);

#define N L*L 
#define M 2*(L-1)*(L-1) 

#define _prodesc(w,u) (w.a*u.a+w.b*u.b+w.c*u.c)

#define _escala(u,x) \
                         u.a *= (x);\
                         u.b *= (x);\
                         u.c *= (x);

#define _norma2(u)  (u.a*u.a + u.b*u.b + u.c*u.c)

#define _prodvec(u,v,w) \
                         u.a = v.b*w.c - v.c*w.b; \
                         u.b = v.c*w.a - v.a*w.c; \
                         u.c = v.a*w.b - v.b*w.a; 
#define _resta(u,v,w) \
                         u.a = v.a - w.a; \
                         u.b = v.b - w.b; \
                         u.c = v.c - w.c; 
#define _suma(u,v,w) \
			u.a = v.a + w.a; \
                        u.b = v.b + w.b; \
                        u.c = v.c + w.c; 	
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
int p[N][6]; // Indices de las plaquetas vecinas

int main(void )
{
  FILE *input;// -> Archivo de lectura = Posiciones 3d de los nodos 
  FILE *ftermal;// -> Archivo de escritura = Termalización
  FILE *ferror;// -> Archivo de escritura = Error vs. Tamaño del bloque Jacknife 
  FILE *fileCv;// -> Archivo de escritura = Valores Rg2 con errores
  FILE *pipegp = popen("gnuplot -persist","w");//Tubería a gnuplot (gráficas)

  char namein[255];
  char nametermal[255];
  char nameerror[255];
  char nameCv[255];

  //OBSERVABLES
  double caloresp;//calor específico
  double energia[NF];// Energía de curvatura correspondiente a cada configuración
  double Energia[NF];// Valor acumulado de la energía de curvatura 
  double media_energia;// media de la energia
  double Energia2[NF];// Valor acumulado de la energía de curvatura al cuadradado

  //Valores-bloque Jacknife
  double energiaJK[NF];
  double sumJK_energia;

  double energia2JK[NF];
  double sum2JK_energia;
  
  double calorespJK;

  //Error
  double error_caloresp;

  int i,k,dir,f,b,n;

  int indtermal;
  int fmax;

  int nbloq;

  //TOPOLOGÍA:

  indice_vecnos_prox(); // Cargamos los índices v1[N][6]

  for(i=0; i<N; i++)  // Cargamos los índices p[N][6]
    {
      for(dir=0; dir<6; dir++)
	{
	  p[i][dir]=index_plqta_prox(i,dir);
	}
    }
 

  // LECTURA ARCHIVOS FUENTE (Correspondientes a las posiciones de los nodos de la membrana)
 
  f=0;
  Energia[NF-1]=0.0F;
  
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
      energia[f]=energia_curvatura();
      Energia[NF-1]+=energia[f];
      Energia[f]=Energia[NF-1];

      Energia2[NF-1]+=(energia[f]*energia[f]);
      Energia2[f]=Energia2[NF-1];
      f++;
      sprintf(namein,"./RUNS/L%d/K%.1f/xpos_L%d_K%.1f-%d.dat",L,K,L,K,f);
    }
  
  if(f!=NF)
    {
      printf("Error:El número de archivos leidos en ./RUNS/L%d/K%.1f/ no es %d\n",L,K,NF); 
      return 0;
    }

  //TERMALIZACION:
  sprintf(nametermal,"./MEDIDAS_CV/L%d/K%.1f/termalizacion_Cv_L%d_K%.1f.dat",L,K,L,K);
  ftermal=fopen(nametermal,"w");
  printf("\n Plot: Termalización L=%d K=%.1f \n",L,K);
  printf(" X=1/(nº archivos conservados)\n");
  printf(" Y=promedio del calor específico\n");
  printf(" Ind.=índice del primer archivo conservado \n");
  printf("\n Ind. X  Y\n");
  for(f=NF/20; f<NF; f++)
    {      
      if((NF%f)==0)
      {
	media_energia=(Energia[NF-1]-Energia[NF-f-1])/(double) f;
	caloresp=(Energia2[NF-1]-Energia2[NF-f-1])/(double)f-media_energia*media_energia;
	caloresp*=(K*K)/((double)N); // promedio calor especifico
	fprintf(ftermal,"%lf %lf\n",1.0F/(double) f,caloresp);
	printf(" %d %lf %lf\n",NF-f,1.0F/(double) f,caloresp);
      }
    }
  media_energia=Energia[NF-1]/(double) NF;
  caloresp=Energia2[NF-1]/(double)NF-media_energia*media_energia;
  caloresp*=(K*K)/((double)N); // promedio calor especifico
  fprintf(ftermal,"%lf %lf\n",1.0F/(double) NF,caloresp);
  printf(" %d %lf %lf\n",0,1.0F/(double) NF,caloresp);
  fclose(ftermal);

  //GRÁFICA GNUPLOT TERMALIZACION:

  fprintf(pipegp,"set terminal push\n");//Guarda la terminal gp por defecto

  //Guarda el gráfico en formato latex en el disco:
  fprintf(pipegp,"set terminal epslatex color colortext\n");
  fprintf(pipegp,"set output \"./MEDIDAS_CV/L%d/K%.1f/Plottermal_Cv-L%d-K%d.tex\" \n",L,K,L,(int)(10*K));
  fprintf(pipegp, "set title \' Gráfico Termalización ($L=%d$ $K=%.1f$) \' \n",L,K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel \'1/ (n\\textdegree de archivos conservados)\' \n");
  fprintf(pipegp, "set ylabel \' $\\langle C_v \\rangle$ \'\n");
  fprintf(pipegp, "plot \"%s\" title \'$C_v$\' w lp\n",nametermal);
  
  //vuelve a la anterior terminal gnuplot:
  fprintf(pipegp,"set output\n");
  fprintf(pipegp,"set terminal pop\n");

  //gráfico en pantalla:
  fprintf(pipegp, "set title \" Termalización L=%d K=%.1f\" \n",L,K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel\" nº de archivos conservados\"\n",TAU);
  fprintf(pipegp, "set ylabel\" promedio Cv \"\n");
  fprintf(pipegp, "plot \"%s\" title \"C_v\" w lp\n",nametermal);
    
  fflush(pipegp);//vacía el buffer de la tubería gnuplot
  

  //Entrada del valor del índice del archivo en el que se produce la termalización:
  printf("\n Escribe el valor del índice del archivo correspondiente a la termalización=");
 
  while(scanf("%d",&indtermal)==0 || indtermal<0)
    {
      while (getchar()!= '\n');// para leer un único dato por línea
      printf("\n Debe ser un número entero positivo="); 
    }

  fmax=NF-indtermal;//fmax es ahora el nº total de archivos para los cálculos
  if(fmax!=NF)
    {
      for(f=0; f<fmax; f++)//redefinimos los observables
	{
	  Energia[f]=Energia[f+indtermal]-Energia[indtermal-1];
	  Energia2[f]=Energia2[f+indtermal]-Energia2[indtermal-1];
	}   
      media_energia=Energia[fmax-1]/(double) fmax;
      caloresp=Energia2[fmax-1]/(double)fmax-media_energia*media_energia;
      caloresp*=(K*K)/((double)N); // promedio calor especifico
    }
  printf("\n -> Media Cv= %lf\n",caloresp);
 
  // Error en función del nº de bloques

  sprintf(nameerror,"./MEDIDAS_CV/L%d/K%.1f/error_Cv_L%d_K%.1f.dat",L,K,L,K);
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

	  error_caloresp=0.0F;

	  energiaJK[0]=Energia[n-1];
	  energia2JK[0]=Energia2[n-1];

	  sumJK_energia=(Energia[fmax-1]-energiaJK[0])/(double)(fmax-n);
	  sum2JK_energia=(Energia2[fmax-1]-energia2JK[0])/(double)(fmax-n);

	  calorespJK=sum2JK_energia-sumJK_energia*sumJK_energia;
	  calorespJK*=(K*K)/((double)N);
	  error_caloresp+=pow(caloresp-calorespJK,2.0);

	  for (k=1; k<b; k++)
	    {
	      energiaJK[k]=Energia[(k+1)*n-1]-Energia[k*n-1];
	      energia2JK[k]=Energia2[(k+1)*n-1]-Energia2[k*n-1];
	      sumJK_energia=(Energia[fmax-1]-energiaJK[k])/(double)(fmax-n);
	      sum2JK_energia=(Energia2[fmax-1]-energia2JK[k])/(double)(fmax-n);

	      calorespJK=sum2JK_energia-sumJK_energia*sumJK_energia;
	      calorespJK*=(K*K)/((double)N);
	      error_caloresp+=pow(caloresp-calorespJK,2.0);
	    }

	  error_caloresp=(double)(b-1.)/((double) b) * error_caloresp;
	  error_caloresp=sqrt(error_caloresp);
	  fprintf(ferror,"%d %lf\n",n,error_caloresp); 
	  printf(" %d %lf\n",n,error_caloresp);  
	}
    }
  fclose(ferror);
      
  //GRAFICA GNUPLOT ERROR: Error en función del tamaño del bloque Jacknife

  fprintf(pipegp,"set terminal push\n");//Almacena el anterior terminal gp
  
  //Guardamos el gráfico en formato latex en el disco:
  fprintf(pipegp,"set terminal epslatex color colortext\n");
  fprintf(pipegp,"set output \"./MEDIDAS_CV/L%d/K%.1f/Ploterror_Cv-L%d-K%d.tex\" \n",L,K,L,(int)(10*K));
  fprintf(pipegp, "set title \' Error vs Tamaño Jacknife ($L=%d$ $K=%.1f$) \' \n",L,K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel \'Tamaño del bloque Jacknife (sweeps/$\\tau_0$)\' \n");
  fprintf(pipegp, "set ylabel \' Error $C_v$ \'\n");
  fprintf(pipegp, "plot \"%s\" title \'Error $C_v$\' w lp\n",nameerror);

  //Vuelve a la anterior terminal gnuplot:
  fprintf(pipegp,"set output\n");
  fprintf(pipegp,"set terminal pop\n");

  //Gráfico en pantalla:  
  fprintf(pipegp, "set title \" Error vs. tamaño del bloque Jacknife L=%d K=%.1f\" \n",L, K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel\" tamaño del bloque Jacknife (sweeps/tau)\"\n");
  fprintf(pipegp, "set ylabel\"Error Cv \"\n");  
  fprintf(pipegp, "plot \"%s\" title \"Cv\" w lp\n",nameerror);

  //Vacía el buffer de la tubería gnuplot y se cierra:
  fflush(pipegp);
  //close(pipegp);
  

  //Lectura valor del donde estabiliza el error

  printf("\n Escribe el tamaño del bloque Jacknife en donde se estabiliza el error=");  
  while(scanf("%d",&nbloq)==0 || nbloq<0)
    {
      while (getchar()!= '\n');// para leer un único dato por línea
      printf("\n Debe ser un número entero positivo="); 
    }
  
  ferror=fopen(nameerror,"r");//Abrimos el archivo de errores para lectura
  while(nbloq!=n)
    {
      if(fscanf(ferror,"%d %lf",&n,&error_caloresp)==EOF)
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
 
  printf("\n ->Cv=%lf +- %lf\n", caloresp, error_caloresp);

  //Escribe el valor de la medida con su error en el archivo Medidas_Cv..
  sprintf(nameCv,"./MEDIDAS_CV/L%d/Medidas_Cv_L%d.dat",L,L);
  fileCv=fopen(nameCv,"a");
  fprintf(fileCv,"%lf %lf %lf\n", K, caloresp, error_caloresp,indtermal,nbloq);
  close(fileCv);
  close(pipegp);
  return 1;
}

double energia_curvatura(void)
{
  
  int i,j,k;
  int s1,s2;
  double energia;
  
  vector r1,r2,r3,ntemp,n[M];
  double norma2,invnorma;
 
  energia=0.0F;

  for(s1=1; s1<L; s1++)
    {
      for(s2=0; s2<L-1; s2++)
	{

	  i=s1+L*s2;

	  //Vectores relativos:

	  _resta(r1,x[v1[i][1]],x[i]);
	  _resta(r2,x[v1[i][2]],x[i]);
	  _resta(r3,x[v1[i][3]],x[i]);

	  //Productos vectoriales:

	   _prodvec(ntemp,r1,r2);
	   norma2=_norma2(ntemp);
	   invnorma=1.0F/sqrt(norma2);
	   _escala(ntemp,invnorma);
	   n[p[i][1]]=ntemp;

	   _prodvec(ntemp,r2,r3);
	   norma2=_norma2(ntemp);
	   invnorma=1.0F/sqrt(norma2);
	   _escala(ntemp,invnorma);
	   n[p[i][2]]=ntemp;
	}
    }

  // Productos dirección 1

  for(s1=1; s1<L-1; s1++)
    {
      for(s2=0; s2<L-1; s2++)
	{
	  i=s1+L*s2;
	  energia+=_prodesc(n[p[i][1]],n[p[i][0]]);
	}      
    }
  // Productos dirección 2

  for(s1=0; s1<L; s1++)
    {
      for(s2=0; s2<L-1; s2++)
	{
	  i=s1+L*s2;
	  energia+=_prodesc(n[p[i][2]],n[p[i][1]]);
	}      
    }
  // Productos dirección 3

  for(s1=1; s1<L; s1++)
    {
      for(s2=1; s2<L-1; s2++)
	{
	  i=s1+L*s2;
	  energia+=_prodesc(n[p[i][3]],n[p[i][2]]);
	}      
    }

    return energia;
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
/*INDEX_PLQTA_DIR1(int i):
*  Devuelve el índice de la plaqueta vecina del nodo i en la dirección 1
*   * Las direcciones de las plaquetas vecinas se muestran en la figura:
*
* 		           x-----x
* 		          / \ 1 / \
* 		         / 2 \ / 0 \
* 		        x-----i-----x
* 		         \ 3 / \ 5 /
* 		          \ / 4 \ /
* 		           x-----x
*
*/
int index_plqta_dir1(int i)
{
  int p;

  if(v1[i][1]!=-1 && v1[i][2]!=-1)
    p=2*i - (2* (i/L) + 1);
  else
    {
      p=-1;
    }
return p;
}

/*INDEX_PLQTA_DIR0(int i):
*  Devuelve el índice de la plaqueta vecina del nodo i en la dirección 0
*   * Las direcciones de las plaquetas vecinas se muestran en la figura:
*
* 		           x-----x
* 		          / \ 1 / \
* 		         / 2 \ / 0 \
* 		        x-----i-----x
* 		         \ 3 / \ 5 /
* 		          \ / 4 \ /
* 		           x-----x
*
*/
int index_plqta_dir0(int i)
{
  int p;

  if(v1[i][0]!=-1 && v1[i][1]!=-1)
    p=2*i - (2* (i/L));
  else
    {
      p=-1;
    }
  return p;
}


/*INDEX_PLQTA_PROX(int i int dir):
*  Devuelve el índice de la plaqueta vecina del nodo i en la dirección dir
*   * Las direcciones de las plaquetas vecinas se muestran en la figura:
*
* 		           x-----x
* 		          / \ 1 / \
* 		         / 2 \ / 0 \
* 		        x-----i-----x
* 		         \ 3 / \ 5 /
* 		          \ / 4 \ /
* 		           x-----x
*
*/
int index_plqta_prox(int i, int dir)
{
  int p;

  switch(dir)
    {
    case 0:
      p=index_plqta_dir0(i);
      break;
    case 1:
      p=index_plqta_dir1(i);
      break;
    case 2:
      p=index_plqta_dir0(v1[i][3]);
      break;
    case 3:
      p=index_plqta_dir1(v1[i][4]);
      break;
    case 4:
      p=index_plqta_dir0(v1[i][4]);
      break;
    case 5:
      p=index_plqta_dir1(v1[i][5]);
      break;
    }
  return p;
}

