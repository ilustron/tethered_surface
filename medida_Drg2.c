 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double radio_giraton(void);


double energia_curvatura(void);
void indice_vecnos_prox();

int sigma_min(int );
int sigma_max(int );

int index_plqta_dir0(int i);
int index_plqta_dir1(int i);
int index_plqta_prox(int i, int dir);

int Drg2_file(void);
void plot_historico(void);
void plot_termalizacion(void);
void error_Jacknife(int );

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

double se[NF];// Energía de curvatura correspondiente a cada configuración
double Se[NF];// Valor acumulado de la energía de curvatura 

double rg2[NF];//observable 
double Rg2[NF];//valor acumulado  


double rg2se[NF];
double Rg2Se[NF];

double drg2;
double error_drg2;
int main(void )
{
  FILE *input;
  FILE *output;

  char drg2data[255];
  char namein[255];
  char nameout[255];

  //Indices y contadores enteros
  int i,dir,k,f,n;
  int indtermal;
  int nbloq;

  //CÁLCULO de los ÍNDICES VÉCINOS
  indice_vecnos_prox(); // Cargamos los índices v1[N][6]

  for(i=0; i<N; i++)  // Cargamos los índices p[N][6]
    {
      for(dir=0; dir<6; dir++)
	{
	  p[i][dir]=index_plqta_prox(i,dir);
	}
    }

  
  //MEDIDA/LECTURA del Radio del giratón al para cada configuración:
  sprintf(drg2data,"./MEDIDAS_Drg2/L%d/K%.1f/Drg2_L%d_K%.1f.dat",L,K,L,K);
  if((input=fopen(drg2data,"r"))==NULL)
    {
      if(Drg2_file()==1)// 
	return 1;
    }
  else
    {
      i=0;      
      while(fscanf(input,"%lf %lf %lf %lf %lf %lf %lf",&rg2[f],&Rg2[f],&se[f],&Se[f],&rg2se[f],&Rg2Se[f],&drg2)!=EOF)
	{
	  i++;
	}
      fclose(input);     
      if(i!=NF)
	{
	  printf("Error lineal: El fichero %s no contiene %d líneas\n",drg2data,NF);
	  return 1;
      	}
    }

  //GRÁFICO Histórico
  plot_historico();

  //TERMALIZACION: 
  plot_termalizacion();

  //LECTURA del índice del archivo en donde comienza la termalización
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
  sprintf(namein,"./MEDIDAS_Drg2/L%d/K%.1f/error_Drg2_L%d_K%.1f.dat",L,K,L,K); //Abrimos el archivo de errores para lectura

  input=fopen(namein,"r");
  n=0;
  while(nbloq!=n)
    {
      if(fscanf(input,"%d %lf",&n,&error_drg2)==EOF)
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
  printf("\n -> Drg2=%lf +- %lf\n", drg2, error_drg2);

 //ESCRITURA RESULTADOS: 
  sprintf(nameout,"./MEDIDAS_Drg2/L%d/Medidas_Drg2_L%d.dat",L,L);
  output=fopen(nameout,"a");
  fprintf(output,"%lf %lf %lf %d %d\n", K, drg2, error_drg2,NF-indtermal,nbloq);
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
  energia/=((double) N);
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
int Drg2_file(void )
{
  int f,i;
  FILE *input;
  FILE *output;
  char namein[255];
  char nameout[255];

  f=0;
  Se[NF-1]=0.0F;
  Rg2[NF-1]=0.0F;
  Rg2Se[NF-1]=0.0F;
  sprintf(nameout,"./MEDIDAS_Drg2/L%d/K%.1f/Drg2_L%d_K%.1f.dat",L,K,L,K);
  output=fopen(nameout,"w");
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
      se[f]=energia_curvatura();
      Se[NF-1]+=se[f];
      Se[f]=Se[NF-1];

      rg2[f]=radio_giraton();
      Rg2[NF-1]+=rg2[f];
      Rg2[f]=Rg2[NF-1];

      rg2se[f]=rg2[f]*se[f];
      Rg2Se[NF-1]+=rg2se[f];
      Rg2Se[f]=Rg2Se[NF-1];

      drg2=Rg2Se[f]/(double)(f+1)-Se[f]*Rg2[f]/(double)(f+1)/(double)(f+1);

      // Escribimos los resultados en el archivo
      fprintf(output,"%lf %lf %lf %lf %lf %lf %lf \n",rg2[f],Rg2[f],se[f],Se[f],rg2se[f],Rg2Se[f],drg2);
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

  sprintf(nameout,"./MEDIDAS_Drg2/L%d/K%.1f/Drg2_L%d_K%.1f.dat",L,K,L,K);
  fprintf(pipegp,"set terminal qt 1\n");//Guarda la terminal gp por defecto
  fprintf(pipegp,"set terminal push\n");//Guarda la terminal gp por defecto

  //Guarda el gráfico en formato latex en el disco:
  fprintf(pipegp,"set terminal epslatex color colortext\n");
  fprintf(pipegp,"set output \"./MEDIDAS_Drg2/L%d/K%.1f/Phistorico_Drg2-L%d-K%d.tex\" \n",L,K,L,(int)(10*K));
  fprintf(pipegp, "set title \' Histórico medidas ($L=%d$ $K=%.1f$) \' \n",L,K);
  fprintf(pipegp, "set xlabel \'Sweeps\' \n");
  fprintf(pipegp, "set ylabel \' $\\langle DR_g^2 \\rangle$ \'\n");
  fprintf(pipegp, "plot \"%s\" u 1 title \'$Drg2$\' w lp\n",nameout);
  
  //vuelve a la anterior terminal gnuplot:
  fprintf(pipegp,"set output\n");
  fprintf(pipegp,"set terminal pop\n");

  //gráfico en pantalla:
  fprintf(pipegp, "set title \" Histórico medidas L=%d K=%.1f\" \n",L,K);
  fprintf(pipegp, "set xlabel\" Sweeps\"\n");
  fprintf(pipegp, "set ylabel\" Drg2 \"\n");
  fprintf(pipegp, "plot \"%s\"  u 7 title \"Se\" w lp \n",nameout);
    
  fflush(pipegp);//vacía el buffer de la tubería gnuplot:
  close(pipegp); 
}
void plot_termalizacion(void)
{
  int f;
  double media_se,media_rg2,media_drg2;
  FILE *pipegp = popen("gnuplot -persist","w");// Tubería a gnuplot (gráficas)  
  FILE *output;
  
  char nameout[255];

  sprintf(nameout,"./MEDIDAS_Drg2/L%d/K%.1f/logtermal_Drg2_L%d_K%.1f.dat",L,K,L,K);
  output=fopen(nameout,"w");
  printf("\n Plot: Termalización L=%d K=%.1f \n",L,K);
  printf(" X=1/(nº archivos conservados)\n");
  printf(" Y=promedio de la energía de curvatura\n");
  printf(" Ind.=índice del primer archivo conservado \n");
  printf("\n Ind. X  Y\n");
  for(f=NF/20; f<NF; f++)
    {      
      if((NF%f)==0)
      {
	media_se=(Se[NF-1]-Se[NF-f-1])/(double) f;
	media_rg2=(Rg2[NF-1]-Rg2[NF-f-1])/(double) f;
	media_drg2=(Rg2Se[NF-1]-Rg2Se[NF-f-1])/(double)f-media_se*media_rg2;
	fprintf(output,"%lf %lf\n",1.0F/(double) f,media_drg2);
	printf(" %d %lf %lf\n",NF-f,1.0F/(double) f,media_drg2);
      }
    }

  media_se=Se[NF-1]/(double)NF;
  media_rg2=Rg2[NF-1]/(double) NF;
  media_drg2=Rg2Se[NF-1]/(double)NF-media_se*media_rg2;
  fprintf(output,"%lf %lf\n",1.0F/(double)NF,media_drg2);
  printf(" %d %lf %lf\n",0,1.0F/(double) f,media_drg2);
  fclose(output);

  //GRÁFICA GNUPLOT TERMALIZACION:

  fprintf(pipegp,"set terminal qt 2\n");
  fprintf(pipegp,"set terminal push\n");//Guarda la terminal gp por defecto

  //Guarda el gráfico en formato latex en el disco:
  fprintf(pipegp,"set terminal epslatex color colortext\n");
  fprintf(pipegp,"set output \"./MEDIDAS_Drg2/L%d/K%.1f/Plogtermal_Drg2-L%d-K%d.tex\" \n",L,K,L,(int)(10*K));
  fprintf(pipegp, "set title \' Gráfico Termalización ($L=%d$ $K=%.1f$) \' \n",L,K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel \'1/ (n\\textdegree de archivos conservados)\' \n");
  fprintf(pipegp, "set ylabel \' $\\langle DR_g^2 \\rangle$ \'\n");
  fprintf(pipegp, "plot \"%s\" title \'$DR^2_g$\' w lp\n",nameout);
  
  //vuelve a la anterior terminal gnuplot:
  fprintf(pipegp,"set output\n");
  fprintf(pipegp,"set terminal pop\n");

  //gráfico en pantalla:
  fprintf(pipegp, "set title \" Termalización L=%d K=%.1f\" \n",L,K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel\" nº de archivos conservados\"\n");
  fprintf(pipegp, "set ylabel\" promedio Drg2 \"\n");
  fprintf(pipegp, "plot \"%s\" title \"Drg2\" w lp\n",nameout);
    
  fflush(pipegp);//vacía el buffer de la tubería gnuplot:
  close(pipegp);
}
void error_Jacknife(int index)
{

  int b,n,f,k;
  int fmax;

  double seJK[NF];
  double sumJK_se;

  double rg2JK[NF];
  double sumJK_rg2;

  double rg2seJK[NF];
  double sumJK_rg2se;

  double drg2JK;

  double media_se;
  double media_rg2;

  double Senew[NF];
  double Rg2new[NF];
  double Rg2Senew[NF];
  double Seref;
  double Rg2ref;
  double Rg2Seref;

  FILE *output;// Archivo de escritura = Error vs. Tamaño del bloque Jacknife 
  char nameout[255];

  FILE *pipegp = popen("gnuplot -persist","w");// Tubería a gnuplot (gráficas)

  fmax=NF-index;//fmax es ahora el nº total de archivos para los cálculos
 
  
  Seref=0.0F;
  Rg2ref=0.0F;
  Rg2Seref=0.0F;   
  if(fmax<NF)
  {
    Seref=Se[index-1];
    Rg2ref=Rg2[index-1];
    Rg2Seref=Rg2Se[index-1];
  }
  

  for(f=0; f<fmax; f++) 
    {
      Senew[f]=Se[f+index]-Seref;
      Rg2new[f]=Rg2[f+index]-Rg2ref;
      Rg2Senew[f]=Rg2Se[f+index]-Rg2Seref;
    }   
  
  media_se=Senew[fmax-1]/(double)fmax;
  media_rg2=Rg2new[fmax-1]/(double)fmax;
  drg2=Rg2Senew[fmax-1]/(double)fmax-media_se*media_rg2;
  sprintf(nameout,"./MEDIDAS_Drg2/L%d/K%.1f/error_Drg2_L%d_K%.1f.dat",L,K,L,K);
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
	  error_drg2=0.0F;

	  seJK[0]=Senew[n-1];
	  rg2JK[0]=Rg2new[n-1];
	  rg2seJK[0]=Rg2Senew[n-1];

	  sumJK_se=(Senew[fmax-1]-seJK[0])/(double)(fmax-n);
	  sumJK_rg2=(Rg2new[fmax-1]-rg2JK[0])/(double)(fmax-n);
	  sumJK_rg2se=(Rg2Senew[fmax-1]-rg2seJK[0])/(double)(fmax-n);
	  
	  drg2JK=sumJK_rg2se-sumJK_se*sumJK_rg2;
	  error_drg2+=pow(drg2-drg2JK,2.0);
	  
	  for (k=1; k<b; k++)
	    {
	      seJK[k]=Senew[(k+1)*n-1]-Senew[k*n-1];
	      rg2JK[k]=Rg2new[(k+1)*n-1]-Rg2new[k*n-1];
	      rg2seJK[k]=Rg2Senew[(k+1)*n-1]-Rg2Senew[k*n-1];

	      sumJK_se=(Senew[fmax-1]-seJK[k])/(double)(fmax-n);
	      sumJK_rg2=(Rg2new[fmax-1]-rg2JK[k])/(double)(fmax-n);
	      sumJK_rg2se=(Rg2Senew[fmax-1]-rg2seJK[k])/(double)(fmax-n);

	      drg2JK=sumJK_rg2se-sumJK_se*sumJK_rg2;
	      error_drg2+=pow(drg2-drg2JK,2.0);
	    }

	  error_drg2=(double)(b-1.0F)/((double) b) * error_drg2;
	  error_drg2=sqrt(error_drg2);
	  fprintf(output,"%d %lf\n",n,error_drg2); 
	  printf(" %d %lf\n",n,error_drg2);  
	}
    }
  fclose(output);
      
  //GRAFICA GNUPLOT ERROR
  fprintf(pipegp,"set terminal qt 3\n");
  fprintf(pipegp,"set terminal push\n");//Almacena el anterior terminal gp

  //Guardamos el gráfico en formato latex en el disco:
  fprintf(pipegp,"set terminal epslatex color colortext\n");
  fprintf(pipegp,"set output \"./MEDIDAS_Drg2/L%d/K%.1f/Perror_Drg2-L%d-K%d.tex\" \n",L,K,L,(int)(10*K));
  fprintf(pipegp, "set title \' Error vs Tamaño Jacknife ($L=%d$ $K=%.1f$) \' \n",L,K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel \'Tamaño del bloque Jacknife (sweeps/$\\tau_0$)\' \n");
  fprintf(pipegp, "set ylabel \' Error $Drg2$ \'\n");
  fprintf(pipegp, "plot \"%s\" title \'Error $Drg2$\' w lp\n",nameout);

  //Vuelve a la anterior terminal gnuplot:
  fprintf(pipegp,"set output\n");
  fprintf(pipegp,"set terminal pop\n");

  //Gráfico en pantalla:  
  fprintf(pipegp, "set title \" Error vs. tamaño del bloque Jacknife L=%d K=%.1f\" \n",L, K);
  fprintf(pipegp, "set logscale x\n");
  fprintf(pipegp, "set xlabel\" tamaño del bloque Jacknife (sweeps/tau)\"\n");
  fprintf(pipegp, "set ylabel\"Error Drg2 \"\n");  
  fprintf(pipegp, "plot \"%s\" title \"Drg2\" w lp\n",nameout);


  //Vacía el buffer de la tubería gnuplot y se cierra:
  fflush(pipegp);
  close(pipegp);
}
