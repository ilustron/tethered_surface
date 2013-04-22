//COMPILAR CON:
//gcc -O2 -DK=0.9 -DL=16 -DNF=1000 -DTAU=16000 -DTERMAL=8000000 medida_Se.c -lm -o medida_Se.out 

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
  FILE *input;
  FILE *output;
  FILE *pipe;
  FILE *fileSe;

  char namein[255];
  char nameout[255];
  char nameSe[255];

  //OBSERVABLES
  double energia[NF];// Energía de curvatura correspondiente a cada configuración
  double Energia[NF];// Valor acumulado de la energía de curvatura 
  double media_energia;// media del radio del giraton al cuadrado

  //Valores -bloque Jacknife
  double energiaJK[NF];
  double sumJK_energia;

  //Error
  double error_energia;

  int i,k,dir,f,b,n;

  int ftermal,fmax;

  indice_vecnos_prox(); // Cargamos los índices v1[N][6]

  for(i=0; i<N; i++)  // Cargamos los índices p[N][6]
    {
      for(dir=0; dir<6; dir++)
	{
	  p[i][dir]=index_plqta_prox(i,dir);
	}
    }
 

  //Archivo LOG
  //sprintf(namelog,"./MEDIDAS_Rg2/L%d/K%.1f/medida_Rg2_L%d_K%.1f.log",L,K,L,K);
  //filelog=fopen(namelog,"w"); 
  //fprintf(filelog,"MEDIDA DEL RADIO DE GIRATÓN AL CUADRADO:\n");
  //fprintf(filelog,"L=%d N=%d M=%d K=%.1f \n",L,N,M,K);
  //fprintf(filelog,"NF=%d TAU=%d NTERMAL=%d\n",NF,TAU,NTERMAL);

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
      f++;
      sprintf(namein,"./RUNS/L%d/K%.1f/xpos_L%d_K%.1f-%d.dat",L,K,L,K,f);
    }
  
  if(f!=NF)
    {
      printf("Error:El número de archivos leidos en ./RUNS/L%d/K%.1f/ no es %d\n",L,K,NF); 
      return 0;
    }

  //TERMALIZACION:
  sprintf(nameout,"./MEDIDAS_Se/L%d/K%.1f/termalizacionSe_L%d_K%.1f.dat",L,K,L,K);
  output=fopen(nameout,"w"); 
  for(f=NF; f>0; f--)
    {      
      if((NF%f)==0)
      {
	      b=NF/f;
	      n=b-1;
	      media_energia=Energia[n]/(double)b;
	      fprintf(output,"%d %lf\n",b,media_energia);
      }
    }
  fclose(output);

  //GRÁFICA GNUPLOT TERMALIZACION:
  pipe = popen("gnuplot -persist","w");
  fprintf(pipe, "set title \" Termalización \" \n");
  fprintf(pipe, "set logscale x\n");
  fprintf(pipe, "set xlabel\" número de sweeps/(tau=%d)\"\n",TAU);
  fprintf(pipe, "set ylabel\" promedio energia acumulado\"\n");
  fprintf(pipe, "plot \"%s\" title \"Se\" w lp\n",nameout);
  fprintf(pipe,"set terminal push\n");
  fprintf(pipe,"set terminal png\n");
  fprintf(pipe,"set output \"./MEDIDAS_Se/L%d/K%.1f/termalizacionSe_L%d_K%.1f.png\" \n",L,K,L,K);
  fprintf(pipe,"replot\n");
  fprintf(pipe,"set output\n");
  fprintf(pipe,"set terminal pop\n");
  fflush(pipe);

  getchar(); //Pausa
  
  //Esta parte corresponde a introducir por teclado el valor observado para la termalización

  ftermal=0;

  /*printf("Escribe el valor de sweep/TAU correspondiente a la termalización:\n");
  
  while(scanf("%d",&ftermal)==0 || ftermal<0)
    {
      while (getchar()!= '\n');// para leer un único dato por línea
      printf("Debe ser un número entero positivo:\n"); 
    }

  //Escribimos el valor de la termalización oobservada en el archivo .log
  fprintf(filelog,"Termalización:\n");
  fprintf(filelog," nfile_termal=%d sweep_obs_termalizacion=%d\n",ftermal,ftermal*TAU+NTERMAL);*/
  

  if((fmax=NF-ftermal)!=NF)
    {
      for(f=0; f<fmax; f++)//redefinimos los observables
	{
	  Energia[f]=Energia[f+ftermal]-Energia[ftermal-1];
	}
      
      sprintf(nameout,"./MEDIDAS_Se/L%d/K%.1f/termalizacionSeNEW_L%d_K%.1f.dat",L,K,L,K);
      output=fopen(nameout,"w"); 
      for(f=fmax; f>0; f--)
	{
	  if((fmax%f)==0)
	    {
	      b=fmax/f;
	      n=b-1;
	      media_energia=Energia[n]/(double)b;
	      fprintf(output,"%d %lf\n",b,media_energia);
	    }
	}//El último valor que toma media_energia es el correcto
      fclose(output);
      
      pipe = popen("gnuplot -persist","w");
      fprintf(pipe, "set title \" Termalización \"\n");
      fprintf(pipe, "set logscale x\n");
      fprintf(pipe, "set xlabel\" número de sweeps/(tau=%d)\"\n",TAU);  
      fprintf(pipe, "set ylabel\" promedio energia \"\n");
      fprintf(pipe, "plot \"%s\" title \"Energia\" w lp\n",nameout);
      fprintf(pipe,"set terminal push\n");
      fprintf(pipe,"set terminal png\n");
      fprintf(pipe,"set output \"./MEDIDAS_Se/L%d/K%.1f/termalizacionSeNEW_L%d_K%.1f.png\" \n",L,K,L,K);
      fprintf(pipe,"replot\n");
      fprintf(pipe,"set output\n");
      fprintf(pipe,"set terminal pop\n");
      fflush(pipe);
      close(pipe);
    }
 

  // Error en función del nº de bloques

  sprintf(nameout,"./MEDIDAS_Se/L%d/K%.1f/error_Se_L%d_K%.1f.dat",L,K,L,K);
  output=fopen(nameout,"w"); 
  for(b=2; b<=fmax; b++)
    {
      if((fmax%b)==0)
	{
	  n=fmax/b;
	  error_energia=0.0F;

	  energiaJK[0]=Energia[n-1];

	  sumJK_energia=(Energia[fmax-1]-energiaJK[0])/(double)(fmax-n);
	  error_energia+=pow(sumJK_energia-media_energia,2.0);
	  
	  for (k=1; k<b; k++)
	    {
	      energiaJK[k]=Energia[(k+1)*n-1]-Energia[k*n-1];

	      sumJK_energia=(Energia[fmax-1]-energiaJK[k])/(double)(fmax-n);
	      error_energia+=pow(sumJK_energia-media_energia,2.0);
	    }

	  error_energia=(double)(b-1)/((double) b) * error_energia;
	  error_energia=sqrt(error_energia);
	  fprintf(output,"%d %lf\n",n,error_energia);  
	}
    }
  fclose(output);
      
  //GRAFICA GNUPLOT ERROR: Error en función del tamaño del bloque Jacknife

  pipe = popen("gnuplot -persist","w");  
  fprintf(pipe, "set title \" Error en función del tamaño del bloque Jacknife\" \n");
  fprintf(pipe, "set logscale x\n");
  fprintf(pipe, "set xlabel\" tamaño del bloque jacknife = sweeps/(tau=%d)\"\n",TAU);
  fprintf(pipe, "set ylabel\"error energia \"\n");  
  fprintf(pipe, "plot \"%s\" title \"Radio2 g\" w lp\n",nameout);
  fprintf(pipe,"set terminal push\n");
  fprintf(pipe,"set terminal png\n");
  fprintf(pipe,"set output \"./MEDIDAS_Rg2/L%d/K%.1f/errorRg2_L%d_K%.1f.png\" \n",L,K,L,K);
  fprintf(pipe,"replot\n");
  fprintf(pipe,"set output\n");
  fprintf(pipe,"set terminal pop\n");
  fflush(pipe);
  close(pipe);

  getchar();//pausa
  /*
  // Valor del error
  printf("Escribe el tamaño del bloque Jacknife en donde se estabiliza el error:\n");  
  while(scanf("%d",&nbloq)==0 || nbloq<0)
    {
      while (getchar()!= '\n');// para leer un único dato por línea
      printf("Debe ser un número entero positivo:\n"); 
    }
  
  input=fopen(nameout,"r");
  while(nbloq!=n)
    {
      if(fscanf(input,"%d %lf",&n,&error_energia)==EOF)
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

  //Escribimos en el archivo log los valores del tamaño del bloque Jacknife al que estabiliza el error
  fprintf(filelog,"Error:\n");
  fprintf(filelog,"tamaño bloque=%d tamaño bloque (sweep)=%d \n",nbloq,nbloq*TAU);
  close(filelog);
  */

  printf("media Se=%lf +- %lf\n", media_energia, error_energia);//El último valor del error corresponde a un tamaño de bloque 1

  //Escribe el valor de la medida con su error en el archivo valores_Rg2
  sprintf(nameSe,"./MEDIDAS_Se/L%d/valores_Se_L%d.dat",L,L);
  fileSe=fopen(nameSe,"a");
  fprintf(fileSe,"%lf %lf %lf\n", K, media_energia, error_energia);
  close(fileSe);

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

