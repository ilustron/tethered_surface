#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void indice_vecnos_prox();

int sigma_min(int );
int sigma_max(int );

#define K 2.0
#define L 16
#define N L*L 
#define M 2*(L-1)*(L-1) 
#define F 1

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
  FILE *input;
  FILE *output;

  char namein[255];
  char nameout[255];

  int s1,s2;
  int i;

  //CÁLCULO de los ÍNDICES VÉCINOS
  indice_vecnos_prox(); // Cargamos los índices v1[N][6]

  //Lectura de los archivos de posicion
  
  sprintf(namein,"./RUNS/L%d/K%.1f/xpos_L%d_K%.1f-%d.dat",L,K,L,K,F);
  if((input=fopen(namein,"r"))==NULL)
    {
      printf("Error existencial: El archivo %s no existe\n",namein); 
      return 0;
    }    
  i=0;      
  while(fscanf(input,"%lf %lf %lf",&x[i].a,&x[i].b,&x[i].c)!=EOF)
    {
      printf("%lf %lf %lf\n",x[i].a,x[i].b,x[i].c);
      i++;
    }
  fclose(input);     
  if(i!=N)
    {
      printf("Error: El fichero %s no contiene %d líneas\n",namein,N);
      return 0;
    }
  
  
  //Escritura del archivo para graficar
  sprintf(nameout,"snapshot_x_L%d_K%.1f-%d.dat",L,K,F);
  output=fopen(nameout,"w");
  for(s1=0; s1<L; s1++)//cuadrícula principal
    {
      for(s2=0; s2<L; s2++)
	{
	  i=s1+L*s2;
	  fprintf(output,"%lf %lf %lf\n",x[i].a,x[i].b,x[i].c);
	}
      fprintf(output,"\n");
    }
  fclose(output);
  return 1;
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
