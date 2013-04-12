/* FREESURF 2.0
*
*
* 		           X-----X-----X-----X-----X---(N-1)
* 		          / \   / \   /	\   / \   / \ M /
* 		         /   \ /   \ /	 \ /   \ /   \ /
* 		        X-----X-----X-----X-----X-----X
* 		       / \   / \   / \	 / \   / \   /
* 		      /   \ /   \ /   \ /   \ /   \ /
* 		     X-----X-----X-----X-----X-----X
*		    / \	  / \   / \   / \   / \   /
*		   /   \ /   \ /   \ /   \ /   \ /
*	          X-----X-----X-----X-----X-----X
*	         / 2L+1/ \   / \   / \	 / \   /
*	        /2L \ /   \ /   \ /   \	/   \ /
*	      (L)--(L+1)---X-----X-----X---(2L-1)
*	      / \ 1 / \	  / \   / \   / \2L-1
*	     / 0 \ /   \ /   \ /   \ /	 \ /
*	   (0)---(1)----X-----X-----X---(L-1)
*
*
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void regparam(void );

int sigma_min(int );
int sigma_max(int );
void indice_vecnos_prox(void);
void indice_vecnos_seg(void);
void indice_vecnos_ord(void);

int index_plqta_dir0(int i);
int index_plqta_dir1(int i);
int index_plqta_prox(int i, int dir);
int index_plqta_seg(int i, int dir);
void indice_plqta_ord(void);
void ini_random(unsigned long long semilla);

void files_topology(void);

void tipos_nodos(void);
void vectores_iniciales(void);
double sweep_spiral(double delta, double kappa);

#define SEMILLA 123421
#define KAPPA 0.5

#define L 16
#define N L*L
#define M 2*(L-1)*(L-1)
#define NC (L-4)*(L-4)
#define NB 2*(L-1)
#define NF 2*(L-3)// Se cumple N = 2*NB + NC + 2*NF

#define TERMALIZACION 10 // número de sweeps para la termalización
#define NTAU 10 // Número de configuraciones que se registran
#define TAU 1 // Se registran las configuraciones cada TAU sweeps 
#define RESET 5 //Periodo de sweeps en el que se resetean las posiciones 

#define DELTA0 0.1 // Lado del cubo inicial de donde se elige aleatoriamente el vector epsilon
#define MIN_Racept 0.40 // Porcentaje mínimo de la razon de aceptacion de configuraiones
#define MAX_Racept 0.60 // Porcentaje máximo de la razon de aceptacion de configuraiones

 
#define SI 1
#define NO 0

// Álgebra vectorial:

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
#define _restigual(u,v) \
                         u.a -= v.a; \
                         u.b -= v.b; \
                         u.c -= v.c;
#define _suma(u,v,w) \
			u.a = v.a + w.a; \
                        u.b = v.b + w.b; \
                        u.c = v.c + w.c;
#define _sumigual(u,v) \
                         u.a += v.a; \
                         u.b += v.b; \
                         u.c += v.c;
#define _inverso(u,v) \
		        u.a = -v.a; \
                        u.b = -v.b; \
                        u.c = -v.c;

// Parámetros números aleatorios:

#define MAX_RAND     2147483647L
#define I_FMAX_RAND     4.656612875e-10
#define NormRAN (1.0F/( (float) RAND_MAX+1.0F))
#define RAND() ( (float) rand() * NormRAN )
typedef unsigned long long randint;

//Variables for pseudo-random number generators :
randint ira[256];        //The Parisi-Rapuano wheel
unsigned char ip,ip1,ip2,ip3; //Variables for Parisis-Rapuano PRNG
randint zseed;            //The random variable for the congruential PRNG

#define TWOBRMINUS1   18446744073709551615ULL //2^64 - 1
#define TWOBR         18446744073709551616.   //2^64
#define TWOBRM1       9223372036854775808ULL  //2^63
#define FNORM (5.4210108624275218e-20)// max double such that RAND_MAX*NORMF<1

#define CGRANDOM  ( zseed=zseed*3202034522624059733LLU+1)
#define PRRANDOM  ( (ira[ip++]=ira[ip1++]+ira[ip2++]) ^ira[ip3++] )
//#define MYRANDOM   PRRANDOM                         // P-R generator
#define MYRANDOM   CGRANDOM                       // Congruential generator
#define HIRANDOM  ((randint) (PRRANDOM + CGRANDOM)) //Addition of both
#define FRANDOM (FNORM*MYRANDOM)   //A floating point, 0<=pseudo-random<1
#define FHIRANDOM (FNORM*HIRANDOM) //A floating point, 0<=pseudo-random<1

// Variables globales

int v1[N][6];
int dirv1_ini[N],dirv1_fin[N],zv1[N];

int v2[N][6];
int dirv2_ini[N],dirv2_fin[N],zv2[N];

int v1ord[N][6];

int p1ord[N][6];
int zp1[N];
int p2ord[N][6];

// Clasificaciónn puntos según vecinos:
int centro[NC];
int borde[2][NB];
int frontera[2][NF];

typedef struct{double a,b,c;} vector;

vector x[N],n[M];

typedef struct{int s1,s2;}vector2D; 

vector2D basev1[6];

typedef struct{vector2D min,max;}rect;

rect rombo[6];

int main(void)
{
  int i,cont;
  int sweep;
  double delta,ratio_acept;
  char nameout[255];
  FILE *output;

  regparam();// Guardamos los parámetros de la ejecución en archivo.
  

  /* 1ª PARTE: PARÁMETROS TOPOLÓGICOS */

  indice_vecnos_prox();  /*Cálcula v[N][6], dirv_ini[i], zv1[i]*/
  indice_vecnos_seg(); /*Cálcula v[N][6], dirv2_ini[i], z2[i]*/
  indice_vecnos_ord();
  indice_plqta_ord();
  // files_topology();  /*Guarda en archivos lo anteriormente calculado */

  // 2º PARTE: Configuración inicial

  tipos_nodos();
  vectores_iniciales();  /*Cálcula x[N] y n[N] para red triangular plana */
  sprintf(nameout,"xpos%d_inicial.dat",L,cont);
  output=fopen(nameout,"w");
  
  for(i=0; i<N; i++)
    fprintf(output,"%lf %lf %lf \n",x[i].a,x[i].b,x[i].c);
  fclose(output);

  // 3ª PARTE: Metropolis

  ini_random(SEMILLA);

  delta=DELTA0;

  for(sweep=0; sweep<TERMALIZACION; sweep++)// Termalizacion
    {
      ratio_acept=sweep_spiral(delta,KAPPA);

      if(sweep%RESET==0)
	{
	  for(i=1; i<N; i++)
	    {
	      _restigual(x[i],x[0]);
	    }
	  x[0].a=0.0F;
	  x[0].b=0.0F;
	  x[0].c=0.0F;
	}

      if(ratio_acept<MIN_Racept)
	delta/=2.0F;	
      if(ratio_acept>MAX_Racept)
	delta*=2.0F;

    }

  for(cont=0; cont<NTAU; cont++)// Registro configuraciones 
    {
      for(sweep=0; sweep<TAU; sweep++)
	{
	  ratio_acept=sweep_spiral(delta,KAPPA);

	  if(sweep%RESET==0)
	    {
	      for(i=1; i<N; i++)
		{
		  _restigual(x[i],x[0]);
		}
	      x[0].a=0.0F;
	      x[0].b=0.0F;
	      x[0].c=0.0F;
	    } 

	  if(ratio_acept<MIN_Racept)
	    delta/=2.0F;	
	  if(ratio_acept>MAX_Racept)
	    delta*=2.0F;
	}

      sprintf(nameout,"xpos_L%d_K%.1f-%d.dat",L,KAPPA,cont);
      output=fopen(nameout,"w");
      
      for(i=0; i<N; i++)
	fprintf(output,"%lf %lf %lf \n",x[i].a,x[i].b,x[i].c);
      fclose(output);

    }
    
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

/* INDICE_VECNOS_PROX()
*  Cálcula para cada nodo i de los N:
*   * zv1[i] : Número de cordinación para primeros vecinos
*   * dirv1_ini[i] : Primera dirección en la que el nodo i tiene vecino (sentido de giro positivo)
*   * dirv1_fin[i] : Última dirección en la que el nodo i tiene vecino (sentido de giro positivo)
*   * v1[i][j] : Índice del primer vecino del nodo i en la dirección j (ver figura)
*
* 		    v1[i][2]-----v1[i][1]
* 		          / \   / \
* 		         /   \ /   \
* 		 v1[i][3]-----i-----v0[i][0]
* 		         \   / \   /
* 		          \ /   \ /
* 		    v1[i][4]-----v[i][5]
*
*
*/
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
      rombo[dir].min.s1=sigma_min(basev1[dir].s1);
      rombo[dir].max.s1=sigma_max(basev1[dir].s1);
      rombo[dir].min.s2=sigma_min(basev1[dir].s2);
      rombo[dir].max.s2=sigma_max(basev1[dir].s2);
    }

  // Cálculo de v[N][6],zv[N]

  for(s2=0; s2<L; s2++)
    {
      for(s1=0; s1<L; s1++)
	{
	  i=s1+L*s2;
	  vec_new=NO;

	  for(dir=0; dir<6; dir++)
	    {
	      vec_old=vec_new;
	      if(s1>=rombo[dir].min.s1 && s1<rombo[dir].max.s1 && s2>=rombo[dir].min.s2 && s2<rombo[dir].max.s2)
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
	    }
	}
    }
}
/* INDICE_VECNOS_SEG()
*  Cálcula para cada nodo i de los N:
*   * zv2[i] : Número de cordinación para primeros vecinos
*   * dirv2_ini[i] : Primera dirección en la que el nodo i tiene segundo vecino (sentido de giro positivo)
*   * dirv2_fin[i] : Última dirección en la que el nodo i tiene segundo vecino (sentido de giro positivo)
*   * v2[i][j] : Índice del segundo vecino del nodo i en la dirección j (ver figura)
*
*			      v2[i][1]
*			     / \
*		 v2[i][2]---X---X---v2[i][0]
*			 \ / \ / \ /
*			  X---X---/
*			 / \ / \ / \
*		 v2[i][3]---X---X---v2[i][5]
* 			     \ /
* 			      v2[i][4]
*
*/
void indice_vecnos_seg(void)
{
  int s1,s2;
  int s1min,s1max,s2min,s2max;
  int dir;
  int s1_dir,s2_dir,s1_dir_old,s2_dir_old;
  int vecino_dir[N];
  int vecino_dir_anterior;
  int i;

  for(i=0; i<N; i++)
    vecino_dir[N]=NO;

  for(s1_dir=1,s2_dir=1,dir=0; dir<6; dir++)
    {
      switch(s1_dir)
	{
	case -2:
	  s1min=2;
	  s1max=L;
	  break;
	case -1:
	  s1min=1;
	  s1max=L;
	  break;
	case 1:
	  s1min=0;
	  s1max=L-1;
	  break;
	case 2:
	  s1min=0;
	  s1max=L-2;
	  break;
	}

      switch(s2_dir)
	{
	case -2:
	  s2min=2;
	  s2max=L;
	  break;
	case -1:
	  s2min=1;
	  s2max=L;
	  break;
	case 1:
	  s2min=0;
	  s2max=L-1;
	  break;
	case 2:
	  s2min=0;
	  s2max=L-2;
	  break;
	}

      for(s2=0; s2<L; s2++)
	{
	for(s1=0; s1<L; s1++)
	  {
	    i=s1+L*s2;
	    vecino_dir_anterior=vecino_dir[i];

	    if(s1>=s1min && s1<s1max && s2>=s2min && s2<s2max)
	      {
		v2[i][dir] = i + s1_dir + L*s2_dir;
		zv2[i]+=1;
		vecino_dir[i]=SI;
		if(vecino_dir_anterior==NO)
		  dirv2_ini[i]=dir;
	      }
	    else
	      {
		v2[i][dir]=-1;
		vecino_dir[i]=NO;
		if(vecino_dir_anterior==SI)
		  {
		    if(dir==0)
		      dirv2_fin[i]=5;
		    else
		      dirv2_fin[i]=dir-1;
		  }
	      }
	  }
	}

      s1_dir_old=s1_dir;
      s2_dir_old=s2_dir;

      s1_dir=-s2_dir_old;
      s2_dir=s1_dir_old + s2_dir_old;
    }

  dirv2_ini[L+1]=5;
  dirv2_fin[L+1]=3; /*Punto excepcional: Vecinos segundos no conexos*/

  dirv2_ini[N-(L+2)]=2;
  dirv2_fin[N-(L+2)]=0; /*Punto excepcional: vecinos segundos no conexos*/

  }

void indice_vecnos_ord(void)
{
  int i,dir,ord;

  for(i=0; i<N; i++)
    {
      for(ord=0; ord<6; ord++)
	v1ord[i][ord]=-1;
    }


  for(i=0; i<N; i++)
    {
      if(zv1[i]==6)
	{
	  for(dir=dirv2_ini[i],ord=0; ord<6; dir++,ord++)
	    {
	      if(dir>5)
		dir=0;

	      v1ord[i][ord]=v1[i][dir];
	    }
	}
      else
	{
	  for(dir=dirv1_ini[i],ord=0; ord<zv1[i]; dir++,ord++)
	    {
	      if(dir>5)
		dir=0;

	      v1ord[i][ord]=v1[i][dir];
	    }
	}
    }
}

void indice_plqta_ord(void)
{
  int i,dir,dir_min,ord;

  for(i=0; i<N; i++)
    {
      for(ord=0; ord<6; ord++)
	{
	  p1ord[i][ord]=-1;
	  p2ord[i][ord]=-1;
	}
    }


  for(i=0; i<N; i++)
    {
      if(zv1[i]==6)
	{
	  zp1[i]=6;
	  dir_min=dirv2_ini[i];
	}
      else
	{
	  zp1[i]=zv1[i]-1;
	  dir_min=dirv1_ini[i];
	}

      for(dir=dir_min, ord=0; ord<zp1[i]; dir++, ord++)
	{
	  if(dir>5)
	    dir=0;

	  p1ord[i][ord]=index_plqta_prox(i,dir);
	}


      for(dir=dirv2_ini[i],ord=0; ord<zv2[i]; dir++,ord++)
	{
	  if(dir>5)
	    dir=0;

	  p2ord[i][ord]=index_plqta_seg(i,dir);
	}
    }
  /*CORECCIONES:*/

  /*Corremos los indices en una unidad de p2ord correspondientes a L  N-(L+1) */
  p2ord[L][2]=p2ord[L][1];
  p2ord[L][1]=p2ord[L][0];
  p2ord[L][0]=-1;

  p2ord[N-(L+1)][2]=p2ord[N-(L+1)][1];
  p2ord[N-(L+1)][1]=p2ord[N-(L+1)][0];
  p2ord[N-(L+1)][0]=-1;

  /*Asignamos el valor que falta en los no conexos*/

  p2ord[L+1][4]=index_plqta_seg(L+1,3);
  p2ord[N-(L+2)][4]=index_plqta_seg(N-(L+2),0);
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

/*INDEX_PLQTA_SEG(int i int dir):
*  Devuelve el índice de la plaqueta vecina segunda del nodo i en la dirección dir
*   * Las direcciones de las plaquetas vecinas segundas se muestran en la figura:
*
* 		           x
* 		          / \
* 		         / 1 \
* 	          x-----x-----x-----x
* 		   \ 2 / \   / \ 0 /
* 		    \ /   \ /   \ /
* 	             x-----i-----X
*		    / \   / \   / \
*		   / 3 \ /   \ / 6 \
*		  x-----x-----x-----x
*		         \ 5 /
*		       	  \ /
*			   x
*
*/
int index_plqta_seg(int i, int dir)
{

  int p;

  switch(dir)
    {
    case 0:
      p=index_plqta_dir1(v1[i][0]);
      break;
    case 1:
      p=index_plqta_dir0(v1[i][2]);
      break;
    case 2:
      p=index_plqta_dir1(v1[i][3]);
      break;
    case 3:
      p=index_plqta_dir0(v2[i][3]);
      break;
    case 4:
      p=index_plqta_dir1(v2[i][4]);
      break;
    case 5:
      p=index_plqta_dir0(v1[i][5]);
      break;
    }

  return p;
}
/* TIPOS_NOD0S()
*  Cálcula:
*   * centro[0 =< i < NC] (Representados por @ en la figura)
*   * frontera[0=< j < 2][0 =< i < NF] (Representados por o  en la figura)    
*   * borde[0=< j < 2][0 =< i < NB] (Representados por *  en la figura)
*
* 		           ****************************(N-1)
* 		          * \   / \   /	\   / \   / \   /
* 		         *   \ /   \ /   \ /   \ /   \ /
* 		        *-----ooooooooooooooooooo-----*
* 		       * \   o \   / \	 / \   / \   *
* 		      *   \ o   \ /   \ /   \ /   \ *
* 		     *-----o-----@-----@-----o-----*
*		    * \	  o \   / \   / \   o \   *
*		   *   \ o   \ /   \ /   \ o   \ *
*	          *-----o-----@-----@-----o-----*
*	         *     / \   / \   / \	 o \   *
*	        *   \ /   \ /   \ /   \ o   \ *
*	      (L)--(L+1)oooooooooooooooo---  *
*	      / \   / \	  / \   / \   / \   *
*	     /   \ /   \ /   \ /   \ /	 \ *
*	   (0)*****************************
*
*
*/
void tipos_nodos(void)
{
  int s1,s2;
  int i,cont;

  //Indices de los nodos centrales, aquellos que tienen todos sus vecinos próximos y segundos

  cont=0;
  for(s2=0; s2<L; s2++)
	{
	  for(s1=0; s1<L; s1++)
	    {
	      i=s1+L*s2;
	      if( zv1[i]==6 && zv2[i]==6)
		{
		  centro[cont]=i;
		  cont++;
		}
	    }
	}
  //Indices de los nodos del contorno

  cont=0;
  s2=0;
  for(s1=0; s1<L; s1++)
	{
	  i=s1+L*s2;
	  borde[0][cont]=i;
	  cont++;
	}
  s1=L-1;
  for(s2=1; s2<L-1; s2++)
	{
	  i=s1+L*s2;
	  borde[0][cont]=i;
	  cont++;
	}

  cont=0;
  s2=L-1;
  for(s1=L-1; s1>0; s1--)
	{
	  i=s1+L*s2;
	  borde[1][cont]=i;
	  cont++;
	}
  s1=0;
  for(s2=L-1; s2>0; s2--)
	{
	  i=s1+L*s2;
	  borde[1][cont]=i;
	  cont++;
	}

  //Indices de los nodos de la frontera

  cont=0;
  s2=1;
  for(s1=1; s1<L-1; s1++)
	{
	  i=s1+L*s2;
	  frontera[0][cont]=i;
	  cont++;
	}
  s1=L-2;
  for(s2=2; s2<L-2; s2++)
	{
	  i=s1+L*s2;
	  frontera[0][cont]=i;
	  cont++;
	}
  cont=0;
  s2=L-2;
  for(s1=L-2; s1>1; s1--)
	{
	  i=s1+L*s2;
	  frontera[1][cont]=i;
	  cont++;
	}
  s1=1;
  for(s2=L-2; s2>1; s2--)
	{
	  i=s1+L*s2;
	  frontera[1][cont]=i;
	  cont++;
	}
}
/*VECTORES INICIALES: Cálculo vectores iniciales correspondientes red triangular regular plana*/
void vectores_iniciales(void)
{
  vector sigma1,sigma2;
  int i,p,s1,s2;
  double l;
  
  l=1.0F; /* Distancia entre vecinos red triangular */
	
  sigma1.a=l;
  sigma1.b=0.0F;
  sigma1.c=0.0F;
  
  sigma2.a=0.5F*l;
  sigma2.b=0.5F*sqrt(3)*l;
  sigma2.c=0.0F;

  /* Vectores de posición iniciales = red triangular regular*/
  

  for(i=0; i<N; i++)
    {
      s1=i % L;
      s2=i/L;

      x[i].a= sigma1.a * (double) s1+ sigma2.a * (double) s2;
      x[i].b= sigma2.b * (double) s2;
      x[i].c= 0.0F;
    }

  /* Vectores normales iniciales = red triangular regular */

  for(p=0; p<M; p++)
    {
      n[p].a=0.0F;
      n[p].b=0.0F;
      n[p].c=1.0F;
    }
}

/*SWEEP_SPIRAL: Aplica algoritmo de metropolis*/
double sweep_spiral(double delta, double kappa)
{
  int sum_acept;
  double ratio_acept;
  int index,frag;
  int i;
  int j;
  int p;
  int k;
  int l;
  double deltaH; // Incremento de energía total
  double deltaHt; // Incremento de energía elástica
  double deltaHb; // Incremento de energía de curvatura

  double prod_ep2_sumr,prod_rad,prod_circ_old,prod_circ_new;

  double norma2,invnorma,norma2ep2;

  vector ep2ilon,xtemp,rtemp[6],ntemp[6],r[6];
  vector sumr,difn;

  sum_acept=0;
  /* 1ª NODOS BORDE:
   * Se aplica Metropolis de forma circular a los nodos del borde exterior
   */
for(frag=0;frag<2;frag++)
    {
      for(index=0; index<NB-1; index++)
	{
	  i=borde[frag][index]; //La variable i es el índice del nodo a modificar

	  /* Inicializamos variables */
	  deltaH=0.0F;
	  deltaHt=0.0F;
	  deltaHb=0.0F;
	  sumr.a=0.0F;
	  sumr.b=0.0F;
	  sumr.c=0.0F;	  
	  /* Cálculo de vectores: */      
	  ep2ilon.a=delta*(2*FRANDOM-1);// Elección aleatoria del vector ep2ilon:
	  ep2ilon.b=delta*(2*FRANDOM-1);//  Elegido aleatoriamente y uniformente en   
	  ep2ilon.c=delta*(2*FRANDOM-1);//  un cubo de lado 2*delta
	  norma2ep2=_norma2(ep2ilon);// Norma al cuadrado del vector ep2ilon (necesario para el 
	  // incremento de energía elástica).
	  _suma(xtemp,x[i],ep2ilon); //Nuevo vector de posición para el punto i
	  for(j=0; j<zv1[i]; j++) // Vectores relativos con los vecinos próximos
	    {
	      _resta(rtemp[j], x[v1ord[i][j]], xtemp);
	      _resta(r[j], x[v1ord[i][j]], x[i]); 
	      _sumigual(sumr,r[j]);// Vector suma de los r[i] (necesario para el incremento de energía elástica)
	    }
	  for(k=0; k<zv1[i]-1; k++) // Nuevos vectores normales de las palquetas adyacentes
	    {
	      _prodvec(ntemp[k],rtemp[k],rtemp[k+1]);
	      norma2=_norma2(ntemp[k]);
	      invnorma=1.0F/sqrt(norma2);
	      _escala(ntemp[k],invnorma);
	    }
	  
	  /* Cálculo incremento energía: */
	  
	  //Incremento de la energía elástica
	  
	  prod_ep2_sumr=_prodesc(ep2ilon,sumr);
	  deltaHt=(((double) zv1[i])*norma2ep2)-(2.0F*prod_ep2_sumr);
	  
	  //Incremento de la energía de curvatura:
	  for(k=0; k<zv1[i]-2; k++)// Productos circulares de las normales
	    {
	      prod_circ_old=_prodesc(n[p1ord[i][k]],n[p1ord[i][k+1]]);
	      prod_circ_new=_prodesc(ntemp[k],ntemp[k+1]);
	      deltaHb+=(prod_circ_new-prod_circ_old);
	    }

	  for(l=0; l<zv2[i]; l++)// Productos radiales de las normales
	    {
	      _resta(difn,ntemp[l],n[p1ord[i][l]]);
	      prod_rad=_prodesc(difn,n[p2ord[i][l]]);
	      deltaHb += prod_rad;
	    }
	 
	  deltaHb*=KAPPA;
	  deltaH=deltaHt-deltaHb;//Incremento total de Energía
	  
	  /* Se decide si se acepta la nueva posición del nodo i*/
	  if(deltaH<=0 || exp(-deltaH)> FRANDOM)
	    {
	      x[i]=xtemp;// Actualizamos el vector de posición x[i]
	      for(k=0; k<zv1[i]-1; k++)
		{
		  n[p1ord[i][k]]=ntemp[k];// Actualizamos las normales de las plaquetas vecinas al nodo i
 		}
	      sum_acept+=1; /*aumenta una unidad la cantidad de nuevas posiciones aceptadas*/
	    }
	}
      /* Uĺtimos puntos del borde: La dirección inicial de v1 no coincide con la de v2
       * 
       */
      i=borde[frag][NB-1]; //La variable i es el índice del nodo a modificar
    
	  /* Inicializamos variables */
	  deltaH=0.0F;
	  deltaHt=0.0F;
	  deltaHb=0.0F;
	  sumr.a=0.0F;
	  sumr.b=0.0F;
	  sumr.c=0.0F;
	  
	  /* Cálculo de vectores: */      
	  ep2ilon.a=delta*(2*FRANDOM-1);// Elección aleatoria del vector ep2ilon:
	  ep2ilon.b=delta*(2*FRANDOM-1);//  Elegido aleatoriamente y uniformente en   
	  ep2ilon.c=delta*(2*FRANDOM-1);//  un cubo de lado 2*delta
	  norma2ep2=_norma2(ep2ilon);// Norma al cuadrado del vector ep2ilon (necesario para el 
	  // incremento de energía elástica).
	  _suma(xtemp,x[i],ep2ilon); //Nuevo vector de posición para el punto i
	  for(j=0; j<zv1[i]; j++) // Vectores relativos con los vecinos próximos
	    {
	      _resta(rtemp[j], x[v1ord[i][j]], xtemp);
	      _resta(r[j], x[v1ord[i][j]], x[i]); 
	      _sumigual(sumr,r[j]);// Vector suma de los r[i] (necesario para el incremento de energía elástica)
	    }
	  for(k=0; k<zv1[i]-1; k++) // Nuevos vectores normales de las palquetas adyacentes
	    {
	      _prodvec(ntemp[k],rtemp[k],rtemp[k+1]);
	      norma2=_norma2(ntemp[k]);
	      invnorma=1.0F/sqrt(norma2);
	      _escala(ntemp[k],invnorma);
	    }
	  
	  /* Cálculo incremento energía: */
	  
	  //Incremento de la energía elástica
	  
	  prod_ep2_sumr=_prodesc(ep2ilon,sumr);
	  deltaHt=(((double) zv1[i])*norma2ep2)-(2.0F*prod_ep2_sumr);
	  
	  //Incremento de la energía de curvatura:
	  for(k=0; k<zv1[i]-2; k++)// Productos circulares de las normales
	    {
	      prod_circ_old=_prodesc(n[p1ord[i][k]],n[p1ord[i][k+1]]);
	      prod_circ_new=_prodesc(ntemp[k],ntemp[k+1]);
	      deltaHb+=(prod_circ_new-prod_circ_old);
	    }

	  for(l=1; l<zv2[i]+1; l++)// Productos radiales de las normales(no empieza en 0!)
	    {
	      _resta(difn,ntemp[l],n[p1ord[i][l]]);
	      prod_rad=_prodesc(difn,n[p2ord[i][l]]);
	      deltaHb += prod_rad;
	    }
	 
	  deltaHb*=KAPPA;
	  deltaH=deltaHt-deltaHb;//Incremento total de Energía
	  
	  /* Se decide si se acepta la nueva posición del nodo i*/
	  if(deltaH<=0 || exp(-deltaH)> FRANDOM)
	    {
	      x[i]=xtemp;// Actualizamos el vector de posición x[i]
	      for(k=0; k<zv1[i]-1; k++)
		{
		  n[p1ord[i][k]]=ntemp[k];// Actualizamos las normales de las plaquetas vecinas al nodo i
 		}
	      sum_acept+=1; /*aumenta una unidad la cantidad de nuevas posiciones aceptadas*/
	    }
	
    }
     
  /* 2ª NODOS FRONTERA:
   * Se aplica Metropolis de forma circular a los nodos del borde interior
   */
  for(frag=0;frag<2;frag++)
    {
      /* Esquinas de la frontera: Las 2as plaquetas no son conexas
       * 
       */
      i=frontera[frag][0]; //La variable i es el índice del nodo a modificar

      /* Inicializamos variables */
      deltaH=0.0F;
      deltaHt=0.0F;
      deltaHb=0.0F;
      sumr.a=0.0F;
      sumr.b=0.0F;
      sumr.c=0.0F;

      /* Cálculo de vectores: */      
      ep2ilon.a=delta*(2*FRANDOM-1);// Elección aleatoria del vector ep2ilon:
      ep2ilon.b=delta*(2*FRANDOM-1);//  Elegido aleatoriamente y uniformente en   
      ep2ilon.c=delta*(2*FRANDOM-1);//  un cubo de lado 2*delta
      norma2ep2=_norma2(ep2ilon);// Norma al cuadrado del vector ep2ilon (necesario para el 
	                         // incremento de energía elástica).
      _suma(xtemp,x[i],ep2ilon); //Nuevo vector de posición para el punto i
      for(j=0; j<6; j++) // Vectores relativos con los vecinos próximos
	{
	  _resta(rtemp[j], x[v1ord[i][j]], xtemp);
	  _resta(r[j], x[v1ord[i][j]], x[i]); 
	  _sumigual(sumr,r[j]);// Vector suma de los r[i] (necesario para el incremento de energía elástica)
	}
      for(k=0; k<5; k++) // Nuevos vectores normales de las palquetas adyacentes
	{
	  _prodvec(ntemp[k],rtemp[k],rtemp[k+1]);
	  norma2=_norma2(ntemp[k]);
	  invnorma=1.0F/sqrt(norma2);
	  _escala(ntemp[k],invnorma);
	}
      _prodvec(ntemp[5],rtemp[5],rtemp[0]);
      norma2=_norma2(ntemp[5]);
      invnorma=1.0F/sqrt(norma2);
      _escala(ntemp[5],invnorma);

      /* Cálculo incremento energía: */

      //Incremento de la energía elástica
     
      prod_ep2_sumr=_prodesc(ep2ilon,sumr);
      deltaHt=(6.0F*norma2ep2)-(2.0F*prod_ep2_sumr);

      //Incremento de la energía de curvatura:
      for(k=0; k<5; k++)// Productos circulares de las normales
	{
	  prod_circ_old=_prodesc(n[p1ord[i][k]],n[p1ord[i][k+1]]);
	  prod_circ_new=_prodesc(ntemp[k],ntemp[k+1]);
	  deltaHb+=(prod_circ_new-prod_circ_old);
	}
      prod_circ_old=_prodesc(n[p1ord[i][5]],n[p1ord[i][0]]);
      prod_circ_new=_prodesc(ntemp[5],ntemp[0]);
      deltaHb+=(prod_circ_new-prod_circ_old);
      for(l=0; l<3; l++)// Productos radiales de las normales
	{
	  _resta(difn,ntemp[l],n[p1ord[i][l]]);
	  prod_rad=_prodesc(difn,n[p2ord[i][l]]);
	  deltaHb += prod_rad;
	}
      _resta(difn,ntemp[4],n[p1ord[i][4]]);
      prod_rad=_prodesc(difn,n[p2ord[i][4]]);
      deltaHb += prod_rad;

      deltaHb*=KAPPA;
      deltaH=deltaHt-deltaHb;//Incremento total de Energía

      /* Se decide si se acepta la nueva posición del nodo i*/
      if(deltaH<=0 || exp(-deltaH)> FRANDOM)
	{
	  x[i]=xtemp;// Actualizamos el vector de posición x[i]
	  for(k=0; k<6; k++)
	    {
	      n[p1ord[i][k]]=ntemp[k];// Actualizamos las normales de las plaquetas vecinas al nodo i
	    }
	  sum_acept+=1; /*aumenta una unidad la cantidad de nuevas posiciones aceptadas*/
	}

      for(index=1; index<NF; index++)
	{
	  i=frontera[frag][index]; //La variable i es el índice del nodo a modificar

	  /* Inicializamos variables */
	  deltaH=0.0F;
	  deltaHt=0.0F;
	  deltaHb=0.0F;
	  sumr.a=0.0F;
	  sumr.b=0.0F;
	  sumr.c=0.0F;
	  
	  /* Cálculo de vectores: */      
	  ep2ilon.a=delta*(2*FRANDOM-1);// Elección aleatoria del vector ep2ilon:
	  ep2ilon.b=delta*(2*FRANDOM-1);//  Elegido aleatoriamente y uniformente en   
	  ep2ilon.c=delta*(2*FRANDOM-1);//  un cubo de lado 2*delta
	  norma2ep2=_norma2(ep2ilon);// Norma al cuadrado del vector ep2ilon (necesario para el 
	  // incremento de energía elástica).
	  _suma(xtemp,x[i],ep2ilon); //Nuevo vector de posición para el punto i
	  for(j=0; j<6; j++) // Vectores relativos con los vecinos próximos
	    {
	      _resta(rtemp[j], x[v1ord[i][j]], xtemp);
	      _resta(r[j], x[v1ord[i][j]], x[i]); 
	      _sumigual(sumr,r[j]);// Vector suma de los r[i] (necesario para el incremento de energía elástica)
	    }
	  for(k=0; k<5; k++) // Nuevos vectores normales de las palquetas adyacentes
	    {
	      _prodvec(ntemp[k],rtemp[k],rtemp[k+1]);
	      norma2=_norma2(ntemp[k]);
	      invnorma=1.0F/sqrt(norma2);
	      _escala(ntemp[k],invnorma);
	    }
	  _prodvec(ntemp[5],rtemp[5],rtemp[0]);
	  norma2=_norma2(ntemp[5]);
	  invnorma=1.0F/sqrt(norma2);
	  _escala(ntemp[5],invnorma);
	  
	  /* Cálculo incremento energía: */
	  
	  //Incremento de la energía elástica
	  
	  prod_ep2_sumr=_prodesc(ep2ilon,sumr);
	  deltaHt=(6.0F*norma2ep2)-(2.0F*prod_ep2_sumr);
	  
	  //Incremento de la energía de curvatura:
	  for(k=0; k<5; k++)// Productos circulares de las normales
	    {
	      prod_circ_old=_prodesc(n[p1ord[i][k]],n[p1ord[i][k+1]]);
	      prod_circ_new=_prodesc(ntemp[k],ntemp[k+1]);
	      deltaHb+=(prod_circ_new-prod_circ_old);
	    }
	  prod_circ_old=_prodesc(n[p1ord[i][5]],n[p1ord[i][0]]);
	  prod_circ_new=_prodesc(ntemp[5],ntemp[0]);
	  deltaHb+=(prod_circ_new-prod_circ_old);
	  for(l=0; l<zv2[i]; l++)// Productos radiales de las normales
	    {
	      _resta(difn,ntemp[l],n[p1ord[i][l]]);
	      prod_rad=_prodesc(difn,n[p2ord[i][l]]);
	      deltaHb += prod_rad;
	    }
	 
	  deltaHb*=KAPPA;
	  deltaH=deltaHt-deltaHb;//Incremento total de Energía
	  
	  /* Se decide si se acepta la nueva posición del nodo i*/
	  if(deltaH<=0 || exp(-deltaH)> FRANDOM)
	    {
	      x[i]=xtemp;// Actualizamos el vector de posición x[i]
	      for(k=0; k<6; k++)
		{
		  n[p1ord[i][k]]=ntemp[k];// Actualizamos las normales de las plaquetas vecinas al nodo i
		}
	      sum_acept+=1; /*aumenta una unidad la cantidad de nuevas posiciones aceptadas*/
	    }
	}
      
    }
  
  /* 3ª NODOS CENTRALES: 
   *  Se aplica el algoritmo de Metropolis para cada nodo central forma secuencial
   */
  for(index=0; index<NC; index++)
    {
      i=centro[index]; //La variable i es el índice del nodo a modificar

      /* Inicializamos variables */
      deltaH=0.0F;
      deltaHt=0.0F;
      deltaHb=0.0F;
      sumr.a=0.0F;
      sumr.b=0.0F;
      sumr.c=0.0F;

      /* Cálculo de vectores: */      
      ep2ilon.a=delta*(2*FRANDOM-1);// Elección aleatoria del vector ep2ilon:
      ep2ilon.b=delta*(2*FRANDOM-1);//  Elegido aleatoriamente y uniformente en   
      ep2ilon.c=delta*(2*FRANDOM-1);//  un cubo de lado 2*delta
      norma2ep2=_norma2(ep2ilon);// Norma al cuadrado del vector ep2ilon (necesario para el 
	                         // incremento de energía elástica).
      _suma(xtemp,x[i],ep2ilon); //Nuevo vector de posición para el punto i
      for(j=0; j<6; j++) // Vectores relativos con los vecinos próximos
	{
	  _resta(rtemp[j], x[v1ord[i][j]], xtemp);
	  _resta(r[j], x[v1ord[i][j]], x[i]); 
	  _sumigual(sumr,r[j]);// Vector suma de los r[i] (necesario para el incremento de energía elástica)
	}
      for(k=0; k<5; k++) // Nuevos vectores normales de las palquetas adyacentes
	{
	  _prodvec(ntemp[k],rtemp[k],rtemp[k+1]);
	  norma2=_norma2(ntemp[k]);
	  invnorma=1.0F/sqrt(norma2);
	  _escala(ntemp[k],invnorma);
	}
      _prodvec(ntemp[5],rtemp[5],rtemp[0]);
      norma2=_norma2(ntemp[5]);
      invnorma=1.0F/sqrt(norma2);
      _escala(ntemp[5],invnorma);

      /* Cálculo incremento energía: */

      //Incremento de la energía elástica
     
      prod_ep2_sumr=_prodesc(ep2ilon,sumr);
      deltaHt=(6.0F*norma2ep2)-(2.0F*prod_ep2_sumr);

      //Incremento de la energía de curvatura:
      for(k=0; k<5; k++)// Productos circulares de las normales
	{
	  prod_circ_old=_prodesc(n[p1ord[i][k]],n[p1ord[i][k+1]]);
	  prod_circ_new=_prodesc(ntemp[k],ntemp[k+1]);
	  deltaHb+=(prod_circ_new-prod_circ_old);
	}
      prod_circ_old=_prodesc(n[p1ord[i][5]],n[p1ord[i][0]]);
      prod_circ_new=_prodesc(ntemp[5],ntemp[0]);
      deltaHb+=(prod_circ_new-prod_circ_old);
      for(l=0; l<6; l++)// Productos radiales de las normales
	{
	  _resta(difn,ntemp[l],n[p1ord[i][l]]);
	  prod_rad=_prodesc(difn,n[p2ord[i][l]]);
	  deltaHb += prod_rad;
	}

      deltaHb*=KAPPA;
      deltaH=deltaHt-deltaHb;//Incremento total de Energía

      /* Se decide si se acepta la nueva posición del nodo i*/
      if(deltaH<=0 || exp(-deltaH)> FRANDOM)
	{
	  x[i]=xtemp;// Actualizamos el vector de posición x[i]
	  for(k=0; k<6; k++)
	    {
	      n[p1ord[i][k]]=ntemp[k];// Actualizamos las normales de las plaquetas vecinas al nodo i
	    }
	  sum_acept+=1; /*aumenta una unidad la cantidad de nuevas posiciones aceptadas*/
	}
    }

  ratio_acept=((double) sum_acept)/ ((double) N);

  return ratio_acept;
}

// Función números aleatorios
void ini_random(unsigned long long semilla)
{
  int i;

  zseed=semilla;
  for (i=0;i<11;i++)         /* Just in case initial randomseed were small */
      CGRANDOM;

  ip=128;
  ip1=ip-24;
  ip2=ip-55;
  ip3=ip-61;
  for (i=ip3; i<ip; i++)
      ira[i] = CGRANDOM;

  for (i=0;i<1111;i++)
    if (!PRRANDOM)
      printf("Found zero in the first P-R random numbers generated\n");
}
// regparam: Parametros de la simulacion

void regparam(void )
{
  FILE *output;
  int i;

  output=fopen("regparams.log","w"); 
  
  fprintf(output,"Tamaño de la red: \n");
  fprintf(output,"L=%d N=%d M=%d \n",L,N,M); 
  fprintf(output,"NC=%d NB=%d NF=%d \n",NC,NB,NF);
  fprintf(output,"\n");
  fprintf(output,"Sweep2: \n");
  fprintf(output,"Termalizacion=%d Ntau=%d tau=%d Reset=%d \n",TERMALIZACION,NTAU,TAU,RESET);
  fprintf(output,"Semilla=%d \n",SEMILLA);
  fprintf(output,"\n");
  fprintf(output,"Parametros Metropolis:\n");
  fprintf(output,"Kappa=%f\n",KAPPA);
  fprintf(output,"Delta inicial=%f mín.razon_acept=%f max.razon_acept=%f \n",DELTA0,MIN_Racept,MAX_Racept);
  fprintf(output,"\n");
  fclose(output);

}

/*FILES_TOPOLOGY: Almacena en archivos los parámetros topológicos*/
void files_topology(void )
{
  FILE *output;  

  int i;


  output=fopen("indicev1.dat","w"); /* dirv1_ini[i]   dirv1_fin[i]   zv1[i] */
  for(i=0; i<N; i++)
    {
      fprintf(output,"%4d %4d %4d \n",dirv1_ini[i], dirv1_fin[i], zv1[i]);
    }
  fclose(output);

  output=fopen("v1.dat","w");
  for(i=0; i<N; i++)
    {
      fprintf(output,"%3d %3d %3d %3d %3d %3d\n",v1[i][0],v1[i][1],v1[i][2],v1[i][3],v1[i][4],v1[i][5]);
    }
  fclose(output);
output=fopen("indicev2.dat","w"); /* dirv1_ini[i]   dirv1_fin[i]   zv1[i] */
  for(i=0; i<N; i++)
    {
      fprintf(output,"%4d %4d %4d \n",dirv2_ini[i], dirv2_fin[i], zv2[i]);
    }
  fclose(output);
output=fopen("index_comparacion.dat","w");
  for(i=0; i<N; i++)
    {
      fprintf(output,"%4d %4d %4d %4d\n",dirv1_ini[i], zv1[i],dirv2_ini[i], zv2[i]);
    }
  fclose(output);
  output=fopen("v2.dat","w");
  for(i=0; i<N; i++)
    {
      fprintf(output,"%3d %3d %3d %3d %3d %3d\n",v2[i][0],v2[i][1],v2[i][2],v2[i][3],v2[i][4],v2[i][5]);
    }
  fclose(output);

output=fopen("v_ord.dat","w");
  for(i=0; i<N; i++)
    {
      fprintf(output,"%3d %3d %3d %3d %3d %3d\n",v1ord[i][0],v1ord[i][1],v1ord[i][2],v1ord[i][3],v1ord[i][4],v1ord[i][5]);
    }
  fclose(output);

output=fopen("pordinal.dat","w");/* 0   1   2   3   4   5   6   7   8   9  10  11 */
  for(i=0; i<N; i++)
    {
      fprintf(output,"%3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d \n",p1ord[i][0],p1ord[i][1],p1ord[i][2],p1ord[i][3],p1ord[i][4],p1ord[i][5],p2ord[i][0],p2ord[i][1],p2ord[i][2],p2ord[i][3],p2ord[i][4],p2ord[i][5]);
    }
  fclose(output);
}
