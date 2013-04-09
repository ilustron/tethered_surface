#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double radio_giraton(void);
double energia_curvatura(void);

void index_vcnos_prox();

int index_plqta_dir0(int i);
int index_plqta_dir1(int i);
int index_plqta_prox(int i, int dir);


#define L 16
#define N L*L 
#define M 2*(L-1)*(L-1) 
#define K 1.1 // Kappa 
#define NF 24000 // nº de archivos

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

typedef struct{double a,b,c;} vector;

vector x[N];
int v[N][6],p[N][6];// indices puntos y plaquetas vecinos

int main(void )
{
  FILE *input;
  FILE *output;
  FILE *pipe = popen("gnuplot -persist","w");
  
  char namein[255];
  char nameout[255];
  
  int i,f,dir;
  int b,n; // nº de elementos del bloque

  int btermal,ntermal;

  // Observables:
  double energia[NF]; 
  double radio2g[NF];
  double caloresp[NF];

  double Energia[NF]; // Acumulados
  double Energia2[NF]; // Acumulados
  double Radio2g[NF];

  double media_radio2g;// media del radio del giraton al cuadrado
  double media_energia;// media del radio de la Energia de curvatura
  double media_energia2;
  double Cv;
  
  double refRadio2g;
  double refEnergia;
  double refEnergia2;

  // Cargamos los índices v[N][6]

  index_vcnos_prox();
  
  // Cargamos los índices p[N][6]

  for(i=0; i<N; i++)
    {
      for(dir=0; dir<6; dir++)
	{
	  p[i][dir]=index_plqta_prox(i,dir);
	}
    } 
  
  // Lectura de archivos de posiciones

  f=0;

  sprintf(namein,"xpos_L%d_K1.1-%d.dat",L,f);

  Radio2g[NF-1]=0.0F;
  Energia2[f]=0.0F;
  Energia2[NF-1]=0.0F;

  while((input=fopen(namein,"r"))!=NULL)
    {
      i=0;
      
      while(fscanf(input,"%lf %lf %lf",&x[i].a,&x[i].b,&x[i].c)!=EOF)
	i++;
   
      if(i!=N)
	printf("Error: El fichero %s no contiene %d líneas\n",namein,N);      
      fclose(input);

      radio2g[f]=radio_giraton();

      Radio2g[NF-1]+=radio2g[f];
      Radio2g[f]=Radio2g[NF-1];

      energia[f]=energia_curvatura();

      Energia[NF-1]+=energia[f];
      Energia[f]=Energia[NF-1];

      Energia2[NF-1]+=(energia[f]*energia[f]);
      Energia2[f]=Energia2[NF-1];

      f++;
      sprintf(namein,"xpos_L%d_K1.1-%d.dat",L,f);
    }
  
  if(f!=NF)
    {
      printf("El número de archivos de posiciones leidos no es %d\n",NF);
      return 0;
    }

  // Escala logaritmica

  output=fopen("termalizacion.dat","w"); 
  for(f=1; f<=NF; f++)
    {      
       if((NF%f)==0)
	 {
	      b=NF/f;
	      n=b-1;
	      media_radio2g=Radio2g[n]/(double)b;
	      media_energia=Energia[n]/(double)b;
	      Cv=Energia2[n]/(double)b-Energia[n]*Energia[n]/(double)b/(double)b;
	      Cv*=(K*K)/((double)N);
	      fprintf(output,"%d %lf %lf %lf\n",b,media_energia,Cv,media_radio2g);
	 }
    }
  fclose(output);
  
  fprintf(pipe, "set logscale x\n");
  fprintf(pipe, "set multiplot layout 3,1 \n");
  fprintf(pipe, "plot \"./termalizacion.dat\" u 1:2 title \"energia\" w lp\n");
  fprintf(pipe, "plot \"./termalizacion.dat\" u 1:3 title \"Cv\" w lp\n");
  fprintf(pipe, "plot \"./termalizacion.dat\" u 1:4 title \"Radio2 g\" w lp\n");
  fprintf(pipe, "unset multiplot\n"); 
  fflush(pipe);


  // Lee el sweep correspondiente a la termalizacion

  printf("Escribe el sweep correspondiente a la termalización:\n");
  ntermal=scanf("%d",&btermal);

  //redefinimos los observables

  refRadio2g=Radio2g[btermal-1];
  refEnergia=Energia[btermal-1];
  refEnergia2=Energia2[btermal-1l];

  for(f=0; f<(NF-btermal); f++)
    {
      Radio2g[f]=Radio2g[f+btermal]-refRadio2g;
      Energia[f]=Energia[f+btermal]-refEnergia;
      Energia2[f]=Energia2[f+btermal]-refEnergia2;
    }
  

  output=fopen("termalizacion2.dat","w"); 
  for(f=1; f<=(NF-btermal); f++)
    {
      if(((NF-btermal)%f)==0)
	 {
	      b=(NF-btermal)/f;
	      n=b-1;
	      media_radio2g=Radio2g[n]/(double)b;
	      media_energia=Energia[n]/(double)b;
	      Cv=Energia2[n]/(double)b-Energia[n]*Energia[n]/(double)b/(double)b;
	      Cv*=(K*K)/((double)N);
	      fprintf(output,"%d %lf %lf %lf\n",b,media_energia,Cv,media_radio2g);
	 }
    }
  fclose(output);

  pipe = popen("gnuplot -persist","w");  
  fprintf(pipe, "set logscale x\n");
  fprintf(pipe, "set multiplot layout 3,1 \n");
  fprintf(pipe, "plot \"./termalizacion2.dat\" u 1:2 title \"energia\" w lp\n");
  fprintf(pipe, "plot \"./termalizacion2.dat\" u 1:3 title \"Cv\" w lp\n");
  fprintf(pipe, "plot \"./termalizacion2.dat\" u 1:4 title \"Radio2 g\" w lp\n");
  fprintf(pipe, "unset multiplot"); 
  close(pipe);

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

	  _resta(r1,x[v[i][1]],x[i]);
	  _resta(r2,x[v[i][2]],x[i]);
	  _resta(r3,x[v[i][3]],x[i]);

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

void index_vcnos_prox()
{
  int s1,s2;
  int s1min,s1max,s2min,s2max;
  int dir;
  int s1_dir,s2_dir,s1_dir_old,s2_dir_old;
  int i;

  for(s1_dir=1,s2_dir=0,dir=0; dir<6; dir++)
    {
      switch(s1_dir)
	{
	case -1:
	  s1min=1;
	  s1max=L;
	  break;
	case 0:
	  s1min=0;
	  s1max=L;
	  break;
	case 1:
	  s1min=0;
	  s1max=L-1;
	  break;
	}

      switch(s2_dir)
	{
	case -1:
	  s2min=1;
	  s2max=L;
	  break;
	case 0:
	  s2min=0;
	  s2max=L;
	  break;
	case 1:
	  s2min=0;
	  s2max=L-1;
	  break;
	}

      for(s2=0; s2<L; s2++)
	{
	for(s1=0; s1<L; s1++)
	  {
	    i=s1+L*s2;
	    if(s1>=s1min && s1<s1max && s2>=s2min && s2<s2max)
	      v[i][dir] = i + s1_dir + L*s2_dir;
	    else
	      v[i][dir]=-1;
	  }
	}

      s1_dir_old=s1_dir;
      s2_dir_old=s2_dir;

      s1_dir=-s2_dir_old;
      s2_dir=s1_dir_old + s2_dir_old;
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

  if(v[i][1]!=-1 && v[i][2]!=-1)
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

  if(v[i][0]!=-1 && v[i][1]!=-1)
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
      p=index_plqta_dir0(v[i][3]);
      break;
    case 3:
      p=index_plqta_dir1(v[i][4]);
      break;
    case 4:
      p=index_plqta_dir0(v[i][4]);
      break;
    case 5:
      p=index_plqta_dir1(v[i][5]);
      break;
    }
  return p;
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

