#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define L 192
#define N L*L 
#define M 2*(L-1)*(L-1) 
#define K 2.0 // Kappa 
#define NF 1000 // nº de archivos
#define B 10 // tamaño del bloque (Jack-Knife)

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

double gmean(int dir1,int dir2);
void index_vcnos_prox();
int sigma_min(int si);
int sigma_max(int si);
double stimateJK(int , int , double []);
double stimateJKhex(int , int , int ,int ,double [][3][NF]);

vector x[N];
int v[N][6]; // indice primeros vecinos

struct vector2D{
int s1,s2;
} basev1[6];

struct rect{
struct vector2D min,max;
} rombo[6];
 
int main(void)
{
  FILE *input;
  FILE *output;
  
  char namein[255];
  char nameout[255];

  int i,j,k,l,b,n,f,dir1,dir2,s1,s2;  

  // observables
  double g[3][3][NF],G[3][3][NF],G2[3][3][NF];
  double gxx[NF],Gxx[NF],Gxx2[NF];
  double gyy[NF],Gyy[NF],Gyy2[NF];;
  double gxy[NF],Gxy[NF],Gxy2[NF];;
  double gxxgyy[NF],GxxGyy[NF];
  double sigma1, sigma2;

  // promedios
  double gc[3][3];
  double gxxc;
  double gyyc;
  double gxyc;
  double gxxgyyc;
  double s;

  //Errores
  double error_gc[3][3];
  double error_gxxc;
  double error_gxyc;
  double error_gyyc;
  double error_gxxgyyc;
  double error_sigma1,error_sigma2;
  double error_s;

  double m1JK_g[3][3];
  double m2JK_g[3][3];
  double m1JK_gxx,m1JK_gxy,m1JK_gyy,m1JK_gxxgyy;
  double m2JK_gxx,m2JK_gxy,m2JK_gyy;

  double gcJK[3][3];
  double gxxcJK,gxycJK,gyycJK,gxxgyycJK;
  double sigma1JK,sigma2JK;
  double sJK;


  // Cargamos los índices v[N][6]

  index_vcnos_prox();

  // Se leen las posiciones de cada archivo

  f=0;
  sprintf(namein,"./RUNS/L%d/K%.1f/xpos_L%d_K%.1f-%d.dat",L,K,L,K,f);
  
  while((input=fopen(namein,"r"))!=NULL)
    {
      i=0;
      while(fscanf(input,"%lf %lf %lf",&x[i].a,&x[i].b,&x[i].c)!=EOF)
	i++;
      if(i!=N)
	printf("Error: El fichero %s no contiene %d líneas\n",namein,N);      
      fclose(input);
      for(dir1=0; dir1<3; dir1++)
	{
	  for(dir2=dir1; dir2<3; dir2++)
	    {
	      g[dir1][dir2][f]=gmean(dir1,dir2);

	      G[dir1][dir2][NF-1]+=g[dir1][dir2][f];
	      G[dir1][dir2][f]=G[dir1][dir2][NF-1];

	      G2[dir1][dir2][NF-1]+=(g[dir1][dir2][f]*g[dir1][dir2][f]);
	      G2[dir1][dir2][f]=G2[dir1][dir2][NF-1];
	    }
	}
      // Componentes cartesianas en función de las configuraciones
    
      gxx[f]=g[0][0][f];//gxx
      Gxx[NF-1]+=gxx[f];
      Gxx[f]=Gxx[NF-1];
      
      Gxx2[NF-1]+=(gxx[f]*gxx[f]);
      Gxx2[f]=Gxx2[NF-1];

      gyy[f]=(g[1][1][f]+g[2][2][f]+2*g[1][2][f])/3.0;//gyy
      Gyy[NF-1]+=gyy[f];
      Gyy[f]=Gyy[NF-1];
      
      Gyy2[NF-1]+=(gyy[f]*gyy[f]);
      Gyy2[f]=Gyy2[NF-1];

      gxy[f]=(g[0][1][f]+g[0][2][f])/sqrt(3.0);//gxy
      Gxy[NF-1]+=gxy[f];
      Gxy[f]=Gxy[NF-1];
      
      Gxy2[NF-1]+=(gxy[f]*gxy[f]);
      Gxy2[f]=Gxy2[NF-1];

      gxxgyy[f]=gxx[f]*gyy[f];//gxx*gyy
      GxxGyy[NF-1]+=gxxgyy[f];
      GxxGyy[f]=GxxGyy[NF-1];

      f++;      
      sprintf(namein,"./RUNS/L%d/K%.1f/xpos_L%d_K%.1f-%d.dat",L,K,L,K,f);
    }

  if(f!=NF)
    {
      printf("El número de archivos de posiciones leidos no es %d\n",NF);
      return 0;
    } 
  
  // Promedios observables
  for(s1=0; s1<3; s1++)
	{
	  for(s2=s1; s2<3; s2++)
	    gc[s1][s2]=G2[s1][s2][NF-1]/(double)NF-G[s1][s2][NF-1]*G[s1][s2][NF-1]/(double)NF/(double)NF;
    	}
  gxxc=Gxx2[NF-1]/(double)NF -Gxx[NF-1]*Gxx[NF-1]/(double)NF/(double)NF;
  gyyc=Gyy2[NF-1]/(double)NF -Gyy[NF-1]*Gyy[NF-1]/(double)NF/(double)NF;
  gxyc=Gxy2[NF-1]/(double)NF -Gxy[NF-1]*Gxy[NF-1]/(double)NF/(double)NF;
  gxxgyyc=GxxGyy[NF-1]/(double)NF -Gxx[NF-1]*Gyy[NF-1]/(double)NF/(double)NF;
  sigma1=-(gxxgyyc/gyyc);
  sigma2=-1+2*(gxyc/gxxc);
  s=gxxc-gxxgyyc-2*gxyc;

  //Error de los observables hexagonales

  for(b=B; b<=B; b++)
    {
      if((NF%b)==0)
	{
	  n=NF/b;

	  for(s1=0; s1<3; s1++)
	    {
	      for(s2=s1; s2<3; s2++)
		{
		  error_gc[s1][s2]=0.0F;
		  for (k=0; k<b; k++)
		    {
		      m1JK_g[s1][s2]=stimateJKhex(n, k, s1, s2, G);
		      m2JK_g[s1][s2]=stimateJKhex(n, k, s1, s2, G2);
		      gcJK[s1][s2]=m2JK_g[s1][s2]-m1JK_g[s1][s2]*m1JK_g[s1][s2];
		      error_gc[s1][s2]+=pow(gc[s1][s2]-gcJK[s1][s2],2.0);
		    }
		   error_gc[s1][s2]=(double)(b-1)/((double) b) * error_gc[s1][s2];
		   error_gc[s1][s2]=sqrt(error_gc[s1][s2]);
		}
	    }
	}
    }
  
  // Errores de los observables cartesianos

  for(b=B; b<=B; b++)
    {
      if((NF%b)==0)
	{
	  n=NF/b;

	  error_gxxc=0.0F;
	  error_gxyc=0.0F;
	  error_gyyc=0.0F;
	  error_gxxgyyc=0.0F;
	  error_sigma1=0.0F;
	  error_sigma2=0.0F;

	  for (k=0; k<b; k++)
	    {
	      m1JK_gxx=stimateJK(n, k, Gxx);
	      m1JK_gyy=stimateJK(n, k, Gyy);
	      m1JK_gxy=stimateJK(n, k, Gxy);
	      m1JK_gxxgyy=stimateJK(n, k, GxxGyy);

	      m2JK_gyy=stimateJK(n, k, Gyy2);
	      m2JK_gxx=stimateJK(n, k, Gxx2);	      
	      m2JK_gxy=stimateJK(n, k, Gxy2);	      

	      gxxcJK=m2JK_gxx-m1JK_gxx*m1JK_gxx;
	      gyycJK=m2JK_gyy-m1JK_gyy*m1JK_gyy;
	      gxycJK=m2JK_gxy-m1JK_gxy*m1JK_gxy;
	      gxxgyycJK=m1JK_gxxgyy-m1JK_gxx*m1JK_gyy;

	      sigma1JK=-gxxgyycJK/gyycJK;
	      sigma2JK=-1+2*(gxycJK/gxxcJK);
	      sJK=gxxcJK-gxxgyycJK-2*gxycJK;
	      
	      error_gxxc+=pow(gxxc-gxxcJK,2.0);
	      error_gyyc+=pow(gyyc-gyycJK,2.0);
	      error_gxyc+=pow(gxyc-gxycJK,2.0);
	      error_gxxgyyc+=pow(gxxgyyc-gxxgyycJK,2.0);
	      error_sigma1+=pow(sigma1-sigma1JK,2.0);
	      error_sigma2+=pow(sigma2-sigma2JK,2.0);
	      error_s+=pow(s-sJK,2.0);
		      
	    }
	  error_gxxc=(double)(b-1)/((double) b) * error_gxxc;
	  error_gxxc=sqrt(error_gxxc);

	  error_gyyc=(double)(b-1)/((double) b) * error_gyyc;
	  error_gyyc=sqrt(error_gyyc);

	  error_gxyc=(double)(b-1)/((double) b) * error_gxyc;
	  error_gxyc=sqrt(error_gxyc);

	  error_gxxgyyc=(double)(b-1)/((double) b) * error_gxxgyyc;
	  error_gxxgyyc=sqrt(error_gxxgyyc);

	  error_sigma1=(double)(b-1)/((double) b) * error_sigma1;
	  error_sigma1=sqrt(error_sigma1);

	  error_sigma2=(double)(b-1)/((double) b) * error_sigma2;
	  error_sigma2=sqrt(error_sigma2);


	  
	}
    }
  
  printf("<gii>c para direciones hexagonales:\n");
  printf("\t g11=%lf+-%lf\n",gc[0][0],error_gc[0][0]);
  printf("\t g12=%lf+-%lf\n",gc[0][1],error_gc[0][1]);
  printf("\t g13=%lf+-%lf\n",gc[0][2],error_gc[0][2]);
  printf("\t g22=%lf+-%lf\n",gc[1][1],error_gc[1][1]);
  printf("\t g23=%lf+-%lf\n",gc[1][2],error_gc[1][2]);
  printf("\t g33=%lf+-%lf\n",gc[2][2],error_gc[2][2]);
  
  printf("Direcciones ortogonales (n=%d)\n",n);
  printf("\t <gxx²>c=%lf+-%lf\n",gxxc,error_gxxc);
  printf("\t <gxy²>c=%lf+-%lf\n",gxyc,error_gxyc);
  printf("\t <gyy²>c=%lf+-%lf\n",gyyc,error_gyyc);
  printf("\t <gxxgyy>c=%lf+-%lf\n",gxxgyyc,error_gxxgyyc);
  printf("\t sigma1=%lf+-%lf\n",sigma1,error_sigma1);
  printf("\t sigma2=%lf+-%lf\n",sigma2,error_sigma2);
  printf("\t s=%lf+-%lf\n",s,error_s);

  return 1;
}

double stimateJK(int n, int k, double Obs[])
{
  double obsJK,sumJK_obs;

  if(k==0)
    {
      obsJK=Obs[n-1];
      sumJK_obs=(Obs[NF-1]-obsJK)/(double)(NF-n);
    }
  else
    {
      obsJK=Obs[(k+1)*n-1]-Obs[k*n-1];
      sumJK_obs=(Obs[NF-1]-obsJK)/(double)(NF-n);
    }
  return sumJK_obs;
}
double stimateJKhex(int n, int k, int s1,int s2,double Obs[][3][NF])
{
  double obsJK,sumJK_obs;

  if(k==0)
    {
      obsJK=Obs[s1][s2][n-1];
      sumJK_obs=(Obs[s1][s2][NF-1]-obsJK)/(double)(NF-n);
    }
  else
    {
      obsJK=Obs[s1][s2][(k+1)*n-1]-Obs[s1][s2][k*n-1];
      sumJK_obs=(Obs[s1][s2][NF-1]-obsJK)/(double)(NF-n);
    }
  return sumJK_obs;
}


double gmean(int dir1,int dir2)
{

  int s1,s2;
  int i,n;
  double g,G; 
  vector r1,r2;

  n=0;
  G=0;
  // Recorremos los puntos de la red que tienen vecinos en dir 0, 1 y 2  
  for(s2=2; s2<L-2; s2++)
    {
      for(s1=2; s1<L-2; s1++)
	{
	  i=s1+L*s2;
      
	  _resta(r1, x[v[i][dir1]], x[i]);
	  _resta(r2, x[v[i][dir2]], x[i]);
	  g=_prodesc(r1,r2);
	  G+=g;
	  n++;
	}
    }

  G/=(double)n;
  return G;
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

void index_vcnos_prox()
{
  int s1,s2;
  int s1min,s1max,s2min,s2max;
  int dir;
  int s1_dir,s2_dir,s1_old,s2_old;
  int i;

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

  for(s2=0; s2<L; s2++)
    {
      for(s1=0; s1<L; s1++)
	{
	  i=s1+L*s2;
	  for(dir=0; dir<6; dir++)
	    {
	      if(s1>=rombo[dir].min.s1 && s1<rombo[dir].max.s1 && s2>=rombo[dir].min.s2 && s2<rombo[dir].max.s2)
		v[i][dir] = i + basev1[dir].s1 + L * basev1[dir].s2;
	      else
		v[i][dir]=-1;
	    }
	}
    }  
}
