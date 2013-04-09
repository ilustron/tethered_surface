#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define L 16
#define N L*L 
#define M 2*(L-1)*(L-1) 
#define K 10 // Kappa 
#define NF 24000 // nº de archivos
#define B 200 // nº de bloques (Jack-Knife)

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
typedef struct{double xx,xy,yy;} metrica;

metrica gmean();
void index_vcnos_prox();
int sigma_min(int );
int sigma_max(int );
double stimateJK(int , int , double []);
double stimateJKhex(int , int , int ,int ,double [][3][NF]);


vector x[N];
int v[N][6],vord[N][6],dirv_ini[N],zv[N]; // indice primeros vecinos
metrica g[NF];

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
  
  double gxx[NF],Gxx[NF],Gxx2[NF];
  double gyy[NF],Gyy[NF],Gyy2[NF];;
  double gxy[NF],Gxy[NF],Gxy2[NF];;
  double gxxgyy[NF],GxxGyy[NF];
  double sigma1, sigma2;

  // promedios

  double gxxc;
  double gyyc;
  double gxyc;
  double gxxgyyc;
  double s;

  //Errores

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

  output=fopen("indicev1.dat","w"); /* dirv1_ini[i]   zv1[i] */
  for(i=0; i<N; i++)
    {
      fprintf(output,"%4d %4d \n",dirv_ini[i], zv[i]);
    }
  fclose(output);
output=fopen("v_ord.dat","w");
  for(i=0; i<N; i++)
    {
      fprintf(output,"%3d %3d %3d %3d %3d %3d\n",vord[i][0],vord[i][1],vord[i][2],vord[i][3],vord[i][4],vord[i][5]);
    }
  fclose(output);
  // Se leen las posiciones de cada archivo

  f=0;
  sprintf(namein,"xpos_L%d_K10-%d.dat",L,f);
  while((input=fopen(namein,"r"))!=NULL)
    {
      i=0;
      while(fscanf(input,"%lf %lf %lf",&x[i].a,&x[i].b,&x[i].c)!=EOF)
	i++;
      if(i!=N)
	printf("Error: El fichero %s no contiene %d líneas\n",namein,N);      
      fclose(input);      

      // Cálculo de los promedios de las componentes de g
      
      g[f]=gmean();

      //printf("componentes métrica:\n");
      //printf("gxx=%lf gxy=%lf gyy=%lf\n",g[f].xx,g[f].xy,g[f].yy);

      //getchar();
      Gxx[NF-1]+=g[f].xx;//gxx
      Gxx[f]=Gxx[NF-1];
      
      Gxx2[NF-1]+=(g[f].xx*g[f].xx);
      Gxx2[f]=Gxx2[NF-1];

      Gyy[NF-1]+=g[f].yy;//gyy
      Gyy[f]=Gyy[NF-1];
      
      Gyy2[NF-1]+=(g[f].yy*g[f].yy);
      Gyy2[f]=Gyy2[NF-1];

      Gxy[NF-1]+=g[f].xy;//gxy
      Gxy[f]=Gxy[NF-1];
      
      Gxy2[NF-1]+=(g[f].xy*g[f].xy);
      Gxy2[f]=Gxy2[NF-1];

      gxxgyy[f]=g[f].xx*g[f].yy;//gxx*gyy
      GxxGyy[NF-1]+=gxxgyy[f];
      GxxGyy[f]=GxxGyy[NF-1];

      f++;      
      sprintf(namein,"xpos_L%d_K10-%d.dat",L,f);
    }

  if(f!=NF)
    {
      printf("El número de archivos de posiciones leidos no es %d\n",NF);
      return 0;
    } 

  // promedios observables

  gxxc=Gxx2[NF-1]/(double)NF -Gxx[NF-1]*Gxx[NF-1]/(double)NF/(double)NF;
  gyyc=Gyy2[NF-1]/(double)NF -Gyy[NF-1]*Gyy[NF-1]/(double)NF/(double)NF;
  gxyc=Gxy2[NF-1]/(double)NF -Gxy[NF-1]*Gxy[NF-1]/(double)NF/(double)NF;
  gxxgyyc=GxxGyy[NF-1]/(double)NF -Gxx[NF-1]*Gyy[NF-1]/(double)NF/(double)NF;
  sigma1=-(gxxgyyc/gyyc);
  sigma2=-1+2*(gxyc/gxxc);
  s=gxxc-gxxgyyc-2*gxyc;

  // Errores de los observables cartesianos

  for(b=200; b<=200; b++)
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
  int vec_new,vec_old;
  int ord;
  int SI,NO;

  SI=1;
  NO=0;
 
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
		  v[i][dir] = i + basev1[dir].s1 + L * basev1[dir].s2;
		  zv[i]+=1;
		  vec_new=SI;
		}
	      else
		{
		  v[i][dir]=-1;
		  vec_new=NO;
		}
	      if(vec_old-vec_new<0)
		dirv_ini[i]=dir;
	    }
	}
    }

  // Cálculo de vord[N][6]
  for(i=0; i<N; i++)
    {
      for(dir=dirv_ini[i],ord=0; ord<zv[i]; dir++,ord++)
	{
	  if(dir>5)
	    dir=0;
	  
	  vord[i][ord]=v[i][dir];
	}
    }
}
metrica gmean()
{

  int s1,s2;
  int i,j,k,m,np;
  double gxx,gxy,gyy,Gxx,Gxy,Gyy;
  double hx,hy,Hx,Hy,hx_mean,hy_mean;
  metrica g; 
  vector r[6],n;

  
  
  // Recorremos todos los puntos de la red

  Gxx=Gxy=Gyy=0.0F;

  np=0; // contador del número de puntos
  
  for(s2=1; s2<2; s2++)
    {
      for(s1=7; s1<8; s1++)
	{
	  i=s1+L*s2;

	  Hx=Hy=0.0F;

	  np++;

	  for(j=0; j<zv[i]; j++) // Vectores relativos con los vecinos próximos
	    {
	      _resta(r[j], x[vord[i][j]], x[i]);
	      //printf("Componentes del vector posicion %d:\n",vord[i][j]);
	      //printf("%lf %lf %lf \n",x[vord[i][j]].a,x[vord[i][j]].b,x[vord[i][j]].c);
	      //printf("Componentes del vector posicion %d:\n",i);
	      //printf("%lf %lf %lf \n",x[i].a,x[i].b,x[i].c);
	      //printf("Componentes del vector relativo %d:\n",j);
	      //printf("%lf %lf %lf \n",r[j].a,r[j].b,r[j].c);
	      //printf("\n");
	    }

	  m=zv[i]-1;//número de plaquetas (para zv<6)

	  for(k=0; k<m; k++) // Vectores normales a las plaquetas
	    {
	      _prodvec(n,r[k],r[k+1]);
	      hx=-(n.a/n.c);
	      hy=-(n.b/n.c);
	      Hx+=hx;
	      Hy+=hy;
	      //printf("Componentes vector normal %d (no unitario):",k);
	      //printf("%lf %lf %lf",n.a,n.b,n.c);
	      //printf(" hx=%lf hy=%lf \n",hx,hy);
	    }
	  if(zv[i]==6)
	    {
	      _prodvec(n,r[5],r[0]);
	      //printf("Componentes vector normal 6 (no unitario):");
	      hx=-(n.a/n.c);
	      hy=-(n.b/n.c);
	      //printf("%lf %lf %lf ",n.a,n.b,n.c);
	      //printf(" hx=%lf hy=%lf \n",hx,hy);
	      Hx+=hx;
	      Hy+=hy;
	      m=zv[i];
	    }

	  hx_mean=Hx/((double)m);
	  hy_mean=Hy/((double)m);
	  //printf("hx=%lf hy=%lf \n",hx_mean,hy_mean);
	  gxx=1.0+hx_mean*hx_mean;
	  gxy=hx_mean*hy_mean;
	  gyy=1.0F+hy_mean*hy_mean;
	  //printf("gxx=%lf gxy=%lf gyy=%lf\n",gxx,gxy,gyy);
	  Gxx+=gxx;
	  Gxy+=gxy;
	  Gyy+=gyy;

	}
    }
  //printf("%lf %lf %lf \n",Gxx,Gxy,Gyy);
  g.xx=Gxx/(double)np;
  g.xy=Gxy/(double)np;
  g.yy=Gyy/(double)np;
  return g;
}
