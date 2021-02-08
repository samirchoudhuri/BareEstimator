 //variance calculation for 1st estimator
//./<executable> visinput100k
# include <stdio.h>
# include <math.h>
# include <stdlib.h>

double A_150=5.13e-4,alpha=2.8,betav=2.34;
double k_B=1.38e3 /* in Jy*/,cspeed=3.0e8,pi=3.141593;

int main(int argc, char *argv[])
{
  char fname[40];
  int ii,Nbin,Nv;
  long gcount,nbasln,maxbasln=8000;
  double signv,*uu,*vv,*uread,*vread,DD,Umin,kv,Umax,Uval,nu,pu,SA,SB,SC;
  double u1,u2,v1,v2,uval;
  double *ubin,*psmean,*sim_var,*wt;
  double sigma,sigmau=16.6;
  unsigned long *uindex,*vindex,total=0;
  float *uufloat,*vvfloat;

  int *mval,ln,lmax,a,b,c,d;
  double **wvar,**rvar,A,theta_0=2.75e-2,deluab2,*variance,fac1;

  FILE *fp,*fp1;

  system ("date");

  if(argv[1]==NULL)
  { fprintf(stderr,"Usage: argv[0] Input file\n");return 1;}
  
  A=M_PI*pow(theta_0,2.)/2.;

  fp=fopen(argv[1],"r");
  fscanf(fp,"%*d%ld%*ld%lf%*lf",&gcount,&nu);
  fscanf(fp,"%*d%*f%s",fname);
  fscanf(fp,"%lf%lf%lf%d%*d%lf",&Umin,&Umax,&DD,&Nbin,&sigma);
  fclose(fp);

  printf("%ld %s %lf %lf %lf %d %lf A=%e\nsigma=%lf\n",gcount,fname,Umin,Umax,DD,Nbin,nu,A,sigma);

  ubin=(double*)calloc(Nbin,sizeof(double));
  psmean=(double*)calloc(Nbin,sizeof(double));
  sim_var=(double*)calloc(Nbin,sizeof(double));
  wt=(double*)calloc(Nbin,sizeof(double));
  
  fp=fopen("psmean_gmrt2","r");
  for(ii=0;ii<Nbin;++ii)
    fscanf(fp,"%lf%lf%lf%lf",&ubin[ii],&psmean[ii],&sim_var[ii],&wt[ii]);
  fclose(fp);

  printf("%e %e %e %e\n",ubin[Nbin-1],psmean[Nbin-1],sim_var[Nbin-1],wt[Nbin-1]);
  
  uufloat=(float*)calloc(gcount,sizeof(float));
  vvfloat=(float*)calloc(gcount,sizeof(float));
  uu =(double*) malloc(sizeof(double)*gcount);
  if(uu==NULL)
    { fprintf(stderr,"malloc failurs33\n"); return 1;}
  vv = (double*)malloc(sizeof(double)*gcount);
  if(vv==NULL)
    { fprintf(stderr,"malloc failurs\n"); return 1;}
  uread =(double*) malloc(sizeof(double)*gcount);
  if(uread==NULL)
    { fprintf(stderr,"malloc failurs33\n"); return 1;}
  vread = (double*)malloc(sizeof(double)*gcount);
  if(vread==NULL)
    { fprintf(stderr,"malloc failurs\n"); return 1;}

  fp=fopen(fname,"r");
  fread(uread,sizeof(double),gcount,fp);
  fread(vread,sizeof(double),gcount,fp);
  fclose(fp);

  nbasln=0;
  for(ii=0;ii<gcount;ii++)
    {
      uval=sqrt(uread[ii]*uread[ii]+vread[ii]*vread[ii]);
      if((uval>=Umin)&&(uval<=Umax))
  	{
  	  signv= (vread[ii]<0.) ? -1. : 1. ;
  	  uu[nbasln]=signv*uread[ii];
  	  vv[nbasln]=signv*vread[ii];
	  ++nbasln;
	}
    }

  /* index data according to u and v */
  total=(unsigned long)nbasln;
  uindex=(unsigned long*)malloc((size_t)(sizeof(unsigned long))*nbasln);
  vindex=(unsigned long*)malloc((size_t)(sizeof(unsigned long))*nbasln);
  
  //copy into float variable
  for(ii=0;ii<nbasln;++ii)
    {
      uufloat[ii]=uu[ii];
      vvfloat[ii]=vv[ii];
    }
  
  indexx(total,uufloat-1,uindex-1);
  indexx(total,vvfloat-1,vindex-1);
 
  variance=(double*)calloc(Nbin,sizeof(double));
  mval=(int*)calloc(maxbasln,sizeof(int));
  wvar=(double**)calloc(maxbasln,sizeof(double*));
  rvar=(double**)calloc(maxbasln,sizeof(double*));
  wvar[0]=(double*)calloc(maxbasln*maxbasln,sizeof(double));
  rvar[0]=(double*)calloc(maxbasln*maxbasln,sizeof(double));

  for(ii=1;ii<maxbasln;++ii)
    {
      wvar[ii]=wvar[0]+ii*maxbasln;
      rvar[ii]=rvar[0]+ii*maxbasln;
    }

  // logaritmic bins
  kv=Nbin/log10(Umax/Umin);

  lmax=0;
  for(ii=0;ii<nbasln;++ii)
    {
      u1=uu[uindex[ii]-1];v1=vv[uindex[ii]-1];
      Uval=sqrt(u1*u1+v1*v1);
      Nv=(int)floor(kv*log10(Uval/Umin));//Nbin
      if(Nv<Nbin) //for appropriate end of the bin
  	{
	  lmax=ii;
	  ln=0;
	  while((lmax<nbasln) && ((uu[uindex[lmax]-1]-u1)<=1.*DD)) 
	    {
	      if((pow(uu[uindex[lmax]-1]-u1,2.)+pow(vv[uindex[lmax]-1]-v1,2.))<1.*DD*DD)
		{
		  mval[ln]=uindex[lmax]-1;
		  ++ln;
		  if(ln>maxbasln) 
		    {
		      printf("increase ln range\n");
		      exit(0);
		    }
		}
	      ++lmax;
	    }
	  //half circle searching complete
	  pu=pow(2.*k_B*nu*nu*1.0e12/(cspeed*cspeed),2.)*A_150*pow(150.0/nu,2.0*alpha)*pow(1000./(2.*pi*Uval),betav);
	  //pu=psmean[Nv];
	  for(a=0;a<ln;++a)
	    {
	      wvar[a][a]=0.0;rvar[a][a]=(pu+(2.*pow(sigma,2.)/(1.*A)));
	      for(b=a+1;b<ln;++b)
		{
		  deluab2=pow(uu[mval[a]]-uu[mval[b]],2.)+pow(vv[mval[a]]-vv[mval[b]],2);
		  fac1=exp(-1.*deluab2/(sigmau*sigmau));
		  wvar[a][b]=1.*fac1;
   		  wvar[b][a]=wvar[a][b];
		  rvar[a][b]=1.*pu*fac1;
		  rvar[b][a]=rvar[a][b];
		}
	    }
	  
	  for(b=1;b<ln;++b) 
	    {
	      variance[Nv]+=2.*wvar[0][b]*wvar[0][b]*(rvar[0][0]*rvar[b][b]+rvar[0][b]*rvar[0][b]);
	      for(c=b+1;c<ln;++c)
		{
		  variance[Nv]+=4.*(wvar[0][b]*wvar[0][c]*(rvar[0][b]*rvar[0][c]+rvar[0][0]*rvar[b][c])+wvar[0][b]*wvar[b][c]*(rvar[0][c]*rvar[b][b]+rvar[0][b]*rvar[b][c])+wvar[0][c]*wvar[b][c]*(rvar[0][c]*rvar[b][c]+rvar[0][b]*rvar[c][c]));
		  for(d=c+1;d<ln;++d)
		    {
		      SA=(wvar[0][c]*wvar[b][d]+wvar[0][d]*wvar[b][c])*rvar[0][b]*rvar[c][d];
		      SB=(wvar[0][b]*wvar[c][d]+wvar[0][d]*wvar[b][c])*rvar[0][c]*rvar[b][d];
		      SC=(wvar[0][c]*wvar[b][d]+wvar[0][b]*wvar[c][d])*rvar[0][d]*rvar[b][c];
		      variance[Nv]+=4.*(SA+SB+SC);
		    }
		}
	    }
	}
    }

fp=fopen("var_gmrt2new","w");
  for(ii=0;ii<Nbin;++ii)
    {
      variance[ii]/=(4.*wt[ii]*wt[ii]);
      fprintf(fp,"%e\t%e\t%e\t%e\n",ubin[ii],psmean[ii],sim_var[ii],sqrt(variance[ii]));
    }
fclose(fp);
  system ("date");

}//main close
