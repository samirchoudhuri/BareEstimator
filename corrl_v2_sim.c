// This program caclulates visibility correlations Bare Estimator
//./<executable> <input FITS file> <ouput .dat file> <input>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <fitsio.h>
# include <unistd.h>
# include <nr.h>
# include <nrutil.h>
# include "read_fits_func.h"

void printerror(int status){
  if (status){
    fits_report_error(stderr, status);
    exit( status ); 
  }
}

extern float *uu,*vv,**RRre,**RRim,**LLre,**LLim,del_chan,nu_chan0,chan0;
extern long nstokes,nchan,ncmplx,gcount;
extern int flagLL;

int main(int argc, char *argv[])
{ 
  char INFITS[128],OUTFITS[128],input[128];
  int chan,chan1,chan2,n_ave,ii;
  int nbasln,Nbin;
  double Umax,Umin;
  double sigman,sigmau,DD,du0,theta_0,A;
  double **corr_re,**corr_im,**u_val,**num_v,kv;
  unsigned long *uindex,*vindex,total=0; 
  int jj,nn,mm,kk,ll,Nv;
  double u1,u2,v1,v2,duv,Uval,wt,Re,Im,Re1,Im1; 
  FILE *fp;

  if(argc!=4)
    {
      printf("Usage: %s <input FITS file> <ouput .dat file> <input>\n", argv[0]);
      return 1;
    }
  
  sscanf(argv[1],"%s",INFITS);
  sscanf(argv[2],"%s",OUTFITS);
  sscanf(argv[3],"%s",input);
 
  fp=fopen(input,"r"); 
  fscanf(fp,"%d%d",&chan1,&chan2);
  fscanf(fp,"%lf%lf%lf",&Umax,&Umin,&theta_0);
  fscanf(fp,"%lf%lf%lf%d",&DD,&sigman,&sigmau,&Nbin);
  fscanf(fp,"%d",&flagLL);
  fclose(fp);
  
  theta_0=(M_PI*theta_0)/(180.*60.);//theta_0 in rad
  A=M_PI*pow(theta_0,2.)/2.;

  chan1=chan1-1; chan2=chan2-1;
  n_ave=chan2-chan1+1;

  printf("\n\nchan1-1=%d\tchan2-1=%d\tn_ave=%d\nUmax=%lf\tUmin=%lf\n\n",chan1,chan2,n_ave,Umax,Umin);
  printf("theta_0=%erad A=%e\n",theta_0,A);
  printf("DD=%e sigman=%e sigmau=%e Nbin=%d\n",DD,sigman,sigmau,Nbin);

 if(access(INFITS, F_OK)!=0){
    printf("Input File %s does not exists\n",INFITS);
    exit(0);
  }

  if(access(argv[4], F_OK)==0){
    printf("Output File %s exists\n",argv[4]);
    exit(0);
  }

  //It will fill RRre[][],RRim[][],LLre[][],LLim[][]
  nbasln=readfits(INFITS,Umax,Umin,chan1,chan2);
  
  printf("\nchan0=%e nu_chan0=%e del_chan=%e\n",chan0,nu_chan0,del_chan);
  printf("nbasln=%d\n",nbasln); 
  
  double *corrfact;
  corrfact = (double*) calloc(nchan, sizeof(double));
  for(ii=0; ii<nchan; ii++)
    {    
      corrfact[ii] = 1.+ (del_chan/nu_chan0)*(ii+1.+0.5-chan0);//lamda[chan0]/lambda[ii]
    }
  
  // logaritmic bins
  kv=Nbin/log10(Umax/Umin);

  corr_re=(double**)calloc(Nbin,sizeof(double*));
  corr_im=(double**)calloc(Nbin,sizeof(double*));
  num_v=(double**)calloc(Nbin,sizeof(double*));
  u_val=(double**)calloc(Nbin,sizeof(double*));
  
  corr_re[0]=(double*)calloc(Nbin*n_ave,sizeof(double));
  corr_im[0]=(double*)calloc(Nbin*n_ave,sizeof(double));
  u_val[0]=(double*)calloc(Nbin*n_ave,sizeof(double));
  num_v[0]=(double*)calloc(Nbin*n_ave,sizeof(double*));

  for(ii=1;ii<Nbin;++ii)
    {
      corr_re[ii]=corr_re[0]+ii*n_ave;
      corr_im[ii]=corr_im[0]+ii*n_ave;
      u_val[ii]=u_val[0]+ii*n_ave;
      num_v[ii]=num_v[0]+ii*n_ave;
    }

  /* index data according to u and v */
  total=(unsigned long)nbasln;
  uindex=(unsigned long*)malloc((size_t)(sizeof(unsigned long))*nbasln);
  vindex=(unsigned long*)malloc((size_t)(sizeof(unsigned long))*nbasln);

  indexx(total,uu-1,uindex-1);
  indexx(total,vv-1,vindex-1);

 

  printf("starting correlation\n");
  for(ii=0;ii<nbasln;++ii)
    {
      du0=0.;jj=1;
      kk=uindex[ii]-1;
      while((du0<=1.1*DD)&&(ii+jj<nbasln))
  	{
  	  mm=uindex[ii+jj]-1;
  	  du0=(uu[mm]-uu[kk]);

	  for(chan=0;chan<n_ave;++chan)
	    {
	      u1=uu[kk]*corrfact[chan+chan1];v1=vv[kk]*corrfact[chan+chan1];
	      u2=uu[mm]*corrfact[chan+chan1]; v2=vv[mm]*corrfact[chan+chan1];
	      duv=sqrt((u1-u2)*(u1-u2)+(v1-v2)*(v1-v2));

	      //add if both stokes are unflagged
	      if((duv<=DD)&&(RRre[kk][chan]>-1.e7)&&(RRre[mm][chan]>-1.e7)&&(LLre[kk][chan]>-1.e7)&&(LLre[mm][chan]>-1.e7))
		{
		  Uval=0.5*(sqrt(u1*u1+v1*v1)+sqrt(u2*u2+v2*v2));
		  Nv=(int)floor(kv*log10(Uval/Umin));
		  
		  if((Nv>=0)&&(Nv<Nbin))
		    {
		      wt=exp(-1.*duv*duv/(sigmau*sigmau));
		      
		      Re=(RRre[kk][chan]+LLre[kk][chan])/2.;
		      Im=(RRim[kk][chan]+LLim[kk][chan])/2.;
		      
		      Re1=(RRre[mm][chan]+LLre[mm][chan])/2.;
		      Im1=(RRim[mm][chan]+LLim[mm][chan])/2.;
		      
		      corr_re[Nv][chan]+=(wt*(Re*Re1+Im*Im1));
		      corr_im[Nv][chan]+=(wt*(Im*Re1-Re*Im1));
		      u_val[Nv][chan]+=(Uval*wt*wt);
		      num_v[Nv][chan]+=(wt*wt);
		    }
		}
	    }
	  ++jj;
	}
    }
  printf("correlation done \n");
  
  fp=fopen(OUTFITS,"w");
  for(chan=0;chan<n_ave;++chan)
    {
      for(ii=0;ii<Nbin;++ii)
	{
	  if(num_v[ii][chan]>0.)
	    {
	      u_val[ii][chan]/=(1.*num_v[ii][chan]);
	      corr_re[ii][chan]/=(1.*num_v[ii][chan]*A);
	      corr_im[ii][chan]/=(1.*num_v[ii][chan]);
	      fprintf(fp,"%e %e %e %e\n",u_val[ii][chan],corr_re[ii][chan],corr_im[ii][chan],num_v[ii][chan]);
	    }
	}
      fprintf(fp,"\n");
    }
  fclose(fp);
 
  free(corrfact);
  free(corr_re);
  free(corr_im);
  free(num_v);
  free(u_val);
  return 0;
}
