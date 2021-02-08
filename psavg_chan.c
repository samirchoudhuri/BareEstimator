//To calculate the mean and variance of the power spectrum from 100 files
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>

int main(int argc, char *argv[])
{
  int ii,jj,N,Nbin;
  double *Uval,*Vval,*ps,*psseq,dummyps,dummyu,*num;
  char in_file[150];
  FILE *fp;

  if(argc!=5)
    {
      printf("Usage: %s <No. of file> <Nbin> <ouput file> <inputfile name>\n", argv[0]);
      return 1;
    }
  
  N=atoi(argv[1]);
  Nbin=atoi(argv[2]);
  
  Uval= (double*)calloc(Nbin,sizeof(double*));
  Vval= (double*)calloc(Nbin,sizeof(double*));
  ps= (double*)calloc(Nbin,sizeof(double*));
  psseq= (double*)calloc(Nbin,sizeof(double*));
  num= (double*)calloc(Nbin,sizeof(double*));

  fp=fopen(argv[4],"r");
  for(ii=0;ii<N;++ii)
    {
      for(jj=0;jj<Nbin;jj++)
	{
	  fscanf(fp,"%lf%lf%*f%lf",&dummyu,&dummyps,&num[jj]);
	  Uval[jj]+=dummyu;
	  ps[jj]+=dummyps;
	  psseq[jj]+=(dummyps*dummyps);
	}
      
    }
  fclose(fp);  
  
  fp = fopen(argv[3],"w");
  for(jj=0;jj<Nbin;++jj)
    {
      Uval[jj]=Uval[jj]/(1.*N);
      ps[jj]=ps[jj]/(1.*N);
      psseq[jj]=(psseq[jj]/(1.*N))-(ps[jj]*ps[jj]);
      //fprintf(fp,"%e %e %e %e %e\n",Uval[jj],Vval[jj],ps[jj],sqrt(psseq[jj]),num[jj]);
      fprintf(fp,"%e %e %e %e\n",Uval[jj],ps[jj],sqrt(psseq[jj]),num[jj]);
    }
  fclose(fp);
}
