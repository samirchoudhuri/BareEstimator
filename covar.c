//it will generate corelation coefficient
# include <stdio.h>
# include <math.h>

main()
{
  int ii,jj,kk,Nbin=10,Nloop=20;
  FILE *fp;
  double ps[Nbin][Nloop],pstot1,pstot2,pstot12,p11,p22,covar[Nbin][Nbin],r12[Nbin][Nbin];

  fp=fopen("ps_rndm2nonoise","r");
  for(ii=0;ii<Nbin;++ii)
    {
      fscanf(fp,"%*lf");
      for(jj=0;jj<Nloop;++jj)
	{
	  fscanf(fp,"%lf",&ps[ii][jj]);
	}
    }
  fclose(fp);

  for(ii=0;ii<Nbin;++ii)
    for(jj=0;jj<Nbin;++jj)
      {
	pstot1=0.;
	pstot2=0.;
	pstot12=0.;
	for(kk=0;kk<Nloop;++kk)
	  {
	    p11=ps[ii][kk];
	    p22=ps[jj][kk];
	    
	    pstot1+=p11;
	    pstot2+=p22;
	    pstot12+=(p11*p22);
  	   }
	pstot1=pstot1/(1.*Nloop);
	pstot2=pstot2/(1.*Nloop);
	pstot12=pstot12/(1.*Nloop);
	covar[ii][jj]=(pstot12-pstot1*pstot2);
	//fprintf(fp,"%d\t%d\t%e\n",ii,jj,(pstot12-pstot1*pstot2));
      }
  
  fp=fopen("correff_rndm2.dat","w");
  for(ii=0;ii<Nbin;++ii)
    for(jj=0;jj<Nbin;++jj)
       {
  	 r12[ii][jj]=covar[ii][jj]/sqrt(covar[ii][ii]*covar[jj][jj]);
  	 if(ii==jj) r12[ii][jj]=0.;
	 fprintf(fp,"%d\t%d\t%e\n",ii,jj,r12[ii][jj]);
       }
  fclose(fp);
  
  
}
