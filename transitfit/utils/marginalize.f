CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine marginalize(npt,time,mag,merr,itime,tmodel,nfit,sol,
     .  serr,Dpvary,npar,sol2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nfit,npar,i,k
      double precision time(npt),mag(npt),merr(npt),itime(npt),
     .  tmodel(npt),sol(nfit),serr(nfit,2),Dpvary(nfit),msum,likelihood,
     .  sol2(nfit),plow,phgh,msumold

     
      k=6
      msum=0.0d0
      if((k.ne.npar).and.(serr(k,2).ne.0.0))then
        msum=msum+likelihood(npt,tmodel,time,mag,merr,
     .      itime,nfit,sol,serr)
        msumold=msum
        do 10 i=1,nfit
            sol2(i)=sol(i)
 10     continue
        plow=sol(k)
        phgh=sol(k)
C       start loop
        do 11 i=1,10000
            plow=plow-Dpvary(k)
            phgh=phgh+Dpvary(k)
            sol2(k)=plow
            msum=likelihood(npt,tmodel,time,mag,merr,
     .        itime,nfit,sol2,serr)
            write(6,*) plow,msum
            sol2(k)=phgh
            msum=likelihood(npt,tmodel,time,mag,merr,
     .        itime,nfit,sol2,serr)
            write(6,*) phgh,msum
c            write(6,*) "msum:",msum/msumold
            msumold=msum
 11     continue
      endif
     
      return
      end