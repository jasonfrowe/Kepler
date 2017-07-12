CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine windowfn(npt,nfit,time,mag,merr,period,per1,per2,avg,a,
     .   nb,nc,steps,panx,pany,nstar,pers)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Attempts to model the window function for the data by using shape
C     parameters from COSINEFIT.  A noiseless light curve is sampled 
C     using the TIME(npt) data.  The model data is then run through 
C     PDM to produce a periodogram that approximates the window 
C     function.
C     
      implicit none
      integer nmax,npt,nfmax,nfit,i,j,nb,nc,steps,panx,pany,nstarmax,
     .   nstar,i2,k,cc,na
      parameter(nmax=2000000,nstarmax=300)
      double precision time(npt),mag(npt),per1,per2,avg,a(nfit),temp,
     .  res(nmax),bper,period,w(nstarmax),pi,btheta,merr(npt),
     .  pers(nstarmax),arg,sn,snlimit

      snlimit=3.6 !S/N for jmfourw (useless for window function)
      pi=3.141592654
      do 5 i=1,nstar
         w(i)=(2.0*pi)/pers(i)
 5    continue
      na=nstar*nfit*2+nstar
      
      j=(na-nstar)/nstar
      
      do 10 i2=1,npt
         temp=0.0
         do 20 k=1,nstar
            cc=(j+1)*(k-1)+2
            do 21 i=cc,cc+j-2,2
               arg=dble((i-cc+2)/2)*a(cc-1)*time(i2)+a(i+1)
               temp=temp+a(i)*cos(arg)
 21         continue
 20      continue
         res(i2)=temp
 10   continue


      call jmfourw(npt,time,res,merr,per1,per2,steps,bper,btheta,sn,
     .   panx,pany,snlimit,1)
c      call pdm(npt,time,res,nb,nc,per1,per2,steps,bper,btheta,panx,pany)
C      call plotph(npt,time,res,period,1,4)

      return
      end
