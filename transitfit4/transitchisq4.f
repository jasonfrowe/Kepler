      program transitchisq
      implicit none
      integer iargc,nmax,nunit,npt,nfit,i,j
      parameter(nmax=2000000,nfit=18)
      integer dtype(nmax)
      double precision time(nmax),flux(nmax),ferr(nmax),itime(nmax),
     .  Keplertime,sol(nfit),serr(nfit,2),err(nfit),tdur,transitdur,
     .  tmodel(nmax),pcut(nmax),epo,per,ph1,ph2,ph,std,stdev,mean,
     .  chisq
      character*80 obsfile,inputsol
      
      if(iargc().lt.2) goto 901

C     Parse the name of the observations data file from the commandline
      call getarg(1,obsfile)
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,flux,ferr,itime,Keplertime)
      call readkeplc(nunit,nmax,npt,time,flux,ferr,itime,Keplertime)
      close(nunit)!release unit number as we are done with the file.
      
      call getarg(2,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
      call getfitpars(nunit,nfit,sol,serr,err)
      close(nunit) !release unit number as we are done with file
      
      do 11 i=1,npt
        dtype(i)=0
 11   continue
      
      call transitmodel(nfit,sol,npt,time,itime,tmodel,dtype)
      
      tdur=transitdur(nfit,sol)/24.0d0 !in days
c      write(0,*) "tdur: ",tdur
      
      epo=sol(1)
      per=sol(2)
      ph1=(tdur)/sol(2)/2.0d0
      ph2=1.0d0+(-tdur)/sol(2)/2.0d0      
      j=0
      do 10 i=1,npt
        ph=(time(i)-epo)/per-int((time(i)-epo)/per)
        if(ph.lt.0.0) ph=ph+1.0
        if((ph.gt.ph1).and.(ph.lt.ph2))then
            j=j+1
            pcut(j)=flux(i)-tmodel(i)
c            write(6,*) time(i),pcut(j)
        endif
 10   continue
      std=stdev(j,pcut,mean)
c      write(0,*) "std: ",j,std

      j=0
      do 12 i=1,npt
        ph=(time(i)-epo)/per-int((time(i)-epo)/per)
        if(ph.lt.0.0) ph=ph+1.0
        if((ph.le.ph1).or.(ph.ge.ph2))then
            j=j+1
            pcut(j)=flux(i)-tmodel(i)
c            write(6,*) time(i),pcut(j)
        endif
 12   continue
 
      chisq=0.0d0
      do 13 i=1,j
        chisq=chisq+pcut(i)*pcut(i)/(std*std)
 13   continue
c      write(0,*) "Chisq: ",chisq,chisq/dble(j)
      write(6,*) chisq/dble(j)
        
      
      goto 999
 901  write(0,*) 'Usage: transitchisq photometry a1.dat'
      goto 999
 902  write(0,*) 'Cannot open ',inputsol
      goto 999
 903  write(0,*) 'Cannot open ',obsfile
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function stdev(npt,pts,mean)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Calculates standard deviation of data set given the mean.
      implicit none

      integer npt,i
      double precision pts(npt),mean,s,ep,adev,p,var,sdev,mean2

      s=0.
      do 11 i=1,npt
         s=s+pts(i)
 11   continue
      mean=s/npt

      ep=0.
      var=0.
      do 10 i=1,npt
         s=pts(i)-mean
         ep=ep+s
         p=s*s
         var=var+p
 10   continue
      var=(var-ep**2/npt)/(npt-1)
      sdev=sqrt(var)

      stdev=sdev

      return
      end      
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function transitdur(nfit,sol)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit
      double precision sol(nfit),Pi,Msun,Rsun,G,aConst,b,Psec,adrs,
     .  cincl,temp(4),tpi,rdr,bb,eccn,Per

      Pi=acos(-1.d0)   !Pi
      tpi=2.0d0*pi
      Msun=1.9891d30 !kg  mass of Sun
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)
      
      eccn=0.0

      bb=sol(3)
      b=sqrt(bb)
      per=sol(2)
      Psec=per*8.64d4 !sec ; period of planet
      eccn=sqrt(sol(7)*sol(7)+sol(8)*sol(8)) !eccentricity
      adrs=sol(5)*per/tpi*sqrt((1.0+sol(4))**2.0-sol(3))
c     .  *(1+sol(8))/sqrt(1-eccn*eccn)
c      write(0,*) sol(5),per,sol(4)
c      write(0,*) sol(3)
c      write(0,*) "adrs:",adrs
      
      rdr=sol(4)
      
      cincl=b/adrs !cos(i)
        
      temp(1)=Psec/Pi
      temp(2)=1.0d0/adrs
      temp(3)=(1+rdr)**2.0-bb
      temp(4)=1-cincl*cincl
C     Transit duration in hours
      transitdur=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))/3600.0
      
      
      return
      end
