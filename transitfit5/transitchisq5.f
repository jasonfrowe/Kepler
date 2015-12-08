      program transitchisq
      implicit none
      integer iargc,nmax,nunit,npt,nfit,i,j,nplanet,np,k,nplanetmax
      parameter(nmax=2000000,nfit=108,nplanetmax=10)
      integer dtype(nmax),ntt(nplanetmax)
      double precision time(nmax),flux(nmax),ferr(nmax),itime(nmax),
     .  Keplertime,sol(nfit),serr(nfit,2),err(nfit),tdur,transitdur,
     .  tmodel(nmax),pcut(nmax),epo,per,ph1,ph2,ph,std,stdev,mean,
     .  chisq,tobs(nplanetmax,nmax),omc(nplanetmax,nmax),ttcor,
     .  tcor(nmax)
      character*80 obsfile,inputsol,cline,ttfile
      
      if(iargc().lt.3) goto 901

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
c      write(0,*) "reading in input solution"
      call getfitpars(nunit,nfit,nplanet,sol,serr,err)
c      write(0,*) "done reading input solution"
c      call getfitpars(nunit,nfit,sol,serr,err)
      close(nunit) !release unit number as we are done with file
      
      call getarg(3,cline) !planet number for chi-sq calc..
      read(cline,*) np
      
      do 22 i=1,nplanet
        if(iargc().ge.3+i)then
            call getarg(3+i,ttfile)

            if(ttfile.eq.'null')then
                ntt(i)=0
            else
                nunit=10
                open(unit=nunit,file=ttfile,status='old',err=905)
                call readttfile(nunit,nplanetmax,nmax,i,ntt,tobs,omc)
                close(nunit)
            endif
            
        else
            ntt(i)=0
        endif
c        write(0,*) "ntt",i,ntt(i)
 22   continue
      
      do 11 i=1,npt
        dtype(i)=0
 11   continue
      
      call transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npt,time,itime,
     .  ntt,tobs,omc,tmodel,dtype)
      
      tdur=transitdur(nfit,sol,np)/24.0d0 !in days
c      write(0,*) "tdur: ",tdur
      
      do 23 i=1,npt
        call lininterp(tobs,omc,nplanetmax,nmax,np,ntt,time(i),ttcor)
        tcor(i)=time(i)-ttcor
 23   continue
      
      epo=sol(8+(np-1)*10+1)      
      per=sol(8+(np-1)*10+2)
      ph1=(tdur)/per/2.0d0
      ph2=1.0d0+(-tdur)/per/2.0d0      
c      write(0,*) ph1,ph2,per

      j=0
      do 10 i=1,npt
        ph=(tcor(i)-epo)/per-int((tcor(i)-epo)/per)
        if(ph.lt.0.0) ph=ph+1.0
        if((ph.gt.ph1).and.(ph.lt.ph2))then
            j=j+1
            pcut(j)=flux(i)-tmodel(i)
c            write(6,*) time(i),pcut(j)
        endif
 10   continue
      std=stdev(j,pcut,mean)
      k=0
      do 14 i=1,j
        if(abs(pcut(i)).lt.3.0*std)then
            k=k+1
            pcut(k)=pcut(i)
        endif
 14   continue
      std=stdev(k,pcut,mean)  !standard deviation of out-of-transit
c      write(0,*) "std: ",k,std
c
      j=0
      do 12 i=1,npt
        ph=(tcor(i)-epo)/per-int((tcor(i)-epo)/per)
        if(ph.lt.0.0) ph=ph+1.0
        if((ph.le.ph1).or.(ph.ge.ph2))then
            j=j+1
            pcut(j)=flux(i)-tmodel(i)
c            write(6,*) tcor(i),flux(i)
        endif
 12   continue
 
      chisq=0.0d0
      do 13 i=1,j
        chisq=chisq+pcut(i)*pcut(i)/(std*std)
 13   continue
c      write(0,*) "Chisq: ",chisq,chisq/dble(j)
      write(6,*) chisq/dble(j)
        
      
      goto 999
 901  write(0,*) 'Usage: transitchisq photometry n1.dat nplanet'
      goto 999
 902  write(0,*) 'Cannot open ',inputsol
      goto 999
 903  write(0,*) 'Cannot open ',obsfile
      goto 999
 905  write(0,*) "Error opening ",ttfile
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
      double precision function transitdur(nfit,sol,np)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,col,np
      double precision sol(nfit),Pi,Msun,Rsun,G,aConst,b,Psec,adrs,
     .  cincl,temp(4),tpi,rdr,bb

      Pi=acos(-1.d0)   !Pi
      tpi=2.0d0*pi
      Msun=1.9891d30 !kg  mass of Sun
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)

      col=8+(np-1)*10
      bb=sol(col+3)
      b=sqrt(bb)
      Psec=sol(col+2)*8.64d4 !sec ; period of planet
      adrs=1000.0*sol(1)*G*(Psec)**2/(3.0d0*Pi)
      adrs=adrs**(1.0d0/3.0d0) !a/R*
      cincl=b/adrs !cos(i)
      rdr=sol(col+4)
        
      temp(1)=Psec/Pi
      temp(2)=1.0d0/adrs
      temp(3)=(1+rdr)**2.0-bb
      temp(4)=1-cincl*cincl
C     Transit duration in hours
      transitdur=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))/3600.0
      
      
      return
      end
