      program mcmcsetup
      implicit none
      integer iargc,nmax,nunit,npt,nfit,i,j,nplanet,np,k,nplanetmax,nf,
     .   col,ii,nptv,npta
      parameter(nmax=2000000,nfit=108,nplanetmax=10)
      integer dtype(nmax),ntt(nplanetmax)
      double precision time(nmax),flux(nmax),ferr(nmax),itime(nmax),
     .  Keplertime,sol(nfit),serr(nfit,2),err(nfit),tdur,transitdur,
     .  tmodel(nmax),pcut(nmax),epo,per,ph1,ph2,ph,std,stdev,mean,
     .  chisq,tobs(nplanetmax,nmax),omc(nplanetmax,nmax),ttcor,
     .  tcor(nmax),chisq2,sol2(nfit),dchi,serr2(nfit,2),dvar,var(18),
     .   vtime(nmax),vel(nmax),verr(nmax),vetime(nmax),aT(nmax),
     .   aM(nmax),aE(nmax),aIT(nmax),ubs(18)
      character*80 obsfile,inputsol,cline,ttfile,rvfile
      character*3 titles(18)
      data var/0.003,0.01,0.01,0.01,0.01,0.01,0.1,1.0e-6,1.0e-6,5.0e-6,
     .   0.003,1.0e-4,0.003,0.003,0.3,1.0,1.0,0.5/
      data ubs/0.003,0.01,0.01,0.01,0.01,0.01,1.0,1.0e-6,1.0e-5,1.0e-5,
     .   0.10,0.003,0.01,0.01,1.0,1.0,1.0,1.0/

      if(iargc().lt.4) goto 901

C     Parse the name of the observations data file from the commandline
      call getarg(1,obsfile)
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,flux,ferr,itime,Keplertime)
      call readkeplc(nunit,nmax,npt,time,flux,ferr,itime,Keplertime)
      close(nunit)!release unit number as we are done with the file.

      call getarg(2,rvfile) !get filename for RV data
      if(rvfile.eq.'null')then
        nptv=0
      else
        nunit=10
        open(unit=nunit,file=rvfile,status='old',err=904)
        call readrv(nunit,nmax,nptv,vtime,vel,verr,vetime)
        close(nunit)
c        write(6,*) "rv:",(vtime(i),i=1,nptv)
      endif

C     Store all the observations in the master file
      do 7 i=1,npt !central database of all data
        dtype(i)=0 !0 marks that we have photometric data
        aT(i)=time(i)
        aM(i)=flux(i)
        aE(i)=ferr(i)
        aIT(i)=itime(i)
c        write(0,*) "itime",itime(i)
 7    continue
      npta=npt
      do 8 i=1,nptv
        dtype(npt+i)=1 !mark data as RV
        aT(npt+i)=vtime(i)
        aM(npt+i)=vel(i)
        aE(npt+i)=verr(i)
        aIT(npt+i)=vetime(i)
c        write(0,*) aT(npt+i),aM(npt+i)
 8    continue
      npta=npta+nptv
c      write(0,*) npta,npt

      call getarg(3,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
c      write(0,*) "reading in input solution"
      call getfitpars(nunit,nfit,nplanet,sol,serr,err)
c      write(0,*) "done reading input solution"
c      call getfitpars(nunit,nfit,sol,serr,err)
      close(nunit) !release unit number as we are done with file

      call getarg(4,cline) !planet number for chi-sq calc..
      read(cline,*) np
      if(np.lt.0) then
         write(0,*) "np must be greater than zero"
         stop
      endif

      do 22 i=1,nplanet !deal with transit-timing variations
        if(iargc().ge.4+i)then
            call getarg(4+i,ttfile)

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

      call transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npta,aT,aIT,
     .  ntt,tobs,omc,tmodel,dtype)

C     calculate global chi-sq
      chisq=0.0d0
      do 10 i=1,npta
         chisq=chisq+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/
     .      (aE(i)*aE(i))
 10   continue

      write(0,*) "chisq: ",chisq



      nf=8+nplanet*10
      do 21 i=1,nf
         serr2(i,1)=serr(i,1)
         serr2(i,2)=serr(i,2)
 21   continue
      serr2(1,2)=0.003

      do 12 i=2,5
         dvar=var(i)
         if(serr(i,2).ne.0.0)then
            do 13 j=1,nf
               sol2(j)=sol(j)
 13         continue
            sol2(i)=sol2(i)+dvar
            call transitmodel(nfit,nplanet,nplanetmax,sol2,nmax,npta,
     .       aT,aIT,ntt,tobs,omc,tmodel,dtype)
            chisq2=0.0d0
            do 14 j=1,npta
               chisq2=chisq2+(aM(j)-tmodel(j))*(aM(j)-tmodel(j))/
     .          (aE(j)*aE(j))
 14         continue
            dchi=chisq2-chisq
            serr2(i,2)=abs(dvar/dchi)
            serr2(i,2)=min(serr2(i,2),ubs(i))
            write(0,*) i,serr2(i,2)
         endif
 12   continue

      i=7
      dvar=var(i)
      if(serr(i,2).ne.0.0)then
         do 15 j=1,nf
            sol2(j)=sol(j)
 15      continue
         sol2(i)=sol2(i)+dvar
         call transitmodel(nfit,nplanet,nplanetmax,sol2,nmax,npta,
     .       aT,aIT,ntt,tobs,omc,tmodel,dtype)
         chisq2=0.0d0
         do 16 j=1,npta
            chisq2=chisq2+(aM(j)-tmodel(j))*(aM(j)-tmodel(j))/
     .          (aE(j)*aE(j))
 16      continue
         dchi=chisq2-chisq
         serr2(i,2)=abs(dvar/dchi)
         serr2(i,2)=min(serr2(i,2),ubs(i))
         write(0,*) i,serr2(i,2)
      endif

      i=8
      dvar=var(i)
      if(serr(i,2).ne.0.0)then
         do 17 j=1,nf
            sol2(j)=sol(j)
 17      continue
         sol2(i)=sol2(i)+dvar
         call transitmodel(nfit,nplanet,nplanetmax,sol2,nmax,npta,
     .       aT,aIT,ntt,tobs,omc,tmodel,dtype)
         chisq2=0.0d0
         do 18 j=1,npta
            chisq2=chisq2+(aM(j)-tmodel(j))*(aM(j)-tmodel(j))/
     .          (aE(j)*aE(j))
 18      continue
         dchi=chisq2-chisq
         serr2(i,2)=abs(dvar/dchi)
         serr2(i,2)=min(serr2(i,2),ubs(i))
         write(0,*) i,serr2(i,2)
      endif

      do 19 k=1,nplanet
         col=10*(k-1)
         do 20 ii=9,18
            i=ii+col
            dvar=var(ii)
            if((ii.eq.11).or.(ii.eq.13).or.(ii.eq.14))then
               if(serr(i,2).ne.0)serr2(i,2)=0.003
            else

               if(serr(i,2).ne.0.0)then
                  do 23 j=1,nf
                     sol2(j)=sol(j)
 23               continue
                  sol2(i)=sol2(i)+dvar
                  call transitmodel(nfit,nplanet,nplanetmax,sol2,nmax,
     .             npta,aT,aIT,ntt,tobs,omc,tmodel,dtype)
                  chisq2=0.0d0
                  do 24 j=1,npta
                     chisq2=chisq2+(aM(j)-tmodel(j))*(aM(j)-tmodel(j))/
     .                (aE(j)*aE(j))
 24               continue
                  dchi=chisq2-chisq
                  serr2(i,2)=abs(dvar/dchi)
                  serr2(i,2)=min(serr2(i,2),ubs(ii))
                  write(0,*) i,serr2(i,2)
               endif
            endif

 20      continue
 19   continue

      if(np.gt.0)then !only export one planet
         col=10*(np-1)
         do 25 i=1,8
            sol2(i)=sol(i)
 25      continue
         do 26 i=9,18
            sol2(i)=sol(i+col)
            serr2(i,1)=serr2(i+col,1)
            serr2(i,2)=serr2(i+col,2)
 26      continue
         call exportfit(nfit,1,sol2,serr2,err,titles)
      else  !export all planets
         call exportfit(nfit,nplanet,sol,serr2,err,titles)
      endif

      goto 999
 901  write(0,*) 'Usage: transitchisq photometry rv n1.dat nplanet'
      goto 999
 902  write(0,*) 'Cannot open ',inputsol
      goto 999
 903  write(0,*) 'Cannot open ',obsfile
      goto 999
 904  write(0,*) "Error opening ",rvfile
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
