      program transitcut5
C     Cut out all data, except the transits!
      implicit none
      integer iargc,nunit,nmax,npt,nfit,nplanet,nplanetmax,i
      parameter(nmax=2000000,nfit=108,nplanetmax=10)
      integer ntt(nplanetmax),tflag(nmax)
      double precision time(nmax),flux(nmax),ferr(nmax),sol(nfit),
     .   serr(nfit,2),err(nfit),tobs(nplanetmax,nmax),
     .   omc(nplanetmax,nmax),phase(nmax),itime(nmax),Keplertime
      character*80 filename,inputsol,ttfile

      if(iargc().lt.2) goto 901

      call getarg(1,filename)

      nunit=10
      open(unit=nunit,file=filename,status='old',err=902)
      call readkeplc(nunit,nmax,npt,time,flux,ferr,itime,
     .  Keplertime)
      close(nunit)
      if(npt.eq.0) goto 999
      flux=flux-1.0d0

      call getarg(2,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=904)
C     We start by reading in solution from input file
c     write(0,*) "reading in input solution"
      call getfitpars(nunit,nfit,nplanet,sol,serr,err)
c     write(0,*) "done reading input solution"
      close(nunit) !release unit number as we are done with file

      do 21 i=1,nplanet
        if(iargc().ge.2+i)then
            call getarg(2+i,ttfile)

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
 21   continue

      write(0,*) "nplanet:",nplanet

      do 22 i=1,nplanet
        call marktransit(i,npt,phase,time,tflag,nfit,sol,ntt,tobs,omc)
 22   continue


 9    do 10 i=1,npt
c         write(0,*) tflag(i)
         if(tflag(i).eq.1) write(500,*) time(i)-0.5d0+54900.0d0,flux(i),
     .      ferr(i),itime(i)
 10   continue
 500  format(4(F17.11,1X))

      goto 999
 901  write(0,*) "Usage: transitcut5 <filein> <n1.dat>"
      write(0,*) "<n1.dat> transit solution"
      goto 999
 902  write(0,*) "Cannot open ",filename
      goto 999
 904  write(6,*) "Cannot open ",inputsol
      goto 999
 905  write(0,*) "Error opening ",ttfile
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine marktransit(np,npt,phase,time,tflag,nfit,sol,ntt,tobs,
     .  omc)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nplanetmax,nmax,nfit
      parameter(nmax=2000000,nplanetmax=10)
      integer tflag(npt),i,np,ntt(nplanetmax),col
      double precision time(npt),sol(nfit),ph1,ph2,transitdur,
     .  toff,tdur,phase(npt),period,tobs(nplanetmax,nmax),
     .  omc(nplanetmax,nmax),tcor(nmax),ttcor,epo,tdurfac

      tdur=2.0d0*transitdur(np,nfit,sol)/86400.0d0+0.03
      tdur=max(tdur,2.55d0/24.0) !make sure duration is at least 2.5 hours
c      write(0,*) 'tdur: ',tdur

      col=10*(np-1)
      epo=sol(9+col)
      period=sol(10+col)

C     tdurfac - 0.5: remove exactly 1 transit duration
C     tdurfac - 1.0: remove +/- 1 transit duration centred on transit
      tdurfac=1.0d0
      ph1=0.75-tdurfac*tdur/period
      if(ph1.lt.0.5)ph1=0.5
      ph2=0.75+tdurfac*tdur/period
      if(ph2.gt.1.0)ph2=1.0
c      write(0,*) "ph1,ph2",ph1,ph2

      toff=0.75-(epo/period-int(epo/period))
      if(toff.lt.0.0)toff=toff+1.0
c      write(0,*) "Toff:",toff

      do 24 i=1,npt
         call lininterp(tobs,omc,nplanetmax,nmax,np,ntt,time(i),ttcor)
         tcor(i)=time(i)-ttcor
c         write(0,*) tcor(i),ttcor
 24   continue

      call phasept(npt,tcor,phase,period,toff)

      do 10 i=1,npt
c        write(0,*) phase(i),ph1,ph2
        if((phase(i).ge.ph1).and.(phase(i).le.ph2))then
            tflag(i)=1
c            write(0,*) time(i)
c            read(5,*)
        endif
 10   continue


      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function transitdur(np,nfit,sol)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer np,nfit
      double precision sol(nfit),b,Psec,G,aConst,Pi,adrs,cincl,temp(4),
     .  bb,rdr

      Pi=acos(-1.d0)   !Pi
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)

      b=sol(11+10*(np-1))
      bb=b*b
      Psec=sol(10+10*(np-1))*8.64d4 !sec ; period of planet
      adrs=1000.0*sol(1)*G*(Psec)**2/(3.0d0*Pi)
      adrs=adrs**(1.0d0/3.0d0) !a/R*
      cincl=b/adrs !cos(i)
      rdr=sol(12+10*(np-1))
c      write(0,*) bb,adrs,rdr

      temp(1)=Psec/Pi
      temp(2)=1.0d0/adrs
      temp(3)=(1+rdr)**2.0-bb
      temp(4)=1-cincl*cincl
C     Transit duration in days
      transitdur=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine phasept(npt,time,phase,period,toff)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     npt - number of data points
C     time - times of observations
C     mag - the observations
C     phase - returned phase of observation
C     period - fixed period for data
      implicit none
      integer npt
      double precision time(npt),phase(npt),period,toff

      integer i
      double precision temp

      do 10 i=1,npt
         temp=time(i)
C        Get the phase
         phase(i)=temp/period-int(temp/period)
C        apply optional phase offset to make plot pretty
         phase(i)=phase(i)+toff
C        make sure phase is between 0 and 1
         if(phase(i).lt.0.0) phase(i)=phase(i)+1.0
         if(phase(i).gt.1.0) phase(i)=phase(i)-1.0
 10   continue
      return
      end
