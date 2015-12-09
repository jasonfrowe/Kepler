      program prdump
C     Updated PRDUMP
      implicit none
      integer today(3),nunit,nmax,npt,nfit,nplanet,nplanetmax,i,nsample
      parameter(nfit=108,nmax=2000000,nplanetmax=10)
      double precision Pi,tPi,time(nmax),flux(nmax),ferr(nmax),
     .   itime(nmax),Keplertime,sol(nfit),serr(nfit,2),err(nfit),
     .   transitdur,tdur,m1,m2,w,dtype(nmax),tmodel(nmax),per,epoch,
     .   meananom,dnintg,dnintgm1,stime
      character*80 obsfile,inputsol,ttfile
C     TT variations
      integer ntt(nplanetmax)
      double precision tobs(nplanetmax,nmax),omc(nplanetmax,nmax)


      Pi=acos(-1.d0)   !Pi
      tPi=2.0d0*Pi

      call idate(today)

      if(iargc().lt.2) goto 901 !check number of command line arguements
      call getarg(1,obsfile) !get filename for input solution
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
      call readkeplc(nunit,nmax,npt,time,flux,ferr,itime,Keplertime)
      do 10 i=1,npt !central database of all data
        dtype(i)=0 !0 marks that we have photometric data
 10   continue
      close(nunit)!release unit number as we are done with the file.

      call getarg(2,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
      call getfitpars(nunit,nfit,nplanet,sol,serr,err)
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


      tdur=transitdur(1,nfit,sol)/86400.0d0
      write(0,*) "Tdur: ",tdur

      w=0.0d0!  assuming circular orbits

      m1=(-tdur)/sol(10)-floor((-tdur)/sol(10))
      m1=m1*tPi+w
      if(m1.gt.tPi) m1=m1-tPi
      if(m1.lt.0.0d0) m1=m1+tPi

      m2=(tdur)/sol(10)-floor((tdur)/sol(10))
      m2=m2*tPi+w
      if(m2.gt.tPi) m2=m2-tPi
      if(m2.lt.0.0d0) m2=m2+tPi

      write(0,*) "m1,m2:",m1,m2

      call transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npt,time,itime,
     .  ntt,tobs,omc,tmodel,dtype)

c      write(6,500) "# Kepler Observations, N=",npt
 500  format(A25,I6)
      write(6,501) today(2),today(1),2000+today(3)
 501  format('<!-- data for Kepler XXXb, generated ',2(i2.2,'/'),i4.4,
     .  ' -->')
      write(6,504) '<dataPoints>'
 504  format(A12)

      epoch=sol(9)
      per=sol(10)
      do 11 i=1,npt
         meananom=((time(i)-epoch)/per-floor((time(i)-epoch)/per))*tPi
         if(meananom.gt.tPi) meananom=meananom-tPi
         if(meananom.lt.0.0d0) meananom=meananom+tPi
         if((meananom.gt.m1).or.(meananom.lt.m2)) then
            write(6,502) meananom,flux(i)
C            write(6,*) meananom,flux(i)
 502        format('<pt><ma>',F10.8,'</ma><i>',F10.8,'</i></pt>')
         endif
 11   continue
      write(6,505) '</dataPoints>'
 505  format(A13)

      nsample=200
      dnintg=dble(nsample)
      dnintgm1=2.0*dnintg-2.0
      stime=tdur*2.0d0
      write(0,*) "stime:",stime,tdur
      do 12 i=1,nsample
        time(i)=epoch+stime*(2.0*dble(i)-dnintg-1.0)/dnintgm1
c        write(0,*) time(i),stime*(2.0*dble(i)-dnintg-1.0)/dnintgm1
 12   continue

      npt=nsample
      call transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npt,time,itime,
     .  ntt,tobs,omc,tmodel,dtype)

c      write(6,500) "# Model Observations,  N=",nsample
      write(6,506) '<curvePoints>'
 506  format(A13)
      do 13 i=1,npt
        meananom=((time(i)-epoch)/per-floor((time(i)-epoch)/per))*tPi
        if((meananom.gt.m1).or.(meananom.lt.m2))
     .      write(6,503) meananom,tmodel(i)
 503        format('<pt><ma>',F10.8,'</ma><i>',F10.8,'</i></pt>')
 13   continue
      write(6,507) '</curvePoints>'
 507  format(A14)

      goto 999
 901  write(0,*) "Usage: modelbuild <obsfile> <modelpars>"
      write(0,*) " <obsfile>: file containing photometry"
      write(0,*) " <modelpars>: model parameters to adjust"
      write(0,*) "              omit to use defaults"
      goto 999
 902  write(0,*) "Cannot open ",inputsol
      goto 999
 903  write(0,*) "Cannot open ",obsfile
      goto 999
 905  write(0,*) "Error opening ",ttfile
      goto 999
 999  end

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

      bb=sol(11+10*(np-1))
      b=sqrt(bb)
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
