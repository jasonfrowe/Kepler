      program datadump5
C     Produces photometry for BLENDER analysis.
      implicit none
      integer iargc,nunit,nfit,nmax,nplanetmax,npt,nplanet,dplanet,i
      parameter(nfit=108,nmax=2000000,nplanetmax=10)
      integer ntt(nplanetmax),dtype(nmax)
      double precision time(nmax),flux(nmax),ferr(nmax),itime(nmax),
     .   Keplertime,sol(nfit),serr(nfit,2),err(nfit),
     .   tobs(nplanetmax,nmax),omc(nplanetmax,nmax),tmodel(nmax),esw,
     .   ecw,eccn,w,Eanom,Manom,phi0,epoch,phi,Pi,tPi,pid2,ttcor,t,per,
     .   rprs,tmodel2(nmax)
      character*80 obsfile,inputsol,cline,ttfile

      Pi=acos(-1.d0)!define Pi and 2*Pi
      tPi=2.0d0*Pi
      pid2=Pi/2.0d0

C     arguments are photometry, n1.dat and planet number
      if(iargc().lt.3) goto 901

C     Parse the name of the observations data file from the commandline
      call getarg(1,obsfile)
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
      call readkeplc(nunit,nmax,npt,time,flux,ferr,itime,Keplertime)
      close(nunit)!release unit number as we are done with the file.

      do 17 i=1,npt !central database of all data
        dtype(i)=0 !0 marks that we have photometric data
 17   continue

      call getarg(2,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
      write(0,*) "reading in input solution"
      call getfitpars(nunit,nfit,nplanet,sol,serr,err)
      write(0,*) "done reading input solution"
      close(nunit) !release unit number as we are done with file

      call getarg(3,cline)
      read(cline,*,err=904) dplanet
      if((dplanet.le.0).or.(dplanet.gt.nplanet)) goto 905

      do 21 i=1,nplanet
        if(iargc().ge.3+i)then
            call getarg(3+i,ttfile)

            if(ttfile.eq.'null')then
                ntt(i)=0
            else
                nunit=10
                open(unit=nunit,file=ttfile,status='old',err=906)
                call readttfile(nunit,nplanetmax,nmax,i,ntt,tobs,omc)
                close(nunit)
            endif

        else
            ntt(i)=0
        endif
c        write(0,*) "ntt",i,ntt(i)
 21   continue

      rprs=sol(10*(dplanet-1)+8+4)
      sol(10*(dplanet-1)+8+4)=0.0d0 !set r/R*=0
      sol(10*(dplanet-1)+8+7)=0.0d0 !set K=0
      sol(10*(dplanet-1)+8+8)=0.0d0 !set TED=0
      sol(10*(dplanet-1)+8+9)=0.0d0 !set ELL=0
      sol(10*(dplanet-1)+8+10)=0.0d0 !set ALB=0
      call transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npt,time,
     .   itime,ntt,tobs,omc,tmodel,dtype)
!      sol(10*(dplanet-1)+8+4)=rprs
!      do 11 i=1,nplanet
!         if(i.ne.dplanet) sol(10*(i-1)+8+4)=0.0d0
! 11   continue
!      call transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npt,time,
!     .   itime,ntt,tobs,omc,tmodel2,dtype)

!      per=sol(10*(dplanet-1)+8+2)
!      ecw=sol(10*(dplanet-1)+8+6)
!      esw=sol(10*(dplanet-1)+8+5)
!      eccn=sqrt(ecw*ecw+esw*esw) !eccentricity
!      if(eccn.ge.1.0) eccn=0.99
!      if(eccn.eq.0.0d0)then
!         w=0.0d0
!      else
!         if(ecw.eq.0.0d0)then
!            w=Pi/2.0d0
!         else
!            w=atan(esw/ecw)
!         endif
!         if((ecw.gt.0.0d0).and.(esw.lt.0.0d0))then
!            w=tPi+w
!         elseif((ecw.lt.0.0d0).and.(esw.ge.0.0d0))then
!            w=Pi+w
!         elseif((ecw.le.0.0d0).and.(esw.lt.0.0d0))then
!            w=Pi+w
!         endif
!      endif
!      epoch=sol(10*(dplanet-1)+8+1)
!      Eanom=tan(w/2.0d0)/sqrt((1.0d0+eccn)/(1.0d0-eccn)) !mean anomaly
!      Eanom=2.0d0*atan(Eanom)
!      phi0=Eanom-eccn*sin(Eanom)

      do 10 i=1,npt
         call lininterp(tobs,omc,nplanetmax,nmax,dplanet,ntt,time(i),
     .          ttcor)
         t=time(i)-epoch-ttcor

         phi=t/per-int(t/per)
         if(phi.lt.0.0)phi=phi+1.0
         Manom=phi*tPi+phi0

         write(6,500) time(i)+54900.0d0-0.5d0,flux(i)-tmodel(i),
     .       ferr(i),itime(i)

c         if(itime(i).gt.0.001)then
c            write(6,501) time(i)+54900.0d0-0.5d0,flux(i)-tmodel(i),
c     .       ferr(i),0
c         else
c            write(6,501) time(i)+54900.0d0-0.5d0,flux(i)-tmodel(i),
c     .       ferr(i),1
c         endif

c         write(6,*) time(i),flux(i)-tmodel(i),tmodel2(i)-1.0
c         write(6,500) time(i)+54900.0d0-0.5d0,phi,
c     .      1.0d0+flux(i)-tmodel(i),tmodel2(i),ferr(i)
 500    format(F17.11,4(1X,F17.11))
 501    format(F17.11,1X,F17.11,1X,F17.11,1X,I1)
 10   continue

      goto 999
 901  write(0,*) "Usage: transitfit5 <photfile> <fitpars> <planetnumber>
     . [ttfiles]"
      goto 999
 902  write(0,*) "Error opening ",inputsol
      goto 999
 903  write(0,*) "Error opening ",obsfile
      goto 999
 904  write(0,*) "Error reading planetnumber: ",cline
      goto 999
 905  write(0,*) "Error planetnumber out of range: ",dplanet, nplanet
      goto 999
 906  write(0,*) "Error opening ",ttfile
      goto 999
 999  end

