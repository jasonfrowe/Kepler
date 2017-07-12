      program transitsn5
C     determine S/N of transit.
      implicit none
      integer iargc,nunit,nfit,nmax,nplanetmax,npt,nplanet,dplanet,i
      parameter(nfit=108,nmax=2000000,nplanetmax=10)
      integer ntt(nplanetmax),dtype(nmax)
      double precision time(nmax),flux(nmax),ferr(nmax),itime(nmax),
     .   Keplertime,sol(nfit),serr(nfit,2),err(nfit),
     .   tobs(nplanetmax,nmax),omc(nplanetmax,nmax),tmodel(nmax),esw,
     .   ecw,eccn,w,Eanom,Manom,phi0,epoch,phi,Pi,tPi,pid2,ttcor,t,per,
     .   rprs,tmodel2(nmax),res(nmax),stdev,std,mean,SNR
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


C     Calculate standard deviation
C     Get transit model
      call transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npt,time,
     .   itime,ntt,tobs,omc,tmodel,dtype)
C     Get residuals
      do 10 i=1,npt
         res(i)=flux(i)-tmodel(i)
 10   continue
C     get standard devation
      std=stdev(npt,res,mean)
c      write(0,*) "Std: ",std


c      rprs=sol(10*(dplanet-1)+8+4)
c      sol(10*(dplanet-1)+8+4)=0.0d0 !set r/R*=0
c      call transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npt,time,
c     .   itime,ntt,tobs,omc,tmodel,dtype)
c      sol(10*(dplanet-1)+8+4)=rprs
      do 11 i=1,nplanet
         if(i.ne.dplanet) sol(10*(i-1)+8+4)=0.0d0
 11   continue
      call transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npt,time,
     .   itime,ntt,tobs,omc,tmodel,dtype)


      SNR=0.d0
      do 12 i=1,npt
!         write(0,*) tmodel(i),sol(8),tmodel(i)-sol(8)-1.0
         SNR=SNR+((tmodel(i)-sol(8)-1.0)/std)**2.0d0
 12   continue
      SNR=sqrt(SNR)
      write(6,500) SNR,STD

 500    format(F17.11,4(1X,F17.11))

      goto 999 !clean exit
 901  write(0,*) "Usage: transitfit5 <photfile> <fitpars> <planetnumber>
     . [ttfiles]"
      stop 11
 902  write(0,*) "Error opening ",inputsol
      stop 12
 903  write(0,*) "Error opening ",obsfile
      stop 13
 904  write(0,*) "Error reading planetnumber: ",cline
      stop 14
 905  write(0,*) "Error planetnumber out of range: ",dplanet, nplanet
      stop 15
 906  write(0,*) "Error opening ",ttfile
      stop 16
 999  end

