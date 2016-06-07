      program transitdepth5
      implicit none
      integer iargc,nmax,nunit,npt,nfit,i,j,nplanet,np,k,nplanetmax,col
      parameter(nmax=2000000,nfit=108,nplanetmax=10)
      integer dtype(1),ntt(nplanetmax)
      double precision time(1),flux(nmax),ferr(nmax),itime(1),
     .  Keplertime,sol(nfit),serr(nfit,2),err(nfit),tdur,transitdur,
     .  tmodel(1),pcut(nmax),epo,per,ph1,ph2,ph,std,stdev,mean,
     .  chisq,tobs(nplanetmax,nmax),omc(nplanetmax,nmax),ttcor,
     .  tcor(nmax),sol2(nfit),tdep
      character*80 obsfile,inputsol,cline,ttfile

      if(iargc().lt.2) goto 901

      call getarg(1,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
c      write(0,*) "reading in input solution"
      call getfitpars(nunit,nfit,nplanet,sol,serr,err)
c      write(0,*) "done reading input solution"
c      call getfitpars(nunit,nfit,sol,serr,err)
      close(nunit) !release unit number as we are done with file

      call getarg(2,cline) !planet number for chi-sq calc..
      read(cline,*) np

      do 22 i=1,nplanet
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
 22   continue

      do 11 i=1,npt
        dtype(i)=0
 11   continue

      col=10*(np-1)
      do 10 i=1,8
         sol2(i)=sol(i)
 10   continue
      do 12 i=9,18
         sol2(i)=sol(i+col)
 12   continue

c      do 13 i=1,18
c         write(0,*) sol2(i)
c 13   continue

      dtype(1)=0
      time(1)=sol2(9)

c      call transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npt,time,itime,
c     .  ntt,tobs,omc,tmodel,dtype)

      call transitmodel(nfit,1,nplanetmax,sol2,nmax,1,time,
     .  itime,ntt,tobs,omc,tmodel,dtype)

      tdep=(1.0d0-tmodel(1)+sol(8))*1.0d6

      write(6,500) tdep
! 500  format(F8.1)
 500  format(F13.5)

      goto 999
 901  write(0,*) 'Usage: transitdepth5 n1.dat nplanet'
      goto 999
 902  write(0,*) 'Cannot open ',inputsol
      goto 999
 903  write(0,*) 'Cannot open ',obsfile
      goto 999
 905  write(0,*) "Error opening ",ttfile
      goto 999
 999  end
