      program plottimeseries
      implicit none
      integer nmax,nunit,nfit,i,npt,nptr,col,j,ii,nplanetmax
      parameter(nmax=2000000,nfit=108,nplanetmax=10)
      integer dtype(nmax),rdtype(nmax),nplanet,order(nplanetmax)
      real rbb(4),px(nmax),py(nmax)
      double precision sol(nfit),time(nmax),
     .  flux(nmax),ferr(nmax),itime(nmax),Keplertime,serr(nfit,2),
     .  err(nfit),rtime(nmax),rflux(nmax),rferr(nmax),
     .  ritime(nmax),tmin,tmax,mmin,mmax,T0,pos
      character*80 inputsol,obsfile,rawobs
c      data order /2,5,1,3,4,6,7,8,9,10/
      data order /1,2,3,4,5,6,7,8,9,10/
     
      pos=0.035 !offset for transit markers
     
      if(iargc().lt.3) goto 901

      call getarg(3,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
      write(0,*) "reading in input solution"
      call getfitpars(nunit,nfit,nplanet,sol,serr,err)
      write(0,*) "done reading input solution"
c      call getfitpars(nunit,nfit,sol,serr,err)
      close(nunit) !release unit number as we are done with file
      
C     Parse the name of the observations data file from the commandline
      call getarg(2,obsfile)
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,flux,ferr,itime,Keplertime)
      call readkeplc(nunit,nmax,npt,time,flux,ferr,itime,Keplertime)
      close(nunit)!release unit number as we are done with the file.
      
C     Parse the name of the observations data file from the commandline
      call getarg(1,rawobs)
      nunit=10
      open(unit=nunit,file=rawobs,status='old',err=904)
c      call readdata(nunit,nmax,nptr,rtime,rflux,rferr,ritime,Keplertime)
      call readkeplc(nunit,nmax,nptr,rtime,rflux,rferr,ritime,
     .  Keplertime)
      close(nunit)!release unit number as we are done with the file.
      
      
      mmin=rflux(1)
      mmax=rflux(1)
      tmin=rtime(1)
      tmax=rtime(1)
      do 15 i=2,nptr
        mmin=min(mmin,rflux(i))
        mmax=max(mmax,rflux(i))
        tmin=min(tmin,rtime(i))
        tmax=max(tmax,rtime(i))
 15   continue
      rbb(1)=real(tmin)
      rbb(2)=real(tmax)
      rbb(3)=real(mmin)
      rbb(4)=real(mmax)
      rbb(3)=rbb(3)-(rbb(4)-rbb(3))*pos*real(nplanet)
      
      do 17 i=1,nptr !central database of all data
        rdtype(i)=0 !0 marks that we have photometric data
 17   continue

      call pgopen('?') !open PGPlot device
      call pgask(.true.) !don't ask for new page.. just do it.
      call PGPAP ( 11.0 ,0.5) !paper size
      call pgsubp(1,3)  !break up plot into grid
      
      call pgpage()
      call pgslw(2)
      call pgsch(4.0)

c      goto 500

      call pgvport(0.15,0.95,0.0,0.8)
      call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))
      call pgptxt(rbb(1)-0.10*(rbb(2)-rbb(1)),(rbb(3)+rbb(4))/2.0,90.0,
     .      0.5,"Relative Flux")
      call pgslw(1)
      
      do 10 i=1,nptr
        px(i)=real(rtime(i))
        py(i)=real(rflux(i))
c        write(6,*) px(i),py(i)
 10   continue
      call pgpt(nptr,px,py,-1)
      call pgslw(2)
      call pgbox('BCTS1',0.0,0,'BCNTSV1',0.0,0)
      
      call pgsch(2.0)
      do 14 ii=1,nplanet
        j=order(ii)
        col=10*(j-1)
        call pgsci(j+1)
        T0=sol(9+col)-int((sol(9+col)-Tmin)/sol(10+col))*sol(10+col)
        do 16 i=1,int((Tmax-T0)/sol(10+col))+1
            px(1)=real(T0)+real(i-1)*real(sol(10+col))
            py(1)=rbb(3)+(rbb(4)-rbb(3))*pos*real(ii)
            call pgpt(1,px,py,13)
 16     continue
 14   continue
      call pgsci(1)
      call pgslw(1)
      call pgsch(4.0)

CCCCCCC PLOT CORRECTED TIME SERIES CCCCCCC

 500  continue

      mmin=flux(1)
      mmax=flux(1)
      do 11 i=2,npt
        mmin=min(mmin,flux(i))
        mmax=max(mmax,flux(i))
 11   continue
      rbb(3)=real(mmin)
      rbb(4)=real(mmax)
      rbb(3)=rbb(3)-(rbb(4)-rbb(3))*pos*real(nplanet)
      
      do 12 i=1,npt !central database of all data
        dtype(i)=0 !0 marks that we have photometric data
 12   continue
      
      call pgpage()
      call pgslw(2)
      call pgsch(4.0)
      call pgvport(0.15,0.95,0.2,1.0)
      call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))
c      call pglabel("BJD-2454900","","")
      call pgptxt((rbb(1)+rbb(2))/2.0,rbb(3)-
     .          0.30*(rbb(4)-rbb(3)),0.0,
     .          0.5,"BJD-2451545")
      call pgptxt(rbb(1)-0.10*(rbb(2)-rbb(1)),(rbb(3)+rbb(4))/2.0,90.0,
     .      0.5,"Relative Flux")
      call pgslw(1)
      
      do 13 i=1,npt
        px(i)=real(time(i))
        py(i)=real(flux(i))
c        write(6,*) px(i),py(i)
 13   continue
      call pgslw(2)
      call pgpt(npt,px,py,-1)
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)

      call pgsch(2.0)
      do 18 ii=1,nplanet
        j=order(ii)
        write(6,*) j
        col=10*(j-1)
        call pgsci(j+1)
        T0=sol(9+col)-int((sol(9+col)-Tmin)/sol(10+col))*sol(10+col)
        do 19 i=1,int((Tmax-T0)/sol(10+col))+1
            px(1)=real(T0)+real(i-1)*real(sol(10+col))
            py(1)=rbb(3)+(rbb(4)-rbb(3))*pos*real(ii)
            call pgpt(1,px,py,13)
 19     continue
 18   continue
      call pgsci(1)
      call pgslw(1)
      call pgsch(4.0)
      
  
      call pgclos()
      goto 999
 901  write(0,*) "Usage: transitplot5 <rawphot> <photfile> <fitpars>"
      goto 999
 902  write(0,*) "Error opening ",inputsol
      goto 999
 903  write(0,*) "Error opening ",obsfile
      goto 999
 904  write(0,*) "Error opening ",rawobs
      goto 999
 999  end
