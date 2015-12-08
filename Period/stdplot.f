      program stdplot
      implicit none
      integer nmax,nunit,npt
      parameter(nmax=52496)
      integer kID(nmax),nfreq(nmax)
      real kmag(nmax),Teff(nmax),logg(nmax),rad(nmax),t(nmax),
     .  std(2,nmax),px(nmax),py(nmax)
      character*80 filename
     
      if(iargc().lt.1) goto 903
      call getarg(1,filename)
      
      nunit=10
      call readdata(nunit,filename,nmax,npt,kID,Kmag,Teff,logg,rad,t,
     .  std,nfreq)
     
      call pgopen('?')
      call pgask(.true.)
      call PGPAP ( 6.0 ,1.0) !square plot
      call pgpage()
     
      call plotTefflogg(npt,Teff,logg,t)
      call pgpage()
      call plotTeffstd1(nmax,npt,Teff,std,py,t)
      call pgpage()
      call plotloggstd1(nmax,npt,logg,std,py,t)
      call pgpage()
      call plotKmagstd1(nmax,npt,Kmag,std,py,t)
 
      call pgclos()
 
      goto 999
 903  write(0,*) "Usage stdplot <std.dat>"
      write(0,*) " "
      write(0,*) "Where <std.dat> is the output from keplerper"
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotKmagstd1(nmax,npt,Kmag,std,py,t)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i,nmax
      real Kmag(npt),std(2,npt),bb(4),py(nmax),t(nmax)
      
      bb(1)= 99.9e30
      bb(2)=-99.9e30
      bb(3)= 99.9e30
      bb(4)=-99.9e30
      do 10 i=1,npt
        if(Kmag(i).gt.0.0)then
            bb(1)=min(Kmag(i),bb(1))
            bb(2)=max(Kmag(i),bb(2))
        endif
        if(std(1,i).gt.0.0)then
            bb(3)=min(std(1,i),bb(3))
            bb(4)=max(std(1,i),bb(4))
        endif
        py(i)=std(1,i)
 10   continue
 
      call pgwindow(bb(1),bb(2),bb(3),1000.0)
      call pgbbuf()
      do 11 i=1,npt
        if(t(i).lt.-0.5)then
            call pgsci(3)
            call pgpt1(Kmag(i),py(i),17)
        elseif(t(i).gt.0.5)then
            call pgsci(2)
            call pgpt1(Kmag(i),py(i),17)
        else
            call pgsci(1)
            call pgpt1(Kmag(i),py(i),1)
        endif
 11   continue 
      call pgebuf()
      call pgsci(1)
      call pgslw(2) !thicker lines
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Kepler magnitude","Standard Deviation (ppm)","")
 
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotloggstd1(nmax,npt,logg,std,py,t)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i,nmax
      real logg(npt),std(2,npt),bb(4),py(nmax),t(nmax)
      
      bb(1)= 99.9e30
      bb(2)=-99.9e30
      bb(3)= 99.9e30
      bb(4)=-99.9e30
      do 10 i=1,npt
        if(logg(i).gt.0.0)then
            bb(1)=min(logg(i),bb(1))
            bb(2)=max(logg(i),bb(2))
        endif
        if(std(1,i).gt.0.0)then
            bb(3)=min(std(1,i),bb(3))
            bb(4)=max(std(1,i),bb(4))
        endif
        py(i)=std(1,i)
 10   continue
 
      call pgwindow(bb(1),bb(2),bb(3),1000.0)
            call pgbbuf()
      do 11 i=1,npt
        if(t(i).lt.-0.5)then
            call pgsci(3)
            call pgpt1(logg(i),py(i),17)
        elseif(t(i).gt.0.5)then
            call pgsci(2)
            call pgpt1(logg(i),py(i),17)
        else
            call pgsci(1)
            call pgpt1(logg(i),py(i),1)
        endif
 11   continue 
      call pgebuf()
      call pgsci(1)
      call pgslw(2) !thicker lines
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("logg (cgs)","Standard Deviation (ppm)","")
 
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotTeffstd1(nmax,npt,Teff,std,py,t)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i,nmax
      real Teff(npt),std(2,npt),bb(4),py(nmax),t(npt)
      
      bb(1)= 99.9e30
      bb(2)=-99.9e30
      bb(3)= 99.9e30
      bb(4)=-99.9e30
      do 10 i=1,npt
        if(Teff(i).gt.0.0)then
            bb(1)=min(Teff(i),bb(1))
            bb(2)=max(Teff(i),bb(2))
        endif
        if(std(1,i).gt.0.0)then
            bb(3)=min(std(1,i),bb(3))
            bb(4)=max(std(1,i),bb(4))
        endif
        py(i)=std(1,i)
 10   continue
 
      call pgwindow(bb(1),bb(2),bb(3),1000.0)
      call pgbbuf()
      do 11 i=1,npt
        if(t(i).lt.-0.5)then
            call pgsci(3)
            call pgpt1(teff(i),py(i),17)
        elseif(t(i).gt.0.5)then
            call pgsci(2)
            call pgpt1(teff(i),py(i),17)
        else
            call pgsci(1)
            call pgpt1(teff(i),py(i),1)
        endif
 11   continue 
      call pgebuf()
      call pgsci(1)
      call pgslw(2) !thicker lines
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Teff (K)","Standard Deviation (ppm)","")
 
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotTefflogg(npt,Teff,logg,t)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i
      real Teff(npt),logg(npt),bb(4),t(npt)
      
      bb(1)= 99.9e30
      bb(2)=-99.9e30
      bb(3)= 99.9e30
      bb(4)=-99.9e30
      do 10 i=1,npt
        if(Teff(i).gt.0.0)then
            bb(1)=min(Teff(i),bb(1))
            bb(2)=max(Teff(i),bb(2))
        endif
        if(logg(i).gt.0.0)then
            bb(3)=min(logg(i),bb(3))
            bb(4)=max(logg(i),bb(4))
        endif
 10   continue
 
      call pgwindow(bb(1),bb(2),bb(3),bb(4))
      call pgbbuf()
      do 11 i=1,npt
        if(t(i).lt.-0.5)then
            call pgsci(3)
            call pgpt1(teff(i),logg(i),17)
        elseif(t(i).gt.0.5)then
            call pgsci(2)
            call pgpt1(teff(i),logg(i),17)
        else
            call pgsci(1)
            call pgpt1(teff(i),logg(i),1)
        endif
 11   continue
      call pgsci(1)
      call pgebuf()
      call pgslw(2) !thicker lines
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Teff (K)","logg (cgs)","")
 
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readdata(nunit,filename,nmax,npt,kID,Kmag,Teff,logg,
     .  rad,t,std,nfreq)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nmax,npt,nfreq(nmax),kID(nmax),i
      real Kmag(nmax),Teff(nmax),logg(nmax),rad(nmax),t(nmax),
     .  std(2,nmax)
      character*80 filename
     
      open(unit=nunit,file=filename,status='old',err=901)
      
      i=1
  9   continue
 10   read(nunit,*,end=11,err=9) kID(i),Kmag(i),Teff(i),logg(i),rad(i),
     .  t(i),std(1,i),std(2,i),nfreq(i)
        i=i+1
        goto 10
 11   continue
      npt=i-1 
      
      write(0,*) "Number of points read: ",npt
      
      close(nunit)
      goto 999
 901  write(0,*) "Cannot open ",filename
      stop
      goto 999
 999  return
      end