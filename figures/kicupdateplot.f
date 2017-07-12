      program KICupdateplot
      implicit none
      integer iargc,nunit,i,npar
      parameter(npar=20)
      real rbb(4),dumr(npar),px(2),py(2),G,Rsun,Msun
      character*80 filename,dumc

      G=6.67259E-8
      Rsun=6.9599E10
      Msun=1.989E33
     
      if(iargc().lt.1) goto 901
      call getarg(1,filename)
      nunit=10
      open(unit=nunit,file=filename,status='old',err=902)
      read(nunit,*) dumc
      read(nunit,*) dumc
      
      call pgopen('?')
c      call pgpage()
      call PGPAP ( 6.0 ,1.0)  
      call pgvport(0.15,0.85,0.2,0.9)
      call pgpage()
      
      rbb(1)=10000.0
      rbb(2)=3000.0
      rbb(3)=5.5
      rbb(4)=2.5
      call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))
      call pgslw(3)
      call pgsch(1.5)
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("T\deff\u (K)","log(g) (cgs)","")
      
      call pgslw(1)
      call pgsch(1.0)
 10   read(nunit,*,end=11) dumc,(dumr(i),i=1,npar)
        px(1)=dumr(17) !Teff KIC
        px(2)=dumr(7) !Teff YY
        py(1)=dumr(19) !logg KIC
c        py(2)=dumr(13) !logg YY
        py(2)=log10(G*dumr(1)*Msun/(Rsun*Rsun*dumr(3)*dumr(3)))

c        write(0,*) px(1),px(2)
c        call pgsci(2)
c        call pgpt1(px(2),py(2),17)
c        call pgsci(1)
        call pgline(2,px,py)
        
      goto 10
 11   continue
 
      rewind(nunit)
      read(nunit,*) dumc
      read(nunit,*) dumc 
      
      call pgsch(0.6)
 12   read(nunit,*,end=13) dumc,(dumr(i),i=1,npar)
        px(1)=dumr(17) !Teff KIC
        px(2)=dumr(7) !Teff YY
        py(1)=dumr(19) !logg KIC
c        py(2)=dumr(13) !logg YY
        py(2)=log10(G*dumr(1)*Msun/(Rsun*Rsun*dumr(3)*dumr(3)))
        
c        write(0,*) px(1),px(2)
        call pgsci(2)
        call pgpt1(px(2),py(2),17)
        call pgsci(1)
c        call pgline(2,px,py)
        
      goto 12
 13   continue
      
      call pgclos()
      close(nunit)
      
      goto 999
 901  write(0,*) "Usage: kicupdateplot stellarpars.txt"
      goto 999
 902  write(0,*) "Cannot open ",filename
      goto 999
 999  end