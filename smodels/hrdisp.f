      program hrdisp
      implicit none
      integer nmax,nplot,npt,iargc,nline,i
      parameter(nmax=10000)
C     nmax defines the maximum number of points in a file that can be 
C     plotted.
      real bounds(4),age(nmax),Teff(nmax),Lum(nmax),Rad(nmax),ageb(2),
     .	temp(4),sf,rconv(nmax),rconvb(2),Radb(2),Mass(nmax),logg(nmax),
     .  rho(nmax)
      character*80 filelist,filename
C     filelist is an ascii file listing the HR datafiles for plotting      
      
C     We read in the filelist from the command line.  If no commandline
C     arguement is given, then we exit with an error.
      if(iargc().lt.1) goto 901
      call getarg(1,filelist)

C     Open up the ascii file listing the HR files      
      open(unit=10,file=filelist,status='old',err=902)
  
      nline=1 !nline defines the current line in filelist    
C     Note that filename is currently constrained to 80 characters.
C     Refine the variable filename above to expand.
 10   read(10,500,end=11) filename !the end statement defines our loop
 500  format(A80)
C       open up the HR file for reading
      	open(unit=11,file=filename,status='old',err=903)
      	read(filename(9:12),*) mass(nline)
C       retrieve data from the file
      	call getdata(nmax,npt,age,Teff,Lum,Rad,rconv,mass(nline),logg,
     .      rho)
C       initialize ploting bounds to first data point we read in.
      	if(nline.eq.1)then
      		bounds(1)=Teff(1)
      		bounds(2)=Teff(1)
      		bounds(3)=rho(1)!logg(1)!Lum(1)
      		bounds(4)=rho(1)!logg(1)!Lum(1)
      		ageb(1)=age(1)
      		ageb(2)=age(1)
      		rconvb(1)=rconv(1)
      		rconvb(2)=rconv(1)
      		Radb(1)=Rad(1)
      		Radb(2)=Rad(2)
        endif
        call getbounds(npt,age,Teff,rho,Rad,bounds,ageb,rconv,rconvb,
     .		Radb)
      	close(11) !closing the file and freeing the unit number
      	nline=nline+1
      	goto 10
 11   continue
      write(0,*) "Temp: ",bounds(1),bounds(2)
      write(0,*) "rho : ",bounds(3),bounds(4)
      write(0,*) "Age: ",ageb(1)**2,ageb(2)**2
 
C     initialize PGPLOT device
      call pgopen('?')
C     initialize multipage (may not be necessary)
      call pgpage()
C     prompt when starting new page
      call pgask(.true.)
      call PGPAP ( 6.0 ,1.0) !square plot
c      call pgsch(2.0) !bigger text

      sf=0.05
      temp(1)=bounds(1)-sf*(bounds(2)-bounds(1))
      temp(2)=bounds(2)+sf*(bounds(2)-bounds(1))
      temp(3)=bounds(3)-sf*(bounds(4)-bounds(3))
      temp(4)=bounds(4)+sf*(bounds(4)-bounds(3))
      
      do 12 i=1,4
        bounds(i)=temp(i)
 12   continue
      
c      call pgwindow(log10(12500.0),log10(3500.),bounds(4),bounds(3))
      call pgwindow(log10(12500.0),log10(3500.),1.2,-3.0)
      call pgslw(3)
      call pgbox('BCLNTS1',0.0,0,'BCLNTS',0.0,0)
      call pgslw(1)
      
      rewind(10) !read in filenames from start of file
      nline=1 !nline defines the current line in filelist    
C     Note that filename is currently constrained to 80 characters.
C     Refine the variable filename above to expand.
 13   read(10,500,end=14) filename !the end statement defines our loop
C       open up the HR file for reading
      	open(unit=11,file=filename,status='old',err=903)
C       retrieve data from the file
      	call getdata(nmax,npt,age,Teff,Lum,Rad,rconv,mass(nline),logg,
     .      rho)
      	call plotdata(npt,age,Teff,rho,Rad,ageb,rconv,rconvb,Radb)
      	nline=nline+1
      	goto 13
 14   continue
      	
      
      call pgclos()
      close(10)
      
      goto 999
 901  write(0,*) "Usage: hrdisp <filelist>"
      write(0,*) "  where <filelist> defines the HR files to plot"
      goto 999
 902  write(0,*) "Cannot open ",filelist
      goto 999
 903  write(0,*) "Error reading ",filelist
      write(0,*) "Line # ",nline
      write(0,*) "Cannot open ",filename
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotdata(npt,age,Teff,Lum,Rad,ageb,rconv,rconvb,Radb)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,CI1,CI2,P,i,icolor
      real age(npt),Teff(npt),Lum(npt),ageb(2),SIGN,CONTRA,BRIGHT,x(2),
     .	y(2),rconv(npt),rconvb(2),Rad(npt),Radb(2),temp,dumr,rhos
      
      call pgline(npt,Teff,Lum)
      
      CI1 = 16 
      CI2 = 86
      P      = 2
      SIGN   = 1.0
      CONTRA = 1.0!0.56
      BRIGHT = 0.5 
      call PGSCIR ( CI1 , CI2 )         
      call Palett ( P , SIGN*CONTRA , BRIGHT )
      call PGSITF ( 0 )
      
      call pgslw(2)
      call pgbbuf()
      do 10 i=2,npt
c      	icolor=int((CI2-CI1)*(age(i)-ageb(1))/(ageb(2)-ageb(1)))+CI1
c      	if(icolor.gt.CI2) icolor=CI2
      	icolor=int((CI2-CI1)*(rconv(i)-rconvb(1))/
     .		(rconvb(2)-rconvb(1)))+CI1
c      	icolor=int((CI2-CI1)*(Rad(i)-Radb(1))/
c     .		(Radb(2)-Radb(1)))+CI1
      	call pgsci(icolor)
c      	call pgpt1(Teff(i),Lum(i),17)
        x(1)=Teff(i-1)
        x(2)=Teff(i)
        y(1)=Lum(i-1)
        y(2)=Lum(i)
        call pgline(2,x,y)
c        write(6,*) age(i),teff(i),lum(i)
 10   continue
      call pgebuf()
 
      call pgsci(1)
      call pgslw(3)
c      call pglabel("Teff","Log(L/Lsun)","")   
c      call pglabel("Teff (K)", "logg (cgs)","")
      call pglabel("Teff (K)", "rho (g/cm^3)","")
      call pgslw(2)
c      call pgpt1(log10(5777.0),0.0,9)

C     Multis
c      call pgpt1(log10(6500.0),log10(0.80),9)!152
c      call pgpt1(log10(5750.0),log10(1.12),9)!157
c      call pgpt1(log10(5500.0),log10(1.91),9)!191
c      call pgpt1(log10(5750.0),log10(0.40),9)!209
c      call pgpt1(log10(5206.0),log10(2.31),9)!896
c      call pgpt1(log10(4211.0),log10(5.10),9)!877

C     KOI-218B
c      call pgpt1(log10(6200.0),log10(0.2732),9)

      open(unit=11,file='koirhos.txt',status='old')
  11  read(11,*,end=12) temp,dumr,dumr,rhos
        call pgpt1(log10(temp),log10(rhos),9)
        goto 11
  12  continue
      close(11)
     
c      call pgpt1(log10(4980.0),log10(3.2),9)!HD189733
c      call pgpt1(log10(5942.0),log10(1.397),9)!HD209458
c      call pgpt1(log10(5545.0),log10(2.103),9)!KOI217

c      call pgpt1(log10(10000.0),4.00,9)
c      call pgpt1(log10(9400.0),4.0,9)
c      call pgpt1(log10(5000.0),log(0.64),9)
      call pgslw(1)
      
      return
      end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getdata(nmax,npt,age,Teff,Lum,Rad,rconv,mass,logg,rho)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nmax,npt,i,j,k,dumi,izone,nelem
      real age(nmax),Teff(nmax),Lum(nmax),Rad(nmax),dumr,rconv(nmax),
     .  mass,logg(nmax),G,Rsun,Msun,rho(nmax),Pi
      character cflag
      
      G=6.67259E-8
      Rsun=6.9599E10
      Msun=1.989E33
      Pi=acos(-1.0)!define Pi and 2*Pi
      
      i=1
 10   read(11,*,end=11,err=11) age(i),nelem,izone,cflag
c
        age(i)=age(i)!+1.0
        age(i)=sqrt(age(i))
        read(11,500) Teff(i),Lum(i),Rad(i),dumr,dumr,rconv(i)
c        Teff(i)=10.0**Teff(i)
 500    format(6E13.6)
c        write(6,*) Teff(i),Lum(i),Rad(i)
        logg(i)=log10(G*mass*Msun/(Rsun*Rsun*(10**Rad(i))**2))
        rho(i)=mass*Msun/(Pi*4.0/3.0*(Rsun*10**Rad(i))**3.0)
        rho(i)=log10(rho(i)+0.0001)
c        write(6,*) mass,10**rad(i),rho(i)
c        read(5,*)
c        if(izone.eq.
        if(izone.eq.2) then
        	read(11,500) dumr,dumr,rconv(i)
        else
        	read(11,500) dumr
        endif
        rconv(i)=10.0**(log10(rconv(i))-Rad(i))
c        write(6,*) age(i)**2,cflag,rconv(i),izone
        if(izone.ge.3) read(11,500) dumr!,dumr,rconv2
        k=nelem
        do 12 j=1,k
      	  read(11,*)
 12     continue
        i=i+1
        if(i.gt.nmax)then
        	write(6,*) "Warning, increase nmax"
        	goto 11
        endif
        goto 10
 11   continue
      npt=i-1
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getbounds(npt,age,Teff,Lum,Rad,bounds,ageb,rconv,
     .	rconvb,Radb)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i
      real age(npt),Teff(npt),Lum(npt),bounds(4),ageb(2),agemax,
     .	rconv(npt),rconvb(2),Rad(npt),Radb(2)
      
      agemax=sqrt(14000.0)
      
      do 10 i=1,npt
      	bounds(1)=min(Teff(i),bounds(1))
      	bounds(2)=max(Teff(i),bounds(2))
      	bounds(3)=min(Lum(i),bounds(3))
      	bounds(4)=max(Lum(i),bounds(4))
      	ageb(1)=min(age(i),ageb(1))
      	ageb(2)=max(age(i),ageb(2))
      	rconvb(1)=min(rconv(i),rconvb(1))
      	rconvb(2)=max(rconv(i),rconvb(2))
      	Radb(1)=min(Rad(i),Radb(1))
      	Radb(2)=max(Rad(i),Radb(2))
 10   continue
 
      if(ageb(2).gt.agemax) ageb(2)=agemax
 
      return
      end

c     --------------------------------------------
      SUBROUTINE PALETT ( TYPE , CONTRA , BRIGHT )
c     --------------------------------------------
c
c     Sets a "palette" of colors in the range of defined color indices 
c     This subroutine is distributed with PGPLOT in one of the demos
c
      INTEGER TYPE
      REAL CONTRA, BRIGHT
C
      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL CL(5), CR(5), CG(5), CB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)
C
      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
C
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
C
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
C
      DATA CL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA CR /0.0, 0.0, 0.0, 0.3, 1.0/
      DATA CG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA CB /0.0, 0.5, 1.0, 1.0, 1.0/
C
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
C
      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
     :         0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
     :         0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
     :         0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
     :         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
C
      IF (TYPE.EQ.1) THEN
C        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.2) THEN
C        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.3) THEN
C        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.4) THEN
C        -- freeze
         CALL PGCTAB(CL, CR, CG, CB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.5) THEN
C        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
      END IF
      END     	
      