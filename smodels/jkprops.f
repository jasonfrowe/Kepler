      program jkprops
      implicit none
C     Given J-K calculate Teff,logg,R*
      integer i,nmax,npt,iargc
      parameter(nmax=16)
      real teffs(nmax),jhs(nmax),hks(nmax),jks(nmax),ms(nmax),rs(nmax),
     .  Jmag,Kmag,JK,G,Rsun,Msun,loggs(nmax),teff,logg,rad
      character*80 filename,dumc,cline
      
      G=6.67259E-8
      Rsun=6.9599E10
      Msun=1.989E33
      
      if(iargc().lt.2) goto 902
      
      call getarg(1,cline)
      read(cline,*,err=902) Jmag
      call getarg(2,cline)
      read(cline,*,err=902) Kmag
      
      JK=Jmag-Kmag
      
      
      filename='/Users/rowe/Documents/transitfit/schmidtkaler.dat'
      
      open(unit=10,file=filename,status='old',err=901)
      read(10,*) dumc
      
      i=1
 10   read(10,*,end=11) dumc,teffs(i),jhs(i),hks(i),ms(i),rs(i)
        jks(i)=jhs(i)+hks(i)
        loggs(i)=log10(G*ms(i)*Msun/(Rsun*Rsun*rs(i)*rs(i)))
        i=i+1
        goto 10
 11   continue
      close(10)
      npt=i-1
      
      if(JK.gt.jks(npt)) write(6,*) teffs(npt),loggs(npt),rs(npt)
      if(JK.lt.jks(1)) write(6,*) teffs(1),loggs(1),rs(1)
      do 12 i=1,npt-1
        if((JK.gt.jks(i)).and.(JK.lt.jks(i+1)))then
            Teff=(JK-jks(i))/(jks(i+1)-jks(i))*(Teffs(i+1)-Teffs(i))+
     .          Teffs(i)
            logg=(JK-jks(i))/(jks(i+1)-jks(i))*(loggs(i+1)-loggs(i))+
     .          loggs(i)
            rad=(JK-jks(i))/(jks(i+1)-jks(i))*(rs(i+1)-rs(i))+rs(i)
            write(6,*) Teff,logg,rad
        endif
 12   continue
        
      
      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 902  write(0,*) "Usage: jkprops Jmag Kmag"
      goto 999
 999  end
