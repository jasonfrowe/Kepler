CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program limbdarkening
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer iargc,nmodel,tabteff,tablogg,tabfeh,i,j,k
      parameter(nmodel=3800)
      double precision logg,teff,temps(76),loggs(11),
     .  fehs(8),temp(7),nl(4),nl1(76,11,8),nl2(76,11,8),nl3(76,11,8),
     .  nl4(76,11,8)
      character*80 file1,file2,cline,pm,dumc
      data fehs /-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.2,0.5/

C     set up grid indices for Teff,logg and [Fe/H]      
      do 12 i=1,39
        temps(i)=3500.0+250.0*(i-1)
 12   continue
      do 13 i=1,37
        j=i+39
        temps(j)=13000.0+1000.0*i
 13   continue
      do 14 i=1,11
        loggs(i)=0.5*(i-1)
 14   continue
      
      if(iargc().lt.2) goto 901
      call getarg(1,cline)
      read(cline,*) teff
      call getarg(2,cline)
      read(cline,*) logg
      if(teff.lt.10.0) teff=5700.0
      if(logg.lt.0.1) logg=4.5
      
c      write(6,*) "Teff: ",teff
c      write(6,*) "logg: ",logg
      
c      file1="/home/rowe/p555/transitfit/kepler.list"
      file1="/home/rowe/p555/transitfit/models.list"
c      file2="/home/rowe/p555/transitfit/kepler.new"
      file2="/home/rowe/p555/transitfit/kepler.ch_25_26.ld"
c      file2="/home/rowe/p555/transitfit/table.most"
      
      open(unit=10,file=file1,status='old',err=902)
      open(unit=11,file=file2,status='old',err=903)
      do 5 i=1,7
        read(11,502) dumc
 5    continue
 502  format(A1)
      
 10   read(10,501,end=11) tabteff,tablogg,tabfeh
 501  format(I5,1X,I2,1X,I3)
      read(11,*,end=11) (temp(i),i=1,7),(nl(j),j=1,4)
C       get array position 
        if(tabteff.le.13000)then
            i=(tabteff-3500)/250+1
c            write(6,*) i,temps(i),tabteff
        elseif(tabteff.gt.13000)then
            i=(tabteff-13000)/1000+39
c            write(6,*) i,temps(i),tabteff
        endif
        j=tablogg/5+1
c        write(6,*) j,loggs(j),tablogg
        if(tabfeh.le.0)then
            k=(tabfeh+25)/5+1
c            write(6,*) k,fehs(k),tabfeh
        elseif(tabfeh.eq.2)then
            k=7
        elseif(tabfeh.eq.5)then
            k=8
        endif
        nl1(i,j,k)=nl(1)
        nl2(i,j,k)=nl(2)
        nl3(i,j,k)=nl(3)
        nl4(i,j,k)=nl(4)
c        write(6,*) nl(1),nl(2),nl(3),nl(4)
c        read(5,*)
      goto 10
 11   continue
     
      close(10)
      
      if(teff.le.13000.0)then
        i=(int(teff)-3500)/250+1
      elseif(teff.gt.13000)then
        i=(int(teff)-13000)/1000+39
      endif
      j=int(10.0*logg)/5+1
      
      k=6
c      write(6,*) i,j,k
      write(6,500) nl1(i,j,k),nl2(i,j,k),nl3(i,j,k),nl4(i,j,k)
 500  format(4(F7.4,1X))
      
      goto 999
 901  write(0,*) "Usage: limbdarkening teff logg"
      write(0,*) "Where Teff is in K and logg is cgs"
      goto 999
 902  write(0,*) "Cannot open ",file1
      goto 999
 903  write(0,*) "Cannot open ",file2
      goto 999
 999  end
      
