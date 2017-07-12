      subroutine limbdarkensetup(nl1,nl2,nl3,nl4)
      implicit none
      integer tabteff,tablogg,tabfeh,i,j,k
      double precision dumc,nl(4),temp(7),nl1(76,11,8),nl2(76,11,8),
     .  nl3(76,11,8),nl4(76,11,8)
      character*80 file1,file2

c      file1="/home/rowe/p555/transitfit/kepler.list"
      file1="./models.list"
c      file2="/home/rowe/p555/transitfit/kepler.new"
      file2="./kepler.ch_25_26.ld"
c      file2="/home/rowe/p555/transitfit/table.most"
      
      open(unit=12,file=file1,status='old',err=902)
      open(unit=13,file=file2,status='old',err=903)
      do 5 i=1,7
        read(13,502) dumc
 5    continue
 502  format(A1)
      
 10   read(12,501,end=11) tabteff,tablogg,tabfeh
 501  format(I5,1X,I2,1X,I3)
      read(13,*,end=11) (temp(i),i=1,7),(nl(j),j=1,4)
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
     
      close(12)
      close(13)
      
      goto 999
 901  write(0,*) "Usage: limbdarkening teff logg"
      write(0,*) "Where Teff is in K and logg is cgs"
      goto 999
 902  write(0,*) "Cannot open ",file1
      goto 999
 903  write(0,*) "Cannot open ",file2
      goto 999
 999  return
      end