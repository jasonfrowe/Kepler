CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine starplot(nplanet,nfit,sol,srad,Teff,logg,ntype,b,order)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer ncolour,i,Teff,j,nplanet,nfit,col,ii
      parameter(ncolour=37)
      integer colourtemp(ncolour),sR(ncolour),sG(ncolour),sB(ncolour),
     .   CI1,CI2,inl(3),icol,ntype(nplanet),order(nplanet)
      real plotsize,bounds(4),srad,In(3),nl(4),Intensity,r,mu,x,y,br,
     .   rp,shiftx,shifty,neg,dx,dy,th,dth,pi,srad2
      double precision nl1(76,11,8),nl2(76,11,8),nl3(76,11,8),
     .  nl4(76,11,8),logg,sol(nfit),b(nplanet)
      character*80 dumc,cfile

      Pi=acos(-1.e0) !real

c      srad=1.241 !real
c      Teff=5794  !int
c      logg=4.2807 !double

      call limbdarkensetup(nl1,nl2,nl3,nl4)

C      cfile="/home/rowe/p555/transitfit/Kepler/starcoloursRGB.txt"
      cfile="/Users/rowe/Documents/transitfit/starcoloursRGB.txt"

      open(unit=10,file=cfile,status='old',err=905)
      do 15 i=1,ncolour
        read(10,*) dumc,colourtemp(i),sR(i),sG(i),sB(i)
 15   continue
      close(10)

      plotsize=4.0
      bounds(1)=0.0
      bounds(2)=plotsize*1.4
      bounds(3)=0.0
      bounds(4)=plotsize

      x=(bounds(1)+bounds(2))/2.0
      y=(bounds(3)+bounds(4))/2.0

      call pgwindow(bounds(1),bounds(2),bounds(3),bounds(4))
      call pgscr(14,0.0,0.0,0.0)
      call pgsci(14)
      call pgrect(0.0,bounds(2),0.0,bounds(4))

      CI1 = 16
      CI2 = 86
      call PGSCIR(CI1,CI2)

      call starcolour(Teff,In,ncolour,colourtemp,sR,sG,sB)

      if(Teff.le.13000.0)then
         inl(1)=(Teff-3500)/250+1
      elseif(Teff.gt.13000)then
         inl(1)=(Teff-13000)/1000+39
      endif
      inl(2)=int(10.0*logg)/5+1
      inl(3)=6 ![Fe/H]=0

      nl(1)=real(nl1(inl(1),inl(2),inl(3)))
      nl(2)=real(nl2(inl(1),inl(2),inl(3)))
      nl(3)=real(nl3(inl(1),inl(2),inl(3)))
      nl(4)=real(nl4(inl(1),inl(2),inl(3)))

      do 12 i=ci1,ci2
         Intensity=real(i-ci1)/real(ci2-ci1)
         call pgscr(i,In(1)*Intensity/255.0,In(2)*Intensity/255.0,
     .    In(3)*Intensity/255.0)
 12   continue

      do 10 i=ci2,ci1,-1
         r=real(i-ci1)/real(ci2-ci1)
         mu=sqrt(1.0-r*r)
         Intensity=1.0
         do 11 j=1,4
            Intensity=Intensity-nl(j)*(1.0-mu**(real(j)/2.0))
 11      continue
c           write(0,*) i,r,Intensity
         icol=Intensity*(ci2-ci1)+ci1
         call pgsci(icol)
         call pgcirc(x,y,r*srad)
c         write(6,*) r*srad
 10   continue

CCCCC Okay.. extra star here...

c      srad2=1.0*srad
c      do 18 i=ci2,ci1,-1
c         r=real(i-ci1)/real(ci2-ci1)
c         mu=sqrt(1.0-r*r)
c         Intensity=1.0
c         do 19 j=1,4
c            Intensity=Intensity-nl(j)*(1.0-mu**(real(j)/2.0))
c 19      continue
cc           write(0,*) i,r,Intensity
c         icol=Intensity*(ci2-ci1)+ci1
c         call pgsci(icol)
c         call pgcirc(x+2.0*srad2,y+y/3.0,r*srad2)
cc         write(6,*) r*srad
c 18   continue



      shiftx=0.0
      shifty=-1.0


c      th=pi/2.0
c      dth=Pi/real(nplanet+2)

c      do 17 i=1,nplanet
      do 17 ii=1,nplanet
         i=order(ii)
         col=10*(i-1)

c         br=min(real(sol(11+col)),0.99)
         br=min(real(b(i)),0.99)
         rp=real(sol(12+col))*srad

c         th=th+dth
c         dx=br*sin(th)*srad
c         dy=br*cos(th)*srad
         dth=(Pi-2.0*asin(br))/real(nplanet+2)
         th=asin(br)+dth*i
         dx= br*srad/tan(th)
c         dx= cos(th)*srad*0.9
         dy=-br*srad
c         write(6,*) th,dth,dx,dy


 16      if(ntype(i).eq.2)then
            call pgscr(14,0.0,0.0,0.0)
         else
            call pgscr(14,0.8,0.8,0.8)
         endif
         call pgsci(14)
cc         if(i.eq.2)then
cc            call pgcirc(x+2.0*srad2+dx,y+y/3.0+dy,rp)
cc         else
            call pgcirc(x+dx,y+dy,rp)
cc         endif
 17   continue

      call pgsci(1)

      goto 999
 905  write(6,*) "Cannot open starcolourRGB.txt"
      pause
      goto 999
 999  return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine starcolour(Teff,In,ncolour,colourtemp,sR,sG,sB)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer Teff,ncolour,colourtemp(ncolour),sR(ncolour),sG(ncolour),
     .  sB(ncolour),i
      real In(3)
      logical loop

      if(Teff.ge.colourtemp(1))then
        In(1)=real(sR(1))
        In(2)=real(sG(1))
        In(3)=real(sB(1))
      elseif(Teff.le.colourtemp(ncolour))then
        In(1)=real(sR(ncolour))
        In(2)=real(sG(ncolour))
        In(3)=real(sB(ncolour))
      else
        loop=.true.
        i=1
        do while(loop)
           if((Teff.lt.colourtemp(i)).and.(Teff.ge.colourtemp(i+1)))then
                In(1)=real(sR(i))
                In(2)=real(sG(i))
                In(3)=real(sB(i))
                loop=.false.
            endif
            i=i+1
        enddo
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine limbdarkensetup(nl1,nl2,nl3,nl4)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer tabteff,tablogg,tabfeh,i,j,k
      double precision dumc,nl(4),temp(7),nl1(76,11,8),nl2(76,11,8),
     .  nl3(76,11,8),nl4(76,11,8)
      character*80 file1,file2

c      file1="/home/rowe/p555/transitfit/models.list"
      file1="/Users/rowe/Documents/transitfit/models.list"
c      file2="/home/rowe/p555/transitfit/kepler.ch_25_26.ld"
      file2="/Users/rowe/Documents/transitfit/kepler.ch_25_26.ld"


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
