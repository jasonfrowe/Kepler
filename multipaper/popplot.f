      program popplot
      implicit none
      integer nmax,n,nunit,ndata,i,j,ndyn,nex,nkex,nsdata,ns,nfp,k,
     .   nsmall,nast,pflag,naused,ncolour,ci1,ci2
      parameter(nmax=5000,ndata=26,nex=17,nkex=3,nsdata=24,ncolour=37)
      integer excl(nmax),nkoi(nmax),astkoi(nmax),colourtemp(ncolour),
     .   sR(ncolour),sG(ncolour),sB(ncolour)
      real data(nmax,ndata),dynkoi(nmax),pexcl(nex),kexcl(nkex),dumr,
     .   sdata(nmax,nsdata),koifp,S,Tsun,pars2(nsdata),koi2,px,py,pe,
     .   hzx(2),hzy(2),seff,In(3),Intensity
      character dynchar
      character*80 multiprops,dumc,dynfile,sfile,fpfile,sfile2,
     .   sparsfile,cfile
      data pexcl/375.01,422.01,435.02,490.02,771.01,1032.01,1096.01,
     .   1192.01,1206.01,1274.01,1356.01,1421.01,1463.01,1174.01,
     .   1208.01,1008.01,99.01/
      data kexcl/1750.01,1750.02,1955.03/

      Tsun=5781.6

      nunit=10
      dynfile="dynamictest.20130905.dat"
      open(unit=nunit,file=dynfile,status='old',err=904)
      read(nunit,*) dumc
      i=1
 10   read(nunit,*,end=11) dynkoi(i),dumr,dumr,dumr,dynchar
         if(dynchar.eq.'C')then
c            write(0,*) dynkoi(i)
            i=i+1
         endif
      goto 10
 11   continue
      ndyn=i-1
      close(nunit)


      multiprops="multiprops.dat" !file that contains multi-analysis

      nunit=10
C     open up multiprops file that has the multi-planet analysis
      open(unit=nunit,file=multiprops,status='old',err=901)
      read(nunit,*) dumc !first line is a header
      i=1
 14   read(nunit,*,end=15) (data(i,j),j=1,ndata)
c         write(0,*) i,koi(i)  !this line was for debugging

         do 70 j=1,nex
            if(pexcl(j).eq.data(i,1))then
c               write(0,*) pexcl(j),data(n,1)
               data(i,4)=-1.0
             endif
 70      continue

         do 125 j=1,nkex
            if(kexcl(j).eq.data(i,1))then
c               write(0,*) kexcl(j),data(n,1)
               excl(i)=1
            endif
 125     continue

         do 128 j=1,ndyn
            if(dynkoi(j).eq.data(i,1))then
c               write(0,*) dynkoi(j),data(n,1)
               excl(i)=1
            endif
 128     continue


         i=i+1
      goto 14
 15   continue
      close(nunit)
      n=i-1

      call pgopen('?')
      call pgpage()
      call PGPAP ( 8.0 ,1.0)
      call pgsch(1.5)
      call pgslw(2)

      call pgvport(0.15,0.85,0.15,0.85)
      call pgwindow(log10(0.1),log10(5000.0),log10(0.2),log10(20.0))
      call pgbox('BCLNTS1',0.0,0,'BCLNTS1',0.0,0)
      call pglabel("S (S\d\(2284)\u)","Radius (R\d\(2284)\u)","")


      cfile="/Users/rowe/Documents/transitfit/starcoloursRGB.txt"
      open(unit=nunit,file=cfile,status='old',err=907)
      do 61 i=1,ncolour
        read(nunit,*) dumc,colourtemp(i),sR(i),sG(i),sB(i)
 61   continue
      close(nunit)
      CI1 = 16
      CI2 = 86
      do 60 j=ci1,ci2
         i=3000+int((7000.0-3000.0)*(real(j)-real(ci1))/real(ci2-ci1))
         call starcolour(i,In,ncolour,colourtemp,sR,sG,sB)
c         write(0,*) in(1),in(2),in(3)
         call pgscr(j,In(1)/255.0,In(2)/255.0,In(3)/255.0)
         hzx(1)=log10(seff(real(i),1))
         hzx(2)=hzx(1)
         hzy(1)=log10(0.2)
         hzy(2)=log10(20.0)
c         write(0,*) i,hzx(1)!,10**hzx(1)
         call pgsci(j)
         call pgline(2,hzx,hzy)
         hzx(1)=log10(seff(real(i),5))
         hzx(2)=hzx(1)
         call pgline(2,hzx,hzy)
 60   continue
      call pgsci(1)
      call pgbox('BCLNTS1',0.0,0,'BCLNTS1',0.0,0)

C     read and plot singles
      sfile="koi_multi.20120926.dat"
      open(unit=nunit,file=sfile,status='old',err=902)
      read(nunit,*,err=902) dumc
      read(nunit,*,err=902) dumc
      i=1 !count number of lines read
 18   read(nunit,*,end=17,err=16) (sdata(i,j),j=1,nsdata)
         if(sdata(i,1).gt.3058.01) goto 18
         if(floor(sdata(i,1)).eq.sdata(i,1)) sdata(i,1)=sdata(i,1)+0.01
         i=i+1
      goto 18
 16   write(0,*) "Error on line: ",i+2
      goto 18
 17   continue
      close(nunit)
      ns=i-1
      write(0,*) "Number of lines read: ",ns

      fpfile="Fergalfile.20130801.q1q8.awk.txt"
      open(unit=nunit,file=fpfile,status='old',err=906)
 19   read(nunit,*,end=21) koifp,nfp
         do 20 i=1,ns
            if(koifp.eq.sdata(i,1))then
               sdata(i,24)=nfp
            endif
 20      continue
      goto 19
 21   continue
      close(nunit)

      sfile2="koi_characteristics.20130509.dat"
      open(unit=nunit,file=sfile2,status='old',err=905)
      read(nunit,*,err=904) dumc
      read(nunit,*,err=904) dumc
      i=1 !count number of lines read
 51   read(nunit,*,end=52,err=53) koi2,(pars2(j),j=1,20)
         if(koi2.gt.3058.01) goto 51
         if(floor(koi2).eq.koi2) koi2=koi2+0.01
         do 54 j=1,ns
            if(sdata(j,1).eq.koi2) then
               sdata(j,21)=pars2(20) !SN replacement
c               write(0,*) koi(j),per(i),pars2(5)
               sdata(j,6)=pars2(5) !period

               do 55 k=1,nex
                  if(pexcl(k).eq.sdata(j,1))then
c                     write(0,*) pexcl(k),koi(j)
                     sdata(j,6)=-1.0
                  endif
 55            continue

            endif
 54      continue
         i=i+1
      goto 51
 53   write(0,*) "Error on line: ",i+2
      goto 51
 52   continue
      close(nunit)

      call pgbbuf()
      do 22 i=1,ns
         if((sdata(i,24).lt.1).and.(sdata(i,21).gt.10.0).and.
     .    (sdata(i,14)+sdata(i,13).lt.0.98).and.
     .    (sdata(i,6).gt.1.6))then
            S=sdata(i,9)**2.0*(sdata(i,7)/Tsun)**4.0*
     .       (sdata(i,15))**-2.0
            call pgpt1(log10(S),log10(sdata(i,4)),1)
         endif
 22   continue
      call pgebuf()

C     FP,SN,P,b cuts

      do 44 i=1,nmax
         nkoi(i)=0
 44   continue

      do 45 i=1,n
         if((data(i,24).gt.1).and.(data(i,21).ge.10.0).and.
     .    (data(i,6)+data(i,8).lt.0.98).and.(data(i,4).gt.1.6))then
           nkoi(int(data(i,1)))=nkoi(int(data(i,1)))+1
         endif
 45   continue

      nsmall=0
      call pgbbuf()
      do 47 i=1,n
c         if((data(i,24).gt.1).and.(data(i,21).ge.10.0).and.
c     .    (data(i,6)+data(i,8).lt.0.98).and.(data(i,4).gt.1.6).and.
c     .    (nkoi(int(data(i,1))).gt.1).and.(data(i,20).gt.0.0))then
         if((data(i,24).gt.1).and.(data(i,21).ge.10.0).and.
     .    (data(i,6)+data(i,7)+data(i,8).lt.1.00).and.
     .    (data(i,4).gt.1.6).and.
     .    (nkoi(int(data(i,1))).gt.1).and.(data(i,20).gt.0.0))then

            if(data(i,23).gt.3)then
               if((data(i,20).gt.50).and.(data(i,18).gt.8.0))
     .            write(0,*) data(i,1),data(i,20),data(i,18)
               call pgpt1(log10(data(i,20)),log10(data(i,18)),17)
               if(data(i,20).lt.2.0) write(0,*) data(i,1),data(i,18),
     .          data(i,20)
               if((data(i,18).lt.1.25).and.(data(i,18).gt.0))
     .          nsmall=nsmall+1
            else
               call pgpt1(log10(data(i,20)),log10(data(i,18)),21)
            endif
         endif
 47   continue
      call pgebuf()
      write(6,*) "R<1.25: ",nsmall

C     Make a plot of rhostar vs rhoc-rhostar

      call pgpage()
      call pgvport(0.15,0.85,0.15,0.85)
      call pgwindow(log10(0.01),log10(10.0),-10.0,10.0)
      call pgbox('BCLNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("\(2143)\d*\u (g/cm\u3\d)",
     .   "(\(2143)\d*\u - \(2143)\dc\u)/\gs","")

      sparsfile="spars.dat"
      open(unit=nunit,file=sparsfile,status='old',err=903)
      i=1
      read(nunit,*) dumc !first line is a comment
 57   read(nunit,*,end=56) astkoi(i),(dumr,j=1,12),pflag
         if(pflag.eq.5) i=i+1
      goto 57
 56   continue
      nast=i-1
      close(nunit)

      naused=0
      do 58 j=1,nast
         do 59 i=1,n
            px=log10(data(i,12))
            py=(data(i,12)-data(i,10))
            pe=sqrt(data(i,11)*data(i,11)+data(i,13)*data(i,13))
            py=py/pe
            if(int(data(i,1)).eq.astkoi(j))then
               call pgpt1(px,py,17)
c               call pgerr1(6,px,py,pe,1.0)
               naused=naused+1
            else
               call pgpt1(px,py,1)
            endif
 59      continue
 58   continue
      write(0,*) "Number of targets with astero: ",naused

      call pgclos()


      goto 999
 901  write(0,*) "Cannot open ",multiprops
      goto 999
 902  write(0,*) "Cannot open ",sfile
      goto 999
 903  write(0,*) "Cannot open ",sparsfile
      goto 999
 904  write(0,*) "Cannot open: ",dynfile
      goto 999
 905  write(0,*) "Cannot open ",sfile2
      goto 999
 906  write(0,*) "Cannot open ",fpfile
      goto 999
 907  write(0,*) "Cannot open ",cfile
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function seff(teff,type)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer ntypemax,type
      parameter(ntypemax=5)
      real seffsun(ntypemax),a(ntypemax),b(ntypemax),c(ntypemax),
     .   d(ntypemax),teff,Ts
      data seffsun/1.7753,1.0512,1.0140,0.3438,0.3179/
      data a/1.4316e-4,1.3242e-4,8.1774e-5,5.8942e-5,5.4513e-5/
      data b/2.9875e-9,1.5418e-8,1.7063e-9,1.6558e-9,1.5313e-9/
      data c/-7.5702e-12,-7.9895e-12,-4.3241e-12,-3.0045e-12,
     .   -2.7786e-12/
      data d/-1.1635e-15,-1.8328e-15,-6.6462e-16,-5.2983e-16,
     .   -4.8997e-16/

C     1: Venus, 2: Runaway Greenhouse, 3:Moist Greenhouse,
C     4:Maximum Greenhouse, 5: Early Mars

      Ts=Teff-5780.0

      if((type.gt.0).and.(type.le.ntypemax))then
         seff=seffsun(type)+a(type)*Ts+b(type)*Ts*Ts+c(type)*Ts**3.0+
     .    d(type)*Ts**4.0
      else
         seff=0.0
      endif

      return
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
