      program tdurplots
      implicit none
      integer iargc,nunit,i,nmax,flag,nc,nmc,nfit,npt,j,k,ii,jj,nbin,
     .   kk
      parameter(nmax=500000,nfit=5)
      integer nbstar(nmax),iorder(nfit)
      real chi(nmax),ans(nfit,nmax),minmax(nfit,2),ncut,bx(nmax),
     .   by(nmax),px(nmax),bmax,lx(2),ly(2)
      character*80 mcmcfile,title(nfit)
      data iorder/1,5,2,3,4/

      ncut=26.9

      title(1)="binfrac"
      title(2)="x"
      title(3)="y"
      title(4)="errfrac"
      title(5)="n"

      if(iargc().lt.1) goto 901

      call getarg(1,mcmcfile)

      nunit=10
      open(unit=nunit,file=mcmcfile,status='old',err=902)
      i=1
      read(nunit,500) flag,nc,nmc,chi(i),(ans(j,i),j=1,nfit),
     . nbstar(i)
      do 12 j=1,nfit
         minmax(j,1)=ans(j,i)
         minmax(j,2)=ans(j,1)
 12   continue
      i=i+1
 10   read(nunit,500,end=11) flag,nc,nmc,chi(i),(ans(j,i),j=1,nfit),
     . nbstar(i)
 500     format(I1,2(1X,I1),6(1X,1PE17.10),1X,I5)
         do 13 j=1,nfit
            minmax(j,1)=min(ans(j,i),minmax(j,1))
            minmax(j,2)=max(ans(j,i),minmax(j,2))
 13      continue
         i=i+1
      goto 10
 11   continue
      close(nunit)
      npt=i-1
      write(0,*) "npt: ",npt

      call pgopen('?')
      call pgpage()
      call PGPAP ( 8.0 ,1.0)
      call pgsubp(nfit+2,nfit+2)
      call pgsch(3.0)
      call pgslw(2)

      nbin=20
      do 14 ii=1,nfit
         i=iorder(ii)
         do 15 jj=ii,nfit
            j=iorder(jj)
            call pgpanl(ii+1,jj+1)
            if(i.eq.j)then
               do 17 k=1,npt
                  px(k)=ans(i,k)
 17            continue
               call bindata(nbin,npt,px,bx,by,minmax(i,1),minmax(i,2),
     .          bmax)

               call pgvport(0.0,1.0,0.0,1.0)
               call pgwindow(minmax(i,1),minmax(i,2),0.0,bmax+0.1*bmax)

               lx(1)=bx(1)
               lx(2)=bx(1)
               ly(1)=0.0
               ly(2)=by(1)
               call pgline(2,lx,ly)
               do 18 kk=1,nbin-1
                  lx(1)=bx(kk)
                  lx(2)=bx(kk+1)
                  ly(1)=by(kk)
                  ly(2)=by(kk)
                  call pgline(2,lx,ly)
                  lx(1)=bx(kk+1)
                  lx(2)=bx(kk+1)
                  ly(1)=by(kk)
                  ly(2)=by(kk+1)
                  call pgline(2,lx,ly)
 18            continue
c               call pgline(nbin,bx,by)
               call pgwindow(minmax(i,1),minmax(i,2),0.0,1.1)
               call pgbox('BCTS',0.0,0,'BCTS',0.0,0)
               if(ii.eq.1) then
c                  call pgbox('BCTS',0.0,0,'BCNTS1',0.0,0)
                  call pglabel("",title(j),"")
               endif
               if(jj.eq.nfit) then
                  call pgbox('BCNTS1',0.0,0,'BCTS',0.0,0)
                  call pglabel(title(i),"","")
               endif
            else
               call pgvport(0.0,1.0,0.0,1.0)
               call pgwindow(minmax(i,1),minmax(i,2),minmax(j,1),
     .          minmax(j,2))
               call pgbox('BCTS',0.0,0,'BCTS',0.0,0)
               if(ii.eq.1) then
                  call pgbox('BCTS',0.0,0,'BCNTS1',0.0,0)
                  call pglabel("",title(j),"")
               endif
               if(jj.eq.nfit) then
                  call pgbox('BCNTS1',0.0,0,'BCTS',0.0,0)
                  call pglabel(title(i),"","")
               endif
               call pgbbuf()
               do 16 k=1,npt
                  if((nbstar(k).lt.26.9).and.(ans(1,k).lt.0.08))then
                     call pgsci(2)
                     call pgpt1(ans(i,k),ans(j,k),17)
                  else
                     call pgsci(1)
                     call pgpt1(ans(i,k),ans(j,k),1)
                  endif
 16            continue
               call pgebuf()
               call pgsci(1)
            endif
 15      continue
 14   continue

      goto 999
 901  write(0,*) "Usage: tdurplot mcmcfile"
      goto 999
 902  write(0,*) "Cannot open ",mcmcfile
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bindata(nbin,npt,pdata,bdatax,bdatay,datamin,datamax,
     .   bmax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nbin,npt,i,bin
      real pdata(npt),bdatax(nbin),bdatay(nbin),datamin,datamax,
     .   binsize,bmax,tsum

C     get size of bin.
      binsize=(datamax-datamin)/real(nbin-1)

C     initalize bdata
      do 20 i=1,nbin
         bdatay(i)=0.
 20   continue

C     loop over data and place data into proper bin.
      do 10 i=1,npt
C        calculate bin number
         bin=int((pdata(i)-datamin)/binsize)+1
         if((bin.gt.0).and.(bin.le.nbin)) bdatay(bin)=bdatay(bin)+1.0
 10   continue

C     get the max value in a bin and assign bin value.
      tsum=0.
      bmax=0.
      do 30 i=1,nbin
         bmax=max(bdatay(i),bmax)
         bdatax(i)=datamin+real(i-1)*binsize!+binsize/2.0 !added -1
         tsum=tsum+bdatay(i)
 30   continue
c      write(6,*) "Sum:",tsum

      return
      end
