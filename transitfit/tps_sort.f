      program tps_sort
      implicit none
      integer nmax,i,j,iargc,nunit,npt
      parameter(nmax=600000)
      integer kid1(nmax),kid2
      real*8 kmag,teff,logg,rad,maghi,maglow
      character*80 filename,cline(nmax)
      
      if(iargc().lt.3) goto 901
      call getarg(1,filename)
      call getarg(2,cline(1))
      read(cline,*)maglow
      call getarg(3,cline(1))
      read(cline,*)maghi
      
      nunit=10
      open(unit=nunit,file=filename,status='old',err=902)
      read(nunit,500) cline(1)
 500  format(A80)
      write(6,500) cline(1)
      
      i=1
 12   read(nunit,500,end=13) cline(i)
        read(cline(i),*) kid1(i)
        i=i+1
        goto 12
 13   continue
      npt=i-1
      close(nunit)
      
c      filename="/home/rowe/Kepler/q1/kictpsq1.txt" !read in KIC info
c      filename="/media/Streetsville/rowe/Kepler/tpsq3/kic_10_11.txt"
c      filename="/media/Etobicoke/Kepler/tpsq3/kic_10_11.txt"
      filename="/home/rowe/p555/transitfit/Kepler/kic_0_16.5.txt"
      open(unit=nunit,file=filename,status='old',err=902)
      
      
      do 14 j=1,npt
        i=0
 10     read(nunit,*,end=11) kid2,kmag,teff,logg,rad
            if(kid1(j).eq.kid2)then
                i=1
                goto 11
            endif
            goto 10
 11     continue
        if(i.eq.0)then
            kmag=0.0
            teff=0.0
            rad=0.0
            logg=0.0
        endif
        if((kmag.gt.maglow).and.(kmag.le.maghi))then
            write(6,500) cline(j)
        endif
c        if((teff.lt.4000).and.(teff.gt.10).and.(logg.gt.4.0))then
c            write(6,500) cline(j)
c        endif
        rewind(nunit)
 14   continue
      
      goto 999
 901  write(0,*) "Usage: tps_sort tps_results maglow maghi"
      goto 999
 902  write(0,*) "Cannot open ",filename
      goto 999
 999  end