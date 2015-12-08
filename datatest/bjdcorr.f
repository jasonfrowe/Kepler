C     Tags on BJD corrections.
      program bjdcorrection
      implicit none
      integer kid,gap
      double precision mjd,flux,ferr,bjdcorr,median,time
      character*80 filename,dir1,dir2,dir3,ext1,ext2,ext3,listfile,
     .  mjdfile,medfile,corfile,newfile,ext4,gapfile
      
      dir1="q10phot_mjd/"               !The MJD time-series files
      dir2="q10phot_mjd/"!"q7phot_bjd/" !where the BJD Corrs are
      dir3="q10phot/"                   !BJD time series is stored here
      ext1=".dat"
      ext2=".bjd"
      ext3=".dat"
      ext4=".gap"
      listfile="q10phot_mjd.list"
      
c      medfile=dir1(1:11)//"median.dat"
c      open(unit=13,file=medfile,status='old',err=904)
      open(unit=10,file=listfile,status='old',err=901)
 10   read(10,*,end=11) filename
c      read(13,*,end=11) kid,median
        write(6,*) filename
        mjdfile=dir1(1:12)//filename(1:11)//ext1(1:4)
        corfile=dir2(1:12)//filename(1:11)//ext2(1:4)
        newfile=dir3(1:8)//filename(1:11)//ext3(1:4)
        gapfile=dir1(1:12)//filename(1:11)//ext4(1:4)
        open(unit=11,file=mjdfile,status='old',err=15)
        goto 16
 15         mjdfile=dir1(1:12)//filename(1:12)//ext1(1:4)
            corfile=dir2(1:12)//filename(1:12)//ext2(1:4)
            newfile=dir3(1:8)//filename(1:12)//ext3(1:4)
            gapfile=dir1(1:12)//filename(1:12)//ext4(1:4)
            open(unit=11,file=mjdfile,status='old',err=902)
 16     continue
        open(unit=12,file=corfile,status='old',err=903)
        open(unit=14,file=newfile)
        open(unit=15,file=gapfile,status='old',err=905)
 12         read(11,*,end=13) mjd,flux,ferr
            read(12,*,end=13) bjdcorr
            read(15,*,end=13) gap
                if((flux.eq.0.0).and.(ferr.eq.0.0)) goto 12
                flux=flux-1.0d0
                time=mjd+bjdcorr              
                if(gap.eq.0)write(14,501) time,flux,ferr
c                write(14,501) time,flux,ferr
 501            format(F14.8,1X,F11.8,1X,F11.8)
            goto 12
 13         continue
c            write(6,*) "done..",kid
c            read(5,*)
        close(11)
        close(12)
        close(14)
        close(15)
        goto 10
 11   continue
      close(10)
c      close(13)
      
      goto 999
 901  write(0,*) "Cannot open ",listfile,"err901"
      goto 999
 902  write(0,*) "Cannot open ",mjdfile,"err902"
      goto 999
 903  write(0,*) "Cannot open ",corfile,"err903"
      goto 999
 904  write(0,*) "Cannot open ",medfile,"err904"
      goto 999
 905  write(0,*) "Cannot open ",gapfile,"err905"
      goto 999
 999  end