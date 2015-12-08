      program tpsclean
      implicit none
      integer nmax,nunit,i,j,flag,nunit2,npt,nunit3,nkoi,nEB,nVar,
     .  nunit4,nKIC,nunit5
      parameter(nmax=3000000)
      integer kid(nmax),koi(nmax),eb(nmax),var(nmax),kic(nmax),
     .  koinum,koiID(nmax)
      double precision mes(nmax),ses(nmax),per(nmax),epo(nmax),tce(5),
     .  dumr,kmag(nmax),mag(nmax)
      character*80 filename,filenameEB,filenameKID,filenameVAR,
     .  filenameKIC
      
      nunit=10
      filename="tps_q1q8.dat"
      nunit2=11
      filenameEB="eb.list"
      nunit3=12
      filenameKID="KOI.list"
      nunit4=13
      filenameVAR="var.list"
      nunit5=14
      filenameKIC="kic_colour_16.5.txt"
      
      
      open(unit=nunit3,file=filenameKID,status='old',err=904)      
      i=1
 19   read(nunit3,*,end=20) koi(i),koiID(i)
        i=i+1
      goto 19
 20   continue
      nkoi=i-1
      close(nunit2)

      open(unit=nunit2,file=filenameEB,status='old',err=902)
      i=1
 24   read(nunit2,*,end=23) eb(i)
        i=i+1
      goto 24
 23   continue
      nEB=i-1
      close(nunit2)
      
      open(unit=nunit4,file=filenameVAR,status='old',err=905)
      i=1
 29   read(nunit4,*,end=28) var(i)
        i=i+1
      goto 29
 28   continue
      nVAR=i-1
      close(nunit4)

      open(unit=nunit5,file=filenameKIC,status='old',err=906)
      i=1
 30   read(nunit5,*,end=31) kic(i),kmag(i)
        i=i+1
      goto 30
 31   continue
      nKIC=i-1
      close(nunit5)
      
      open(unit=nunit,file=filename,status='old',err=901)
      
      i=1      
      read(nunit,*) (tce(j),j=1,5)
      kid(i)=int(tce(1))
      per(i)=tce(2)
      ses(i)=tce(3)
      mes(i)=tce(4)
      epo(i)=tce(5)
      
 11   read(nunit,*,end=12) (tce(j),j=1,5)
        flag=0
        if(tce(4).lt.7.0) flag=1
        if(tce(3).lt.sqrt(2.0)) flag=1
        if(flag.eq.0)then
            do 10 j=1,i  !scan through existing KOIs
                if(kid(j).eq.int(tce(1)))then !does KOI already exist?
                    flag=1
                    if(tce(4).gt.mes(j))then !if better, take it
                        kid(j)=int(tce(1))
                        per(j)=tce(2)
                        ses(j)=tce(3)
                        mes(j)=tce(4)
                        epo(j)=tce(5)
                    endif
                endif
 10         continue
        endif
        if(flag.eq.0)then !if it's a new KOI...
            i=i+1
            kid(i)=int(tce(1))
            per(i)=tce(2)
            ses(i)=tce(3)
            mes(i)=tce(4)
            epo(i)=tce(5)
c            write(6,500) kid(i),per(i),ses(i),mes(i),epo(i)
        endif
      goto 11
 12   continue
      npt=i
      
      do 25 i=1,npt
        flag=1
        do 26 j=1,nEB
            if(eb(j).eq.kid(i)) then
               flag=1
               write(6,*) kid(i),per(i)
            endif
 26     continue
        do 27 j=1,nKOI
            if(koi(j).eq.kid(i)) then
                flag=0
                koinum=koiID(j)               
            endif
 27     continue
c        do 32 j=1,nVAR
c            if(var(j).eq.kid(i)) flag=1
c 32     continue
        do 33 j=1,nKIC
            if(KIC(j).eq.kid(i)) mag(i)=kmag(j)
 33     continue
ccc        if(flag.eq.0) write(6,500) kid(i),per(i),ses(i),mes(i),epo(i),
ccc     .      mag(i),koinum
 25   continue
 500  format(I8,5(1X,F12.6),1X,I4)
      
      close(nunit)
      
      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 902  write(0,*) "Cannot open ",filenameEB
      goto 999     
 903  write(0,*) "error on line",j
      goto 999
 904  write(0,*) "Cannot open ",filenameKID
      goto 999
 905  write(0,*) "Cannot open ",filenameVAR
      goto 999
 906  write(0,*) "Cannot open ",filenameKIC
      goto 999    
 999  end
