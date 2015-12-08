      program notesmerge
      implicit none
      integer nmax,i,nunit,nmaster,j,np,k,nfit
      parameter(nmax=200000,np=8,nfit=18)
      integer kepid(nmax),fitn(nmax),class(nmax,np),nkepid,nfitn,nclass,
     .  nfound(nmax),tceID(nmax),ntceID,ntcemin,ntce
      double precision sol(nmax,nfit),err(nmax,nfit),rchi(nmax),
     .  lchi(nmax),per,epo,messes(nmax),mes(nmax),nmesses,nmesmepo,rs,
     .  dur,depthest,kepmag,perdif,epodif,mindiff,diff
      character*80 filein(np),filem,txt1,txt2,txt3,filet

C     Useful AWK line to make sure lines have 3 columns
C     cat notes_dh.txt | awk '{if (NF>2) print $0; else print $1,$2,-1}' > test.dat
      
      filem="tps_stats.cut.dat" !master list
      filet="tps_q1q10.pt1.dat" !tps input
      filein(1)="notes_CJB.all.txt" !chris burke
      filein(2)="notes_dh.txt" !dan huber
      filein(3)="notes_fm.all.txt" !fergal mullally
      filein(4)="notes_jessie.txt" !jessie christiansen
      filein(5)="notes_mh.txt" !mathieu havel
      filein(6)="notes_set.all.txt" !susan mullally
      filein(7)="notes_eq.txt" !elisa quintana
      filein(8)="notes_jr.txt" !jason rowe
      
      do 10 i=1,nmax
        nfound(i)=0
        do 9 j=1,np
            class(i,j)=-1
  9     continue
 10   continue
 
      nunit=10
      open(unit=nunit,file=filem,status='old',err=901)
      i=1
 11   read(nunit,*,end=12) kepid(i),fitn(i),(sol(i,j),j=1,nfit),
     .  (err(i,j),j=1,nfit),rchi(i),lchi(i)
        if(isnan(rchi(i))) rchi(i)=-1.0
        if(isnan(lchi(i))) rchi(i)=-1.0
        do 18 j=1,nfit
            if(isnan(sol(i,j))) sol(i,j)=-1.0
            if(isnan(err(i,j))) err(i,j)=-1.0
 18     continue
        if(sol(i,6).lt.0.0) sol(i,6)=abs(sol(i,6))
        if(sol(i,6).gt.90.0) sol(i,6)=sol(i,6)-90*int(sol(i,6)/90.0)
        err(i,4)=min(sol(i,4),err(i,4))
        err(i,6)=min(sol(i,6),err(i,6))
        if(sol(i,7).lt.0.0) sol(i,7)=sol(i,7)+sol(i,4)*(abs(int(
     .      sol(i,7)/sol(i,4)))+1)
        err(i,7)=min(sol(i,7),err(i,7))
        i=i+1
      goto 11
 12   continue
      nmaster=i-1
      close(nunit)
      
      open(unit=nunit,file=filet,status='old',err=904)
 21   read(nunit,*,end=19) nkepid,per,nmesses,nmesmepo,epo,rs,dur,
     .  depthest,ntceID,kepmag
      
        mindiff=99.9e30
        do 20 i=1,nmaster
            if(nkepid.eq.kepid(i))then
                perdif=abs(sol(i,5)-per)
                epodif=abs(sol(i,7)+54900-epo)
                diff=perdif+epodif
                if(diff.lt.mindiff)then
                    ntce=i
                    ntcemin=ntceID
                    mindiff=diff
                endif
            endif
 20     continue
        tceid(ntce)=ntcemin
      goto 21
 19   continue
      close(nunit)
      
      do 13 i=1,np
        k=1
        open(unit=nunit,file=filein(i),status='old',err=902)
 14     read(nunit,*,end=15,err=903) txt1,txt2,txt3
c            write(6,*) txt1,txt2,txt3
c            write(6,*) 'bump'
            read(txt1,*) nkepid
            read(txt2,*) nfitn
            if(txt3.eq.' ')then
                nclass=-1
            else
                read(txt3,*,err=903) nclass
            endif
            k=k+1 !count line number
            do 16 j=1,nmaster
                if((kepid(j).eq.nkepid).and.(fitn(j).eq.nfitn))then
                    class(j,i)=nclass
                    nfound(j)=1
                endif
 16         continue
        goto 14
 15     continue
        close(nunit)
 13   continue
      
      do 17 i=1,nmaster
        if(nfound(i).eq.1)then
            write(6,500) kepid(i),fitn(i),(class(i,j),j=1,np),
     .          sol(i,1),(sol(i,j),j=3,7),(err(i,j),j=4,7),rchi(i),
     .          lchi(i),tceID(i)
        endif
 500    format(I9,1X,I2,8(1X,I2),3(1X,F7.3),1X,F12.6,1X,F7.2,1X,F13.6,
     .      1X,F7.3,1X,F12.6,1X,F7.2,1X,F12.6,2(1X,1PE13.6),1X,I6)
 17   continue
      
      
      goto 999
 901  write(0,*) "Cannot open ",filem
      goto 999
 902  write(0,*) "Cannot open ",filein(i)
      goto 999
 903  write(0,*) "Error: line ",k
      write(0,*) txt1,txt2,txt3
      write(0,*) "file: ",filein(i)
      goto 999
 904  write(0,*) "Cannot open ",filet
      goto 999
 999  end