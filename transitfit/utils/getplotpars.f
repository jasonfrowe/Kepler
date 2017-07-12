CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getplotpars(nunit,filename,ntype,fitparsfile,id,id2,
     .  kID,Kmag,Teff,logg,rad,eoff,escale,xwidth,bins)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,ntype,id,kID,Teff,id2,bins
      double precision Kmag,logg,rad,eoff,escale,xwidth
      character*80 filename,fitparsfile,cline
      
      read(nunit,*) filename  !1
      read(nunit,*) ntype     !2
      read(nunit,*) fitparsfile !3
      read(nunit,500) cline 
 500  format(A80)
      read(cline,*,err=10,end=10) id,id2 !4
      goto 11
 10     read(cline,*) id
        id2=1
 11   continue
      read(nunit,*) kID  !5
      read(nunit,*) Kmag !6
      read(nunit,*) Teff !7
      read(nunit,*) logg !8
      read(nunit,*) rad  !9
      read(nunit,*) eoff !10
      read(nunit,*) escale  !11
      read(nunit,*) xwidth !12
      bins=0
      read(nunit,*,end=12) bins !13
      write(6,*)"BINS:",bins
 12   continue
      
      return
      end     