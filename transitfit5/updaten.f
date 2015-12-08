      program updaten
      implicit none
      integer nunit,nfit,nplanet,iargc,i,flag
      parameter(nfit=108)
      double precision sol(nfit),serr(nfit,2),err(nfit),rhostar,nl(4)
      character*3 titles(nfit)
      character*80 filename,cline
      
      if(iargc().lt.6) goto 901
      call getarg(1,filename)
      nunit=10
      open(unit=nunit,file=filename,status='old',err=902)
      call getfitpars(nunit,nfit,nplanet,sol,serr,err)
      close(nunit)
      
      call getarg(2,cline)
      read(cline,*) rhostar
      sol(1)=rhostar
      serr(1,2)=0.0d0
      do 10 i=1,4
        call getarg(2+i,cline)
        read(cline,*) nl(i)
        sol(i+1)=nl(i)
        serr(i+1,2)=0.0d0
 10   continue

C     Check to see if we need RVs
c      flag=0
c      do 11 i=1,nplanet
c        if(sol(15+10*(i-1)).ne.0.0) flag=1
c
c 11   continue
      flag=1
 
      if(flag.eq.1)then
        serr(7,2)=-1.0
        do 12 i=1,nplanet
            serr(15+10*(i-1),2)=-1.0 !K
c            serr(13+10*(i-1),2)=-1.0 !ecosw
c            serr(14+10*(i-1),2)=-1.0 !esinw
c            write(0,*) 15+10*(i-1)
 12     continue
      endif
      
      do 13 i=1,nplanet
        sol(11+10*(i-1))=0.5

c        write(0,*) 15+10*(i-1)      
 13   continue
      
      call exportfit(nfit,nplanet,sol,serr,err,titles)
      
      goto 999
 901  write(0,*) "Usage: updaten n1.dat rhostar n1 n2 n3 n4"
      write(0,*) " n1.dat - transitfit parameter file"
      write(0,*) " rhostar - mean stellar density g/cm^3"
      write(0,*) " n1-4 - non-linear limb-darkening coefficients"
      goto 999
 902  write(0,*) "Cannot open ",filename
      goto 999
 999  end
