      program mcmccorr
C     Fit correlations in MCMC burn-in for a second-run
      implicit none
      integer nplanet,nmax,iargc,nunit,i,nfit,j,dumi,npt,k,ii
      parameter(nmax=2000000) !maximum planets
      double precision rho(nmax),b(nmax),rprs(nmax),dumr,ans(2),mx,
     .  sx,s,ep,var,ecw(nmax),esw(nmax)
      character*80 filename
      
      nunit=10
      call getarg(1,filename)
      
      open(unit=nunit,file=filename,status='old',err=901)
      read(nunit,*) nfit
      nplanet=(nfit-8)/10
      
      do 10 i=1,nplanet
        j=1
        ii=0
c        write(0,*) 4+9+10*(i-1)
 11     read(nunit,*,end=12) dumr,dumi,dumi,rho(j),
     .    (dumr,k=1,9+10*(i-1)),b(j),rprs(j),ecw(j),esw(j)
            if((rho(j).eq.rho(j-1)).and.(b(j).eq.b(j-1)).and.
     .          (rprs(j).eq.rprs(j-1)).and.(j.gt.1)) goto 11
            ii=ii+1
            if(ii.ge.10)then
                j=j+1
                ii=0
            endif
            if(j.gt.nmax) goto 902
        goto 11
 12     continue
        npt=j-1
                
        write(0,*) "npt: ",npt
        call lfit(npt,rho,b,ans)
        mx=0.0d0
        do 13 j=1,npt
            mx=mx+b(j)-(ans(1)+ans(2)*rho(j))
 13     continue
        mx=mx/dble(npt)
        do 14 j=1,npt
            s=(b(j)-(ans(1)+ans(2)*rho(j)))-mx
            ep=ep+s
            var=var+s*s
 14     continue
        var=(var-ep**2/npt)/(npt-1)
        sx=sqrt(var)
        write(0,500) ans(2),sx,ans(1)
 500    format(3(1X,1PE17.10))
 
        call lfit(npt,rho,rprs,ans)
        do 15 j=1,npt
            mx=mx+rprs(j)-(ans(1)+ans(2)*rho(j))
 15     continue
        mx=mx/dble(npt)
        do 16 j=1,npt
            s=(rprs(j)-(ans(1)+ans(2)*rho(j)))-mx
            ep=ep+s
            var=var+s*s
 16     continue
        var=(var-ep**2/npt)/(npt-1)
        sx=sqrt(var)
        write(0,500) ans(2),sx,ans(1)

c        call lfit(npt,rho,ecw,ans)
c        do 17 j=1,npt
c            mx=mx+ecw(j)-(ans(1)+ans(2)*rho(j))
c 17     continue
c        mx=mx/dble(npt)
c        do 18 j=1,npt
c            s=(ecw(j)-(ans(1)+ans(2)*rho(j)))-mx
c            ep=ep+s
c            var=var+s*s
c 18     continue
c        var=(var-ep**2/npt)/(npt-1)
c        sx=sqrt(var)
c        write(0,500) ans(2),sx,ans(1)
 
c        call lfit(npt,rho,esw,ans)
c        do 19 j=1,npt
c            mx=mx+esw(j)-(ans(1)+ans(2)*rho(j))
c 19     continue
c        mx=mx/dble(npt)
c        do 20 j=1,npt
c            s=(esw(j)-(ans(1)+ans(2)*rho(j)))-mx
c            ep=ep+s
c            var=var+s*s
c 20     continue
c        var=(var-ep**2/npt)/(npt-1)
c        sx=sqrt(var)
c        write(0,500) ans(2),sx,ans(1)
 
        rewind(nunit)
        read(nunit,*)dumi
 10   continue
      
      close(nunit)
      
      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999    
 902  write(0,*) "Segfault detected.  Increase nmax"
      goto 999
 999  end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lfit(npt,x,y,ans)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i
      double precision x(npt),y(npt),ans(2),Sx,Sy,Sxx,Sxy,d,S
      
      S=dble(npt)
      
      Sx=0.0d0
      do 10 i=1,npt
        Sx=Sx+x(i)
 10   continue
      
      Sy=0.0d0
      do 11 i=1,npt
        Sy=Sy+y(i)
 11   continue     
      
      Sxx=0.0d0
      do 12 i=1,npt
        Sxx=Sxx+x(i)*x(i)
 12   continue
      
      Sxy=0.0d0
      do 13 i=1,npt
        Sxy=Sxy+x(i)*y(i)
 13   continue
      
      d=S*Sxx-Sx*Sx
      
      ans(1)=(Sxx*Sy-Sx*Sxy)/d
      ans(2)=(S*Sxy-Sx*Sy)/d
      
      return
      end
      
