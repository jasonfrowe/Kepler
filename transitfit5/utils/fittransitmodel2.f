      subroutine fittransitmodel2(npt,nfit,sol,sol2,serr,fvec,iwa,lwa,
     .  wa)
      implicit none
      integer npt,n,nfit,i,info,iwa(nfit),lwa,m
      double precision sol(nfit),serr(nfit,2),fvec(npt),tol,wa(lwa),
     .  sol2(nfit)
      external fcn
      
      n=0
      do 10 i=1,nfit
        if(serr(i,2).ne.0.0)then
            n=n+1
            sol2(n)=sol(i)
        endif
 10   continue
      
      tol=1.0d-8

      write(0,503) sol(1),sol(8),sol(9),sol(10),sol(11),sol(12),
     .      sol(13),sol(14),sol(15),sol(16)      
c      write(0,503) (sol(i),i=1,26)
      
      m=npt
c      call lmdif1(fcn,m,n,x,fvec,tol,info,iwa,wa,lwa)
      call lmdif1(fcn,m,n,sol2,fvec,tol,info,iwa,wa,lwa)
      write(0,*) "info: ",info

      n=0
      do 11 i=1,nfit
        if(serr(i,2).ne.0.0)then
            n=n+1
            sol(i)=sol2(n)
        endif
 11   continue

c      write(0,503) (sol(i),i=1,26)
      write(0,503) sol(1),sol(8),sol(9),sol(10),sol(11),sol(12),
     .      sol(13),sol(14),sol(15),sol(16) 
 503  format(28(1PE10.3,1X))


      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine fcn(m,n,x,fvec,iflag)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer m,n,iflag,nmax,nrad,ntheta,nfit,i,j,nplanet,nfrho,
     .  nplanetmax
      parameter(nmax=2000000,nfit=108,nplanetmax=10)
      integer dtype(nmax),ntt(nplanetmax)
      double precision x(n),fvec(m),aT(nmax),aM(nmax),aE(nmax),
     .  aIT(nmax),sol(nfit),serr(nfit,2),sol2(nfit),rhoin(9),rhoierr(9),
     .  y,yy,yp,rhoi,drho,dsig,chifac,tobs(nplanetmax,nmax),
     .  omc(nplanetmax,nmax),c1,c2,c3,c4
      common /Fitting2/ nfrho,nplanet,dtype,aT,aM,aE,aIT,sol,
     .  serr,rhoi,rhoierr,chifac,ntt,tobs,omc
      data rhoin/-1.0d2,-3.0d0,-2.0d0,-1.0d0,0.0d0,1.0d0,2.0d0,3.0d0,
     .  1.0d2/  
      
      do 10 i=1,nfit
        sol2(i)=sol(i)
 10   continue

      j=0
      do 11 i=1,nfit
        if(serr(i,2).ne.0.0)then
            j=j+1
            sol2(i)=x(j)
        endif
        if(j.gt.n)write(0,*) "whoops.. j>n"
 11   continue

C     Check that limb-darkening is valid
      c1=sol2(2)      !limb-darkening
      c2=sol2(3)
      c3=sol2(4)
      c4=sol2(5)
      if((c1.eq.0.0).and.(c2.eq.0.0))then
        if((c3.lt.0.0).or.(c3.gt.1.0).or.(c4.lt.0.0).or.(c4.gt.1.0))then
            fvec=9.9e30
            return
         endif
      elseif((c3.eq.0.0).and.(c4.eq.0.0))then
         if((c1+c2.gt.1.0).or.(c1.lt.0).or.(c1+2.0d0*c2.lt.0))then
            fvec=9.9e30
            return
         endif
      endif


c      call transitmodel(nfit,nplanet,sol2,m,aT,aIT,fvec,dtype)
      call transitmodel(nfit,nplanet,nplanetmax,sol2,nmax,m,aT,aIT,ntt,
     .  tobs,omc,fvec,dtype)
    
      if(nfrho.eq.0)then
        y=0.0d0
        yy=0.0d0
        do 13 i=1,m
            y=y+(fvec(i)-aM(i))/aE(i)
            yy=yy+y*y
 13     continue
        drho=1.0d3*sol2(1)-rhoi
        call getrhosig(rhoierr,rhoin,9,drho,dsig)
c        write(0,*) 1.0d3*sol2(1),rhoi,dsig
        yp=sqrt( (chifac*yy+dsig*dsig)/yy )
      else
        yp=1.0d0
      endif
c      write(0,*) y,yy,yp
      
      do 12 i=1,m
        fvec(i)=(fvec(i)-aM(i))/aE(i)*yp
c        write(6,*) i,fvec(i)
 12   continue
c      read(5,*) 

      return
      end
      
