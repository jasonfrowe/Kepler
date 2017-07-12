      subroutine fittransitmodel2(npt,nfit,sol,sol2,serr,fvec,iwa,lwa,
     .  wa)
      implicit none
      integer npt,n,nfit,i,info,iwa(nfit),lwa,m
      double precision sol(nfit),serr(nfit,2),fvec(npt),tol,wa(lwa),
     .  sol2(nfit)
      external fcn
      
      n=0
      do 10 i=1,nfit
c        write(0,*) i,serr(i,2)
        if(serr(i,2).ne.0.0)then
            n=n+1
            sol2(n)=sol(i)
        endif
 10   continue
c      write(0,*) "n:",n
      
      tol=1.0d-8

      write(0,503) sol(1),sol(2),sol(10),sol(11),sol(12),sol(13),
     .      sol(14),sol(15),sol(16),sol(17),sol(18),sol(19),sol(20)
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
      write(0,503) sol(1),sol(2),sol(10),sol(11),sol(12),sol(13),
     .      sol(14),sol(15),sol(16),sol(17),sol(18),sol(19),sol(20)
 503  format(28(1PE10.3,1X))


      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine fcn(m,n,x,fvec,iflag)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer m,n,iflag,nmax,nrad,ntheta,nfit,i,j,nplanet,nfrho,
     .  nplanetmax
      parameter(nmax=500000,nfit=101,nplanetmax=10)
      integer dtype(nmax),ntt(nplanetmax)
      double precision x(n),fvec(m),aT(nmax),aM(nmax),aE(nmax),
     .  aIT(nmax),sol(nfit),serr(nfit,2),sol2(nfit),
     .  tobs(nplanetmax,nmax),
     .  omc(nplanetmax,nmax)
      common /Fitting2/ nplanet,dtype,aT,aM,aE,aIT,sol,
     .  serr,ntt,tobs,omc
      
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

c      write(0,503) (sol2(i),i=1,20)
c      write(0,503) (serr(i,2),i=1,20)
c      read(5,*)
c 503  format(28(1PE10.3,1X))

      call transitmodel(nfit,nplanet,nplanetmax,sol2,nmax,m,aT,aIT,ntt,
     .  tobs,omc,fvec,dtype)
c      do 13 i=1,m
c      write(0,*) fvec(i)
c      read(5,*)
c 13   continue
      
      do 12 i=1,m
        fvec(i)=(fvec(i)-aM(i))/aE(i)
 12   continue


      return
      end
      