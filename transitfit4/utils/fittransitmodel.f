CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine fittransitmodel(npta,aT,aM,aE,dtype,nfit,sol,serr,ia,
     .  covar,alpha,dchistop,err)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npta,nfit,i,it,itmax
      integer ia(nfit),dtype(npta)
      double precision aT(npta),aM(npta),aE(npta),sol(nfit),
     .  alpha(nfit,nfit),covar(nfit,nfit),chisq,alamda,alamhi,
     .  serr(nfit,2),bchi,ochi,dchi,dchistop,err(nfit)
      logical loop
      
      do 10 i=1,nfit
        if(serr(i,2).ne.0.0d0)then
            ia(i)=1 !mark parameter to fit
        else
            ia(i)=0 !mark parameter as fixed
        endif
 10   continue
      ia(18)=0 !never fit for dilution

      alamhi=100000.
      alamda=-1.0
      loop=.true.
      itmax=50  !takes about 5 iterations to get minimum.
      it=1
      do 30 while(loop)
        call mrqmin(aT,aM,aE,dtype,npta,sol,ia,nfit,covar,alpha,nfit,
     .      chisq,alamda)
C        inclination bounds check
         if(sol(3).lt.0.0d0)then
            sol(3)=abs(sol(3))
            ia(3)=0
         endif

         write(0,504) "Best fit for iteration # ",it
 504     format(A25,I4)
         bchi=chisq
         if(it.gt.1)then
            dchi=(ochi-bchi)!/ochi
         else
            dchi=1.0e30
         endif
         write(0,*) "alamda",alamda,dchi
         write(0,503) (sol(i),i=1,17),bchi,dchi
         ochi=bchi
 503     format(19(1PE10.3,1X))
         it=it+1
c         if(it.gt.itmax) loop=.false.
c         if((abs(dchi).lt.dchistop).and.(abs(dchi).gt.0.0)) loop=.false.
c         if(abs(dchi).lt.dchistop) then
          if((it.gt.itmax).or.
     .    ((dchi.lt.dchistop).and.(dchi.gt.0.0d0)).or.
     .    (alamda.gt.alamhi))then
            loop=.false.
            alamda=0.0
            call mrqmin(aT,aM,aE,dtype,npta,sol,ia,nfit,covar,alpha,
     .          nfit,chisq,alamda)
         endif
 30   enddo
 
      do 33 i=1,nfit
        err(i)=0.0
 33   continue
      do 31 i=1,nfit
         err(i)=sqrt(abs(covar(i,i))*bchi/dble(npta-1))
 31   continue
     
      return
      end
