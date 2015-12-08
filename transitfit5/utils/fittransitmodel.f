CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine fittransitmodel(npta,aT,aM,aE,dtype,nfitm,nplanet,sol,
     .  serr,ia,covar,alpha,dchistop,err)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npta,nfitm,i,it,itmax,npt,nplanet,nfit,ii
      integer ia(nfitm),dtype(npta)
      double precision aT(npta),aM(npta),aE(npta),sol(nfitm),
     .  alpha(nfitm,nfitm),covar(nfitm,nfitm),chisq,alamda,alamhi,
     .  serr(nfitm,2),bchi,ochi,dchi,dchistop,err(nfitm)
      logical loop
      
      nfit=nplanet*10+8
c      write(0,*) "nplanet",nplanet,nfit
      
      do 10 i=1,nfit
        if(serr(i,2).ne.0.0d0)then
            ia(i)=1 !mark parameter to fit
        else
            ia(i)=0 !mark parameter as fixed
        endif
 10   continue
      ia(6)=0 !never fit for dilution

      alamhi=100000.
      alamda=-1.0
      loop=.true.
      itmax=50  !takes about 5 iterations to get minimum.
      it=1
      do 30 while(loop)
        call mrqmin(aT,aM,aE,dtype,npta,sol,ia,nfit,covar,alpha,nfitm,
     .      chisq,alamda)
C        inclination bounds check
         do 11 ii=1,nplanet
            if(sol(10*(ii-1)+8+3).lt.0.0d0)then
                sol(10*(ii-1)+8+3)=abs(sol(10*(ii-1)+8+3))
                ia(10*(ii-1)+8+3)=0
            endif
 11      continue
 
         write(0,504) "Best fit for iteration # ",it
 504     format(A25,I4)
         bchi=chisq
         if(it.gt.1)then
            dchi=(ochi-bchi)!/ochi
         else
            dchi=1.0e30
         endif
         write(0,*) "alamda",alamda,dchi
         write(0,503) sol(1),sol(8),sol(9),sol(10),sol(11),sol(12),
     .      bchi,dchi
         ochi=bchi
 503     format(17(1PE10.3,1X))
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
     .          nfitm,chisq,alamda)
         endif
c         write(0,*) "Chi-squared: ",chisq
 30   enddo
 
 
      do 33 i=1,nfit
        err(i)=0.0
 33   continue
      do 31 i=1,nfit
         err(i)=sqrt(covar(i,i)*bchi/dble(npta-1))
         write(0,*) i,ia(i),err(i)
 31   continue
     
      return
      end
