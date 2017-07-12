CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function fittransitmodel(nfit,sol,serr,npars,
     .  ipars,p,pvary,Dpvary,x,y,tp,dchistop)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,i,j,ipars(nfit),npars,it,itmax,iter,tp(nfit+1)
      double precision sol(nfit),serr(nfit,2),p(nfit+1,nfit),
     .  pvary(nfit),y(nfit+1),x(nfit),cfunk,Dpvary(nfit),ftol,bchi,ochi,
     .  dchi,dchistop
      logical loop
      
      j=0
      do 10 i=1,nfit
        if(serr(i,2).ne.0.0)then
            j=j+1
            ipars(j)=i
        endif
 10   continue
      npars=j
 
      write(0,500) "To Fit: ",(ipars(i),i=1,npars)
 500  format(A8,13(1X,I2))

      do 22 i=1,npars
        p(npars+1,i)=sol(ipars(i))
 22   continue

      loop=.true.
      itmax=25  !takes about 5 iterations to get minimum.
      it=1
      do 30 while(loop)
C       Going to put parameters in order of what I want fitted.
C       Initial guesses go into element nfit of p
        do 11 i=1,npars
            pvary(i)=Dpvary(ipars(i))
 11     continue
      
         do 14 i=1,npars
            do 15 j=1,npars
               p(i,j)=p(npars+1,j) !copy row 11 into the other rows
 15         continue
            p(i,i)=p(npars+1,i)+pvary(i)!increase diagonals 
 14      continue

c         if(it.eq.1) then
         write(0,*) "p:"
         do 12 j=1,npars
            write(0,501) (p(i,j),i=1,npars+1)
 501        format(13(1PE10.3,1X))
 12      continue
c         endif
c         write(0,*)
c         write(0,501) (p(i,i),i=1,npars)

         do 18 i=1,npars+1
            do 19 j=1,npars
               x(j)=p(i,j)  !put each row into a temp array
 19         continue
            y(i)=cfunk(x)   !now we initialize y with cfunk
 18      continue
         if(it.eq.1) ochi=y(npars+1)
         write(0,*) "y:"
         write(0,501) ((y(i)-y(npars+1))/y(npars+1),i=1,npars+1)
         write(0,*) "Chi-Sq:",y(npars+1)
         do 21 i=1,npars
            if ((y(i)-y(npars+1))/y(npars+1).lt.0) pvary(i)=-pvary(i)
 21      continue
         
         ftol=1.0d-16  !set tolerance to something small.
         call amoeba(p,y,nfit+1,nfit,npars,ftol,cfunk,iter)
         write(0,*) 
         
C        Now we sort the chi-squares to find the minimum.
         call rqsort(npars+1,y,tp)!does not destroy arrays.returns index

c         write(6,*) y(tp(1)),iter
c         write(6,*) (p(tp(1),i),i=1,nfit)
C        move best chi-square into 11 position.
         do 20 i=1,npars
            p(npars+1,i)=p(tp(1),i)
 20      continue
 
         write(0,504) "Best fit for iteration # ",it
 504     format(A25,I4)
         bchi=y(tp(1))
         dchi=(bchi-ochi)/ochi
         write(0,503) (p(npars+1,i),i=1,npars),bchi,dchi
         ochi=bchi
 503     format(7(1PE10.3,1X))
         it=it+1
         if(it.gt.itmax) loop=.false.
c         if((abs(dchi).lt.dchistop).and.(abs(dchi).gt.0.0)) loop=.false.
         if(abs(dchi).lt.dchistop) loop=.false.
 30   enddo
 
      bchi=y(tp(1))
      do 31 i=1,npars
         sol(ipars(i))=p(tp(1),i)
         Dpvary(ipars(i))=pvary(i)
 31   continue
      write(0,502) (sol(i),i=1,nfit)
 502  format(13(F7.3,1X))
      write(0,*) "Chi-Sq:",bchi
      fittransitmodel=bchi
 
      return
      end