CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE mrqmin(x,y,sig,dtype,ndata,a,ia,ma,covar,alpha,nca, 
     *     chisq,alamda) 
ccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER ma,nca,ndata,ia(ma),MMAX,dtype(ndata)
      REAL*8 alamda,chisq,a(ma),alpha(nca,nca),covar(nca,nca), 
     *     sig(ndata),x(ndata),y(ndata)
      PARAMETER (MMAX=1000) 
      INTEGER j,k,l,mfit 
      REAL*8 ochisq,atry(MMAX),beta(MMAX),da(MMAX) 
      SAVE ochisq,atry,beta,da,mfit 
      if(alamda.lt.0.)then 
         mfit=0 
         do 11 j=1,ma 
            if (ia(j).ne.0) mfit=mfit+1 
 11      enddo  
         alamda=0.001 
         call mrqcof(x,y,sig,dtype,ndata,a,ia,ma,alpha,beta,nca,chisq)
         ochisq=chisq 
         do 12 j=1,ma 
            atry(j)=a(j) 
 12      enddo 
      endif 
      do 14 j=1,mfit 
         do 13 k=1,mfit
            covar(j,k)=alpha(j,k) 
 13      enddo 
         covar(j,j)=alpha(j,j)*(1.+alamda) 
         da(j)=beta(j) 
 14   enddo 
      call gaussj(covar,mfit,nca,da,1,1)  
      if(alamda.eq.0.)then 
         call covsrt(covar,nca,ma,ia,mfit) 
         call covsrt(alpha,nca,ma,ia,mfit) 
         return 
      endif 
      j=0 
      do 15 l=1,ma 
         if(ia(l).ne.0) then 
            j=j+1 
            atry(l)=a(l)+da(j) 
         endif 
 15   enddo 
      call mrqcof(x,y,sig,dtype,ndata,atry,ia,ma,covar,da,nca,chisq)
      if(chisq.lt.ochisq)then 
         alamda=0.1*alamda 
         ochisq=chisq 
         do 17 j=1,mfit 
            do 16 k=1,mfit 
               alpha(j,k)=covar(j,k) 
 16         enddo 
            beta(j)=da(j) 
 17      enddo 
         do 18 l=1,ma 
            a(l)=atry(l) 
 18      enddo 
      else 
         alamda=10.*alamda 
         chisq=ochisq 
      endif 
      return 
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE mrqcof(x,y,sig,dtype,ndata,a,ia,ma,alpha,beta,nalp, 
     *     chisq) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER ma,nalp,ndata,ia(ma),MMAX,dtype(ndata)
      REAL*8 chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata), 
     *     y(ndata)
      PARAMETER (MMAX=1000) 
      INTEGER mfit,i,j,k,l,m 
      REAL*8 dy,sig2i,wt,ymod,dyda(MMAX) 
      mfit=0 
      do 11 j=1,ma 
         if (ia(j).ne.0) mfit=mfit+1 
 11   enddo 
      do 13 j=1,mfit 
         do 12 k=1,j
            alpha(j,k)=0. 
 12      enddo 
         beta(j)=0. 
 13   enddo 
      chisq=0. 
      do 16 i=1,ndata 
c         write(0,*) "input",x(i),a,ymod,ma,i
         call funcs(x(i),a,ymod,dyda,ma,i,dtype(i)) 
         sig2i=1./(sig(i)*sig(i)) 
         dy=y(i)-ymod
c         write(6,*) "dy:",(dyda(j),j=1,ma)
         j=0 
         do 15 l=1,ma 
            if(ia(l).ne.0) then 
               j=j+1 
               wt=dyda(l)*sig2i 
               k=0 
               do 14 m=1,l 
                  if(ia(m).ne.0) then 
                     k=k+1 
                     alpha(j,k)=alpha(j,k)+wt*dyda(m) 
                  endif 
 14            enddo 
               beta(j)=beta(j)+dy*wt 
            endif 
 15      enddo 
         chisq=chisq+dy*dy*sig2i  
 16   enddo
      do 18 j=2,mfit 
         do 17 k=1,j-1 
            alpha(k,j)=alpha(j,k) 
 17      enddo
 18   enddo
      return 
      END 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE funcs(x,a,y,dyda,na,n,dtype) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER na,n,nplanet,dtype
      REAL*8 x,y,a(na),dyda(na)
      integer nfit,npt,nmax,i,j,nplanetmax
      parameter(nmax=2000000,nfit=108,nplanetmax=10)
      integer ipars(nfit),ntt(nplanetmax)
      real*8 exptime,itime(nmax),sol2(nfit),vary(nfit),
     .  y1,y2,s0,s1,s2,tobs(nplanetmax,nmax),omc(nplanetmax,nmax)
      common /Fitting/ npt,nplanet,itime,ntt,tobs,omc
c      data vary/2.0d-3,1.0d-4,1.0d-3,1.0d-4,1.0d-1,5.0d-5,1.0d-2,1.0d-2,
c     .  1.0,0.01,
c     .  0.01,0.01,0.01,0.01/
      vary(1)=0.1 !RHO
      vary(2)=0.01 !NL1
      vary(3)=0.01 !NL2
      vary(4)=0.01 !NL3
      vary(5)=0.01 !NL4
      vary(6)=0.01 !DIL
      vary(7)=0.01 !VOF
      vary(8)=5.0d-5 !ZPT
      do 10 i=1,nplanet
        vary(9+(i-1)*10)=2.0d-3 !EPO
        vary(10+(i-1)*10)=1.0d-4 !PER
        vary(11+(i-1)*10)=1.0d-2 !BBB
        vary(12+(i-1)*10)=1.0d-4 !RDR
        vary(13+(i-1)*10)=1.0d-2 !ECW
        vary(14+(i-1)*10)=1.0d-2 !ESW
        vary(15+(i-1)*10)=100.0   !KRV
        vary(16+(i-1)*10)=10.0   !TED
        vary(17+(i-1)*10)=10.0   !ELL
        vary(18+(i-1)*10)=10.0   !ALB
 10   continue

      exptime=itime(n)
c      call transitmodel(na,nplanet,a,1,x,exptime,y,dtype)
      call transitmodel(na,nplanet,nplanetmax,a,nmax,1,x,exptime,ntt,
     .  tobs,omc,y,dtype)

      do 12 i=1,na
        do 13 j=1,na
            sol2(j)=a(j)
 13     continue
        s0=a(i)
        s1=s0+vary(i)
        s2=s0-vary(i)
        if(ipars(i)-10*((ipars(i)-9)/10).eq.12)then
            s1=abs(s1)
            s2=abs(s2)
        endif
        if(ipars(i)-10*((ipars(i)-9)/10).eq.11)then
            if(s1.lt.0.0d0) s1=abs(s1)    
            if(s1.gt.1.0d0) s1=1.0d0
            if(s2.lt.0.0d0) s2=abs(s2)    
            if(s2.gt.1.0d0) s2=1.0d0
c            if(s1.eq.s2)then
c                s1=1.0
c                s2=0.0
c            endif
        endif
        sol2(i)=s1
        call transitmodel(na,nplanet,nplanetmax,sol2,nmax,1,x,exptime,
     .      ntt,tobs,omc,y1,dtype)
        sol2(i)=s2
c        call transitmodel(na,nplanet,sol2,1,x,exptime,y2,dtype)
        call transitmodel(na,nplanet,nplanetmax,sol2,nmax,1,x,exptime,
     .      ntt,tobs,omc,y2,dtype) 
        dyda(i)=(y1-y2)/(s1-s2)
c        write(6,*) x,y,dyda(i)
c        write(6,*) i,a(i),dyda(i)
c        write(6,*) x,y1,y2
 12   continue
c      read(5,*)

      return 
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE gaussj(a,n,np,b,m,mp) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER m,mp,n,np,NMAX 
      REAL*8 a(np,np),b(np,mp) 
      PARAMETER (NMAX=1000) 
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX), 
     *     ipiv(NMAX) 
      REAL*8 big,dum,pivinv 
      do 11 j=1,n 
         ipiv(j)=0 
 11   enddo 
      do 22 i=1,n 
         big=0. 
         do 13 j=1,n 
            if(ipiv(j).ne.1)then 
               do 12 k=1,n 
                  if (ipiv(k).eq.0) then 
                     if (abs(a(j,k)).ge.big)then 
                        big=abs(a(j,k)) 
                        irow=j 
                        icol=k 
                     endif 
                  else if (ipiv(k).gt.1) then 
                     write(0,*) "sing:",k
                     return
c                     stop!pause 'singular matrix in gaussj'
                  endif 
 12            enddo 
            endif 
 13      enddo 
         ipiv(icol)=ipiv(icol)+1
         if (irow.ne.icol) then 
            do 14 l=1,n 
               dum=a(irow,l) 
               a(irow,l)=a(icol,l) 
               a(icol,l)=dum 
 14         enddo
            do 15 l=1,m 
               dum=b(irow,l) 
               b(irow,l)=b(icol,l) 
               b(icol,l)=dum 
 15         enddo 
         endif 
         indxr(i)=irow 
         indxc(i)=icol 
         if (a(icol,icol).eq.0.) then
            write(0,*) "sing2:",icol
            return
c            stop!pause 'singular matrix in gaussj'
         endif 
         pivinv=1./a(icol,icol) 
         a(icol,icol)=1. 
         do 16 l=1,n 
            a(icol,l)=a(icol,l)*pivinv 
 16      enddo 
         do 17 l=1,m 
            b(icol,l)=b(icol,l)*pivinv 
 17      enddo 
         do 21 ll=1,n 
            if(ll.ne.icol)then 
               dum=a(ll,icol) 
               a(ll,icol)=0. 
               do 18 l=1,n 
                  a(ll,l)=a(ll,l)-a(icol,l)*dum 
 18            enddo
               do 19 l=1,m 
                  b(ll,l)=b(ll,l)-b(icol,l)*dum 
 19            enddo 
            endif 
 21      enddo
 22   enddo 
      do 24 l=n,1,-1 
         if(indxr(l).ne.indxc(l))then 
            do 23 k=1,n 
               dum=a(k,indxr(l)) 
               a(k,indxr(l))=a(k,indxc(l)) 
               a(k,indxc(l))=dum 
 23         enddo 
         endif 
 24   enddo 
      return 
      END 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE covsrt(covar,npc,ma,ia,mfit) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER ma,mfit,npc,ia(ma) 
      REAL*8 covar(npc,npc) 
      INTEGER i,j,k 
      REAL*8 swap
      do 12 i=mfit+1,ma 
         do 11 j=1,i 
            covar(i,j)=0. 
            covar(j,i)=0. 
 11      enddo 
 12   enddo 
      k=mfit 
      do 15 j=ma,1,-1 
         if(ia(j).ne.0)then 
            do 13 i=1,ma 
               swap=covar(i,k) 
               covar(i,k)=covar(i,j) 
               covar(i,j)=swap 
 13         enddo 
            do 14 i=1,ma 
               swap=covar(k,i) 
               covar(k,i)=covar(j,i) 
               covar(j,i)=swap 
 14         enddo 
            k=k-1 
         endif 
 15   enddo 
      return 
      END
