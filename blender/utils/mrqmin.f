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
         call funcs(x(i),dtype(i),a,ymod,dyda,ma,i) 
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
      SUBROUTINE funcs(x,dtype,a,y,dyda,na,n) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER na,n,dtype
      REAL*8 x,y,a(na),dyda(na)
      integer nfit,npt,nmax,i
      parameter(nmax=600000,nfit=18)
      integer ipars(nfit)
      real*8 sol(nfit),exptime,sol2(nfit),itime(nmax),dfridr,
     .  pvary(nfit),err,vary(nfit),y1,y2,inclmin,tin,s0,s1,s2
      common /Fitting/ ipars,npt,itime,sol,pvary,inclmin
      data vary/0.1,0.1,0.01,0.01,1.0e-5,0.01,0.002,5.0e-5,0.05,
     .  0.1,0.1,0.1,0.1,0.05,0.1,20.0,1.0,0.01/  


      do 10 i=1,nfit
         sol2(i)=sol(i)
 10   continue
      do 11 i=1,na
         sol2(ipars(i))=a(i)
 11   continue
      exptime=itime(n)
      call transitmodel(1,x,exptime,dtype,y,nfit,sol2)

      do 12 i=1,na
        s0=sol2(ipars(i))
c        if(pvary(i).gt.0.0)then
c            s1=s0+pvary(i)
c            s2=s0+pvary(i)
c        else
            s1=s0+vary(ipars(i))
            s2=s0-vary(ipars(i))
c        endif
        if(ipars(i).eq.4)then
            s1=abs(s1)
            s2=abs(s2)
        endif
        if(ipars(i).eq.6)then
            if(s1.gt.90.0d0)s1=180.0-s1    
            if(s1.lt.inclmin) s1=inclmin
            if(s2.gt.90.0d0)s2=180.0-s2    
            if(s2.lt.inclmin) s2=inclmin
            if(s1.eq.s2)then
                s1=inclmin
                s2=90.0
            endif
        endif
        sol2(ipars(i))=s1
        call transitmodel(1,x,exptime,dtype,y1,nfit,sol2)
        sol2(ipars(i))=s2
        call transitmodel(1,x,exptime,dtype,y2,nfit,sol2)
        dyda(i)=(y1-y2)/(s1-s2)
        sol2(ipars(i))=s0
c        write(6,*) i,a(i),dyda(i)
c        write(6,*) x,y1,y2
 12   continue
c      write(6,*) "what the fuck man"
        
c      do 12 i=1,na
c       dyda(i)=dfridr(a(i),vary(ipars(i)),err,ipars(i),sol2,nfit,x,
c     .  exptime)
c        if(ipars(i).eq.6)then
c            if(dyda(i).eq.0.0)then
c                tin=sol2(ipars(i))
c                call transitmodel(1,x,exptime,y1,nfit,sol2)
c                sol2(ipars(i))=(90+inclmin)/2.0
c                call transitmodel(1,x,exptime,y2,nfit,sol2)
c                dyda(i)=(y1-y2)/(tin-sol2(ipars(i)))
c            endif
c        endif    
cc       dyda(i)=dfridr(a(i),pvary(i),err,ipars(i),sol2,nfit,x,
cc     .  exptime)
cc        write(6,*) i,a(i),dyda(i),pvary(i)
c 12   continue 
 
c      write(6,*) "func:",x,y
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
                     goto 25!stop!pause 'singular matrix in gaussj' 
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
            goto 25!stop!pause 'singular matrix in gaussj'
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
 25   return 
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