      subroutine gettrack(mass,Z,nmodelmax,nmodel,tage,tTeff,tlogL,trad,
     .  trho,tdrhodt,tmcore,tdloggdt,Tefferr,loggerr)
      implicit none
      integer nline,inda,indz,indm,i,k,j,nmodel,nmodelmax,noff
      double precision mass,Z,tage(nmodelmax),tTeff(nmodelmax),
     .  tlogL(nmodelmax),Tsun,trad(nmodelmax),trho(nmodelmax),
     .  tdrhodt(nmodelmax),Pi,Msun,Rsun,Rs,Ms,tmcore(nmodelmax),logg,
     .  dist,tdloggdt(nmodelmax),Tefferr,loggerr,G,logg1,logg2
      parameter (nline=150)
      parameter (inda=3)
      parameter (indz=11)
      parameter (indm=35)
      real*8 avalue(inda)
      real*8 zvalue(indz)
      real*8 xmass(indm)
      common /grid/avalue,zvalue,xmass
      real*8 work(5,nline,3),pwork(5,nline)
      data avalue/0.0,0.3,0.6/
      data zvalue/0.00001,0.0001,0.0004,0.001,0.004,0.007,
     +0.01,0.02,0.04,0.06,0.08/
      data xmass/ 0.4,0.5,0.6,0.7,0.8,0.9,1.0,
     +1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
     +2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,
     +3.2,3.4,3.6,3.8,4.0,4.2,4.5,5.0/
      integer ka,kz,km,kf,ntrack,m,kvar
      integer a_one,z_one,m_one
      common /single/a_one,z_one,m_one
      real*8 valpha,vz,vmass,y(3)
c      character*72 fout(2)

      Tsun=5781.6d0 !Solar temperature (K) 
      Pi=acos(-1.d0)!define Pi and 2*Pi
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      Msun=1.9891d30 !kg  mass of Sun
      G=6.67259E-11
      
      valpha=0.0
      vz=Z
      vmass=mass
c      write(0,34)valpha, vz, vmass
   34 format(1x,'Input : [a/Fe]=',f6.3,1x,'Z=',f8.5,1x,
     +       'M/Msol=',f6.3)
      
C  search the nearby tracks
      call findex(avalue,valpha,inda,ka,a_one)
      call findex(zvalue,vz,indz,kz,z_one)
      call findex(xmass,vmass,indm,km,m_one)

      ntrack=nline!24
c      do 100 kf=1,2
      kf=2
      kvar=ka
      if(a_one.eq.1)then
       call interp_z(vz,vmass,pwork,ntrack,kvar,kz,km,kf)
      else
       call interp_z(vz,vmass,work(1,1,1),ntrack,1,kz,km,kf)
       call interp_z(vz,vmass,work(1,1,2),ntrack,2,kz,km,kf)
       call interp_z(vz,vmass,work(1,1,3),ntrack,3,kz,km,kf)
         do 10 i=1,ntrack
         do 20 k=1,5
         do 1 m=1,3
          y(m)=work(k,i,m)
    1    continue
         call newton(avalue,y,3,valpha,pwork(k,i))
   20 continue
   10 continue
      endif
c      filename=fout(kf)
c      call write_track(filename,pwork,ntrack)
      
      noff=10 !removing PMS crap.
      ntrack=ntrack-noff
      do 30 j=1,ntrack
        tage(j)=pwork(1,j+noff)
        tTeff(j)=10.0d0**(pwork(2,j+noff))
        tlogL(j)=pwork(3,j+noff)
        trad(j)=sqrt( (Tsun/tTeff(j))**4.0d0*10.0d0**tlogL(j) )
        Ms=vmass*Msun
        Rs=trad(j)*Rsun
        trho(j)=Ms/(4.0d0/3.0d0*pi*Rs*Rs*Rs)
        tmcore(j)=pwork(5,j)
c        write(6,*) tlogL(j),trad(j),tteff(j)
 30   continue
      do 31 j=2,ntrack-1
        Ms=vmass*Msun
        Rs=trad(j-1)*Rsun
        logg1=log10(G*Ms/(Rs*Rs))+2.0d0
        Rs=trad(j)*Rsun
        logg2=log10(G*Ms/(Rs*Rs))+2.0d0
        if((loggerr.gt.0.0d0).and.(Tefferr.gt.0.0d0))then
            dist=sqrt(((logg2-logg1)/loggerr)**2.0+
     .          ((tTeff(j)-tTeff(j-1))/Tefferr)**2.0)
            tdloggdt(j)=dist/(tage(j)-tage(j-1))
        else
            tdloggdt(j)=1.0d0
        endif
        tdrhodt(j)=(trho(j)-trho(j-1))/(tage(j)-tage(j-1))
 31   continue
      tdrhodt(1)=tdrhodt(2)
      tdloggdt(1)=tdloggdt(2)
      nmodel=ntrack      
c      read(5,*)
      
c      ntrack=nline
c  100 continue
      
      return
      end