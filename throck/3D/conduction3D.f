CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine conduction(nlat,nlon,nR,tcond,rad,drad,dtime,dtheta,
     .  dlam,theta,dEn,T,sarea,Eint)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nlat,nlon,nR,i,j,k,ii
      double precision tcond(nlat,nlon,nR),dtime,rad(nR),drad,
     .  theta(nlat),dtheta,dlam,
     .  dEn(nlat,nlon,nR),T(nlat,nlon,nR),cdist,area,dthetad2,cond,
     .  sarea(nlat),Eint,dradd2,Tn
      
      dthetad2=dtheta/2.0d0
      dradd2=drad/2.0d0
      
      do 18 k=1,nR

C       first we sweep in theta direction
      
C       start with corner pieces
        i=1 !nlat
        j=1 !nlon
        ii=nlon/2+j
        area=drad*dlam*sin(theta(i)+dthetad2)*rad(k)
        cdist=dtheta*rad(k)
        dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i+1,j,k),
     .    T(i,j,k),cdist,dtime) 
        area=drad*dlam*sin(theta(i)-dthetad2)*rad(k)
        dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(nlat,ii,k),
     .    T(i,j,k),cdist,dtime)
c      write(6,*) "area:",cond(tcond,area,T(i+1,j),T(i,j),cdist,dtime)
      
        i=nlat
        j=1
        ii=nlon/2+j
        area=drad*dlam*sin(theta(i)+dthetad2)*rad(k)
        cdist=dtheta*rad(k)
        dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(1,ii,k),
     .    T(i,j,k),cdist,dtime) 
        area=drad*dlam*sin(theta(i)-dthetad2)*rad(k)
        dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i-1,j,k),
     .    T(i,j,k),cdist,dtime)
      
        i=1
        j=nlon
        ii=nlon/2
        area=drad*dlam*sin(theta(i)+dthetad2)*rad(k)
        cdist=dtheta*rad(k)
        dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i+1,j,k),
     .    T(i,j,k),cdist,dtime) 
        area=drad*dlam*sin(theta(i)-dthetad2)*rad(k)
        dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(nlat,ii,k),
     .    T(i,j,k),cdist,dtime)
      
        i=nlat
        j=nlon
        ii=nlon/2
        area=drad*dlam*sin(theta(i)+dthetad2)*rad(k)
        cdist=dtheta*rad(k)
        dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(1,ii,k),
     .    T(i,j,k),cdist,dtime) 
        area=drad*dlam*sin(theta(i)-dthetad2)*rad(k)
        dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i-1,j,k),
     .    T(i,j,k),cdist,dtime)
      
C       bottom row
        i=1
        do 10 j=2,nlon-1
            ii=nlon/2+j
            if(ii.gt.nlon) ii=ii-nlon/2
            area=drad*dlam*sin(theta(i)+dthetad2)*rad(k)
            cdist=dtheta*rad(k)
            dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i+1,j,k),
     .        T(i,j,k),cdist,dtime) 
            area=drad*dlam*sin(theta(i)-dthetad2)*rad(k)
            dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(nlat,ii,k),
     .        T(i,j,k),cdist,dtime)
 10     continue
 
C       top row
        i=nlat
        do 11 j=2,nlon-1
            ii=nlon/2+j
            if(ii.gt.nlon) ii=ii-nlon/2
            area=drad*dlam*sin(theta(i)+dthetad2)*rad(k)
            cdist=dtheta*rad(k)
            dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(1,ii,k),
     .        T(i,j,k),cdist,dtime) 
            area=drad*dlam*sin(theta(i)-dthetad2)*rad(k)
            dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i-1,j,k),
     .        T(i,j,k),cdist,dtime)
 11     continue
      
C       the rest
        do 12 i=2,nlat-1
            do 13 j=1,nlon
                area=drad*dlam*sin(theta(i)+dthetad2)*rad(k)
                cdist=dtheta*rad(k)
                dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,
     .            T(i+1,j,k),T(i,j,k),cdist,dtime) 
                area=drad*dlam*sin(theta(i)-dthetad2)*rad(k)
                dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,
     .            T(i-1,j,k),T(i,j,k),cdist,dtime)
 13         continue
 12     continue
 
C       Now we sweep in the lambda direction

C       do the corner pieces first
        i=1
        j=1
        area=drad*dtheta*rad(k)
        cdist=dlam*sin(theta(i)+dthetad2)*rad(k)
        dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i,j+1,k),
     .    T(i,j,k),cdist,dtime) 
        cdist=dlam*sin(theta(i)-dthetad2)*rad(k)
        dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i,nlon,k),
     .    T(i,j,k),cdist,dtime)
      
        i=nlat
        j=1
        area=drad*dtheta*rad(k)
        cdist=dlam*sin(theta(i)+dthetad2)*rad(k)
        dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i,j+1,k),
     .    T(i,j,k),cdist,dtime) 
        cdist=dlam*sin(theta(i)-dthetad2)*rad(k)
        dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i,nlon,k),
     .    T(i,j,k),cdist,dtime)
      
        i=1
        j=nlon
        area=drad*dtheta*rad(k)
        cdist=dlam*sin(theta(i)+dthetad2)*rad(k)
        dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i,1,k),
     .    T(i,j,k),cdist,dtime) 
        cdist=dlam*sin(theta(i)-dthetad2)*rad(k)
        dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i,j-1,k),
     .    T(i,j,k),cdist,dtime)
      
        i=nlat
        j=nlon
        area=drad*dtheta*rad(k)
        cdist=dlam*sin(theta(i)+dthetad2)*rad(k)
        dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i,1,k),
     .    T(i,j,k),cdist,dtime) 
        cdist=dlam*sin(theta(i)-dthetad2)*rad(k)
        dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i,j-1,k),
     .    T(i,j,k),cdist,dtime)
      
C       left side
        j=1
        do 14 i=2,nlat-1
            area=drad*dtheta*rad(k)
            cdist=dlam*sin(theta(i)+dthetad2)*rad(k)
            dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i,j+1,k),
     .        T(i,j,k),cdist,dtime) 
            cdist=dlam*sin(theta(i)-dthetad2)*rad(k)
            dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i,nlon,k),
     .        T(i,j,k),cdist,dtime)
 14     continue
 
C       right side
        j=nlon
        do 15 i=2,nlat-1
            area=drad*dtheta*rad(k)
            cdist=dlam*sin(theta(i)+dthetad2)*rad(k)
            dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i,1,k),
     .        T(i,j,k),cdist,dtime) 
            cdist=dlam*sin(theta(i)-dthetad2)*rad(k)
            dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i,j-1,k),
     .        T(i,j,k),cdist,dtime)
 15     continue
      
C       the rest
        do 16 i=1,nlat
            do 17 j=2,nlon-1
                area=drad*dtheta*rad(k)
                cdist=dlam*sin(theta(i)+dthetad2)*rad(k)
                dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i,j+1,k),
     .              T(i,j,k),cdist,dtime) 
                cdist=dlam*sin(theta(i)-dthetad2)*rad(k)
                dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i,j-1,k),
     .              T(i,j,k),cdist,dtime)
 17         continue
 16     continue
     
 18   continue    
     
C     Now we handle the radial direction

C     deal with surface layer
      k=1
      do 23 j=1,nlon
        do 24 i=1,nlat
            area=sarea(i)*(rad(k)-dradd2)*(rad(k)-dradd2)
            cdist=drad
            dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i,j,k+1),
     .        T(i,j,k),cdist,dtime)
 24     continue
 23   continue

C     Inner most layer temperature stays fixed
C 
C     deal with inner most layer
      k=nR
      do 25 j=1,nlon
        do 26 i=1,nlat
            area=sarea(i)*(rad(k)+dradd2)*(rad(k)+dradd2)
            cdist=drad
            dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i,j,k-1),
     .        T(i,j,k),cdist,dtime)
            area=sarea(i)*(rad(k)-dradd2)*(rad(k)-dradd2)
c            Tn=2.0d0*T(i,j,k)-T(i,j,k-1)
c            if(Tn.gt.0.0d0)
c     .       dEn(i,j,k)=dEn(i,j,k)+cond(tcond,area,Tn,T(i,j,k)
c     .          ,cdist,dtime)
            dEn(i,j,k)=dEn(i,j,k)+Eint*area*dtime
 26     continue
 25   continue
 
C     deal with rest of layers
      do 27 k=2,nR-1
        do 28 j=1,nlon
            do 29 i=1,nlat
                area=sarea(i)*(rad(k)-dradd2)*(rad(k)-dradd2)
                cdist=drad
                dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i,j,k+1),
     .              T(i,j,k),cdist,dtime)
                area=sarea(i)*(rad(k)+dradd2)*(rad(k)+dradd2)
                dEn(i,j,k)=dEn(i,j,k)+cond(tcond(i,j,k),area,T(i,j,k-1),
     .              T(i,j,k),cdist,dtime)
 29         continue
 28     continue
 27   continue
     
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function cond(tcond,area,T2,T1,cdist,dtime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision tcond,area,T1,T2,cdist,dtime
      
      cond=tcond*area*(T2-T1)/cdist*dtime
      
      return
      end

