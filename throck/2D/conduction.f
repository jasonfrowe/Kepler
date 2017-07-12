CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine conduction(nlat,nlon,tcond,rad,drad,dtime,dtheta,dlam,
     .  theta,dEn,T)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nlat,nlon,i,j,k
      double precision tcond,dtime,rad,drad,theta(nlat),dtheta,dlam,
     .  dEn(nlat,nlon),T(nlat,nlon),cdist,area,dthetad2,cond
      
      dthetad2=dtheta/2.0d0
      
C     first we sweep in theta direction
      
C     start with corner pieces
      i=1 !nlat
      j=1 !nlon
      k=nlon/2+j
      area=drad*dlam*sin(theta(i)+dthetad2)*rad
      cdist=dtheta*rad
      dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i+1,j),T(i,j),cdist,dtime) 
      area=drad*dlam*sin(theta(i)-dthetad2)*rad
      dEn(i,j)=dEn(i,j)+cond(tcond,area,T(nlat,k),T(i,j),cdist,dtime)
c      write(6,*) "area:",cond(tcond,area,T(i+1,j),T(i,j),cdist,dtime)
      
      i=nlat
      j=1
      k=nlon/2+j
      area=drad*dlam*sin(theta(i)+dthetad2)*rad
      cdist=dtheta*rad
      dEn(i,j)=dEn(i,j)+cond(tcond,area,T(1,k),T(i,j),cdist,dtime) 
      area=drad*dlam*sin(theta(i)-dthetad2)*rad
      dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i-1,j),T(i,j),cdist,dtime)
      
      i=1
      j=nlon
      k=nlon/2
      area=drad*dlam*sin(theta(i)+dthetad2)*rad
      cdist=dtheta*rad
      dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i+1,j),T(i,j),cdist,dtime) 
      area=drad*dlam*sin(theta(i)-dthetad2)*rad
      dEn(i,j)=dEn(i,j)+cond(tcond,area,T(nlat,k),T(i,j),cdist,dtime)
      
      i=nlat
      j=nlon
      k=nlon/2
      area=drad*dlam*sin(theta(i)+dthetad2)*rad
      cdist=dtheta*rad
      dEn(i,j)=dEn(i,j)+cond(tcond,area,T(1,k),T(i,j),cdist,dtime) 
      area=drad*dlam*sin(theta(i)-dthetad2)*rad
      dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i-1,j),T(i,j),cdist,dtime)
      
C     bottom row
      i=1
      do 10 j=2,nlon-1
        k=nlon/2+j
        if(k.gt.nlon) k=k-nlon/2
        area=drad*dlam*sin(theta(i)+dthetad2)*rad
        cdist=dtheta*rad
        dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i+1,j),T(i,j),cdist,dtime) 
        area=drad*dlam*sin(theta(i)-dthetad2)*rad
        dEn(i,j)=dEn(i,j)+cond(tcond,area,T(nlat,k),T(i,j),cdist,dtime)
 10   continue
 
C     top row
      i=nlat
      do 11 j=2,nlon-1
        k=nlon/2+j
        if(k.gt.nlon) k=k-nlon/2
        area=drad*dlam*sin(theta(i)+dthetad2)*rad
        cdist=dtheta*rad
        dEn(i,j)=dEn(i,j)+cond(tcond,area,T(1,k),T(i,j),cdist,dtime) 
        area=drad*dlam*sin(theta(i)-dthetad2)*rad
        dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i-1,j),T(i,j),cdist,dtime)
 11   continue
      
C     the rest
      do 12 i=2,nlat-1
        do 13 j=1,nlon
          area=drad*dlam*sin(theta(i)+dthetad2)*rad
          cdist=dtheta*rad
          dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i+1,j),T(i,j),cdist,dtime) 
          area=drad*dlam*sin(theta(i)-dthetad2)*rad
          dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i-1,j),T(i,j),cdist,dtime)
 13     continue
 12   continue
 
C     Now we sweep in the lambda direction

C     do the corner pieces first
      i=1
      j=1
      area=drad*dtheta*rad
      cdist=dlam*sin(theta(i)+dthetad2)*rad
      dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i,j+1),T(i,j),cdist,dtime) 
      cdist=dlam*sin(theta(i)-dthetad2)*rad
      dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i,nlon),T(i,j),cdist,dtime)
      
      i=nlat
      j=1
      area=drad*dtheta*rad
      cdist=dlam*sin(theta(i)+dthetad2)*rad
      dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i,j+1),T(i,j),cdist,dtime) 
      cdist=dlam*sin(theta(i)-dthetad2)*rad
      dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i,nlon),T(i,j),cdist,dtime)
      
      i=1
      j=nlon
      area=drad*dtheta*rad
      cdist=dlam*sin(theta(i)+dthetad2)*rad
      dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i,1),T(i,j),cdist,dtime) 
      cdist=dlam*sin(theta(i)-dthetad2)*rad
      dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i,j-1),T(i,j),cdist,dtime)
      
      i=nlat
      j=nlon
      area=drad*dtheta*rad
      cdist=dlam*sin(theta(i)+dthetad2)*rad
      dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i,1),T(i,j),cdist,dtime) 
      cdist=dlam*sin(theta(i)-dthetad2)*rad
      dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i,j-1),T(i,j),cdist,dtime)
      
C     left side
      j=1
      do 14 i=2,nlat-1
        area=drad*dtheta*rad
        cdist=dlam*sin(theta(i)+dthetad2)*rad
        dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i,j+1),T(i,j),cdist,dtime) 
        cdist=dlam*sin(theta(i)-dthetad2)*rad
        dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i,nlon),T(i,j),cdist,dtime)
 14   continue
 
C     right side
      j=nlon
      do 15 i=2,nlat-1
        area=drad*dtheta*rad
        cdist=dlam*sin(theta(i)+dthetad2)*rad
        dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i,1),T(i,j),cdist,dtime) 
        cdist=dlam*sin(theta(i)-dthetad2)*rad
        dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i,j-1),T(i,j),cdist,dtime)
 15   continue
      
C     the rest
      do 16 i=1,nlat
        do 17 j=2,nlon-1
          area=drad*dtheta*rad
          cdist=dlam*sin(theta(i)+dthetad2)*rad
          dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i,j+1),T(i,j),cdist,dtime) 
          cdist=dlam*sin(theta(i)-dthetad2)*rad
          dEn(i,j)=dEn(i,j)+cond(tcond,area,T(i,j-1),T(i,j),cdist,dtime)
 17     continue
 16   continue
     
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

