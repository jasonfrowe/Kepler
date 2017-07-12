      program corotscatterplot
      implicit none
      integer nunit,nmax,npt,nerr,iargc,nj
      parameter(nmax=104)
      real mag(nmax),rms(3,nmax),expt(3,nmax),ub,bper(nmax),noise(nmax)
      character*2 sptype(nmax)
      character*80 filename,cline
      
      ub=-1.0
      if(iargc().lt.2) goto 900
      call getarg(1,filename) !read in filename from command line
      call getarg(2,cline)
      read(cline,*) nj
      if((nj.lt.1).or.(nj.gt.3)) goto 900
      if(iargc().ge.3) then
        call getarg(3,cline)
        read(cline,*) ub
      endif
      
      nunit=10 !set unit number for file input
      open(unit=nunit,file=filename,status='old',err=901)
      call getdata(nunit,nmax,npt,mag,rms,expt,sptype,bper,nerr)
      close(nunit)
      if(nerr.eq.1) goto 902
      
      call rmsplot(nmax,npt,mag,rms,expt,sptype,bper,nj,ub,noise)
     
C     Error parsing 
      goto 999
 900  write(0,*) "Usage: corotscatterplot <filename> <nj>"
      write(0,*) " "
      write(0,*) "<filename> : datefile contain corot scatter stats"
      write(0,*) "<nj> : stats to use."
      write(0,*) "        1:  3 hour"
      write(0,*) "        2:  6 hour"
      write(0,*) "        3: 12 hour"
      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 902  write(0,*) "nmax must be increased to match length of ",filename
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rmsplot(nmax,npt,mag,rms,expt,sptype,bper,nj,ub,noise)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nmax,npt,i,j,nplot,sym,nj
      real mag(nmax),rms(3,nmax),expt(3,nmax),xp,yp,
     .  bounds(4),ub,temp(4),sf,bper(nmax),noise(nmax),s3,s12,p3,p12,n2
      character*2 sptype(nmax)
      character*80 title
       
      call pgopen('?')
      call pgpage()
      call pgask(.true.)
      call PGPAP ( 7.0 ,1.0) !square plot
c      call pgsch(2.0) !bigger text

      bounds(1)=mag(1)
      bounds(2)=mag(2)
      j=nj
      bounds(3)=rms(j,1)/expt(j,1)
      bounds(4)=rms(j,1)/expt(j,1)
      do 10 i=2,npt
        bounds(1)=min(mag(i),bounds(1))
        bounds(2)=max(mag(i),bounds(2))
        bounds(3)=min(rms(j,i)/expt(j,i),bounds(3))
        bounds(4)=max(rms(j,i)/expt(j,i),bounds(4))
 10   continue

      if(ub.gt.0) bounds(4)=ub
 
      sf=0.05
      temp(1)=bounds(1)-sf*(bounds(2)-bounds(1))
      temp(2)=bounds(2)+sf*(bounds(2)-bounds(1))
      temp(3)=bounds(3)-sf*(bounds(4)-bounds(3))
      temp(4)=bounds(4)+sf*(bounds(4)-bounds(3))
      
      do 12 i=1,4
        bounds(i)=temp(i)
 12   continue
      
      call pgwindow(bounds(1),bounds(2),bounds(3),bounds(4))
      call pgslw(3)
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      if(j.eq.1) then
        title="3 Hour"
      elseif(j.eq.2) then
        title="6 Hour"
      elseif(j.eq.3) then
        title="12 Hour"
      endif
      call pglabel("Magnitude","RMS/EXPT",title)
      call pgslw(1)
      
      call pgsch(2.0)
      call pgbbuf()
      do 11 i=1,npt
        if(sptype(i)(1:1).eq.'A') then
        	call pgsci(5) !cyan
        	if(bper(i).gt.0)then
        		sym=17
            else
        	    sym=4
        	endif
        elseif(sptype(i)(1:1).eq.'F') then
            call pgsci(4) !yellow
            if(bper(i).gt.0)then
            	sym=16
            else
            	sym=6
            endif
        elseif(sptype(i)(1:1).eq.'G') then
            call pgsci(3) !green
            if(bper(i).gt.0)then
            	sym=13
        	else
        		sym=7
            endif
        elseif(sptype(i)(1:1).eq.'K') then
            call pgsci(6) !magenta
            if(bper(i).gt.0)then
            	sym=-4
            else
                sym=11
            endif
        elseif(sptype(i)(1:1).eq.'M') then
            call pgsci(2) !red
            if(bper(i).gt.0)then
            	sym=18
            else
            	sym=12
            endif
        endif
        xp=mag(i)
        yp=rms(j,i)/expt(j,i)
        call pgpt1(xp,yp,sym)
c        write(0,*) i,mag(i),bper(i),sym
 11   continue
      call pgebuf()
      call pgsch(1.0)
      call pgsci(1)

      call pgpage()
      bounds(3)=rms(1,1)/rms(3,1)
      bounds(4)=bounds(3)
      do 14 i=2,npt
      	bounds(3)=min(rms(1,i)/rms(3,i),bounds(3))
      	bounds(4)=max(rms(1,i)/rms(3,i),bounds(4))
 14   continue
      temp(3)=bounds(3)-sf*(bounds(4)-bounds(3))
      temp(4)=bounds(4)+sf*(bounds(4)-bounds(3))
      do 15 i=3,4
        bounds(i)=temp(i)
 15   continue
      call pgslw(3)
      call pgwindow(bounds(1),bounds(2),bounds(3),bounds(4))
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Magnitude","RMS3/RMS12"," ")
      call pgslw(1)
      call pgsch(2.0)
      call pgbbuf()
      do 13 i=1,npt
        if(sptype(i)(1:1).eq.'A') then
        	call pgsci(5) !cyan
        	if(bper(i).gt.0)then
        		sym=17
            else
        	    sym=4
        	endif
        elseif(sptype(i)(1:1).eq.'F') then
            call pgsci(4) !blue
            if(bper(i).gt.0)then
            	sym=16
            else
            	sym=6
            endif
        elseif(sptype(i)(1:1).eq.'G') then
            call pgsci(3) !green
            if(bper(i).gt.0)then
            	sym=13
        	else
        		sym=7
            endif
        elseif(sptype(i)(1:1).eq.'K') then
            call pgsci(6) !magenta
            if(bper(i).gt.0)then
            	sym=-4
            else
                sym=11
            endif
        elseif(sptype(i)(1:1).eq.'M') then
            call pgsci(2) !red
            if(bper(i).gt.0)then
            	sym=18
            else
            	sym=12
            endif
        endif
      	xp=mag(i)
        yp=rms(1,i)/rms(3,i)
        call pgpt1(xp,yp,sym)
 13   continue
      call pgebuf()
      call pgsch(1.0)
      call pgsci(1)
 
      call pgpage()
	  
	  do 16 i=1,npt
c	    s3=rms(1,i)*rms(1,i)
c	    s12=rms(3,i)*rms(3,i)
c	    p3=expt(1,i)*expt(1,i)
c	    p12=expt(3,i)*expt(3,i)
c	    n2=(s3*p12 - s12*p3)/(s12-s3)
        n2=rms(3,i)*rms(3,i)-expt(3,i)*expt(3,i)	    
	    if(n2.ge.0) then
	  		noise(i)=sqrt(n2)
	    else
	        noise(i)=-1.0 !mark errors (invalid assumptions)
	    endif
c	  	write(0,*) noise(i),s3,s12,p3,p12
 16   continue
      bounds(3)= 99.9e30
      bounds(4)=-99.9e30
      do 17 i=1,npt
      	if(noise(i).gt.0)bounds(3)=min(noise(i),bounds(3))
      	if(noise(i).gt.0)bounds(4)=max(noise(i),bounds(4))
 17   continue
      temp(3)=bounds(3)-sf*(bounds(4)-bounds(3))
      temp(4)=bounds(4)+sf*(bounds(4)-bounds(3))
      do 18 i=3,4
        bounds(i)=temp(i)
 18   continue
      call pgslw(3)
      call pgwindow(bounds(1),bounds(2),0.0,1000.0)
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Magnitude","Noise (umag)"," ")
      call pgslw(1)
      call pgsch(2.0)
      call pgbbuf()
      do 19 i=1,npt
        if(sptype(i)(1:1).eq.'A') then
        	call pgsci(5) !cyan
        	if(bper(i).gt.0)then
        		sym=17
            else
        	    sym=4
        	endif
        elseif(sptype(i)(1:1).eq.'F') then
            call pgsci(4) !yellow
            if(bper(i).gt.0)then
            	sym=16
            else
            	sym=6
            endif
        elseif(sptype(i)(1:1).eq.'G') then
            call pgsci(3) !green
            if(bper(i).gt.0)then
            	sym=13
        	else
        		sym=7
            endif
        elseif(sptype(i)(1:1).eq.'K') then
            call pgsci(6) !magenta
            if(bper(i).gt.0)then
            	sym=-4
            else
                sym=11
            endif
        elseif(sptype(i)(1:1).eq.'M') then
            call pgsci(2) !red
            if(bper(i).gt.0)then
            	sym=18
            else
            	sym=12
            endif
        endif
      	xp=mag(i)
        yp=noise(i)
        if((yp.ge.0).and.(bper(i).le.0))call pgpt1(xp,yp,sym)
        if((yp.ge.0).and.(bper(i).le.0).and.(xp.gt.11.0).and.
     .		(xp.lt.13.0))then
        	write(0,*) xp,yp,sptype(i)
        endif
 19   continue
      call pgebuf()
      call pgsch(1.0)
      call pgsci(1)
 
      call pgclos()
 
      return
      end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getdata(nunit,nmax,npt,mag,rms,expt,sptype,bper,nerr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nmax,i,npt,id,j,nerr
      real mag(nmax),rms(3,nmax),expt(3,nmax),Teff,bper(nmax)
      character dumc
      character*2 sptype(nmax)
      character*3 cl
      
      do 10 i=1,3 !3 first lines are comments
      	read(nunit,*) dumc
 10   continue
 
      i=1
 11   read(nunit,*,end=12) id,mag(i),Teff,sptype(i),cl,
     .	  (rms(j,i),j=1,3),(expt(j,i),j=1,3),bper(i)
c         do 13 j=1,3
c         	expt(j,i)=(expt(j,i)/sqrt(8.0))
c 13      continue
c        write(6,*) i,id,mag(i),sptype(i),cl,(rms(j,i),j=1,3)
      	i=i+1
      	if(i.gt.nmax+1) goto 901
      	goto 11
 12   continue
      npt=i-1
      
C     Error parsing
      nerr=0
      goto 999
 901  nerr=1
 999  return
      end