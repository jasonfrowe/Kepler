      program boothisto
      implicit none
      integer nmax,npt,col,nbin,i,j,k,nfit,nunit,npt2,
     .   nbinmax
      parameter(nmax=3500000,nfit=17+12,nbinmax=500)
      integer nd(nmax)
      real dd(nmax),datamin,datamax,bdatax(nmax),
     .   bdatay(nmax),bmax,fc,ave,var,pi,dif,ans(2,nmax),
     .	 errs(6),sol(nfit),serr(nfit,2),Dpvary(nfit),err(nfit),doe,toff
      real Psec(nmax),M1(nmax),M2(nmax),R1(nmax),G,Msun,Mearth,Mjup,
     .  Rsun,aConst,asemi,incl(nmax),out1(nbinmax,nfit),AU,tPi,
     .  out2(nbinmax,nfit),bsize,Rjup,R2(nmax),med(nfit),mode(nfit),
     .  momax,pid2,eccn,ecw(nmax),esw(nmax),adrs,zrs,b,temp(4),inc
      double precision Kamp,dtmp(5)
      character*80 filename,inputsol,cline
      character*80 titles(nfit)
      titles(1)='Stellar Mass (M\d\(2281)\u)'
      titles(2)='Planet Mass (M\dJ\u)'
      titles(3)='Stellar Radius (R\d\(2281)\u)'
      titles(4)='Planet Radius (R\dJ\u)'
      titles(5)='Period (days)'
      titles(6)='Inclination (deg)'
      titles(7)='Epoch (HJD-2454900)'
      titles(8)='Photometric Zeropoint (ppm)'
      titles(9)='Geometric Albedo'
      titles(10)='Non-Linear Limbdarkening(1)'
      titles(11)='Non-Linear Limbdarkening(2)'
      titles(12)='Non-Linear Limbdarkening(3)'
      titles(13)='Non-Linear Limbdarkening(4)'
      titles(14)="e cos w"!"Eccentricity"
      titles(15)="e sin w"!"Argument of Pericenter (deg)"
      titles(16)='Secondary Eclipse Depth (ppm)'
      titles(17)="Velocity Offset (m/s)"
      titles(18)="a/R*"
      titles(19)="Rp/R*"
      titles(20)="b"
      titles(21)="Planet Density (g/cm\u3\d)"
      titles(22)="a (AU)"
      titles(23)="Planet logg (cgs)"
      titles(24)="Stellar Density (g/cm\u3\d)"
      titles(25)="Stellar logg (cgs)"
      titles(26)="K (m/s)"
      titles(27)="e"
      titles(28)="w (radians)"
      titles(29)="Transit Duration"
      
      Pi=acos(-1.d0)   !Pi
      pid2=Pi/2.0 !Pi/2
      tpi=2.0d0*Pi
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      Msun=1.9891d30 !kg  mass of Sun
      Mearth=5.974d24 !kg mass of Earth
      Mjup=317.833d0*Mearth !kg  mass of Jupiter  
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      Rjup=142980.0d0*1000.0d0/2.0d0 !m  radius of Jupiter 
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)
      AU=1.49598e11
      
      nbin=30 !default binning for histograms
      if(iargc().lt.2) goto 903 !need one commandline arguement
     
      write(0,*) "Reading input solution"
      call getarg(1,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
      call getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff)
      close(nunit) !release unit number as we are done with file
      
      call getarg(2,filename)
      open(unit=10,file=filename,status='old',err=901)
      if(iargc().ge.3)then
      	call getarg(3,cline)
      	read(cline,*) nbin
      endif
      if(nbin.lt.2) goto 902

      if(nbin.gt.nbinmax) goto 904
            
      call pgopen('?') !open PGPlot device
      call pgask(.true.) !don't ask for new page.. just do it.
      call PGPAP ( 7.0 ,1.0)
      call pgsubp(5,5)
c      call pgpage() !start a multi-page document

C     Now we loop over all Frequencies.
      do 10 i=1,nfit !there are 3 variables for each frequency.
c      do 10 i=1,17
c         if(i.eq.8) goto 10
c         if(i.eq.17) goto 10
c         write(6,*) i
         col=i
         npt=nmax
         if(i.eq.18)then !a/R*
            open(unit=12,file="adrs.dat")
            do 17 j=1,npt2
                asemi=(M1(j)+M2(j))**(1.0e0/3.0e0)*
     .            Psec(j)**(2.0d0/3.0d0)*aConst
                dd(j)=asemi/R1(j)!*tan(Pi*(90.0-incl(j))/180.0)
                write(12,*) M1(j)/Msun,R1(j)/Rsun,incl(j),dd(j)
 17         continue
            close(12)
            asemi=(sol(1)*Msun+sol(2)*Mjup)**(1.0e0/3.0e0)*
     .        (sol(5)*8.64e4)**(2.0d0/3.0d0)*aConst
            sol(18)=asemi/(sol(3)*Rsun)!*tan(Pi*(90.0-incl(j))/180.0)
            npt=npt2
         elseif(i.eq.19)then !Rp/R*
            do 21 j=1,npt2
                dd(j)=R2(j)/R1(j)
 21         continue
            sol(19)=sol(4)*Rjup/(sol(3)*Rsun)
            npt=npt2
         elseif(i.eq.20)then !b
            do 23 j=1,npt2
                asemi=(M1(j)+M2(j))**(1.0e0/3.0e0)*
     .            Psec(j)**(2.0d0/3.0d0)*aConst
                dd(j)=asemi/R1(j)*tan(Pi*(90.0-incl(j))/180.0)
 23         continue
            asemi=(sol(1)*Msun+sol(2)*Mjup)**(1.0e0/3.0e0)*
     .        (sol(5)*8.64e4)**(2.0d0/3.0d0)*aConst
            sol(20)=asemi/(sol(3)*Rsun)*tan(Pi*(90.0-sol(6))/180.0)
            npt=npt2
         elseif(i.eq.21)then !planet density 
            do 24 j=1,npt2
                dd(j)=M2(j)/(4.0*Pi*R2(j)**3.0/3.0)/1000.0
 24         continue
            sol(21)=Mjup*sol(2)/(4.0*Pi*(Rjup*sol(4))**3.0/3.0)/1000.0
            npt=npt2
         elseif(i.eq.22)then !asemi
            do 25 j=1,npt2
                dd(j)=(M1(j)+M2(j))**(1.0e0/3.0e0)*
     .            Psec(j)**(2.0d0/3.0d0)*aConst/AU
 25         continue
            sol(22)=(sol(1)*Msun+sol(2)*Mjup)**(1.0e0/3.0e0)*
     .        (sol(5)*8.64e4)**(2.0d0/3.0d0)*aConst/AU
            npt=npt2
         elseif(i.eq.23)then !logg_p
            sol(23)=log10(1.0e2*G*Mjup*sol(2)/(Rjup*Rjup*sol(4)*sol(4)))
            k=0
            do 26 j=1,npt2
                if((M2(j).gt.0.0).and.(R2(j).gt.0.0))then
                    k=k+1
                    dd(k)=log10(1.0e2*G*M2(j)/(R2(j)*R2(j)))
                endif
 26         continue
            npt=k
         elseif(i.eq.24)then !stellar density 
            do 28 j=1,npt2
                dd(j)=M1(j)/(4.0*Pi*R1(j)**3.0/3.0)/1000.0
 28         continue
            sol(24)=Msun*sol(1)/(4.0*Pi*(Rsun*sol(3))**3.0/3.0)/1000.0
            npt=npt2
         elseif(i.eq.25)then !stellar logg
            do 29 j=1,npt2
                dd(j)=log10(1.0e2*G*M1(j)/(R1(j)*R1(j)))
 29         continue
            sol(25)=log10(1.0e2*G*Msun*sol(1)/(Rsun*Rsun*sol(3)*sol(3)))
            npt=npt2
         elseif(i.eq.26)then !K
c            write(6,*) dtmp(1),dtmp(2),dtmp(3),dtmp(4)
            eccn=sqrt(ecw(j)*ecw(j)+esw(j)*esw(j))
            dtmp(1)=2.0d0*dble(pi*G)*dble(Mjup*sol(2))**3
            dtmp(2)=dble((sin(Pi*(90.0-sol(6))/180.0+Pid2))**3)
            dtmp(3)=dble(sol(5)*8.64e4*(1.0d0-eccn*
     .          eccn)**(3.0e0/2.0e0))
            dtmp(4)=dble(msun*sol(1)+mjup*sol(2))*
     .          dble(msun*sol(1)+mjup*sol(2))
            Kamp=dtmp(1)*dtmp(2)/(dtmp(3)*dtmp(4))
            Kamp=Kamp**(1.0d0/3.0d0)
            sol(26)=real(Kamp)
c            write(6,*) "K",sol(26)
c            write(6,*) dtmp(1),dtmp(2),dtmp(3),dtmp(4)
            do 30 j=1,npt2
                eccn=sqrt(ecw(j)*ecw(j)+esw(j)*esw(j))                
                dtmp(1)=2.0d0*dble(pi*G)*dble(M2(j))**3
                dtmp(2)=dble((sin(Pi*(90.0-incl(j))/180.0+Pid2))**3)
                dtmp(3)=dble(Psec(j)*(1.0d0-eccn*
     .              eccn)**(3.0d0/2.0d0))
                dtmp(4)=dble(M1(j)+M2(j))*dble(M1(j)+M2(j))
c                write(6,*) dble(M1(j)+M2(j))*dble(M1(j)+M2(j))
c                write(6,*) M1(j)+M2(j)
                Kamp=dtmp(1)*dtmp(2)/(dtmp(3)*dtmp(4))
                Kamp=Kamp**(1.0d0/3.0d0)
                dd(j)=real(Kamp)
c                write(6,*) dtmp(1),dtmp(2),dtmp(3),dtmp(4)
c                write(6,*) "K",dd(j)
c                read(5,*)
 30         continue
            npt=npt2          
         elseif(i.eq.27) then !eccentricity 
            do 33 j=i,npt2
                dd(j)=sqrt(ecw(j)*ecw(j)+esw(j)*esw(j))
c                if(dd(j).ne.0.0d0) then
c                    write(6,*) dd(j),j
c                endif
 33         continue
c            do 40 j=1,npt2
c                if(dd(j).ne.0.0d0) then
c                    write(6,*) dd(j),j
c                endif
c 40         continue
            sol(27)=sqrt(sol(14)*sol(14)+sol(15)*sol(15)) 
            npt=npt2
         elseif(i.eq.28) then !w (omega)
            do 34 j=1,npt2
                dtmp(1)=sqrt(ecw(j)*ecw(j)+esw(j)*esw(j))
                if(dtmp(1).eq.0.0d0)then
                    dd(j)=0.0d0
                else
                    dd(j)=atan(esw(j)/ecw(j))
                    if((ecw(j).gt.0.0d0).and.(esw(j).lt.0.0d0))then
                        dd(j)=tPi+dd(j)
                    elseif((ecw(j).lt.0.0d0).and.(esw(j).gt.0.0d0))then 
                        dd(j)=Pi+dd(j)
                    elseif((ecw(j).lt.0.0d0).and.(esw(j).lt.0.0d0))then
                        dd(j)=Pi+dd(j)
                    endif
                endif
c                write(6,*) dtmp(1),dd(j),ecw(j)
c                read(5,*)
 34         continue
            dtmp(1)=sqrt(sol(14)*sol(14)+sol(15)*sol(15))
            if(dtmp(1).eq.0.0d0)then
                sol(28)=0.0d0
            else
                sol(28)=atan(sol(15)/sol(14))
                if((sol(14).gt.0.0d0).and.(sol(15).lt.0.0d0))then
                    sol(28)=tPi+sol(28)
                elseif((sol(14).lt.0.0d0).and.(sol(15).gt.0.0d0))then 
                    sol(28)=Pi+sol(28)
                elseif((sol(14).lt.0.0d0).and.(sol(15).lt.0.0d0))then
                    sol(28)=Pi+sol(28)
                endif
            endif
            npt=npt2
         elseif(i.eq.29) then !Transit duration.
c            do 36 j=1,npt2 !transit duration
c                asemi=(M1(j)+M2(j))**(1.0e0/3.0e0)*
c     .            Psec(j)**(2.0d0/3.0d0)*aConst
c                adrs=asemi/R1(j)
c                b=asemi/R1(j)*tan(Pi*(90.0-incl(j))/180.0)
c                eccn=sqrt(ecw(j)*ecw(j)+esw(j)*esw(j)) 
c                zrs=adrs*(tPi/Psec(j))/sqrt(1.0-b*b)*
c     .              sqrt(1.0-eccn*eccn)/(1.0+esw(j))
c                dd(j)=2.0/zrs/3600.0
cc                write(6,*) dd(j),zrs,adrs
cc                write(6,*) b,eccn,asemi
cc                read(5,*)
c 36         continue
c            npt=npt2
            do 38 j=1,npt2
                asemi=(M1(j)+M2(j))**(1.0e0/3.0e0)*
     .            Psec(j)**(2.0d0/3.0d0)*aConst
                inc=Pi*incl(j)/180.0d0
                temp(1)=Psec(j)/Pi
                temp(2)=R1(j)/asemi
                temp(3)=(1+(R2(j)/R1(j)))**2.0-
     .              ((asemi/R1(j))*cos(inc))**2.0
                temp(4)=1-cos(inc)*cos(inc)     
                dd(j)=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))/3600.0
c                write(0,*) temp(1),temp(2),temp(3),temp(4)
c                write(0,*) dd(j),inc,asemi
c                read(5,*)
 38         continue
            npt=npt2
            asemi=(msun*sol(1)+mjup*sol(2))**(1.0e0/3.0e0)*
     .          (sol(5)*8.64e4)**(2.0d0/3.0d0)*aConst
            inc=Pi*sol(6)/180.0d0
            temp(1)=sol(5)*8.64e4/Pi
            temp(2)=rsun*sol(3)/asemi
            temp(3)=(1+(rjup*sol(4)/(rsun*sol(3))))**2.0-
     .          ((asemi/(sol(3)*rsun))*cos(inc))**2.0
            temp(4)=1-cos(inc)*cos(inc)     
            sol(29)=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))/3600.0
         else
            call getdata(col,npt,dd) !read column of bootstrapped values
c            write(6,*) "npt2",npt
            if(i.eq.5)then
                do 13 j=1,npt
                    Psec(j)=dd(j)*8.64d4
 13             continue
                npt2=npt
            elseif(i.eq.1)then
                do 14 j=1,npt
                    M1(j)=dd(j)*Msun
 14             continue
            elseif(i.eq.2)then
                do 15 j=1,npt
                    M2(j)=dd(j)*Mjup
 15             continue
            elseif(i.eq.3)then
                do 16 j=1,npt
                    R1(j)=dd(j)*Rsun
 16             continue
            elseif(i.eq.4)then
                do 22 j=1,npt
                    R2(j)=dd(j)*Rjup
 22             continue
            elseif(i.eq.6)then
                do 18 j=1,npt
                    incl(j)=dd(j)
 18             continue
            elseif(i.eq.14)then
                do 31 j=1,npt
                    ecw(j)=dd(j)
 31             continue
            elseif(i.eq.15)then
                do 32 j=1,npt
                    esw(j)=dd(j)
 32             continue              
            endif
         endif

c      if(i.eq.27)then
c        do 39 j=1,npt
c            if(dd(j).ne.0.0d0) then
c                write(6,*) dd(j),j
c            endif
c 39     continue
c      endif         
         
C        Find median          
         call rqsort(npt,dd,nd) !changed npt to k
         if(npt/2.gt.0) then
            med(i)=dd(nd(npt/2))
         else
            med(i)=dd(1)
         endif


C     These lines will sigclip at 4 sigma if required.  (probably not)
         call avevar(dd,npt,med(i),var)
c         write(6,*) i,med(i),sqrt(var)
         k=0
         do 12 j=1,npt
            if(abs(dd(j)-med(i)).le.4.0*sqrt(var))then
               k=k+1
               dd(k)=dd(j)
            endif
 12      continue
         npt=k
C     End of sigmaclipping

C        Find median          
         call rqsort(npt,dd,nd) !changed npt to k
         if(npt/2.gt.0) then
            med(i)=dd(nd(npt/2))
         else
            med(i)=dd(1)
         endif

c         if(i.eq.7)then
c            do 37 j=1,npt
c                dd(j)=dd(j)-med(i)
c 37         continue
c            med(i)=0.0d0
c         endif
         
         call avevar(dd,npt,med(i),var)
C     If we have a phase parameter, then we need to check for cyclic 
C     numbers around Pi

         ans(1,i)=med(i)
         ans(2,i)=sqrt(var)
      
         call findrange(npt,dd,datamin,datamax)
         if(datamin.eq.datamax)goto 11
c         write(6,*) "datamin,datamax",datamin,datamax
         call bindata(nbin,npt,dd,bdatax,bdatay,datamin,datamax,bmax)

         bsize=(bdatax(2)-bdatax(1))/2.0
         do 19 j=1,nbin
            out1(j,i)=bdatax(j)+bsize
            out2(j,i)=bdatay(j)/bmax
 19      continue
 
C        find mode
         mode(i)=bdatax(1)+bsize
         momax=bdatay(1)
         do 27 j=2,nbin
            if(bdatay(j).gt.momax)then
                momax=bdatay(j)
                mode(i)=bdatax(j)+bsize
            endif
 27      continue           
      
         call pgpage() !fresh plotting surface
         call pgslw(1)
         call pgsch(2.0)
c         call windowsetup(xb1,xb2,yb1,yb2) !make a square plotting surface
c         call pgvport(xb1,xb2,yb1,yb2)
         call pgvport(0.2,1.0,0.2,0.9)
         fc=2.*(datamax-datamin)/real(nbin) !center histograms
         call pgwindow(datamin,datamax,0.,bmax+0.1*bmax) !set size
         
C        Add axis labels
         call pglabel(titles(i),"Relative Probability","")
         call pgbin(nbin,bdatax,bdatay,.false.) !plot the histogram
         
         
         if(i.eq.8)then
            call pgwindow(datamin*1.0e6,datamax*1.0e6,0.,1.1)
c         elseif(i.eq.7)then
c            call pgwindow(datamin-60.0,datamax-60.0,0.,1.1)
         else
            call pgwindow(datamin,datamax,0.,1.1) !set size
         endif
         call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0) !add boarders
          
c         if(i.eq.6) sol(i)=ans(1,i)
         call errorest(npt,dd,nbin,bdatax,bdatay,bmax,med(i),errs)

c         write(6,*) "npt:",npt
         if(i.eq.8)then
            write(6,500) med(i)*1.0e6,mode(i)*1.0e6,ans(2,i)*1.0e6,
     .          (errs(j)*1.0e6,j=1,6)
         else
            write(6,500) sol(i),med(i),ans(2,i),(errs(j),j=1,6)
         endif
c 500     format(9(1PE16.9,1X))
c 500     format(9(F11.7,1X))
 500     format(9(F11.6,1X))
c 500     format(9(F15.12,1X))
         
C        this line is heavy on disk-io but easy on memory.
 11      rewind(10) ! rewind the bootstrap file so we can read it again.
c         read(10,*) dumi !not needed, but no harm. Just skips first line
 10   continue
      call pgslw(1)
      call pgclos() !close plotting device
      
c      do 12 i=0,nfreq-1
c         write(6,500) (ans(1,j),j=3*i+1,3*i+3),(ans(2,j),j=3*i+1,3*i+3)
c 12   continue
c 500  format(200(E14.7,1X))

c      open(unit=13,file="probs.dat")
c      write(13,502) "M*","M*_P","Mp","Mp_P","R*","R*_P","Rp","Rp_P",
c     .  "Per","Per_P","i","i_P","E","E_p","Ecl","Ecl_P","a/R*","a/R*_P"
c 502  format(18(A9,3X))
c      do 20 j=1,nbin
c	  write(13,501) out1(j,1),out2(j,1),out1(j,2),out2(j,2),
c     .    out1(j,3),out2(j,3),out1(j,4),out2(j,4),out1(j,5),out2(j,5),
c     .    out1(j,6),out2(j,6),out1(j,7),out2(j,7),out1(j,16),out2(j,16),
c     .    out1(j,18),out2(j,18)
cc        write(6,501) out1(j,1),out2(j,1)
c 20   continue
c      close(13)
c 501  format(18(F11.7,1X))
      
      
      close(10)
      goto 999
 901  write(6,*) "Cannot open ",filename
      goto 999
 902  write(6,*) "Error: Nbin must be at least 2"
      goto 999
 903  write(6,*) "Usage: boothist <freqname> <bootstats> [nbin]"
      write(6,*) "<freqname> : Frequency output list from MOSTPer"
      write(6,*) "<bootstats>: Bootstrap statistics from COSBoot"
      write(6,*) "[nbin] (optional) : binned parameter for histrograms, 
     .default=30"
      goto 999
 904  write(6,*) "Increase nbinmax to at least",nbin
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine errorest(npt,dd,nbin,bdatax,bdatay,bmax,frsol,errs)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nbin,nfreq,i,j,mdpnt,k,nfac
      real dd(npt),bdatax(nbin),bdatay(nbin),bmax,frsol,errs(6),shgt,
     .	intl,inth,per,perold,e1,e2,inthold,intlold
      logical loop
      
C     find out which bin contains the middle.      
      do 5 i=1,nbin-1
c        write(0,*) frsol,bdatax(i),bdatax(i+1)
      	if((bdatax(i).le.frsol).and.(bdatax(i+1).gt.frsol))then
        	mdpnt=i
      	endif
 5    continue
C     case when "bestfit" is outside histogram (i.e. eccentricity)
      if(frsol.le.bdatax(1)) mdpnt=1
      if(frsol.ge.bdatax(nbin)) mdpnt=nbin 
c      write(0,*) "mdpnt:",mdpnt
      
      perold=0 !initalization of variables
      intlold=0
      inthold=0
      do 14 i=1,6  !initialize variable to zero.
        errs(i)=0.0
 14   continue
C     If prog is slow, increase nfac, but lose accuracy 
      nfac=10
C     bmax is the value of the largest histogram bin
      do 13 k=1,bmax-1,nfac
      	shgt=real(bmax-k)
C     	first we move in the forward direction      
      	loop=.true. !loop until loop is false
      	i=mdpnt  !removed -1
      	do 10 while(loop)
      		if((bdatay(i).ge.shgt).and.(bdatay(i+1).lt.shgt))then
      			loop=.false.
c      			inth=bdatax(i)
      			inth=(bdatay(i+1)-shgt)/(bdatay(i+1)-bdatay(i))
     .				*(bdatax(i+1)-bdatax(i))+
     .				bdatax(i)
      		endif
        	if(i.ge.nbin-1) then
        		loop=.false.
        		inth=bdatax(nbin)
    		endif
    		i=i+1
c     		write(6,*) "i:",i
 10   	enddo
       
C     	now we move in the reverse direction
      	loop=.true. !loop until loop is false
        i=mdpnt+1
      	do 11 while(loop)
      		if((bdatay(i).ge.shgt).and.(bdatay(i-1).lt.shgt))then
      			loop=.false.
c      			intl=bdatax(i)
				intl=bdatax(i)-
     .				(bdatay(i)-shgt)/(bdatay(i)-bdatay(i-1))
     .				*(bdatax(i)-bdatax(i-1))
      		endif
        	if(i.le.2) then
        		loop=.false.
        		intl=bdatax(1)
    		endif
    		i=i-1
 11   	enddo

        j=0
      	do 12 i=1,npt,nfac
      		if((dd(i).gt.intl).and.(dd(i).lt.inth))then
      			j=j+1    
      		endif
 12	  	continue
        j=j*nfac
 
 		per=real(j)/real(npt)
 		if((per.ge.0.683).and.(perold.lt.0.683))then
      		e1=inth-(per-0.683)/(per-perold)*(inth-inthold)
      		e2=intl+(per-0.683)/(per-perold)*(intlold-intl)
      		errs(1)=e1-frsol
      		errs(2)=e2-frsol
c			write(6,*) "1 sigma:",frsol,errs(1),errs(2)
 		endif
 		if((per.ge.0.954).and.(perold.lt.0.954))then
      		e1=inth-(per-0.954)/(per-perold)*(inth-inthold)
      		e2=intl+(per-0.954)/(per-perold)*(intlold-intl)
      		errs(3)=e1-frsol
      		errs(4)=e2-frsol
c			write(6,*) "2 sigma:",frsol,errs(3),errs(4)
 		endif 
 		if((per.ge.0.9973).and.(perold.lt.0.9973))then
      		e1=inth-(per-0.9973)/(per-perold)*(inth-inthold)
      		e2=intl+(per-0.9973)/(per-perold)*(intlold-intl)
      		errs(5)=e1-frsol
      		errs(6)=e2-frsol
c			write(6,*) "3 sigma:",frsol,errs(5),errs(6)
 		endif     	
      	perold=per
      	intlold=intl
      	inthold=inth
c      	write(6,*) per,intl,inth,shgt
 13	  continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,nunit,i
      real sol(nfit),serr(nfit,2),Dpvary(nfit),toff,
     .  err(nfit),doe
      character*160 command
      
c      write(0,*) "----------------------------"
c      write(0,*) "Reading in starting solution"
c      write(0,*) "----------------------------"
      
      i=0 !initialize line counter
      toff=0.0 !initialize toff to zero
C     Start of loop to read in file
  10  read(nunit,500,end=11) command
 500  format(A160) !command much be contained within first 80 characters
        i=i+1 !increase line counter
        if(command(1:1).eq."#") then
            write(0,*) "Comment line on line ",i
        elseif(command(1:3).eq."SMA") then
            read(command(5:160),*) sol(1),serr(1,1),serr(1,2),Dpvary(1),
     .          err(1)
c           write(0,501) "SMA: ",sol(1),serr(1,1),serr(1,2),Dpvary(1),
c     .          err(1)
        elseif(command(1:3).eq."PMA") then
            read(command(5:160),*) sol(2),serr(2,1),serr(2,2),Dpvary(2),
     .          err(2)
c           write(0,501) "PMA: ",sol(2),serr(2,1),serr(2,2),Dpvary(2),
c     .          err(2)
        elseif(command(1:3).eq."SRA") then
            read(command(5:160),*) sol(3),serr(3,1),serr(3,2),Dpvary(3),
     .          err(3)
c           write(0,501) "SRA: ",sol(3),serr(3,1),serr(3,2),Dpvary(3),
c     .          err(3)
        elseif(command(1:3).eq."PRA") then
            read(command(5:160),*) sol(4),serr(4,1),serr(4,2),Dpvary(4),
     .          err(4)
c           write(0,501) "PRA: ",sol(4),serr(4,1),serr(4,2),Dpvary(4),
c     .          err(4)
        elseif(command(1:3).eq."PER") then
            read(command(5:160),*) sol(5),serr(5,1),serr(5,2),Dpvary(5),
     .          err(5)
c           write(0,501) "PER: ",sol(5),serr(5,1),serr(5,2),Dpvary(5),
c     .          err(5)
        elseif(command(1:3).eq."INC") then
            read(command(5:160),*) sol(6),serr(6,1),serr(6,2),Dpvary(6),
     .          err(6)
c           write(0,501) "INC: ",sol(6),serr(6,1),serr(6,2),Dpvary(6),
c     .          err(6)
        elseif(command(1:3).eq."EPO") then
            read(command(5:160),*) sol(7),serr(7,1),serr(7,2),Dpvary(7),
     .          err(7)
c           write(0,501) "EPO: ",sol(7),serr(7,1),serr(7,2),Dpvary(7),
c     .          err(7)
        elseif(command(1:3).eq."ZPT") then
            read(command(5:160),*) sol(8),serr(8,1),serr(8,2),Dpvary(8),
     .          err(8)
c           write(0,501) "ZPT: ",sol(8),serr(8,1),serr(8,2),Dpvary(8),
c     .          err(8)
        elseif(command(1:3).eq."ALB") then
            read(command(5:160),*) sol(9),serr(9,1),serr(9,2),Dpvary(9),
     .          err(9)
c           write(0,501) "ALB: ",sol(9),serr(9,1),serr(9,2),Dpvary(9),
c     .          err(9)
        elseif(command(1:3).eq."NL1") then
            read(command(5:160),*) sol(10),serr(10,1),serr(10,2),
     .          Dpvary(10),err(10)
c         write(0,501) "NL1: ",sol(10),serr(10,1),serr(10,2),Dpvary(10),
c     .          err(10)
        elseif(command(1:3).eq."NL2") then
            read(command(5:160),*) sol(11),serr(11,1),serr(11,2),
     .          Dpvary(11),err(11)
c         write(0,501) "NL2: ",sol(11),serr(11,1),serr(11,2),Dpvary(11),
c     .          err(11)
        elseif(command(1:3).eq."NL3") then
            read(command(5:160),*) sol(12),serr(12,1),serr(12,2),
     .          Dpvary(12),err(12)
c         write(0,501) "NL3: ",sol(12),serr(12,1),serr(12,2),Dpvary(12),
c     .          err(12)
        elseif(command(1:3).eq."NL4") then
            read(command(5:160),*) sol(13),serr(13,1),serr(13,2),
     .          Dpvary(13),err(13)
c         write(0,501) "NL4: ",sol(13),serr(13,1),serr(13,2),Dpvary(13),
c     .          err(13)
        elseif(command(1:3).eq."ECW") then
            read(command(5:160),*) sol(14),serr(14,1),serr(14,2),
     .          Dpvary(14),err(14)
c         write(0,501) "ECN: ",sol(14),serr(14,1),serr(14,2),Dpvary(14),
c     .          err(14)
        elseif(command(1:3).eq."ESW") then
            read(command(5:160),*) sol(15),serr(15,1),serr(15,2),
     .          Dpvary(15),err(15)
c         write(0,501) "WWW: ",sol(15),serr(15,1),serr(15,2),Dpvary(15),
c     .          err(15)
        elseif(command(1:3).eq."TED") then
            read(command(5:160),*) sol(16),serr(16,1),serr(16,2),
     .          Dpvary(16),err(16)
c         write(0,501) "TED: ",sol(16),serr(16,1),serr(16,2),Dpvary(16),
c     .          err(16)
        elseif(command(1:3).eq."VOF") then
            read(command(5:160),*) sol(17),serr(17,1),serr(17,2),
     .          Dpvary(17),err(17)
        elseif(command(1:3).eq."OFF") then
            read(command(5:160),*,err=12,end=12) toff,doe
            goto 13
 12         doe=0.0
            read(command(5:160),*) toff
 13     continue
c            write(0,501) "OFF: ",toff,doe
        endif
 501    format(A5,5(1X,PE17.10))
      goto 10 !loop back to read statement
 11   continue !end loop when EOF is reached
c      write(0,*) "----------------------------"
             
      return
      end 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readfreqs(nfreq,nmax,frsol)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfreq,nmax,nvar,i,j
      real frsol(nmax)
      double precision ztime
      
      read(11,*) ztime
      read(11,*) nvar,nfreq
      
      do 10 i=1,nfreq
        j=3*(i-1)+1 !read all variables in order.
      	read(11,*) frsol(j),frsol(j+1),frsol(j+2)
 10   continue
      return
      end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE avevar(data,n,ave,var)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n,m
      REAL ave,var,data(n)
C Given array data(1:n), returns its mean as ave and its variance as var.
      INTEGER j
      REAL s,ep
c      ave=0.0
c      m=min(1000,n)
c      do 11 j=1,m
c         ave=ave+data(j)
c 11   continue
c      ave=ave/m
      var=0.0
      ep=0.0
      do 12 j=1,n
         s=data(j)-ave
         ep=ep+s
         var=var+s*s
 12   continue
      var=(var-ep**2/n)/(n-1) !Corrected two-pass formula (14.1.8).
      return
      END 
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bindata(nbin,npt,pdata,bdatax,bdatay,datamin,datamax,
     .   bmax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nbin,npt,i,bin
      real pdata(npt),bdatax(nbin),bdatay(nbin),datamin,datamax,
     .   binsize,bmax,tsum
c      write(6,*) 'start'
C     get size of bin.
      binsize=(datamax-datamin)/real(nbin-1)
      
C     initalize bdata
      do 20 i=1,nbin
         bdatay(i)=0.
 20   continue

C     loop over data and place data into proper bin.
      do 10 i=1,npt
C        calculate bin number
         bin=int((pdata(i)-datamin)/binsize)+1
         if((bin.gt.0).and.(bin.le.nbin)) bdatay(bin)=bdatay(bin)+1.0
 10   continue

C     get the max value in a bin and assign bin value.
      tsum=0.
      bmax=0.
      do 30 i=1,nbin
         bmax=max(bdatay(i),bmax)
         bdatax(i)=datamin+real(i-1)*binsize !added -1
         tsum=tsum+bdatay(i)
 30   continue
c      write(6,*) "Sum:",tsum
c      write(6,*) 'end'
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine findrange(npt,pdata,datamin,datamax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i
      real pdata(npt),datamin,datamax
      
C     initialize datamin and datamax
      datamin=pdata(1)
      datamax=pdata(1)
      
      do 10 i=2,npt
         datamin=min(pdata(i),datamin)
         datamax=max(pdata(i),datamax)
 10   continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine windowsetup(xb1,xb2,yb1,yb2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      setup plotting area
      real xb1,xb2,yb1,yb2,xs1,xs2,
     .   ys1,ys2,pratio
      
C     set text size      
      call pgsch(1.0)
C     inquire window size
      call pgqvsz(2,xs1,xs2,ys1,ys2)
      pratio=(xs2-xs1)/(ys2-ys1)
C     if we have a landscape surface
      if(xs2-xs1.gt.ys2-ys1)then
C        start 10% from the edge
         yb1=(ys2-ys1)*0.10+ys1
C        use 60% of surface available
         yb2=(ys2-ys1)*0.70+yb1
C        now normalize to 1
         yb1=yb1/(ys2-ys1)         
         yb2=yb2/(ys2-ys1)
C        now use same size in other dimension      
         xb1=yb1
         xb2=yb1+(yb2-yb1)/pratio
      else
C        start 10% from the edge
         xb1=(xs2-xs1)*0.10+xs1
C        use 60% of surface available
         xb2=(xs2-xs1)*0.60+xb1
C        now normalize to 1
         xb1=xb1/(xs2-xs1)         
         xb2=xb2/(xs2-xs1)
C        now use same size in other dimension      
         yb1=xb1
         yb2=xb1+(xb2-xb1)*pratio
      endif
      
      return
      end 
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getdata(col,npt,dd)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer col,npt,i,j,k
      real dd(npt),dumr
      double precision dp,dpn
      
c      read(10,*,end=11) (dumr,j=1,col-1),dp
c      if(col.eq.7)then
c        dd(1)=0.0
c      else
c        dd(1)=real(dp)
c      endif
      
      k=0
      i=1
  9   continue
 10   read(10,*,end=11) (dumr,j=1,col-1),dpn
        if(col.eq.7)then
            dd(i)=real(dpn-dp)
        else
            dd(i)=real(dpn)
        endif
 
c         write(6,*) dd(i)
c         if((col.eq.6).and.(dd(i).gt.89.95)) goto 10
c         if(col.eq.8) dd(i)=dd(i)*1.0e6

c         k=k+1
c         if(k.lt.10) goto 10
c         k=0

         i=i+1
         goto 10
 11   continue
      
      npt=i-1
      
 500  format(200(E14.7,1X))
      goto 999
 901  write(6,*) "Error on line: ",i
      goto 999
      
 999  return
      end
      
c**********************************************************************
      subroutine rqsort(n,a,p)
c======================================================================
c     Return integer array p which indexes array a in increasing order.
c     Array a is not disturbed.  The Quicksort algorithm is used.
c
c     B. G. Knapp, 86/12/23
c
c     Reference: N. Wirth, Algorithms and Data Structures,
c     Prentice-Hall, 1986
c======================================================================
      implicit none

c     Input:
      integer   n
      real      a(n)

c     Output:
      integer   p(n)

c     Constants
      integer   LGN, Q
      parameter (LGN=32, Q=11)
c        (LGN = log base 2 of maximum n;
c         Q = smallest subfile to use quicksort on)

c     Local:
      real      x
      integer   stackl(LGN),stackr(LGN),s,t,l,m,r,i,j

c     Initialize the stack
      stackl(1)=1
      stackr(1)=n
      s=1

c     Initialize the pointer array
      do 1 i=1,n
         p(i)=i
    1 continue

    2 if (s.gt.0) then
         l=stackl(s)
         r=stackr(s)
         s=s-1

    3    if ((r-l).lt.Q) then

c           Use straight insertion
            do 6 i=l+1,r
               t = p(i)
               x = a(t)
               do 4 j=i-1,l,-1
                  if (a(p(j)).le.x) goto 5
                  p(j+1) = p(j)
    4          continue
               j=l-1
    5          p(j+1) = t
    6       continue
         else

c           Use quicksort, with pivot as median of a(l), a(m), a(r)
            m=(l+r)/2
            t=p(m)
            if (a(t).lt.a(p(l))) then
               p(m)=p(l)
               p(l)=t
               t=p(m)
            endif
            if (a(t).gt.a(p(r))) then
               p(m)=p(r)
               p(r)=t
               t=p(m)
               if (a(t).lt.a(p(l))) then
                  p(m)=p(l)
                  p(l)=t
                  t=p(m)
               endif
            endif

c           Partition
            x=a(t)
            i=l+1
            j=r-1
    7       if (i.le.j) then
    8          if (a(p(i)).lt.x) then
                  i=i+1
                  goto 8
               endif
    9          if (x.lt.a(p(j))) then
                  j=j-1
                  goto 9
               endif
               if (i.le.j) then
                  t=p(i)
                  p(i)=p(j)
                  p(j)=t
                  i=i+1
                  j=j-1
               endif
               goto 7
            endif

c           Stack the larger subfile
            s=s+1
            if ((j-l).gt.(r-i)) then
               stackl(s)=l
               stackr(s)=j
               l=i
            else
               stackl(s)=i
               stackr(s)=r
               r=j
            endif
            goto 3
         endif
         goto 2
      endif
      return
      end
