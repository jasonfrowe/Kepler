CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine autoexplore(nmax,npt,time,mag,merr,etime,tmodel,nfit,
     .  sol,serr,Dpvary,soltemp)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     A very simple and inefficient method to find the scale length 
C     vectors to find to amoeba
C     This subroutine scans sol+/-Dpvary to find the maximum 
C     likelihood.  If it is not bounded, then Dpvary is increased and
C     the search continues.  Once maximum is found, then Dpvary is
C     adjusted accordingly 
      implicit none
      integer npt,nfit,i,j,nplot,nmax
      parameter(nplot=100)
      double precision time(npt),mag(npt),merr(npt),
     .  etime(npt),sol(nfit),serr(nfit,2),Dpvary(nfit),
     .  soltemp(nfit),likelihood,vary,x1,x2,dnptm1,dnplotm1,chimin,
     .  chimax,chi(nplot),solmin,dsol,tmodel(nmax)
          
      dnptm1=dble(npt-1)
      dnplotm1=dble(nplot-1)    
          
      do 10 i=1,nfit !explore each parameter
        if(serr(i,2).ne.0.0d0)then !only explore fitted parameters
            do 11 j=1,nfit  !make a copy of the solution array
                soltemp(j)=sol(j)
 11         continue
            vary=Dpvary(i) !the exploration space
            x1=sol(i)-vary !bounds to explore
            x2=sol(i)+vary
            
            soltemp(i)=x1
            chi(1)=likelihood(npt,tmodel,time,mag,merr,etime,nfit,
     .          soltemp,serr)
            chimin=chi(1)
            chimax=chi(1)
            solmin=soltemp(i)
            do 12 j=2,nplot
                soltemp(i)=x1+(x2-x1)*dble(j-1)/dnplotm1
                chi(j)=likelihood(npt,tmodel,time,mag,merr,etime,nfit,
     .              soltemp,serr)
c                write(0,*) soltemp(i),chi(j)
c                read(5,*)
                if(chi(j).lt.chimin) then
                    chimin=chi(j)
                    solmin=soltemp(i)
                endif
                chimax=max(chi(j),chimax)
 12         continue
            dsol=sol(i)-solmin
            if((solmin.eq.x1).or.(solmin.eq.x2)) then
                write(0,*) "Edge! Increase Dpvary",i
                Dpvary(i)=2.0*Dpvary(i)
                write(0,*) sol(i),Dpvary(i),chimin
                if(Dpvary(i).eq.0.0d0)then
                    write(6,*)"Ahh.. zero. New Dpvary:"
                    read(5,*) Dpvary(i)
                endif
c                read(5,*)
                goto 11  !ugh.. horrible use of goto here
            endif
            Dpvary(i)=-dsol
        endif
 10   continue    
      
      return
      end