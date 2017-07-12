      double precision function transn(nfit,solg,npt,time,flux,ferr,
     .  pcut)
      implicit none
      integer nfit,npt1,i,j,npt
      double precision solg(nfit),time1(1),exptime1(1),dtype1(1),
     .  tmodel1(1),tdepth,epo,per,ph1,ph2,std,stdev,flux(npt),ferr(npt),
     .  pcut(npt),tdur,transitdur,ph,time(npt),pmean

      tdur=transitdur(nfit,solg)/86400.0d0 !transit duration in days
      
      npt1=1
      time1(1)=solg(7)
      exptime1(1)=1765.5/86400.0d0
      dtype1(1)=0
      call transitmodel(npt1,time1,exptime1,dtype1,tmodel1,nfit,solg)
      tdepth=(1.0d0-10**((tmodel1(1)-solg(8))/-2.5d0))
c      write(0,*) "TDEPTH:",tdepth
      
c      write(6,*) int(tdepth+0.5),tdur
      
      epo=solg(7)
      per=solg(5)
      ph1=(tdur)/solg(5)/2.0d0
      ph2=1.0d0+(-tdur)/solg(5)/2.0d0
c      write(6,*) ph1,ph2
      
      j=0
      do 10 i=1,npt
        ph=(time(i)-epo)/per-int((time(i)-epo)/per)
        if(ph.lt.0.0) ph=ph+1.0
        if((ph.gt.ph1).and.(ph.lt.ph2))then
            j=j+1
            pcut(j)=flux(i)!(1.0d0-10**((flux(i)-solg(8))/-2.5d0))
c            write(6,*) time(i),pcut(j)
        endif
 10   continue
 
      std=stdev(j,pcut,pmean)
      
      transn=tdepth/std*sqrt(dble(npt-j))
      write(0,*) "STATS:", npt,std,j
      write(0,*) "S/N:",transn
      write(0,*) "STD: ",std
      write(0,*) "Tdepth: ",tdepth
      
      return
      end