      program claretlimb
      implicit none
      integer nunit,i,iargc,nTeff,nFeH
      real teff,logg,feh,temps(61),fehs(19),diff,diffmin,c(4),pars(4),
     .  cmin(4)
      character*80 filename,dumc,cline
      data temps/3500.,3750.,4000.,4250.,4500.,4750.,5000.,5250.,5500.,
     .  5750.,6000.,6250.,6500.,6750.,7000.,7250.,7500.,7750.,8000.,
     .  8250.,8500.,8750.,9000.,9250.,9500.,9750.,10000.,10500.,11000.,
     .  11500.,12000.,12500.,13000.,14000.,15000.,16000.,17000.,18000.,
     .  19000.,20000.,21000.,22000.,23000.,24000.,25000.,26000.,27000.,
     .  28000.,29000.,30000.,31000.,32000.,33000.,34000.,35000.,37500.,
     .  40000.,42500.,45000.,47500.,50000./
      data fehs/-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,-0.3,
     .  -0.2,-0.1,0.0,0.1,0.2,0.3,0.5,1.0/

      cmin(1)=0.5118
      cmin(2)=0.0525
      cmin(3)=0.4590
      cmin(4)=-0.2727
      
c      filename="/Users/rowe/Documents/transitfit/Claret-limb.txt"
      filename="/Users/rowe/Documents/transitfit/Claret-limbquad.txt"
      nunit=10
      
      open(unit=nunit,file=filename,status='old',err=901)
      
      if(iargc().lt.2) goto 902
      call getarg(1,cline)
      read(cline,*) teff
      call getarg(2,cline)
      read(cline,*) logg
      if(iargc().ge.3) then
        call getarg(3,cline)
        read(cline,*) feh
      endif
      
      diffmin=abs(Teff-temps(1))
      nTeff=1
      do 11 i=2,61
        diff=abs(Teff-temps(i))
        if(diff.lt.diffmin)then
            nTeff=i
            diffmin=diff
        endif
 11   continue
      
      diffmin=abs(FeH-fehs(1))
      nFeH=1
      do 12 i=2,19
        diff=abs(FeH-fehs(i))
        if(diff.lt.diffmin)then
            nFeH=i
            diffmin=diff
        endif
 12   continue

c      write(0,*) Teff,Temps(nTeff),FeH,fehs(nFeH)   
      
c      do 10 i=1,34
c        read(nunit,*) dumc !skip header
c 10   continue
 
           
      diffmin=10.0
      cmin(1)=0.5118
      cmin(2)=0.0525
      cmin(3)=0.4590
      cmin(4)=-0.2727

c      cmin(1)=0.4035
c      cmin(2)=0.2622
 13   read(nunit,*,end=14) (pars(i),i=1,4),(c(i),i=1,2)
        if((pars(2).eq.Temps(nTeff)).and.(pars(3).eq.fehs(nFeH)))then
            diff=abs(logg-pars(1))
            if(diff.lt.diffmin)then
                diffmin=diff
                do 15 i=1,2
                    cmin(i)=c(i)
 15             continue
            endif
        endif
      goto 13
 14   continue
 
      write(6,500) (cmin(i),i=1,2)
 500  format(4(F7.4,1X))
      close(nunit)
      
      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 902  write(0,*) "Usage: claretlimb Teff log(g) [Fe/H]"
      write(0,*) "  Teff is in K"
      write(0,*) "  log(g) is cgs"
      write(0,*) "  [Fe/H] is [m/H] and is optional (assumes solar)"
      write(6,500) (cmin(i),i=1,4)
      goto 999
 999  end
