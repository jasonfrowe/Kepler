      program mcmcphys
      implicit none
      integer nunit,iargc,nfit,nfitm,nplanet,flag,ng,i,nmax,j,npmax,k
      parameter(nfitm=108,nmax=1e6,npmax=10)
      double precision mstar,mstarerr,rstar,rstarerr,rhostar,nl(4),rchi,
     .  dil,voff,zpt,T0(npmax),Per(npmax),b(npmax),rprs(npmax),
     .  ecw(npmax),esw(npmax),Krv(npmax),ted(npmax),ell(npmax),
     .  alb(npmax),adrs,Psec,Pi,G,bb,incl,Mearth,Rearth,Mp,Rp
      character*80 filename,command
      
      Pi=acos(-1.d0)!define Pi and 2*Pi
      G=6.674d-11 !N m^2 kg^-2  Gravitation constamt
      Mearth=5.974d24 !kg Earth Mass
      Rearth=12756.0d0/2.0d0*1000.0d0
      
      if(iargc().lt.5) goto 901
      
      call getarg(1,filename)
      nunit=10
      open(unit=nunit,file=filename,status='old',err=902)
      read(nunit,*) nfit
      nplanet=(nfit-8)/10
      write(0,*) "nPlanet: ",nPlanet
      
      call getarg(2,command)
      read(command,*) mstar
      call getarg(3,command)
      read(command,*) mstarerr
      call getarg(4,command)
      read(command,*) rstar
      call getarg(5,command)
      read(command,*) rstarerr
      
      k=1
 10   read(nunit,*,end=11) rchi,flag,ng,rhostar,(nl(i),i=1,4),dil,voff,
     .  zpt,(T0(j),Per(j),b(j),rprs(j),ecw(j),esw(j),Krv(j),ted(j),
     .  ell(j),alb(j),j=1,nplanet)
     
        do 12 i=1,nplanet
            Psec=86400.0d0*Per(i)
            adrs=(1.0d3*rhostar*Psec*Psec*G/
     .          (3.0d0*Pi))**(1.0d0/3.0d0)
            bb=b(i)*b(i)
            incl=acos(b(i)/adrs)
 12     continue
 
        k=k+1
        if(k.gt.nmax) goto 903
      goto 10
 11   continue
      close(nunit)
      
      goto 999
 901  write(0,*) "Usage: mcmcphys <filename> <M> <Merr> <R> <Rerr>"
      write(0,*) "  filename = transitmcmc5 output"
      write(0,*) "  M, Merr = star mass and uncertainity (Msun)"
      write(0,*) "  R, Rerr = star radius and uncertainity (Rsun)"     
      goto 999
 902  write(0,*) "Cannot open ",filename
      goto 999
 903  write(0,*) "Segfault catch: increase nmax"
      goto 999
 999  end
