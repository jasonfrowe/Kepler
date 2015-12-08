      program kicpick
      implicit none
      integer nmax,kid,Teff,i,nunit,npt,nflag
      parameter(nmax=300000)
      integer kidobs(nmax)
      real kepmag,logg,rstar
      character*80 kicdatabase,kicobs
      kicdatabase="kic_colour.format.txt"
      kicobs="kicobs.list"

      nunit=10
      open(unit=nunit,file=kicobs,status='old',err=901)
      i=1
 10   read(nunit,*,end=11) kidobs(i)
         i=i+1
      goto 10
 11   continue
      npt=i-1
      write(0,*) "Number of KIDs: ",npt
      close(nunit)

      open(unit=nunit,file=kicdatabase,status='old',err=902)
 13      read(nunit,*,end=14) kid,kepmag,teff,logg,rstar
            if((kepmag.gt.0.0).and.(teff.lt.5000).and.(teff.gt.0).and.
     .       (logg.gt.0.0).and.(rstar.gt.0.0).and.(kepmag.lt.16.0)
     .       .and.(logg.gt.4.0))then
               nflag=0
               do 12 i=1,npt
                  if(kid.eq.kidobs(i))then
                     nflag=1
                  endif
 12            continue
               if(nflag.eq.0)then
                  write(6,*) kid,kepmag,teff,logg,rstar
               endif
            endif
         goto 13
 14      continue
      close(nunit)
      goto 999
901   write(0,*) "Cannot open ",kicobs
      goto 999
902   write(0,*) "Cannot open ",kicdatabase
      goto 999
999   end
