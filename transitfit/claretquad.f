      program claretquad
      implicit none
      integer nunit,i,iargc,nTeff,nFeH,teffidx,loggidx,fehidx
      real teff,logg,feh,teffs(79),loggs(11),fehs(19),diff,diffmin,c(2),
     .  pars(4),cmin(2),c1(4),c2(4),c3(4),c4(4),c5(4),c6(4),c7(4),c8(4)
      real teffmin,teffmax,loggmin,loggmax,fehmin,fehmax
      character*80 filename,dumc,cline
      data teffs/3500.,3750.,4000.,4250.,4500.,4750.,5000.,5250.,5500.,
     .  5750.,6000.,6250.,6500.,6750.,7000.,7250.,7500.,7750.,8000.,
     .  8250.,8500.,8750.,9000.,9250.,9500.,9750.,10000.,10250.,10500.,
     .  10750.,11000.,11250.,11500.,11750.,12000.,12250.,12500.,12750.,
     .  13000.,14000.,15000.,16000.,17000.,18000.,19000.,20000.,21000.,
     .  22000.,23000.,24000.,25000.,26000.,27000.,28000.,29000.,30000.,
     .  31000.,32000.,33000.,34000.,35000.,36000.,37000.,37500.,38000.,
     .  39000.,40000.,41000.,42000.,42500.,43000.,44000.,45000.,46000.,
     .  47000.,47500.,48000.,49000.,50000./
      data loggs/0.,0.5,1,1.5,2.,2.5,3,3.5,4.,4.5,5./
      data fehs/-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,-0.3,
     .  -0.2,-0.1,0.0,0.1,0.2,0.3,0.5,1.0/

      filename="claret-limbquad.txt"  !name of Claret tables
      nunit=10

      open(unit=nunit,file=filename,status='old',err=901) !open file

      if(iargc().lt.2) goto 902  !check command line arguments
      call getarg(1,cline)
      read(cline,*) teff  !read in Teff
      call getarg(2,cline)
      read(cline,*) logg  !read in log(g)
      if(iargc().ge.3) then
        call getarg(3,cline)
        read(cline,*) feh !read in [Fe/H]
      endif


      call locate(teffs,79,teff,teffidx)
      teffmin=teffs(teffidx)
      teffmax=teffs(teffidx+1)

      call locate(loggs,11,logg,loggidx)
      loggmin=loggs(loggidx)
      loggmax=loggs(loggidx+1)

      call locate(fehs,19,feh,fehidx)
      fehmin=fehs(fehidx)
      fehmax=fehs(fehidx+1)

c      write(0,*) teffmin,teffmax
c      write(0,*) loggmin,loggmax
c      write(0,*) fehmin,fehmax

      do 10 i=1,13
        read(nunit,*) dumc !skip header
 10   continue

      cmin(1)=0.5118 !set default output
      cmin(2)=0.0525
c      cmin(3)=0.4590
c      cmin(4)=-0.2727

 13   read(nunit,*,end=14) (pars(i),i=1,4),(c(i),i=1,2)

      if((pars(2).eq.teffmin).and.(pars(1).eq.loggmin)
     .  .and.(pars(3).eq.fehmin))then
               do 15 i=1,2
                    c1(i)=c(i)
 15             continue
      endif

      if((pars(2).eq.teffmin).and.(pars(1).eq.loggmin)
     .  .and.(pars(3).eq.fehmax))then
                 do 16 i=1,2
                     c2(i)=c(i)
  16             continue
      endif

      if((pars(2).eq.teffmax).and.(pars(1).eq.loggmin)
     .  .and.(pars(3).eq.fehmax))then
                 do 17 i=1,2
                     c3(i)=c(i)
  17             continue
      endif

      if((pars(2).eq.teffmax).and.(pars(1).eq.loggmin)
     .  .and.(pars(3).eq.fehmin))then
                 do 18 i=1,2
                     c4(i)=c(i)
  18             continue
      endif

      if((pars(2).eq.teffmin).and.(pars(1).eq.loggmax)
     .  .and.(pars(3).eq.fehmin))then
                 do 19 i=1,2
                     c5(i)=c(i)
  19             continue
      endif

      if((pars(2).eq.teffmin).and.(pars(1).eq.loggmax)
     .  .and.(pars(3).eq.fehmax))then
                 do 20 i=1,2
                     c6(i)=c(i)
  20             continue
      endif

      if((pars(2).eq.teffmax).and.(pars(1).eq.loggmax)
     .  .and.(pars(3).eq.fehmax))then
                 do 21 i=1,2
                     c7(i)=c(i)
  21             continue
      endif

      if((pars(2).eq.teffmax).and.(pars(1).eq.loggmax)
     .  .and.(pars(3).eq.fehmin))then
                 do 22 i=1,2
                     c8(i)=c(i)
  22             continue
      endif

      goto 13
 14   continue

        call trilinear(c1(1),c2(1),c3(1),c4(1),c5(1),c6(1),c7(1),c8(1),
     .      cmin(1),logg,teff,feh,loggmin,loggmax,teffmin,teffmax,
     .      fehmin,fehmax)

        call trilinear(c1(2),c2(2),c3(2),c4(2),c5(2),c6(2),c7(2),c8(2),
     .      cmin(2),logg,teff,feh,loggmin,loggmax,teffmin,teffmax,
     .      fehmin,fehmax)

c        call trilinear(c1(3),c2(3),c3(3),c4(3),c5(3),c6(3),c7(3),c8(3),
c     .      cmin(3),logg,teff,feh,loggmin,loggmax,teffmin,teffmax,
c     .      fehmin,fehmax)

c        call trilinear(c1(4),c2(4),c3(4),c4(4),c5(4),c6(4),c7(4),c8(4),
c     .      cmin(4),logg,teff,feh,loggmin,loggmax,teffmin,teffmax,
c     .      fehmin,fehmax)

      write(6,500) cmin !write answer to stdout
 500  format(4(F7.4,1X))
      close(nunit)  !close input file

      goto 999 !exit
 901  write(0,*) "Cannot open ",filename
      goto 999
 902  write(0,*) "Usage: claretquad Teff log(g) [Fe/H]"
      write(0,*) "  Teff is in K"
      write(0,*) "  log(g) is cgs"
      write(0,*) "  [Fe/H] is [m/H] and is optional (assumes solar)"
      goto 999
 999  end



      subroutine locate(array,length,x,idx)
      integer idx,length
      real x,array(length)
c     Given an array and a value x, returns idx
c     such that x is between array(idx) and array(idx+1).

      integer low,mid,up
      low=0                             !initialize limits
      up=length+1
 10   if(up-low.gt.1)then
         mid=(up+low)/2                 !compute midpoint
         if((array(length).ge.array(1)).eqv.(x.ge.array(mid)))then
            low=mid             !replace either lower or upper limit
         else
            up=mid
         endif
         goto 10
      endif
      if(x.eq.array(1))then
         idx=1
      else if(x.eq.array(length))then
         idx=length-1
      else
         idx=low
      endif
      return
      end


      subroutine trilinear(c1,c2,c3,c4,c5,c6,c7,c8,cout,logg,teff,feh,
     . loggmin,loggmax,teffmin,teffmax,fehmin,fehmax)
      implicit none
      real c1,c2,c3,c4,c5,c6,c7,c8,cout,t,u,v
      real loggmin,loggmax,teffmin,teffmax,fehmin,fehmax
      real logg,teff,feh

      t=(logg-loggmin)/(loggmax-loggmin) !normalize
      u=(teff-teffmin)/(teffmax-teffmin)
      v=(feh-fehmin)/(fehmax-fehmin)

      cout=(1-t)*(1-u)*(1-v)*c1+
     . t*(1-u)*(1-v)*c2+
     . t*u*(1-v)*c3+
     . (1-t)*u*(1-v)*c4+
     . (1-t)*(1-u)*v*c5+
     . t*(1-u)*v*c6+
     . t*u*v*c7+
     . (1-t)*u*v*c8


      return
      end
