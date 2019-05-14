c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCE_INIT.FOR    (ErikSoft   28 February 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Calculates close-approach limits RCE (in AU) and physical radii RPHYS
c (in AU) for all objects, given their masses M, coordinates X, velocities
c V, densities RHO, and close-approach limits RCEH (in Hill radii).
c
c Also calculates the changeover distance RCRIT, used by the hybrid
c symplectic integrator. RCRIT is defined to be the larger of N1*HILL and
c N2*H*VMAX, where HILL is the Hill radius, H is the timestep, VMAX is the
c largest expected velocity of any body, and N1, N2 are parameters (see
c section 4.2 of Chambers 1999, Monthly Notices, vol 304, p793-799).
c
c N.B. Designed to use heliocentric coordinates, but should be adequate using
c ===  barycentric coordinates.
c
c------------------------------------------------------------------------------
c
c      subroutine mce_init (tstart,algor,h,jcen,rcen,rmax,cefac,nbod,
c     %  nbig,m,x,v,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile,rcritflag)
      subroutine mce_init (tstart,algor,h,jcen,rcen,rmax,cefac,nbod,
     %  nbig,m,x,v,s,rho,rceh,rphys,rce,rcrit,rcritflag)
c
      implicit none
      include 'mercury.inc'
c
      real*8 N2,THIRD
      parameter (N2=.4d0,THIRD=.3333333333333333d0)
c
c Input/Output
      integer nbod,nbig,algor,opt(8),rcritflag
      real*8 tstart,h,jcen(3),rcen,rmax,cefac,m(nbod),x(3,nbod)
      real*8 v(3,nbod),s(3,nbod),rho(nbod),rceh(nbod),rphys(nbod)
      real*8 rce(nbod),rcrit(nbod)
      character*8 id(nbod)
      character*80 outfile
c
c Local
      integer j
      real*8 a(NMAX),hill(NMAX),temp,amin,vmax,k_2,rhocgs,rcen_2
      character*80 header,c(NMAX)
      character*8 mio_re2c, mio_fl2c
c
c------------------------------------------------------------------------------
c
      rhocgs = AU * AU * AU * K2 / MSUN
      k_2 = 1.d0 / K2
      rcen_2 = 1.d0 / (rcen * rcen)
      amin = HUGE
c
c Calculate the Hill radii
      call mce_hill (nbod,m,x,v,hill,a)
c
c Determine the maximum close-encounter distances, and the physical radii
      temp = 2.25d0 * m(1) / PI
      do j = 2, nbod
        rce(j)   = hill(j) * rceh(j)
        rphys(j) = hill(j) / a(j) * (temp/rho(j))**THIRD
        amin = min (a(j), amin)
      end do
c
c If required, calculate the changeover distance used by hybrid algorithm
      if (rcritflag.eq.1) then
        vmax = sqrt (m(1) / amin)
        temp = N2 * h * vmax
        do j = 2, nbod
          rcrit(j) = max(hill(j) * cefac, temp)
        end do
      end if
c JR commented out.. don't need file-output
cc
cc Write list of object's identities to close-encounter output file
c      header(1:8)   = mio_fl2c (tstart)
c      header(9:16)  = mio_re2c (dble(nbig - 1),   0.d0, 11239423.99d0)
c      header(12:19) = mio_re2c (dble(nbod - nbig),0.d0, 11239423.99d0)
c      header(15:22) = mio_fl2c (m(1) * k_2)
c      header(23:30) = mio_fl2c (jcen(1) * rcen_2)
c      header(31:38) = mio_fl2c (jcen(2) * rcen_2 * rcen_2)
c      header(39:46) = mio_fl2c (jcen(3) * rcen_2 * rcen_2 * rcen_2)
c      header(47:54) = mio_fl2c (rcen)
c      header(55:62) = mio_fl2c (rmax)
cc
c      do j = 2, nbod
c        c(j)(1:8) = mio_re2c (dble(j - 1), 0.d0, 11239423.99d0)
c        c(j)(4:11) = id(j)
c        c(j)(12:19) = mio_fl2c (m(j) * k_2)
c        c(j)(20:27) = mio_fl2c (s(1,j) * k_2)
c        c(j)(28:35) = mio_fl2c (s(2,j) * k_2)
c        c(j)(36:43) = mio_fl2c (s(3,j) * k_2)
c        c(j)(44:51) = mio_fl2c (rho(j) / rhocgs)
c      end do
cc
cc Write compressed output to file
c  50  open (22, file=outfile, status='old', access='append', err=50)
c      write (22,'(a1,a2,i2,a62,i1)') char(12),'6a',algor,header(1:62),
c     %  opt(4)
c      do j = 2, nbod
c        write (22,'(a51)') c(j)(1:51)
c      end do
c      close (22)
cc
cc------------------------------------------------------------------------------
cc
      return
      end

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCE_HILL.FOR    (ErikSoft   4 October 2000)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Calculates the Hill radii for all objects given their masses, M,
c coordinates, X, and velocities, V; plus the mass of the central body, M(1)
c Where HILL = a * (m/3*m(1))^(1/3)
c
c If the orbit is hyperbolic or parabolic, the Hill radius is calculated using:
c       HILL = r * (m/3*m(1))^(1/3)
c where R is the current distance from the central body.
c
c The routine also gives the semi-major axis, A, of each object's orbit.
c
c N.B. Designed to use heliocentric coordinates, but should be adequate using
c ===  barycentric coordinates.
c
c------------------------------------------------------------------------------
c
      subroutine mce_hill (nbod,m,x,v,hill,a)
c
      implicit none
      include 'mercury.inc'
      real*8 THIRD
      parameter (THIRD = .3333333333333333d0)
c
c Input/Output
      integer nbod
      real*8 m(nbod),x(3,nbod),v(3,nbod),hill(nbod),a(nbod)
c
c Local
      integer j
      real*8 r, v2, gm
c
c------------------------------------------------------------------------------
c
      do j = 2, nbod
        gm = m(1) + m(j)
        call mco_x2a (gm,x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),a(j),
     %    r,v2)
c If orbit is hyperbolic, use the distance rather than the semi-major axis
        if (a(j).le.0) a(j) = r
        hill(j) = a(j) * (THIRD * m(j) / m(1))**THIRD
      end do
c
c------------------------------------------------------------------------------
c
      return
      end

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_X2A.FOR    (ErikSoft   4 October 2000)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Calculates an object's orbital semi-major axis given its Cartesian coords.
c
c------------------------------------------------------------------------------
c
      subroutine mco_x2a (gm,x,y,z,u,v,w,a,r,v2)
c
      implicit none
c
c Input/Output
      real*8 gm,x,y,z,u,v,w,a,r,v2
c
c------------------------------------------------------------------------------
c
      r  = sqrt(x * x  +  y * y  +  z * z)
      v2 =      u * u  +  v * v  +  w * w
      a  = gm * r / (2.d0 * gm  -  r * v2)
c
c------------------------------------------------------------------------------
c
      return
      end

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_H2DH.FOR    (ErikSoft   2 March 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Convert coordinates with respect to the central body to democratic
c heliocentric coordinates.
c
c------------------------------------------------------------------------------
c
C      subroutine mco_h2dh (time,jcen,nbod,nbig,h,m,xh,vh,x,v,ngf,ngflag,
C     %  opt)
      subroutine mco_h2dh (time,jcen,nbod,nbig,h,m,xh,vh,x,v)
C      JR remove ngf,ngflag and opt, because they are not used
c
      implicit none
c
c Input/Output
      integer nbod,nbig
      real*8 time,jcen(3),h,m(nbod),xh(3,nbod),vh(3,nbod),x(3,nbod)
      real*8 v(3,nbod)
c
c Local
      integer j
      real*8 mtot,temp,mvsum(3)
c
c------------------------------------------------------------------------------
c
      mtot = 0.d0
      mvsum(1) = 0.d0
      mvsum(2) = 0.d0
      mvsum(3) = 0.d0
c
      do j = 2, nbod
        x(1,j) = xh(1,j)
        x(2,j) = xh(2,j)
        x(3,j) = xh(3,j)
        mtot = mtot + m(j)
        mvsum(1) = mvsum(1)  +  m(j) * vh(1,j)
        mvsum(2) = mvsum(2)  +  m(j) * vh(2,j)
        mvsum(3) = mvsum(3)  +  m(j) * vh(3,j)
      end do
c
      temp = 1.d0 / (m(1) + mtot)
c
      mvsum(1) = temp * mvsum(1)
      mvsum(2) = temp * mvsum(2)
      mvsum(3) = temp * mvsum(3)
c
      do j = 2, nbod
        v(1,j) = vh(1,j) - mvsum(1)
        v(2,j) = vh(2,j) - mvsum(2)
        v(3,j) = vh(3,j) - mvsum(3)
      end do
c
c------------------------------------------------------------------------------
c
      return
      end

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MDT_HY.FOR    (ErikSoft   2 March 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Integrates NBOD bodies (of which NBIG are Big) for one timestep H
c using a second-order hybrid-symplectic integrator algorithm
c
c DTFLAG = 0 implies first ever call to this subroutine,
c        = 1 implies first call since number/masses of objects changed.
c        = 2 normal call
c
c N.B. Input/output must be in democratic heliocentric coordinates.
c ===
c
c------------------------------------------------------------------------------
c
c      subroutine mdt_hy (time,tstart,h0,tol,rmax,en,am,jcen,rcen,nbod,
c     %  nbig,m,x,v,s,rphys,rcrit,rce,stat,id,ngf,algor,opt,dtflag,
c     %  ngflag,opflag,colflag,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,
c     %  outfile,mem,lmem)
      subroutine mdt_hy (time,tstart,h0,tol,rmax,en,am,jcen,rcen,nbod,
     %  nbig,m,x,v,s,rphys,rcrit,rce,stat,algor,opt,dtflag,
     %  ngflag,opflag,colflag,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo)
c
      implicit none
      include 'mercury.inc'
c
c Input/Output
      integer nbod,nbig,stat(nbod),algor,opt(8),dtflag,ngflag,opflag
      integer colflag,nclo,iclo(CMAX),jclo(CMAX)
      real*8 time,tstart,h0,tol,rmax,en(3),am(3),jcen(3),rcen
      real*8 m(nbod),x(3,nbod),v(3,nbod),s(3,nbod),rphys(nbod)
      real*8 rce(nbod),rcrit(nbod),tclo(CMAX),dclo(CMAX)
      real*8 ixvclo(6,CMAX),jxvclo(6,CMAX)
C      character*80 outfile(3),mem(NMESS)
c
c Local
      integer j,nce,ice(NMAX),jce(NMAX),ce(NMAX),iflag
      real*8 a(3,NMAX),hby2,hrec,x0(3,NMAX),v0(3,NMAX),mvsum(3),temp
      real*8 angf(3,NMAX),ausr(3,NMAX)
      external mfo_hkce
c
c------------------------------------------------------------------------------
c
      save a, hrec, angf, ausr
      hby2 = h0 * .5d0
      nclo = 0
      colflag = 0
c
c If accelerations from previous call are not valid, calculate them now
      if (dtflag.ne.2) then
        if (dtflag.eq.0) hrec = h0
        call mfo_hy (jcen,nbod,nbig,m,x,rcrit,a,stat)
        dtflag = 2
        do j = 2, nbod
          angf(1,j) = 0.d0
          angf(2,j) = 0.d0
          angf(3,j) = 0.d0
          ausr(1,j) = 0.d0
          ausr(2,j) = 0.d0
          ausr(3,j) = 0.d0
        end do
c If required, apply non-gravitational and user-defined forces
c        if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
c        if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,x,v,angf,ngf)
      end if
c
c Advance interaction Hamiltonian for H/2
      do j = 2, nbod
        v(1,j) = v(1,j)  +  hby2 * (angf(1,j) + ausr(1,j) + a(1,j))
        v(2,j) = v(2,j)  +  hby2 * (angf(2,j) + ausr(2,j) + a(2,j))
        v(3,j) = v(3,j)  +  hby2 * (angf(3,j) + ausr(3,j) + a(3,j))
      end do
c
c Advance solar Hamiltonian for H/2
      mvsum(1) = 0.d0
      mvsum(2) = 0.d0
      mvsum(3) = 0.d0
      do j = 2, nbod
        mvsum(1) = mvsum(1)  +  m(j) * v(1,j)
        mvsum(2) = mvsum(2)  +  m(j) * v(2,j)
        mvsum(3) = mvsum(3)  +  m(j) * v(3,j)
      end do
c
      temp = hby2 / m(1)
      mvsum(1) = temp * mvsum(1)
      mvsum(2) = temp * mvsum(2)
      mvsum(3) = temp * mvsum(3)
      do j = 2, nbod
        x(1,j) = x(1,j)  +  mvsum(1)
        x(2,j) = x(2,j)  +  mvsum(2)
        x(3,j) = x(3,j)  +  mvsum(3)
      end do
c
c Save the current coordinates and velocities
      call mco_iden (time,jcen,nbod,nbig,h0,m,x,v,x0,v0)
c
c Advance H_K for H
      do j = 2, nbod
        iflag = 0
        call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),
     %    v(3,j),h0,iflag)
      end do
c
c Check whether any object separations were < R_CRIT whilst advancing H_K
      call mce_snif (h0,2,nbod,nbig,x0,v0,x,v,rcrit,ce,nce,ice,jce)
c
c If objects had close encounters, advance H_K using Bulirsch-Stoer instead
      if (nce.gt.0) then
c         write(0,*) "close: mdt_hy"
        do j = 2, nbod
          if (ce(j).ne.0) then
            x(1,j) = x0(1,j)
            x(2,j) = x0(2,j)
            x(3,j) = x0(3,j)
            v(1,j) = v0(1,j)
            v(2,j) = v0(2,j)
            v(3,j) = v0(3,j)
          end if
        end do
        call mdt_hkce (time,tstart,h0,hrec,tol,rmax,en(3),jcen,rcen,
     %    nbod,nbig,m,x,v,s,rphys,rcrit,rce,stat,algor,opt,
     %    ngflag,colflag,ce,nce,ice,jce,nclo,iclo,jclo,dclo,tclo,ixvclo,
     %    jxvclo,mfo_hkce)
      end if
c
c Advance solar Hamiltonian for H/2
      mvsum(1) = 0.d0
      mvsum(2) = 0.d0
      mvsum(3) = 0.d0
      do j = 2, nbod
        mvsum(1) = mvsum(1)  +  m(j) * v(1,j)
        mvsum(2) = mvsum(2)  +  m(j) * v(2,j)
        mvsum(3) = mvsum(3)  +  m(j) * v(3,j)
      end do
c
      temp = hby2 / m(1)
      mvsum(1) = temp * mvsum(1)
      mvsum(2) = temp * mvsum(2)
      mvsum(3) = temp * mvsum(3)
      do j = 2, nbod
        x(1,j) = x(1,j)  +  mvsum(1)
        x(2,j) = x(2,j)  +  mvsum(2)
        x(3,j) = x(3,j)  +  mvsum(3)
      end do
c
c Advance interaction Hamiltonian for H/2
      call mfo_hy (jcen,nbod,nbig,m,x,rcrit,a,stat)
c      if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
C Removed by JR
C      if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,x,v,angf,ngf)
c
      do j = 2, nbod
        v(1,j) = v(1,j)  +  hby2 * (angf(1,j) + ausr(1,j) + a(1,j))
        v(2,j) = v(2,j)  +  hby2 * (angf(2,j) + ausr(2,j) + a(2,j))
        v(3,j) = v(3,j)  +  hby2 * (angf(3,j) + ausr(3,j) + a(3,j))
      end do
c
c------------------------------------------------------------------------------
c
      return
      end

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_IDEN.FOR    (ErikSoft   2 November 2000)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Makes a new copy of a set of coordinates.
c
c------------------------------------------------------------------------------
c
      subroutine mco_iden (time,jcen,nbod,nbig,h,m,xh,vh,x,v)
c
      implicit none
c
c Input/Output
      integer nbod,nbig
      real*8 time,jcen(3),h,m(nbod),x(3,nbod),v(3,nbod),xh(3,nbod)
      real*8 vh(3,nbod)
c
c Local
      integer j
c
c------------------------------------------------------------------------------
c
      do j = 1, nbod
        x(1,j) = xh(1,j)
        x(2,j) = xh(2,j)
        x(3,j) = xh(3,j)
        v(1,j) = vh(1,j)
        v(2,j) = vh(2,j)
        v(3,j) = vh(3,j)
      enddo
c
c------------------------------------------------------------------------------
c
      return
      end

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MDT_HKCE.FOR    (ErikSoft   1 March 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Integrates NBOD bodies (of which NBIG are Big) for one timestep H under
c the Hamiltonian H_K, including close-encounter terms.
c
c------------------------------------------------------------------------------
c
c      subroutine mdt_hkce (time,tstart,h0,hrec,tol,rmax,elost,jcen,
c     %  rcen,nbod,nbig,m,x,v,s,rphy,rcrit,rce,stat,id,ngf,algor,opt,
c     %  ngflag,colflag,ce,nce,ice,jce,nclo,iclo,jclo,dclo,tclo,ixvclo,
c     %  jxvclo,outfile,mem,lmem)
      subroutine mdt_hkce (time,tstart,h0,hrec,tol,rmax,elost,jcen,
     %  rcen,nbod,nbig,m,x,v,s,rphy,rcrit,rce,stat,algor,opt,
     %  ngflag,colflag,ce,nce,ice,jce,nclo,iclo,jclo,dclo,tclo,ixvclo,
     %  jxvclo)
c
      implicit none
      include 'mercury.inc'
c
c Input/Output
      integer nbod,nbig,nce,ice(nce),jce(nce),stat(nbod),ngflag,ce(nbod)
      integer algor,opt(8),colflag,nclo,iclo(CMAX)
      integer jclo(CMAX)
      real*8 time,tstart,h0,hrec,tol,rmax,elost,jcen(3),rcen
      real*8 m(nbod),x(3,nbod),v(3,nbod),s(3,nbod)
      real*8 rce(nbod),rphy(nbod),rcrit(nbod)
      real*8 tclo(CMAX),dclo(CMAX),ixvclo(6,CMAX),jxvclo(6,CMAX)
C      character*80 outfile(3),mem(NMESS)
C      character*8 id(nbod)

c
c Local
      integer iback(NMAX),index(NMAX),ibs(NMAX),jbs(NMAX),nclo_old
      integer i,j,k,nbs,nbsbig,statbs(NMAX)
      integer nhit,ihit(CMAX),jhit(CMAX),chit(CMAX),nowflag,dtflag
      real*8 tlocal,hlocal,hdid,tmp0
      real*8 mbs(NMAX),xbs(3,NMAX),vbs(3,NMAX),sbs(3,NMAX)
      real*8 rcritbs(NMAX),rcebs(NMAX),rphybs(NMAX)
      real*8 ngfbs(4,NMAX),x0(3,NMAX),v0(3,NMAX)
      real*8 thit(CMAX),dhit(CMAX),thit1,temp
      character*8 idbs(NMAX)
c
c------------------------------------------------------------------------------
c
c N.B. Don't set nclo to zero!!
      nbs = 1
      nbsbig = 0
      mbs(1) = m(1)
      if (algor.eq.11) mbs(1) = m(1) + m(2)
      sbs(1,1) = s(1,1)
      sbs(2,1) = s(2,1)
      sbs(3,1) = s(3,1)
c
c Put data for close-encounter bodies into local arrays for use with BS routine
      do j = 2, nbod
        if (ce(j).ne.0) then
          nbs = nbs + 1
          if (j.le.nbig) nbsbig = nbs
          mbs(nbs)   = m(j)
          xbs(1,nbs) = x(1,j)
          xbs(2,nbs) = x(2,j)
          xbs(3,nbs) = x(3,j)
          vbs(1,nbs) = v(1,j)
          vbs(2,nbs) = v(2,j)
          vbs(3,nbs) = v(3,j)
          sbs(1,nbs) = s(1,j)
          sbs(2,nbs) = s(2,j)
          sbs(3,nbs) = s(3,j)
          rcebs(nbs) = rce(j)
          rphybs(nbs) = rphy(j)
          statbs(nbs) = stat(j)
          rcritbs(nbs) = rcrit(j)
cJR          idbs(nbs) = id(j)
          index(nbs) = j
          iback(j) = nbs
        end if
      end do
c
      do k = 1, nce
        ibs(k) = iback(ice(k))
        jbs(k) = iback(jce(k))
      end do
c
      tlocal = 0.d0
      hlocal = sign(hrec,h0)
c
c Begin the Bulirsch-Stoer integration
  50  continue
        tmp0 = abs(h0) - abs(tlocal)
        hrec = hlocal
        if (abs(hlocal).gt.tmp0) hlocal = sign (tmp0, h0)
c
c Save old coordinates and integrate
C edit by JR to remove ngf,ngflag,opt to mco_iden
        call mco_iden (time,jcen,nbs,0,h0,mbs,xbs,vbs,x0,v0)
        call mdt_bs2 (time,hlocal,hdid,tol,jcen,nbs,nbsbig,mbs,xbs,vbs,
     %    sbs,rphybs,rcritbs,ngfbs,statbs,dtflag,ngflag,opt,nce,
     %    ibs,jbs)
        tlocal = tlocal + hdid
c
c Check for close-encounter minima
        nclo_old = nclo
        temp = time + tlocal
        call mce_stat (temp,hdid,rcen,nbs,nbsbig,mbs,x0,v0,xbs,vbs,
     %    rcebs,rphybs,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,nhit,ihit,
     %    jhit,chit,dhit,thit,thit1,nowflag,statbs)
c
c If collisions occurred, resolve the collision and return a flag
        if (nhit.gt.0.and.opt(2).ne.0) then
          do k = 1, nhit
            if (chit(k).eq.1) then
              i = ihit(k)
              j = jhit(k)
              call mce_coll (thit(k),tstart,elost,jcen,i,j,nbs,nbsbig,
     %          mbs,xbs,vbs,sbs,rphybs,statbs)
              colflag = colflag + 1
            end if
          end do
        end if
c
c If necessary, continue integrating objects undergoing close encounters
      if ((tlocal - h0)*h0.lt.0) goto 50
c
c Return data for the close-encounter objects to global arrays
      do k = 2, nbs
        j = index(k)
        m(j)   = mbs(k)
        x(1,j) = xbs(1,k)
        x(2,j) = xbs(2,k)
        x(3,j) = xbs(3,k)
        v(1,j) = vbs(1,k)
        v(2,j) = vbs(2,k)
        v(3,j) = vbs(3,k)
        s(1,j) = sbs(1,k)
        s(2,j) = sbs(2,k)
        s(3,j) = sbs(3,k)
        stat(j) = statbs(k)
      end do
      do k = 1, nclo
        iclo(k) = index(iclo(k))
        jclo(k) = index(jclo(k))
      end do
c
c------------------------------------------------------------------------------
c
      return
      end

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCE_COLL.FOR    (ErikSoft   2 October 2000)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Resolves a collision between two objects, using the collision model chosen
c by the user. Also writes a message to the information file, and updates the
c value of ELOST, the change in energy due to collisions and ejections.
c
c N.B. All coordinates and velocities must be with respect to central body.
c ===
c
c------------------------------------------------------------------------------
c
      subroutine mce_coll (time,tstart,elost,jcen,i,j,nbod,nbig,m,xh,
     %  vh,s,rphys,stat)
c
      implicit none
      include 'mercury.inc'
c
c Input/Output
      integer i,j,nbod,nbig,stat(nbod)
      real*8 time,tstart,elost,jcen(3)
      real*8 m(nbod),xh(3,nbod),vh(3,nbod),s(3,nbod),rphys(nbod)
c      character*80 outfile,mem(NMESS)
c      character*8 id(nbod)
c
c Local
      integer year,month,itmp
      real*8 t1
      character*38 flost,fcol
      character*6 tstring
c
c------------------------------------------------------------------------------
c
c If two bodies collided, check that the less massive one is removed
c (unless the more massive one is a Small body)
      if (i.ne.0) then
        if (m(j).gt.m(i).and.j.le.nbig) then
          itmp = i
          i = j
          j = itmp
        end if
      end if
Cc
Cc Write message to info file (I=0 implies collision with the central body)
C  10  open (23, file=outfile, status='old', access='append', err=10)
Cc
C      if (opt(3).eq.1) then
C        call mio_jd2y (time,year,month,t1)
C        if (i.eq.0) then
C          flost = '(1x,a8,a,i10,1x,i2,1x,f8.5)'
C          write (23,flost) id(j),mem(67)(1:lmem(67)),year,month,t1
C        else
C          fcol  = '(1x,a8,a,a8,a,i10,1x,i2,1x,f4.1)'
C          write (23,fcol) id(i),mem(69)(1:lmem(69)),id(j),
C     %      mem(71)(1:lmem(71)),year,month,t1
C        end if
C      else
C        if (opt(3).eq.3) then
C          t1 = (time - tstart) / 365.25d0
C          tstring = mem(2)
C          flost = '(1x,a8,a,f18.7,a)'
C          fcol  = '(1x,a8,a,a8,a,1x,f14.3,a)'
C        else
C          if (opt(3).eq.0) t1 = time
C          if (opt(3).eq.2) t1 = time - tstart
C          tstring = mem(1)
C          flost = '(1x,a8,a,f18.5,a)'
C          fcol  = '(1x,a8,a,a8,a,1x,f14.1,a)'
C        end if
C        if (i.eq.0.or.i.eq.1) then
C          write (23,flost) id(j),mem(67)(1:lmem(67)),t1,tstring
C        else
C          write (23,fcol) id(i),mem(69)(1:lmem(69)),id(j),
C     %      mem(71)(1:lmem(71)),t1,tstring
C        end if
C      end if
C      close (23)
c
c Do the collision (inelastic merger)
      call mce_merg (jcen,i,j,nbod,nbig,m,xh,vh,s,stat,elost)
c
c------------------------------------------------------------------------------
c
      return
      end

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c    MCE_STAT.FOR    (ErikSoft   1 March 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Returns details of all close-encounter minima involving at least one Big
c body during a timestep. The routine estimates minima using the initial
c and final coordinates X(0),X(1) and velocities V(0),V(1) of the step, and
c the stepsize H.
c
c  ICLO, JCLO contain the indices of the objects
c  DCLO is their minimum separation
c  TCLO is the time of closest approach relative to current time
c
c The routine also checks for collisions/near misses given the physical radii
c RPHYS, and returns the time THIT of the collision/near miss closest to the
c start of the timestep, and the identities IHIT and JHIT of the objects
c involved.
c
c  NHIT = +1 implies a collision
c         -1    "    a near miss
c
c N.B. All coordinates & velocities must be with respect to the central body!!
c ===
c------------------------------------------------------------------------------
c
      subroutine mce_stat (time,h,rcen,nbod,nbig,m,x0,v0,x1,v1,rce,
     %  rphys,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,nhit,ihit,jhit,
     %  chit,dhit,thit,thit1,nowflag,stat)
c
      implicit none
      include 'mercury.inc'
c
c Input/Output
      integer nbod,nbig,stat(nbod),nowflag
      integer nclo,iclo(CMAX),jclo(CMAX)
      integer nhit,ihit(CMAX),jhit(CMAX),chit(CMAX)
      real*8 time,h,rcen,m(nbod),x0(3,nbod),v0(3,nbod)
      real*8 x1(3,nbod),v1(3,nbod),rce(nbod),rphys(nbod)
      real*8 dclo(CMAX),tclo(CMAX),thit(CMAX),dhit(CMAX),thit1
      real*8 ixvclo(6,CMAX),jxvclo(6,CMAX)
c      character*80 outfile,mem(NMESS)
c
c Local
      integer i,j
      real*8 d0,d1,d0t,d1t,hm1,tmp0,tmp1
      real*8 dx0,dy0,dz0,du0,dv0,dw0,dx1,dy1,dz1,du1,dv1,dw1
      real*8 xmin(NMAX),xmax(NMAX),ymin(NMAX),ymax(NMAX)
      real*8 d2min,d2ce,d2near,d2hit,temp,tmin
c
c------------------------------------------------------------------------------
c
      nhit = 0
      thit1 = sign(HUGE, h)
      hm1 = 1.d0 / h
c
c Calculate maximum and minimum values of x and y coords for each object
      call mce_box (nbod,h,x0,v0,x1,v1,xmin,xmax,ymin,ymax)
c
c Adjust values by the maximum close-encounter radius plus a fudge factor
      do j = 2, nbod
        temp = rce(j) * 1.2d0
        xmin(j) = xmin(j)  -  temp
        xmax(j) = xmax(j)  +  temp
        ymin(j) = ymin(j)  -  temp
        ymax(j) = ymax(j)  +  temp
      end do
c
c Check for close encounters between each pair of objects
      do i = 2, nbig
        do j = i + 1, nbod
          if (   xmax(i).ge.xmin(j).and.xmax(j).ge.xmin(i)
     %      .and.ymax(i).ge.ymin(j).and.ymax(j).ge.ymin(i)
     %      .and.stat(i).ge.0.and.stat(j).ge.0) then
c
c If the X-Y boxes for this pair overlap, check circumstances more closely
            dx0 = x0(1,i) - x0(1,j)
            dy0 = x0(2,i) - x0(2,j)
            dz0 = x0(3,i) - x0(3,j)
            du0 = v0(1,i) - v0(1,j)
            dv0 = v0(2,i) - v0(2,j)
            dw0 = v0(3,i) - v0(3,j)
            d0t = (dx0*du0 + dy0*dv0 + dz0*dw0) * 2.d0
c
            dx1 = x1(1,i) - x1(1,j)
            dy1 = x1(2,i) - x1(2,j)
            dz1 = x1(3,i) - x1(3,j)
            du1 = v1(1,i) - v1(1,j)
            dv1 = v1(2,i) - v1(2,j)
            dw1 = v1(3,i) - v1(3,j)
            d1t = (dx1*du1 + dy1*dv1 + dz1*dw1) * 2.d0
c
c Estimate minimum separation during the time interval, using interpolation
            d0 = dx0*dx0 + dy0*dy0 + dz0*dz0
            d1 = dx1*dx1 + dy1*dy1 + dz1*dz1
            call mce_min (d0,d1,d0t,d1t,h,d2min,tmin)
            d2ce  = max (rce(i), rce(j))
            d2hit = rphys(i) + rphys(j)
            d2ce   = d2ce  * d2ce
            d2hit  = d2hit * d2hit
            d2near = d2hit * 4.d0
c
c If the minimum separation qualifies as an encounter or if a collision
c is in progress, store details
            if ((d2min.le.d2ce.and.d0t*h.le.0.and.d1t*h.ge.0)
     %        .or.(d2min.le.d2hit)) then
              nclo = nclo + 1
              if (nclo.gt.CMAX) then
                write(0,*) "nclo>nmax"
c 230            open (23,file=outfile,status='old',access='append',
c     %            err=230)
c                write (23,'(/,2a,/,a)') mem(121)(1:lmem(121)),
c     %            mem(132)(1:lmem(132)),mem(82)(1:lmem(82))
c                close (23)
              else
                tclo(nclo) = tmin + time
                dclo(nclo) = sqrt (max(0.d0,d2min))
                iclo(nclo) = i
                jclo(nclo) = j
c
c Make sure the more massive body is listed first
                if (m(j).gt.m(i).and.j.le.nbig) then
                  iclo(nclo) = j
                  jclo(nclo) = i
                end if
c
c Make linear interpolation to get coordinates at time of closest approach
                tmp0 = 1.d0 + tmin*hm1
                tmp1 = -tmin*hm1
                ixvclo(1,nclo) = tmp0 * x0(1,i)  +  tmp1 * x1(1,i)
                ixvclo(2,nclo) = tmp0 * x0(2,i)  +  tmp1 * x1(2,i)
                ixvclo(3,nclo) = tmp0 * x0(3,i)  +  tmp1 * x1(3,i)
                ixvclo(4,nclo) = tmp0 * v0(1,i)  +  tmp1 * v1(1,i)
                ixvclo(5,nclo) = tmp0 * v0(2,i)  +  tmp1 * v1(2,i)
                ixvclo(6,nclo) = tmp0 * v0(3,i)  +  tmp1 * v1(3,i)
                jxvclo(1,nclo) = tmp0 * x0(1,j)  +  tmp1 * x1(1,j)
                jxvclo(2,nclo) = tmp0 * x0(2,j)  +  tmp1 * x1(2,j)
                jxvclo(3,nclo) = tmp0 * x0(3,j)  +  tmp1 * x1(3,j)
                jxvclo(4,nclo) = tmp0 * v0(1,j)  +  tmp1 * v1(1,j)
                jxvclo(5,nclo) = tmp0 * v0(2,j)  +  tmp1 * v1(2,j)
                jxvclo(6,nclo) = tmp0 * v0(3,j)  +  tmp1 * v1(3,j)
              end if
            end if
c
c Check for a near miss or collision
            if (d2min.le.d2near) then
              nhit = nhit + 1
              ihit(nhit) = i
              jhit(nhit) = j
              thit(nhit) = tmin + time
              dhit(nhit) = sqrt(d2min)
              chit(nhit) = -1
              if (d2min.le.d2hit) chit(nhit) = 1
c
c Make sure the more massive body is listed first
              if (m(jhit(nhit)).gt.m(ihit(nhit)).and.j.le.nbig) then
                ihit(nhit) = j
                jhit(nhit) = i
              end if
c
c Is this the collision closest to the start of the time step?
              if ((tmin-thit1)*h.lt.0) then
                thit1 = tmin
                nowflag = 0
                if (d1.le.d2hit) nowflag = 1
              end if
            end if
          end if
c
c Move on to the next pair of objects
        end do
      end do
c
c------------------------------------------------------------------------------
c
      return
      end

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCE_MIN.FOR    (ErikSoft  1 December 1998)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Calculates minimum value of a quantity D, within an interval H, given initial
c and final values D0, D1, and their derivatives D0T, D1T, using third-order
c (i.e. cubic) interpolation.
c
c Also calculates the value of the independent variable T at which D is a
c minimum, with respect to the epoch of D1.
c
c N.B. The routine assumes that only one minimum is present in the interval H.
c ===
c------------------------------------------------------------------------------
c
      subroutine mce_min (d0,d1,d0t,d1t,h,d2min,tmin)
c
      implicit none
c
c Input/Output
      real*8 d0,d1,d0t,d1t,h,d2min,tmin
c
c Local
      real*8 a,b,c,temp,tau
c
c------------------------------------------------------------------------------
c
      if (d0t*h.gt.0.or.d1t*h.lt.0) then
        if (d0.le.d1) then
          d2min = d0
          tmin = -h
        else
          d2min = d1
          tmin = 0.d0
        end if
      else
        temp = 6.d0*(d0 - d1)
        a = temp + 3.d0*h*(d0t + d1t)
        b = temp + 2.d0*h*(d0t + 2.d0*d1t)
        c = h * d1t
c
        temp =-.5d0*(b + sign (sqrt(max(b*b - 4.d0*a*c,0.d0)), b) )
        if (temp.eq.0) then
          tau = 0.d0
        else
          tau = c / temp
        end if
c
c Make sure TAU falls in the interval -1 < TAU < 0
        tau = min(tau, 0.d0)
        tau = max(tau, -1.d0)
c
c Calculate TMIN and D2MIN
        tmin = tau * h
        temp = 1.d0 + tau
        d2min = tau*tau*((3.d0+2.d0*tau)*d0 + temp*h*d0t)
     %        + temp*temp*((1.d0-2.d0*tau)*d1 + tau*h*d1t)
c
c Make sure D2MIN is not negative
        d2min = max(d2min, 0.d0)
      end if
c
c------------------------------------------------------------------------------
c
      return
      end

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c    MCE_BOX.FOR    (ErikSoft   30 September 2000)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Given initial and final coordinates and velocities, the routine returns
c the X and Y coordinates of a box bounding the motion in between the
c end points.
c
c If the X or Y velocity changes sign, the routine performs a quadratic
c interpolation to estimate the corresponding extreme value of X or Y.
c
c------------------------------------------------------------------------------
c
      subroutine mce_box (nbod,h,x0,v0,x1,v1,xmin,xmax,ymin,ymax)
c
      implicit none
      include 'mercury.inc'
c
c Input/Output
      integer nbod
      real*8 h,x0(3,nbod), x1(3,nbod), v0(3,nbod),v1(3,nbod)
      real*8   xmin(nbod), xmax(nbod), ymin(nbod),ymax(nbod)
c
c Local
      integer j
      real*8 temp
c
c------------------------------------------------------------------------------
c
      do j = 2, nbod
        xmin(j) = min (x0(1,j), x1(1,j))
        xmax(j) = max (x0(1,j), x1(1,j))
        ymin(j) = min (x0(2,j), x1(2,j))
        ymax(j) = max (x0(2,j), x1(2,j))
c
c If velocity changes sign, do an interpolation
        if ((v0(1,j).lt.0.and.v1(1,j).gt.0).or.
     %      (v0(1,j).gt.0.and.v1(1,j).lt.0)) then
          temp = (v0(1,j)*x1(1,j) - v1(1,j)*x0(1,j)
     %            - .5d0*h*v0(1,j)*v1(1,j)) / (v0(1,j) - v1(1,j))
          xmin(j) = min (xmin(j),temp)
          xmax(j) = max (xmax(j),temp)
        end if
c
        if ((v0(2,j).lt.0.and.v1(2,j).gt.0).or.
     %      (v0(2,j).gt.0.and.v1(2,j).lt.0)) then
          temp = (v0(2,j)*x1(2,j) - v1(2,j)*x0(2,j)
     %            - .5d0*h*v0(2,j)*v1(2,j)) / (v0(2,j) - v1(2,j))
          ymin(j) = min (ymin(j),temp)
          ymax(j) = max (ymax(j),temp)
        end if
      end do
c
c------------------------------------------------------------------------------
c
      return
      end

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCE_MERG.FOR    (ErikSoft   2 October 2000)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% c
c Author: John E. Chambers
c
c Merges objects I and J inelastically to produce a single new body by
c conserving mass and linear momentum.
c   If J <= NBIG, then J is a Big body
c   If J >  NBIG, then J is a Small body
c   If I = 0, then I is the central body
c
c N.B. All coordinates and velocities must be with respect to central body.
c ===
c
c------------------------------------------------------------------------------
c
      subroutine mce_merg (jcen,i,j,nbod,nbig,m,xh,vh,s,stat,elost)
c
      implicit none
      include 'mercury.inc'
c
c Input/Output
      integer i, j, nbod, nbig, stat(nbod)
      real*8 jcen(3),m(nbod),xh(3,nbod),vh(3,nbod),s(3,nbod),elost
c
c Local
      integer k
      real*8 tmp1, tmp2, dx, dy, dz, du, dv, dw, msum, mredu, msum_1
      real*8 e0, e1, l2
c
c------------------------------------------------------------------------------
c
c If a body hits the central body
      if (i.le.1) then
        call mxx_en (jcen,nbod,nbig,m,xh,vh,s,e0,l2)
c
c If a body hit the central body...
        msum   = m(1) + m(j)
        msum_1 = 1.d0 / msum
        mredu  = m(1) * m(j) * msum_1
        dx = xh(1,j)
        dy = xh(2,j)
        dz = xh(3,j)
        du = vh(1,j)
        dv = vh(2,j)
        dw = vh(3,j)
c
c Calculate new spin angular momentum of the central body
        s(1,1) = s(1,1)  +  s(1,j)  +  mredu * (dy * dw  -  dz * dv)
        s(2,1) = s(2,1)  +  s(2,j)  +  mredu * (dz * du  -  dx * dw)
        s(3,1) = s(3,1)  +  s(3,j)  +  mredu * (dx * dv  -  dy * du)
c
c Calculate shift in barycentric coords and velocities of central body
        tmp2 = m(j) * msum_1
        xh(1,1) = tmp2 * xh(1,j)
        xh(2,1) = tmp2 * xh(2,j)
        xh(3,1) = tmp2 * xh(3,j)
        vh(1,1) = tmp2 * vh(1,j)
        vh(2,1) = tmp2 * vh(2,j)
        vh(3,1) = tmp2 * vh(3,j)
        m(1) = msum
        m(j) = 0.d0
        s(1,j) = 0.d0
        s(2,j) = 0.d0
        s(3,j) = 0.d0
c
c Shift the heliocentric coordinates and velocities of all bodies
        do k = 2, nbod
          xh(1,k) = xh(1,k) - xh(1,1)
          xh(2,k) = xh(2,k) - xh(2,1)
          xh(3,k) = xh(3,k) - xh(3,1)
          vh(1,k) = vh(1,k) - vh(1,1)
          vh(2,k) = vh(2,k) - vh(2,1)
          vh(3,k) = vh(3,k) - vh(3,1)
        end do
c
c Calculate energy loss due to the collision
        call mxx_en (jcen,nbod,nbig,m,xh,vh,s,e1,l2)
        elost = elost + (e0 - e1)
      else
c
c If two bodies collided...
        msum   = m(i) + m(j)
        msum_1 = 1.d0 / msum
        mredu  = m(i) * m(j) * msum_1
        dx = xh(1,i) - xh(1,j)
        dy = xh(2,i) - xh(2,j)
        dz = xh(3,i) - xh(3,j)
        du = vh(1,i) - vh(1,j)
        dv = vh(2,i) - vh(2,j)
        dw = vh(3,i) - vh(3,j)
c
c Calculate energy loss due to the collision
        elost = elost  +  .5d0 * mredu * (du*du + dv*dv + dw*dw)
     %        -  m(i) * m(j) / sqrt(dx*dx + dy*dy + dz*dz)
c
c Calculate spin angular momentum of the new body
        s(1,i) = s(1,i)  +  s(1,j)  +  mredu * (dy * dw  -  dz * dv)
        s(2,i) = s(2,i)  +  s(2,j)  +  mredu * (dz * du  -  dx * dw)
        s(3,i) = s(3,i)  +  s(3,j)  +  mredu * (dx * dv  -  dy * du)
c
c Calculate new coords and velocities by conserving centre of mass & momentum
        tmp1 = m(i) * msum_1
        tmp2 = m(j) * msum_1
        xh(1,i) = xh(1,i) * tmp1  +  xh(1,j) * tmp2
        xh(2,i) = xh(2,i) * tmp1  +  xh(2,j) * tmp2
        xh(3,i) = xh(3,i) * tmp1  +  xh(3,j) * tmp2
        vh(1,i) = vh(1,i) * tmp1  +  vh(1,j) * tmp2
        vh(2,i) = vh(2,i) * tmp1  +  vh(2,j) * tmp2
        vh(3,i) = vh(3,i) * tmp1  +  vh(3,j) * tmp2
        m(i) = msum
      end if
c
c Flag the lost body for removal, and move it away from the new body
      stat(j) = -2
      xh(1,j) = -xh(1,j)
      xh(2,j) = -xh(2,j)
      xh(3,j) = -xh(3,j)
      vh(1,j) = -vh(1,j)
      vh(2,j) = -vh(2,j)
      vh(3,j) = -vh(3,j)
      m(j)   = 0.d0
      s(1,j) = 0.d0
      s(2,j) = 0.d0
      s(3,j) = 0.d0
c
c------------------------------------------------------------------------------
c
      return
      end

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MDT_BS2.FOR    (ErikSoft   2 March 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Integrates NBOD bodies (of which NBIG are Big) for one timestep H0
c using the Bulirsch-Stoer method. The accelerations are calculated using the
c subroutine FORCE. The accuracy of the step is approximately determined
c by the tolerance parameter TOL.
c
c N.B. This version only works for conservative systems (i.e. force is a
c ===  function of position only) !!!! Hence, non-gravitational forces
c      and post-Newtonian corrections cannot be used.
c
c N.B. Input/output must be in coordinates with respect to the central body.
c ===
c
c------------------------------------------------------------------------------
c
      subroutine mdt_bs2 (time,h0,hdid,tol,jcen,nbod,nbig,mass,x0,v0,s,
     %  rphys,rcrit,ngf,stat,dtflag,ngflag,opt,nce,ice,jce)
c
      implicit none
      include 'mercury.inc'
c
      real*8 SHRINK,GROW
      parameter (SHRINK=.55d0,GROW=1.3d0)
c
c Input/Output
      integer nbod, nbig, opt(8), stat(nbod), dtflag, ngflag
      real*8 time,h0,hdid,tol,jcen(3),mass(nbod),x0(3,nbod),v0(3,nbod)
      real*8 s(3,nbod),ngf(4,nbod),rphys(nbod),rcrit(nbod)
      integer nce,ice(nce),jce(nce)
c
c Local
      integer j,j1,k,n
      real*8 tmp0,tmp1,tmp2,errmax,tol2,h,h2(12),hby2,h2by2
      real*8 xend(3,NMAX),b(3,NMAX),c(3,NMAX)
      real*8 a(3,NMAX),a0(3,NMAX),d(6,NMAX,12),xscal(NMAX),vscal(NMAX)
c
c------------------------------------------------------------------------------
c
      tol2 = tol * tol
c
c Calculate arrays used to scale the relative error (R^2 for position and
c V^2 for velocity).
      do k = 2, nbod
        tmp1 = x0(1,k)*x0(1,k) + x0(2,k)*x0(2,k) + x0(3,k)*x0(3,k)
        tmp2 = v0(1,k)*v0(1,k) + v0(2,k)*v0(2,k) + v0(3,k)*v0(3,k)
        xscal(k) = 1.d0 / tmp1
        vscal(k) = 1.d0 / tmp2
      end do
c
c Calculate accelerations at the start of the step
      call mfo_hkce (time,jcen,nbod,nbig,mass,x0,v0,s,rcrit,a0,stat,ngf,
     %  ngflag,opt,nce,ice,jce)
c
 100  continue
c
c For each value of N, do a modified-midpoint integration with N substeps
      do n = 1, 12
        h = h0 / (dble(n))
        hby2  = .5d0 * h
        h2(n) = h * h
        h2by2 = .5d0 * h2(n)
c
        do k = 2, nbod
          b(1,k) = .5d0*a0(1,k)
          b(2,k) = .5d0*a0(2,k)
          b(3,k) = .5d0*a0(3,k)
          c(1,k) = 0.d0
          c(2,k) = 0.d0
          c(3,k) = 0.d0
          xend(1,k) = h2by2 * a0(1,k)  +  h * v0(1,k)  +  x0(1,k)
          xend(2,k) = h2by2 * a0(2,k)  +  h * v0(2,k)  +  x0(2,k)
          xend(3,k) = h2by2 * a0(3,k)  +  h * v0(3,k)  +  x0(3,k)
        end do
c
        do j = 2, n
          call mfo_hkce (time,jcen,nbod,nbig,mass,xend,v0,s,rcrit,a,
     %      stat,ngf,ngflag,opt,nce,ice,jce)
          tmp0 = h * dble(j)
          do k = 2, nbod
            b(1,k) = b(1,k) + a(1,k)
            b(2,k) = b(2,k) + a(2,k)
            b(3,k) = b(3,k) + a(3,k)
            c(1,k) = c(1,k) + b(1,k)
            c(2,k) = c(2,k) + b(2,k)
            c(3,k) = c(3,k) + b(3,k)
            xend(1,k) = h2(n)*c(1,k) + h2by2*a0(1,k) + tmp0*v0(1,k)
     %                + x0(1,k)
            xend(2,k) = h2(n)*c(2,k) + h2by2*a0(2,k) + tmp0*v0(2,k)
     %                + x0(2,k)
            xend(3,k) = h2(n)*c(3,k) + h2by2*a0(3,k) + tmp0*v0(3,k)
     %                + x0(3,k)
          end do
        end do
c
        call mfo_hkce (time,jcen,nbod,nbig,mass,xend,v0,s,rcrit,a,stat,
     %    ngf,ngflag,opt,nce,ice,jce)
c
        do k = 2, nbod
          d(1,k,n) = xend(1,k)
          d(2,k,n) = xend(2,k)
          d(3,k,n) = xend(3,k)
          d(4,k,n) = h*b(1,k) + hby2*a(1,k) + v0(1,k)
          d(5,k,n) = h*b(2,k) + hby2*a(2,k) + v0(2,k)
          d(6,k,n) = h*b(3,k) + hby2*a(3,k) + v0(3,k)
        end do
c
c Update the D array, used for polynomial extrapolation
        do j = n - 1, 1, -1
          j1 = j + 1
          tmp0 = 1.d0 / (h2(j) - h2(n))
          tmp1 = tmp0 * h2(j1)
          tmp2 = tmp0 * h2(n)
          do k = 2, nbod
            d(1,k,j) = tmp1 * d(1,k,j1)  -  tmp2 * d(1,k,j)
            d(2,k,j) = tmp1 * d(2,k,j1)  -  tmp2 * d(2,k,j)
            d(3,k,j) = tmp1 * d(3,k,j1)  -  tmp2 * d(3,k,j)
            d(4,k,j) = tmp1 * d(4,k,j1)  -  tmp2 * d(4,k,j)
            d(5,k,j) = tmp1 * d(5,k,j1)  -  tmp2 * d(5,k,j)
            d(6,k,j) = tmp1 * d(6,k,j1)  -  tmp2 * d(6,k,j)
          end do
        end do
c
c After several integrations, test the relative error on extrapolated values
        if (n.gt.3) then
          errmax = 0.d0
c
c Maximum relative position and velocity errors (last D term added)
          do k = 2, nbod
            tmp1 = max( d(1,k,1)*d(1,k,1), d(2,k,1)*d(2,k,1),
     %                  d(3,k,1)*d(3,k,1) )
            tmp2 = max( d(4,k,1)*d(4,k,1), d(5,k,1)*d(2,k,1),
     %                  d(6,k,1)*d(6,k,1) )
            errmax = max( errmax, tmp1*xscal(k), tmp2*vscal(k) )
          end do
c
c If error is smaller than TOL, update position and velocity arrays and exit
          if (errmax.le.tol2) then
            do k = 2, nbod
              x0(1,k) = d(1,k,1)
              x0(2,k) = d(2,k,1)
              x0(3,k) = d(3,k,1)
              v0(1,k) = d(4,k,1)
              v0(2,k) = d(5,k,1)
              v0(3,k) = d(6,k,1)
            end do
c
            do j = 2, n
              do k = 2, nbod
                x0(1,k) = x0(1,k) + d(1,k,j)
                x0(2,k) = x0(2,k) + d(2,k,j)
                x0(3,k) = x0(3,k) + d(3,k,j)
                v0(1,k) = v0(1,k) + d(4,k,j)
                v0(2,k) = v0(2,k) + d(5,k,j)
                v0(3,k) = v0(3,k) + d(6,k,j)
              end do
            end do
c
c Save the actual stepsize used
            hdid = h0
c
c Recommend a new stepsize for the next call to this subroutine
            if (n.ge.8) h0 = h0 * SHRINK
            if (n.lt.7) h0 = h0 * GROW
            return
          end if
        end if
c
      end do
c
c If errors were too large, redo the step with half the previous step size.
      h0 = h0 * .5d0
      goto 100
c
c------------------------------------------------------------------------------
c
      end

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MFO_USER.FOR    (ErikSoft   2 March 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Applies an arbitrary force, defined by the user.
c
c If using with the symplectic algorithm MAL_MVS, the force should be
c small compared with the force from the central object.
c If using with the conservative Bulirsch-Stoer algorithm MAL_BS2, the
c force should not be a function of the velocities.
c
c N.B. All coordinates and velocities must be with respect to central body
c ===
c------------------------------------------------------------------------------
c
      subroutine mfo_user (time,jcen,nbod,nbig,m,x,v,a)
c
      implicit none
      include 'mercury.inc'
c
c Input/Output
      integer nbod, nbig
      real*8 time,jcen(3),m(nbod),x(3,nbod),v(3,nbod),a(3,nbod)
c
c Local
      integer j
c
c------------------------------------------------------------------------------
c
      do j = 1, nbod
        a(1,j) = 0.d0
        a(2,j) = 0.d0
        a(3,j) = 0.d0
      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     MXX_EN.FOR    (ErikSoft   21 February 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Calculates the total energy and angular-momentum for a system of objects
c with masses M, coordinates X, velocities V and spin angular momenta S.
c
c N.B. All coordinates and velocities must be with respect to the central
c ===  body.
c
c------------------------------------------------------------------------------
c
      subroutine mxx_en  (jcen,nbod,nbig,m,xh,vh,s,e,l2)
c
      implicit none
      include 'mercury.inc'
c
c Input/Output
      integer nbod,nbig
      real*8 jcen(3),m(nbod),xh(3,nbod),vh(3,nbod),s(3,nbod),e,l2
c
c Local
      integer j,k,iflag,itmp(8)
      real*8 x(3,NMAX),v(3,NMAX),temp,dx,dy,dz,r2,tmp,ke,pe,l(3)
      real*8 r_1,r_2,r_4,r_6,u2,u4,u6,tmp2(4,NMAX)
c
c------------------------------------------------------------------------------
c
      ke = 0.d0
      pe = 0.d0
      l(1) = 0.d0
      l(2) = 0.d0
      l(3) = 0.d0
c
c Convert to barycentric coordinates and velocities
      call mco_h2b(temp,jcen,nbod,nbig,temp,m,xh,vh,x,v,tmp2,iflag,itmp)
c
c Do the spin angular momenta first (probably the smallest terms)
      do j = 1, nbod
        l(1) = l(1) + s(1,j)
        l(2) = l(2) + s(2,j)
        l(3) = l(3) + s(3,j)
      end do
c
c Orbital angular momentum and kinetic energy terms
      do j = 1, nbod
        l(1) = l(1)  +  m(j)*(x(2,j) * v(3,j)  -  x(3,j) * v(2,j))
        l(2) = l(2)  +  m(j)*(x(3,j) * v(1,j)  -  x(1,j) * v(3,j))
        l(3) = l(3)  +  m(j)*(x(1,j) * v(2,j)  -  x(2,j) * v(1,j))
        ke = ke + m(j)*(v(1,j)*v(1,j)+v(2,j)*v(2,j)+v(3,j)*v(3,j))
      end do
c
c Potential energy terms due to pairs of bodies
      do j = 2, nbod
        tmp = 0.d0
        do k = j + 1, nbod
          dx = x(1,k) - x(1,j)
          dy = x(2,k) - x(2,j)
          dz = x(3,k) - x(3,j)
          r2 = dx*dx + dy*dy + dz*dz
          if (r2.ne.0) tmp = tmp + m(k) / sqrt(r2)
        end do
        pe = pe  -  tmp * m(j)
      end do
c
c Potential energy terms involving the central body
      do j = 2, nbod
        dx = x(1,j) - x(1,1)
        dy = x(2,j) - x(2,1)
        dz = x(3,j) - x(3,1)
        r2 = dx*dx + dy*dy + dz*dz
        if (r2.ne.0) pe = pe  -  m(1) * m(j) / sqrt(r2)
      end do
c
c Corrections for oblateness
      if (jcen(1).ne.0.or.jcen(2).ne.0.or.jcen(3).ne.0) then
        do j = 2, nbod
          r2 = xh(1,j)*xh(1,j) + xh(2,j)*xh(2,j) + xh(3,j)*xh(3,j)
          r_1 = 1.d0 / sqrt(r2)
          r_2 = r_1 * r_1
          r_4 = r_2 * r_2
          r_6 = r_4 * r_2
          u2 = xh(3,j) * xh(3,j) * r_2
          u4 = u2 * u2
          u6 = u4 * u2
          pe = pe + m(1) * m(j) * r_1
     %       * (jcen(1) * r_2 * (1.5d0*u2 - 0.5d0)
     %       +  jcen(2) * r_4 * (4.375d0*u4 - 3.75d0*u2 + .375d0)
     %       +  jcen(3) * r_6
     %       *(14.4375d0*u6 - 19.6875d0*u4 + 6.5625d0*u2 - .3125d0))
        end do
      end if
c
      e = .5d0 * ke  +  pe
      l2 = sqrt(l(1)*l(1) + l(2)*l(2) + l(3)*l(3))
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     MFO_HY.FOR    (ErikSoft   2 October 2000)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Calculates accelerations due to the Interaction part of the Hamiltonian
c of a hybrid symplectic integrator for a set of NBOD bodies (NBIG of which
c are Big), where Small bodies do not interact with one another.
c
c------------------------------------------------------------------------------
c
      subroutine mfo_hy (jcen,nbod,nbig,m,x,rcrit,a,stat)
c
      implicit none
      include 'mercury.inc'
c
c Input/Output
      integer nbod, nbig, stat(nbod)
      real*8 jcen(3), m(nbod), x(3,nbod), a(3,nbod), rcrit(nbod)
c
c Local
      integer k
      real*8 aobl(3,NMAX),acen(3)
c
c------------------------------------------------------------------------------
c
c Initialize accelerations to zero
      do k = 1, nbod
        a(1,k) = 0.d0
        a(2,k) = 0.d0
        a(3,k) = 0.d0
      end do
c
c Calculate direct terms
      call mfo_drct (2,nbod,nbig,m,x,rcrit,a,stat)
c
c Add accelerations due to oblateness of the central body
      if (jcen(1).ne.0.or.jcen(2).ne.0.or.jcen(3).ne.0) then
        call mfo_obl (jcen,nbod,m,x,aobl,acen)
        do k = 2, nbod
          a(1,k) = a(1,k) + aobl(1,k) - acen(1)
          a(2,k) = a(2,k) + aobl(2,k) - acen(2)
          a(3,k) = a(3,k) + aobl(3,k) - acen(3)
        end do
      end if
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     MFO_HKCE.FOR    (ErikSoft   27 February 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Calculates accelerations due to the Keplerian part of the Hamiltonian
c of a hybrid symplectic integrator, when close encounters are taking place,
c for a set of NBOD bodies (NBIG of which are Big). Note that Small bodies
c do not interact with one another.
c
c------------------------------------------------------------------------------
c
      subroutine mfo_hkce (time,jcen,nbod,nbig,m,x,v,spin,rcrit,a,stat,
     %  ngf,ngflag,opt,nce,ice,jce)
c
      implicit none
      include 'mercury.inc'
c
c Input/Output
      integer nbod,nbig,stat(nbod),ngflag,opt(8),nce,ice(nce),jce(nce)
      real*8 time,jcen(3),rcrit(nbod),ngf(4,nbod),m(nbod)
      real*8 x(3,nbod),v(3,nbod),a(3,nbod),spin(3,nbod)
c
c Local
      integer i, j, k
      real*8 tmp2,dx,dy,dz,s,s_1,s2,s_3,faci,facj,rc,rc2,q,q2,q3,q4,q5
c
c------------------------------------------------------------------------------
c
c Initialize accelerations
      do j = 1, nbod
        a(1,j) = 0.d0
        a(2,j) = 0.d0
        a(3,j) = 0.d0
      end do
c
c Direct terms
      do k = 1, nce
        i = ice(k)
        j = jce(k)
        dx = x(1,j) - x(1,i)
        dy = x(2,j) - x(2,i)
        dz = x(3,j) - x(3,i)
        s2 = dx * dx  +  dy * dy  +  dz * dz
        rc = max (rcrit(i), rcrit(j))
        rc2 = rc * rc
c
        if (s2.lt.rc2) then
          s_1 = 1.d0 / sqrt(s2)
          s_3 = s_1 * s_1 * s_1
          if (s2.le.0.01*rc2) then
            tmp2 = s_3
          else
            s = 1.d0 / s_1
            q = (s - 0.1d0*rc) / (0.9d0 * rc)
            q2 = q * q
            q3 = q * q2
            q4 = q2 * q2
            q5 = q2 * q3
            tmp2 = (1.d0 - 10.d0*q3 + 15.d0*q4 - 6.d0*q5) * s_3
          end if
c
          faci = tmp2 * m(i)
          facj = tmp2 * m(j)
          a(1,j) = a(1,j)  -  faci * dx
          a(2,j) = a(2,j)  -  faci * dy
          a(3,j) = a(3,j)  -  faci * dz
          a(1,i) = a(1,i)  +  facj * dx
          a(2,i) = a(2,i)  +  facj * dy
          a(3,i) = a(3,i)  +  facj * dz
        end if
      end do
c
c Solar terms
      do i = 2, nbod
        s2 = x(1,i)*x(1,i) + x(2,i)*x(2,i) + x(3,i)*x(3,i)
        s_1 = 1.d0 / sqrt(s2)
        tmp2 = m(1) * s_1 * s_1 * s_1
        a(1,i) = a(1,i)  -  tmp2 * x(1,i)
        a(2,i) = a(2,i)  -  tmp2 * x(2,i)
        a(3,i) = a(3,i)  -  tmp2 * x(3,i)
      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c    MCE_SNIF.FOR    (ErikSoft   3 October 2000)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Given initial and final coordinates and velocities X and V, and a timestep
c H, the routine estimates which objects were involved in a close encounter
c during the timestep. The routine examines all objects with indices I >= I0.
c
c Returns an array CE, which for each object is:
c                           0 if it will undergo no encounters
c                           2 if it will pass within RCRIT of a Big body
c
c Also returns arrays ICE and JCE, containing the indices of each pair of
c objects estimated to have undergone an encounter.
c
c N.B. All coordinates must be with respect to the central body!!!!
c ===
c
c------------------------------------------------------------------------------
c
      subroutine mce_snif (h,i0,nbod,nbig,x0,v0,x1,v1,rcrit,ce,nce,ice,
     %  jce)
c
      implicit none
      include 'mercury.inc'
c
c Input/Output
      integer i0,nbod,nbig,ce(nbod),nce,ice(NMAX),jce(NMAX)
      real*8 x0(3,nbod),v0(3,nbod),x1(3,nbod),v1(3,nbod),h,rcrit(nbod)
c
c Local
      integer i,j
      real*8 d0,d1,d0t,d1t,d2min,temp,tmin,rc,rc2
      real*8 dx0,dy0,dz0,du0,dv0,dw0,dx1,dy1,dz1,du1,dv1,dw1
      real*8 xmin(NMAX),xmax(NMAX),ymin(NMAX),ymax(NMAX)
c
c------------------------------------------------------------------------------
c
      if (i0.le.0) i0 = 2
      nce = 0
      do j = 2, nbod
        ce(j) = 0
      end do
c
c Calculate maximum and minimum values of x and y coordinates
      call mce_box (nbod,h,x0,v0,x1,v1,xmin,xmax,ymin,ymax)
c
c Adjust values for the Big bodies by symplectic close-encounter distance
      do j = i0, nbig
        xmin(j) = xmin(j) - rcrit(j)
        xmax(j) = xmax(j) + rcrit(j)
        ymin(j) = ymin(j) - rcrit(j)
        ymax(j) = ymax(j) + rcrit(j)
      end do
c
c Identify pairs whose X-Y boxes overlap, and calculate minimum separation
      do i = i0, nbig
        do j = i + 1, nbod
          if (xmax(i).ge.xmin(j).and.xmax(j).ge.xmin(i)
     %      .and.ymax(i).ge.ymin(j).and.ymax(j).ge.ymin(i)) then
c
c Determine the maximum separation that would qualify as an encounter
            rc = max(rcrit(i), rcrit(j))
            rc2 = rc * rc
c
c Calculate initial and final separations
            dx0 = x0(1,i) - x0(1,j)
            dy0 = x0(2,i) - x0(2,j)
            dz0 = x0(3,i) - x0(3,j)
            dx1 = x1(1,i) - x1(1,j)
            dy1 = x1(2,i) - x1(2,j)
            dz1 = x1(3,i) - x1(3,j)
            d0 = dx0*dx0 + dy0*dy0 + dz0*dz0
            d1 = dx1*dx1 + dy1*dy1 + dz1*dz1
c
c Check for a possible minimum in between
            du0 = v0(1,i) - v0(1,j)
            dv0 = v0(2,i) - v0(2,j)
            dw0 = v0(3,i) - v0(3,j)
            du1 = v1(1,i) - v1(1,j)
            dv1 = v1(2,i) - v1(2,j)
            dw1 = v1(3,i) - v1(3,j)
            d0t = (dx0*du0 + dy0*dv0 + dz0*dw0) * 2.d0
            d1t = (dx1*du1 + dy1*dv1 + dz1*dw1) * 2.d0
c
c If separation derivative changes sign, find the minimum separation
            d2min = HUGE
            if (d0t*h.le.0.and.d1t*h.ge.0) call mce_min (d0,d1,d0t,d1t,
     %        h,d2min,tmin)
c
c If minimum separation is small enough, flag this as a possible encounter
            temp = min (d0,d1,d2min)
            if (temp.le.rc2) then
              ce(i) = 2
              ce(j) = 2
              nce = nce + 1
              ice(nce) = i
              jce(nce) = j
            end if
          end if
        end do
      end do
c
c------------------------------------------------------------------------------
c
      return
      end

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_SINE.FOR    (ErikSoft  17 April 1997)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Calculates sin and cos of an angle X (in radians).
c
c------------------------------------------------------------------------------
c
      subroutine mco_sine (x,sx,cx)
c
      implicit none
c
c Input/Output
      real*8 x,sx,cx
c
c Local
      real*8 pi,twopi
c
c------------------------------------------------------------------------------
c
      pi = 3.141592653589793d0
      twopi = 2.d0 * pi
c
      if (x.gt.0) then
        x = mod(x,twopi)
      else
        x = mod(x,twopi) + twopi
      end if
c
      cx = cos(x)
c
      if (x.gt.pi) then
        sx = -sqrt(1.d0 - cx*cx)
      else
        sx =  sqrt(1.d0 - cx*cx)
      end if
c
c------------------------------------------------------------------------------
c
      return
      end

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     MFO_OBL.FOR    (ErikSoft   2 October 2000)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Calculates barycentric accelerations of NBOD bodies due to oblateness of
c the central body. Also returns the corresponding barycentric acceleration
c of the central body.
c
c N.B. All coordinates must be with respect to the central body!!!!
c ===
c------------------------------------------------------------------------------
c
      subroutine mfo_obl (jcen,nbod,m,x,a,acen)
c
      implicit none
      include 'mercury.inc'
c
c Input/Output
      integer nbod
      real*8 jcen(3), m(nbod), x(3,nbod), a(3,nbod), acen(3)
c
c Local
      integer i
      real*8 jr2,jr4,jr6,r2,r_1,r_2,r_3,u2,u4,u6,tmp1,tmp2,tmp3,tmp4
c
c------------------------------------------------------------------------------
c
      acen(1) = 0.d0
      acen(2) = 0.d0
      acen(3) = 0.d0
c
      do i = 2, nbod
c
c Calculate barycentric accelerations on the objects
        r2 = x(1,i)*x(1,i) + x(2,i)*x(2,i) + x(3,i)*x(3,i)
        r_1 = 1.d0 / sqrt(r2)
        r_2 = r_1 * r_1
        r_3 = r_2 * r_1
        jr2 = jcen(1) * r_2
        jr4 = jcen(2) * r_2 * r_2
        jr6 = jcen(3) * r_2 * r_2 * r_2
        u2 = x(3,i) * x(3,i) * r_2
        u4 = u2 * u2
        u6 = u4 * u2
c
        tmp1 = m(1) * r_3
        tmp2 =jr2*(7.5d0*u2 - 1.5d0)
     %       +jr4*(39.375d0*u4 - 26.25d0*u2 + 1.875d0)
     %       +jr6*(187.6875d0*u6 -216.5625d0*u4 +59.0625d0*u2 -2.1875d0)
        tmp3 = jr2*3.d0 + jr4*(17.5d0*u2 - 7.5d0)
     %       + jr6*(86.625d0*u4 - 78.75d0*u2 + 13.125d0)
c
        a(1,i) = x(1,i) * tmp1 * tmp2
        a(2,i) = x(2,i) * tmp1 * tmp2
        a(3,i) = x(3,i) * tmp1 * (tmp2 - tmp3)
c
c Calculate barycentric accelerations on the central body
        tmp4 = m(i) / m(1)
        acen(1) = acen(1)  -  tmp4 * a(1,i)
        acen(2) = acen(2)  -  tmp4 * a(2,i)
        acen(3) = acen(3)  -  tmp4 * a(3,i)
      end do
c
c------------------------------------------------------------------------------
c
      return
      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     MFO_DRCT.FOR    (ErikSoft   27 February 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Calculates direct accelerations between bodies in the interaction part
c of the Hamiltonian of a symplectic integrator that partitions close
c encounter terms (e.g. hybrid symplectic algorithms or SyMBA).
c The routine calculates accelerations between all pairs of bodies with
c indices I >= I0.
c
c------------------------------------------------------------------------------
c
      subroutine mfo_drct (i0,nbod,nbig,m,x,rcrit,a,stat)
c
      implicit none
      include 'mercury.inc'
c
c Input/Output
      integer i0, nbod, nbig, stat(nbod)
      real*8 m(nbod), x(3,nbod), a(3,nbod), rcrit(nbod)
c
c Local
      integer i,j
      real*8 dx,dy,dz,s,s_1,s2,s_3,rc,rc2,q,q2,q3,q4,q5,tmp2,faci,facj
c
c------------------------------------------------------------------------------
c
      if (i0.le.0) i0 = 2
c
      do i = i0, nbig
        do j = i + 1, nbod
          dx = x(1,j) - x(1,i)
          dy = x(2,j) - x(2,i)
          dz = x(3,j) - x(3,i)
          s2 = dx * dx  +  dy * dy  +  dz * dz
          rc = max(rcrit(i), rcrit(j))
          rc2 = rc * rc
c
          if (s2.ge.rc2) then
            s_1 = 1.d0 / sqrt(s2)
            tmp2 = s_1 * s_1 * s_1
          else if (s2.le.0.01*rc2) then
            tmp2 = 0.d0
          else
            s_1 = 1.d0 / sqrt(s2)
            s   = 1.d0 / s_1
            s_3 = s_1 * s_1 * s_1
            q = (s - 0.1d0*rc) / (0.9d0 * rc)
            q2 = q  * q
            q3 = q  * q2
            q4 = q2 * q2
            q5 = q2 * q3
            tmp2 = (10.d0*q3 - 15.d0*q4 + 6.d0*q5) * s_3
          end if
c
          faci = tmp2 * m(i)
          facj = tmp2 * m(j)
          a(1,j) = a(1,j)  -  faci * dx
          a(2,j) = a(2,j)  -  faci * dy
          a(3,j) = a(3,j)  -  faci * dz
          a(1,i) = a(1,i)  +  facj * dx
          a(2,i) = a(2,i)  +  facj * dy
          a(3,i) = a(3,i)  +  facj * dz
        end do
      end do
c
c------------------------------------------------------------------------------
c
      return
      end

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_H2B.FOR    (ErikSoft   2 March 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Converts coordinates with respect to the central body to barycentric
c coordinates.
c
c------------------------------------------------------------------------------
c
      subroutine mco_h2b (time,jcen,nbod,nbig,h,m,xh,vh,x,v,ngf,ngflag,
     %  opt)
c
      implicit none
c
c Input/Output
      integer nbod,nbig,ngflag,opt(8)
      real*8 time,jcen(3),h,m(nbod),xh(3,nbod),vh(3,nbod),x(3,nbod)
      real*8 v(3,nbod),ngf(4,nbod)
c
c Local
      integer j
      real*8 mtot,temp
c
c------------------------------------------------------------------------------
c
      mtot = 0.d0
      x(1,1) = 0.d0
      x(2,1) = 0.d0
      x(3,1) = 0.d0
      v(1,1) = 0.d0
      v(2,1) = 0.d0
      v(3,1) = 0.d0
c
c Calculate coordinates and velocities of the central body
      do j = 2, nbod
        mtot = mtot  +  m(j)
        x(1,1) = x(1,1)  +  m(j) * xh(1,j)
        x(2,1) = x(2,1)  +  m(j) * xh(2,j)
        x(3,1) = x(3,1)  +  m(j) * xh(3,j)
        v(1,1) = v(1,1)  +  m(j) * vh(1,j)
        v(2,1) = v(2,1)  +  m(j) * vh(2,j)
        v(3,1) = v(3,1)  +  m(j) * vh(3,j)
      enddo
c
      temp = -1.d0 / (mtot + m(1))
      x(1,1) = temp * x(1,1)
      x(2,1) = temp * x(2,1)
      x(3,1) = temp * x(3,1)
      v(1,1) = temp * v(1,1)
      v(2,1) = temp * v(2,1)
      v(3,1) = temp * v(3,1)
c
c Calculate the barycentric coordinates and velocities
      do j = 2, nbod
        x(1,j) = xh(1,j) + x(1,1)
        x(2,j) = xh(2,j) + x(2,1)
        x(3,j) = xh(3,j) + x(3,1)
        v(1,j) = vh(1,j) + v(1,1)
        v(2,j) = vh(2,j) + v(2,1)
        v(3,j) = vh(3,j) + v(3,1)
      enddo
c
c------------------------------------------------------------------------------
c
      return
      end

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_DH2H.FOR    (ErikSoft   2 March 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Converts democratic heliocentric coordinates to coordinates with respect
c to the central body.
c
c------------------------------------------------------------------------------
c
      subroutine mco_dh2h (time,jcen,nbod,nbig,h,m,x,v,xh,vh)
c
      implicit none
c
c Input/Output
      integer nbod,nbig
      real*8 time,jcen(3),h,m(nbod),x(3,nbod),v(3,nbod),xh(3,nbod)
      real*8 vh(3,nbod)
c
c Local
      integer j
      real*8 mvsum(3),temp
c
c------------------------------------------------------------------------------
c
      mvsum(1) = 0.d0
      mvsum(2) = 0.d0
      mvsum(3) = 0.d0
c
      do j = 2, nbod
        xh(1,j) = x(1,j)
        xh(2,j) = x(2,j)
        xh(3,j) = x(3,j)
        mvsum(1) = mvsum(1)  +  m(j) * v(1,j)
        mvsum(2) = mvsum(2)  +  m(j) * v(2,j)
        mvsum(3) = mvsum(3)  +  m(j) * v(3,j)
      end do
c
      temp = 1.d0 / m(1)
      mvsum(1) = temp * mvsum(1)
      mvsum(2) = temp * mvsum(2)
      mvsum(3) = temp * mvsum(3)
c
      do j = 2, nbod
        vh(1,j) = v(1,j) + mvsum(1)
        vh(2,j) = v(2,j) + mvsum(2)
        vh(3,j) = v(3,j) + mvsum(3)
      end do
c
c------------------------------------------------------------------------------
c
      return
      end

c------------------------------------------------------------------
c
c*************************************************************************
c                        DRIFT_ONE.F
c*************************************************************************
c This subroutine does the danby-type drift for one particle, using
c appropriate vbles and redoing a drift if the accuracy is too poor
c (as flagged by the integer iflg).
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 x,y,z         ==>  initial position in jacobi coord
c                                    (real scalar)
c                 vx,vy,vz      ==>  initial position in jacobi coord
c                                    (real scalar)
c                 dt            ==>  time step
c             Output:
c                 x,y,z         ==>  final position in jacobi coord
c                                       (real scalars)
c                 vx,vy,vz      ==>  final position in jacobi coord
c                                       (real scalars)
c                 iflg          ==>  integer (zero for successful step)
c
c Authors:  Hal Levison & Martin Duncan
c Date:    2/10/93
c Last revision: 2/10/93
c

      subroutine drift_one(mu,x,y,z,vx,vy,vz,dt,iflg)

      include 'swift.inc'

c...  Inputs Only:
      real*8 mu,dt

c...  Inputs and Outputs:
      real*8 x,y,z
      real*8 vx,vy,vz

c...  Output
      integer iflg

c...  Internals:
      integer i
      real*8 dttmp

c----
c...  Executable code

           call drift_dan(mu,x,y,z,vx,vy,vz,dt,iflg)

      if(iflg .ne. 0) then

        do i = 1,10
          dttmp = dt/10.d0
               call drift_dan(mu,x,y,z,vx,vy,vz,dttmp,iflg)
          if(iflg .ne. 0) return
        enddo

      endif

        return
        end    
c
c*************************************************************************
c                        DRIFT_DAN.F
c*************************************************************************
c This subroutine does the Danby and decides which vbles to use
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 x0,y0,z0         ==>  initial position in jacobi coord
c                                    (real scalar)
c                 vx0,vy0,vz0      ==>  initial position in jacobi coord
c                                    (real scalar)
c                 dt0            ==>  time step
c             Output:
c                 x0,y0,z0         ==>  final position in jacobi coord
c                                       (real scalars)
c                 vx0,vy0,vz0      ==>  final position in jacobi coord
c                                       (real scalars)
c                 iflg             ==>  integer flag (zero if satisfactory)
c                    (non-zero if nonconvergence)
c
c Authors:  Hal Levison & Martin Duncan
c Date:    2/10/93
c Last revision: April 6/93 - MD adds dt and keeps dt0 unchanged

      subroutine drift_dan(mu,x0,y0,z0,vx0,vy0,vz0,dt0,iflg)

      include 'swift.inc'

c...  Inputs Only:
      real*8 mu,dt0

c...  Inputs and Outputs:
      real*8 x0,y0,z0
      real*8 vx0,vy0,vz0

c...  Output
      integer iflg

c...  Internals:
      real*8 x,y,z,vx,vy,vz,dt
      real*8 f,g,fdot,c1,c2
      real*8 c3,gdot
      real*8 u,alpha,fp,r0,v0s
      real*8 a,asq,en
      real*8 dm,ec,es,esq,xkep
      real*8 fchk,s,c

c----
c...  Executable code

c...  Set dt = dt0 to be sure timestep is not altered while solving
c...  for new coords.
      dt = dt0
      iflg = 0
        r0 = sqrt(x0*x0 + y0*y0 + z0*z0)
        v0s = vx0*vx0 + vy0*vy0 + vz0*vz0
        u = x0*vx0 + y0*vy0 + z0*vz0
        alpha = 2.0*mu/r0 - v0s

        if (alpha.gt.0.d0) then
           a = mu/alpha
           asq = a*a
           en = sqrt(mu/(a*asq))
           ec = 1.0d0 - r0/a
           es = u/(en*asq)
           esq = ec*ec + es*es
           dm = dt*en - int(dt*en/TWOPI)*TWOPI
           dt = dm/en
           if((dm*dm .gt. 0.16d0) .or. (esq.gt.0.36d0)) goto 100

           if(esq*dm*dm .lt. 0.0016) then

           call drift_kepmd(dm,es,ec,xkep,s,c)
           fchk = (xkep - ec*s +es*(1.-c) - dm)

          if(fchk*fchk .gt. DANBYB) then
        iflg = 1
        return
          endif

               fp = 1. - ec*c + es*s
               f = (a/r0) * (c-1.) + 1.
               g = dt + (s-xkep)/en
               fdot = - (a/(r0*fp))*en*s
               gdot = (c-1.)/fp + 1.

               x = x0*f + vx0*g
               y = y0*f + vy0*g
               z = z0*f + vz0*g
               vx = x0*fdot + vx0*gdot
               vy = y0*fdot + vy0*gdot
               vz = z0*fdot + vz0*gdot

               x0 = x
               y0 = y
               z0 = z
               vx0 = vx
               vy0 = vy
               vz0 = vz

          iflg = 0
          return

      endif

         endif

 100     call drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)

         if(iflg .eq.0) then
           f = 1.0 - (mu/r0)*c2
           g = dt - mu*c3
           fdot = -(mu/(fp*r0))*c1
           gdot = 1. - (mu/fp)*c2

           x = x0*f + vx0*g
           y = y0*f + vy0*g
           z = z0*f + vz0*g
           vx = x0*fdot + vx0*gdot
           vy = y0*fdot + vy0*gdot
           vz = z0*fdot + vz0*gdot

           x0 = x
           y0 = y
           z0 = z
           vx0 = vx
           vy0 = vy
           vz0 = vz
      endif

        return
        end   

c-----------------------------------------------------------------------------
c
c*************************************************************************
c                        DRIFT_KEPU.F
c*************************************************************************
c subroutine for solving kepler's equation using universal variables.
c
c             Input:
c                 dt            ==>  time step (real scalor)
c                 r0            ==>  Distance between `Sun' and paritcle
c                                     (real scalor)
c                 mu            ==>  Reduced mass of system (real scalor)
c                 alpha         ==>  energy (real scalor)
c                 u             ==>  angular momentun  (real scalor)
c             Output:
c                 fp            ==>  f' from p170
c                                       (real scalors)
c                 c1,c2,c3      ==>  c's from p171-172
c                                       (real scalors)
c                 iflg          ==>  =0 if converged; !=0 if not
c
c Author:  Hal Levison
c Date:    2/3/93
c Last revision: 2/3/93

      subroutine drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)

      include 'swift.inc'

c...  Inputs:
      real*8 dt,r0,mu,alpha,u

c...  Outputs:
      real*8 fp,c1,c2,c3
      integer iflg

c...  Internals:
      real*8 s,st,fo,fn

c----
c...  Executable code

        call drift_kepu_guess(dt,r0,mu,alpha,u,s)

        st = s
c..     store initial guess for possible use later in
c..     laguerre's method, in case newton's method fails.

        call drift_kepu_new(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
        if(iflg.ne.0) then
           call drift_kepu_fchk(dt,r0,mu,alpha,u,st,fo)
           call drift_kepu_fchk(dt,r0,mu,alpha,u,s,fn)
           if(abs(fo).lt.abs(fn)) then
               s = st
           endif
           call drift_kepu_lag(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
        endif

        return
        end    

c
c********************************************************************#
c                  DRIFT_KEPMD
c********************************************************************#
c  Subroutine for solving kepler's equation in difference form for an
c  ellipse, given SMALL dm and SMALL eccentricity.  See DRIFT_DAN.F
c  for the criteria.
c  WARNING - BUILT FOR SPEED : DOES NOT CHECK HOW WELL THE ORIGINAL
c  EQUATION IS SOLVED! (CAN DO THAT IN THE CALLING ROUTINE BY
c  CHECKING HOW CLOSE (x - ec*s +es*(1.-c) - dm) IS TO ZERO.
c
c  Input:
c      dm      ==> increment in mean anomaly M (real*8 scalar)
c      es,ec       ==> ecc. times sin and cos of E_0 (real*8 scalars)
c
c       Output:
c            x          ==> solution to Kepler's difference eqn (real*8 scalar)
c            s,c        ==> sin and cosine of x (real*8 scalars)
c

        subroutine drift_kepmd(dm,es,ec,x,s,c)

      implicit none

c...    Inputs
      real*8 dm,es,ec

c...  Outputs
      real*8 x,s,c

c...    Internals
      real*8 A0, A1, A2, A3, A4
        parameter(A0 = 39916800.d0, A1 = 6652800.d0, A2 = 332640.d0)
      parameter(A3 = 7920.d0, A4 = 110.d0)
      real*8 dx
      real*8 fac1,fac2,q,y
      real*8 f,fp,fpp,fppp


c...    calc initial guess for root
      fac1 = 1.d0/(1.d0 - ec)
      q = fac1*dm
      fac2 = es*es*fac1 - ec/3.d0
      x = q*(1.d0 -0.5d0*fac1*q*(es -q*fac2))

c...  excellent approx. to sin and cos of x for small x.
      y = x*x
      s = x*(A0-y*(A1-y*(A2-y*(A3-y*(A4-y)))))/A0
      c = sqrt(1.d0 - s*s)

c...    Compute better value for the root using quartic Newton method
        f = x - ec*s + es*(1.-c) - dm
        fp = 1. - ec*c + es*s
        fpp = ec*s + es*c
        fppp = ec*c - es*s
        dx = -f/fp
        dx = -f/(fp + 0.5*dx*fpp)
        dx = -f/(fp + 0.5*dx*fpp + 0.16666666666666666*dx*dx*fppp)
        x = x + dx

c...  excellent approx. to sin and cos of x for small x.
        y = x*x
        s = x*(A0-y*(A1-y*(A2-y*(A3-y*(A4-y)))))/A0
        c = sqrt(1.d0 - s*s)

      return
      end

c-------------------------------------------------------------------
c
c*************************************************************************
c                        DRIFT_KEPU_LAG.F
c*************************************************************************
c subroutine for solving kepler's equation in universal variables.
c using LAGUERRE'S METHOD
c
c             Input:
c                 s             ==>  inital value of universal variable
c                 dt            ==>  time step (real scalor)
c                 r0            ==>  Distance between `Sun' and paritcle
c                                     (real scalor)
c                 mu            ==>  Reduced mass of system (real scalor)
c                 alpha         ==>  energy (real scalor)
c                 u             ==>  angular momentun  (real scalor)
c             Output:
c                 s             ==>  final value of universal variable
c                 fp            ==>  f' from p170
c                                       (real scalors)
c                 c1,c2,c3      ==>  c's from p171-172
c                                       (real scalors)
c                 iflgn          ==>  =0 if converged; !=0 if not
c
c Author:  Hal Levison
c Date:    2/3/93
c Last revision: 4/21/93

      subroutine drift_kepu_lag(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)

      include 'swift.inc'

c...  Inputs:
      real*8 s,dt,r0,mu,alpha,u

c...  Outputs:
      real*8 fp,c1,c2,c3
      integer iflg

c...  Internals:
      integer nc,ncmax
      real*8 ln
      real*8 x,fpp,ds,c0,f
      real*8 fdt

      integer NTMP
      parameter(NTMP=NLAG2+1)

c----
c...  Executable code

c...    To get close approch needed to take lots of iterations if alpha<0
        if(alpha.lt.0.0) then
           ncmax = NLAG2
        else
           ncmax = NLAG2
        endif

        ln = 5.0
c...    start laguere's method
        do nc =0,ncmax
           x = s*s*alpha
           call drift_kepu_stumpff(x,c0,c1,c2,c3)
           c1 = c1*s
           c2 = c2*s*s
           c3 = c3*s*s*s
           f = r0*c1 + u*c2 + mu*c3 - dt
           fp = r0*c0 + u*c1 + mu*c2
           fpp = (-40.0*alpha + mu)*c1 + u*c0
           ds = - ln*f/(fp + dsign(1.d0,fp)*sqrt(abs((ln - 1.0)*
     &       (ln - 1.0)*fp*fp - (ln - 1.0)*ln*f*fpp)))
           s = s + ds

           fdt = f/dt

c..        quartic convergence
           if( fdt*fdt.lt.DANBYB*DANBYB) then
             iflg = 0
             return
           endif
c...      Laguerre's method succeeded
        enddo

        iflg = 2

        return

        end    

c-------------------------------------------------------------------
c
c*************************************************************************
c                        DRIFT_KEPU_STUMPFF.F
c*************************************************************************
c subroutine for the calculation of stumpff functions
c see Danby p.172  equations 6.9.15
c
c             Input:
c                 x             ==>  argument
c             Output:
c                 c0,c1,c2,c3   ==>  c's from p171-172
c                                       (real scalors)
c Author:  Hal Levison
c Date:    2/3/93
c Last revision: 2/3/93
c Modified by JEC: 31/3/98
c
      subroutine drift_kepu_stumpff(x,c0,c1,c2,c3)

      include 'swift.inc'

c...  Inputs:
      real*8 x

c...  Outputs:
      real*8 c0,c1,c2,c3

c...  Internals:
      integer n,i
      real*8 xm,x2,x3,x4,x5,x6

c----
c...  Executable code

      n = 0
      xm = 0.1
      do while(abs(x).ge.xm)
         n = n + 1
         x = x * .25d0
      enddo
c
      x2 = x  * x
      x3 = x  * x2
      x4 = x2 * x2
      x5 = x2 * x3
      x6 = x3 * x3
c
      c2 = 1.147074559772972d-11*x6 - 2.087675698786810d-9*x5
     %   + 2.755731922398589d-7*x4  - 2.480158730158730d-5*x3
     %   + 1.388888888888889d-3*x2  - 4.166666666666667d-2*x + .5d0
c
      c3 = 7.647163731819816d-13*x6 - 1.605904383682161d-10*x5
     %   + 2.505210838544172d-8*x4  - 2.755731922398589d-6*x3
     %   + 1.984126984126984d-4*x2  - 8.333333333333333d-3*x
     %   + 1.666666666666667d-1
c
      c1 = 1. - x*c3
      c0 = 1. - x*c2
c
      if(n.ne.0) then
         do i=n,1,-1
            c3 = (c2 + c0*c3)*.25d0
            c2 = c1*c1*.5d0
            c1 = c0*c1
            c0 = 2.*c0*c0 - 1.
            x = x * 4.
          enddo
       endif

       return
       end     

c----------------------------------------------------------------------
c
c*************************************************************************
c                        DRIFT_KEPU_FCHK.F
c*************************************************************************
c Returns the value of the function f of which we are trying to find the root
c in universal variables.
c
c             Input:
c                 dt            ==>  time step (real scalar)
c                 r0            ==>  Distance between `Sun' and particle
c                                     (real scalar)
c                 mu            ==>  Reduced mass of system (real scalar)
c                 alpha         ==>  Twice the binding energy (real scalar)
c                 u             ==>  Vel. dot radial vector (real scalar)
c                 s             ==>  Approx. root of f
c             Output:
c                 f             ==>  function value ( = 0 if O.K.) (integer)
c
c Author:  Martin Duncan
c Date:    March 12/93
c Last revision: March 12/93

      subroutine drift_kepu_fchk(dt,r0,mu,alpha,u,s,f)

c...  Inputs:
      real*8 dt,r0,mu,alpha,u,s

c...  Outputs:
      real*8 f

c...  Internals:
      real*8  x,c0,c1,c2,c3

c----
c...  Executable code

        x=s*s*alpha
        call drift_kepu_stumpff(x,c0,c1,c2,c3)
        c1=c1*s
        c2 = c2*s*s
        c3 = c3*s*s*s
        f = r0*c1 + u*c2 + mu*c3 - dt

        return
        end     

c-----------------------------------------------------------------------
c
c*************************************************************************
c                        DRIFT_KEPU_NEW.F
c*************************************************************************
c subroutine for solving kepler's equation in universal variables.
c using NEWTON'S METHOD
c
c             Input:
c                 s             ==>  inital value of universal variable
c                 dt            ==>  time step (real scalor)
c                 r0            ==>  Distance between `Sun' and paritcle
c                                     (real scalor)
c                 mu            ==>  Reduced mass of system (real scalor)
c                 alpha         ==>  energy (real scalor)
c                 u             ==>  angular momentun  (real scalor)
c             Output:
c                 s             ==>  final value of universal variable
c                 fp            ==>  f' from p170
c                                       (real scalors)
c                 c1,c2,c3      ==>  c's from p171-172
c                                       (real scalors)
c                 iflgn          ==>  =0 if converged; !=0 if not
c
c Author:  Hal Levison
c Date:    2/3/93
c Last revision: 4/21/93
c Modified by JEC: 31/3/98

      subroutine drift_kepu_new(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflgn)

      include 'swift.inc'

c...  Inputs:
      real*8 s,dt,r0,mu,alpha,u

c...  Outputs:
      real*8 fp,c1,c2,c3
      integer iflgn

c...  Internals:
      integer nc
      real*8 x,c0,ds,s2
      real*8 f,fpp,fppp,fdt

c----
c...  Executable code

      do nc=0,6
         s2 = s * s
         x = s2*alpha
         call drift_kepu_stumpff(x,c0,c1,c2,c3)
         c1 = c1*s
         c2 = c2*s2
         c3 = c3*s*s2
         f = r0*c1 + u*c2 + mu*c3 - dt
         fp = r0*c0 + u*c1 + mu*c2
         fpp = (mu - r0*alpha)*c1 + u*c0
         fppp = (mu - r0*alpha)*c0 - u*alpha*c1
         ds = - f/fp
         ds = - f/(fp + .5d0*ds*fpp)
         ds = -f/(fp + .5d0*ds*fpp + ds*ds*fppp*.1666666666666667d0)
         s = s + ds
         fdt = f/dt

c..      quartic convergence
         if( fdt*fdt.lt.DANBYB*DANBYB) then
             iflgn = 0
             return
         endif
c...     newton's method succeeded

        enddo

c..     newton's method failed
        iflgn = 1
        return

        end  

c-------------------------------------------------------------------
c
c*************************************************************************
c                        DRIFT_KEPU_GUESS.F
c*************************************************************************
c Initial guess for solving kepler's equation using universal variables.
c
c             Input:
c                 dt            ==>  time step (real scalor)
c                 r0            ==>  Distance between `Sun' and paritcle
c                                     (real scalor)
c                 mu            ==>  Reduced mass of system (real scalor)
c                 alpha         ==>  energy (real scalor)
c                 u             ==>  angular momentun  (real scalor)
c             Output:
c                 s             ==>  initial guess for the value of
c                                    universal variable
c
c Author:  Hal Levison & Martin Duncan
c Date:    3/12/93
c Last revision: April 6/93
c Modified by JEC: 8/6/98

      subroutine drift_kepu_guess(dt,r0,mu,alpha,u,s)

      include 'swift.inc'

c...  Inputs:
      real*8 dt,r0,mu,alpha,u

c...  Inputs and Outputs:
      real*8 s

c...  Internals:
      integer iflg
      real*8 y,sy,cy,sigma,es
      real*8 x,a
      real*8 en,ec,e

c----
c...  Executable code

        if (alpha.gt.0.0) then
c...       find initial guess for elliptic motion

            if( dt/r0 .le. 0.4)  then
              s = dt/r0 - (dt*dt*u)/(2.0*r0*r0*r0)
         return
            else
              a = mu/alpha
              en = sqrt(mu/(a*a*a))
              ec = 1.0 - r0/a
              es = u/(en*a*a)
              e = sqrt(ec*ec + es*es)
              y = en*dt - es
c
              call mco_sine (y,sy,cy)
c
              sigma = dsign(1.d0,(es*cy + ec*sy))
              x = y + sigma*.85*e
              s = x/sqrt(alpha)
       endif

        else
c...       find initial guess for hyperbolic motion.
      call drift_kepu_p3solve(dt,r0,mu,alpha,u,s,iflg)
      if(iflg.ne.0) then
         s = dt/r0
      endif
        endif

        return
        end     

c----------------------------------------------------------------------
c
c*************************************************************************
c                        DRIFT_KEPU_P3SOLVE.F
c*************************************************************************
c Returns the real root of cubic often found in solving kepler
c problem in universal variables.
c
c             Input:
c                 dt            ==>  time step (real scalar)
c                 r0            ==>  Distance between `Sun' and paritcle
c                                     (real scalar)
c                 mu            ==>  Reduced mass of system (real scalar)
c                 alpha         ==>  Twice the binding energy (real scalar)
c                 u             ==>  Vel. dot radial vector (real scalar)
c             Output:
c                 s             ==>  solution of cubic eqn for the
c                                    universal variable
c                 iflg          ==>  success flag ( = 0 if O.K.) (integer)
c
c Author:  Martin Duncan
c Date:    March 12/93
c Last revision: March 12/93

      subroutine drift_kepu_p3solve(dt,r0,mu,alpha,u,s,iflg)

c...  Inputs:
      real*8 dt,r0,mu,alpha,u

c...  Outputs:
      integer iflg
      real*8 s

c...  Internals:
      real*8 denom,a0,a1,a2,q,r,sq2,sq,p1,p2

c----
c...  Executable code

      denom = (mu - alpha*r0)/6.d0
      a2 = 0.5*u/denom
      a1 = r0/denom
      a0 =-dt/denom

      q = (a1 - a2*a2/3.d0)/3.d0
      r = (a1*a2 -3.d0*a0)/6.d0 - (a2**3)/27.d0
      sq2 = q**3 + r**2

      if( sq2 .ge. 0.d0) then
      sq = sqrt(sq2)

      if ((r+sq) .le. 0.d0) then
         p1 =  -(-(r + sq))**(1.d0/3.d0)
      else
         p1 = (r + sq)**(1.d0/3.d0)
      endif
      if ((r-sq) .le. 0.d0) then
         p2 =  -(-(r - sq))**(1.d0/3.d0)
      else
         p2 = (r - sq)**(1.d0/3.d0)
      endif

      iflg = 0
      s = p1 + p2 - a2/3.d0

      else
      iflg = 1
      s = 0
      endif

      return
      end     

