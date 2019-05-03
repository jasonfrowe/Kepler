c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MAL_HCON.FOR    (ErikSoft   28 March 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Does an integration using an integrator with a constant stepsize H.
c Input and output to this routine use coordinates XH, and velocities VH,
c with respect to the central body, but the integration algorithm uses
c its own internal coordinates X, and velocities V.
c
c The programme uses the transformation routines COORD and BCOORD to change
c to and from the internal coordinates, respectively.
c
c------------------------------------------------------------------------------
c
      subroutine mal_hcon (time,tstart,tstop,dtout,algor,h0,tol,jcen,
     %  rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,
     %  stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,lmem,onestep,
     %  coord,bcoord)
c
      implicit none
      include 'mercury.inc'
c
c Input/Output
      integer algor,nbod,nbig,stat(nbod),opt(8),opflag,ngflag
      integer lmem(NMESS),ndump,nfun
      real*8 time,tstart,tstop,dtout,h0,tol,jcen(3),rcen,rmax
      real*8 en(3),am(3),cefac,m(nbod),xh(3,nbod),vh(3,nbod)
      real*8 s(3,nbod),rho(nbod),rceh(nbod),ngf(4,nbod)
      character*8 id(nbod)
      character*80 outfile(3),dumpfile(4),mem(NMESS)
c
c Local
      integer i,j,k,n,itmp,nclo,nhit,jhit(CMAX),iclo(CMAX),jclo(CMAX)
      integer dtflag,ejflag,stopflag,colflag,nstored
      real*8 x(3,NMAX),v(3,NMAX),xh0(3,NMAX),vh0(3,NMAX)
      real*8 rce(NMAX),rphys(NMAX),rcrit(NMAX),epoch(NMAX)
      real*8 hby2,tout,tmp0,tdump,tfun,tlog,dtdump,dtfun
      real*8 dclo(CMAX),tclo(CMAX),dhit(CMAX),thit(CMAX)
      real*8 ixvclo(6,CMAX),jxvclo(6,CMAX),a(NMAX)
      external onestep,coord,bcoord
c
c------------------------------------------------------------------------------
c
c Initialize variables. DTFLAG = 0/2: first call ever/normal
      dtout  = abs(dtout)
      dtdump = abs(h0) * ndump
      dtfun  = abs(h0) * nfun
      dtflag = 0
      nstored = 0
      hby2 = 0.500001d0 * abs(h0)
c
c Calculate close-encounter limits and physical radii
      call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig,
     %  m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),1)
c
c Set up time of next output, times of previous dump, log and periodic effect
      if (opflag.eq.-1) then
        tout = tstart
      else
        n = int (abs (time-tstart) / dtout) + 1
        tout = tstart  +  dtout * sign (dble(n), tstop - tstart)
        if ((tstop-tstart)*(tout-tstop).gt.0) tout = tstop
      end if
      tdump = time
      tfun  = time
      tlog  = time
c
c Convert to internal coordinates and velocities
      call coord (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag,opt)
c
c------------------------------------------------------------------------------
c
c  MAIN  LOOP  STARTS  HERE
c
 100  continue
c
c Is it time for output ?
      if (abs(tout-time).le.hby2.and.opflag.ge.-1) then
c
c Beware: the integration may change direction at this point!!!!
        if (opflag.eq.-1.and.dtflag.ne.0) dtflag = 1
c
c Convert to heliocentric coordinates and output data for all bodies
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        call mio_out (time,jcen,rcen,rmax,nbod,nbig,m,xh,vh,s,rho,
     %    stat,id,opt,opflag,algor,outfile(1))
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id,
     %    0,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem,
     %    outfile,nstored,0)
        tmp0 = tstop - tout
        tout = tout + sign( min( abs(tmp0), abs(dtout) ), tmp0 )
c
c Update the data dump files
        do j = 2, nbod
          epoch(j) = time
        end do
        call mio_dump (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen,
     %    rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat,
     %    id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)
        tdump = time
      end if
c
c If integration has finished, convert to heliocentric coords and return
      if (abs(tstop-time).le.hby2.and.opflag.ge.0) then
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        return
      end if
c
c Make sure the integration is heading in the right direction
 150  continue
      tmp0 = tstop - time
      if (opflag.eq.-1) tmp0 = tstart - time
      h0 = sign (h0, tmp0)
c
c Save the current heliocentric coordinates and velocities
      if (algor.eq.1) then
        call mco_iden (time,jcen,nbod,nbig,h0,m,x,v,xh0,vh0,ngf,ngflag,
     %    opt)
      else
        call bcoord(time,jcen,nbod,nbig,h0,m,x,v,xh0,vh0,ngf,ngflag,opt)
      end if
      call onestep (time,tstart,h0,tol,rmax,en,am,jcen,rcen,nbod,nbig,
     %  m,x,v,s,rphys,rcrit,rce,stat,id,ngf,algor,opt,dtflag,ngflag,
     %  opflag,colflag,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,outfile,
     %  mem,lmem)
      time = time + h0
c
c------------------------------------------------------------------------------
c
c  CLOSE  ENCOUNTERS
c
c If encounter minima occurred, output details and decide whether to stop
      if (nclo.gt.0.and.opflag.ge.-1) then
        itmp = 1
        if (colflag.ne.0) itmp = 0
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id,nclo,
     %    iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem,
     %    outfile,nstored,itmp)
        if (stopflag.eq.1) return
      end if
c
c------------------------------------------------------------------------------
c
c  COLLISIONS
c
c If collisions occurred, output details and remove lost objects
      if (colflag.ne.0) then
c
c Reindex the surviving objects
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat,
     %    id,mem,lmem,outfile(3),itmp)
c
c Reset flags, and calculate new Hill radii and physical radii
        dtflag = 1
        if (opflag.ge.0) opflag = 1
        call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig,
     %    m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),1)
        call coord (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag,opt)
      end if
c
c------------------------------------------------------------------------------
c
c  COLLISIONS  WITH  CENTRAL  BODY
c
c Check for collisions with the central body
      if (algor.eq.1) then
        call mco_iden(time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
      else
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
      end if
      itmp = 2
      if (algor.eq.11.or.algor.eq.12) itmp = 3
      call mce_cent (time,h0,rcen,jcen,itmp,nbod,nbig,m,xh0,vh0,xh,vh,
     %  nhit,jhit,thit,dhit,algor,ngf,ngflag)
c
c If something hit the central body, restore the coords prior to this step
      if (nhit.gt.0) then
        call mco_iden (time,jcen,nbod,nbig,h0,m,xh0,vh0,xh,vh,ngf,
     %    ngflag,opt)
        time = time - h0
c
c Merge the object(s) with the central body
        do k = 1, nhit
          i = 1
          j = jhit(k)
          call mce_coll (thit(k),tstart,en(3),jcen,i,j,nbod,nbig,m,xh,
     %      vh,s,rphys,stat,id,opt,mem,lmem,outfile(3))
        end do
c
c Remove lost objects, reset flags and recompute Hill and physical radii
        call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat,
     %    id,mem,lmem,outfile(3),itmp)
        if (opflag.ge.0) opflag = 1
        dtflag = 1
        call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig,
     %    m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),0)
        if (algor.eq.1) then
          call mco_iden (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag,
     %      opt)
        else
          call coord (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag,opt)
        end if
c
c Redo that integration time step
        goto 150
      end if
c
c------------------------------------------------------------------------------
c
c  DATA  DUMP  AND  PROGRESS  REPORT
c
c Convert to heliocentric coords and do the data dump
      if (abs(time-tdump).ge.abs(dtdump).and.opflag.ge.-1) then
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        do j = 2, nbod
          epoch(j) = time
        end do
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id,
     %    0,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem,
     %    outfile,nstored,0)
        call mio_dump (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen,
     %    rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat,
     %    id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)
        tdump = time
      end if
c
c Convert to heliocentric coords and write a progress report to the log file
      if (abs(time-tlog).ge.abs(dtdump).and.opflag.ge.0) then
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        call mxx_en (jcen,nbod,nbig,m,xh,vh,s,en(2),am(2))
        call mio_log (time,tstart,en,am,opt,mem,lmem)
        tlog = time
      end if
c
c------------------------------------------------------------------------------
c
c  CHECK  FOR  EJECTIONS  AND  DO  OTHER  PERIODIC  EFFECTS
c
      if (abs(time-tfun).ge.abs(dtfun).and.opflag.ge.-1) then
        if (algor.eq.1) then
          call mco_iden (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,
     %      opt)
        else
          call bcoord(time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        end if
c
c Recompute close encounter limits, to allow for changes in Hill radii
        call mce_hill (nbod,m,xh,vh,rce,a)
        do j = 2, nbod
          rce(j) = rce(j) * rceh(j)
        end do
c
c Check for ejections
        itmp = 2
        if (algor.eq.11.or.algor.eq.12) itmp = 3
        call mxx_ejec (time,tstart,rmax,en,am,jcen,itmp,nbod,nbig,m,xh,
     %    vh,s,stat,id,opt,ejflag,outfile(3),mem,lmem)
c
c Remove ejected objects, reset flags, calculate new Hill and physical radii
        if (ejflag.ne.0) then
          call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat,
     %      id,mem,lmem,outfile(3),itmp)
          if (opflag.ge.0) opflag = 1
          dtflag = 1
          call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig,
     %      m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),0)
          if (algor.eq.1) then
            call mco_iden (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,
     %        ngflag,opt)
          else
            call coord (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag,
     %        opt)
          end if
        end if
        tfun = time
      end if
c
c Go on to the next time step
      goto 100
c
c------------------------------------------------------------------------------
c
      end
