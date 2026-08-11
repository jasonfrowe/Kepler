      subroutine occultquad(z0,u1,u2,p,muo1,mu0,nz)
C  This routine computes the lightcurve for occultation
C  of a quadratically limb-darkened source without microlensing.
C  Please cite Mandel & Agol (2002) if you make use of this routine
C  in your research.  Please report errors or bugs to agol@tapir.caltech.edu
      implicit none
      integer i,nz
      double precision z0(nz),u1,u2,p,muo1(nz),mu0(nz),
     &       mu(nz),lambdad(nz),etad(nz),lambdae(nz),lam,
     &       pi,x1,x2,x3,z,omega,kap0,kap1,q,kc,Kk,Ek,Pk,n,cel
      if(abs(p-0.5d0).lt.1.d-3) p=0.5d0
C
C Input:
C
C rs   radius of the source (set to unity)
C z0   impact parameter in units of rs
C p    occulting star size in units of rs
C u1   linear    limb-darkening coefficient (gamma_1 in paper)
C u2   quadratic limb-darkening coefficient (gamma_2 in paper)
C
C Output:
C
C muo1 fraction of flux at each z0 for a limb-darkened source
C mu0  fraction of flux at each z0 for a uniform source
C
C Limb darkening has the form:
C  I(r)=[1-u1*(1-sqrt(1-(r/rs)^2))-u2*(1-sqrt(1-(r/rs)^2))^2]/(1-u1/3-u2/6)/pi
C 
C To use this routine
C
C Now, compute pure occultation curve:
      omega=1.d0-u1/3.d0-u2/6.d0
      pi=acos(-1.d0)
C Loop over each impact parameter:
      do i=1,nz
C substitute z=z0(i) to shorten expressions
        z=z0(i)
        x1=(p-z)**2
        x2=(p+z)**2
        x3=p**2-z**2
C the source is unocculted:
C Table 3, I.
        if(z.ge.1.d0+p) then
          lambdad(i)=0.d0
          etad(i)=0.d0
          lambdae(i)=0.d0
          goto 10
        endif
C the  source is completely occulted:
C Table 3, II.
        if(p.ge.1.d0.and.z.le.p-1.d0) then
          lambdad(i)=1.d0
          etad(i)=1.d0
          lambdae(i)=1.d0
          goto 10
        endif
C the source is partly occulted and the occulting object crosses the limb:
C Equation (26):
        if(z.ge.abs(1.d0-p).and.z.le.1.d0+p) then
          kap1=acos(min((1.d0-p*p+z*z)/2.d0/z,1.d0))
          kap0=acos(min((p*p+z*z-1.d0)/2.d0/p/z,1.d0))
          lambdae(i)=p*p*kap0+kap1
          lambdae(i)=(lambdae(i)-0.5d0*sqrt(max(4.d0*z*z-
     &               (1.d0+z*z-p*p)**2,0.d0)))/pi
        endif
C the occulting object transits the source star (but doesn't
C completely cover it):
        if(z.le.1.d0-p) lambdae(i)=p*p
C the edge of the occulting star lies at the origin- special 
C expressions in this case:
        if(abs(z-p).lt.1.d-4*(z+p)) then
C Table 3, Case V.:
          if(z.ge.0.5d0) then
            lam=0.5d0*pi
            q=0.5d0/p
            kc=sqrt(max(0.d0,1.d0-q*q))
            Kk=cel(kc,1.d0,1.d0,1.d0)
            Ek=cel(kc,1.d0,1.d0,kc*kc)
C Equation 34: lambda_3
            lambdad(i)=1.d0/3.d0+16.d0*p/9.d0/pi*(2.d0*p*p-1.d0)*Ek-
     &                 (32.d0*p**4-20.d0*p*p+3.d0)/9.d0/pi/p*Kk
C Equation 34: eta_1
            etad(i)=1.d0/2.d0/pi*(kap1+p*p*(p*p+2.d0*z*z)*kap0-
     &              (1.d0+5.d0*p*p+z*z)/4.d0*sqrt((1.d0-x1)*(x2-1.d0)))
            if(p.eq.0.5d0) then
C Case VIII: p=1/2, z=1/2
              lambdad(i)=1.d0/3.d0-4.d0/pi/9.d0
              etad(i)=3.d0/32.d0
            endif
            goto 10
          else
C Table 3, Case VI.:
            lam=0.5d0*pi
            q=2.d0*p
            kc=sqrt(max(0.d0,1.d0-q*q))
            Kk=cel(kc,1.d0,1.d0,1.d0)
            Ek=cel(kc,1.d0,1.d0,kc*kc)
C Equation 34: lambda_4
            lambdad(i)=1.d0/3.d0+2.d0/9.d0/pi*(4.d0*(2.d0*p*p-1.d0)*Ek+
     &                 (1.d0-4.d0*p*p)*Kk)
C Equation 34: eta_2
            etad(i)=p*p/2.d0*(p*p+2.d0*z*z)
            goto 10
          endif
        endif
C the occulting star partly occults the source and crosses the limb:
C Table 3, Case III:
        if((z.gt.0.5d0+abs(p-0.5d0).and.z.lt.1.d0+p).or.(p.gt.0.5d0.
     &      and.z.gt.abs(1.d0-p)*1.0001d0.and.z.lt.p)) then
          lam=0.5d0*pi
          q=sqrt((1.d0-(p-z)**2)/4.d0/z/p)
          if(q.gt.1.0d0) q=0.99999d0
          kc=sqrt(max(0.d0,1.d0-q*q))
          Kk=cel(kc,1.d0,1.d0,1.d0)
          Ek=cel(kc,1.d0,1.d0,kc*kc)
          n=1.d0/x1-1.d0
          Pk=cel(kc,1.d0+n,1.d0,1.d0)
C Equation 34, lambda_1:
          lambdad(i)=1.d0/9.d0/pi/sqrt(p*z)*(((1.d0-x2)*(2.d0*x2+
     &        x1-3.d0)-3.d0*x3*(x2-2.d0))*Kk+4.d0*p*z*(z*z+
     &        7.d0*p*p-4.d0)*Ek-3.d0*x3/x1*Pk)
          if(z.lt.p) lambdad(i)=lambdad(i)+2.d0/3.d0
C Equation 34, eta_1:
          etad(i)=1.d0/2.d0/pi*(kap1+p*p*(p*p+2.d0*z*z)*kap0-
     &          (1.d0+5.d0*p*p+z*z)/4.d0*sqrt((1.d0-x1)*(x2-1.d0)))
          goto 10
        endif
C the occulting star transits the source:
C Table 3, Case IV.:
        if(p.le.1.d0.and.z.le.(1.d0-p)*1.0001d0) then
          lam=0.5d0*pi
          q=sqrt((x2-x1)/(1.d0-x1))
          if(q.gt.1.0d0) q=0.99999d0
          kc=sqrt(max(0.d0,1.d0-q*q))
          Kk=cel(kc,1.d0,1.d0,1.d0)
          Ek=cel(kc,1.d0,1.d0,kc*kc)
          n=x2/x1-1.d0
          Pk=cel(kc,1.d0+n,1.d0,1.d0)
C Equation 34, lambda_2:
          lambdad(i)=2.d0/9.d0/pi/sqrt(1.d0-x1)*((1.d0-5.d0*z*z+p*p+
     &         x3*x3)*Kk+(1.d0-x1)*(z*z+7.d0*p*p-4.d0)*Ek-3.d0*x3/x1*Pk)
          if(z.lt.p) lambdad(i)=lambdad(i)+2.d0/3.d0
          if(abs(p+z-1.d0).le.1.d-4) then
            lambdad(i)=2/3.d0/pi*acos(1.d0-2.d0*p)-4.d0/9.d0/pi*
     &            sqrt(p*(1.d0-p))*(3.d0+2.d0*p-8.d0*p*p)
          endif
C Equation 34, eta_2:
          etad(i)=p*p/2.d0*(p*p+2.d0*z*z)
        endif
 10     continue
C Now, using equation (33):
        muo1(i)=1.d0-((1.d0-u1-2.d0*u2)*lambdae(i)+(u1+2.d0*u2)*
     &      lambdad(i)+u2*etad(i))/omega
C Equation 25:
        mu0(i)=1.d0-lambdae(i)
      enddo
      return
      end

      function cel(kc, p, a, b)
      implicit none
      double precision cel, kc, p, a, b
      double precision ca, e, f, g, h, m, p_loc, a_loc, b_loc, kc_loc
      parameter (ca=1.0d-12)

      if (kc .eq. 0.0d0) then
         cel = 0.0d0
         return
      endif

      kc_loc = abs(kc)
      p_loc = p
      a_loc = a
      b_loc = b
      e = kc_loc
      m = 1.0d0

      if (p_loc .gt. 0.0d0) then
         p_loc = sqrt(p_loc)
         b_loc = b_loc / p_loc
      else
         f = kc_loc * kc_loc
         g = 1.0d0 - p_loc
         e = (f - p_loc) / g
         h = (b_loc - a_loc * p_loc) / g
         p_loc = sqrt(e / g)
         a_loc = (a_loc - b_loc) / g
         b_loc = -h / (p_loc * g) + a_loc * p_loc
      endif

 10   continue
      f = a_loc
      a_loc = a_loc + b_loc / p_loc
      g = e / p_loc
      b_loc = b_loc + f * g
      b_loc = b_loc + b_loc
      p_loc = g + p_loc
      g = m
      m = m + kc_loc
      if (abs(g - kc_loc) .gt. g * ca) then
         kc_loc = 2.0d0 * sqrt(e)
         e = kc_loc * m
         goto 10
      endif
      cel = (a_loc * m + b_loc) / (m * (m + p_loc)) * 1.5707963267948966d0
      return
      end
