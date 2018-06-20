      subroutine trilinear(c1,c2,c3,c4,c5,c6,c7,c8,cout,logg,teff,feh,
     . loggmin,loggmax,teffmin,teffmax,fehmin,fehmax)
      implicit none
      real*8 c1,c2,c3,c4,c5,c6,c7,c8,cout,t,u,v
      real*8 loggmin,loggmax,teffmin,teffmax,fehmin,fehmax
      real*8 logg,teff,feh

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
