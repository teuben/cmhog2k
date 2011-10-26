c
c	on some machines this arithmetic statement function
c	will work, on strict ansi (e.g. the f2c compiler) you
c	will need the cvmgt.src function, on the cray CVMGT is
c	an inline vector instruction
c
#if defined(SUN) && !defined(LINUX)
      real aa,bb,cvmgt
      integer ii
      cvmgt(aa,bb,ii) = ii*aa + (1-ii)*bb
#endif SUN
#if defined(CONVEX)
      real aa,bb,cvmgt
      integer ii
      cvmgt(aa,bb,ii) = -ii*aa + (1+ii)*bb
#endif CONVEX
