c-----------------------------------------------------------------------
c  FIELD VARIABLES
c  d = mass density [gm/cm**3]
c  u,v,w = x,y,z components of velocity [cm/sec]
c  p = pressure [dynes]
c
c  d_init = initial surface density 
      real d(jn,kn),v(jn,kn),w(jn,kn)
#ifdef SELFGRAV
      real d_init(jn,kn), d_st(jn,kn)
#endif
      common /fieldr/ d, v, w
      common /selfinit/ d_init, d_st
  
