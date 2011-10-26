#ifdef SELFGRAV
      real H0,Ro_self
      real phi(jn,kn),coeff,tsgst,tsgfl
      real selfgyed(jn,kn),selfgycen(jn,kn)
     &    ,selfgzed(jn,kn),selfgzcen(jn,kn)
#ifdef MILLERPOT
      complex*8 FTDpot
c      real       DDD
      dimension FTDpot(N1,NCOPY*N2/2+1,N1)
c      dimension FTDpot(N1,NCOPY*N2/2+1,N1),DDD(NCOPY*N2)

      integer*8 planf,plani
      real massn,VV
      complex*8 FTmassn,FTV
      dimension massn(N2),FTmassn(hN2+1),VV(N2),FTV(hN2+1)
#endif /* MILLERPOT */
#ifdef COORDTRANS 
      real amp_k(NCOPY*hN1+1, N2)
      integer*8 planfnm,planinm
      real mass
      dimension mass(NCOPY*N1,N2)
      complex*8 out
      dimension out(NCOPY*hN1+1,N2)
#endif /* COORDTRANS */
#ifdef CARTESIAN
      real selfyrat, selfymin, selfdymin, selfzmin, selfdzmin
      real Lx,Ly,Lx2,Ly2
      real din(Nxx,Nyy), pin(Nxx,Nyy)
      real xc(Nxx), yc(Nyy)
      real xy_rdd(Nxx,Nyy), xy_thd(Nxx,Nyy)
     &    ,xy_frd(Nxx,Nyy), xy_fth(Nxx,Nyy)

      real rth_xd(jn,kn),rth_yd(jn,kn),rth_fx(jn,kn),rth_fy(jn,kn)

      integer inNxx, inNyy

      integer*8 p_fwd, p_rvs
      real indata
      dimension indata(Nxx,Nyy)
      complex*8 out
      dimension out(Nxx/2+1,Nyy)
#endif /* CARTESIAN */

      common /diskthick/ H0,Ro_self
      common /selfgrav/ phi,coeff,tsgst,tsgfl,selfgyed,selfgycen
     &                 ,selfgzed,selfgzcen
#ifdef MILLERPOT
      common /poten/ FTDpot,planf,plani,massn,FTmassn,VV,FTV
c      common /poten/ FTDpot,planf,plani,massn,FTmassn,VV,FTV,DDD
#endif /*  MILLERPOT */
#ifdef COORDTRANS 
      common /pureF/ amp_k,planfnm,planinm,mass,out
#endif /* COORDTRANS */

#ifdef CARTESIAN
      common /gridconf/ selfyrat, selfymin, selfdymin
     &                          , selfzmin, selfdzmin
      common /boxsize/Lx,Ly,Lx2,Ly2
      common /selfint/ din,pin
      common /xyGrid/ xc,yc,xy_rdd,xy_thd,xy_frd,xy_fth,inNxx,inNyy
      common /rthGrid/ rth_xd, rth_yd, rth_fx, rth_fy

      common /fftplan/ p_fwd, p_rvs, indata, out 
#endif /* CARTESIAN */


#endif /* SELFGRAV */
