c-----------------------------------------------------------------------
c  PARAMETERS - Always set by CPP macros
c  in,jn,kn = number of array elements in i,j,k direction
c  tiny[huge] = smallest[biggest] number allowed (machine dependent)
c
      integer jn, kn, ijkn,icartx,icarty
#ifdef SELFGRAV
      integer N1,N2,hN1,hN2,NCOPY
#ifdef CARTESIAN
      integer Nxx,Nyy,Nxy,dims,npad
#endif
#endif
      real tiny,huge
      parameter(jn= NO_OF_J_ZONES, kn= NO_OF_K_ZONES)
      parameter(tiny= A_SMALL_NUMBER, huge= A_BIG_NUMBER)
      parameter(ijkn=MAX_OF_IJK)
      parameter(icartx=IMG_X, icarty=IMG_Y)
#ifdef SELFGRAV
      parameter(N1=jn-6,N2=kn-6,hN1=N1/2,hN2=N2/2,NCOPY=2)
#ifdef CARTESIAN
c      parameter(dims=2,npad=4,Nxx=npad*(jn-6),Nyy=Nxx,Nxy=Nxx*Nyy)
      parameter(dims=2,npad=2 ,Nxx=2**12,Nyy=Nxx,Nxy=Nxx*Nyy)
#endif
#endif
