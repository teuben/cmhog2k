












c=======================================================================
c/////////////////////////  SUBROUTINE SOLVER  \\\\\\\\\\\\\\\\\\\\\\\\c
      subroutine solver
c
c  PERFORMS PPM UPDATE FOR ONE TIMESTEP
c   MODIFIED VERSION FOR ADVECTION TESTS
c
c  written by: Jim Stone
c  date:       January, 1991
c  modified1:
c
c  PURPOSE:  Advances the fluid equations by one timestep using the
c    Lagrange plus remap PPM algorithm in 3-D.
c
c  LOCALS:
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c  PARAMETERS - Always set by CPP macros
c  in,jn,kn = number of array elements in i,j,k direction
c  tiny[huge] = smallest[biggest] number allowed (machine dependent)
c
      integer jn, kn, ijkn,icartx,icarty
      real tiny,huge
      parameter(jn= 262, kn= 144)
      parameter(tiny= 1.0e-10, huge= 1.0e+10)
      parameter(ijkn=1030)
      parameter(icartx=256, icarty=256)
c-----------------------------------------------------------------------
c  GRID VARIABLES
c  is,ie = array index of first and last ACTIVE zone in i-direction
c  itp,ibtm = array index of entire computational volume; may be
c    different from is,ie if only computing on a sub-volume
c  nxz = # of ACTIVE zones in x-direction = ie-is+1; same for nyz,nzz
c
      integer        js  , je , ks  , ke , nyz, nzz
     &              ,jbtm, jtp, kbtm, ktp
      common /gridi/ js  , je , ks  , ke , nyz, nzz
     &              ,jbtm, jtp, kbtm, ktp
c
c  x ,y ,z  = coordinates of Eulerian zone EDGES
c  dx,dy,dx = distance between Eulerian zone EDGES in x,y,z directions
c    e.g. dx(i)=x(i+1)-x(i)
c  xn,yn,zn = coordinates of Lagrangean zone EDGES
c  dxn,dyn,dzn = distance between Lagrangean zone EDGES
c  vgx,vgy,vgz = grid velocity applied during remap in x,y,z directions
c
      real
     &   y (jn+1),           z (kn+1)
     & , dy (jn  ),          dz (kn  )
     & , rvoleul (jn),   radius (jn  ),zcenteul(kn)
     & , yn(jn+1,kn+1), zn(jn+1,kn+1)
     & , dyn(jn  ,kn  ),dzn(jn  ,kn  )
     & , rvollag (jn,kn)
     & , diskyed(jn,kn),diskycen(jn,kn)
     & , baryed(jn,kn)
     & , barycen(jn,kn),barzed(jn,kn),barzcen(jn,kn)
     & , spyed(jn,kn),spycen(jn,kn)
     & , spzed(jn,kn),spzcen(jn,kn)
     & , bhycen(jn,kn),bhyed(jn,kn)
     & ,massflux(kn),omassflux(kn),xyz(jn+1,kn+1),dxyz(jn,kn),vgy,vgz
      equivalence (yn,zn,xyz) , (dyn,dzn,dxyz)
      common /gridr/  y,z,dy,dz,rvoleul,rvollag,radius
     &		     ,zcenteul,xyz,dxyz,vgy,vgz
     &               ,diskyed,diskycen,baryed,barycen
     &               ,barzed,barzcen,spyed,spycen
     &               ,spzed,spzcen,massflux,omassflux
     &               ,bhycen,bhyed
c-----------------------------------------------------------------------
c  FIELD VARIABLES
c  d = mass density [gm/cm**3]
c  u,v,w = x,y,z components of velocity [cm/sec]
c  p = pressure [dynes]
c
c  d_init = initial surface density 
      real d(jn,kn),v(jn,kn),w(jn,kn)
      common /fieldr/ d, v, w
      common /selfinit/ d_init, d_st
  
c-----------------------------------------------------------------------
c  ROOT VARIABLES
c  idiff = diffusive flux in Lagrangean step switch (0=off)
c  ifltn = flattener switch (0=off) 
c  ifsen = stop flag [code only runs while ifsen=0]
c  istp  = steepener switch for d interpolation (0=off)
c  nlim  = cycle limit
c  nhy   = number of hydro cycles executed
c  nwarn = number of warnings during execution
c
      integer
     &  idiff   ,ifltn   ,ifsen,istp,nlim,nhy,nwarn,nwritten
      common /rooti/
     &  idiff   ,ifltn   ,ifsen,istp,nlim,nhy,nwarn,nwritten
c
c  co    = Courant #
c  dfloor= default value of density, same for e,u,v,w
c  dt    = timestep
c  dtcs  = CFL limit for sound waves
c  dtdump= cpu time between restart dumps
c  dthdf = problem time between hdf dumps
c  dthist= problem time between history dumps
c  dtmin = minimum timestep below which warnings are issued
c  dtu   = CFL limit for u-velocity, same for v,w
c  ciso  = isothermal sound speed
c  pmin  = minimum allowed pressure
c  tdump = cpu time of last restart dump
c  time  = problem time
c  thdf  = problem time of last hdf dump
c  thist = problem time of last history dump
c  tlim  = problem time limit
c  trem  = cpu time remaining before auto-stop due to cpu time limit
c  tsave = cpu reserve time for final data dumps during auto-stop
c  ttotal= total cpu time allowed
c  tused = cpu time used
c  
      real
     &  co      ,dfloor  ,dt      ,dtcs    ,dtdump  ,dthdf   ,dthist
     & ,dtmin   ,dtmovie ,dtu     ,dtv     ,dtw     ,ciso     ,pmin
     & ,tdump   ,time    ,thdf    ,thist   ,tlim    ,trem    ,tsave
     & ,ttotal  ,tused   ,tmovie  ,vfloor  ,wfloor
     & ,bartime ,barfract,rl,vspiral, spifract
     & ,cma1,cmi1,cma2,cmi2,cma3,cmi3
     & ,spamp,spang,sppat,spsc,sr,pc
      common /rootr/
     &  co      ,dfloor  ,dt      ,dtcs    ,dtdump  ,dthdf   ,dthist
     & ,dtmin   ,dtmovie ,dtu     ,dtv     ,dtw     ,ciso     ,pmin
     & ,tdump   ,time    ,thdf    ,thist   ,tlim    ,trem    ,tsave
     & ,ttotal  ,tused   ,tmovie  ,vfloor  ,wfloor,cmax,cmin
     & ,bartime ,barfract,rl,vspiral, spifract
     & ,cma1,cmi1,cma2,cmi2,cma3,cmi3
     & ,spamp,spang,sppat,spsc,sr,pc
c
c  hdffile = hdf     file name
c  resfile = restart file name
c  id      = two letter problem tag appended to file names
c
      character*9 hdffile, resfile, movfile1,movfile2,movfile3,strfile
      character*2 id
      common /rootch/ hdffile, movfile1,movfile2,movfile3, resfile, id
     &                ,strfile
c-----------------------------------------------------------------------
c  SCRATCH SPACE
c
      real  wijk0(jn,kn), wijk1(jn,kn), wijk2(jn,kn)
     &     ,wijk3(jn,kn), wijk4(jn,kn), wijk5(jn,kn)
c
      common /scr1/ wijk0, wijk1, wijk2, wijk3, wijk4, wijk5
      real psdymin,psdymax,psdyrat,psddymin,psdzmin,psdzmax,psddzmin
      integer npsd
      parameter (npsd=5000)
      real psdy(npsd),psdz(npsd)
      integer psdyind(npsd),psdzind(npsd)
      common /pseudo/ psdymin,psdymax,psdyrat,psddymin,psdzmin
     &               ,psdzmax,psddzmin,psdy,psdz,psdyind,psdzind
 

c-----------------------------------------------------------------------
      integer j,k,ixyz,n
      real  dl(jn,kn), dr(jn,kn) 
     &     ,pb(jn,kn)
     &     ,vl(jn,kn), vr(jn,kn), vb(jn,kn)
     &     ,wl(jn,kn), wr(jn,kn), wb(jn,kn)
      equivalence (dl,wijk0),(dr,wijk1)
     &           ,(vl,wl,pb,wijk4),(vr,wr,vb,wb,wijk5)
c
      real  df(jn,kn)
     &     ,vf(jn,kn), wf(jn,kn)
      equivalence (df,wijk0)
     &           ,(vf,wijk3),(wf,wijk4)

      real*8 pi
      PARAMETER (pi = 3.141592654)

c
      external yintlgrg, ylgrg, yintrmp, yremap
     &        ,zintlgrg, zlgrg, zintrmp, zremap, isoshockquad
     &        ,bvaldyeul, bvaldylag, bvalvy, bvalwy
     &        ,bvaldz, bvalvz, bvalwz

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////
c=======================================================================
c
      call cut(d,v)



c#ifdef SELFGRAV
c      do 10 k=ks,ke
c         do 10 j=js,je+1
c            v(j,k) = v(j,k) + 0.5*dt*selfgyed (j,k)
c            w(j,k) = w(j,k) + 0.5*dt*selfgzcen(j,k)
c   10 continue
c#endif
      ixyz = mod(nhy,2)
      do 130 n=ixyz,ixyz+1
c
c  Update in y-direction
c
      if (mod(n,2) .eq. 0 .and. nyz .gt. 1) then
        call bvaldyeul
        call bvalvy
        call bvalwy
        call yintlgrg          (d,v,w,  dy,dl,dr,vl,vr)
        call isoshockquad (js,je+1,ks,ke,    dl,dr,vl,vr,pb,vb)
        do 50 k=ks,ke
          vb(js,k)=min(vb(js,k),0.0)
50      continue
        call ylgrg             (d,  v,w,                     pb,vb)
c
        call bvaldyeul
        call bvalvy
        call bvalwy
        call yintrmp(d,v,w,dyn,df,vf,wf)
        call yremap (d,v,w,    df,vf,wf)
      endif
c
c  Update in z-direction
c
      if (mod(n,2) .eq. 1 .and. nzz .gt. 1) then
        call bvaldz
        call bvalwz
        call bvalvz
        call zintlgrg          (d,  v,w,dz,dl,dr,wl,wr)
        call isoshockquad (js,je,ks,ke+1,    dl,dr,wl,wr,pb,wb)
        call zlgrg             (d,  v,w,                     pb,wb)
c
        call bvaldz
        call bvalvz
        call bvalwz
        call zintrmp(d,v,w,dzn,df,vf,wf)
        call zremap (d,v,w,    df,vf,wf)
      endif
130   continue
c
      return
      end
