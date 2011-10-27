












c=======================================================================
c//////////////////////////  SUBROUTINE SETUP  \\\\\\\\\\\\\\\\\\\\\\\\c
      subroutine setup
c
c  SETS UP CONTROL PARAMETERS, GRID, BC, EOS, AND IC FOR NEW RUN
c
c  written by: Jim Stone
c  date:       January, 1991
c  modified1:
c
c  PURPOSE:
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
c  BOUNDARY VARIABLES
c  niib = integer flag for every zone at Inner I Boundary
c  diib = constant value of d at IIB (used for niib=3); same for e,u,v,w
c    same for OIB, IJB, OJB, IKB, OKB
c
      integer  nijb(kn), nikb(jn)
     &        ,nojb(kn), nokb(jn)
      common /bndryi/ nijb, nikb, nojb, nokb
c
      real
     &   dijb(kn), vijb(kn), wijb(kn)
     &  ,dojb(kn), vojb(kn), wojb(kn)
     &  ,dikb(jn), vikb(jn), wikb(jn)
     &  ,dokb(jn), vokb(jn), wokb(jn)
c
      common /bndryr/  dijb, vijb, wijb
     &                ,dojb, vojb, wojb
     &                ,dikb, vikb, wikb
     &                ,dokb, vokb, wokb
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
      integer j,k
      external ggen,galaxy,nudt
      namelist /hycon/ nlim,tlim,tsave,ttotal,co,pmin,idiff,ifltn,istp
     &   ,dfloor,vfloor,wfloor
      namelist /ijb/ nijb, dijb, vijb, wijb
      namelist /ojb/ nojb, dojb, vojb, wojb
      namelist /ikb/ nikb, dikb, vikb, wikb
      namelist /okb/ nokb, dokb, vokb, wokb
      namelist /eos/ ciso
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
      nlim   = 1 000 000
      tlim   = 0.0
      tsave  = 60.0
      ttotal = 100.0*(3600.0)
      co     = 0.5
      pmin   = tiny
      idiff  = 1
      ifltn  = 0
      istp   = 0
      dfloor = tiny
      vfloor = tiny
      wfloor = tiny
      read (1,hycon)
      write(2,hycon)
c
      call ggen
c
      do 30 k=1,kn
          nijb(k) = 2
          dijb(k) = dfloor
          vijb(k) = vfloor
          wijb(k) = wfloor
30    continue
      read (1,ijb)
      write(2,ijb)
c
      do 40 k=1,kn
          nojb(k) = 2
          dojb(k) = dfloor
          vojb(k) = vfloor
          wojb(k) = wfloor
40    continue
      read (1,ojb)
      write(2,ojb)
c
      do 50 j=1,jn
          nikb(j) = 4
          dikb(j) = dfloor
          vikb(j) = vfloor
          wikb(j) = wfloor
50    continue
      read (1,ikb)
      write(2,ikb)
c
      do 60 j=1,jn
          nokb(j) = 4
          dokb(j) = dfloor
          vokb(j) = vfloor
          wokb(j) = wfloor
60    continue
      read (1,okb)
      write(2,okb)
c
      read (1,eos)
      write(2,eos)
c
c  Initialize all field variables to floor values and call problem
c
      do 70 k=1,kn
      do 70 j=1,jn
        d(j,k) = dfloor
        v(j,k) = vfloor
        w(j,k) = wfloor
70    continue
      call galaxy
c
c  Initial timestep
c
      dt = huge
      call nudt
      dtmin = dt/100.0
c
      return
      end
