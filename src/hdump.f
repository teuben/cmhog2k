












c=======================================================================
c////////////////////////  SUBROUTINE hdump  \\\\\\\\\\\\\\\\\\\\\\\\c
      subroutine hdump
c
c  FORMATTED WRITE OF SELECTED VARIABLES
c
c     written by: Jim Stone
c     date:       February,1993
c     modified1: 
c
c  PURPOSE:  Dumps scalar "history" variables in a formatted write
c  for analysis.  Currently implemented variables are:
c   scal( 1) = time
c   scal( 2) = dt
c   scal( 3) = mass
c   scal( 4) = 0.5*d*v**2
c   scal( 5) = 0.5*d*w**2
c  More variables can be added by increasing nscal, adding the
c  appropriate lines to the do 10 loop, and changing the 2001 format
c  statement
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
      integer j,k,nscal,irl
      parameter(nscal=11)
      real scal(nscal),dvol,darea,qa,massrl,grav
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
      qa =(z(ke+1)-z(ks)) * 0.5*(y(je+1)**2-y(js)**2)
      do 1 i=1,nscal
        scal(i) = 0.0
1     continue
      massrl=0.0
      scal(1) = time
      scal(2) = dt
c
c  Integrate quantities
c
      do 10 k=ks,ke
      do 10 j=js,je
        dvol = (dz(k)*0.5*(y(j+1)**2-y(j)**2))/qa
        scal(3) = scal(3) + dvol*d(j,k)
        scal(4) = scal(4) + dvol*d(j,k)*0.5*v(j,k)**2
        scal(5) = scal(5) + dvol*d(j,k)*0.5*w(j,k)**2
        grav=(1-barfract)*diskycen(j,k)
     &       +barfract*barycen(j,k) + spifract*spycen(j,k)
        scal(6) = scal(6) + dvol*d(j,k)*v(j,k)*
     &            (w(j,k)-sqrt(abs(radius(j)*grav)))
10    continue
      do 11 j=js,je
        if (y(j+1).gt.rl) then
          irl=j-1
          goto 12
        endif
11    continue
12    qa = (z(ke+1)-z(ks))*0.5*(y(irl+1)**2-y(js)**2)
      do 13 k=ks,ke
      do 13 j=js,irl
        dvol = (dz(k)*0.5*(y(j+1)**2-y(j)**2))/qa
        massrl=massrl+dvol*d(j,k)
      scal(7)=scal(7)+dvol*d(j,k)*v(j,k)
13    continue
      scal(7)=scal(7)/massrl
c    
      qa = (z(ke+1)-z(ks))
      do 15 k=ks,ke
	darea = dz(k)/qa
      scal(8) = scal(8) + darea*massflux(k)
      scal(9) = scal(9) + darea*omassflux(k)
15    continue
c
c  Write out variables to file connected to unit 3 opened in MAIN
c  program unit (zeus2d.src)
c
      write(3,2001) (scal(i),i=1,nscal)
2001  format(1x,1pe12.5,10e14.5)
c
      return
      end

