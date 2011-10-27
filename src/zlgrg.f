












c=======================================================================
c///////////////////////////  SUBROUTINE ZLGRG  \\\\\\\\\\\\\\\\\\\\\\\c
      subroutine zlgrg(d,v,w,pb,wb)
c
c  SOLVES THE LAGRANGEAN CONSERVATION LAWS USING FLUXES FROM R.S.
c
c  written by: Jim Stone
c  date:       January, 1991
c  modified1:
c
c  PURPOSE:  Updates the conservation laws in Lagrangean form using
c    fluxes in the z-direction computed by the Riemann solver.
c
c  INPUT ARGUMENTS: d,e,u,v,w,pb,wb
c
c  OUTPUT ARGUMENTS: d,e,w
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
      real d(jn,kn),v(jn,kn),w(jn,kn)
     &   ,wb(jn,kn),pb(jn,kn)
      integer j,k
      real qa,qb,qc,qd,qawt,qbwt,nu,di(kn),forcez1,forcez2
     &	   ,vlag,wlag,gravlag
      real grav,gravp,gravm,zcentlag(kn)
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
      do 100 j=js,je
c
c  Compute diffusive fluxes
c
      if (idiff .ne. 0) then
      di(ks-1) = 1.0/d(j,ks-1)
      do 10 k=ks,ke+1
	if (nyz.gt.1) then
	  qbwt=0.25*radius(j)*(dz(k)+dz(k-1))
     &	       /(0.5*(dy(j+1)+dy(j-1))+dy(j))
	  qb=(radius(j-1)*(v(j-1,k)+v(j-1,k-1))
     &	     -radius(j+1)*(v(j+1,k)+v(j+1,k-1)))/radius(j)
	  qb=qb*qbwt
        else
	  qb=0.0
        endif
	qc=w(j,k-1)-w(j,k)
        qd = max(0.0,(qb+qc))
        nu =0.1*min(d(j,k),d(j,k-1))*qd
        di(k) = 1.0/d(j,k)
        wb (j,k) = wb(j,k)           + nu*(di(    k)-di(    k-1))
        pb (j,k) = pb(j,k)           - nu*(w (j,k)-w (j,k-1))
10    continue
      endif
c
c  Compute new Lagrangean grid positions
c
      do 30 k=ks,ke+1
        zn(j,k) = z(k) + dt*wb(j,k)/radius(j)
30    continue
c     inner and outer ghost zone edges assigned using periodic symmetry
      zn(j,ks-1) = zn(j,ks  ) - (zn(j,ke+1)-zn(j,ke  ))
      zn(j,ks-2) = zn(j,ks-1) - (zn(j,ke  )-zn(j,ke-1))
      zn(j,ks-3) = zn(j,ks-2) - (zn(j,ke-1)-zn(j,ke-2))
      zn(j,ke+2) = zn(j,ke+1) + (zn(j,ks+1)-zn(j,ks  ))
      zn(j,ke+3) = zn(j,ke+2) + (zn(j,ks+2)-zn(j,ks+1))
      zn(j,ke+4) = zn(j,ke+3) + (zn(j,ks+3)-zn(j,ks+2))
      do 40 k=ks-3,ke+3
        dzn(j,k) = zn(j,k+1) - zn(j,k)
40    continue
c
      do 45 k=ks,ke
	zcentlag(k)=(zn(j,k)+zn(j,k+1))/2.0
45    continue
c
c  Update conservation laws
c
      do 50 k=ks,ke
        qa = dt/(d(j,k)*dz(k)*radius(j))
        w(j,k) = w(j,k) - qa*( pb(j,k+1)- pb(j,k))
	grav=barfract*barzcen(j,k)+spifract*spzcen(j,k)


	forcez1=grav-(w(j,k)*v(j,k)/radius(j))
c       interpolate vlag,wlag,and gravlag
	gravp=barfract*barzcen(j,k+1)+spifract*spzcen(j,k+1)


	gravm=barfract*barzcen(j,k-1)+spifract*spzcen(j,k-1)


	if (zcentlag(k).gt.zcenteul(k)) then
	  vlag=v(j,k)+((v(j,k+1)-v(j,k))*
     &         (zcentlag(k)-zcenteul(k))/(zcenteul(k+1)-zcenteul(k)))
	  wlag=w(j,k)+((w(j,k+1)-w(j,k))*
     &         (zcentlag(k)-zcenteul(k))/(zcenteul(k+1)-zcenteul(k)))
	  gravlag=grav+((gravp-grav)*
     &         (zcentlag(k)-zcenteul(k))/(zcenteul(k+1)-zcenteul(k)))
	else
	  vlag=v(j,k)+((v(j,k)-v(j,k-1))*
     &	       (zcentlag(k)-zcenteul(k))/(zcenteul(k)-zcenteul(k-1)))
	  wlag=w(j,k)+((w(j,k)-w(j,k-1))*
     &	       (zcentlag(k)-zcenteul(k))/(zcenteul(k)-zcenteul(k-1)))
	  gravlag=grav+((grav-gravm)*
     &	       (zcentlag(k)-zcenteul(k))/(zcenteul(k)-zcenteul(k-1)))
        endif
	forcez2=gravlag-(wlag*vlag/radius(j))
	w(j,k)=w(j,k)+0.5*dt*(forcez1+forcez2)
        d(j,k) = d(j,k)*dz(k)/dzn(j,k)
50    continue
c
100   continue
      return
      end
