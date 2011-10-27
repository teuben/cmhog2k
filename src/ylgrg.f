












c=======================================================================
c///////////////////////////  SUBROUTINE YLGRG  \\\\\\\\\\\\\\\\\\\\\\\c
      subroutine ylgrg(d,v,w,pb,vb)
c
c  SOLVES THE LAGRANGEAN CONSERVATION LAWS USING FLUXES FROM R.S.
c
c  written by: Jim Stone
c  date:       January, 1991
c  modified1:
c
c  PURPOSE:  Updates the conservation laws in Lagrangean form using
c    fluxes in the y-direction computed by the Riemann solver.
c
c  INPUT ARGUMENTS: d,e,u,v,w,pb,vb
c
c  OUTPUT ARGUMENTS: d,e,v
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
     &   ,vb(jn,kn),pb(jn,kn)
      integer j,k
      real qa,qb,qc,qd,qawt,qcwt,nu,di(jn),ycentlag(jn),
     &     forcey1,forcey2,gravlag,wlag,grav,gravp,gravm
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
      do 100 k=ks,ke
c
c  Compute diffusive fluxes
c
      if (idiff .ne. 0) then
      di(js-1) = 1.0/d(js-1,k)
      do 10 j=js,je+1
        qb=(radius(j-1)*v(j-1,k)-radius(j)*v(j,k))/y(j)
        if (nzz.gt.1) then
	  qcwt=0.25*(dy(j)+dy(j-1))/(y(j)*(0.5*(dz(k+1)+dz(k-1))+dz(k)))
          qc=(w(j,k-1)+w(j-1,k-1))-(w(j,k+1)+w(j-1,k+1))
	  qc=qc*qcwt
        else
          qc=0.0
        endif                 
        qd = max(0.0,(qb+qc))
        nu =0.1*min(d(j,k),d(j-1,k))*qd
        di(j) = 1.0/d(j,k)
        vb (j,k) = vb(j,k)           + nu*(di(  j  ) -di(  j-1  ))
        pb (j,k) = pb(j,k)           - nu*(v (j,k) -v (j-1,k))
10    continue
      endif
c
c  Compute new lagrangean grid positions
c
      do 30 j=js,je+1
        yn(j,k) = y(j) + dt*vb(j,k)
30    continue
c     inner ghost zone EDGES assigned by reflection
      yn(js-1,k) = yn(js  ,k) - (yn(js+1,k)-yn(js  ,k))
      yn(js-2,k) = yn(js-1,k) - (yn(js+2,k)-yn(js+1,k))
      yn(js-3,k) = yn(js-2,k) - (yn(js+3,k)-yn(js+2,k))
c     outer ghost zone EDGES assigned by outflow
      yn(je+2,k) = yn(je+1,k) + (yn(je+1,k)-yn(je  ,k))
      yn(je+3,k) = yn(je+2,k) + (yn(je+1,k)-yn(je  ,k))
      yn(je+4,k) = yn(je+3,k) + (yn(je+1,k)-yn(je  ,k))
      do 40 j=js-3,je+3
        dyn(j,k) = yn(j+1,k) - yn(j,k)
        rvollag(j,k)=abs(0.5*((yn(j+1,k)**2)-(yn(j,k)**2)))
40    continue
c
      do 45 j=js,je
	ycentlag(j)=(2.0*(yn(j+1,k)**3-yn(j,k)**3))
     &		    /(3.0*(yn(j+1,k)**2-yn(j,k)**2))
45    continue
c
c  Update conservation laws
c
      do 50 j=js,je
        qa = dt/(d(j,k)*rvoleul(j))
	qb = (y(j)+y(j+1)+yn(j,k)+yn(j+1,k))/4.0
        v(j,k) = v(j,k) - qa*qb*( pb(j+1,k)- pb(j,k))
	grav=(1-barfract)*diskycen(j,k)+barfract*barycen(j,k)
     &       +spifract*spycen(j,k)


	forcey1=grav+(w(j,k)**2/radius(j))
c       interpolate gravity & w in center of lagrangian zone
        gravp=(1-barfract)*diskycen(j+1,k)+barfract*barycen(j+1,k)
     &        +spifract*spycen(j+1,k)


	gravm=(1-barfract)*diskycen(j-1,k)+barfract*barycen(j-1,k)
     &        +spifract*spycen(j-1,k)


        if (ycentlag(j).gt.radius(j)) then
          wlag=w(j,k)+((w(j+1,k)-w(j,k))*
     &        (ycentlag(j)-radius(j))/(radius(j+1)-radius(j)))
          gravlag=grav+((gravp-grav)*
     &        (ycentlag(j)-radius(j))/(radius(j+1)-radius(j)))
	else
	  wlag=w(j,k)+((w(j,k)-w(j-1,k))*
     &        (ycentlag(j)-radius(j))/(radius(j)-radius(j-1)))
	  gravlag=grav+((grav-gravm)*
     &        (ycentlag(j)-radius(j))/(radius(j)-radius(j-1)))
	endif
	forcey2=gravlag+(wlag**2/ycentlag(j))
	v(j,k)=v(j,k)+(0.5*dt*(forcey1+forcey2))
        d(j,k)=d(j,k)*rvoleul(j)/rvollag(j,k)
50    continue
c
100   continue
      return
      end
