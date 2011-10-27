

c -*- fortran -*-













c=======================================================================
c////////////////////////  CMHOG MAIN PROGRAM  \\\\\\\\\\\\\\\\\\\\\\\\c
c  written by: Jim Stone
c  date:       January, 1991
c  modified1:  -- see CVS logs of module cmhog2 --
c
c  PURPOSE:  Implementation of Colella and Woodwards Lagrange plus
c  remap PPM scheme for gas dynamics.  Designed to work on vector and
c  massively parallel machines.
c
c  HISTORY: none
c
c  DOCUMENTATION: see cmhog2/doc/hydro.tex
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
      integer iswres,iswhdf,iswhst,iswmov,j,k
      real zcs,arg

c doosu cpu time estimate
      real etime
      real elapsed(2)
      real time_step
      real t_start, t_end
      real time_self
      integer n_estimate
      common /selftime/ time_self
c doosu end

      external mstart,dataio,solver,nudt,intchk,empty
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
      write(6,"('Welcome to CMHOG version 0.9 [teuben 28-sep-2003]')")
      iswres = 0
      iswhdf = 1
      iswhst = 1
      iswmov = 1
      call mstart
      open(unit=51,file='time_'//id,status='unknown')
c      open (unit=778,file='star_'//id//'.dat',status='unknown',
c     1      form='unformatted')
      call dataio(iswres,iswhdf,iswhst,iswmov)
      write(6,"('Setup complete with ',i2,' warning(s):')") nwarn
      write(6,"(' entering main loop...')") 
c
c doosu cpu time estimate
c      n_estimate = 1
c      open(unit=52,file='cpuTime'//id//'.dat',status='unknown')
c for restart
c      time = 0.3
1000  continue
c      t_start = etime(elapsed)

        barfract=time/bartime
	barfract=min(1.0,barfract)

        spifract=barfract

c        barfract =0.

        call solver
        nhy  = nhy + 1
        time = time + dt
c       re-calculate spiral potential, if needed
	if (spamp.gt.0.0) then
        do 1010 k=ks-1,ke+1
        do 1010 j=js-3,je+3
          if (y(j).gt.sr) then
            arg=2.0*((log(y(j)-sr)/tan(spang))-
     &          (zcenteul(k)-vspiral*time))+pc
            spyed(j,k)=-spamp*(y(j)-sr)**2*exp(-spsc*(y(j)-sr))*
     &                ((3.0-spsc*(y(j)-sr))*cos(arg)-
     &                (2.0/tan(spang))*sin(arg))
          else
            spyed(j,k)=0.0
          endif
          if (radius(j).gt.sr) then
            arg=2.0*((log(radius(j)-sr)/tan(spang))-
     &          (zcenteul(k)-vspiral*time))+pc
            spycen(j,k)=-spamp*(radius(j)-sr)**2*
     &                 exp(-spsc*(radius(j)-sr))*
     &                 ((3.0-spsc*(radius(j)-sr))*cos(arg)-
     &                 (2.0/tan(spang))*sin(arg))
            spzcen(j,k)=-2.0*spamp*(radius(j)-sr)**2*
     &                 exp(-spsc*(radius(j)-sr))*sin(arg)
            arg=2.0*((log(radius(j)-sr)/tan(spang))-
     &		(z(k)-vspiral*time))+pc
            spzed(j,k)=-2.0*spamp*(radius(j)-sr)**2*
     &                 exp(-spsc*(radius(j)-sr))*sin(arg)
          else
            spycen(j,k)=0.0
            spzed(j,k)=0.0
            spzcen(j,k)=0.0
          endif
1010    continue
	endif
        call intchk(iswres,iswhdf,iswhst)
        if (ifsen .eq. 1) goto 2000
        call dataio(iswres,iswhdf,iswhst,iswmov)
        call empty
        call nudt
   
c      if (n_estimate .gt. 50) goto 1500
c      t_end = etime(elapsed)
c      time_step = t_end - t_start
c      write(52,*) n_estimate, time_step, time_self
c      n_estimate = n_estimate + 1
      goto 1000

c1500  continue
c      close(unit=52)
c      print*,'Estimation of CPU TIME complete!'
c      stop

2000  continue
      iswres = 1
      iswhdf = 1
      iswhst = 1
      iswmov = 1
      call dataio(iswres,iswhdf,iswhst,iswmov)
      zcs = float(nhy*nyz*nzz)/tused
      write(6,"('Execution terminated with ',i4,' warning(s)')") nwarn
      write(6,"('zone-cycles per cpu second =',1pe12.5)") zcs
      close(unit=51)
c      close(unit=778)
      end
