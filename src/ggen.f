












c=======================================================================
c/////////////////////////  SUBROUTINE GGEN  \\\\\\\\\\\\\\\\\\\\\\\\\\c
      subroutine ggen
c
c  INITIALIZES GRID ZONES
c
c     written by: Jim Stone
c     date:       January, 1991
c     modified1: 
c
c  PURPOSE:  Initializes the grid in a new run according to the control
c  parameters in the input deck namelists "ggen1","ggen2","ggen3".
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
      real psdymin,psdymax,psdyrat,psddymin,psdzmin,psdzmax,psddzmin
      integer npsd
      parameter (npsd=5000)
      real psdy(npsd),psdz(npsd)
      integer psdyind(npsd),psdzind(npsd)
      common /pseudo/ psdymin,psdymax,psdyrat,psddymin,psdzmin
     &               ,psdzmax,psddzmin,psdy,psdz,psdyind,psdzind
 

      integer  nbl,igrid,imin,imax,iter,i,j,k
      real  ymin, ymax, yrat, dymin, dfndyr, yr, deltyr, erryr
     &     ,zmin, zmax, zrat, dzmin, dfndzr, zr, deltzr, errzr, fn
      logical  lgrid
      namelist /ggen2/ nbl,ymin,ymax,igrid,yrat,dymin,vgy,lgrid
      namelist /ggen3/ nbl,zmin,zmax,igrid,zrat,dzmin,vgz,lgrid
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////
c=======================================================================
c-----------  GENERATE Y GRID  -----------------------------------------
c  Read in blocks of y-grid zones.  Variables same as x-grid blocks.
c
      js   = 4
      imax = js
c      namelist /ggen2/ nbl,ymin,ymax,igrid,yrat,dymin,vgy,lgrid
      nbl   = 1
      ymin  = 0.0
      ymax  = 0.0
      igrid = 1
      yrat  = 1.0
      dymin = 0.0
      vgy   = 0.0
      lgrid = .false.
c
210   continue
      read (1,ggen2)
      write(2,ggen2)
c
      imin = imax
      imax = imax + nbl
      if (imin.eq.js .and. nbl.eq.1 .and. lgrid.eqv..true.) then
        js   = 1
        imin = 1
        imax = 2
      else
        if (imax-1 .gt. jn-3) then
      write(6,"('ERROR: number of zones in y-direction exceeds array')")
      write(6,"(' bounds',/1x,'imax = ',i4,'  jn = ',i4)") imax,jn
        stop
        endif
      endif
c
c  1). Compute dy(imin) from given value of yrat.
c
      if (igrid .eq. 1) then
        if (yrat .eq. 1.0)  then
          dy(imin) = (ymax-ymin)/float(nbl)
        else
          dy(imin) = (ymax-ymin)*(yrat-1.0)/(yrat**nbl - 1.0)
        endif
      endif
c
c  2). Compute yrat from given value of dymin.  Newton Raphson iteration
c  is required to find the root (yrat) of the function:
c     fn(yr) = (ymax-ymin) - dy(imin)*[(yr)**nbl - 1]/[yr-1] = 0
c  apparently the initial guess for yrat is rather essential for 
c  large nbl. In old code you may find yr=1.01, which doesn't work.
c
      if (igrid .eq. 2) then
        dy(imin) = dymin
        yr = (ymax/ymin)**(1.0/nbl)
        do 220 iter=1,20
          fn = (ymax - ymin) - dymin*(yr**nbl - 1.0)/(yr - 1.0)
          dfn dyr =  -nbl*dymin*yr**(nbl - 1)/(yr - 1.0)
     &               + dymin*(yr**nbl - 1.0)/(yr - 1.0)**2
          deltyr  = -fn/dfndyr
          err yr  = abs(deltyr/yr)
          yr = yr + deltyr
          write(*,*) 'iter,yrat,yrat_err ',iter,yr,erryr
          if (erryr .lt. 1.0e-6) goto 230
220     continue
        write(6,"('ERROR from GGEN: Newton-Raphson did not converge')")
        write(6,"('for yrat',/1x,'imin=',i3,' yr=',1pe12.5)") imin,yr
        write(6,"('deltyr = '1e12.5,' fn = ',1e12.5)") deltyr,fn
        stop
c
230     continue
        yrat = yr
      endif
c
c  Set up Eulerian grid EDGES from i=imin to imax.  Go back and read
c  another block of y-grid zones, if needed.
c
      y(imin  ) = ymin
      y(imin+1) = ymin + dy(imin)
      do 300 i=imin+2,imax
        dy(i-1) = dy(i-2) * yrat
        y (i  ) = y (i-1) + dy(i-1)
300   continue
      if (.not. lgrid) go to 210
c
c  Setup boundary grid zones
c
      je   = imax-1
      jbtm = js
      jtp  = je
      nyz  = je-js+1
c
      if (nyz .gt. 1) then
        dy(js-1) = dy(js  )
        dy(js-2) = dy(js+1)
        dy(js-3) = dy(js+2)
        dy(je+1) = dy(je  )
        dy(je+2) = dy(je-1)
        dy(je+3) = dy(je-2)
c
        y(js-1) = y(js  ) - dy(js-1)
        y(js-2) = y(js-1) - dy(js-2)
        y(js-3) = y(js-2) - dy(js-3)
        y(je+2) = y(je+1) + dy(je+1)
        y(je+3) = y(je+2) + dy(je+2)
        y(je+4) = y(je+3) + dy(je+3)

      do 305 j=js-3,je+3
        rvoleul(j)=abs (0.5 * ((y(j+1)**2)-(y(j)**2)))
305   continue
      endif
c
      do 310 j=js,je
          radius(j)=abs((2.0*(y(j+1)**3-y(j)**3))
     & 		    /(3.0*(y(j+1)**2-y(j)**2)))
310   continue
      if (nyz.gt.1) then
        radius(js-1)=abs((2.0*(y(js)**3-y(js-1)**3))
     &		     /(3.0*(y(js)**2-y(js-1)**2)))
        radius(js-2)=abs((2.0*(y(js-1)**3-y(js-2)**3))
     &		     /(3.0*(y(js-1)**2-y(js-2)**2)))
        radius(js-3)=abs((2.0*(y(js-2)**3-y(js-3)**3))
     &		     /(3.0*(y(js-2)**2-y(js-3)**2)))
        radius(je+1)=abs((2.0*(y(je+2)**3-y(je+1)**3))
     &		     /(3.0*(y(je+2)**2-y(je+1)**2)))
        radius(je+2)=abs((2.0*(y(je+3)**3-y(je+2)**3))
     &		     /(3.0*(y(je+3)**2-y(je+2)**2)))
        radius(je+3)=abs((2.0*(y(je+4)**3-y(je+3)**3))
     &		     /(3.0*(y(je+4)**2-y(je+3)**2)))
      endif
c
c-----------  GENERATE Z GRID  -----------------------------------------
c  Read in blocks of x-grid zones.  Variables same as x-grid blocks.
c
      ks   = 4
      imax = ks
c      namelist /ggen3/ nbl,zmin,zmax,igrid,zrat,dzmin,vgz,lgrid
      nbl   = 1
      zmin  = 0.0
      zmax  = 0.0
      igrid = 1
      zrat  = 1.0
      dzmin = 0.0
      vgz   = 0.0
      lgrid = .false.
c
410   continue
      read (1,ggen3)
      write(2,ggen3)
c
      imin = imax
      imax = imax + nbl
      if (imin.eq.ks .and. nbl.eq.1 .and. lgrid.eqv..true.) then
        ks   = 1
        imin = 1
        imax = 2
      else
        if (imax-1 .gt. kn-3) then
      write(6,"('ERROR: number of zones in z-direction exceeds array')")
      write(6,"(' bounds',/1x,'imax = ',i4,'  kn = ',i4)") imax,kn
        stop
        endif
      endif
c
c  1). Compute dz(imin) from given value of zrat.
c
      if (igrid .eq. 1) then
        if (zrat .eq. 1.0)  then
          dz(imin) = (zmax-zmin)/float(nbl)
        else
          dz(imin) = (zmax-zmin)*(zrat-1.0)/(zrat**nbl - 1.0)
        endif
      endif
c
c  2). Compute zrat from given value of dzmin.  Newton Raphson iteration
c  is required to find the root (zrat) of the function:
c     fn(zr) = (zmax-zmin) - dz(imin)*[(zr)**nbl - 1]/[zr-1] = 0
c
      if (igrid .eq. 2) then
        dz(imin) = dzmin
        zr = 1.01
        do 420 iter=1,20
          fn = (zmax - zmin) - dzmin*(zr**nbl - 1.0)/(zr - 1.0)
          dfn dzr =  -nbl*dzmin*zr**(nbl - 1)/(zr - 1.0)
     &               + dzmin*(zr**nbl - 1.0)/(zr - 1.0)**2
          deltzr  = -fn/dfndzr
          err zr  = abs(deltzr/zr)
          zr = zr + deltzr
          if (errzr .lt. 1.0e-6) goto 430
420     continue
        write(6,"('ERROR from GGEN: Newton-Raphson did not converge')")
        write(6,"('for zrat',/1x,'imin=',i3,' zr=',1pe12.5)") imin,zr
        write(6,"('deltzr = ',1e12.5,' fn = ',1e12.5)") deltzr,fn
        stop
c
430     continue
        zrat = zr
      endif
c
c  Set up Eulerian grid EDGES from i=imin to imax.  Go back and read
c  another block of z-grid zones, if needed.
c
      z(imin  ) = zmin
      z(imin+1) = zmin + dz(imin)
      do 500 i=imin+2,imax
        dz(i-1) = dz(i-2) * zrat
        z (i  ) = z (i-1) + dz(i-1)
500   continue
      if (.not. lgrid) go to 410
c
c  Setup boundary grid zones
c
      ke   = imax-1
      kbtm = ks
      ktp  = ke
      nzz  = ke-ks+1
c
      if (nzz .gt. 1) then
        dz(ks-1) = dz(ks  )
        dz(ks-2) = dz(ks+1)
        dz(ks-3) = dz(ks+2)
        dz(ke+1) = dz(ke  )
        dz(ke+2) = dz(ke-1)
        dz(ke+3) = dz(ke-2)
c
        z(ks-1) = z(ks  ) - dz(ks-1)
        z(ks-2) = z(ks-1) - dz(ks-2)
        z(ks-3) = z(ks-2) - dz(ks-3)
        z(ke+2) = z(ke+1) + dz(ke+1)
        z(ke+3) = z(ke+2) + dz(ke+2)
        z(ke+4) = z(ke+3) + dz(ke+3)
      endif
c
	do 510 k=ks,ke
	  zcenteul(k)=(z(k)+z(k+1))/2.0
510     continue
        if (nzz.gt.1) then
	  zcenteul(ks-1) = (z(ks-1)+z(ks))/2.0
	  zcenteul(ke+1) = (z(ke+1)+z(ke+2))/2.0
          zcenteul(ks-2) = (z(ks-2)+z(ks-1))/2.0
          zcenteul(ke+2) = (z(ke+2)+z(ke+3))/2.0
        endif
c
      return
      end
