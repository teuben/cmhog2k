

c -*- fortran -*-













c	 will turn on HDF output, else ascii output ala Glenn Piner
c=======================================================================
c/////////////////////////  SUBROUTINE HDFALL  \\\\\\\\\\\\\\\\\\\\\\\\c
      subroutine hdfall(filename)
c
c  MAKES AN HDF DUMP CONTAINING ALL ACTIVE FIELD ARRAYS
c
c     written by: Jim Stone
c     date:       January,1989
c     modified1:  PJT - april 1995: 3d -> 2d for the isothermal bar project
c
c  PURPOSE: Makes an hdf dump of all the active field variables.
c    Data is written in the Scientific Data Set to "filename".  Note
c    that data must be stored contiguously in order to interface
c    correctly to the C hdf routines.
c
c  EXTERNALS: HDF library routines, PGAS (only if non-isothermal)
c
c  LOCALS:
c-----------------------------------------------------------------------
c      implicit NONE
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
      real psdymin,psdymax,psdyrat,psddymin,psdzmin,psdzmax,psddzmin
      integer npsd
      parameter (npsd=5000)
      real psdy(npsd),psdz(npsd)
      integer psdyind(npsd),psdzind(npsd)
      common /pseudo/ psdymin,psdymax,psdyrat,psddymin,psdzmin
     &               ,psdzmax,psddzmin,psdy,psdz,psdyind,psdzind
 

c #include "scratch.h"
      character*9  filename
c
      integer  i,j,k,rank,shape(3),ret,nyd,nzd
      real data(jn*kn),yscale(jn),zscale(kn)
      character*32 string
      equivalence (data,wijk0)
c
      integer  dssdims,dssdast,dssdisc,dsadata,dspdata
      external dssdims,dssdast,dssdisc,dsadata,dspdata
c     &  ,pgas
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
      do 20 j=jbtm-3,jtp+3
        yscale(j-jbtm+4) = y(j) + 0.5*dy(j)
20    continue
      do 30 k=kbtm-3,ktp+3
        zscale(k-kbtm+4) = z(k) + 0.5*dz(k)
30    continue
      nyd = jtp-jbtm+7
      nzd = ktp-kbtm+7
c
c  double check, since the two versions (stone vs. piner) seem to use 
c  different indices
c
      if (nyd .ne. je-js+7) then
         write(*,*) 'Problem with Radial grid',nyd,js,je         
      endif
      if (nzd .ne. ke-ks+7) then
         write(*,*) 'Problem with Angular grid',nzd,ks,ke
      endif


c      rank     = 3
      rank     = 2
      shape(1) = nyd
      shape(2) = nzd
      ret = dssdims(rank,shape)
c      ret = dssdisc(1,shape(1),xscale)
c      ret = dssdisc(2,shape(2),yscale)
c      ret = dssdisc(3,shape(3),zscale)
      ret = dssdisc(1,shape(1),yscale)
      ret = dssdisc(2,shape(2),zscale)

c
c  y-velocity
c
      do 200 k=ks-3,ke+3
      do 200 j=js-3,je+3
        data((k-ks+3)*nyd + (j-js+3)+1) = v(j,k)
200   continue   
      write(string,"('R-VELOCITY AT TIME=',1pe8.2)") time
      ret = dssdast(string,'km/sec   ','1pe8.2','Polar')
      ret = dsadata(filename,rank,shape,data)
c
c  z-velocity
c
      do 300 k=ks-3,ke+3
      do 300 j=js-3,je+3
        data((k-ks+3)*nyd + (j-js+3)+1) = w(j,k)
300   continue   
      write(string,"('PHI-VELOCITY AT TIME=',1pe8.2)") time
      ret = dssdast(string,'km/sec   ','1pe8.2','Polar')
      ret = dsadata(filename,rank,shape,data)
c
c  density
c
      do 400 k=ks-3,ke+3
      do 400 j=js-3,je+3
        data((k-ks+3)*nyd + (j-js+3)+1) = d(j,k)
400   continue
      write(string,"('DENSITY AT TIME=',1pe8.2)") time
      ret = dssdast(string,'Msolar/pc**2','1pe8.2','Polar')
      ret = dsadata(filename,rank,shape,data)

c
c  neutral fraction
c
c

c     forces
      do 810 k=ks,ke
      do 810 j=js,je
        data((k-ks)*nyd + (j-js)+1) = diskyed(j,k)
810   continue
      write(string,"('DISKYED AT TIME=',1pe8.2)") time
      ret = dssdast(string,'force','1pe8.2','Polar')
      ret = dsadata(filename,rank,shape,data)

      do 820 k=ks,ke
      do 820 j=js,je
        data((k-ks)*nyd + (j-js)+1) = diskycen(j,k)
820   continue
      write(string,"('DISKYCEN AT TIME=',1pe8.2)") time
      ret = dssdast(string,'force','1pe8.2','Polar')
      ret = dsadata(filename,rank,shape,data)

      do 830 k=ks,ke
      do 830 j=js,je
        data((k-ks)*nyd + (j-js)+1) = baryed(j,k)
830   continue
      write(string,"('BARYED AT TIME=',1pe8.2)") time
      ret = dssdast(string,'force','1pe8.2','Polar')
      ret = dsadata(filename,rank,shape,data)

      do 840 k=ks,ke
      do 840 j=js,je
        data((k-ks)*nyd + (j-js)+1) = barycen(j,k)
840   continue
      write(string,"('BARYCEN AT TIME=',1pe8.2)") time
      ret = dssdast(string,'force','1pe8.2','Polar')
      ret = dsadata(filename,rank,shape,data)

      do 850 k=ks,ke
      do 850 j=js,je
        data((k-ks)*nyd + (j-js)+1) = barzed(j,k)
850   continue
      write(string,"('BARZED AT TIME=',1pe8.2)") time
      ret = dssdast(string,'force','1pe8.2','Polar')
      ret = dsadata(filename,rank,shape,data)

      do 860 k=ks,ke
      do 860 j=js,je
        data((k-ks)*nyd + (j-js)+1) = barzcen(j,k)
860   continue
      write(string,"('BARZCEN AT TIME=',1pe8.2)") time
      ret = dssdast(string,'force','1pe8.2','Polar')
      ret = dsadata(filename,rank,shape,data)
      return
      end
