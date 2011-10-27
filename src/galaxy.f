

c -*- fortran -*-













c=======================================================================
c////////////////////////  SUBROUTINE GALAXY  \\\\\\\\\\\\\\\\\\\\\\\\\c  defines the potential from either an analytical subroutine
c  (it should define inipotential() and potential() as defined
c   by NEMO), or via up to two grid (FITS) files defining the
c   axisymmetric and  non-axisymmetric parts of the potential
c   on an X-Y grid.
c
c
c
      subroutine galaxy
c
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
c  FIELD VARIABLES
c  d = mass density [gm/cm**3]
c  u,v,w = x,y,z components of velocity [cm/sec]
c  p = pressure [dynes]
c
c  d_init = initial surface density 
      real d(jn,kn),v(jn,kn),w(jn,kn)
      common /fieldr/ d, v, w
      common /selfinit/ d_init, d_st
  
      real psdymin,psdymax,psdyrat,psddymin,psdzmin,psdzmax,psddzmin
      integer npsd
      parameter (npsd=5000)
      real psdy(npsd),psdz(npsd)
      integer psdyind(npsd),psdzind(npsd)
      common /pseudo/ psdymin,psdymax,psdyrat,psddymin,psdzmin
     &               ,psdzmax,psddzmin,psdy,psdz,psdyind,psdzind
 

      INTEGER FFTW_R2HC
      PARAMETER (FFTW_R2HC=0)
      INTEGER FFTW_HC2R
      PARAMETER (FFTW_HC2R=1)
      INTEGER FFTW_DHT
      PARAMETER (FFTW_DHT=2)
      INTEGER FFTW_REDFT00
      PARAMETER (FFTW_REDFT00=3)
      INTEGER FFTW_REDFT01
      PARAMETER (FFTW_REDFT01=4)
      INTEGER FFTW_REDFT10
      PARAMETER (FFTW_REDFT10=5)
      INTEGER FFTW_REDFT11
      PARAMETER (FFTW_REDFT11=6)
      INTEGER FFTW_RODFT00
      PARAMETER (FFTW_RODFT00=7)
      INTEGER FFTW_RODFT01
      PARAMETER (FFTW_RODFT01=8)
      INTEGER FFTW_RODFT10
      PARAMETER (FFTW_RODFT10=9)
      INTEGER FFTW_RODFT11
      PARAMETER (FFTW_RODFT11=10)
      INTEGER FFTW_FORWARD
      PARAMETER (FFTW_FORWARD=-1)
      INTEGER FFTW_BACKWARD
      PARAMETER (FFTW_BACKWARD=+1)
      INTEGER FFTW_MEASURE
      PARAMETER (FFTW_MEASURE=0)
      INTEGER FFTW_DESTROY_INPUT
      PARAMETER (FFTW_DESTROY_INPUT=1)
      INTEGER FFTW_UNALIGNED
      PARAMETER (FFTW_UNALIGNED=2)
      INTEGER FFTW_CONSERVE_MEMORY
      PARAMETER (FFTW_CONSERVE_MEMORY=4)
      INTEGER FFTW_EXHAUSTIVE
      PARAMETER (FFTW_EXHAUSTIVE=8)
      INTEGER FFTW_PRESERVE_INPUT
      PARAMETER (FFTW_PRESERVE_INPUT=16)
      INTEGER FFTW_PATIENT
      PARAMETER (FFTW_PATIENT=32)
      INTEGER FFTW_ESTIMATE
      PARAMETER (FFTW_ESTIMATE=64)
      INTEGER FFTW_ESTIMATE_PATIENT
      PARAMETER (FFTW_ESTIMATE_PATIENT=128)
      INTEGER FFTW_BELIEVE_PCOST
      PARAMETER (FFTW_BELIEVE_PCOST=256)
      INTEGER FFTW_DFT_R2HC_ICKY
      PARAMETER (FFTW_DFT_R2HC_ICKY=512)
      INTEGER FFTW_NONTHREADED_ICKY
      PARAMETER (FFTW_NONTHREADED_ICKY=1024)
      INTEGER FFTW_NO_BUFFERING
      PARAMETER (FFTW_NO_BUFFERING=2048)
      INTEGER FFTW_NO_INDIRECT_OP
      PARAMETER (FFTW_NO_INDIRECT_OP=4096)
      INTEGER FFTW_ALLOW_LARGE_GENERIC
      PARAMETER (FFTW_ALLOW_LARGE_GENERIC=8192)
      INTEGER FFTW_NO_RANK_SPLITS
      PARAMETER (FFTW_NO_RANK_SPLITS=16384)
      INTEGER FFTW_NO_VRANK_SPLITS
      PARAMETER (FFTW_NO_VRANK_SPLITS=32768)
      INTEGER FFTW_NO_VRECURSE
      PARAMETER (FFTW_NO_VRECURSE=65536)
      INTEGER FFTW_NO_SIMD
      PARAMETER (FFTW_NO_SIMD=131072)

      integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
      parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)

c#include "fftw_f77.i"
c-----------------------------------------------------------------------
      integer i,j,k,ndim,npar
      real arg
      character name*5, pot0*80, pot1*80, errtxt*40,comment*80
      real crval1,crval2,crpix1,crpix2,cdelt1,cdelt2
      real*8 amp,amode,n,aob,qm,rhoc,gheight,
     &                r_th,gam,bhmass
      real*8 dlimit,rcloud,starff
      integer istar
      real*8 par(11),pos(3),acc(3),pot,dumtime
      integer flun,bsize,status,rwmode,naxes(2)
      logical anyf
c doosu thick
c      real gheight

      real*8 pi
      PARAMETER (pi = 3.141592654)

      INTEGER nx,ny
      REAL gpot0,gpot1,gcrpix1,gcrpix2,gcdelt1,gcdelt2,gcrval1,gcrval2
      REAL gomega,grtstart,grtscale
      COMMON/ggrid/nx,ny,gpot0(1030,1030),gpot1(1030,1030),
     *          grmax,gomega,
     * 		gcrpix1,gcrpix2,gcdelt1,gcdelt2,gcrval1,gcrval2,
     *          grhalo,gvhalo,ggamma,grtstart,grtscale

      external inipotential,potential
c doosu spiral r position
      integer sr_ind
c doosu end

      real*8 pot_out(jn,kn),diskf(jn,kn),bulf(jn,kn),barf(jn,kn)
      REAL*8 FXBUL,FYBUL,FXDISK,FYDISK
      COMMON /CHECKINIT/ FXBUL,FYBUL,FXDISK,FYDISK
      REAL*8 FXBAR,FYBAR
      COMMON /CHECKINIT2/ FXBAR,FYBAR
      integer idum,idum2

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////
c=======================================================================
c
c     assumed that igrav=1 and igeom=1
      namelist /pgen/ amp,npar,amode,n,aob,rl,qm,rhoc,bartime,gheight,
     &                r_th,gam,bhmass
      common /diskthick / gheight
      namelist /spiral/ spamp,spsc,spang,sppat,sr,pc
      namelist /pgrid/pot0,pot1,naxis1,naxis2,rmax,omega,gamma,
     *                crval1,crval2,crpix1,crpix2,cdelt1,cdelt2,
     *                rhalo,vhalo,rtstart,rtscale,rnstart,rnscale,
     *                nrh,nth
      ndim=3
      npar=10
      name='dummy'
      pot0=' '
      pot1=' '
      naxis1=0
      naxis2=0
      crval1=0.0
      crval2=0.0
      cdelt1=1.0
      cdelt2=1.0
      crpix1=1.0
      crpix2=1.0
      rmax=0.0
      omega=0.0
      gamma=1.0
      rhalo=1.0
      vhalo=0.0
      rtstart=1000.0
      rtscale=1000.0

      rnstart=-0.3
      rnscale=0.05

      nrh=0
      nth=0

      write(2,*) 'PINER MODEL'
      read (1,pgen)
      write(2,pgen)
      read (1,spiral)
      write(2,spiral)
      read (1,pgrid)
      write(2,pgrid)

c wooyoung test
c      open(unit=777,file='test_Nm_a2',status='unknown')
c doosu start
      open(unit=50,file='coord_'//id,status='unknown')
c      open(unit=55,file='square',status='unknown')
      do j=1,jn
c         write(50,*) j, y(j)
         write(50,*) j, radius(j)
c         write(55,*) j, y(j),dy(j),y(j)*dz(ks),dy(j)-y(j)*dz(ks)
      enddo
c       print*,js,je,y(js),y(je)
      close(unit=50)
c      close(unit=55)
      open(unit=51,file='coordall_'//id,status='unknown'
     &    ,form='unformatted')
      write(51) (j,j=1,jn),(y(j),j=1,jn),(radius(j),j=1,jn)
      write(51) (k,k=1,kn),(z(k),k=1,kn),(zcenteul(k),k=1,kn)
      close(unit=51)
c      stop
c doosu end


C
C  there are two modes of using PGRID:    when pot0 given, it will 
C  use the external file (a fits file), else it will use an internal
C  grid, i.e. pot0/pot1 is filled on the fly from the analytical
C  potential
C
      if (pot0 .ne. ' ') then
	write(*,*) 'Support for FITSIO was not enabled - '
	write(*,*) 'You probably should edit cmhog.def and Makefile'
	STOP
C
C  the send mode of using 'pgrid' is setting the number of pixels of the grid,
C  and an internal grid is computed exactly from the analytical potential
C  NEMO routines can also generate those externally on a grid, but this method
C  makes the program independant of NEMO
C
      else if (naxis1.gt.0 .and. naxis2.gt.0) then
					
	write(*,*) 'Using an internal grid; size=',naxis1,naxis2,
     *             ' rmax=',rmax
	nx = naxis1
	ny = naxis2
	gcrpix1 = crpix1
	gcrpix2 = crpix2
	gcrval1 = crval1
	gcrval2 = crval2
	gcdelt1 = cdelt1
	gcdelt2 = cdelt2
	grmax = rmax
	gomega = omega
        ggamma = gamma
	write(*,*) 'Pattern speed omega = ',omega
      endif
C
      do 5 i=1,npar
	if (i.eq.1) par(i)=omega
	if (i.eq.2) par(i)=amode
	if (i.eq.3) par(i)=n
	if (i.eq.4) par(i)=aob
	if (i.eq.5) par(i)=rl
	if (i.eq.6) par(i)=0.0
	if (i.eq.7) par(i)=rhoc
	if (i.eq.8) par(i)=r_th
	if (i.eq.9) par(i)=gam
	if (i.eq.10) par(i)=bhmass
5     continue
c        print out the grid values, both center and edges
      if (.false.) then
         do  101 j=js-3,je+3
            write(*,*) 'J: ',j,radius(j),y(j)
 101     continue
         do 102 k=ks-1,ke+1
            write(*,*) 'K: ',k,zcenteul(k),z(k)
 102     continue

      endif

c   	first set par(6) = 0, i.e. determine the radial (Y) forces 
c 	when no bar is present.

      call inipotential(npar,par,name)
C
      if (nyz.gt.1) then
	do 7 j=js-3,je+3
	 do 6 k=ks,ke+1
          pos(1)=y(j)*cos(zcenteul(k))
	  pos(2)=y(j)*sin(zcenteul(k))
	  pos(3)=0.0
	  call potential (ndim,pos,acc,pot,dumtime)
	  diskyed(j,k)=(acc(1)*cos(zcenteul(k)))
     &			 +(acc(2)*sin(zcenteul(k)))
	  pos(1)=radius(j)*cos(zcenteul(k))
	  pos(2)=radius(j)*sin(zcenteul(k))
	  call potential (ndim,pos,acc,pot,dumtime)
	  diskycen(j,k)=(acc(1)*cos(zcenteul(k)))
     &			  +(acc(2)*sin(zcenteul(k)))
6        continue
        
         if (y(j).lt.rnstart) 
     &    call taper(diskyed,jn,j,ks,ke+1,y(j),rnstart,rnscale)
         if (radius(j).lt.rnstart) 
     &    call taper(diskycen,jn,j,ks,ke+1,radius(j),rnstart,rnscale)

         if (y(j).gt.rtstart) 
     &    call taper(diskyed,jn,j,ks,ke+1,y(j),rtstart,rtscale)
         if (radius(j).gt.rtstart) 
     &    call taper(diskycen,jn,j,ks,ke+1,radius(j),rtstart,rtscale)

7       continue
      endif
c
c   	second, set par(6) = qm, i.e. determine the radial (Y) and
c	tangential (Z) forces when a bar is present.

      par(6)=qm
      call inipotential(npar,par,name)
      if (nyz.gt.1) then
	do 9 j=js-3,je+3
	 do 8 k=ks,ke+1
          pos(1)=y(j)*cos(zcenteul(k))
	  pos(2)=y(j)*sin(zcenteul(k))
	  pos(3)=0.0
	  call potential (ndim,pos,acc,pot,dumtime)
c doosu print potential
          pot_out(j,k) = pot
          barf(j,k) = (FYBAR*cos(zcenteul(k)))
     &                  +(FXBAR*sin(zcenteul(k)))
          diskf(j,k) = (FYDISK*cos(zcenteul(k)))
     &                    +(FXDISK*sin(zcenteul(k)))
          bulf(j,k)  = (FYBUL*cos(zcenteul(k)))
     &                    +(FXBUL*sin(zcenteul(k)))
c doosu end
	  baryed(j,k)=(acc(1)*cos(zcenteul(k)))
     &			+(acc(2)*sin(zcenteul(k)))
	  pos(1)=radius(j)*cos(zcenteul(k))
	  pos(2)=radius(j)*sin(zcenteul(k))
	  call potential (ndim,pos,acc,pot,dumtime)
	  barycen(j,k)=(acc(1)*cos(zcenteul(k)))
     &			 +(acc(2)*sin(zcenteul(k)))
          if (y(j).gt.sr) then
            arg=2.0*((log(y(j)-sr)/tan(spang))-zcenteul(k))+pc
            spyed(j,k)=-spamp*(y(j)-sr)**2*exp(-spsc*(y(j)-sr))*
     &                ((3.0-spsc*(y(j)-sr))*cos(arg)-
     &                (2.0/tan(spang))*sin(arg))
          else
            spyed(j,k)=0.0
          endif
c          if (y(j).gt.16.) then
c            spyed(j,k)=spyed(j,k)*exp(16.0-y(j))
c          endif
          if (radius(j).gt.sr) then
            arg=2.0*((log(radius(j)-sr)/tan(spang))-zcenteul(k))+pc
            spycen(j,k)=-spamp*(radius(j)-sr)**2*
     &                 exp(-spsc*(radius(j)-sr))*
     &                 ((3.0-spsc*(radius(j)-sr))*cos(arg)-
     &                 (2.0/tan(spang))*sin(arg))
          else
            spycen(j,k)=0.0
          endif
c bhmass increase
            bhycen(j,k)=4.29569*radius(j)
     *                  *(radius(j)**2.0+0.001**2.0)**(-1.5)
            bhyed(j,k)=4.29569*y(j)
     *                  *(y(j)**2.0+0.001**2.0)**(-1.5)
c bhmass increase end          
c          if (radius(j).gt.16.) then
c            spycen(j,k)=spycen(j,k)*exp(16.0-radius(j))
c          endif
8        continue
         if (y(j).lt.rnstart) 
     &     call taper(baryed,jn,j,ks,ke+1,y(j),rnstart,rnscale)
         if (radius(j).lt.rnstart) 
     &     call taper(barycen,jn,j,ks,ke+1,radius(j),rnstart,rnscale)

         if (y(j).gt.rtstart) 
     &     call taper(baryed,jn,j,ks,ke+1,y(j),rtstart,rtscale)
         if (radius(j).gt.rtstart) 
     &     call taper(barycen,jn,j,ks,ke+1,radius(j),rtstart,rtscale)
9       continue
      endif
c
c doosu extract information
      open(unit=12,file='forces_'//id//'.dat',status='unknown',
     &     form='unformatted')
      write(12) ((diskf(j,k),j=js-3,je+3),k=ks,ke)
      write(12) ((bulf(j,k),j=js-3,je+3),k=ks,ke)
      write(12) ((barf(j,k),j=js-3,je+3),k=ks,ke)
      write(12) ((dble(baryed(j,k)),j=js-3,je+3),k=ks,ke) 
      close(unit=12)
      open(unit=902,file='total_pot_'//id//'.dat',status='unknown',
     &     form='unformatted')
      write(902) ((pot_out(j,k),j=js-3,je+3),k=ks,ke)
      close(unit=902)
c doosu end
c doosu spiral potential
c      open (unit=900,file='spy.dat',form='unformatted')
c      write(900) spyed
c      close(unit=900)
c      open (unit=901,file='coord_th.dat',status='unknown')
c      do k=ks,ke
c         write(901,*) k,z(k)
c      enddo
c      close(unit=901)
c      stop
c doosu end
c
      if (nzz.gt.1) then
	do 11 j=js,je
	 do 10 k=ks-1,ke+1
          pos(1)=radius(j)*cos(z(k))
	  pos(2)=radius(j)*sin(z(k))
	  pos(3)=0.0
	  call potential (ndim,pos,acc,pot,dumtime)
c       wrong part!!
c	  barzed(j,k)=(acc(1)*sin(z(k)))-(acc(2)*cos(z(k)))
          barzed(j,k)=-(acc(1)*sin(z(k)))+(acc(2)*cos(z(k)))   
	  pos(1)=radius(j)*cos(zcenteul(k))
	  pos(2)=radius(j)*sin(zcenteul(k))
	  call potential (ndim,pos,acc,pot,dumtime)
c  	  barzcen(j,k)=(acc(1)*sin(zcenteul(k)))
c     &			 - (acc(2)*cos(zcenteul(k)))
          barzcen(j,k)=-(acc(1)*sin(zcenteul(k)))
     &                   + (acc(2)*cos(zcenteul(k)))   
          if (radius(j).gt.sr) then
            arg=2.0*((log(radius(j)-sr)/tan(spang))-z(k))+pc
            spzed(j,k)=-2.0*spamp*(radius(j)-sr)**2*
     &                 exp(-spsc*(radius(j)-sr))*sin(arg)
            arg=2.0*((log(radius(j)-sr)/tan(spang))-zcenteul(k))+pc
            spzcen(j,k)=-2.0*spamp*(radius(j)-sr)**2*
     &                 exp(-spsc*(radius(j)-sr))*sin(arg)
          else
            spzed(j,k)=0.0
            spzcen(j,k)=0.0
          endif
c          if (radius(j).gt.16.) then
c            spzed(j,k)=spzed(j,k)*exp(16.0-radius(j))
c            spzcen(j,k)=spzcen(j,k)*exp(16.0-radius(j))
c          endif
10     continue
       if (radius(j).lt.rnstart)
     &   call taper(barzed,jn,j,ks-1,ke+1,radius(j),rnstart,rnscale)
       if (radius(j).lt.rnstart)
     &   call taper(barzcen,jn,j,ks-1,ke+1,radius(j),rnstart,rnscale)

       if (radius(j).gt.rtstart)
     &   call taper(barzed,jn,j,ks-1,ke+1,radius(j),rtstart,rtscale)
       if (radius(j).gt.rtstart)
     &   call taper(barzcen,jn,j,ks-1,ke+1,radius(j),rtstart,rtscale)
  
11    continue
      endif
c
c      do 20 k=ks,ke
c      do 20 j=js,je


c for reatart
c      open(111,file='restart.dat')



      do 20 k=1,kn
      do 20 j=1,jn

c        read(111,*) dd,vv,ww
c Toomre Test
        d(j,k)=amp!/(radius(j)+0.5)*(1.0+0.001*ran2(idum))
c        d(j,k)=dd

       w(j,k)=sqrt(abs(radius(j)*diskycen(j,k)))
       v(j,k)=0.
c        w(j,k)=ww
c        v(j,k)=vv
20    continue
      vgz=par(1)
      vspiral=sppat-vgz


C
      if (.false.) then
         open (unit=7,file='spforces.dat')
         do 12 k=ks,ke
            do 12 j=js,je
               write (7,13) j,k,diskycen(j,k),spycen(j,k),spzcen(j,k) 
 13            format (2i5,3e14.5) 
 12      continue
         close (unit=7)
      endif
      return
      end
c
c  the potential near the center can be tapered to something more axisymmetric
c  below radius of rnstart a gaussian taper of scale length rnscale is applied
c  NOTE:: it is assumed that upon  rad < r0
c
      subroutine taper(a,ny,j,k1,k2,rad,r0,rs)
      implicit none
      integer ny,j,k1,k2
      real a(ny,1),rad,r0,rs
c
      integer k
      real t
      double precision sum1, sum2
      
      sum1 = 0.0d0
      sum2 = 0.0d0

      do 100 k=k1,k2
         sum1 = sum1 + a(j,k)
         sum2 = sum2 + a(j,k)*a(j,k)
        write(*,*) 'k=',k,a(j,k),sum1,sum2
 100  continue
      sum1 = sum1/dble(k2-k1+1)
      sum2 = sum2/dble(k2-k1+1) - sum1*sum1
      if (sum2 .lt. 0) then
         sum2 = -sqrt(-sum2)
      else
         sum2 = sqrt(sum2)
      endif

      if (rs.gt.0) then
         t = exp(- (rad-r0)**2/(2*rs**2))
         write(*,*) 'taper_1 ',sum1,sum2,j,k1,k2,rad,t
         do 101 k=k1,k2
            a(j,k) = t*a(j,k) + (1-t)*sum1
 101     continue
      else
         write(*,*) 'taper_0 ',sum1,sum2,j,k1,k2,rad,0
         do 102 k=k1,k2
            a(j,k) = sum1
 102     continue
      endif

      end

c-------------------------------------------------------------------
c Tapering inner few fixels.
c-------------------------------------------------------------------
      subroutine taper1(a,ny,nz,j,k1,k2)
      implicit none
      integer j,k,k1,k2,ny,nz
      real a(ny,nz)
      double precision sum

      sum = 0.d0

      do 200 k=k1,k2
         sum = sum + a(j,k)
  200 continue
      sum = sum/dble(k2-k1+1)
      do 201 k=k1,k2
         a(j,k) = sum
  201 continue
      return
      end


c----------------------------------------------------------------
c function of random generation [ran2: numerical recipe]
c----------------------------------------------------------------
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *1a311+)-.

