












c=======================================================================
c/////////////////////////  SUBROUTINE MOVIE  \\\\\\\\\\\\\\\\\\\\\\\\c
      subroutine movie(filename1,filename2,filename3)
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
      character*9  filename,filename1,filename2,filename3
      character*1  image(icartx,icarty)
c
c
      integer i,j,k,jx,kx,jpn,jmn,kpn,kmn,ival
      integer irecl,mcount
      real dcon,cartxr,cartxl,cartytop,cartybot,delx,dely
      real cartx(icartx),carty(icarty),cartrad,cartthe,scrtch
      real c1,c2,c3,c4,a1,a2
      real varint(jn,kn),grav
      external bvaldz,bvalvz,bvalwz
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
      irecl=icartx*icarty
      call bvaldz
      call bvalvz
      call bvalwz
      do 100 mcount=1,3
        if (mcount.eq.1) then
          do 11 j=js,je 
          do 11 k=ks,ke+1
cpjt            varint(j,k)=d(j,k)
            varint(j,k)=d(j,k)-10.0
11        continue
          filename=filename1
          cmax=cma1
          cmin=cmi1 
          cmax=-0.1
          cmin=0.1 
        endif
c
        if (mcount.eq.2) then
          do 12 j=js,je 
          do 12 k=ks,ke+1
            varint(j,k)=v(j,k)
12        continue
          filename=filename2
          cmax=cma2
          cmin=cmi2 
        endif
c
        if (mcount.eq.3) then
          do 13 j=js,je 
          do 13 k=ks,ke+1
            grav=(1-barfract)*diskycen(j,k)+barfract*barycen(j,k)
     &           +spifract*spycen(j,k) 
            varint(j,k)=w(j,k)-sqrt(abs(radius(j)*grav))
13        continue
          filename=filename3
          cmax=cma3
          cmin=cmi3 
        endif
c
      dcon=253./(cmax-cmin)
c
      cartxr=y(je+1)
      cartxl=-y(je+1)
      cartytop=cartxr
      cartybot=cartxl
c
      delx=(cartxr-cartxl)/icartx
      dely=(cartytop-cartybot)/icarty
      do 20 i=1,icartx
	cartx(i)=cartxl+real(i-1)*delx
20    continue
      do 30 j=1,icarty
	carty(j)=cartybot+real(j-1)*dely
30    continue
c
      do 40 i=(icartx/2)+1,icartx
      do 40 j=1,icarty
	cartrad=sqrt(cartx(i)**2+carty(j)**2)
	cartthe=atan(carty(j)/(cartx(i)+tiny))
	if ((cartrad.gt.y(je)).or.(cartrad.lt.y(js))) then
	  scrtch=0.0
        else
          do 50 jx=js+1,je
	    if (y(jx).ge.cartrad) goto 60
50        continue
60        jpn=jx
	  jmn=jx-1
	  do 70 kx=ks+1,ke+1
	    if (z(kx).ge.cartthe) goto 80
70        continue
80        kpn=kx
	  kmn=kx-1
	  c1=(cartrad-y(jpn))/(y(jmn)-y(jpn))
	  c2=(cartrad-y(jmn))/(y(jpn)-y(jmn))
	  c3=(cartthe-z(kpn))/(z(kmn)-z(kpn))
	  c4=(cartthe-z(kmn))/(z(kpn)-z(kmn))
	  a1=varint(jmn,kmn)*c1+varint(jpn,kmn)*c2
	  a2=varint(jmn,kpn)*c1+varint(jpn,kpn)*c2
          if (mcount.eq.1) then
cpjt	    scrtch=(log10(a1*c3+a2*c4)-cmin)*dcon
            scrtch=((a1*c3+a2*c4)-cmin)*dcon
          else
	    scrtch=((a1*c3+a2*c4)-cmin)*dcon
          endif
        endif
	ival=int(scrtch)
        image(i,icarty+2-j)=char(min(max(1,ival),253))
        image(icartx+2-i,j)=char(min(max(1,ival),253))
40    continue
      open(unit=15,file=filename,status='unknown',access='direct',
     &  form='unformatted',recl=irecl)
      write(15,rec=1) ((image(i,j),i=1,icartx),j=1,icarty)
      close(unit=15)
100   continue
c
      return
      end
