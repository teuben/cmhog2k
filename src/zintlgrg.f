












c=======================================================================
c///////////////////////  SUBROUTINE ZINTLGRG  \\\\\\\\\\\\\\\\\\\\\\\\c
      subroutine zintlgrg(d,v,w,dzi,dla,dra,wla,wra)
c
c  COMPUTES LEFT AND RIGHT Z-INTERFACE VALUES FOR RIEMANN SOLVER
c
c  written by: Jim Stone
c  date:       January,1991
c  modified1:
c
c  PURPOSE:  Uses piecewise parabolic interpolation to compute left-
c    and right interface values to be fed into Riemann solver during
c    Z-direction sweeps.
c
c  INPUT: d,p,w,dzi
c
c  OUTPUT: dla,dra,pla,pra,wla,wra
c
c  EXTERNALS: CVMGT
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
      real   d(jn,kn),v(jn,kn), w(jn,kn), dzi(kn)
      real dla(jn,kn), dra(jn,kn)
     &    ,wla(jn,kn), wra(jn,kn)

      integer j,k
      real qa,qb,qc,qd,qe,dplus,wplus,dmnus,wmnus
     &	   ,vplus,vmnus,ft,forcez,grav
      real c1(ijkn),c2(ijkn),c3(ijkn), c4(ijkn),  c5(ijkn), c6(ijkn)
     &    ,dd(ijkn),dw(ijkn),dph(ijkn),wph(ijkn),dv(ijkn),vph(ijkn)
     &    ,dl(ijkn),dr(ijkn),wl(ijkn), wr(ijkn),vl(ijkn),vr(ijkn)
     &    ,d6(ijkn),w6(ijkn), cs(ijkn),ftwd(ijkn)
      parameter(ft = 4.0/3.0)
c
c	on some machines this arithmetic statement function
c	will work, on strict ansi (e.g. the f2c compiler) you
c	will need the cvmgt.src function, on the cray CVMGT is
c	an inline vector instruction
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////
c=======================================================================
c  Compute coefficients used in interpolation formulae
c
      do 100 k=ks-2,ke+2
        qa    = dzi(k)/(dzi(k-1) + dzi(k) + dzi(k+1))
        c1(k) = qa*(2.0*dzi(k-1) + dzi(k))/(dzi(k+1) + dzi(k))
        c2(k) = qa*(2.0*dzi(k+1) + dzi(k))/(dzi(k-1) + dzi(k))
100   continue
      do 110 k=ks-1,ke+2
        qa    = dzi(k-2) + dzi(k-1) + dzi(k) + dzi(k+1)
        qb    = dzi(k-1)/(dzi(k-1) + dzi(k))
        qc    = (dzi(k-2) + dzi(k-1))/(2.0*dzi(k-1) + dzi(k  ))
        qd    = (dzi(k+1) + dzi(k  ))/(2.0*dzi(k  ) + dzi(k-1))
        qb    = qb + 2.0*dzi(k)*qb/qa*(qc-qd)
        c3(k) = 1.0 - qb
        c4(k) = qb
        c5(k) =  dzi(k  )/qa*qd
        c6(k) = -dzi(k-1)/qa*qc
110   continue
c
      do 300 j=js,je
c
c  Compute average linear slopes (eqn 1.7)
c
      do 210 k=ks-2,ke+2
        dplus = d(j,k+1)-d(j,k  )
        wplus = w(j,k+1)-w(j,k  )
        vplus = v(j,k+1)-v(j,k  )
c
        dmnus = d(j,k  )-d(j,k-1)
        wmnus = w(j,k  )-w(j,k-1)
        vmnus = v(j,k  )-v(j,k-1)
c
        dd(k) = c1(k)*dplus + c2(k)*dmnus
        dw(k) = c1(k)*wplus + c2(k)*wmnus
        dv(k) = c1(k)*vplus + c2(k)*vmnus
c
c  Monotonize (eqn 1.8)
c
        qa = min(abs(dd(k)),2.0*abs(dmnus),2.0*abs(dplus))
        qb = min(abs(dv(k)),2.0*abs(vmnus),2.0*abs(vplus))
        qc = min(abs(dw(k)),2.0*abs(wmnus),2.0*abs(wplus))
c
        dd(k) = cvmgt(qa, 0.0, dplus*dmnus.gt.0.0)*sign(1.0, dd(k))
        dv(k) = cvmgt(qb, 0.0, vplus*vmnus.gt.0.0)*sign(1.0, dv(k))
        dw(k) = cvmgt(qc, 0.0, wplus*wmnus.gt.0.0)*sign(1.0, dw(k))
210   continue
c
c  construct left and right values (eqn 1.6)
c
      do 220 k=ks-1,ke+2
       dph(k)= c3(k)*d(j,k-1)+c4(k)*d(j,k)+c5(k)*dd(k-1)+c6(k)*dd(k)
       wph(k)= c3(k)*w(j,k-1)+c4(k)*w(j,k)+c5(k)*dw(k-1)+c6(k)*dw(k)
       vph(k)= c3(k)*v(j,k-1)+c4(k)*v(j,k)+c5(k)*dv(k-1)+c6(k)*dv(k)
220   continue
c
      do 230 k=ks-1,ke+1
        dl(k) = dph(k  )
        wl(k) = wph(k  )
        vl(k) = vph(k  )
        dr(k) = dph(k+1)
        wr(k) = wph(k+1)
        vr(k) = vph(k+1)
230   continue
c
c  Monotonize again (eqn 1.10)
c
      do 240 k=ks-1,ke+1
        qa = (dr(k)-d(j,k))*(d(j,k)-dl(k))
        qd = dr(k)-dl(k)
        qe = 6.0*(d(j,k)-0.5*(dr(k)+dl(k)))
        dl(k) = cvmgt(d(j,k), dl(k), qa.le.0.0)
        dr(k) = cvmgt(d(j,k), dr(k), qa.le.0.0)
        dl(k) = cvmgt(3.0*d(j,k)-2.0*dr(k), dl(k), qd**2-qd*qe.lt.0.0)
        dr(k) = cvmgt(3.0*d(j,k)-2.0*dl(k), dr(k), qd**2+qd*qe.lt.0.0)
c
        qa = (wr(k)-w(j,k))*(w(j,k)-wl(k))
        qd = wr(k)-wl(k)
        qe = 6.0*(w(j,k)-0.5*(wr(k)+wl(k)))
        wl(k) = cvmgt(w(j,k), wl(k), qa.le.0.0)
        wr(k) = cvmgt(w(j,k), wr(k), qa.le.0.0)
        wl(k) = cvmgt(3.0*w(j,k)-2.0*wr(k), wl(k), qd**2-qd*qe.lt.0.0)
        wr(k) = cvmgt(3.0*w(j,k)-2.0*wl(k), wr(k), qd**2+qd*qe.lt.0.0)
c
        qa = (vr(k)-v(j,k))*(v(j,k)-vl(k))
        qd = vr(k)-vl(k)
        qe = 6.0*(v(j,k)-0.5*(vr(k)+vl(k)))
        vl(k) = cvmgt(v(j,k), vl(k), qa.le.0.0)
        vr(k) = cvmgt(v(j,k), vr(k), qa.le.0.0)
        vl(k) = cvmgt(3.0*v(j,k)-2.0*vr(k), vl(k), qd**2-qd*qe.lt.0.0)
        vr(k) = cvmgt(3.0*v(j,k)-2.0*vl(k), vr(k), qd**2+qd*qe.lt.0.0)
240   continue
c
c  Construct flattening parameter
c
      if (ifltn .ne. 0) then
      do 250 k=ks-1,ke+1
	qa=w(j,k+1)-w(j,k-1)
	qc=cvmgt((radius(j)*(0.5*(dzi(k-1)+dzi(k+1))+dzi(k)))
     &    /(0.5*(dy(j-1)+dy(j+1))+dy(j))
     &    *(radius(j+1)*v(j+1,k)-radius(j-1)*v(j-1,k))/radius(j),
     &    0.0,nyz.gt.1)
        qd=qa+qc
        qe=cvmgt(1.0,0.0,((qd**2/ciso**2).gt.0.082).and.(qd.lt.0.))
        ftwd(k)=qe*max(0.0,min(0.5,10.0*((d(j,k+1)-d(j,k-1))/
     &          (d(j,k+2)-d(j,k-2)+tiny)-0.75)))
250   continue
      ftwd(ks-2) = ftwd(ks-1)
      ftwd(ke+2) = ftwd(ke+1)
c
c  Flatten (eqn 4.1)
c
      do 260 k=ks-1,ke+1
        qa = cvmgt(max(ftwd(k),ftwd(k+1))
     &            ,max(ftwd(k),ftwd(k-1)), d(j,k+1)-d(j,k-1).lt.0.0)
        dl(k) = d(j,k)*qa + dl(k)*(1.0-qa)
        dr(k) = d(j,k)*qa + dr(k)*(1.0-qa)
        wl(k) = w(j,k)*qa + wl(k)*(1.0-qa)
        wr(k) = w(j,k)*qa + wr(k)*(1.0-qa)
        vl(k) = v(j,k)*qa + vl(k)*(1.0-qa)
        vr(k) = v(j,k)*qa + vr(k)*(1.0-qa)
260   continue
      endif
c
c  Now construct left and right interface values (eqn 1.12)
c
      do 270 k=ks-1,ke+1
        d6(k) = 6.0*(d(j,k)-0.5*(dl(k)+dr(k)))
        w6(k) = 6.0*(w(j,k)-0.5*(wl(k)+wr(k)))
        dd(k) = dr(k) - dl(k)
        dw(k) = wr(k) - wl(k)
        cs(k) = (dt*ciso)/(2.0*dzi(k)*radius(j))
270   continue
      do 280 k=ks,ke+1
        dla(j,k)= dr(k-1)-cs(k-1)*(dd(k-1)-(1.0-ft*cs(k-1))*d6(k-1))
        dra(j,k)= dl(k  )+cs(k  )*(dd(k  )+(1.0-ft*cs(k  ))*d6(k  ))
c
	grav=barfract*barzed(j,k)+spifract*spzed(j,k)


	forcez=grav-(wl(k)*vl(k)/radius(j))
        wla(j,k)= wr(k-1)-cs(k-1)*(dw(k-1)-(1.0-ft*cs(k-1))*w6(k-1))
	wla(j,k)= wla(j,k)+(0.5*dt*forcez)
        wra(j,k)= wl(k  )+cs(k  )*(dw(k  )+(1.0-ft*cs(k  ))*w6(k  ))
	wra(j,k)= wra(j,k)+(0.5*dt*forcez)
280   continue
c
300   continue
      return
      end
