












c=======================================================================
c///////////////////////  SUBROUTINE YINTLGRG  \\\\\\\\\\\\\\\\\\\\\\\\c
      subroutine yintlgrg(d,v,w,dyi,dla,dra,vla,vra)
c
c  COMPUTES LEFT AND RIGHT Y-INTERFACE VALUES FOR RIEMANN SOLVER
c
c  written by: Jim Stone
c  date:       January,1991
c  modified1:
c
c  PURPOSE:  Uses piecewise parabolic interpolation to compute left-
c    and right interface values to be fed into Riemann solver during
c    Y-direction sweeps.
c
c  INPUT: d,p,v,dyi
c
c  OUTPUT: dla,dra,pla,pra,yla,yra
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
      real   d(jn,kn),v(jn,kn),w(jn,kn), dyi(jn)
      real dla(jn,kn), dra(jn,kn)
     &    ,vla(jn,kn), vra(jn,kn)

      integer j,k
      real qa,qb,qc,qd,qe,dplus,vplus,wplus
     &	   ,dmnus,vmnus,wmnus,ft,forcey,grav
      real c1(ijkn),c2(ijkn),c3(ijkn), c4(ijkn),  c5(ijkn), c6(ijkn)
     &    ,dd(ijkn),dv(ijkn),dw(ijkn),dph(ijkn),vph(ijkn),wph(ijkn)
     &    ,dl(ijkn),dr(ijkn),vl(ijkn), vr(ijkn),wl(ijkn),wr(ijkn)
     &    ,d6(ijkn),v6(ijkn), cs(ijkn),ftwd(ijkn)
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
      do 100 j=js-2,je+2
        qa    = dyi(j)/(dyi(j-1) + dyi(j) + dyi(j+1))
        c1(j) = qa*(2.0*dyi(j-1) + dyi(j))/(dyi(j+1) + dyi(j))
        c2(j) = qa*(2.0*dyi(j+1) + dyi(j))/(dyi(j-1) + dyi(j))
100   continue
      do 110 j=js-1,je+2
        qa    = dyi(j-2) + dyi(j-1) + dyi(j) + dyi(j+1)
        qb    = dyi(j-1)/(dyi(j-1) + dyi(j))
        qc    = (dyi(j-2) + dyi(j-1))/(2.0*dyi(j-1) + dyi(j  ))
        qd    = (dyi(j+1) + dyi(j  ))/(2.0*dyi(j  ) + dyi(j-1))
        qb    = qb + 2.0*dyi(j)*qb/qa*(qc-qd)
        c3(j) = 1.0 - qb
        c4(j) = qb
        c5(j) =  dyi(j  )/qa*qd
        c6(j) = -dyi(j-1)/qa*qc
110   continue
c
      do 300 k=ks,ke
c
c  Compute average linear slopes (eqn 1.7)
c
      do 210 j=js-2,je+2
        dplus = d(j+1,k)-d(j  ,k)
        vplus = v(j+1,k)-v(j  ,k)
        wplus = w(j+1,k)-w(j  ,k)
c
        dmnus = d(j  ,k)-d(j-1,k)
        vmnus = v(j  ,k)-v(j-1,k)
        wmnus = w(j  ,k)-w(j-1,k)
c
        dd(j) = c1(j)*dplus + c2(j)*dmnus
        dv(j) = c1(j)*vplus + c2(j)*vmnus
        dw(j) = c1(j)*wplus + c2(j)*wmnus
c
c  Monotonize (eqn 1.8)
c
        qa = min(abs(dd(j)),2.0*abs(dmnus),2.0*abs(dplus))
        qb = min(abs(dw(j)),2.0*abs(wmnus),2.0*abs(wplus))
        qc = min(abs(dv(j)),2.0*abs(vmnus),2.0*abs(vplus))
c
        dd(j) = cvmgt(qa, 0.0, dplus*dmnus.gt.0.0)*sign(1.0, dd(j))
        dw(j) = cvmgt(qb, 0.0, wplus*wmnus.gt.0.0)*sign(1.0, dw(j))
        dv(j) = cvmgt(qc, 0.0, vplus*vmnus.gt.0.0)*sign(1.0, dv(j))
210   continue
c
c  construct left and right values (eqn 1.6)
c
      do 220 j=js-1,je+2
       dph(j)= c3(j)*d(j-1,k)+c4(j)*d(j,k)+c5(j)*dd(j-1)+c6(j)*dd(j)
       vph(j)= c3(j)*v(j-1,k)+c4(j)*v(j,k)+c5(j)*dv(j-1)+c6(j)*dv(j)
       wph(j)= c3(j)*w(j-1,k)+c4(j)*w(j,k)+c5(j)*dw(j-1)+c6(j)*dw(j)
220   continue
c
      do 230 j=js-1,je+1
        dl(j) = dph(j  )
        vl(j) = vph(j  )
        wl(j) = wph(j  )
        dr(j) = dph(j+1)
        vr(j) = vph(j+1)
        wr(j) = wph(j+1)
230   continue
c
c  Monotonize again (eqn 1.10)
c
      do 240 j=js-1,je+1
        qa = (dr(j)-d(j,k))*(d(j,k)-dl(j))
        qd = dr(j)-dl(j)
        qe = 6.0*(d(j,k)-0.5*(dr(j)+dl(j)))
        dl(j) = cvmgt(d(j,k), dl(j), qa.le.0.0)
        dr(j) = cvmgt(d(j,k), dr(j), qa.le.0.0)
        dl(j) = cvmgt(3.0*d(j,k)-2.0*dr(j), dl(j), qd**2-qd*qe.lt.0.0)
        dr(j) = cvmgt(3.0*d(j,k)-2.0*dl(j), dr(j), qd**2+qd*qe.lt.0.0)
c
        qa = (vr(j)-v(j,k))*(v(j,k)-vl(j))
        qd = vr(j)-vl(j)
        qe = 6.0*(v(j,k)-0.5*(vr(j)+vl(j)))
        vl(j) = cvmgt(v(j,k), vl(j), qa.le.0.0)
        vr(j) = cvmgt(v(j,k), vr(j), qa.le.0.0)
        vl(j) = cvmgt(3.0*v(j,k)-2.0*vr(j), vl(j), qd**2-qd*qe.lt.0.0)
        vr(j) = cvmgt(3.0*v(j,k)-2.0*vl(j), vr(j), qd**2+qd*qe.lt.0.0)
c
        qa = (wr(j)-w(j,k))*(w(j,k)-wl(j))
        qd = wr(j)-wl(j)
        qe = 6.0*(w(j,k)-0.5*(wr(j)+wl(j)))
        wl(j) = cvmgt(w(j,k), wl(j), qa.le.0.0)
        wr(j) = cvmgt(w(j,k), wr(j), qa.le.0.0)
        wl(j) = cvmgt(3.0*w(j,k)-2.0*wr(j), wl(j), qd**2-qd*qe.lt.0.0)
        wr(j) = cvmgt(3.0*w(j,k)-2.0*wl(j), wr(j), qd**2+qd*qe.lt.0.0)
240   continue
c
c  Construct flattening parameter 
c
      if (ifltn .ne. 0) then
      do 250 j=js-1,je+1
	qa=(radius(j+1)*v(j+1,k)-radius(j-1)*v(j-1,k))/radius(j)
        qc=cvmgt((0.5*(dyi(j-1)+dyi(j+1))+dyi(j))
     &	  /(radius(j)*(0.5*(dz(k-1)+dz(k+1))+dz(k)))
     &	  *(w(j,k+1)-w(j,k-1)),0.0,nzz.gt.1)
	qd=qa+qc
	qe=cvmgt(1.0,0.0,((qd**2/ciso**2).gt.0.082).and.(qd.lt.0.))
	ftwd(j)=qe*max(0.0,min(0.5,10.0*((d(j+1,k)-d(j-1,k))/
     &		(d(j+2,k)-d(j-2,k)+tiny)-0.75)))
250   continue
      ftwd(js-2) = ftwd(js-1)
      ftwd(je+2) = ftwd(je+1)
c
c  Flatten (eqn 4.1)
c
      do 260 j=js-1,je+1
        qa = cvmgt(max(ftwd(j),ftwd(j+1))
     &            ,max(ftwd(j),ftwd(j-1)), d(j+1,k)-d(j-1,k).lt.0.0)
        dl(j) = d(j,k)*qa + dl(j)*(1.0-qa)
        dr(j) = d(j,k)*qa + dr(j)*(1.0-qa)
        vl(j) = v(j,k)*qa + vl(j)*(1.0-qa)
        vr(j) = v(j,k)*qa + vr(j)*(1.0-qa)
        wl(j) = w(j,k)*qa + wl(j)*(1.0-qa)
        wr(j) = w(j,k)*qa + wr(j)*(1.0-qa)
260   continue
      endif
c
c  Now construct left and right interface values (eqn 1.12)
c
      do 270 j=js-1,je+1
        d6(j) = 6.0*(d(j,k)-0.5*(dl(j)+dr(j)))
        v6(j) = 6.0*(v(j,k)-0.5*(vl(j)+vr(j)))
        dd(j) = dr(j) - dl(j)
        dv(j) = vr(j) - vl(j)
        cs(j) = (dt*ciso)/(2.0*dyi(j))
270   continue
      do 280 j=js,je+1
        dla(j,k)= dr(j-1)-cs(j-1)*(dd(j-1)-(1.0-ft*cs(j-1))*d6(j-1))
        dra(j,k)= dl(j  )+cs(j  )*(dd(j  )+(1.0-ft*cs(j  ))*d6(j  ))
c
	grav=(1-barfract)*diskyed(j,k)+barfract*baryed(j,k)
     &       +spifract*spyed(j,k)


	forcey=grav+(wl(j)**2/y(j))
        vla(j,k)= vr(j-1)-cs(j-1)*(dv(j-1)-(1.0-ft*cs(j-1))*v6(j-1))
	vla(j,k)= vla(j,k)+(0.5*dt*forcey)
        vra(j,k)= vl(j  )+cs(j  )*(dv(j  )+(1.0-ft*cs(j  ))*v6(j  ))
	vra(j,k)= vra(j,k)+(0.5*dt*forcey)
280   continue
c
300   continue
      return
      end
