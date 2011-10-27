












c=======================================================================
c/////////////////////////  SUBROUTINE YINTRMP  \\\\\\\\\\\\\\\\\\\\\\\c
      subroutine yintrmp(d,v,w,dyi,df,vf,wf)
c
c  COMPUTES INTERFACE FLUXES IN Y-DIRECTION FOR REMAP
c
c  written by: Jim Stone
c  date:       January, 1991
c  modified1:
c
c  PURPOSE:  Interpolates the fundamental variables on the Lagrangean
c    mesh for remap.  The interpolated parabolae are then integrated
c    to compute the "effective flux" for the remap at each interface.
c
c  INPUT ARGUMENTS: d=density; p=pressure; e=total energy; u,v,w=x,y,z
c    components of velocity.  These are all zone centered (defined over
c    nxz,nyz,nzz).
c
c  OUTPUT ARGUMENTS: df=density flux; ef=total energy flux; uf,vf,wf=
c    x,y,z-components of velocity fluxes.  These are all centered at
c    interfaces in y-direction (defined over nxz,nyz+1,nzz).
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
c-----------------------------------------------------------------------
c  BOUNDARY VARIABLES
c  niib = integer flag for every zone at Inner I Boundary
c  diib = constant value of d at IIB (used for niib=3); same for e,u,v,w
c    same for OIB, IJB, OJB, IKB, OKB
c
      integer  nijb(kn), nikb(jn)
     &        ,nojb(kn), nokb(jn)
      common /bndryi/ nijb, nikb, nojb, nokb
c
      real
     &   dijb(kn), vijb(kn), wijb(kn)
     &  ,dojb(kn), vojb(kn), wojb(kn)
     &  ,dikb(jn), vikb(jn), wikb(jn)
     &  ,dokb(jn), vokb(jn), wokb(jn)
c
      common /bndryr/  dijb, vijb, wijb
     &                ,dojb, vojb, wojb
     &                ,dikb, vikb, wikb
     &                ,dokb, vokb, wokb
      real  d(jn,kn)
     &    , v(jn,kn), w(jn,kn), dyi(jn,kn)
      real df(jn,kn),vf(jn,kn)
     &    ,wf(jn,kn)

      integer j,k
      real qa, qb, qc, qd, qe, s1, s2, udtim1, udti, ft
     & ,dplus, vplus, wplus 
     & ,dmnus, vmnus, wmnus, yy
      real  c1(ijkn), c2(ijkn), c3(ijkn), c4(ijkn), c5(ijkn), c6(ijkn)
     &    , dd(ijkn), dv(ijkn), dw(ijkn), f(ijkn)
     &    ,dph(ijkn),vph(ijkn),wph(ijkn)
     &    ,d2d(ijkn),dyb(ijkn), dl(ijkn), vl(ijkn)
     &    , wl(ijkn), dr(ijkn), vr(ijkn), wr(ijkn)
     &    , d6(ijkn), v6(ijkn), w6(ijkn), vcorr(ijkn)
      real delta
c      real cs
      parameter(ft=4.0/3.0)
c
c	on some machines this arithmetic statement function
c	will work, on strict ansi (e.g. the f2c compiler) you
c	will need the cvmgt.src function, on the cray CVMGT is
c	an inline vector instruction
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////
c=======================================================================
c
      do 200 k=ks,ke
c  compute volume element correction factor
      do 90 j=js-3,je+3
	  vcorr(j)=rvollag(j,k)/dyi(j,k)
90    continue
c
c  Compute coefficients used in interpolation formulae
c
      do 100 j=js-2,je+2
        qa    = dyi(j,k)/(dyi(j-1,k) + dyi(j,k) + dyi(j+1,k))
        c1(j)=qa*(2.0*dyi(j-1,k)+dyi(j,k))/(dyi(j+1,k)+dyi(j,k))
        c2(j)=qa*(2.0*dyi(j+1,k)+dyi(j,k))/(dyi(j-1,k)+dyi(j,k))
100   continue
      do 110 j=js-1,je+2
        qa = dyi(j-2,k) + dyi(j-1,k) + dyi(j,k) + dyi(j+1,k)
        qb = dyi(j-1,k)/(dyi(j-1,k) + dyi(j,k))
        qc = (dyi(j-2,k)+dyi(j-1,k))/(2.0*dyi(j-1,k)+dyi(j  ,k))
        qd = (dyi(j+1,k)+dyi(j  ,k))/(2.0*dyi(j  ,k)+dyi(j-1,k))
        qb = qb + 2.0*dyi(j,k)*qb/qa*(qc-qd)
        c3(j) = 1.0 - qb 
        c4(j) = qb
        c5(j) =  dyi(j  ,k)/qa*qd
        c6(j) = -dyi(j-1,k)/qa*qc
110   continue
c
c  Compute average linear slopes (eqn 1.7)
c
      do 210 j=js-2,je+2
        dplus = (d(j+1,k)*vcorr(j+1))-(d(j  ,k)*vcorr(j))
        vplus = (v(j+1,k)*vcorr(j+1))-(v(j  ,k)*vcorr(j))
        wplus = (w(j+1,k)*vcorr(j+1))-(w(j  ,k)*vcorr(j))
c
        dmnus = (d(j  ,k)*vcorr(j))-(d(j-1,k)*vcorr(j-1))
        vmnus = (v(j  ,k)*vcorr(j))-(v(j-1,k)*vcorr(j-1))
        wmnus = (w(j  ,k)*vcorr(j))-(w(j-1,k)*vcorr(j-1))
c
        dd(j) = c1(j)*dplus + c2(j)*dmnus
        dv(j) = c1(j)*vplus + c2(j)*vmnus
        dw(j) = c1(j)*wplus + c2(j)*wmnus
210   continue
c
c  construct interface values (eqn 1.6)
c
      do 220 j=js-1,je+2
       dph(j)= (c3(j)*d(j-1,k)*vcorr(j-1))+(c4(j)*d(j,k)*vcorr(j))
     &	+(c5(j)*dd(j-1))+(c6(j)*dd(j))
       vph(j)= (c3(j)*v(j-1,k)*vcorr(j-1))+(c4(j)*v(j,k)*vcorr(j))
     &	+(c5(j)*dv(j-1))+(c6(j)*dv(j))
       wph(j)= (c3(j)*w(j-1,k)*vcorr(j-1))+(c4(j)*w(j,k)*vcorr(j))
     &	+(c5(j)*dw(j-1))+(c6(j)*dw(j))
220   continue
c
c  left and right values
c
      do 250 j=js-1,je+1
        dl(j) = dph(j  )/abs(yn(j,k))
        vl(j) = vph(j  )/abs(yn(j,k))
        wl(j) = wph(j  )/abs(yn(j,k))
        dr(j) = dph(j+1)/abs(yn(j+1,k))
        vr(j) = vph(j+1)/abs(yn(j+1,k))
        wr(j) = wph(j+1)/abs(yn(j+1,k))
c
c  monotonize
c
      dl(j)=cvmgt(max(d(j,k),d(j-1,k)),dl(j),
     &     (dl(j).gt.d(j,k)).and.(dl(j).gt.d(j-1,k)))
      vl(j)=cvmgt(max(v(j,k),v(j-1,k)),vl(j),
     &     (vl(j).gt.v(j,k)).and.(vl(j).gt.v(j-1,k)))
      wl(j)=cvmgt(max(w(j,k),w(j-1,k)),wl(j),
     &     (wl(j).gt.w(j,k)).and.(wl(j).gt.w(j-1,k)))
      dl(j)=cvmgt(min(d(j,k),d(j-1,k)),dl(j),
     &     (dl(j).lt.d(j,k)).and.(dl(j).lt.d(j-1,k)))
      vl(j)=cvmgt(min(v(j,k),v(j-1,k)),vl(j),
     &     (vl(j).lt.v(j,k)).and.(vl(j).lt.v(j-1,k)))
      wl(j)=cvmgt(min(w(j,k),w(j-1,k)),wl(j),
     &     (wl(j).lt.w(j,k)).and.(wl(j).lt.w(j-1,k)))
      dr(j)=cvmgt(max(d(j,k),d(j+1,k)),dr(j),
     &     (dr(j).gt.d(j,k)).and.(dr(j).gt.d(j+1,k)))
      vr(j)=cvmgt(max(v(j,k),v(j+1,k)),vr(j),
     &     (vr(j).gt.v(j,k)).and.(vr(j).gt.v(j+1,k)))
      wr(j)=cvmgt(max(w(j,k),w(j+1,k)),wr(j),
     &     (wr(j).gt.w(j,k)).and.(wr(j).gt.w(j+1,k)))
      dr(j)=cvmgt(min(d(j,k),d(j+1,k)),dr(j),
     &     (dr(j).lt.d(j,k)).and.(dr(j).lt.d(j+1,k)))
      vr(j)=cvmgt(min(v(j,k),v(j+1,k)),vr(j),
     &     (vr(j).lt.v(j,k)).and.(vr(j).lt.v(j+1,k)))
      wr(j)=cvmgt(min(w(j,k),w(j+1,k)),wr(j),
     &     (wr(j).lt.w(j,k)).and.(wr(j).lt.w(j+1,k)))
250   continue
c
c  Monotonize again (eqn 1.10)
c
      do 280 j=js-1,je+1
	qa = dyi(j,k)/yn(j,k)
	f(j)=qa/(6.+3.*qa)
        qa = (dr(j)-d(j,k))*(d(j,k)-dl(j))
        qd = dr(j)-dl(j)
        qe = 6.0*(d(j,k)-0.5*(dr(j)+dl(j)))
        dl(j) = cvmgt(d(j,k), dl(j), qa.le.0.0)
        dr(j) = cvmgt(d(j,k), dr(j), qa.le.0.0)
        dl(j) = cvmgt((6.0*d(j,k)-((4+3*f(j))*dr(j)))
     &	  /(2-3*f(j)),dl(j),qd**2-qd*qe.lt.0.0)
        dr(j) = cvmgt((6.0*d(j,k)-((4-3*f(j))*dl(j)))
     &	  /(2+3*f(j)),dr(j),qd**2+qd*qe.lt.0.0)
c
        qa = (vr(j)-v(j,k))*(v(j,k)-vl(j))
        qd = vr(j)-vl(j)
        qe = 6.0*(v(j,k)-0.5*(vr(j)+vl(j)))
        vl(j) = cvmgt(v(j,k), vl(j), qa.le.0.0)
        vr(j) = cvmgt(v(j,k), vr(j), qa.le.0.0)
        vl(j) = cvmgt((6.0*v(j,k)-((4+3*f(j))*vr(j)))
     &	  /(2-3*f(j)),vl(j),qd**2-qd*qe.lt.0.0)
        vr(j) = cvmgt((6.0*v(j,k)-((4-3*f(j))*vl(j)))
     &	  /(2+3*f(j)),vr(j),qd**2+qd*qe.lt.0.0)
c
        qa = (wr(j)-w(j,k))*(w(j,k)-wl(j))
        qd = wr(j)-wl(j)
        qe = 6.0*(w(j,k)-0.5*(wr(j)+wl(j)))
        wl(j) = cvmgt(w(j,k), wl(j), qa.le.0.0)
        wr(j) = cvmgt(w(j,k), wr(j), qa.le.0.0)
        wl(j) = cvmgt((6.0*w(j,k)-((4+3*f(j))*wr(j)))
     &	  /(2-3*f(j)),wl(j),qd**2-qd*qe.lt.0.0)
        wr(j) = cvmgt((6.0*w(j,k)-((4-3*f(j))*wl(j)))
     &	  /(2+3*f(j)),wr(j),qd**2+qd*qe.lt.0.0)
280   continue
c
c  Now construct interface fluxes (eqn 1.12)
c
      do 290 j=js-1,je+1
        d6(j) = 6.0*(d(j,k)-0.5*((dl(j)*(1-f(j)))+(dr(j)*(1+f(j)))))
        v6(j) = 6.0*(v(j,k)-0.5*((vl(j)*(1-f(j)))+(vr(j)*(1+f(j)))))
        w6(j) = 6.0*(w(j,k)-0.5*((wl(j)*(1-f(j)))+(wr(j)*(1+f(j)))))
290   continue
      do 295 j=js,je+1
	delta = abs(0.5 * ((yn(j,k)**2)-(y(j)**2)))
        delta=cvmgt(-delta,delta,yn(j,k).lt.y(j))
        yy     = (yn(j,k)-dt*vgy) - y(j)
        udtim1 = abs(yy/(2.0*dyi(j-1,k)))
        udti   = abs(yy/(2.0*dyi(j  ,k)))
c
        qa = dr(j-1)-udtim1*(dr(j-1)-dl(j-1)-(1.0-ft*udtim1)*d6(j-1))
        qb = dl(j  )+udti  *(dr(j  )-dl(j  )+(1.0-ft*udti  )*d6(j  ))
        df(j,k) = delta*cvmgt(qa, qb, yy .ge. 0.0)
c
        yy = df(j,k)
        udtim1 = abs(yy/(2.0*d(j-1,k)*rvollag(j-1,k)))
        udti   = abs(yy/(2.0*d(j  ,k)*rvollag(j,k)))
c
        qc = vr(j-1)-udtim1*(vr(j-1)-vl(j-1)-(1.0-ft*udtim1)*v6(j-1))
        qd = vl(j  )+udti  *(vr(j  )-vl(j  )+(1.0-ft*udti  )*v6(j  ))
        vf(j,k) = yy*cvmgt(qc, qd, yy .ge. 0.0)
c
        qc = wr(j-1)-udtim1*(wr(j-1)-wl(j-1)-(1.0-ft*udtim1)*w6(j-1))
        qd = wl(j  )+udti  *(wr(j  )-wl(j  )+(1.0-ft*udti  )*w6(j  ))
        wf(j,k) = yy*cvmgt(qc, qd, yy .ge. 0.0)
295   continue
c
        massflux(k)=df(js,k)
        omassflux(k)=df(je+1,k)

      qa = (z(ke+1)-z(ks))

200   continue

      return
      end
