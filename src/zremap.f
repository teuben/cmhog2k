












c=======================================================================
c///////////////////////////  SUBROUTINE ZREMAP  \\\\\\\\\\\\\\\\\\\\\\c
      subroutine zremap(d,v,w,df,vf,wf)
c
c  PERFORMS REMAP IN Z-DIRECTION
c
c  written by: Jim Stone
c  date:       January, 1991
c
c  PURPOSE:  Remaps the fundamental variables d,e,u,v,w back to the
c    Eulerian mesh after Lagrangean update.  Fluxes for each variable
c    computed in XINTRMP are used.
c
c  INPUT ARGUMENTS: d=density; e=total energy; u,v,w=x,y,z components of
c    velocity.  All are zone centered (defined over nxz,nyz,nzz).
c    df=density flux; ef=total energy flux; uf,vf,wf=x,y,z-components of
c    velocity fluxes.  These are all centered at interfaces in
c    z-direction (defined over nxz,nyz,nzz+1).
c
c  OUTPUT ARGUMENTS: Updated values for d,e,u,v,w.
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
      real d(jn,kn),v(jn,kn),w(jn,kn)
      real df(jn,kn),vf(jn,kn)
     &    ,wf(jn,kn),dm,dmnu
c
      integer j,k
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////
c=======================================================================
c
      do 10 k=ks,ke
      do 10 j=js,je
        dm = d(j,k)*dzn(j,k)*radius(j)
        dmnu = dm - (df(j,k+1)-df(j,k))
        d(j,k) = dmnu/(dz(k)*radius(j))
        v(j,k) = (dm*v(j,k) - (vf(j,k+1)-vf(j,k)))/dmnu
        w(j,k) = (dm*w(j,k) - (wf(j,k+1)-wf(j,k)))/dmnu
10    continue
c
      return
      end
