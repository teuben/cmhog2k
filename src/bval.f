












c=======================================================================
c//////////////  SUBROUTINES BVAL* (multiple routines)  \\\\\\\\\\\\\\\c
c  BOUNDARY VALUE ROUTINES FOR ALL DEPENDENT VARIABLES ON EACH BOUNDARY
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////////
c=======================================================================
c-------------------------  density boundary values  -------------------
c
      subroutine bvaldyeul
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
      integer k
c  include volumes because radial zones will not have 
c  equal volumes- conserve mass
c
c  Inner J Boundary
c
      do 30 k=ks,ke
        if (nijb(k) .eq. 1) d(js-1,k)=d(js,k)
        if (nijb(k) .eq. 2) d(js-1,k)=d(js,k)
        if (nijb(k) .eq. 3) d(js-1,k)=dijb(k)
        if (nijb(k) .eq. 4) d(js-1,k)=d(je,k)
c
        if (nijb(k) .eq. 1) d(js-2,k)=d(js+1,k)
        if (nijb(k) .eq. 2) d(js-2,k)=d(js,k)
        if (nijb(k) .eq. 3) d(js-2,k)=dijb(k)
        if (nijb(k) .eq. 4) d(js-2,k)=d(je-1,k)
c
        if (nijb(k) .eq. 1) d(js-3,k)=d(js+2,k)
        if (nijb(k) .eq. 2) d(js-3,k)=d(js,k)
        if (nijb(k) .eq. 3) d(js-3,k)=dijb(k)
        if (nijb(k) .eq. 4) d(js-3,k)=d(je-2,k)
30    continue
c
c  Outer J Boundary
c
      do 40 k=ks,ke
        if (nojb(k) .eq. 1) d(je+1,k)=d(je,k)
        if (nojb(k) .eq. 2) d(je+1,k)=d(je,k)
        if (nojb(k) .eq. 3) d(je+1,k)=dojb(k)
        if (nojb(k) .eq. 4) d(je+1,k)=d(js,k)
c
        if (nojb(k) .eq. 1) d(je+2,k)=d(je-1,k)
        if (nojb(k) .eq. 2) d(je+2,k)=d(je,k)
        if (nojb(k) .eq. 3) d(je+2,k)=dojb(k)
        if (nojb(k) .eq. 4) d(je+2,k)=d(js+1,k)
c
        if (nojb(k) .eq. 1) d(je+3,k)=d(je-2,k)
        if (nojb(k) .eq. 2) d(je+3,k)=d(je,k)
        if (nojb(k) .eq. 3) d(je+3,k)=dojb(k)
        if (nojb(k) .eq. 4) d(je+3,k)=d(js+2,k)
40    continue
      return
      end
c
      subroutine bvaldylag
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
      integer k
c  include volumes because radial zones will not have 
c  equal volumes- conserve mass
c
c  Inner J Boundary
c
      do 30 k=ks,ke
        if (nijb(k) .eq. 1) d(js-1,k)=d(js,k)
     &	    *rvollag(js,k)/rvollag(js-1,k)
        if (nijb(k) .eq. 2) d(js-1,k)=d(js,k)
     &	    *rvollag(js,k)/rvollag(js-1,k)
        if (nijb(k) .eq. 3) d(js-1,k)=dijb(k)
     &	    *rvollag(js,k)/rvollag(js-1,k)
        if (nijb(k) .eq. 4) d(js-1,k)=d(je,k)
     &	    *rvollag(je,k)/rvollag(js-1,k)
c
        if (nijb(k) .eq. 1) d(js-2,k)=d(js+1,k)
     &	    *rvollag(js+1,k)/rvollag(js-2,k)
        if (nijb(k) .eq. 2) d(js-2,k)=d(js,k)
     &	    *rvollag(js,k)/rvollag(js-2,k)
        if (nijb(k) .eq. 3) d(js-2,k)=dijb(k)
     &	    *rvollag(js,k)/rvollag(js-2,k)
        if (nijb(k) .eq. 4) d(js-2,k)=d(je-1,k)
     &	    *rvollag(je-1,k)/rvollag(js-2,k)
c
        if (nijb(k) .eq. 1) d(js-3,k)=d(js+2,k)
     &	    *rvollag(js+2,k)/rvollag(js-3,k)
        if (nijb(k) .eq. 2) d(js-3,k)=d(js,k)
     &	    *rvollag(js,k)/rvollag(js-3,k)
        if (nijb(k) .eq. 3) d(js-3,k)=dijb(k)
     &	    *rvollag(js,k)/rvollag(js-3,k)
        if (nijb(k) .eq. 4) d(js-3,k)=d(je-2,k)
     &	    *rvollag(je-2,k)/rvollag(js-3,k)
30    continue
c
c  Outer J Boundary
c
      do 40 k=ks,ke
        if (nojb(k) .eq. 1) d(je+1,k)=d(je,k)
     &	    *rvollag(je,k)/rvollag(je+1,k)
        if (nojb(k) .eq. 2) d(je+1,k)=d(je,k)
     &	    *rvollag(je,k)/rvollag(je+1,k)
        if (nojb(k) .eq. 3) d(je+1,k)=dojb(k)
     &	    *rvollag(je,k)/rvollag(je+1,k)
        if (nojb(k) .eq. 4) d(je+1,k)=d(js,k)
     &	    *rvollag(js,k)/rvollag(je+1,k)
c
        if (nojb(k) .eq. 1) d(je+2,k)=d(je-1,k)
     &	    *rvollag(je-1,k)/rvollag(je+2,k)
        if (nojb(k) .eq. 2) d(je+2,k)=d(je,k)
     &	    *rvollag(je,k)/rvollag(je+2,k)
        if (nojb(k) .eq. 3) d(je+2,k)=dojb(k)
     &	    *rvollag(je,k)/rvollag(je+2,k)
        if (nojb(k) .eq. 4) d(je+2,k)=d(js+1,k)
     &	    *rvollag(js+1,k)/rvollag(je+2,k)
c
        if (nojb(k) .eq. 1) d(je+3,k)=d(je-2,k)
     &	    *rvollag(je-2,k)/rvollag(je+3,k)
        if (nojb(k) .eq. 2) d(je+3,k)=d(je,k)
     &	    *rvollag(je,k)/rvollag(je+3,k)
        if (nojb(k) .eq. 3) d(je+3,k)=dojb(k)
     &	    *rvollag(je,k)/rvollag(je+3,k)
        if (nojb(k) .eq. 4) d(je+3,k)=d(js+2,k)
     &	    *rvollag(js+2,k)/rvollag(je+3,k)
40    continue
      return
      end
c
      subroutine bvaldz
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
      integer j
c
c  Inner K Boundary
c
      do 50 j=js,je
        if (nikb(j) .eq. 1) d(j,ks-1) = d   (j,ks  )
        if (nikb(j) .eq. 2) d(j,ks-1) = d   (j,ks  )
        if (nikb(j) .eq. 3) d(j,ks-1) = dikb(j     )
        if (nikb(j) .eq. 4) d(j,ks-1) = d   (j,ke  )
c
        if (nikb(j) .eq. 1) d(j,ks-2) = d   (j,ks+1)
        if (nikb(j) .eq. 2) d(j,ks-2) = d   (j,ks  )
        if (nikb(j) .eq. 3) d(j,ks-2) = dikb(j     )
        if (nikb(j) .eq. 4) d(j,ks-2) = d   (j,ke-1)
c
        if (nikb(j) .eq. 1) d(j,ks-3) = d   (j,ks+2)
        if (nikb(j) .eq. 2) d(j,ks-3) = d   (j,ks  )
        if (nikb(j) .eq. 3) d(j,ks-3) = dikb(j     )
        if (nikb(j) .eq. 4) d(j,ks-3) = d   (j,ke-2)
50    continue
c
c  Outer K Boundary
c
      do 60 j=js,je
        if (nokb(j) .eq. 1) d(j,ke+1) = d   (j,ke  )
        if (nokb(j) .eq. 2) d(j,ke+1) = d   (j,ke  )
        if (nokb(j) .eq. 3) d(j,ke+1) = dokb(j     )
        if (nokb(j) .eq. 4) d(j,ke+1) = d   (j,ks  )
c
        if (nokb(j) .eq. 1) d(j,ke+2) = d   (j,ke-1)
        if (nokb(j) .eq. 2) d(j,ke+2) = d   (j,ke  )
        if (nokb(j) .eq. 3) d(j,ke+2) = dokb(j     )
        if (nokb(j) .eq. 4) d(j,ke+2) = d   (j,ks+1)
c
        if (nokb(j) .eq. 1) d(j,ke+3) = d   (j,ke-2)
        if (nokb(j) .eq. 2) d(j,ke+3) = d   (j,ke  )
        if (nokb(j) .eq. 3) d(j,ke+3) = dokb(j     )
        if (nokb(j) .eq. 4) d(j,ke+3) = d   (j,ks+2)
60    continue
      return
      end
c---------------------  y-velocity boundary values  --------------------
c
      subroutine bvalvy
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
      integer k
      real q1
c
c  Inner J Boundary
c
      do 30 k=ks,ke
        if (nijb(k) .eq. 1) v(js-1,k) =-v   (js  ,k)
        if (nijb(k) .eq. 2) then
           q1 = sign(0.5, v(js+1,k))
           v(js,k)   = v(js+1,k) * (0.5 - q1)
           v(js-1,k) = v(js,k)
        endif
        if (nijb(k) .eq. 3) v(js-1,k) = vijb(     k)
        if (nijb(k) .eq. 4) v(js-1,k) = v   (je  ,k)
c
        if (nijb(k) .eq. 1) v(js-2,k) =-v   (js+1,k)
        if (nijb(k) .eq. 2) v(js-2,k) = v   (js  ,k)
        if (nijb(k) .eq. 3) v(js-2,k) = vijb(     k)
        if (nijb(k) .eq. 4) v(js-2,k) = v   (je-1,k)
c
        if (nijb(k) .eq. 1) v(js-3,k) =-v   (js+2,k)
        if (nijb(k) .eq. 2) v(js-3,k) = v   (js  ,k)
        if (nijb(k) .eq. 3) v(js-3,k) = vijb(     k)
        if (nijb(k) .eq. 4) v(js-3,k) = v   (je-2,k)
30    continue
c
c  Outer J Boundary
c
      do 40 k=ks,ke
        if (nojb(k) .eq. 1) v(je+1,k) =-v   (je  ,k)
        if (nojb(k) .eq. 2) v(je+1,k) = v   (je  ,k)
        if (nojb(k) .eq. 3) v(je+1,k) = vojb(     k)
        if (nojb(k) .eq. 4) v(je+1,k) = v   (js  ,k)
c
        if (nojb(k) .eq. 1) v(je+2,k) =-v   (je-1,k)
        if (nojb(k) .eq. 2) v(je+2,k) = v   (je  ,k)
        if (nojb(k) .eq. 3) v(je+2,k) = vojb(     k)
        if (nojb(k) .eq. 4) v(je+2,k) = v   (js+1,k)
c
        if (nojb(k) .eq. 1) v(je+3,k) =-v   (je-2,k)
        if (nojb(k) .eq. 2) v(je+3,k) = v   (je  ,k)
        if (nojb(k) .eq. 3) v(je+3,k) = vojb(     k)
        if (nojb(k) .eq. 4) v(je+3,k) = v   (js+2,k)
40    continue
      return
      end
c
      subroutine bvalvz
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
      integer j
c
c  Inner K Boundary
c
      do 50 j=js,je
        if (nikb(j) .eq. 1) v(j,ks-1) = v   (j,ks  )
        if (nikb(j) .eq. 2) v(j,ks-1) = v   (j,ks  )
        if (nikb(j) .eq. 3) v(j,ks-1) = vikb(j     )
        if (nikb(j) .eq. 4) v(j,ks-1) = v   (j,ke  )
c
        if (nikb(j) .eq. 1) v(j,ks-2) = v   (j,ks+1)
        if (nikb(j) .eq. 2) v(j,ks-2) = v   (j,ks  )
        if (nikb(j) .eq. 3) v(j,ks-2) = vikb(j     )
        if (nikb(j) .eq. 4) v(j,ks-2) = v   (j,ke-1)
c
        if (nikb(j) .eq. 1) v(j,ks-3) = v   (j,ks+2)
        if (nikb(j) .eq. 2) v(j,ks-3) = v   (j,ks  )
        if (nikb(j) .eq. 3) v(j,ks-3) = vikb(j     )
        if (nikb(j) .eq. 4) v(j,ks-3) = v   (j,ke-2)
50    continue
c
c  Outer K Boundary
c
      do 60 j=js,je
        if (nokb(j) .eq. 1) v(j,ke+1) = v   (j,ke  )
        if (nokb(j) .eq. 2) v(j,ke+1) = v   (j,ke  )
        if (nokb(j) .eq. 3) v(j,ke+1) = vokb(j     )
        if (nokb(j) .eq. 4) v(j,ke+1) = v   (j,ks  )
c
        if (nokb(j) .eq. 1) v(j,ke+2) = v   (j,ke-1)
        if (nokb(j) .eq. 2) v(j,ke+2) = v   (j,ke  )
        if (nokb(j) .eq. 3) v(j,ke+2) = vokb(j     )
        if (nokb(j) .eq. 4) v(j,ke+2) = v   (j,ks+1)
c
        if (nokb(j) .eq. 1) v(j,ke+3) = v   (j,ke-2)
        if (nokb(j) .eq. 2) v(j,ke+3) = v   (j,ke  )
        if (nokb(j) .eq. 3) v(j,ke+3) = vokb(j     )
        if (nokb(j) .eq. 4) v(j,ke+3) = v   (j,ks+2)
60    continue
      return
      end
c---------------------  z-velocity boundary values  --------------------
c
      subroutine bvalwy
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
      integer k
c
c  Inner J Boundary
c
      do 30 k=ks,ke
        if (nijb(k) .eq. 1) w(js-1,k) = w   (js  ,k)
        if (nijb(k) .eq. 2) w(js-1,k) = w   (js  ,k)
        if (nijb(k) .eq. 3) w(js-1,k) = wijb(     k)
        if (nijb(k) .eq. 4) w(js-1,k) = w   (je  ,k)
c
        if (nijb(k) .eq. 1) w(js-2,k) = w   (js+1,k)
        if (nijb(k) .eq. 2) w(js-2,k) = w   (js  ,k)
        if (nijb(k) .eq. 3) w(js-2,k) = wijb(     k)
        if (nijb(k) .eq. 4) w(js-2,k) = w   (je-1,k)
c
        if (nijb(k) .eq. 1) w(js-3,k) = w   (js+2,k)
        if (nijb(k) .eq. 2) w(js-3,k) = w   (js  ,k)
        if (nijb(k) .eq. 3) w(js-3,k) = wijb(     k)
        if (nijb(k) .eq. 4) w(js-3,k) = w   (je-2,k)
30    continue
c
c  Outer J Boundary
c
      do 40 k=ks,ke
        if (nojb(k) .eq. 1) w(je+1,k) = w   (je  ,k)
        if (nojb(k) .eq. 2) w(je+1,k) = w   (je  ,k)
        if (nojb(k) .eq. 3) w(je+1,k) = wojb(     k)
        if (nojb(k) .eq. 4) w(je+1,k) = w   (js  ,k)
c
        if (nojb(k) .eq. 1) w(je+2,k) = w   (je-1,k)
        if (nojb(k) .eq. 2) w(je+2,k) = w   (je  ,k)
        if (nojb(k) .eq. 3) w(je+2,k) = wojb(     k)
        if (nojb(k) .eq. 4) w(je+2,k) = w   (js+1,k)
c
        if (nojb(k) .eq. 1) w(je+3,k) = w   (je-2,k)
        if (nojb(k) .eq. 2) w(je+3,k) = w   (je  ,k)
        if (nojb(k) .eq. 3) w(je+3,k) = wojb(     k)
        if (nojb(k) .eq. 4) w(je+3,k) = w   (js+2,k)
40    continue
      return
      end
c
      subroutine bvalwz
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
      integer j
c
c  Inner K Boundary
c
      do 50 j=js,je
        if (nikb(j) .eq. 1) w(j,ks-1) =-w   (j,ks  )
        if (nikb(j) .eq. 2) w(j,ks-1) = w   (j,ks  )
        if (nikb(j) .eq. 3) w(j,ks-1) = wikb(j     )
        if (nikb(j) .eq. 4) w(j,ks-1) = w   (j,ke  )
c
        if (nikb(j) .eq. 1) w(j,ks-2) =-w   (j,ks+1)
        if (nikb(j) .eq. 2) w(j,ks-2) = w   (j,ks  )
        if (nikb(j) .eq. 3) w(j,ks-2) = wikb(j     )
        if (nikb(j) .eq. 4) w(j,ks-2) = w   (j,ke-1)
c
        if (nikb(j) .eq. 1) w(j,ks-3) =-w   (j,ks+2)
        if (nikb(j) .eq. 2) w(j,ks-3) = w   (j,ks  )
        if (nikb(j) .eq. 3) w(j,ks-3) = wikb(j     )
        if (nikb(j) .eq. 4) w(j,ks-3) = w   (j,ke-2)
50    continue
c
c  Outer K Boundary
c
      do 60 j=js,je
        if (nokb(j) .eq. 1) w(j,ke+1) =-w   (j,ke  )
        if (nokb(j) .eq. 2) w(j,ke+1) = w   (j,ke  )
        if (nokb(j) .eq. 3) w(j,ke+1) = wokb(j     )
        if (nokb(j) .eq. 4) w(j,ke+1) = w   (j,ks  )
c
        if (nokb(j) .eq. 1) w(j,ke+2) =-w   (j,ke-1)
        if (nokb(j) .eq. 2) w(j,ke+2) = w   (j,ke  )
        if (nokb(j) .eq. 3) w(j,ke+2) = wokb(j     )
        if (nokb(j) .eq. 4) w(j,ke+2) = w   (j,ks+1)
c
        if (nokb(j) .eq. 1) w(j,ke+3) =-w   (j,ke-2)
        if (nokb(j) .eq. 2) w(j,ke+3) = w   (j,ke  )
        if (nokb(j) .eq. 3) w(j,ke+3) = wokb(j     )
        if (nokb(j) .eq. 4) w(j,ke+3) = w   (j,ks+2)
60    continue
      return
      end
