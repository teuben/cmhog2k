












c=======================================================================
c///////////////////////  SUBROUTINE ISOSHOCK  \\\\\\\\\\\\\\\\\\\\\\\\c
      subroutine isoshockquad(j1,j2,k1,k2,dl,dr,ul,ur,ps,us)
c
c  RIEMANN SOLVER FOR ISOTHERMAL SHOCKS
c
c  written by: Jim Stone
c  date:       July, 1993
c  modified1:
c
c  LOCALS
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
c  SCRATCH SPACE
c
      real  wijk0(jn,kn), wijk1(jn,kn), wijk2(jn,kn)
     &     ,wijk3(jn,kn), wijk4(jn,kn), wijk5(jn,kn)
c
      common /scr1/ wijk0, wijk1, wijk2, wijk3, wijk4, wijk5
c-----------------------------------------------------------------------
      integer j1,j2,k1,k2
      real  dl(jn,kn), ul(jn,kn), ps(jn,kn)
     &     ,dr(jn,kn), ur(jn,kn), us(jn,kn)
c-----------------------------------------------------------------------
      integer j,k
      real qa,qb,qc,rs
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////
c=======================================================================
c
      do 300 k=k1,k2
      do 300 j=j1,j2
c
c     compute quadratic equation parameters
	qa=ciso*(1./sqrt(dr(j,k))+1./sqrt(dl(j,k)))
	qb=ur(j,k)-ul(j,k)
	qc=-ciso*(sqrt(dr(j,k))+sqrt(dl(j,k)))
        rs = ((-qb+sqrt(qb**2-4.*qa*qc))/(2.*qa))**2
        ps(j,k) = ciso**2*rs
        us(j,k) = ur(j,k)+ciso*(rs-dr(j,k))
     &		    /sqrt(rs*dr(j,k))
c
300   continue
      return
      end

