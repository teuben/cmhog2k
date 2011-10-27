

c -*- fortran -*-













c=======================================================================
c/////////////////////////////  SUBROUTINE MSTART  \\\\\\\\\\\\\\\\\\\\c
      subroutine mstart
c
c  INITIATES SIMULATION
c
c  written by: Jim Stone
c  date:       January, 1991
c  modified1:
c
c  PURPOSE:  Starts a run.
c
c  EXTERNALS: SETUP, {empty}
c
c  LOCALS:
c-----------------------------------------------------------------------
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
      integer irestart
c
      integer incr,strtoi
      external setup,  empty , strtoi
      namelist /rescon/ irestart,tdump,dtdump,id,resfile
      namelist /iocon/ thdf,dthdf,thist,dthist,tmovie,dtmovie
      namelist /mlims/ cma1,cmi1,cma2,cmi2,cma3,cmi3
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
      open(unit=1,file='cmhogin' ,status='old')
c
c     namelist /rescon/ irestart,tdump,dtdump,id,resfile
      irestart = 0
       tdump = 0.0
      dtdump = 0.0
      id     = 'bg'
      resfile = 'res000uc'
      read (1,rescon)
      write(2,rescon)
c
      open(unit=2,file='cmhogout_'//id,status='unknown')

      if (irestart .eq. 0) then
        call setup
      else
c        call mget
c        call empty
      endif
c
       thdf  = 0.0
      dthdf  = 0.0
       thist = 0.0
      dthist = 0.0
       tmovie= 0.0
      dtmovie= 0.0
      read (1,iocon)
      write(2,iocon)
      read (1,mlims)
      write(2,mlims)
c
      open (unit=3,file='history_'//id)
      open (unit=777,file='star_'//id)
      open (unit=776,file='errchk1_'//id)
      open (unit=775,file='errchk2_'//id)
      open (unit=774,file='errchk3_'//id)
      open (unit=773,file='errchk4_'//id)
c      open (unit=772,file='errchk5_'//id)
c      open (unit=771,file='errchk6_'//id)
c      open (unit=770,file='errchk7_'//id)
c      open (unit=778,file='star_'//id//'.dat',status='unknown',
c     1      form='unformatted')
      if (irestart .eq. 0) then
        write(resfile,"(a3,i4.4,a2)") 'res',0,id
        write(hdffile,"(a3,i4.4,a2)") 'hdf',0,id
        write(movfile1,"(a3,i4.4,a2)") 'mde',0,id
        write(movfile2,"(a3,i4.4,a2)") 'mvr',0,id
        write(movfile3,"(a3,i4.4,a2)") 'mvt',0,id
      else
        incr = strtoi(resfile,4,7) + 1
        write(resfile,"(a3,i4.4,a2)") 'res',incr,id
        incr = strtoi(hdffile,4,7) + 1
        write(hdffile,"(a3,i4.4,a2)") 'hdf',incr,id
        incr = strtoi(movfile1,4,7) + 1
        write(movfile1,"(a3,i4.4,a2)") 'mde',incr,id
        incr = strtoi(movfile2,4,7) + 1
        write(movfile2,"(a3,i4.4,a2)") 'mvr',incr,id
        incr = strtoi(movfile3,4,7) + 1
        write(movfile3,"(a3,i4.4,a2)") 'mvt',incr,id
      endif
c
      close(unit=1)
      close(unit=2)
      return
      end
