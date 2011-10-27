












c=======================================================================
c////////////////////////  SUBROUTINE DATAIO  \\\\\\\\\\\\\\\\\\\\\\\\\c
      subroutine dataio(iswres,iswhdf,iswhst,iswmov)
c
c  CONTROLS DATA OUTPUT
c
c  written by: Jim Stone
c  date        January, 1991
c  modified1:
c
c  PURPOSE:  Controls data I/O for restart, HDF, and history dumps.
c
c  INPUT ARGUMENTS: iswres,iswhdf,iswhst=switches for restart,hdf, and
c    history dumps.  Values of 1 ensure dumps will be made.
c
c  OUTPUT ARGUMENTS: none
c
c  EXTERNALS: MSAVE, HDFALL, STRTOI, HDUMP
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
c-----------------------------------------------------------------------
      integer iswres,iswhdf,iswhst,iswmov,i
c
      integer incr,strtoi
      external hdfall, strtoi,  hdump, movie
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
      if (dtmovie.le.0.0) iswmov = 0
      if (dtdump .gt. 0.0 .and. tused .ge. (tdump+dtdump)) then
        tdump = tdump + dtdump
        iswres = 1
      endif
      if (dthdf  .gt. 0.0 .and. time  .ge. (thdf+dthdf)) then
        thdf = thdf + dthdf
        iswhdf = 1
      endif
      if (dthist .gt. 0.0 .and. time  .ge. (thist+dthist)) then
        thist = thist + dthist
        iswhst = 1
      endif
      if (dtmovie .gt. 0.0 .and. time  .ge. (tmovie+dtmovie)) then
        tmovie = tmovie + dtmovie
        iswmov = 1
      endif
c
c  restart dump
c
c      if (iswres .eq. 1) then
c        call msave
c        incr = strtoi(resfile,4,7) + 1
c        write(resfile,"(a3,i4.4,a2)") 'res',incr,id
c        iswres=0
c      endif
c
c  HDF dump
c
      if (iswhdf .eq. 1) then
        call hdfall(hdffile)

c doosu writing time when putting hdf
c retime is real time (yr)
        retime = 1.d3*3.0857d18/(1.d5)/(365.d0*24.d0*60.d0*60.d0)
     &           *time
        write(51,*) hdffile,time,retime
c doosu end

        incr = strtoi(hdffile,4,7) + 1
        write(hdffile,"(a3,i4.4,a2)") 'hdf',incr,id
        write(strfile,"(a3,i4.4,a2)") 'str',incr,id
        write(*,"('Wrote HDF dump ',i4.4,' at time:',1pg15.8,1pg15.8)")
     1       incr-1,time,dt

        iswhdf=0
      endif
c
      if (iswmov .eq. 1) then
        call movie(movfile1,movfile2,movfile3)
        incr = strtoi(movfile1,4,7) + 1
        write(movfile1,"(a3,i4.4,a2)") 'mde',incr,id
        incr = strtoi(movfile2,4,7) + 1
        write(movfile2,"(a3,i4.4,a2)") 'mvr',incr,id
        incr = strtoi(movfile3,4,7) + 1
        write(movfile3,"(a3,i4.4,a2)") 'mvt',incr,id
        iswmov=0
      endif
c
c  history dump
c
      if (iswhst .eq. 1) then
        call hdump
        iswhst=0
      endif
      return
      end
