












c=======================================================================
c///////////////////////////  SUBROUTINE INTCHK  \\\\\\\\\\\\\\\\\\\\\\c
      subroutine intchk(iswres,iswhdf,iswhst)
c
c  CHECKS FOR INTERRUPTS OR STOPPING CRITERION
c
c  written by: Jim Stone
c  date:       January, 1991
c
c  PURPOSE:  Reads the buffer for valid interrupt messages and takes
c    appropriate action. Also checks stopping criteria, and sets ifsen=1
c    if stop condition is detected.
c
c  INPUT ARGUMENTS: none
c
c  OUTPUT ARGUMENTS: iswres,iswhdf,iswhst=switches for restart, hdf, and
c    history dumps; set to 1 if dump to be made in DATAIO
c
c  EXTERNALS:  CHECKIN, BCDFLT, FINDNO, [ETIME,SECOND]
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
      integer iswres,iswhdf,iswhst
c
      integer i,nchar,istrt,iend
      real valnew
      character*80 msg
c
      integer checkin
      external checkin,bcdflt,findno
      real tarray(2),etime
      external etime
c  List of valid interrupt messages
      character*3 intmsg(13)
      data intmsg /  'sto','?  ','pau','abo','tli','nli','tsa','dum'
     &  ,'dtd','hdf','dtf','hst','dth' /
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////
c=======================================================================
c  Check stopping criteria
c
      tused = etime(tarray)
      trem   = ttotal - tused
      if (tlim .gt. 0.0 .and. time .ge. tlim) then
        write(6,2000) tlim,nlim,tsave,time,nhy,trem,tused
2000    format(/1x,'terminating on time limit',
     .         /1x,'tlim=',1pe12.5,'   nlim=',i7,'  tsave=',1e12.5,
     .         /1x,'time=',1e12.5,'  cycle=',i7,'   trem=',1e12.5,
     .           '  tused=',1e12.5)
        ifsen = 1
      endif
c      if (nlim .gt. 0    .and. nhy .ge. nlim) then
c        write(6 ,2010) tlim,nlim,tsave,time,nhy,trem,tused
c2010    format(/1x,'terminating on cycle limit',
c     .         /1x,'tlim=',1pe12.5,'   nlim=',i7,'  tsave=',1e12.5,
c     .         /1x,'time=',1e12.5,'  cycle=',i7,'   trem=',1e12.5,
c     .           '  tused=',1e12.5)
c        ifsen = 1
c      endif
c      if (tsave .gt. 0.0 .and. trem .le. tsave) then
c        write(6,2020) tlim,nlim,tsave,time,nhy,trem,tused
c2020    format(/1x,'terminating on reserve time limit',
c     .         /1x,'tlim=',1pe12.5,'   nlim=',i7,'  tsave=',1e12.5,
c     .         /1x,'time=',1e12.5,'  cycle=',i7,'   trem=',1e12.5,
c     .           '  tused=',1e12.5)
c        ifsen = 1
c      endif
      if (ifsen .eq. 1) return

c  Check for interrupt messages.  If none or illegal message found
c  then return
c
      nchar = checkin(msg,1)
      if (nchar .eq. 0) return
      do i=1,13
        if (msg(1:3) .eq. intmsg(i)) goto 20
      enddo
      write(6,2030) msg(1:3)
2030  format(1x,a3,' is not an interrupt message.  Legal messages are:'
     . ,/1x,'sto ? pau abo tli nli tsa dum dtd hdf dtf hst dth')
      return
20    continue
c
c  Legal interrupt message found, process it
c
c  stop command
c
      if (msg(1:3) .eq. 'sto') then
        write(6,2040) msg,tlim,nlim,tsave,time,nhy,trem,tused
2040    format(1x,a3,': execution stopped with',
     .        /1x,'tlim=',1pe12.5,'   nlim=',i7,'  tsave=',1e12.5,
     .        /1x,'time=',1e12.5,'  cycle=',i7,'   trem=',1e12.5,
     .          '  tused=',1e12.5)
        ifsen = 1
        return
      endif
c
c  status command
c
      if (msg(1:3) .eq. '?  ') then
        write(6,2050) msg,tlim,nlim,tsave,time,nhy,trem,tused
2050    format(1x,a3,': execution continuing with',
     .        /1x,'tlim=',1pe12.5,'   nlim=',i7,'  tsave=',1e12.5,
     .        /1x,'time=',1e12.5,'  cycle=',i7,'   trem=',1e12.5,
     .          '  tused=',1e12.5)
        return
      endif
c
c  pause command
c
      if (msg(1:3) .eq. 'pau') then
        write(6,2060) msg,tlim,nlim,tsave,time,nhy,trem,tused
2060    format(1x,a3,': execution halted with',
     &       /1x,'tlim=',1pe12.5,'   nlim=',i7,'  tsave=',1e12.5,
     &       /1x,'time=',1e12.5,'  cycle=',i7,'   trem=',1e12.5,
     &         '  tused=',1e12.5,/1x,'Hit any key to restart execution')
        nchar = checkin(msg,0)
        return
      endif
c
c  abort command
c
      if (msg(1:3) .eq. 'abo') then
        write(6,"(a3,': ABORT! do you want to abort execution?')") msg
        write(6,"(' (type yes or no)')")
        nchar = checkin(msg,0)
        if (nchar .eq. 0) return
        if (msg(1:3) .ne. 'yes') then
          write(6,"('Abort cancelled, continuing execution ...')")
          return
        else
          write(6,"('ABORT , CRASH , BOOM!!!')")
          stop
        endif
      endif
c
c  reset physical time limit (tlim) command
c
      if (msg(1:3) .eq. 'tli') then
        call findno(msg,istrt,iend)
        if (istrt .lt. 0) then
          write(6,2130) msg
          return
        endif
        call bcdflt(msg,valnew,(istrt-1),(iend-istrt+1))
        if (valnew .lt. 0.0 .or. valnew .ge. huge) then
          write(6,2130) msg
2130      format(1x,a3,': could not read reset number; execution ',
     .                 'continuing')
          return
        endif
        tlim = valnew
        write(6,"(a3,': tlim reset to ',1pe12.5)") msg,tlim
        return
      endif
c
c  reset cycle limit (nlim) command
c
      if (msg(1:3) .eq. 'nli') then
        call findno(msg,istrt,iend)
        if (istrt .lt. 0) then
          write(6,2130) msg
          return
        endif
        call bcdflt(msg,valnew,(istrt-1),(iend-istrt+1))
        if (valnew .lt. 0.0 .or. valnew .ge. huge) then
          write(6,2130) msg
          return
        endif
        nlim = nint(valnew)
        write(6,"(a3,': nlim reset to ',i12)") msg,nlim
      endif
c
c  reset reserve time (tsave) command
c
      if (msg(1:3) .eq. 'tsa') then
        call findno(msg,istrt,iend)
        if (istrt .lt. 0) then
          write(6,2130) msg
          return
        endif
        call bcdflt(msg,valnew,(istrt-1),(iend-istrt+1))
        if (valnew .lt. 0.0 .or. valnew .ge. huge) then
          write(6,2130) msg
          return
        endif
        tsave = valnew
        write(6,"(a3,': tsave reset to ',1pe12.5)") msg,tsave
      endif
c
c  turn restart dump switch on
c
      if (msg(1:3) .eq. 'dum') then
        write(6,"(a3,': restart dump switch on')") msg
        iswres = 1
        return
      endif
c
c  reset dump frequency (dtdump) command
c
      if (msg(1:3) .eq. 'dtd') then
        call findno(msg,istrt,iend)
        if (istrt .lt. 0) then
          write(6,2130) msg
          return
        endif
        call bcdflt(msg,valnew,(istrt-1),(iend-istrt+1))
        if (valnew .lt. 0.0 .or. valnew .ge. huge) then
          write(6,2130) msg
          return
        endif
        dtdump = valnew
        write(6,"(a3,': dtdump reset to ',1pe12.5)") msg,dtdump
      endif
c
c  turn hdf dumps on
c
      if (msg(1:3) .eq. 'hdf') then
        write(6,"(a3,': hdf dump switch on')") msg
        iswhdf = 1
        return
      endif
c
c  reset hdf dump frequency (dthdf) command
c
      if (msg(1:3) .eq. 'dtf') then
        call findno(msg,istrt,iend)
        if (istrt .lt. 0) then
          write(6,2130) msg
          return
        endif
        call bcdflt(msg,valnew,(istrt-1),(iend-istrt+1))
        if (valnew .lt. 0.0 .or. valnew .ge. huge) then
          write(6,2130) msg
          return
        endif
        dthdf = valnew
        write(6,"(a3,': dthdf reset to ',1pe12.5)") msg,dthdf
      endif
c
c  turn history dumps on
c
      if (msg(1:3) .eq. 'hst') then
        write(6,"(a3,': hdf dump switch on')") msg
        iswhst = 1
        return
      endif
c
c  reset hdf dump frequency (dthdf) command
c
      if (msg(1:3) .eq. 'dth') then
        call findno(msg,istrt,iend)
        if (istrt .lt. 0) then  
          write(6,2130) msg  
          return 
        endif  
        call bcdflt(msg,valnew,(istrt-1),(iend-istrt+1))
        if (valnew .lt. 0.0 .or. valnew .ge. huge) then
          write(6,2130) msg
          return 
        endif  
        dthist = valnew   
        write(6,"(a3,': dthist reset to ',1pe12.5)") msg,dthdf
      endif
c
      return
      end
