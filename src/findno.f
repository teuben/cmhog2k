












c=======================================================================
c///////////////////////  SUBROUTINE FINDNO  \\\\\\\\\\\\\\\\\\\\\\\\\\c
      subroutine findno(msg,istrt,iend)
c
c  FINDS START AND END OF RESET NUMBER IN INTERRUPT MESSAGE STRING
c
c     written by: Jim Stone
c     date:       June, 1988
c     modified1:
c
c  INPUT ARGUMENTS:
c     msg   = character*80 string containing the message
c
c  OUTPUT ARGUMENTS:
c     istrt = number of first character of reset number in message
c             If no number is found, a value of -1 is returned
c     iend  = number of last character of reset number in message
c
c  LOCALS:
c-----------------------------------------------------------------------
      character*80 msg
      integer istrt,iend,iblnk
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
      istrt=0
      iend =0
c
      do 10 iblnk=1,80
        if (msg(iblnk:iblnk) .eq. ' ') goto 20
10    continue
20    continue
      do 30 istrt=iblnk,80
        if (msg(istrt:istrt) .ne. ' ') goto 100
30    continue
      istrt = -1
      return
c
100   continue
      do 110 iend=istrt,80
        if (msg(iend:iend) .eq. ' ') goto 120
110   continue
120   iend=iend-1
c
      return
      end
