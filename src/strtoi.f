












c=======================================================================
c////////////////////////  FUNCTION STRTOI  \\\\\\\\\\\\\\\\\\\\\\\\\\\c
      integer function strtoi(string,istrt,iend)
c
c  CONVERTS A SEGMENT OF A CHARACTER*8 STRING INTO AN INTEGER NUMBER
c  grumpf, why not call this str8toi then....
c
c-----------------------------------------------------------------------
c
      character string*(*)
      integer istrt,iend,asciic(48:57),ishift,ival,i
      data asciic / 0,1,2,3,4,5,6,7,8,9 /
      save asciic
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
      ishift = 1
      ival   = 0
      do 10 i=iend,istrt,-1
        ival = ival + ishift*asciic(ichar(string(i:i)))
        ishift = ishift*10
10    continue
      strtoi = ival
c
      return
      end
