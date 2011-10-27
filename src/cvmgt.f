

c
c only used for machines which cannot handle conversion
c from logical or integer (0,1); see "cvmgt.h"
c










































      real function cvmgt (qaa,qbb,ilog)
      real qaa,qbb
      logical ilog

      if (ilog) then
         cvmgt = qaa
      else
         cvmgt = qbb
      endif

      return
      end

