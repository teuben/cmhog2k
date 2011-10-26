      real psdymin,psdymax,psdyrat,psddymin,psdzmin,psdzmax,psddzmin
      integer npsd
      parameter (npsd=5000)
      real psdy(npsd),psdz(npsd)
      integer psdyind(npsd),psdzind(npsd)
      common /pseudo/ psdymin,psdymax,psdyrat,psddymin,psdzmin
     &               ,psdzmax,psddzmin,psdy,psdz,psdyind,psdzind
 

