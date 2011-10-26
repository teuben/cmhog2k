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
