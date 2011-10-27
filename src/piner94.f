

c -*- fortran -*-













C
C   History:
C
C   Taken from Lia:
C       Date: Thu, 21 Oct 1993 08:22:15 EDT
C       From: lia@OBMARA.CNRS-MRS.FR
C       Subject: RE: Bar Potentials
C
C   21-oct-93   Put into potential(5NEMO) format                PJT
C		Special case for (0,0,0) -- doesnt work because of LOG()
C    4-mar-94   Adding the K=2 bar, renaming to piner94.f       PJT
C               (renamed FORBA to FORBA1 - added FORBA2)
C               The potentials of the bars have been messed up a bit
C               to correct for bars being along the Y axis in some
C               simulations. Better was to have this handled in
C               the caller routines (e.g. potential)
C   12-apr-01   better warning when QM exceeds the max for given AOB PJT
C   21-sep-02   internal grid patch (set nx=ny=0 and no grid used)
C    4-oct-02   conversion to two pot0,pot1 grid for axisym + 'bar'   PJT
C   31-jul-03   added rtstart,rtscale                                 PJT
C   11-aug-03   taper at r0 and r1 now                                PJT
C   28-sep-03   bomb when interpolations are going off the grid       PJT
C
C   ToDO:
C	code DACSH()
C
C   NOTE:   it''s not possible to use single quotes in .src files !!!
C 	    (due to out usage of the cpp preprocessor)
C=======================================================================
c  warning: can only have one 'BLOCK DATA' per program
      BLOCK DATA
      IMPLICIT NONE
      INTEGER nx,ny
      REAL gpot0,gpot1,crpix1,crpix2,cdelt1,cdelt2,crval1,crval2,grmax
      REAL gomegp,rhalo,vhalo,gamma,rtstart,rtscale
      COMMON/ggrid/nx,ny,gpot0(1030,1030),gpot1(1030,1030),
     *             grmax,gomegp,
     *             crpix1,crpix2,cdelt1,cdelt2,crval1,crval2,
     *             rhalo,vhalo,gamma,rtstart,rtscale
      DATA nx/0/
      DATA ny/0/
      END
C=======================================================================
      SUBROUTINE inipotential(npar, par, name)
C
      IMPLICIT NONE
      INTEGER npar
      REAL*8 par(*)
      CHARACTER name*(*)
C-----------------------------------------------------------------------
C INITIALIZE THE POTENTIAL AND FORCE CALCULATIONS (after Lias INPOT)
C
      INTEGER nx,ny
      REAL gpot0,gpot1,crpix1,crpix2,cdelt1,cdelt2,crval1,crval2,grmax
      REAL gomegp,rhalo,vhalo,gamma,rtstart,rtscale
      COMMON/ggrid/nx,ny,gpot0(1030,1030),gpot1(1030,1030),
     *             grmax,gomegp,
     *             crpix1,crpix2,cdelt1,cdelt2,crval1,crval2,
     *             rhalo,vhalo,gamma,rtstart,rtscale
      COMMON/ AREAD/FCMASS,AXIRAT,QUAD,CENDEN,RADLAG, INDEX
      REAL*8 FCMASS, AXIRAT, QUAD, CENDEN, RADLAG
      INTEGER INDEX
      REAL*8 BARSTR
C 
      COMMON/AXISY/C1BUL,CONST, CORAD,VT2,VTDR2, RTI2,Mbh,EPSI,rb,gamm,
     &             BULMAS
      REAL*8 C1BUL, CONST, CORAD, VT2, VTDR2, RTI2,Mbh,EPSI,rb,gamm
C 
      COMMON /COM1/ OMEGS, A2, B2, C2, E, NEQ
      INTEGER NEQ
      REAL*8 OMEGS, A2, B2, C2, E
C
      COMMON /COM2/ H, EPS, IOUT, IORB, ISIGN, ICHOOSE
      INTEGER IOUT, IORB, ISIGN, ICHOOSE
      REAL*8 H, EPS
C
      COMMON/PROLA/ELIPM,A,B,C,AAIN,CCIN,AAI,CCO,ACONST
      REAL*8 ELIPM, A, B, C, ACONST, CCO, CCIN, AAI, AAIN
C
      COMMON/PROLA2/COEFF,COEFP,EP,EP2,SHO,CHO,THO,PSIO,CSHO,CTHO
      REAL*8 COEFF, COEFP, EP, EP2, SHO, CHO, THO, PSIO, CSHO, CTHO
C
      INTEGER NITER,IX,IY
      REAL*8 FUN1, ASINH, Z, DACSH
      REAL*8 RTMAX, VTMAX
      REAL*8 AAXIS, BAXIS, BAXIS2, QUADM, BARMAS, BARDEN, BULMAS, BULDEN
      REAL*8 CENMAS
      REAL*8 RQ, FMAS, AMAS, DMAS, DDMAS, QQ, AT, AX, AY, AZ
      REAL*8 PIGROC
      REAL*8 xpos,ypos,rpos,pot,taper
C
      REAL*8 PI, FOURPI, GRAVC
      PARAMETER ( PI = 3.141592654, FOURPI = 4. * PI, GRAVC = 4.29569)

c doosu Englmaier
      REAL*8 Md,Rd
      COMMON /ENGL/ Md,Rd
c doosu end

C
      FUN1 (Z) = SQRT( Z * Z + 1.)
      ASINH (Z) = LOG( Z + SQRT( 1.0 + Z * Z))
c
c set constants from the initpotential par-list
c   1 = pattern speed (returned)
c   2 = ignored in this routine  (should be 1 though)
c   3 = index (0,1,2)
c   4 = radlag
c
c      COMMON/ AREAD/FCMASS,AXIRAT,QUAD,CENDEN,RADLAG, INDEX
c      REAL*8 FCMASS, AXIRAT, QUAD, CENDEN, RADLAG
c      INTEGER INDEX

c first, set defaults:
      index = 1
      axirat = 2.5
      radlag = 6.0
      quad = 45000.0
      cenden = 24000.0

c override the first npar:
c Note: input values for omega (par(1)) and amode (par(2)) 
c       not used in this version)
      IF (npar.GE.3) index = par(3)
      IF (npar.GE.4) axirat = par(4)
      IF (npar.GE.5) radlag = par(5)
      IF (npar.GE.6) quad = par(6)
      IF (npar.GE.7) cenden = par(7)
      IF (npar.GE.8) rb = par(8)
      IF (npar.GE.9) gamm = par(9)
      IF (npar.GE.10) Mbh = par(10)
      EPSI = 0.001
C
C
C  GRAVC IN UNITS OF km/s, kpc and 1e6MSUN
C  not 'KPC**3 GYR**-2 1.E6MSUN**-1' as was quoted in earlier versions
c
c englmaier parameter
      Md = 2.3e5
      Rd = 14.1
C
C  This is the normalization of the model (at 20 kpc Vmx = 164.204 km/s)
C
      RTMAX = 14.1
c      RTMAX = 50.
      RTI2 = 1. / (RTMAX * RTMAX)
c      VTMAX = 164.204
      VTMAX = 260.
c      VTMAX = 0.0
      VT2 = VTMAX * VTMAX
      VTDR2= VT2 * RTI2
C
      FCMASS = 1.0
      AAXIS = 5.0
      CENMAS = 4.87333E4
      CENMAS = CENMAS * FCMASS
      BAXIS = AAXIS / AXIRAT
      BAXIS2 = BAXIS ** 2
      IF( INDEX .EQ. 0) QUADM = CENMAS * (AAXIS ** 2 - BAXIS2) / 5.0
      IF( INDEX .EQ. 1) QUADM = CENMAS * (AAXIS ** 2 - BAXIS2) / 7.0
      IF( INDEX .EQ. 2) QUADM = CENMAS * (AAXIS ** 2 - BAXIS2) / 9.0
      QUAD = MIN(QUAD, QUADM)
      BARMAS = CENMAS * QUAD / QUADM
      IF(INDEX .EQ. 0) BARDEN = BARMAS * 0.2387324    / (AAXIS * BAXIS2)
      IF(INDEX .EQ. 1) BARDEN = BARMAS * 0.5968310366 / (AAXIS * BAXIS2)
      IF(INDEX .EQ. 2) BARDEN = BARMAS * 1.044454314  / (AAXIS * BAXIS2)
      BULMAS = CENMAS - BARMAS
      BULDEN = 3.0E-3 * BULMAS / FOURPI
      CENDEN = MAX(CENDEN, BARDEN + BULDEN)
      BULDEN = CENDEN - BARDEN
      CORAD = 1.0
      BULMAS = CENMAS - BARMAS

      BARSTR = BARMAS / BULMAS

      NITER = 0
   10 RQ = 10.0 / CORAD
      AMAS = FOURPI * BULDEN * CORAD * CORAD * CORAD *
     1       (ASINH(RQ) - RQ / FUN1(RQ))
      DMAS = BULMAS - AMAS
      FMAS = DMAS / BULMAS
      IF (ABS(FMAS) .LE. 1.0E-6) GO TO 11
      IF (NITER .GE. 20) GO TO 13
      NITER = NITER + 1
      DDMAS = FOURPI * BULDEN * CORAD * CORAD * (3.0 * ASINH(RQ) -
     1        RQ * (3.0 + 4.0 * RQ * RQ) / (FUN1(RQ) ** 3))
      CORAD = MAX(0.01 * CORAD, CORAD + DMAS / DDMAS)
      GO TO 10
   11 CONTINUE
      C1BUL = GRAVC * BULDEN * FOURPI * CORAD * CORAD
c      C1BUL = 0.0
c      CONST=0.0
      CONST=12.56637061*GRAVC*CORAD*CORAD*CORAD*BULDEN
C
C     BAR
C
      A=AAXIS
      B=BAXIS
      AAIN = A * A
      AAI = 1. / AAIN
      C = BAXIS
      CCIN = C * C
      CCO = 1. / CCIN
      ELIPM = BARMAS
      ACONST = ELIPM * 4.027209375 / (A * C * C)
C
      IF( INDEX .EQ. 0) THEN
         PIGROC = PI * GRAVC * BARDEN
         EP = SQRT( AAXIS**2 - BAXIS2)
         EP2 = EP * EP
c oops, need to define DACSH here for INDEX=0 case
cpjt         PSIO = DACSH( AAXIS / BAXIS)
c
	 STOP
c
         SHO = SINH (PSIO) 
         CHO = AAXIS / BAXIS
         THO = SHO / CHO 
         CSHO = SHO ** 2 / EP ** 2
         CTHO = THO ** 2 / EP ** 2
         COEFF = 2. * PIGROC / (THO * SHO * SHO)
         COEFP = PIGROC *BAXIS2 / THO
      END IF
C
C
C     CALCULATE PATTERN SPEED
C
      RQ = RADLAG / CORAD
      AMAS = FOURPI * BULDEN * CORAD * CORAD * CORAD *
     1       (ASINH(RQ) - RQ / FUN1(RQ))
      QQ = GRAVC * AMAS / RADLAG
C  BAR NOW LIES ON THE Y AXIS. SO INVERT X AND Y TO HAVE THE SAME LAGRAN
C  POINT AS IN THE HYDRO CODE
      IF( INDEX .EQ. 0) CALL FORBA0( 0.D0, RADLAG, 0.D0,AY,AX,AZ)
      IF( INDEX .EQ. 1) CALL FORBA1( 0.D0, RADLAG, 0.D0,AY,AX,AZ)
      IF( INDEX .EQ. 2) CALL FORBA2( 0.D0, RADLAG, 0.D0,AY,AX,AZ)
      QQ = QQ - AX * RADLAG
      AT = RADLAG * VT2 *
     1     (1. / (1. + 1. * RADLAG * RADLAG * RTI2)) ** 1.5
      AT = AT * RTI2
      QQ = (QQ + AT * RADLAG) / (RADLAG * RADLAG)
      OMEGS = SQRT(QQ)
c wooyoung
      omegs = 33.067387747
c      omegs = 0.0
      par(1) = omegs
c                   set up the grid if nx>0 (rmax > 0 internal, else external)
c 	            can also choose a pattern speed overriding pgen::rl value
	IF (nx.GT.0 .and. grmax.gt.0.0) then
	   crpix1 = 0.5*(nx+1)
	   crpix2 = 0.5*(ny+1)
	   crval1 = 0.0
	   crval2 = 0.0
	   cdelt1 = 2.0*grmax/(nx-1)
	   cdelt2 = 2.0*grmax/(ny-1)
	   write(*,*) 'Internal grid: crpix,crval,cdelt:', 
     *                crpix1,crval1,cdelt1
	   DO 511 iy=1,ny
	      DO 511 ix=1,nx
		 xpos = (ix-crpix1)*cdelt1 + crval1
		 ypos = (iy-crpix2)*cdelt2 + crval2
		 CALL PISOB(ypos, xpos, 0.0d0, pot)
		 gpot0(ix,iy) = pot
		 gpot1(ix,iy) = 0.0
 511	      CONTINUE
	ELSE IF (nx.GT.0) THEN
	   write(*,*) 'Using an external grid of size ',nx,ny
	   IF(quad.GT.0.0) THEN
	     write(*,*) 'Now setting pot0 = pot0 + pot1; quad=',quad
             write(*,*) 'PATCHED: using new averaging in galaxy.src'
	     DO 513 iy=1,ny
	       DO 513 ix=1,nx
		 xpos = (ix-crpix1)*cdelt1 + crval1
		 ypos = (iy-crpix2)*cdelt2 + crval2
		 rpos = sqrt(xpos**2+ypos**2)
		 if (rpos.le.rtstart) then
                   taper = 1.0
		 else
                   taper = exp(-(rpos-rtstart)**2/(2*rtscale**2))
		 endif
c                write(*,*) 'taper: ',ix,iy,xpos,ypos,rpos,taper
c                setting taper=1 is needed to taper in the r-t  plane
c                see new code in galaxy.src (pjt 19-aug-2003)
                 taper = 1.0
                 gpot0(ix,iy) = gpot0(ix,iy) + gpot1(ix,iy) * taper
 513	     CONTINUE
	   ELSE
	     write(*,*) 'Using pot0, quad=',quad
	   ENDIF
	   if (gomegp.GT.0) then
		omegs = gomegp
		par(1) = gomegp
	   endif
	ENDIF
	write(*,*) 'Using pattern speed=',omegs
	write(*,*) 'rb=',CORAD
      print*, 'bar strength : ',BARSTR
c doosu write parameter at cmhogout
      if (QUAD .eq. 0) then
         write(2,*) 'RTMAX    =    ',RTMAX 
         write(2,*) 'VTMAX    =    ',VTMAX 
         write(2,*) 'PATTERN1 =    ',omegs
      else
         write(2,*) 'PATTERN2 =    ',omegs
         write(2,*) 'BULMASS  =    ',BULMAS
         write(2,*) 'BARMASS  =    ',BARMAS
         write(2,*) 'BULDEN   =    ',BULDEN
         write(2,*) 'BARDEN   =    ',BARDEN
         write(2,*) 'Rb       =    ',CORAD
         write(2,*) 'BARSTR   =    ',BARSTR
      endif 
      print*,'Rb     =  ', CORAD
      print*,'BULDEN =  ', BULDEN 
      print*,'BARDEN =  ', BARDEN 
c doosu end

      RETURN
      IF(QUAD.GT.QUADM) THEN
	PRINT 612,QUAD,QUADM
      ENDIF
13    PRINT 613, CORAD,FMAS,NITER
      PRINT 612, QUAD,QUADM
612   FORMAT(' Warning: QUAD too large: ',1E12.5,' > ',1E12.5)
613   FORMAT(' ** NO CONVERGENCE FOR CORE RADIUS; RADIUS=',1PE12.5,
     &       ' FRACTION MASS DIFFERENCE =',E12.5,' NITER =',I5)

      STOP
      END
c----------------------------------------------------------------------
      SUBROUTINE potential(ndim, pos, acc, pot, time)
      INTEGER ndim
      REAL*8 pos(*), acc(*), pot, time

      INTEGER nx,ny
      REAL gpot0,gpot1,crpix1,crpix2,cdelt1,cdelt2,crval1,crval2,grmax
      REAL gomegp,rhalo,vhalo,gamma,rtstart,rtscale
      COMMON/ggrid/nx,ny,gpot0(1030,1030),gpot1(1030,1030),
     *             grmax,gomegp,
     *             crpix1,crpix2,cdelt1,cdelt2,crval1,crval2,
     *             rhalo,vhalo,gamma,rtstart,rtscale

      INTEGER ix,iy
      REAL x,y,dx,dy,x3,y3,x4,y4,x8,y8,tmp

      IF(nx.gt.0) THEN
	 x = (pos(1)-crval1)/cdelt1 + crpix1 - 0.5
	 y = (pos(2)-crval2)/cdelt2 + crpix2 - 0.5
	 ix = nint(x-0.5)
	 iy = nint(y-0.5)
         if (ix.le.1 .or. ix.ge.nx .or. iy.le.1 .or. iy.ge.ny) then
            write(*,*) 'Grid not large enough for ',
     *                  pos(1),pos(2),ix,iy,nx,ny
            stop
         endif
	 dx=x-ix-0.5
	 dy=y-iy-0.5
	 x3=3*dx
	 y3=3*dy
	 x4=4*dx
	 y4=4*dy
	 x8=2*x4
	 y8=2*y4
	 ix=ix+1
	 iy=iy+1
	 pot = gpot0(ix,iy)
         acc(2) = -( (-2+x4+y3)*gpot0(ix-1,iy-1)
     *              +(  -x8   )*gpot0(ix  ,iy-1)
     *              +( 2+x4-y3)*gpot0(ix+1,iy-1)
     *              +(-2+x4   )*gpot0(ix-1,iy)
     *              +(  -x8   )*gpot0(ix  ,iy)
     *              +( 2+x4   )*gpot0(ix+1,iy)
     *              +(-2+x4-y3)*gpot0(ix-1,iy+1)
     *              +(  -x8   )*gpot0(ix  ,iy+1)
     *              +( 2+x4+y3)*gpot0(ix+1,iy+1) )/(12.0*cdelt1)
	 acc(1) = -( (-2+x3+y4)*gpot0(ix-1,iy-1)
     *              +(-2   +y4)*gpot0(ix  ,iy-1)
     *              +(-2-x3+y4)*gpot0(ix+1,iy-1)
     *              +(     -y8)*gpot0(ix-1,iy)
     *              +(     -y8)*gpot0(ix  ,iy)
     *              +(     -y8)*gpot0(ix+1,iy)
     *              +( 2-x3+y4)*gpot0(ix-1,iy+1)
     *              +( 2   +y4)*gpot0(ix  ,iy+1)
     *              +( 2+x3+y4)*gpot0(ix+1,iy+1) )/(12.0*cdelt2)

	 acc(1) = acc(1)*gamma
	 acc(2) = acc(2)*gamma
	 if (vhalo.gt.0.0) then
	    tmp = -vhalo*vhalo/(pos(1)**2 + pos(2)**2 + rhalo*rhalo)
	    acc(1) = acc(1) + pos(1)*tmp
	    acc(2) = acc(2) + pos(2)*tmp
	 endif
c	 write(*,*) 'grid: ',pos(1),pos(2),pot,acc(1),acc(2)
c	 write(*,*) ' c  : ',ix,iy,x,y
c        write(*,*) ' d  : ',dx,dy,x3,x4,x8,y3,y4,y8
c-----------------------------------------------------------------------
c        write(*,*) 'p: ',gpot(ix-1,iy+1),gpot(ix,iy+1),gpot(ix+1,iy+1)
c        write(*,*) 'p: ',gpot(ix-1,iy  ),gpot(ix,iy  ),gpot(ix+1,iy )
c        write(*,*) 'p: ',gpot(ix-1,iy-1),gpot(ix,iy-1),gpot(ix+1,iy-1)
c        write(*,101) gpot(ix-1,iy+1),gpot(ix,iy+1),gpot(ix+1,iy+1)
c        write(*,101) gpot(ix-1,iy  ),gpot(ix,iy  ),gpot(ix+1,iy )
c        write(*,101) gpot(ix-1,iy-1),gpot(ix,iy-1),gpot(ix+1,iy-1)
 101     format('p: ',3f11.2)
	 RETURN
      ENDIF
      acc(1) = 0.0
      acc(2) = 0.0
      acc(3) = 0.0
      pot = 0.0
c			note exchanged X and Y, since their bar
c			is along the Y axis, ours along the X
      CALL FISOB(pos(1), pos(2), pos(3), acc(1), acc(2), acc(3))
      CALL PISOB(pos(1), pos(2), pos(3), pot)
c      write(*,*) 'function: ',pos(1),pos(2),pot,acc(1),acc(2)
      END
C-----------------------------------------------------------------------
      SUBROUTINE FISOB( X, Y, Z, FX, FY, FZ)
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     CALCUL DES FORCES DANS LE POTENTIEL (PISOB)
C     ROUTINE APPELLEE PAR DISOB (DERIVEES POUR RK78)
C
      IMPLICIT NONE
      REAL*8 X, Y, Z, FX, FY, FZ
C
      COMMON/ AREAD/FCMASS,AXIRAT,QUAD,CENDEN,RADLAG, INDEX
      REAL*8 FCMASS, AXIRAT, QUAD, CENDEN, RADLAG
      INTEGER INDEX

      REAL*8 FXBAR,FYBAR
      COMMON /CHECKINIT2/ FXBAR,FYBAR


C 
      REAL*8 FAX, FAY, FAZ, FBX, FBY, FBZ
      CALL FORAX( X, Y, Z, FAX, FAY, FAZ) 
      IF( QUAD .GT. 1.E-9) THEN
         IF( INDEX .EQ. 0) CALL FORBA0( X, Y, Z, FBX, FBY, FBZ) 
         IF( INDEX .EQ. 1) CALL FORBA1( X, Y, Z, FBX, FBY, FBZ) 
         IF( INDEX .EQ. 2) CALL FORBA2( X, Y, Z, FBX, FBY, FBZ) 
         FX = FAX + FBX
         FY = FAY + FBY
         FZ = 0.

         FXBAR = FBX
         FYBAR = FBY

      ELSE
         FX = FAX
         FY = FAY
         FZ = 0.

         FXBAR = 0.
         FYBAR = 0.

      END IF
C 

      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FORAX( X, Y, Z, FX, FY, FZ)
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C FORCES DUES A UN POTENTIEL AXISYMETRIQUE
C
      IMPLICIT NONE
      REAL*8 X, Y, Z, FX, FY, FZ
C
      COMMON/AXISY/C1BUL,CONST,CORAD,VT2,VTDR2,RTI2,Mbh,EPSI,rb,gamm
     &             ,BULMAS
      REAL*8 C1BUL, CONST, CORAD, VT2, VTDR2, RTI2,Mbh,EPSI,rb,gamm
      REAL*8 BULMAS
C
      REAL*8 FUN1, ASINH, ZZ
      REAL*8 R2, R, RQ, QQ

      REAL*8 PI, FOURPI, GRAVC
      PARAMETER ( PI = 3.141592654, FOURPI = 4. * PI, GRAVC = 4.29569)
c doosu englmaier
      REAL*8 Md,Rd
      COMMON /ENGL/ Md,Rd
c doosu end

      REAL*8 QQBUL,QQDISK
      REAL*8 FXBUL,FYBUL,FXDISK,FYDISK
      COMMON /CHECKINIT/ FXBUL,FYBUL,FXDISK,FYDISK 

c wooyoung central black hole start
      REAL*8 QQBH,FXBH,FYBH
c wooyoung end
      

      FUN1 (ZZ) = SQRT( ZZ * ZZ + 1.)
      ASINH (ZZ) = LOG( ZZ + FUN1 (ZZ))
      R2 = X * X + Y * Y
      IF (R2.GT.0.0D0) THEN
         R = SQRT (R2)
         RQ = R / CORAD
         QQ = CONST * (ASINH(RQ) - RQ / FUN1(RQ))/(R * R * R) +
c wooyoung central black hole start
     1        GRAVC * Mbh * ( R2 + EPSI*EPSI )**(-1.5) +
c wooyoung end              
c for Toomre Test

c    1        VTDR2 * (3./ (1. + 2. * R * R * RTI2)) ** 1.5
     1        VTDR2 * (1./ (1. + R * R * RTI2)) ** 1.5
c     *        + 191.*191./R**2.
c     *        + 10.0*10.0/(R+0.1)/R**2.
c         print*,'piner disk force'

c         QQDISK = VTDR2 * (3./ (1. + 2. * R * R * RTI2)) ** 1.5
          QQDISK = VTDR2 * (1./ (1. + R * R * RTI2)) ** 1.5
         FX = - X * QQ
         FY = - Y * QQ
         QQBUL  = CONST * (ASINH(RQ) - RQ / FUN1(RQ))/(R * R * R)
         FXBUL  = - X * QQBUL
         FYBUL  = - Y * QQBUL
         FXDISK = - X * QQDISK
         FYDISK = - Y * QQDISK
      ELSE
         FX = 0.
         FY = 0.
         FXBUL  = 0.
         FYBUL  = 0.
         FXDISK = 0.
         FYDISK = 0.
      ENDIF
      FZ = 0.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FORBA1( XINPUT, YINPUT, Z, AXOUT, AYOUT, AZ)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C FORCES FROM A FERRERS BAR WITH INDEX = 1
C     SUBROUTINE PROL6 (ELIPM,A,B,C,X,Y,Z,AX,AY,AZ,EPT)
C      IMPLICIT real*8 (A-Z)
      IMPLICIT NONE
      REAL*8 XINPUT, YINPUT, Z, AXOUT, AYOUT, AZ
C
      COMMON/PROLA/ELIPM,A,B,C,AAIN,CCIN,AAI,CCO,ACONST
      REAL*8 ELIPM, A, B, C, ACONST, CCO, CCIN, AAI, AAIN
C
      REAL*8 X, Y, AA, CC, XX, YY, RR, GI, SIDE, BK, CK
      REAL*8 KAPPA, EM, EE, E, LNE, A3, A1, I, A33, A13, A11
      REAL*8 X4, Y4, AX, AY
C     CAREFUL : IN HYDRO PROGRAM THE BAR IS ALONG THE X AXIS
C     IN THE ORBIT CALCULATIONS ALONG THE Y AXIS
C
c      print*,'ELIPM, A, B, C, ACONST, CCO, CCIN, AAI, AAIN'
c      print*, ELIPM, A, B, C, ACONST, CCO, CCIN, AAI, AAIN

      X = XINPUT
      Y = YINPUT
      AA = AAIN
      CC = CCIN
      XX = X*X
      YY = Y*Y
      RR = XX+YY
      GI = 1.
      SIDE = XX/CC + (YY)/AA - 1.
      IF (SIDE .LE. 0.) GO TO 10
         BK = RR-AA-CC
         CK = AA*CC*SIDE
c         CK = XX*AA+YY*CC-AA*CC 
         KAPPA = (BK+SQRT(BK*BK+4.*CK))/2.
         AA = AA + KAPPA
         CC = CC + KAPPA
         GI = A*C*C/SQRT(AA)/CC
   10 EM = CC/AA
      EE = 1. - EM
      E = SQRT(EE)
c      LNE = LOG((1.+E)/(1.-E))
C
C     CALCULATE INDEX SYMBOLS IN EASIEST WAY
C
c      A3 = GI*(1.-EM*LNE/(E + E))/EE
c      A1 = 2.*(GI - A3)
c      I = AA*A1 + 2.*CC*A3
 
c      A33 = GI*(4.*EE/EM - 6. + 3.*EM*LNE/E)/(8.*AA*EE*EE)
c      A13 = 2.*(GI/CC - 2.*A33)
c      A11 = 2.*(GI/AA - A13)/3.
 
c      X4 = A11*XX + A13*(YY)
c      Y4 = A13*XX + A33*(YY)
C     Z4 = Y4
c     wooyoung

      LNE = LOG((1.+E)/(sqrt(EM)))
      A3 = 2.0*(LNE-E)/(A*A-C*C)**1.5
      A1 = EE/EM*E/(A*A-C*C)**1.5 - A3/2.0
      I = 2.0*LNE/sqrt(A*A-C*C)
 
      A13 = (A1 - A3)/(A*A-C*C)
      A33 = 2.0/3.0*(1.0/(sqrt(AA)*CC*AA) - A13)
      A11 = 0.25*(2.0/(sqrt(AA)*CC*CC) - A13)
 
      X4 = A11*XX + A13*(YY)
      Y4 = A13*XX + A33*(YY)

c     wooyoung end 


C     EPT = I - 2.*(A1*XX + A3*(YY)) +(XX*X4+YY*Y4)
 
      AX = 4. * X * (X4 - A1)
 
      AY = 4. * Y * (Y4 - A3)
 
      ACONST = ELIPM*4.027209375/(A*C*C)
C     EPT = -EPT * ACONST
      AX = AX * ACONST * A * C * C
      AY = AY * ACONST * A * C * C
      AZ = 0.
      AXOUT=AX
      AYOUT=AY
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FORBA0( X, Y, Z, FX, FY, FZ)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X, Y, Z, FX, FY, FZ
C
      COMMON/PROLA2/COEFF,COEFP,EP,EP2,SHO,CHO,THO,PSIO,CSHO,CTHO
      REAL*8 COEFF, COEFP, EP, EP2, SHO, CHO, THO, PSIO, CSHO, CTHO
C
      REAL*8 R2, Y2, XM, PSI, SH, CH, TH, AA, W10, W11, W20
      R2 = X * X
      Y2 = Y * Y
      XM = SQRT( R2 * CSHO + Y2 * CTHO)
      IF( XM .LE. 1.) THEN
         PSI = PSIO
         SH = SHO
         CH = CHO
         TH = THO
      ELSE
         AA = Y2 + R2 - EP2 
         IF( R2 * EP2 .LE. 1.D-10 * AA) THEN
            SH = EP / SQRT( Y2 - EP2)
         ELSE
            SH = SQRT(( SQRT( AA * AA + 4. * R2 * EP2) - AA) / (R2+R2))
         END IF
         PSI = LOG( SH + SQRT( SH ** 2 + 1.)) 
         CH = COSH (PSI) 
         TH = SH / CH
      END IF
      W10 = PSI + PSI 
      W20 = .5 * W10 - CH * SH
      W11 = TH + TH - W10 
      FX = COEFF * X * W20
      FY = COEFF * Y * W11
      FZ = 0.
      RETURN
      END 
C-----------------------------------------------------------------------
      SUBROUTINE PISOB( X, Y, Z, PHI)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     POTENTIEL AXISYMETRIQUE + BARRE PROLATE
      IMPLICIT NONE
      REAL*8 X, Y, Z, PHI
C
      COMMON/ AREAD/FCMASS,AXIRAT,QUAD,CENDEN,RADLAG, INDEX
      REAL*8 FCMASS, AXIRAT, QUAD, CENDEN, RADLAG
      INTEGER INDEX
C
      REAL*8 PAX, PBA
C
      CALL POTAX( X, Y, Z, PAX) 
      IF ( INDEX .EQ. 0) CALL POTBA0 ( X, Y, Z, PBA)
      IF ( INDEX .EQ. 1) CALL POTBA1 ( X, Y, Z, PBA)
      IF ( INDEX .EQ. 2) CALL POTBA2 ( X, Y, Z, PBA)
      PHI = PAX + PBA 
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE POTAX( X, Y, Z, POT)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C POTENTIEL AXISYMETRIQUE DU BULBE ET DU DISQUE
      IMPLICIT NONE
      REAL*8 X, Y, Z, POT
C
      COMMON/AXISY/C1BUL,CONST,CORAD,VT2,VTDR2,RTI2,Mbh,EPSI,rb,gamm,
     &       BULMAS
      REAL*8 C1BUL, CONST, CORAD, VT2, VTDR2, RTI2,Mbh,EPSI,rb,gamm
      REAL*8 BULMAS
C
      REAL*8 R2, R, RREL, RREL2, POTBU, POTDI, POTBH

      REAL*8 PI, FOURPI, GRAVC
      PARAMETER ( PI = 3.141592654, FOURPI = 4. * PI, GRAVC = 4.29569)
c doosu englmaier
      REAL*8 Md,Rd
      COMMON /ENGL/ Md,Rd
c doosu end

      R2 = X * X + Y * Y
      R = SQRT (R2)
C     BULGE
      RREL = R / CORAD
      RREL2 = RREL * RREL
      POTBU = - C1BUL * LOG( RREL + SQRT( RREL2 + 1.)) / RREL
      POTBH = - GRAVC * Mbh  / SQRT( R2 + EPSI*EPSI )

c for Toomre Test
C     DISK
c  piner(1995) disk
c      POTDI = - 1.5 * VT2 * SQRT( 1.5 / ( 0.5 + R2 * RTI2))
      POTDI = - VT2 * SQRT( 1. / ( 1. + R2 * RTI2))
c       POTDI = 191.*191.*LOG(R)
c     *        - 10.0*10.0*LOG(10.0/(R+0.2))
c      print*,'piner disk potential'
C     TOTAL
      POT = POTBU + POTDI + POTBH

      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE POTBA1( XINPUT, YINPUT, Z, EPT)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C POTENTIAL OF A FERRERS BAR WITH INDEX 1
C     SUBROUTINE PROL6 (ELIPM,A,B,C,X,Y,Z,AX,AY,AZ,EPT)
C      IMPLICIT real*8 (A-Z)
      IMPLICIT NONE
      real*8 XINPUT, YINPUT, Z, EPT
C
      COMMON/PROLA/ELIPM,A,B,C,AAIN,CCIN,AAI,CCO,ACONST
      REAL*8 ELIPM, A, B, C, ACONST, CCO, CCIN, AAI, AAIN
C
      REAL*8 X, Y, AA, CC, XX, YY, RR, GI, SIDE, BK, CK
      REAL*8 KAPPA, EM, EE, E, LNE, A3, A1, I, A33, A13, A11
      REAL*8 X4, Y4
C
C     CAREFUL : IN HYDRO PROGRAM THE BAR IS ALONG THE X AXIS
C     IN THE ORBIT CALCULATIONS ALONG THE Y
C
      X = XINPUT
      Y = YINPUT
      AA = AAIN
      CC = CCIN
      XX = X*X
      YY = Y*Y
      RR = XX+YY
      GI = 1.
      SIDE = XX/CC + (YY)/AA - 1.
      IF (SIDE .LE. 0.) GO TO 10
         BK = RR-AA-CC
         CK = AA*CC*SIDE
c         CK = XX*AA+YY*CC-AA*CC 
         KAPPA = (BK+SQRT(BK*BK+4.*CK))/2.
         AA = AA + KAPPA
         CC = CC + KAPPA
         GI = A*C*C/SQRT(AA)/CC
   10 EM = CC/AA
      EE = 1. - EM
      E = SQRT(EE)
c      LNE = LOG((1.+E)/(1.-E))
C
C     CALCULATE INDEX SYMBOLS IN EASIEST WAY
C
c      A3 = GI*(1.-EM*LNE/(E + E))/EE3 = GI*(1.-EM*LNE/(E + E))/EE
c      A1 = 2.*(GI - A3)
c      I = AA*A1 + 2.*CC*A3

c      A33 = GI*(4.*EE/EM - 6. + 3.*EM*LNE/E)/(8.*AA*EE*EE)
c      A13 = 2.*(GI/CC - 2.*A33)
c      A11 = 2.*(GI/AA - A13)/3.

c      X4 = A11*XX + A13*(YY)
c      Y4 = A13*XX + A33*(YY)

c      EPT = I - 2.*(A1*XX + A3*(YY)) +(XX*X4+YY*Y4)

c     wooyoung

      LNE = LOG((1.+E)/(sqrt(EM)))
      A3 = 2.0*(LNE-E)/(A*A-C*C)**1.5
      A1 = EE/EM*E/(A*A-C*C)**1.5 - A3/2.0
      I = 2.0*LNE/sqrt(A*A-C*C)
 
      A13 = (A1 - A3)/(A*A-C*C)
      A33 = 2.0/3.0*(1.0/(sqrt(AA)*CC*CC) - A13)
      A11 = 0.25*(2.0/(sqrt(AA)*CC*AA) - A13)
 
      X4 = A11*XX + A13*(YY)
      Y4 = A13*XX + A33*(YY)

c     wooyoung end 

      EPT = I - 2.*(A1*XX + A3*(YY)) +(XX*X4+YY*Y4)
 
C     ACONST = ELIPM*4.027209375/(A*C*C)
      EPT = -EPT * ACONST*A*C*C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE POTBA0( X, Y, Z, POT) 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     POTENTIEL DUN SPHEROIDE HOMOGENE PROLAT 
C     PARAMETRES
C     A : PETIT AXE 
C     C : GRAND AXE 
C     EP : EXCENTRICITE   EP = SQRT(C * C - A * A)
C     SHO = SINH(PSIO) .........
      IMPLICIT NONE
      REAL*8 X, Y, Z, POT
C
      COMMON/PROLA2/COEFF,COEFP,EP,EP2,SHO,CHO,THO,PSIO,CSHO,CTHO
      REAL*8 COEFF, COEFP, EP, EP2, SHO, CHO, THO, PSIO, CSHO, CTHO
C
      REAL*8 R2, Y2, XM, PSI, SH, CH, TH, AA, W10, W11, W20
      R2 = X * X
      Y2 = Y * Y
      XM = SQRT( R2 * CSHO + Y2 * CTHO)
      IF( XM .LE. 1.) THEN
         PSI = PSIO
         SH = SHO
         CH = CHO
         TH = THO
      ELSE
         AA = Y2 + R2 - EP2 
         IF( R2 * EP2 .LE. 1.D-30*AA) THEN
            SH = EP / SQRT( Y2 - EP2)
         ELSE
            SH = SQRT(( SQRT( AA * AA + 4. * R2 * EP2) - AA) / (R2+R2))
         END IF
         PSI = LOG( SH + SQRT( SH ** 2 + 1.)) 
         CH = COSH (PSI) 
         TH = SH / CH
      END IF
      W10 = PSI + PSI 
      W20 = .5 * W10 - CH * SH
      W11 = TH + TH - W10 
      POT= ((W20 * R2 + W11 * Y2) / EP2 + W10) * COEFP
      POT = - POT
      RETURN
      END 
c-----------------------------------------------------------------------
      SUBROUTINE forba2 (xin,yin,z,axout,ayout,az)
      IMPLICIT DOUBLE PRECISION (a-z)

      COMMON/PROLA/ELIPM,A,B,C,pad1,pad2,pad3,pad4,pad5
      REAL*8 ELIPM, A, B, C,pad1,pad2,pad3,pad4,pad5

	x=yin
	y=xin
      aa=a*a
      cc=c*c
      xx=x*x
      yy=y*y
      zz=z*z
      rr=xx+yy+zz
      side=xx/aa+(yy+zz)/cc - 1.0
      gi=1.0
      IF (side.LE.0.0) GOTO 100
           kappa=(rr-aa-cc+DSQRT((rr-aa-cc)**2+4*aa*cc*side))*0.5
           aa=aa+kappa
           cc=cc+kappa
           gi=a*c*c/(DSQRT(aa)*cc)
 100  CONTINUE
      em=cc/aa
      ee=1.0-em
      E=dsqrt(ee)
      lne=DLOG((1.0+e)/(1.0-e))
C                              calculate index symbols
      a3=gi*(1.0-em*lne/e*0.5)/ee
      a1=2.0*(gi-a3)
      i=aa*a1+2.0*cc*a3

      a33=gi*(4*ee/em-6+3*em*lne/e)/(8*aa*ee**2)
      a13=2.0*(gi/cc-2.0*a33)
      a11=2.0*(gi/aa-a13)/3.0

      a333=gi*(30-20*ee/em+16*ee**2/em**2-15*em*lne/e)
      a333=a333/(ee**3*aa**2*48)
      a133=2.0*(gi/cc**2-3.0*a333)
      a113=2.0*(gi/(aa*cc)-2.0*a133)/3.0
      a111=2.0*(gi/aa**2-a113)*0.2

      x6=a111*xx+3.0*a113*(yy+zz)
      y6=3.0*a133*xx+a333*(yy+3.0*zz)
      z6=3.0*a133*xx+a333*(3.0*yy+zz)

      x4=a11*xx+a13*(yy+zz)
      y4=a13*xx+a33*(yy+zz)
      z4=y4

      ax=-6.0*x*a1+12.0*x*x4-2.0*x*(2.0*xx*x6+6.0*a133*yy*zz)
      ax=ax-2.0*x*(xx**2*a111+3.0*yy**2*a133+3.0*zz**2*a133)

      ay=-6.0*y*a3+12.0*y*y4-2.0*y*(2.0*yy*y6+6.0*a133*xx*zz)
      ay=ay-2.0*y*(3.0*xx**2*a113+yy**2*a333+3.0*zz**2*a333)

      az=-6.0*z*a3+12.0*z*z4-2.0*z*(2.0*zz*z6+6.0*a133*xx*yy)
      az=az-2.0*z*(3.0*xx**2*a113+3.0*yy**2*a333+zz**2*a333)

      aconst=elipm*35/(32*0.23259*a*c*c)

c flip them forces too !!!!

      axout=ay*aconst
      ayout=ax*aconst
      az=az*aconst
      return
      END


      SUBROUTINE potba2(xin,yin,z,ept)
      IMPLICIT DOUBLE PRECISION (a-z)

      COMMON/PROLA/ELIPM,A,B,C,pad1,pad2,pad3,pad4,pad5
      REAL*8 ELIPM, A, B, C, pad1,pad2,pad3,pad4,pad5

	x = yin
	y = xin
      aa=a*a
      cc=c*c
      xx=x*x
      yy=y*y
      zz=z*z
      rr=xx+yy+zz
      side=xx/aa+(yy+zz)/cc - 1.0
      gi=1.0
      IF (side.LE.0.0) GOTO 100
           kappa=(rr-aa-cc+DSQRT((rr-aa-cc)**2+4*aa*cc*side))*0.5
           aa=aa+kappa
           cc=cc+kappa
           gi=a*c*c/(DSQRT(aa)*cc)
 100  CONTINUE
      em=cc/aa
      ee=1.0-em
      E=dsqrt(ee)
      lne=DLOG((1.0+e)/(1.0-e))
C                               calculate index symbols
      a3=gi*(1.0-em*lne/e*0.5)/ee
      a1=2.0*(gi-a3)
      i=aa*a1+2.0*cc*a3

      a33=gi*(4*ee/em-6+3*em*lne/e)/(8*aa*ee**2)
      a13=2.0*(gi/cc-2.0*a33)
      a11=2.0*(gi/aa-a13)/3.0

      a333=gi*(30-20*ee/em+16*ee**2/em**2-15*em*lne/e)
      a333=a333/(ee**3*aa**2*48)
      a133=2.0*(gi/cc**2-3.0*a333)
      a113=2.0*(gi/(aa*cc)-2.0*a133)/3.0
      a111=2.0*(gi/aa**2-a113)*0.2

      x6=a111*xx+3.0*a113*(yy+zz)
      y6=3.0*a133*xx+a333*(yy+3.0*zz)
      z6=3.0*a133*xx+a333*(3.0*yy+zz)

      x4=a11*xx+a13*(yy+zz)
      y4=a13*xx+a33*(yy+zz)
      z4=y4

      ept = i - 3.0*(a1*xx+a3*(yy+zz)) + 3.0*(xx*x4+yy*y4+zz*z4)
     *         - (x6*xx**2+y6*yy**2+z6*zz**2+6.0*a133*xx*yy*zz)

      aconst=elipm*35/(32*0.23259*a*c*c)
      ept=-ept*aconst
      END



