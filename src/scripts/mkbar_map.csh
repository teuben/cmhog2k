#! /bin/csh -f
#
#  $Id$
#
#  This scripts takes an HDF output snapshot file from cmhog
#  (bar hydro, polar coordinates), projects it to a requested
#  sky view as to be able to compare it with an existing maps
#  (mostly meant for a velocity field). It will need a reference
#  map in order for pixel by pixel comparison, and then computes
#  a difference map between projected model and reference map
#  (the 'observation') with the intent to minimize sigma.
#
#  See also mkbar_cube.csh for a 3D version making cubes.
#
#  In the model the bar is assumed oriented along Y axis (PA=0), and 
#  flows CCW as seen from the positive Z axis.
#  For CW rotating galaxies the code adds 180 to 'inc' and/or 'pa'
#
#  Observations are assumed to have their reference pixel to
#  coincide with the center (ra0,dec0,vsys), although we expect
#  to relax this condition in a future version of this script.
#
#  The new convention for rot=1 means CCW, rot=-1 means CW, denoting
#  the sign of the galaxy angular momentum vector, where the positive
#  Z axis points to the observer (hence doppler recession is -vz).
#
#   9-may-03  1.4 Derived from mkbar_cube.csh with better WCS and refmap
#     sep-03  1.5 added documentation
#  10-nov-04  1.6 merged two seemingly different versions where a $ring file was added
#  14-nov-04  1.7 optional map renormalization for more proper chi2 computation

set version=15-nov-04

if ($#argv == 0) then
  echo Usage: $0 in=HDF_DATASET out=BASENAME refmap=FITSFILE...
  echo Version: $version
  echo Gridding and projecting 2D CMHOG hydro models to specified bar viewing angles
  echo Creates maps to be compared to a refman
  echo Optional parameters:  
  echo "   pa, inc, phi, rot (1=ccw)"
  echo "   rmax, n, beam, color, clean, cube, denlog"
  echo "   wcs, pscale, vscale, vsys"
  echo "   par, inden"
  echo "You also need the NEMO environment"
  exit 0
endif

# 			Required Keywords
unset in
unset out
unset refmap
# 			Geometry (defaults are for some reasonable galaxy)
set pa=30
set inc=60
set phi=30
set rot=1
#                       Spatial gridding (n cells from -rmax : rmax)
set rmax=6
set n=200
set beam=0.25
set color=1
set clean=1
set denlog=0

#                       WCS definition (and mapping) if derived from an observation (cube)
set wcs=1

set pscale=1
set vscale=1
set vsys=0

#                       extra scaling for hydro data into the snapshot (usefull if no wcs)
set hpscale=1
set hvscale=1

#                       optional sigma scaling (set to 1 if none) of velocity difference map
set sigma=5
set npar=1
set nppb=1


#
set par=""
set inden=""
#			Parse commandline to (re)set keywords
foreach a ($*)
  set $a
  if (-e "$par") then
    source $par
  else if (X != X$par) then
    echo Warning, par=$par does not exist
  endif
  set par=""
end

#               Fixed constants (p: arcsec -> degrees    v: km/s -> m/s)
#               probably should not be changed
set puscale=2.77777777777e-4
set vuscale=1e3

#
#  fix inc/pa for ccw(rot=1) or cw(rot=-1) cases for NEMO's euler angles
#
if ($rot == -1) then
   set inc=$inc+180
else if ($rot == 1) then
   set pa=$pa+180
else
   echo "Bad rotation, must be 1 (ccw) or -1 (cw)"
   exit 1
endif
#		Report
echo     Files: in=$in out=$out 
echo -n "Grid: rmax=$rmax n=$n beam=$beam "
if ($wcs) then
   echo \(`nemoinp "$beam*$pscale"` arcsec\)
else
   echo ""
endif
echo Projection: phi=$phi inc=$inc pa=$pa
#               Derived quantities
set cell=`nemoinp "2*$rmax/$n"`
set range="-${rmax}:${rmax}"
echo -n "      Derived: cell=$cell"
if ($wcs) then
   echo \(`nemoinp "2*$rmax/$n*$pscale"` arcsec\)
else
   echo ""
endif

if ($wcs) then
    #   get the wcs from the reference map
    set nx=(`fitshead $refmap | grep ^NAXIS1 | awk '{print $3}'`)
    set ny=(`fitshead $refmap | grep ^NAXIS2 | awk '{print $3}'`)
    set px=(`fitshead $refmap | grep ^CRPIX1 | awk '{print $3}'`)
    set py=(`fitshead $refmap | grep ^CRPIX2 | awk '{print $3}'`)
    set dx=(`fitshead $refmap | grep ^CDELT1 | awk '{print $3}'`)
    set dy=(`fitshead $refmap | grep ^CDELT2 | awk '{print $3}'`)
    set rx=(`fitshead $refmap | grep ^CRVAL1 | awk '{print $3}'`)
    set ry=(`fitshead $refmap | grep ^CRVAL2 | awk '{print $3}'`)

    # set for output fits file
    set crpix=$px,$py
    set cdelt=$dx,$dy
    set crval=$rx,$ry

    echo CRPIX:   $crpix
    echo CRVAL:   $crval
    echo CDELT:   $cdelt
    set wcspars=(crpix=$crpix,1 crval=$crval,0 cdelt=$cdelt,1 radecvel=t)
    set wcspars=(refmap=$refmap)

    set xmin=`nemoinp "$dx*($px-0.5)*3600/$pscale"`
    set xmax=`nemoinp "$dx*($px-$nx-0.5)*3600/$pscale"`
    set ymin=`nemoinp "$dy*(0.5-$py)*3600/$pscale"`
    set ymax=`nemoinp "$dy*($ny+0.5-$py)*3600/$pscale"`

    echo Xrange:  $xmin  $xmax
    echo Yrange:  $ymin  $ymax

    set xrange=${xmin}:${xmax}
    set yrange=${ymin}:${ymax}
else
    echo "### Fatal error: this script does not work without a refmap yet"
    exit
    set wcspars=()
endif


set comment="$0 $in $out $pa,$inc,$phi,$range,$n,$beam"

echo WCSPARS: $wcspars
echo COMMENT: $comment

#> nemo.need tabtos tabmath snaptrans snaprotate snapadd snapgrid ccdsmooth ccdmath ccdfits fitshead

set tmp=tmp$$
if (! -e $out.den.fits) then

    # convert the half-plane HDF file to a full plane snapshot file

    tsd in=$in out=$tmp.tab coord=t
    if (-e "$inden") then
      echo "Taking densities from $inden instead in $in"
      tsd in=$inden out=$tmp.tab.den coord=t
      mv $tmp.tab  $tmp.tab.vel
      tabmath $tmp.tab.vel,$tmp.tab.den $tmp.tab %1,%2,%3,%4,%10 all
    endif
    if ($status) goto cleanup
    tabtos in=$tmp.tab out=$tmp.s0 block1=x,y,vx,vy,mass
    snaptrans in=$tmp.s0 out=$tmp.s1 ctypei=cyl ctypeo=cart
    snaprotate in=$tmp.s1 out=$tmp.s2 theta=180 order=z
    snapadd $tmp.s1,$tmp.s2 $tmp.s3

    # project for skyview, and create a intensity and velocity field

    snaprotate $tmp.s3 - \
        "atand(tand($phi)/cosd($inc)),$inc,$pa" zyz |\
	snapscale - - rscale=$hpscale vscale=$hvscale |\
	snapshift - $tmp.snap vshift=0,0,$vsys mode=sub

    echo -n "Projected model velocities:"
    snapprint $tmp.snap -vz | tabhist - tab=t |& grep min

    foreach mom (0 1 2)
         snapgrid in=$tmp.snap out=$tmp.$mom \
                xrange=$xrange yrange=$yrange nx=$nx ny=$ny moment=$mom mean=t
         ccdsmooth in=$tmp.$mom out=$tmp.$mom.s gauss=$beam
    end
    ccdmath $tmp.1.s,$tmp.0.s $tmp.vel %1/%2
    ## BUG: ifgt() doesn't work
    ##    ccdmath $tmp.0.s - "ifgt(%1,0,log(%1),-10)" | ccdfits - $out.den.fits
    ccdmath $tmp.0.s - "log(%1)" | ccdmath - - "ifeq(%1,0,-10,%1)" |\
        ccdfits - $out.den.fits  \
        object=$in comment="$comment" $wcspars 
    ccdmath $tmp.2.s,$tmp.0.s,$tmp.vel - "sqrt(%1/%2-%3*%3)" |\
        ccdfits - $out.sig.fits \
        object=$in comment="$comment" $wcspars 
    ccdfits $tmp.vel $out.vel.fits \
        object=$in comment="$comment" $wcspars 
    if ($?ring) then
     if (-e $ring) then
      fitsccd $refmap - | ccdmath -,$tmp.vel,$ring $tmp.diff "ifeq(%1,0,0,ifeq(%2,0,0,%1-%2))*%3/$sigma"
     else
      fitsccd $refmap - | ccdmath -,$tmp.vel $tmp.diff "ifeq(%1,0,0,ifeq(%2,0,0,%1-%2))/$sigma"
     endif
    endif
    ccdstat $tmp.diff bad=0 npar=$npar nppb=$nppb
    # ccdstat $tmp.diff bad=0 win=$weight
    # ccdstat $tmp.diff
    # ccdmath $tmp.diff - 'ifeq(%1,0,-99999,%1)' | ccdstat -  min=-1000 max=1000
    ccdfits $tmp.diff $out.diff.fits \
	object=$in comment="$comment" $wcspars 
    

    if ($clean) rm -fr $tmp.*
else
    echo Warning: skipping gridding and projecting
endif

exit 0

cleanup:
    echo Some error occured
    if ($clean) rm $tmp.*



