#! /bin/csh -f
#
#  This scripts takes an HDF output snapshot file from cmhog
#  (bar hydro, polar coordinates), projects it to a requested
#  sky view, and using WIP summarizes the results
#  Bar is conveniently oriented along Y axis (PA=0), and flows CCW.
#  For CW rotating galaxies you may have to add 180 to 'inc' and/or 'pa'
#
#
#  23-sep-95	Created				Peter Teuben
#  15-nov-95    beam=0.2 frang=45
#  24-nov-95    defaults for more central region, fixed dependancies
#   9-jul-96    hacking for N5383 
#  24-jul-00    BIMA proposal N4303 et al.
#   5-sep-00    modified to write cubes instead of moment maps
#  12-mar-01    radecvel=t to make karma swallow these fits files
#  13-mar-01    use phi,inc,pa (no more +/- 180) and documented geometry
#               (notice that earlier versions had sign of radial vel wrong)
#  23-mar-01    added a refmap and fixed refscale; this assumes that the refmap
#		(often a cube) has the reference pixel defined to the be 
#		center of the galaxy (bimasong data often don't do the VELO axis correct)
#  17-apr-01    added velocity referencing using 'vsys' to be at v=0 in the model cube
#               (bugs when model and data have different delta-V)
#  11-apr-02    generic version with new geometry definitions

if ($#argv == 0) then
  echo Usage: $0 in=HDF_DATASET out=BASENAME ...
  echo Gridding and projecting 2D CMHOG hydro models to given bar viewing angles
  exit 0
endif

# 			Required Keywords
unset in
unset out
# 			Defaulted Keywords
set pa=-42
set inc=27
set phi=44
set rot=1
#
set rmax=6
set n=200
set beam=0.05
set color=1
set clean=1
#
set refmap=""
set pscale=0.5
set vscale=1
set vsys=0
#			Velocity gridding for cube  (dv=2*vmax/nvel)
set nvel=50
set vmax=250

set par=""

#			Parse commandline to (re)set keywords, special case for par=
foreach a ($*)
  set $a
  if (-e "$par") then
    source $par
    set par=""
  else
    echo Warning, par=$par does not exist
    set par=""
  endif
end

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

#			Show i'm happy
echo Files: in=$in out=$out 
echo Grid: rmax=$rmax n=$n beam=$beam
echo Projection: phi=$phi inc=$inc pa=$pa
#                       Derived quantities
set cell=`nemoinp "2*$rmax/$n"`
set range="-${rmax}:${rmax}"
echo "      Derived: cell=$cell"

if (-e "$refmap") then
    #   referencing 
    set nz=(`fitshead $refmap | grep ^NAXIS3 | awk '{print $3}'`)
    set pz=(`fitshead $refmap | grep ^CRPIX3 | awk '{print $3}'`)
    set vz=(`fitshead $refmap | grep ^CRVAL3 | awk '{print $3}'`)
    set dz=(`fitshead $refmap | grep ^CDELT3 | awk '{print $3}'`)
    
    set dz1=`nemoinp "2000*$vmax/$nz"`
    set vref=`nemoinp "($vz-1000*$vsys)/($dz1)+$nvel/2+0.5"`
    #set vscale=`nemoinp "$vscale*(2*$vmax/$nvel)/($dz/1000)"`
    set refscale=$pscale,$pscale,$vscale
    #                                    CHECK : is this -0.5 or +0.5   ?????
    set refcen=`nemoinp $n/2-0.5`
    #set refpix=$refcen,$refcen,$vref

    
    #   now assuming model is centered, as well as data cube
    set vref=`nemoinp $nvel/2+0.5`
    ###set vref=`nemoinp $nvel/2-0.5`
    set refpix=$refcen,$refcen,$vref
    
    
    echo $nz $pz $vz $dz 
    echo Vsys at OBS pixel: `nemoinp "(1000*$vsys-$vz)/$dz+$pz"` 
    echo REFPIX:   $refpix
    echo REFSCALE: $refscale
else
    echo BUG: need to rewrite this section for when no refmap given..... since you did not
    exit 1
endif    

#> nemo.need tabtos snaptrans snaprotate snapadd snapgrid ccdsmooth ccdmath ccdfits fitshead

set tmp=tmp$$
if (! -e $out.den.fits) then

    # convert the half-plane HDF file to a full plane snapshot file

    tsd in=$in out=$tmp.tab coord=t
    if ($status) goto cleanup
    tabtos in=$tmp.tab out=$tmp.s0 block1=x,y,vx,vy,mass
    snaptrans in=$tmp.s0 out=$tmp.s1 ctypei=cyl ctypeo=cart
    snaprotate in=$tmp.s1 out=$tmp.s2 theta=180 order=z
    snapadd $tmp.s1,$tmp.s2 $tmp.s3

    # project for skyview, and create a intensity and velocity field

    snaprotate $tmp.s3 $tmp.snap \
        "atand(tand($phi)/cosd($inc)),$inc,$pa" zyz

    foreach mom (0 1 2)
         snapgrid in=$tmp.snap out=$tmp.$mom \
                xrange=$range yrange=$range nx=$n ny=$n moment=$mom mean=t
         ccdsmooth in=$tmp.$mom out=$tmp.$mom.s gauss=$beam
    end
    ccdmath $tmp.1.s,$tmp.0.s $tmp.vel %1/%2
## BUG: ifgt() doesn't work
##    ccdmath $tmp.0.s - "ifgt(%1,0,log(%1),-10)" | ccdfits - $out.den.fits
    ccdmath $tmp.0.s - "log(%1)" | ccdmath - - "ifeq(%1,0,-10,%1)" | ccdfits - $out.den.fits  \
        object=$in comment="$0 $in $out $pa,$inc,$phi,$range,$n,$beam"  \
	refmap=$refmap scale=$refscale refpix=$refpix
    ccdmath $tmp.2.s,$tmp.0.s,$tmp.vel - "sqrt(%1/%2-%3*%3)" | ccdfits - $out.sig.fits \
        object=$in comment="$0 $in $out $pa,$inc,$phi,$range,$n,$beam" 	\
	refmap=$refmap scale=$refscale refpix=$refpix
    ccdfits $tmp.vel $out.vel.fits \
        object=$in comment="$0 $in $out $pa,$inc,$phi,$range,$n,$beam" 	\
	refmap=$refmap scale=$refscale refpix=$refpix


### BUG:      ccdmath - - "ifgt(%1,0,log(%1),-10)" |\


    # now also create the (smoothed) cube
    snapgrid in=$tmp.snap out=- \
          xrange=$range yrange=$range zrange=-${vmax}:${vmax} \
	  xvar=x yvar=y zvar=-vz \
	  nx=$n ny=$n nz=$nvel moment=0 mean=t |\
      ccdsmooth - - gauss=$beam |\
      ccdmath - - "log(%1)" |\
      ccdmath - - "ifeq(%1,0,-10,%1)" |\
      ccdfits - $out.cube.fits \
        object=$in comment="$0 $in $out $pa,$inc,$phi,$range,$n,$beam" \
	refmap=$refmap scale=$refscale refpix=$refpix

    rm -fr $tmp.*
else
    echo Warning: skipping gridding and projecting
endif

exit 0

cleanup:
    if ($clean) rm $tmp.*



