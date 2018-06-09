#!/bin/tcsh

set bin = /data/Kepler/bin

if($#argv != 1) then
    goto error
endif

set kidlc = `echo $1 | awk '{printf "%08d",$1}'`
set kidname = `echo $1 | awk '{printf "%09d",$1}'`
set kidXXXX = `echo $1 | awk '{printf "%04d",$1/100000}'`

echo ftp://archive.stsci.edu/pub/kepler/lightcurves/$kidXXXX/$kidname/kplr"*"_llc.fits
wget --quiet ftp://archive.stsci.edu/pub/kepler/lightcurves/$kidXXXX/$kidname/kplr"*"_llc.fits

$bin/kfitsread kplr$kidname* | awk '{print $0,0}' | sort -bn > klc$kidlc.lc.dat
/bin/rm kplr$kidname*

exit(0)

error:
cat <<EOF
 Usage: $0 KID
  KID - Kepler ID (integer)

EOF
exit(-1)

