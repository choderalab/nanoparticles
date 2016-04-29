#!/bin/bash

cat << EOF > leap.in
source leaprc.ff99SB
source leaprc.gaff
loadamberprep GTX.prepi
loadamberparams GTX.frcmod
x=loadpdb 1k3y_wats.pdb
saveamberparm x prmtop prmcrd
quit
EOF

tleap -f leap.in > leap.out

cat << EOF > min.in

&cntrl
irest=0,ntx=1,
imin=1,maxcyc=100,drms=0.0001,ntmin=2,
ntc=1,ntf=1,
cut=20.0,
ntpr=100,ntwx=0,ntwv=0,ntwe=0,
ipol=0,igb=0,ntb=0,
&end
EOF
sander -O -i min.in -o min.out -p prmtop -c prmcrd -r mincrd

cat << EOF > ptraj.in
trajin mincrd
trajout mincrd.pdb pdb
EOF
ptraj prmtop ptraj.in >& ptraj.out

mv mincrd.pdb.1 prot_xwat_gtx_min.pdb

echo "HEADER GTX" > gtx_bnd.pdb
grep GTX prot_xwat_gtx_min.pdb >> gtx_bnd.pdb

grep WAT prot_xwat_gtx.pdb > xwat_min.pdb
