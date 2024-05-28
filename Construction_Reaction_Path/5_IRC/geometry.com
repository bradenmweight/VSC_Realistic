%chk=geometry.chk
%oldchk=../4_NM_at_TS/geometry.chk
%mem=56GB
%nprocshared=24

#P IRC=(DVV,RCFC,Phase=(16,14),MaxPoints=1000) wB97XD/6-31G* NoSymm
#P guess=read geom=check

Autonmatically read TS geometry from Normal Mode Calculation

0 1








 