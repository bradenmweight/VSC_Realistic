%chk=geometry.chk
%mem=56GB
%nprocshared=24

#P FREQ=SaveNormalModes wB97XD/6-31G* NoSymm

TS Geometry (OPTIMIZED FROM QST3)

0 1 
C         -0.71456       -1.26431        0.12826
C          0.66184       -1.09565        0.03460
C          1.18972        0.18535       -0.14600
C          0.33995        1.28797       -0.23176
C         -1.03494        1.10698       -0.13673
C         -1.56752       -0.16689        0.04333
H         -1.12050       -2.26145        0.26863
H          1.32877       -1.95030        0.10057
H          0.76892        2.27428       -0.37206
H         -1.69270        1.96804       -0.20397
H         -2.64159       -0.30435        0.11714
N          2.56698        0.40598       -0.24594
C          3.59898       -0.21758       -0.22655
O          4.66564       -0.70552       -0.22528
H         34.94801761    1.09688745   -0.32180716
O         36.03482910    0.30720779   -0.31083935
C         36.98764417    0.35903847    0.75422263
C         36.95113534    1.41367935    1.67635151
C         37.95791177   -0.64480835    0.87585448
C         37.88489514    1.46447429    2.72011129
H         36.21031489    2.18013901    1.58348252
C         38.89167227   -0.59401281    1.91961360
H         37.98578689   -1.45005022    0.17178879
C         38.85516495    0.46062936    2.84174108
H         37.85702092    2.26971693    3.42417614
H         39.63249226   -1.36047286    2.01248303
H         39.56811144    0.49941256    3.63867468






