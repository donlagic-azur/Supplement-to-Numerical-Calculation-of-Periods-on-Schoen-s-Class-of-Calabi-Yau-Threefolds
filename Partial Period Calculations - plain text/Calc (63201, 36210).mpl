read "ECHeader.mpl":
read "DEHeader.mpl":
read "PSHeader.mpl":
with(ECHeader);
with(DEHeader);
with(PSHeader);
"DOCUMENT START";
"Inputting data about the surfaces:";
Digits := 1000;
Y := 1;
g2 := `*`(12,X^4-4*X^3*Y+2*X*Y^3+Y^4);
g3 := `*`(4,2*X^6-12*X^5*Y+12*X^4*Y^2+14*X^3*Y^3+3*X^2*Y^4+6*X*Y^5+2*Y^6);
g2_1 := expand(X^4*subs(X = -1/X,g2));
g3_1 := expand(X^6*subs(X = -1/X,g3));
C1 := 4*x^3-g2_1*x-g3_1;
tempS := solve(1/J(C1) = 0):
tempS;
S := {tempS}:
g2_2 := expand(X^4*subs(X = -1/X,X^4*subs(X = `/`(1,4*X),g2)));
g3_2 := expand(X^6*subs(X = -1/X,X^6*subs(X = `/`(1,4*X),g3)));
C2 := 4*x^3-g2_2*x-g3_2;
tempS := solve(1/J(C2) = 0):
tempS;
S := S union {tempS}:
S := sort([seq(s,s in S)]):
print('[C, D, E, A]' = S);
YC := S[1];
YD := S[2];
YE := S[3];
YA := S[4];
txOFF := [`/`(1,16), `/`(3,64), `/`(5,64), `/`(1,32)];
"Plots of zeroes between singularities:";
with(plots):
p[1] := S[1]:
p[2] := S[2]:
p[3] := S[3]:
p[4] := S[4]:
p[0] := p[1]-`/`(1,2):
p[5] := p[4]+`/`(1,2):
N := 30:
q := ptmp:
R1 := [seq([seq([solve(0 = subs(X = q(i,k),C1))],k = 1 .. N-1)],i = 1 .. 5)]:
display(Array(1 .. 5,[seq(pointplot([seq(seq([q(i,k),
    Re(R1[i][k][j])],k = 1 .. N-1),j = 1 .. 3)],color = ifelse(Im(R1[i][floor(.5*N)][2]) = 0,black,red)),i = 1 .. 5)]),size = [1000,
    1200]);
R2 := [seq([seq([solve(0 = subs(X = q(i,k),C2))],k = 1 .. N-1)],i = 1 .. 5)]:
display(Array(1 .. 5,[seq(pointplot([seq(seq([q(i,k),
    Re(R2[i][k][j])],k = 1 .. N-1),j = 1 .. 3)],color = ifelse(Im(R2[i][floor(.5*N)][2]) = 0,black,red)),i = 1 .. 5)]),size = [1000,
    1200]);
"Setting up the differential equation and checking its points of degeneracy:";
with(DETools);
equ := (128*X^11+2960*X^10-2288*X^9-244699*X^8+472822*X^7+945644*X^6-1957592*X^5-73216*X^4+378880*X^3+65536*X^2)*D^4+(1920*X^10+36672*X^9-85816*X^8-2721556*X^7+3603204*X^6+11706472*X^5-17379392*X^4+1281792*X^3+2883584*X^2+327680*X)*D^3+(8192*X^9+121568*X^8-577976*X^7-8304502*X^6+5245660*X^5+37491704*X^4-40077680*X^3+5672192*X^2+4698112*X+262144)*D^2+(10752*X^8+109760*X^7-1011592*X^6-7447296*X^5-524900*X^4+31521592*X^3-25051680*X^2+2687744*X+1261568)*D+3072*X^7+15040*X^6-345408*X^5-1371960*X^4-680812*X^3+4209384*X^2-2807328*X-161792;
sort([solve(lcoeff(equ,D) = 0)],(x, y) -> is(x < y));
evalf(`/`(59,16)-`/`(3*sqrt(273),16),5);
evalf(`/`(59,16)+`/`(3*sqrt(273),16),5);
"Getting expansions around C for the four products of periods (by solving the differential equation):";
SerSol := nformal_sol(equ,YC,500):
Solutions[1] := convert(SerSol[1],polynom):
Solutions[2] := convert(SerSol[2],polynom):
Solutions[3] := convert(SerSol[3],polynom):
Solutions[4] := convert(SerSol[4],polynom):
txC := `+`~(YC,txOFF);
PerProdC := ExplEvalAndComb(C1,C2,txC,Solutions,1,1):
IntC[1] := subs(Z = X-YC,int(subs(X = Z+YC,PerProdC[1]),Z)):
IntC[2] := subs(Z = X-YC,int(subs(X = Z+YC,PerProdC[2]),Z)):
IntC[3] := subs(Z = X-YC,int(subs(X = Z+YC,PerProdC[3]),Z)):
IntC[4] := subs(Z = X-YC,int(subs(X = Z+YC,PerProdC[4]),Z)):
"=======================================";
"===== Main calculations for the segment CD: =====";
"=======================================";
OutputCoincidences(PerProdC,YC):
# Identifying the coincidences at point -16:
# The 1st curve doesn't have a singularity at this point.
# Zeroes 3 and 1 of the 2nd curve coincide at this point.
ResCD := CalcToPoint4(equ,IntC,[YC, -10, -8, -5, -4, `+`(-3,I), -2+`*`(2,I), `/`(-6+`*`(6,I),5),
    `/`(-4+`*`(4,I),5), `/`(-7+`*`(6,I),15), `/`(-6+`*`(4,I),15), `/`(-5+`*`(2,I),15), YD],350):
# The integral (between points -16 and -0.25) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -16 and -0.25) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points -16 and -0.25) of the product of period 2, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -16 and -0.25) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResCD[1]);
OutputCoincidences(ResCD[3],YD):
# Identifying the coincidences at point -0.25:
# Zeroes 3 and 1 of the 1st curve coincide at this point.
# The 2nd curve doesn't have a singularity at this point.
"Getting expansions around D for the four products of periods (by computing a transition matrix):";
txD := `+`~(YD+`/`(1,16),`/`(txOFF,16));
MatD := Matrix(4,4,ExplEvalAndComb(C1,C2,txD,ResCD[3],0,1)):
evalf(MatD,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatD := MTM:-round(MatD);
with(LinearAlgebra):
PerProdD := MatrixVectorMultiply(MatD,ResCD[3]):
IntD := MatrixVectorMultiply(MatD,ResCD[2]-ResCD[1]):
"=======================================";
"===== Main calculations for the segment DE: =====";
"=======================================";
OutputCoincidences(PerProdD,YD):
# Identifying the coincidences at point -0.25:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# The 2nd curve doesn't have a singularity at this point.
ResDE := CalcToPoint4(equ,IntD,[YD, -`/`(1,8), YE],350):
# The integral (between points -0.25 and 0) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -0.25 and 0) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points -0.25 and 0) of the product of period 2, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -0.25 and 0) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResDE[1]);
OutputCoincidences(ResDE[3],YE):
# Identifying the coincidences at point 0:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 3 and 1 of the 2nd curve coincide at this point.
"Getting expansions around E for the four products of periods (by computing a transition matrix):";
txE := `+`~(YE,txOFF);
MatE := Matrix(4,4,ExplEvalAndComb(C1,C2,txE,ResDE[3],0,1)):
evalf(MatE,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatE := MTM:-round(MatE);
with(LinearAlgebra):
PerProdE := MatrixVectorMultiply(MatE,ResDE[3]):
IntE := MatrixVectorMultiply(MatE,ResDE[2]-ResDE[1]):
"=======================================";
"===== Main calculations for the segment EA: =====";
"=======================================";
OutputCoincidences(PerProdE,YE):
# Identifying the coincidences at point 0:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
ResEA := CalcToPoint41(equ,IntE,[YE, `/`(1+I,10), `/`(1+I,4), `/`(2+I,4), `/`(3+`*`(2,I),4),
    `/`(7+`*`(2,I),8), 1, `/`(6,5), YA],350):
# The integral (between points 0 and 2) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0 and 2) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 0 and 2) of the product of period 2, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0 and 2) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResEA[1]);
OutputCoincidences(ResEA[3],YA):
# Identifying the coincidences at point 2:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
"Computing a transition matrix at A (for future reference -- checking solutions):";
txA := `+`~(YA,txOFF);
MatA := Matrix(4,4,ExplEvalAndComb(C1,C2,txA,ResEA[3],0,1)):
evalf(MatA,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatA := MTM:-round(MatA);
#checking if the product (from E to A) of the transition matrices we got is Identity - this reflects the correctness of the solutions, as the matrices were not only instrumental to the calculation of the results, but were also calculated using the information (analytically continued solutions) available at that moment
(((Matrix(4,4,[[1, 1, 0, 0], [-1, -2, 0, 0], [-3, -3, 1, 1], [3, 6, -1, -2]]) . Matrix(4,4,[[1, 0, 0, 0], [0,
    1, 0, 0], [1, 0, -1, 0], [0, 1, 0, -1]])) . Matrix(4,4,[[1, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 1], [0, 0, 1,
    0]])) . Matrix(4,4,[[-1, 0, 1, 0], [3, -1, -3, 1], [2, 0, -1, 0], [-6, 2, 3, -1]])) . Matrix(4,4,[[1, 1, 1,
    1], [0, 1, 0, 1], [0, 0, 1, 1], [0, 0, 0, 1]]);
