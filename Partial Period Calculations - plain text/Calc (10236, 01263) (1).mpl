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
g2_1 := `*`(12,X^4-4*X^3*Y+2*X*Y^3+Y^4);
g3_1 := `*`(4,2*X^6-12*X^5*Y+12*X^4*Y^2+14*X^3*Y^3+3*X^2*Y^4+6*X*Y^5+2*Y^6);
C1 := 4*x^3-g2_1*x-g3_1;
tempS := solve(1/J(C1) = 0):
tempS;
S := {tempS}:
g2_2 := expand(X^4*subs(X = `/`(1,4*X),g2_1));
g3_2 := expand(X^6*subs(X = `/`(1,4*X),g3_1));
C2 := 4*x^3-g2_2*x-g3_2;
tempS := solve(1/J(C2) = 0):
tempS;
S := S union {tempS}:
S := sort([seq(s,s in S)]):
print('[A, B, C, D]' = S);
XA := S[1];
XB := S[2];
XC := S[3];
XD := S[4];
txOFF := [`/`(1,16), `/`(3,64), `/`(5,64), `/`(1,32)];
"Plots of zeroes between singularities: (-inf, A, B, C, D, inf)";
with(plots):
p[1] := S[1]:
p[2] := S[2]:
p[3] := S[3]:
p[4] := S[4]:
p[0] := p[1]-`/`(1,2):
p[5] := p[4]+`/`(1,2):
N := 30:
q := ptmp:
R1 := [seq([seq([fsolve(0 = subs(X = q(i,k),C1),complex)],k = 1 .. N-1)],i = 1 .. 5)]:
display(Array(1 .. 5,[seq(pointplot([seq(seq([q(i,k),
    Re(R1[i][k][j])],k = 1 .. N-1),j = 1 .. 3)],color = ifelse(Im(R1[i][floor(.5*N)][2]) = 0,black,red)),i = 1 .. 5)]),size = [1000,
    1200]);
R2 := [seq([seq([fsolve(0 = subs(X = q(i,k),C2),complex)],k = 1 .. N-1)],i = 1 .. 5)]:
display(Array(1 .. 5,[seq(pointplot([seq(seq([q(i,k),
    Re(R2[i][k][j])],k = 1 .. N-1),j = 1 .. 3)],color = ifelse(Im(R2[i][floor(.5*N)][2]) = 0,black,red)),i = 1 .. 5)]),size = [1000,
    1200]);
"Setting up the differential equation and checking its points of degeneracy:";
with(DETools):
equ := (65536*X^11-378880*X^10-73216*X^9+1957592*X^8+945644*X^7-472822*X^6-244699*X^5+2288*X^4+2960*X^3-128*X^2)*D^4+(983040*X^10-4694016*X^9-2746112*X^8+21772448*X^7+7206408*X^6-5853236*X^5-2172424*X^4-40056*X^3+22528*X^2-640*X)*D^3+(4194304*X^9-15560704*X^8-18495232*X^7+66436016*X^6+10491320*X^5-18745852*X^4-5009710*X^3-177256*X^2+36704*X-512)*D^2+(5505024*X^8-14049280*X^7-32370944*X^6+59578368*X^5-1049800*X^4-15760796*X^3-3131460*X^2-83992*X+9856)*D+1572864*X^7-1925120*X^6-11053056*X^5+10975680*X^4-1361624*X^3-2104692*X^2-350916*X+5056;
sort([solve(lcoeff(equ,D) = 0)],(x, y) -> is(x < y));
evalf(-`/`(59,64)+`/`(3*sqrt(273),64),5);
"Getting expansions around A for the four products of periods (by solving the differential equation):";
SerSol := nformal_sol(equ,XA,500):
Solutions[1] := convert(SerSol[1],polynom):
Solutions[2] := convert(SerSol[2],polynom):
Solutions[3] := convert(SerSol[3],polynom):
Solutions[4] := convert(SerSol[4],polynom):
txA := `+`~(XA,txOFF);
PerProdA := ExplEvalAndComb(C1,C2,txA,Solutions,1,2):
IntA[1] := subs(Z = X-XA,int(subs(X = Z+XA,PerProdA[1]),Z)):
IntA[2] := subs(Z = X-XA,int(subs(X = Z+XA,PerProdA[2]),Z)):
IntA[3] := subs(Z = X-XA,int(subs(X = Z+XA,PerProdA[3]),Z)):
IntA[4] := subs(Z = X-XA,int(subs(X = Z+XA,PerProdA[4]),Z)):
"=======================================";
"===== Main calculations for the segment AB: =====";
"=======================================";
OutputCoincidences(PerProdA,XA):
# Identifying the coincidences at point -0.5:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
ResAB := CalcToPoint4(equ,IntA,[XA, `/`(-2+I,6), `/`(-1+I,4), `/`(`*`(2,I)-1,10), `/`(I,5), `/`(I,12),
    `/`(I,16), `/`(I,32), XB],350):
# The integral (between points -0.5 and 0) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -0.5 and 0) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points -0.5 and 0) of the product of period 2, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -0.5 and 0) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResAB[1]);
OutputCoincidences(ResAB[3],XB):
# Identifying the coincidences at point 0:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
"Getting expansions around B for the four products of periods (by computing a transition matrix):";
txB := `+`~(XB+`/`(1,64),`/`(txOFF,64));
MatB := Matrix(4,4,ExplEvalAndComb(C1,C2,txB,ResAB[3],0,2,0,sqrt(-I) +~ [0, 0, 0, 0],[2, 0])):
evalf(MatB,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatB := MTM:-round(MatB);
with(LinearAlgebra):
PerProdB := MatrixVectorMultiply(MatB,ResAB[3]):
IntB := MatrixVectorMultiply(MatB,ResAB[2]-ResAB[1]):
"=======================================";
"===== Main calculations for the segment BC: =====";
"=======================================";
OutputCoincidences(PerProdB,XB):
# Identifying the coincidences at point 0:
# Zeroes 3 and 1 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
ResBC := CalcToPoint4(equ,IntB,[XB, `/`(1,32), XC],350):
# The integral (between points 0 and 0.0625) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0 and 0.0625) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 0 and 0.0625) of the product of period 2, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0 and 0.0625) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResBC[1]);
OutputCoincidences(ResBC[3],XC):
# Identifying the coincidences at point 0.0625:
# The 1st curve doesn't have a singularity at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
"Getting expansions around C for the four products of periods (by computing a transition matrix):";
txC := `+`~(XC+`/`(1,64),`/`(txOFF,16));
MatC := Matrix(4,4,ExplEvalAndComb(C1,C2,txC,ResBC[3],0,2,0,`+`~(-sqrt(I),[0, 0, 0, 0]),[0, -1])):
evalf(MatC,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatC := MTM:-round(MatC);
with(LinearAlgebra):
PerProdC := MatrixVectorMultiply(MatC,ResBC[3]):
IntC := MatrixVectorMultiply(MatC,ResBC[2]-ResBC[1]):
"=======================================";
"===== Main calculations for the segment CD: =====";
"=======================================";
OutputCoincidences(PerProdC,XC):
# Identifying the coincidences at point 0.0625:
# The 1st curve doesn't have a singularity at this point.
# Zeroes 3 and 1 of the 2nd curve coincide at this point.
ResCD := CalcToPoint4(equ,IntC,[XC, `/`(5+`*`(2,I),60), `/`(3+`*`(2,I),30), `/`(7+`*`(6,I),60), `/`(1+I,5),
    `/`(3+`*`(3,I),10), `/`(1+I,2), `/`(3+I,4), 1, `/`(6,5), 2, `/`(5,2), XD],350):
# The integral (between points 0.0625 and 4) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0.0625 and 4) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 0.0625 and 4) of the product of period 2, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0.0625 and 4) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResCD[1]);
OutputCoincidences(ResCD[3],XD):
# Identifying the coincidences at point 4:
# Zeroes 3 and 1 of the 1st curve coincide at this point.
# The 2nd curve doesn't have a singularity at this point.
"Computing a transition matrix at D (for future reference -- checking solutions):";
txD := `+`~(XD+`/`(1,16),`/`(txOFF,16));
MatD := Matrix(4,4,ExplEvalAndComb(C1,C2,txD,ResCD[3],0,2,0,[sqrt(-I), sqrt(I), sqrt(-I), sqrt(-I)],[0, -2])):
evalf(MatD,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatD := MTM:-round(MatD);
evalf(ResAB[1],40);
evalf(ResBC[1],40);
evalf(ResCD[1],40);
nops;
op(1 .. 5,ResAB[3][1]);
