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
g2_2 := expand(simplify((4*X+2)^4*subs(X = -1/(4*X+2),g2_1)));
g3_2 := expand(simplify((4*X+2)^6*subs(X = -1/(4*X+2),g3_1)));
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
R1 := [seq([seq([solve(0 = subs(X = q(i,k),C1))],k = 1 .. N-1)],i = 1 .. 5)]:
display(Array(1 .. 5,[seq(pointplot([seq(seq([q(i,k),
    Re(R1[i][k][j])],k = 1 .. N-1),j = 1 .. 3)],color = ifelse(Im(R1[i][floor(.5*N)][2]) = 0,black,red)),i = 1 .. 5)]),size = [1000,
    1200]);
R2 := [seq([seq([solve(0 = subs(X = q(i,k),C2))],k = 1 .. N-1)],i = 1 .. 5)]:
display(Array(1 .. 5,[seq(pointplot([seq(seq([q(i,k),
    Re(R2[i][k][j])],k = 1 .. N-1),j = 1 .. 3)],color = ifelse(Im(R2[i][floor(.5*N)][2]) = 0,black,red)),i = 1 .. 5)]),size = [1000,
    1200]);
"Setting up the differential equation and checking its points of degeneracy:";
with(DETools):
equ := (1900544*X^11-7510016*X^10-16322048*X^9+20957640*X^8+116533604*X^7+180292814*X^6+141127495*X^5+60357456*X^4+13420080*X^3+1213056*X^2)*D^4+(28508160*X^10-87613440*X^9-216697600*X^8+75477472*X^7+998525016*X^6+1586638660*X^5+1153968680*X^4+434547288*X^3+81782784*X^2+6065280*X)*D^3+(121634816*X^9-260001792*X^8-812235520*X^7-481784560*X^6+1966231464*X^5+3546231148*X^4+2360301926*X^3+740670984*X^2+105528096*X+4852224)*D^2+(159645696*X^8-177414144*X^7-904080640*X^6-1391123712*X^5+549764200*X^4+1971179084*X^3+1200047188*X^2+288918648*X+24706944)*D+45613056*X^7+3661824*X^6-192325632*X^5-518676288*X^4-200519304*X^3+117866244*X^2+80286804*X+12633408;
sort(identify(evalf([solve(lcoeff(equ,D) = 0)],8)),(x, y) -> is(Re(x) < Re(y)));
"Getting expansions around A for the four products of periods (by solving the differential equation):";
SerSol := nformal_sol(equ,XA,500):
Solutions[1] := convert(SerSol[1],polynom):
Solutions[2] := convert(SerSol[2],polynom):
Solutions[3] := convert(SerSol[3],polynom):
Solutions[4] := convert(SerSol[4],polynom):
txA := `+`~(XA+`/`(1,64),`/`(txOFF,64));
PerProdA := ExplEvalAndComb(C1,C2,txA,Solutions,1,1):
IntA[1] := subs(Z = X-XA,int(subs(X = Z+XA,PerProdA[1]),Z)):
IntA[2] := subs(Z = X-XA,int(subs(X = Z+XA,PerProdA[2]),Z)):
IntA[3] := subs(Z = X-XA,int(subs(X = Z+XA,PerProdA[3]),Z)):
IntA[4] := subs(Z = X-XA,int(subs(X = Z+XA,PerProdA[4]),Z)):
"=======================================";
"===== Main calculations for the segment AB: =====";
"=======================================";
OutputCoincidences(PerProdA,XA):
# Identifying the coincidences at point -0.5625:
# The 1st curve doesn't have a singularity at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
ResAB := CalcToPoint4(equ,IntA,[XA, -`/`(17,32), XB],350):
# The integral (between points -0.5625 and -0.5) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -0.5625 and -0.5) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points -0.5625 and -0.5) of the product of period 2, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -0.5625 and -0.5) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResAB[1]);
OutputCoincidences(ResAB[3],XB):
# Identifying the coincidences at point -0.5:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
"Getting expansions around B for the four products of periods (by computing a transition matrix):";
txB := `+`~(XB+`/`(1,64),`/`(txOFF,16));
MatB := Matrix(4,4,ExplEvalAndComb(C1,C2,txB,ResAB[3],0,1)):
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
# Identifying the coincidences at point -0.5:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
#ResBC&#8788; CalcToPoint4:
ResBC := CalcToPoint4(equ,IntB,[XB, -`/`(15,32), -`/`(7,16), -`/`(13,32), -`/`(5,16), -`/`(1,4), XC],350):
# The integral (between points -0.5 and 0) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -0.5 and 0) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points -0.5 and 0) of the product of period 2, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -0.5 and 0) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResBC[1]);
OutputCoincidences(ResBC[3],XC):
# Identifying the coincidences at point 0:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
"Getting expansions around C for the four products of periods (by computing a transition matrix):";
txC := `+`~(XC+`/`(1,64),`/`(txOFF,16));
MatC := Matrix(4,4,ExplEvalAndComb(C1,C2,txC,ResBC[3],0,1)):
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
# Identifying the coincidences at point 0:
# Zeroes 3 and 1 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
ResCD := CalcToPoint4(equ,IntC,[XC, `/`(1,4), `/`(1,2), `/`(3,4), `/`(3,2), 1, 2, 2, XD],350):
# The integral (between points 0 and 4) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0 and 4) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 0 and 4) of the product of period 2, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0 and 4) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResCD[1]);
OutputCoincidences(ResCD[3],XD):
# Identifying the coincidences at point 4:
# Zeroes 3 and 1 of the 1st curve coincide at this point.
# The 2nd curve doesn't have a singularity at this point.
"Computing a transition matrix at D (for future reference -- checking solutions):";
txD := `+`~(XD+`/`(1,16),`/`(txOFF,16));
MatD := Matrix(4,4,ExplEvalAndComb(C1,C2,txD,ResCD[3],0,1)):
evalf(MatD,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatD := MTM:-round(MatD);
