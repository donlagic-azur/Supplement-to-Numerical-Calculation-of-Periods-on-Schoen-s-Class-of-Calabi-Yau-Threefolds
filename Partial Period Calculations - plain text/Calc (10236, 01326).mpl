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
g2_2 := expand(subs(X = -`/`(1,2)-X,g2_1));
g3_2 := expand(subs(X = -`/`(1,2)-X,g3_1));
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
equ := (896*X^11+2464*X^10-6416*X^9-19056*X^8-557008*X^7-951454*X^6+6856711*X^5+9162576*X^4+3715632*X^3+466560*X^2)*D^4+(13440*X^10+33600*X^9+71360*X^8+92320*X^7-8688024*X^6-13103996*X^5+60044744*X^4+65510424*X^3+21316608*X^2+2332800*X)*D^3+(57344*X^9+129024*X^8+988352*X^7+1579088*X^6-35509800*X^5-45355364*X^4+122801126*X^3+103487688*X^2+25093152*X+1866240)*D^2+(75264*X^8+150528*X^7+2196896*X^6+3163632*X^5-41373080*X^4-42680284*X^3+47470612*X^2+29102904*X+5298048)*D+21504*X^7+37632*X^6+914880*X^5+1120080*X^4-9783456*X^3-7616436*X^2-1921500*X-160704;
solve(factor(lcoeff(equ,D)) = 0);
evalf(-`/`(1,4)+`/`(`*`(13,I)*sqrt(119),28),5);
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
# Identifying the coincidences at point -4.5:
# The 1st curve doesn't have a singularity at this point.
# Zeroes 3 and 1 of the 2nd curve coincide at this point.
ResAB := CalcToPoint4(equ,IntA,[XA, -`/`(5,2), -2, -`/`(6,5), -1, -`/`(4,5), -`/`(3,4), -`/`(3,5), XB],350):
# The integral (between points -4.5 and -0.5) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -4.5 and -0.5) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points -4.5 and -0.5) of the product of period 2, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -4.5 and -0.5) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResAB[1]);
OutputCoincidences(ResAB[3],XB):
# Identifying the coincidences at point -0.5:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 3 and 1 of the 2nd curve coincide at this point.
"Getting expansions around B for the four products of periods (by computing a transition matrix):";
txB := `+`~(XB,txOFF);
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
ResBC := CalcToPoint4(equ,IntB,[XB, `/`(-3+I,8), `/`(-1+I,4), `/`(-1+I,8), XC],350):
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
txC := `+`~(XC,txOFF);
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
ResCD := CalcToPoint4(equ,IntC,[XC, `/`(1,8), `/`(1,4), `/`(1,4), `/`(1,2), `/`(1,2), 1, 1, 2, 2, XD],350):
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
txD := `+`~(XD,txOFF);
MatD := Matrix(4,4,ExplEvalAndComb(C1,C2,txD,ResCD[3],0,1)):
evalf(MatD,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatD := MTM:-round(MatD);
ResAB[3]-subs(log(X-XB) = log(X-XB)+`*`(2*pi,I),ResAB[3]);
evalf(ResAB[3],5);
MatB2 := Matrix(4,4,ExplEvalAndComb(C1,C2,txB,subs(log(X-XB) = log(X-XB)+evalf(`*`(2*Pi,I)),ResAB[3]),0,1)):
#MatB2 ~Matrix:
evalf(MatB2,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatB2 := MTM:-round(MatB2);
MatB;
1/Matrix(4,4,[[-1, 1, 1, -1], [-2, 1, 2, -1], [0, 0, -1, 1], [0, 0, -2, 1]]) . Matrix(4,4,[[-1, 1, -1, 1], [1,
    -2, 1, -2], [0, 0, -1, 1], [0, 0, 1, -2]]);
ExplEvalAndComb1 := proc(C1,C2,tx,ToBeCombined, MatrixIfZeroValuesIfOne,debug := 0, weak := 0,
    tofs := [-sqrt(I),-sqrt(I),-sqrt(I),-sqrt(I)],conj := [0,0],Me := [6,4])
    local b,R1,R2,I1_1,I1_2,I2_1,I2_2,J11,J12,J21,J22,Res;
    for b from 1 to 4 do
        R1[b] := sort([solve(0 = subs(X = tx[b],C1))],(x, y) -> is(Im(x) < Im(y)) or is(x < y));
        R2[b] := sort([solve(0 = subs(X = tx[b],C2))],(x, y) -> is(Im(x) < Im(y)) or is(x < y));
    end do:
    if debug=1 then
        print("Roots of first curve at the chosen points"= Vector(4, R1));
    end if:
    if debug=1 then
        print("Roots of second curve at the chosen points"= Vector(4, R2));
    end if:
    #Ia_c[b] stands for a-th curve, b-th sample point, c-th principal integral
    #Jcd[b] stands for b-th DE sample point, c-th principal integral on 1 st curve, d-principal integral on 2 nd curve
    if weak = 1 then
        for b from 1 to 4 do
            I1_1[b] := 2*WeakEIntegral(subs(X = tx[b],C1),R1[b][1],R1[b][2],tofs[1],Me[1],Me[2]):
            if debug =2 then
                print(evalf(I1_1[b],15)):
            end if:
            I1_2[b] := 2*WeakEIntegral(subs(X = tx[b],C1),R1[b][2],R1[b][3],tofs[2],Me[1],Me[2]):
            if debug =2 then
                print(evalf(I1_2[b],15)):
            end if:
            I2_1[b] := 2*WeakEIntegral(subs(X = tx[b],C2),R2[b][1],R2[b][2],tofs[3],Me[1],Me[2]):
            if debug =2 then
                print(evalf(I2_1[b],15)):
            end if:
            I2_2[b] := 2*WeakEIntegral(subs(X = tx[b],C2),R2[b][2],R2[b][3],tofs[4],Me[1],Me[2]):
            if debug =2 then
                print(evalf(I2_2[b],15)):
            end if:
            J11[b] := evalf(I1_1[b]*I2_1[b]):
            J12[b] := evalf(I1_1[b]*I2_2[b]):
            J21[b] := evalf(I1_2[b]*I2_1[b]):
            J22[b] := evalf(I1_2[b]*I2_2[b]):
        end do:
    else
        for b from 1 to 4 do
            I1_1[b] := 2*evalf(int(sqrt(tofs[1])/sqrt(tofs[1]*subs(X = tx[b],C1)),x = R1[b][1] .. R1[b][2])):
            if debug =2 then
                print(evalf(I1_1[b],15)):
            end if:
            I1_2[b] := 2*evalf(int(sqrt(tofs[2])/sqrt(tofs[2]*subs(X = tx[b],C1)),x = R1[b][2] .. R1[b][3])):
            if debug =2 then
                print(evalf(I1_2[b],15)):
            end if:
            I2_1[b] := 2*evalf(int(sqrt(tofs[3])/sqrt(tofs[3]*subs(X = tx[b],C2)),x = R2[b][1] .. R2[b][2])):
            if debug =2 then
                print(evalf(I2_1[b],15)):
            end if:
            I2_2[b] := 2*evalf(int(sqrt(tofs[4])/sqrt(tofs[4]*subs(X = tx[b],C2)),x = R2[b][2] .. R2[b][3])):
            if debug =2 then
                print(evalf(I2_2[b],15)):
            end if:
            if abs(conj[1])=1 then
                I1_1[b] := conjugate(I1_2[b])*sign(conj[1]):
            end if:
            if abs(conj[1])=2 then
                I1_2[b] := conjugate(I1_1[b])*sign(conj[1]):
            end if:
            if abs(conj[2])=1 then
                I2_1[b] := conjugate(I2_2[b])*sign(conj[2]):
            end if:
            if abs(conj[2])=2 then
                I2_2[b] := conjugate(I2_1[b])*sign(conj[2]):
            end if:
            J11[b] := evalf(I1_1[b]*I2_1[b]):
            J12[b] := evalf(I1_1[b]*I2_2[b]):
            J21[b] := evalf(I1_2[b]*I2_1[b]):
            J22[b] := evalf(I1_2[b]*I2_2[b]):
        end do:
    end if:
    Res[1] := PolynomialCombination(ToBeCombined,4,(-`*`(.1e-4,I)) +~ evalf(tx),J11,MatrixIfZeroValuesIfOne):
    Res[2] := PolynomialCombination(ToBeCombined,4,(-`*`(.1e-4,I)) +~ evalf(tx),J12,MatrixIfZeroValuesIfOne):
    Res[3] := PolynomialCombination(ToBeCombined,4,(-`*`(.1e-4,I)) +~ evalf(tx),J21,MatrixIfZeroValuesIfOne):
    Res[4] := PolynomialCombination(ToBeCombined,4,(-`*`(.1e-4,I)) +~ evalf(tx),J22,MatrixIfZeroValuesIfOne):
    Res := [Res[1], Res[2], Res[3], Res[4]]:
    #if  MatrixIfZeroValuesIfOne = 0 then return Matrix:
    ##if  MatrixIfZeroValuesIfOne = 1 then return Vector:
end proc:
series(Period(C1)*Period(C2),X,10);
2*Pi*(12*G2(f))^(-`/`(1,4))*hypergeom([`/`(1,12), `/`(5,12)],[1],1/J(f));
con := series(G2(C1)^(-`/`(1,4))*hypergeom([`/`(1,12),
    `/`(5,12)],[1],1/J(C1))*G2(C2)^(-`/`(1,4))*hypergeom([`/`(1,12), `/`(5,12)],[1],1/J(C2)),X,10);
expand(108/(12^`/`(3,4))/(3^`/`(3,4))/(4^`/`(1,4))/(9^`/`(3,4))*convert(con,polynom));
evalf(12^`/`(3,4)*3^`/`(3,4)*4^`/`(1,4)*9^`/`(3,4));
