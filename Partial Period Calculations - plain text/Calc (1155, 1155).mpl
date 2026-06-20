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
g2_1 := `*`(3,X^4-12*X^3*Y+14*X^2*Y^2+12*X*Y^3+Y^4);
g3_1 := X^6-18*X^5*Y+75*X^4*Y^2+75*X^2*Y^4+18*X*Y^5+Y^6;
C1 := 4*x^3-g2_1*x-g3_1;
tempS := solve(1/J(C1) = 0):
tempS;
S := {tempS}:
g2_2 := g2_1;
g3_2 := g3_1;
C2 := C1;
tempS := solve(1/J(C2) = 0):
tempS;
S := S union {tempS}:
S := sort([seq(s,s in S)],(a, b) -> is(a < b)):
print('[A, B, C]' = S);
XA := S[1];
XB := S[2];
XC := S[3];
txOFF := [`/`(3,64), `/`(1,16), `/`(5,64)];
"Plots of zeroes between singularities:";
with(plots):
p[1] := S[1]:
p[2] := S[2]:
p[3] := S[3]:
p[0] := p[1]-`/`(1,2):
p[4] := p[3]+`/`(1,2):
N := 30:
q := ptmp:
R1 := [seq([seq([solve(0 = subs(X = q(i,k),C1))],k = 1 .. N-1)],i = 1 .. 4)]:
display(Array(1 .. 4,[seq(pointplot([seq(seq([q(i,k),
    Re(R1[i][k][j])],k = 1 .. N-1),j = 1 .. 3)],color = ifelse(Im(R1[i][floor(.5*N)][2]) = 0,black,red)),i = 1 .. 4)]),size = [1000,
    1200]);
"Setting up the differential equation and checking its points of degeneracy:";
with(DETools):
equ := (X^6-22*X^5+119*X^4+22*X^3+X^2)*D^3+(9*X^5-165*X^4+714*X^3+99*X^2+3*X)*D^2+(19*X^4-276*X^3+842*X^2+78*X+1)*D+8*X^3-84*X^2+128*X+6;
sort([solve(lcoeff(equ,D) = 0)],(x, y) -> is(Re(x) < Re(y)));
"Getting expansions around A for the four products of periods (by solving the differential equation):";
SerSol := nformal_sol(equ,XA,500):
Solutions[1] := convert(SerSol[1],polynom):
Solutions[2] := convert(SerSol[2],polynom):
Solutions[3] := convert(SerSol[3],polynom):
#txA ~;
txA := `+`~(XA,`/`(txOFF,2));
PerProdA := ExplEvalAndComb31(C1,C2,txA,Solutions,1,2,0,[1, 1]):
IntA[1] := subs(Z = X-XA,int(subs(X = Z+XA,PerProdA[1]),Z)):
IntA[2] := subs(Z = X-XA,int(subs(X = Z+XA,PerProdA[2]),Z)):
IntA[3] := subs(Z = X-XA,int(subs(X = Z+XA,PerProdA[3]),Z)):
Digits := 1000;
evalf(PerProdA,5);
"=======================================";
"===== Main calculations for the segment AB: =====";
"=======================================";
OutputCoincidences(PerProdA,XA):
# Identifying the coincidences at point -0.0901699:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
ResAB := CalcToPoint3(equ,IntA,[XA, -`/`(1,20), XB],350):
# The integral (between points -0.0901699 and 0) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -0.0901699 and 0) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points -0.0901699 and 0) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResAB[1]);
OutputCoincidences(ResAB[3],XB):
# Identifying the coincidences at point 0:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
"Getting expansions around B for the four products of periods (by computing a transition matrix):";
txB := `+`~(XB,`/`(txOFF,4));
MatB := Matrix(3,3,ExplEvalAndComb31(C1,C2,txB,ResAB[3],0,1)):
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
# Zeroes 3 and 1 of the 2nd curve coincide at this point.
ResBC := CalcToPoint3(equ,IntB,[XB, `/`(1,20), `/`(3,20), `/`(3,20), `/`(7,20), `/`(7,20), `/`(3,4), `/`(3,4),
    `/`(10,7), `/`(10,7), `/`(5,2), `/`(5,2), 4, 4, 7, 7, XC],350):
# The integral (between points 0 and 11.0902) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0 and 11.0902) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 0 and 11.0902) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResBC[1]);
OutputCoincidences(ResBC[3],XC):
# Identifying the coincidences at point 11.0902:
# Zeroes 3 and 1 of the 1st curve coincide at this point.
# Zeroes 3 and 1 of the 2nd curve coincide at this point.
"Computing a transition matrix at C (for future reference -- checking solutions):";
txC := `+`~(XC,txOFF);
MatC := Matrix(3,3,ExplEvalAndComb31(C1,C2,txC,ResBC[3],0,1)):
evalf(MatC,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatC := MTM:-round(MatC);
ExplEvalAndComb31 := proc(C1,C2,tx,ToBeCombined, MatrixIfZeroValuesIfOne,debug := 0, weak := 0,
    tofs := [-sqrt(I),-sqrt(I)],conj := 0,Me := [6,4], eps := 10^(-100))
    local b,R1,R2,I1_1,I1_2,J11,J12,J22,temp,Res;
    for b from 1 to 3 do
        #print:
        R1[b] := [evalf(solve(0 = subs(X = tx[b],C1)))];
        R2[b] := [evalf(solve(0 = subs(X = tx[b],C2)))];
        if (abs(Im(R1[b][1]))>10^(-100) or abs(Im(R1[b][2]))>10^(-100)) then
            R1[b] := sort([seq(R1[b][temp],temp = 1 .. 3)],(x, y) -> is(Im(x) < Im(y))):
        else
            R1[b] := sort([seq(Re(R1[b][temp]),temp = 1 .. 3)],(x, y) -> is(x < y));
        end if:
        if (abs(Im(R2[b][1]))>10^(-100) or abs(Im(R2[b][2]))>10^(-100)) then
            R2[b] := sort([seq(R2[b][temp],temp = 1 .. 3)],(x, y) -> is(Im(x) < Im(y))):
        else
            R2[b] := sort([seq(Re(R2[b][temp]),temp = 1 .. 3)],(x, y) -> is(x < y));
        end if:
    end do:
    if debug=1 then
        print("Roots of first curve at the chosen points"= evalf(Vector(3, R1),10));
    end if:
    if debug=1 then
        print("Roots of second curve at the chosen points"= evalf(Vector(3, R2),10));
    end if:
    #Ia_c[b] stands for a-th curve, b-th sample point, c-th principal integral
    #Jcd[b] stands for b-th DE sample point, c-th principal integral on 1 st curve, d-principal integral on 2 nd curve
    if weak = 1 then
        for b from 1 to 3 do
            I1_1[b] := 2*WeakEIntegral(subs(X = tx[b],C1),R1[b][1],R1[b][2],tofs[1],Me[1],Me[2]):
            if debug =2 then
                print(evalf(I1_1[b],15)):
            end if:
            I1_2[b] := 2*WeakEIntegral(subs(X = tx[b],C1),R1[b][2],R1[b][3],tofs[2],Me[1],Me[2]):
            if debug =2 then
                print(evalf(I1_2[b],15)):
            end if:
            J11[b] := evalf(`*`(I1_1[b],I1_1[b])):
            J12[b] := evalf(I1_1[b]*I1_2[b]):
            J22[b] := evalf(`*`(I1_2[b],I1_2[b])):
        end do:
    else
        for b from 1 to 3 do
            temp := 2*evalf(int(sqrt(tofs[1])/sqrt(tofs[1]*subs(X = tx[b],C1)),x)):
            I1_1[b] := evalf(subs(x = R1[b][2]-eps,temp)-subs(x = R1[b][1]+eps,temp)):
            if debug =2 then
                print(evalf(I1_1[b],15)):
            end if:
            #temp ~2 evalf:
            I1_2[b] := evalf(subs(x = R1[b][3]-eps,temp)-subs(x = R1[b][2]+eps,temp)):
            if debug =2 then
                print(evalf(I1_2[b],15)):
            end if:
            #I1_2:=2evalf:
            if debug =2 then
                print(evalf(I1_2[b],15)):
            end if:
            if abs(conj)=1 then
                I1_1[b] := conjugate(I1_2[b])*sign(conj):
            end if:
            if abs(conj)=2 then
                I1_2[b] := conjugate(I1_1[b])*sign(conj):
            end if:
            J11[b] := evalf(`*`(I1_1[b],I1_1[b])):
            J12[b] := evalf(I1_1[b]*I1_2[b]):
            J22[b] := evalf(`*`(I1_2[b],I1_2[b])):
        end do:
    end if:
    print(0):
    Res[1] := PolynomialCombination1(ToBeCombined,3,evalf(tx),J11,MatrixIfZeroValuesIfOne):
    print(1):
    Res[2] := PolynomialCombination1(ToBeCombined,3,evalf(tx),J12,MatrixIfZeroValuesIfOne):
    print(2):
    Res[3] := PolynomialCombination1(ToBeCombined,3,evalf(tx),J22,MatrixIfZeroValuesIfOne):
    print(3):
    Res := [Res[1], Res[2], Res[3]]:
    #if  MatrixIfZeroValuesIfOne = 0 then return Matrix:
    ##if  MatrixIfZeroValuesIfOne = 1 then return Vector:
end proc:
PolynomialCombination1 := proc(solutions,n, tx, tv,MatrixIfZeroValuesIfOne)
    local i, j, a,b, equation, E;
    b := [seq([seq(subs(X = tx[i],solutions[j]),j = 1 .. n)],i = 1 .. n)];
    for i from 1 to n do
        equation[i] := evalf(tv[i] = sum(a[j]*b[i][j],j = 1 .. n));
    end do:
    E := solve({seq(equation[i],i = 1 .. n)},[seq(a[j],j = 1 .. n)]);
    print(55):
    if MatrixIfZeroValuesIfOne = 0 then
        return subs(E[1],[seq(a[j],j=1.. n)]);
    end if:
    if MatrixIfZeroValuesIfOne = 1 then
        return subs(E[1],sum(a[j]*solutions[j],j=1.. n));
    end if:
end proc:
locX := `/`(707,128)-`/`(5*sqrt(5),2):
R := [evalf(solve(0 = subs(X = locX,C2)))];
R := sort([seq(Re(R[temp]),temp = 1 .. 3)],(x, y) -> is(x < y));
temp := 2*evalf(int(1/sqrt(subs(X = locX,C1)),x)):
evalf(subs(x = R[2],temp)-subs(x = R[1],temp));
plot(1/sqrt(subs(X = locX,C1)),x = -.5 .. .7);
plot(temp,x = -.5 .. .7);
evalf(subs(x = R[2]+`^`(10,-100),temp));
