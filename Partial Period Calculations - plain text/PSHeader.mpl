# Recovered from the PSHeader package embedded in the Maple workbook.
read "DEHeader.mpl":
with(DEHeader):
with(ArrayTools):
with(IntegerRelations):
PSHeader := module()
    option package;
    description "Shortcuts specific to this project";
    export CalcToPoint3, CalcToPoint4, GetDiffEq4, GetDiffEq5, GetMonodromy3, GetMonodromy, ZRank, ZRankRes,
        OutputCoincidences, ExplEvalAndComb3, ExplEvalAndComb, WeakEIntegral;
    CalcToPoint3 := proc (equ, IntRes, Points, order, prec := 100)
        local ProRes, PREval, k, j, Derv;
        if numelems(Points) mod 2 = 0 then
            print("even number of points");
            return
        end if;
        PREval := Vector(3);
        ProRes := Vector(3);
        Derv := Vector(3);
        for k to 3 do
            ProRes[k] := IntRes[k];
            for j from 3 to numelems(Points) do
                ProRes[k] := prolongint3(equ,ProRes[k],Points[j-2],Points[j],Points[j-1],order);
                j := j+1
            end do;
            j := numelems(Points);
            PREval[k] := evalf(subs(X = Revalf(Points[j]),ProRes[k]));
            printf("The integral (between points %g and %g) of the product of period %d, 1st curve, and period %d, 2nd curve, is:",Revalf(Points[1]),Revalf(Points[j]),ceil(1/2*k),ceil(1/2*k+1/2));
            print(evalf(PREval[k],prec));
            Derv[k] := diff(ProRes[k],X)
        end do;
        [PREval, ProRes, Derv]
    end proc;
    CalcToPoint4 := proc (equ, IntRes, Points, order, prec := 100)
        local ProRes, PREval, k, j, Derv;
        if numelems(Points) mod 2 = 0 then
            print("even number of points");
            return
        end if;
        PREval := Vector(4);
        ProRes := Vector(4);
        Derv := Vector(4);
        for k to numelems(IntRes) do
            ProRes[k] := IntRes[k];
            for j from 3 to numelems(Points) do
                ProRes[k] := prolongint4(equ,ProRes[k],Points[j-2],Points[j],Points[j-1],order);
                j := j+1
            end do;
            j := numelems(Points);
            PREval[k] := evalf(subs(X = Revalf(Points[j]),ProRes[k]));
            printf("The integral (between points %g and %g) of the product of period %d, 1st curve, and period %d, 2nd curve, is:",Revalf(Points[1]),Revalf(Points[j]),ceil(1/2*k),1+(k-1 mod 2));
            print(evalf(PREval[k],prec));
            Derv[k] := diff(ProRes[k],X)
        end do;
        [PREval, ProRes, Derv]
    end proc;
    GetDiffEq4 := proc (S, SerLen, PolyLen, see := 0)
        local k, d0S, d1S, d2S, d3S, d4S, AP, BP, CP, DP, EP, P, Leq, Lvar, R, equ;
        d0S := convert(S,polynom);
        d1S := diff(d0S,X);
        d2S := diff(d1S,X);
        d3S := diff(d2S,X);
        d4S := diff(d3S,X);
        AP := sum(a(k)*X^k,k = 0 .. PolyLen);
        BP := sum(b(k)*X^k,k = 0 .. PolyLen);
        CP := sum(c(k)*X^k,k = 0 .. PolyLen);
        DP := sum(d(k)*X^k,k = 0 .. PolyLen);
        EP := sum(e(k)*X^k,k = 0 .. PolyLen);
        P := AP*d0S+BP*d1S+CP*d2S+DP*d3S;
        P := sort(collect(expand(P),X),X,ascending);
        Leq := seq(coeff(P,X,k) = 0,k = 0 .. 4*PolyLen+14);
        Lvar := seq(a(k),k = 0 .. PolyLen), seq(b(k),k = 0 .. PolyLen), seq(c(k),k = 0 .. PolyLen),
            seq(d(k),k = 0 .. PolyLen);
        R := solve({Leq},[Lvar]);
        if see = 1 then
            print(subs(R[1],AP));
            print(subs(R[1],BP));
            print(subs(R[1],CP));
            print(subs(R[1],DP))
        end if;
        equ := factor(subs(R[1],AP)+D*subs(R[1],BP)+D^2*subs(R[1],CP)+D^3*subs(R[1],DP));
        equ := subs([seq(Lvar[k] = 1,k = 1 .. numelems([Lvar]))],equ)
    end proc;
    GetDiffEq5 := proc (S, SerLen, PolyLen, see := 0)
        local k, d0S, d1S, d2S, d3S, d4S, AP, BP, CP, DP, EP, P, Leq, Lvar, R, equ;
        d0S := convert(S,polynom);
        d1S := diff(d0S,X);
        d2S := diff(d1S,X);
        d3S := diff(d2S,X);
        d4S := diff(d3S,X);
        AP := sum(a(k)*X^k,k = 0 .. PolyLen);
        BP := sum(b(k)*X^k,k = 0 .. PolyLen);
        CP := sum(c(k)*X^k,k = 0 .. PolyLen);
        DP := sum(d(k)*X^k,k = 0 .. PolyLen);
        EP := sum(e(k)*X^k,k = 0 .. PolyLen);
        P := AP*d0S+BP*d1S+CP*d2S+DP*d3S+EP*d4S;
        P := sort(collect(expand(P),X),X,ascending);
        Leq := seq(coeff(P,X,k) = 0,k = 0 .. 5*PolyLen+15);
        Lvar := seq(a(k),k = 0 .. PolyLen), seq(b(k),k = 0 .. PolyLen), seq(c(k),k = 0 .. PolyLen),
            seq(d(k),k = 0 .. PolyLen), seq(e(k),k = 0 .. PolyLen);
        R := solve({Leq},[Lvar]);
        if see = 1 then
            print(subs(R[1],AP));
            print(subs(R[1],BP));
            print(subs(R[1],CP));
            print(subs(R[1],DP));
            print(subs(R[1],EP))
        end if;
        equ := factor(subs(R[1],AP)+D*subs(R[1],BP)+D^2*subs(R[1],CP)+D^3*subs(R[1],DP)+D^4*subs(R[1],EP));
        equ := subs([seq(Lvar[k] = 1,k = 1 .. numelems([Lvar]))],equ)
    end proc;
    GetMonodromy3 := proc (tx, ToBeCombined, p := 15)
        local b, J11, J12, J21, J22, Res;
        for b to 3 do
            J11[b] := evalf(subs(X = evalf(tx[b]),ToBeCombined[1]));
            J12[b] := evalf(subs(X = evalf(tx[b]),ToBeCombined[2]));
            J22[b] := evalf(subs(X = evalf(tx[b]),ToBeCombined[3]))
        end do;
        Res[1] := PolynomialCombination(ToBeCombined,3,evalf(`+`~(-I*10^(-p),tx)),J11,0);
        Res[2] := PolynomialCombination(ToBeCombined,3,evalf(`+`~(-I*10^(-p),tx)),J12,0);
        Res[3] := PolynomialCombination(ToBeCombined,3,evalf(`+`~(-I*10^(-p),tx)),J22,0);
        Res := [Res[1], Res[2], Res[3]];
        return Matrix(3,3,Res)
    end proc;
    GetMonodromy := proc (tx, ToBeCombined, p := 15)
        local b, J11, J12, J21, J22, Res;
        for b to 4 do
            J11[b] := evalf(subs(X = evalf(tx[b]),ToBeCombined[1]));
            J12[b] := evalf(subs(X = evalf(tx[b]),ToBeCombined[2]));
            J21[b] := evalf(subs(X = evalf(tx[b]),ToBeCombined[3]));
            J22[b] := evalf(subs(X = evalf(tx[b]),ToBeCombined[4]))
        end do;
        Res[1] := PolynomialCombination(ToBeCombined,4,evalf(`+`~(-I*10^(-p),tx)),J11,0);
        Res[2] := PolynomialCombination(ToBeCombined,4,evalf(`+`~(-I*10^(-p),tx)),J12,0);
        Res[3] := PolynomialCombination(ToBeCombined,4,evalf(`+`~(-I*10^(-p),tx)),J21,0);
        Res[4] := PolynomialCombination(ToBeCombined,4,evalf(`+`~(-I*10^(-p),tx)),J22,0);
        Res := [Res[1], Res[2], Res[3], Res[4]];
        return Matrix(4,4,Res)
    end proc;
    ZRank := proc (V, prec := 100, ts := 1000, coin := 0)
        local R, N, n, k, j, s, d, ss, sss;
        d := Digits;
        Digits := prec;
        n := numelems(V);
        N := Array([seq(k,k = 1 .. n)]);
        for k to n do
            if abs(V[n+1-k]) < 1/ts then
                Remove(N,n+1-k);
                ss := v[n+1-k];
                if coin = 1 then
                    ss := subs(v[1] = alpha[`&mdash;`(1,2)],v[2] = alpha[`&mdash;`(2,3)],v[3] = alpha[`&mdash;`(3,1)],ss)
                end if;
                print(ss = 0)
            end if
        end do;
        while 1 < numelems(N) do
            n := numelems(N);
            R := PSLQ(V[convert(N,list)]);
            s := 0;
            for j to numelems(R) do
                s := s+R[j]*v[N[j]]
            end do;
            for j to numelems(R) do
                if R[j] <> 0 then
                    if abs(R[j]) < ts then
                        Remove(N,j)
                    end if;
                    break
                end if
            end do;
            if n = numelems(N) then
                break
            end if;
            ss := s;
            sss := s;
            if coin = 1 then
                ss := subs(v[1] = alpha[`&mdash;`(1,2)],v[2] = alpha[`&mdash;`(2,3)],v[3] = alpha[`&mdash;`(3,1)],ss)
            end if;
            print(ss = 0)
        end do;
        Digits := d;
        if coin = 1 then
            return [sss, numelems(N)]
        end if;
        numelems(N)
    end proc;
    ZRankRes := proc (V, prec := 30)
        local W, k;
        W := Vector(2*numelems(V));
        for k to numelems(V) do
            W[-1+2*k] := evalf(Re(V[k]),prec);
            W[2*k] := evalf(Im(V[k]),prec)
        end do;
        ZRank(W,prec)
    end proc;
    OutputCoincidences := proc (tDerv, offset)
        local Derv, rndc, act, M, N, n1, n2, i, j, s, t, d;
        printf("Identifying the coincidences at point %g:",offset);
        Derv[1] := subs(X = X+offset,tDerv[1]);
        Derv[2] := subs(X = X+offset,tDerv[2]);
        Derv[3] := subs(X = X+offset,tDerv[3]);
        if numelems(tDerv) = 4 then
            Derv[4] := subs(X = X+offset,tDerv[4])
        else
            Derv[4] := Derv[3];
            Derv[3] := Derv[2]
        end if;
        rndc := 1.320548+.932582*I;
        act := 0;
        M := [[Derv[1], Derv[2]], [Derv[3], Derv[4]]];
        N := M;
        for i to 2 do
            for j to 2 do
                M[i][j] := subs(ln(X) = W,convert(series(M[i][j],X,1),polynom));
                N[i][j] := coeff(M[i][j],W,2);
                if 1/100000000000000000000000000000000000000000000000000 < abs(N[i][j]) then
                    act := 1
                end if
            end do
        end do;
        if act = 0 then
            for i to 2 do
                for j to 2 do
                    N[i][j] := coeff(M[i][j],W,1)
                end do
            end do
        end if;
        n1 := Re(rndc*N[1][1]+rndc^2*N[1][2]);
        n2 := Re(rndc*N[2][1]+rndc^2*N[2][2]);
        s := ZRank([n1, n2, -n1-n2],30,1000,1);
        t := subs(v[1] = 1,v[2] = 1,v[3] = 1,s[1]);
        if coeff(s[1],v[1])-t mod 2 = 0 then
            d := 1
        end if;
        if coeff(s[1],v[2])-t mod 2 = 0 then
            d := 2
        end if;
        if coeff(s[1],v[3])-t mod 2 = 0 then
            d := 3
        end if;
        printf(ifelse(1 < s[2],"The 1st curve doesn't have a singularity at this point.","Zeroes %d and %d of the 1st curve coincide at this point."),d,1+(d mod 3));
        n1 := Re(rndc*N[1][1]+rndc^2*N[2][1]);
        n2 := Re(rndc*N[1][2]+rndc^2*N[2][2]);
        s := ZRank([n1, n2, -n1-n2],30,1000,1);
        t := subs(v[1] = 1,v[2] = 1,v[3] = 1,s[1]);
        if coeff(s[1],v[1])-t mod 2 = 0 then
            d := 1
        end if;
        if coeff(s[1],v[2])-t mod 2 = 0 then
            d := 2
        end if;
        if coeff(s[1],v[3])-t mod 2 = 0 then
            d := 3
        end if;
        printf(ifelse(1 < s[2],"The 2nd curve doesn't have a singularity at this point.","Zeroes %d and %d of the 2nd curve coincide at this point."),d,1+(d mod 3));
        return
    end proc;
    ExplEvalAndComb3 := proc (C1, C2, tx, ToBeCombined, MatrixIfZeroValuesIfOne, debug := 0, bindex := 0,
        tofs := [1, 1], conj := 0, Me := [6, 4], eps := `^`(10,-100))
        local b, R1, R2, I1_1, I1_2, J11, J12, J22, temp, Res;
        for b to 3 do
            R1[b] := [evalf(solve(0 = subs(X = tx[b],C1)))];
            R2[b] := [evalf(solve(0 = subs(X = tx[b],C2)))];
            if 1/10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 < abs(Im(R1[b][1])) or 1/10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 < abs(Im(R1[b][2])) then
                R1[b] := sort([seq(R1[b][temp],temp = 1 .. 3)],(x, y) -> is(Im(x) < Im(y)))
            else
                R1[b] := sort([seq(Re(R1[b][temp]),temp = 1 .. 3)],(x, y) -> is(x < y))
            end if;
            if 1/10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 < abs(Im(R2[b][1])) or 1/10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 < abs(Im(R2[b][2])) then
                R2[b] := sort([seq(R2[b][temp],temp = 1 .. 3)],(x, y) -> is(Im(x) < Im(y)))
            else
                R2[b] := sort([seq(Re(R2[b][temp]),temp = 1 .. 3)],(x, y) -> is(x < y))
            end if
        end do;
        if debug = 1 then
            print('Roots*of*first*curve*at*the*chosen*points' = evalf(Vector(3,R1),10))
        end if;
        if debug = 1 then
            print('Roots*of*second*curve*at*the*chosen*points' = evalf(Vector(3,R2),10))
        end if;
        if bindex = -1 then
            for b to 3 do
                I1_1[b] := 2*WeakEIntegral(subs(X = tx[b],C1),R1[b][1],R1[b][2],tofs[1],Me[1],Me[2]);
                if debug = 2 then
                    print(evalf(I1_1[b],15))
                end if;
                I1_2[b] := 2*WeakEIntegral(subs(X = tx[b],C1),R1[b][2],R1[b][3],tofs[2],Me[1],Me[2]);
                if debug = 2 then
                    print(evalf(I1_2[b],15))
                end if;
                J11[b] := evalf(I1_1[b]*I1_1[b]);
                J12[b] := evalf(I1_1[b]*I1_2[b]);
                J22[b] := evalf(I1_2[b]*I1_2[b])
            end do
        else
            for b to 3 do
                temp := 2*evalf(int(sqrt(tofs[1])/sqrt(tofs[1]*subs(X = tx[b],C1)),x));
                I1_1[b] := evalf(subs(x = R1[b][2]-eps,temp)-subs(x = R1[b][1]+eps,temp));
                if debug = 2 then
                    print(evalf(I1_1[b],15))
                end if;
                I1_2[b] := evalf(subs(x = R1[b][3]-eps,temp)-subs(x = R1[b][2]+eps,temp));
                if debug = 2 then
                    print(evalf(I1_2[b],15))
                end if;
                if abs(conj) = 1 then
                    I1_1[b] := conjugate(I1_2[b])*sign(conj)
                end if;
                if abs(conj) = 2 then
                    I1_2[b] := conjugate(I1_1[b])*sign(conj)
                end if;
                J11[b] := evalf(tx[b]^bindex*I1_1[b]*I1_1[b]);
                J12[b] := evalf(tx[b]^bindex*I1_1[b]*I1_2[b]);
                J22[b] := evalf(tx[b]^bindex*I1_2[b]*I1_2[b])
            end do
        end if;
        Res[1] := PolynomialCombination(ToBeCombined,3,evalf(tx),J11,MatrixIfZeroValuesIfOne);
        Res[2] := PolynomialCombination(ToBeCombined,3,evalf(tx),J12,MatrixIfZeroValuesIfOne);
        Res[3] := PolynomialCombination(ToBeCombined,3,evalf(tx),J22,MatrixIfZeroValuesIfOne);
        Res := [Res[1], Res[2], Res[3]]
    end proc;
    ExplEvalAndComb := proc (C1, C2, tx, ToBeCombined, MatrixIfZeroValuesIfOne, debug := 0, weak := 0, tofs := [1,
        1, 1, 1], conj := [0, 0], Me := [6, 4])
        local b, R1, R2, I1_1, I1_2, I2_1, I2_2, J11, J12, J21, J22, Res;
        for b to 4 do
            R1[b] := sort([solve(0 = subs(X = tx[b],C1))],(x, y) -> is(Im(x) < Im(y)) or is(x < y));
            R2[b] := sort([solve(0 = subs(X = tx[b],C2))],(x, y) -> is(Im(x) < Im(y)) or is(x < y))
        end do;
        if debug = 1 then
            print('Roots*of*first*curve*at*the*chosen*points' = Vector(4,R1))
        end if;
        if debug = 1 then
            print('Roots*of*second*curve*at*the*chosen*points' = Vector(4,R2))
        end if;
        if weak = 1 then
            for b to 4 do
                I1_1[b] := 2*WeakEIntegral(subs(X = tx[b],C1),R1[b][1],R1[b][2],tofs[1],Me[1],Me[2]);
                if debug = 2 then
                    print(evalf(I1_1[b],15))
                end if;
                I1_2[b] := 2*WeakEIntegral(subs(X = tx[b],C1),R1[b][2],R1[b][3],tofs[2],Me[1],Me[2]);
                if debug = 2 then
                    print(evalf(I1_2[b],15))
                end if;
                I2_1[b] := 2*WeakEIntegral(subs(X = tx[b],C2),R2[b][1],R2[b][2],tofs[3],Me[1],Me[2]);
                if debug = 2 then
                    print(evalf(I2_1[b],15))
                end if;
                I2_2[b] := 2*WeakEIntegral(subs(X = tx[b],C2),R2[b][2],R2[b][3],tofs[4],Me[1],Me[2]);
                if debug = 2 then
                    print(evalf(I2_2[b],15))
                end if;
                J11[b] := evalf(I1_1[b]*I2_1[b]);
                J12[b] := evalf(I1_1[b]*I2_2[b]);
                J21[b] := evalf(I1_2[b]*I2_1[b]);
                J22[b] := evalf(I1_2[b]*I2_2[b])
            end do
        else
            for b to 4 do
                I1_1[b] := 2*evalf(int(sqrt(tofs[1])/sqrt(tofs[1]*subs(X = tx[b],C1)),x = R1[b][1] .. R1[b][2]));
                if debug = 2 then
                    print(evalf(I1_1[b],15))
                end if;
                I1_2[b] := 2*evalf(int(sqrt(tofs[2])/sqrt(tofs[2]*subs(X = tx[b],C1)),x = R1[b][2] .. R1[b][3]));
                if debug = 2 then
                    print(evalf(I1_2[b],15))
                end if;
                I2_1[b] := 2*evalf(int(sqrt(tofs[3])/sqrt(tofs[3]*subs(X = tx[b],C2)),x = R2[b][1] .. R2[b][2]));
                if debug = 2 then
                    print(evalf(I2_1[b],15))
                end if;
                I2_2[b] := 2*evalf(int(sqrt(tofs[4])/sqrt(tofs[4]*subs(X = tx[b],C2)),x = R2[b][2] .. R2[b][3]));
                if debug = 2 then
                    print(evalf(I2_2[b],15))
                end if;
                if abs(conj[1]) = 1 then
                    I1_1[b] := conjugate(I1_2[b])*sign(conj[1])
                end if;
                if abs(conj[1]) = 2 then
                    I1_2[b] := conjugate(I1_1[b])*sign(conj[1])
                end if;
                if abs(conj[2]) = 1 then
                    I2_1[b] := conjugate(I2_2[b])*sign(conj[2])
                end if;
                if abs(conj[2]) = 2 then
                    I2_2[b] := conjugate(I2_1[b])*sign(conj[2])
                end if;
                J11[b] := evalf(I1_1[b]*I2_1[b]);
                J12[b] := evalf(I1_1[b]*I2_2[b]);
                J21[b] := evalf(I1_2[b]*I2_1[b]);
                J22[b] := evalf(I1_2[b]*I2_2[b])
            end do
        end if;
        Res[1] := PolynomialCombination(ToBeCombined,4,evalf(tx),J11,MatrixIfZeroValuesIfOne);
        Res[2] := PolynomialCombination(ToBeCombined,4,evalf(tx),J12,MatrixIfZeroValuesIfOne);
        Res[3] := PolynomialCombination(ToBeCombined,4,evalf(tx),J21,MatrixIfZeroValuesIfOne);
        Res[4] := PolynomialCombination(ToBeCombined,4,evalf(tx),J22,MatrixIfZeroValuesIfOne);
        Res := [Res[1], Res[2], Res[3], Res[4]]
    end proc;
    WeakEIntegral := proc (P, a, b, tofs := -sqrt(I), M := 6, e := 4)
        local sum, k, n, l, ll, T, N, d, u;
        N := M^e;
        d := Digits;
        Digits := 200;
        for k to N do
            T[k] := a+1/2*(b-a)/(1+(N-k)/M)^2
        end do;
        for k from N to 2*N-1 do
            T[k] := b-1/2*(b-a)/(1+(k-N)/M)^2
        end do;
        sum := 0;
        l := 1/tofs;
        ll := 1;
        for k from 2 to 2*N-1 do
            n := evalf(subs(x = T[k],P));
            ll := ll*sqrt(l/n);
            sum := sum+ll*(T[k]-T[k-1]);
            l := n
        end do;
        sum := sqrt(tofs)*sum;
        Digits := d;
        sum
    end proc;
end module:
