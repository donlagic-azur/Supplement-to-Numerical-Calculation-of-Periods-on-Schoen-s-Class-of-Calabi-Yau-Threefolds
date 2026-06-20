# Recovered from the DEHeader package embedded in the Maple workbook.
DEHeader := module()
    option package;
    description "Some useful tools about differential equations";
    export Revalf, nformal_solT, nformal_sol, prolong4, prolong3, prolongint3, prolongint4, PolynomialCombination,
        ExtensionByValues, ExtensionBySeries;
    with(DETools):
    Revalf := proc (a)
        local aa;
        aa := evalf(a);
        if is(abs(Im(aa)) < 1/10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000) then
            aa := Re(aa)
        end if;
        return aa
    end proc;
    nformal_solT := proc (oper, n, s, sigma)
        local result, opert, m, i, qq;
        option cache;
        m := degree(oper,X);
        for i from 0 to m do
            qq[i] := coeff(oper,X,i)
        end do;
        if n < 0 then
            result := 0
            elif n = 0 then
                result := mul(subs(D = T+i,qq[0]),i = 1 .. s)
            else
                result := -add(subs(D = T+n-i,qq[i])*nformal_solT(oper,n-i,s,sigma),i = 1 .. m)/subs(D = n+T,qq[0])
            end if;
            result := convert(series(result,T = sigma,degree(oper,D)),polynom);
            return evalc(result)
        end proc;
        nformal_sol := proc (oper, pt, ord)
            local j, i, result, opert, opers, QQ, QQQ, mm, rr, Mr, aa, sigma, multiplicityss, s, ss, multiplicity,
                rresult, ii, cc;
            option remember;
            forget(nformal_solT);
            result := [];
            opert := subs(X = X+pt,oper);
            opert := expand(add(coeff(opert,D,i)*mul(D-j,j = 0 .. i-1)*X^(-i),i = 0 .. degree(opert,D)));
            opert := simplify(expand(opert));
            opert := simplify(opert*denom(opert));
            QQ := coeff(opert,X,0);
            rr := sort(roots(QQ,D),(x, y) -> x[1] <= y[1]);
            Mr := Re(abs(max(solve(QQ))));
            for aa in rr do
                sigma := aa[1];
                multiplicity := aa[2];
                for ss from 0 to Mr-sigma do
                    if subs(D = sigma+ss,QQ) = 0 then
                        s := ss
                    end if
                end do;
                rresult := add(nformal_solT(evalf(opert),i,s,sigma)*X^(i+T),i = 0 .. ord);
                QQQ := mul(subs(D = i+T,coeff(opert,X,0)),i = 1 .. s);
                mm := ldegree(expand(subs(T = T+sigma,QQQ)));
                if 0 < mm then
                    rresult := diff(rresult,T $ mm)
                end if;
                result := [op(result), subs(T = sigma,rresult/mm!)];
                for ii to multiplicity-1 do
                    result := [op(result), subs(T = sigma,diff(rresult,T $ ii))/(mm+ii)!]
                end do
            end do;
            cc := tcoeff(result[1]);
            result := [seq(result[i]/cc,i = 1 .. nops(result))];
            result := subs(X = X-Revalf(pt),result);
            return result
        end proc;
        prolong4 := proc (equ, f, pta, ptb, ptc, ord := 100)
            local i, C, g, E;
            C := subs(X = ptc,[evalf(f), seq(diff(f,X $ i),i = 1 .. 3)]);
            g := nformal_sol(equ,ptb,ord);
            g := evalf(convert(g,polynom));
            g := X1*g[1]+X2*g[2]+X3*g[3]+X4*g[4];
            E := subs(X = ptc,[g, seq(diff(g,X $ i),i = 1 .. 3)]);
            g := subs(solve(E-C,{X1, X2, X3, X4}),g);
            return evalf(g)
        end proc;
        prolong3 := proc (equ, f, pta, ptb, ptc, ord := 100)
            local i, C, g, E;
            C := subs(X = ptc,[evalf(f), seq(diff(f,X $ i),i = 1 .. 2)]);
            g := nformal_sol(equ,ptb,ord);
            g := evalf(convert(g,polynom));
            g := X1*g[1]+X2*g[2]+X3*g[3];
            E := subs(X = ptc,[g, seq(diff(g,X $ i),i = 1 .. 2)]);
            g := subs(solve(E-C,{X1, X2, X3}),g);
            return evalf(g)
        end proc;
        prolongint3 := proc (equ, f, pta, ptb, ptcc, ord := 100)
            local i, C, g, E, ptc;
            ptc := evalf(ptcc);
            C := subs(X = ptc,[evalf(f), seq(diff(f,X $ i),i = 1 .. 3)]);
            g := nformal_sol(equ,ptb,ord);
            g := subs(X = Z+Revalf(ptb),X1*g[1]+X2*g[2]+X3*g[3]);
            g := subs(Z = X-Revalf(ptb),int(g,Z)+X4);
            E := subs(X = ptc,[g, seq(diff(g,X $ i),i = 1 .. 3)]);
            g := subs(solve(evalf(E-C),{X1, X2, X3, X4}),g);
            return g
        end proc;
        prolongint4 := proc (equ, f, pta, ptb, ptcc, ord := 100)
            local i, C, g, E, ptc;
            ptc := evalf(ptcc);
            C := subs(X = ptc,[evalf(f), seq(diff(f,X $ i),i = 1 .. 4)]);
            g := nformal_sol(equ,ptb,ord);
            g := subs(X = Z+Revalf(ptb),X1*g[1]+X2*g[2]+X3*g[3]+X4*g[4]);
            g := subs(Z = X-Revalf(ptb),int(g,Z)+X5);
            E := subs(X = ptc,[g, seq(diff(g,X $ i),i = 1 .. 4)]);
            g := subs(solve(evalf(E-C),{X1, X2, X3, X4, X5}),g);
            return g
        end proc;
        PolynomialCombination := proc (solutions, n, tx, tv, MatrixIfZeroValuesIfOne)
            local i, j, a, b, equation, E;
            b := [seq([seq(subs(X = tx[i],solutions[j]),j = 1 .. n)],i = 1 .. n)];
            for i to n do
                equation[i] := tv[i] = sum(a[j]*evalf(b[i][j]),j = 1 .. n)
            end do;
            E := solve({seq(equation[i],i = 1 .. n)},[seq(a[j],j = 1 .. n)]);
            if MatrixIfZeroValuesIfOne = 0 then
                return subs(E[1],[seq(a[j],j = 1 .. n)])
            end if;
            if MatrixIfZeroValuesIfOne = 1 then
                return subs(E[1],sum(a[j]*solutions[j],j = 1 .. n))
            end if
        end proc;
        ExtensionByValues := proc (equ, offset, ord, n, tx, tv)
            local sersol, solutions, i;
            sersol := nformal_sol(equ,offset,ord);
            for i to n do
                solutions[i] := convert(sersol[i],polynom)
            end do;
            PolynomialCombination(solutions,n,tx,tv)
        end proc;
        ExtensionBySeries := proc (equ, offset, order, n, tx, S)
            local j;
            ExtensionByValues(equ,offset,order,n,tx,[seq(subs(X = tx[j],S),j = 1 .. n)])
        end proc;
    end module:
