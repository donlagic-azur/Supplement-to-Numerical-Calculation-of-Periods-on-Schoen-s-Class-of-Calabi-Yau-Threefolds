# Recovered from the ECHeader package embedded in the Maple workbook.
ECHeader := module()
    option package;
    description "Some useful tools about elliptic curves";
    export Wei, G2, G3, Delta, J, Period;
    Wei := proc (F)
        local FF;
        FF := -lhs(F)+rhs(F);
        FF := 4*FF/coeff(FF,x,3);
        FF := subs(y = y/sqrt(-coeff(FF,y,2)),FF);
        FF := subs(y = y+1/2*coeff(FF,y,1),FF);
        FF := subs(x = x-1/12*coeff(FF,x,2),FF)+y^2
    end proc;
    G2 := proc (f)
        local ff, fff;
        ff := f/coeff(f,x,3);
        fff := subs(x = x-1/3*coeff(ff,x,2),ff);
        -factor(expand(4*coeff(fff,x,1)))
    end proc;
    G3 := proc (f)
        local ff, fff;
        ff := f/coeff(f,x,3);
        fff := subs(x = x-1/3*coeff(ff,x,2),ff);
        -factor(expand(4*coeff(fff,x,0)))
    end proc;
    Delta := proc (f)
        factor(G2(f)^3-27*G3(f)^2)
    end proc;
    J := proc (f)
        G2(f)^3/Delta(f)
    end proc;
    Period := proc (f)
        local g2n;
        1/6*Pi*12^(3/4)/G2(f)^(1/4)*hypergeom([1/12, 5/12],[1],1/J(f))
    end proc;
end module:
