# Recovered from the PeriodHeader package embedded in the Maple workbook.
with(ArrayTools):
with(ListTools):
with(IntegerRelations):
with(LinearAlgebra):
PeriodHeader := module()
    option package;
    description "All tools used for calculating periods";
    export EnDim, RedVec, Basis2MaxSearchCount, LLLrepPSLQ, QBasisFast, BasisLLL;
    local NearProj, BinInsCount;
    EnDim := proc (A)
        local t;
        if type(A,Vector) = true then
            return Vector([A[1], A[2], A[2], A[3]])
        end if;
        t := (2000000000000000000000000000000*A[1,1]*A[2,2]+100000000000000000000*A[1,1]*A[3,3]-1000000000000000000000000000000*A[1,2]*A[2,1]+A[1,2]*A[2,3]+100000000000000000000*A[2,1]*A[2,3]+10000000000*A[2,1]*A[3,2])/(2*A[1,3]+20000000000*A[3,1]+100000000000000000000*A[2,2]+2000000000000000000000000000000*A[1,1]);
        return rtable(1 .. 4,1 .. 4,[[A[1,1], 1/2*A[1,2], 1/2*A[1,2], A[1,3]], [A[2,1], t, A[2,2]-t, A[2,3]], [A[2,1],
            A[2,2]-t, t, A[2,3]], [A[3,1], 1/2*A[3,2], 1/2*A[3,2], A[3,3]]],subtype = Matrix)
    end proc;
    RedVec := proc (V)
        local v;
        v := igcd(seq(v,v in V));
        if v = 0 then
            return V
        else
            return V/v
        end if
    end proc;
    Basis2MaxSearchCount := proc (A)
        local L, t, n, c;
        L := Array([]);
        for t in A do
            BinInsCount(L,[evalf(t), 1])
        end do;
        while 1 = 1 do
            n := NearProj(L[1][1],L[2][1]);
            if n <> 0 then
                t := L[1][1]-n*L[2][1];
                c := L[1][2]+abs(n)*L[2][2];
                Remove(L,1);
                BinInsCount(L,[t, c]);
                next
            end if;
            if evalf(abs(L[3][1])) < 1/100000000000000000000 then
                break
            end if;
            n := NearProj(L[1][1],L[3][1]);
            if n <> 0 then
                t := L[1][1]-n*L[3][1];
                c := L[1][2]+abs(n)*L[3][2];
                Remove(L,1);
                BinInsCount(L,[t, c]);
                next
            end if;
            n := NearProj(L[2][1],L[3][1]);
            if n <> 0 then
                t := L[2][1]-n*L[3][1];
                c := L[2][2]+abs(n)*L[3][2];
                Remove(L,2);
                BinInsCount(L,[t, c]);
                next
            end if
        end do;
        printf("order of magnitude of error is about 10^%g",evalf(log(max(L[1][2],L[2][2]),10)-35,5));
        return [L[1][1], L[2][1]]
    end proc;
    LLLrepPSLQ := proc (V, prec := 30, debug := 0)
        local M, T, n, i, j;
        n := numelems(V);
        M := convert(IdentityMatrix(n),list);
        M := Transpose(Matrix([seq([seq(M[i*n+j],j = 1 .. n)],i = 0 .. n-1),
            convert(10^prec*convert(V,Vector),list)]));
        M := LLL(Re(M)+1.2325251289365235239852705270275272475765759527352*Im(M))[1];
        Remove(M,n+1);
        T := ListTools:-DotProduct(convert(V,list),convert(M,list));
        if 1/100000000000000000000 < abs(T) and debug = 0 then
            M := 0*convert(M,Vector)
        end if;
        M
    end proc;
    QBasisFast := proc (V, prec := 30, ts := 100000)
        local d, A, B, C, k, c, b, fact;
        fact := 1.23786285372652468467964946453637489038;
        A := [seq(Re(V[k])+I*Im(V[k])*fact,k = 1 .. numelems(V))];
        B := Array([]);
        for k in A do
            if evalf(abs(k)) < 1/ts then
                next
            end if;
            Append(B,k);
            if numelems(B) < 2 then
                next
            end if;
            b := 0;
            C := Re(Vector(LLLrepPSLQ(convert(B,list),prec)));
            for c in C do
                if c <> 0 then
                    b := 1
                end if;
                if ts < evalf(abs(c)) then
                    b := 0;
                    break
                end if
            end do;
            if b = 1 then
                Remove(B,numelems(B))
            end if
        end do;
        [seq(Re(B[k])+I*Im(B[k])/fact,k = 1 .. numelems(B))]
    end proc;
    BasisLLL := proc (V, csize1 := 100, QB1 := [], prec := 40, ts := 10000000, tsL := 100000)
        local qb, QB, A, N, L, t, n, c, C, p, q, k, fact, d, l, csize;
        N := numelems(V);
        A := [];
        if type(csize1,integer) then
            csize := max(csize1,floor(N/csize1));
            if csize < N then
                for k from 0 to floor((N-2)/csize) do
                    C := BasisLLL(convert(V,list)[csize*k+1 .. min(csize*(k+1),N)],csize1,QB1,prec,ts,tsL);
                    A := [op(A), op(C)];
                    print(min(csize*(k+1),N))
                end do
            else
                A := V
            end if
        else
            for k to numelems(csize1)-1 do
                C := BasisLLL(convert(V,list)[csize1[k]+1 .. csize1[k+1]],csize1[k+1],QB1,prec,ts,tsL);
                A := [op(A), op(C)];
                print(csize1[k+1])
            end do
        end if;
        if numelems(QB1) = 0 then
            QB := QBasisFast(sort(A,(x, y) -> is(evalf(abs(x)-abs(y)) < 0)),prec,ts)
        else
            QB := QB1
        end if;
        if numelems(QB) = 0 then
            return [0]
        end if;
        fact := 1.23786285372654642468467964946453637489038;
        QB := [seq(Re(QB[k])+I*Im(QB[k])*fact,k = 1 .. numelems(QB))];
        A := [seq(Re(A[k])+I*Im(A[k])*fact,k = 1 .. numelems(A))];
        A := convert(A,Array);
        Append(A,0);
        QB := convert(QB,Array);
        qb := numelems(QB);
        N := numelems(A);
        L := Array([]);
        d := Digits;
        assume(c,integer);
        for t in A do
            C := Re(Vector(LLLrepPSLQ(convert(Append(QB,t,inplace = false),list),prec)));
            c := C[qb+1];
            if c = 0 then
                print(QB,t,C)
            end if;
            C := -identify(1/c*Vector([C[1 .. qb]]));
            Append(L,C)
        end do;
        QB := [seq(Re(QB[k])+I*Im(QB[k])/fact,k = 1 .. qb)];
        L := LLL(Matrix([seq([seq(IdentityMatrix(N)[p,q]/tsL,q = 1 .. N), seq(L[p][q],q = 1 .. qb)],p = 1 .. N)]));
        L := [seq([seq(L[N-qb+p][N+q],q = 1 .. qb)],p = 1 .. qb)];
        return [seq(DotProduct(L[p],QB,conjugate = false),p = 1 .. qb)]
    end proc;
    NearProj := proc (x, y)
        return round(999999999999999999999999999999/1000000000000000000000000000000*Re(x*conjugate(y))/Re(y*conjugate(y)))
    end proc;
    BinInsCount := proc (L, z)
        Insert(L,1+BinaryPlace(L,z,(x, y) -> is(0 < evalf(abs(x[1])-abs(y[1])))),z)
    end proc;
end module:
