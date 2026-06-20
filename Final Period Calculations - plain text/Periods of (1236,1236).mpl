read "PeriodHeader.mpl":
with(PeriodHeader);
with(ArrayTools):
with(ListTools):
with(IntegerRelations):
with(LinearAlgebra):
interface(rtablesize = 70):
Digits := 100;
ResAB := EnDim(Vector([2.916453830267157038058284863779351617717455618094668166308200840436252296562899540668707737748834932,
    -`*`(3.202605043014076189267141559548412400856144570699652954727039197729866118287923027004791231922093004,I),
    -4.374680745400735557087427295669027426576183427142002249462301260654378444844349311003061606623252398])):
ResBC := EnDim(Vector([-4.374680745400735557087427295669027426576183427142002249462301260654378444844349311003061606623252398-`*`(6.405210086028152378534283119096824801712289141399305909454078395459732236575846054009582463844186007,I),
    8.749361490801471114174854591338054853152366854284004498924602521308756889688698622006123213246504796,
    -4.374680745400735557087427295669027426576183427142002249462301260654378444844349311003061606623252398+`*`(6.405210086028152378534283119096824801712289141399305909454078395459732236575846054009582463844186007,I)])):
ResCD := EnDim(Vector([1.458226915133578519029142431889675808858727809047334083154100420218126148281449770334353868874417466,
    -`*`(3.202605043014076189267141559548412400856144570699652954727039197729866118287923027004791231922093004,I),
    -8.749361490801471114174854591338054853152366854284004498924602521308756889688698622006123213246504796])):
ResDA := EnDim(Vector([4.374680745400735557087427295669027426576183427142002249462301260654378444844349311003061606623252398,
    -`*`(6.405210086028152378534283119096824801712289141399305909454078395459732236575846054009582463844186007,I),
    -13.12404223620220667126228188700708227972855028142600674838690378196313533453304793300918481986975719])):
MatA := EnDim(Matrix(3,3,[[1, 2, 1], [0, 1, 1], [0, 0, 1]]));
MatB := EnDim(Matrix(3,3,[[1, -2, 1], [-2, 3, -1], [4, -4, 1]]));
MatC := EnDim(Matrix(3,3,[[1, 0, 0], [1, -1, 0], [1, -2, 1]]));
MatD := EnDim(Matrix(3,3,[[1, 0, 0], [-3, 1, 0], [9, -6, 1]]));
ResDA+MatD . ResCD+(MatD . MatC) . ResBC+((MatD . MatC) . MatB) . ResAB;
MonodA := EnDim(Matrix(3,3,[[1, 4, 4], [0, 1, 2], [0, 0, 1]]));
MonodB := (1/MatB . EnDim(Matrix(3,3,[[4, 12, 9], [-6, -17, -12], [9, 24, 16]]))) . MatB;
MonodC := (((1/MatB . (1/MatC)) . EnDim(Matrix(3,3,[[1, 2, 1], [0, 1, 1], [0, 0, 1]]))) . MatC) . MatB;
MonodD := (((((1/MatB . (1/MatC)) . (1/MatD)) . EnDim(Matrix(3,3,[[1, 0, 0], [-6, 1, 0], [36, -12,
    1]]))) . MatD) . MatC) . MatB;
qAB := ResAB:
qBC := 1/MatB . ResBC:
qCD := (1/MatB . (1/MatC)) . ResCD:
qDA := ((1/MatB . (1/MatC)) . (1/MatD)) . ResDA:
wA := Vector([wtf_1, wtf_2, wtf_3, wtf_4]):
wB := wA+qAB:
wC := wB+qBC:
wD := wC+qCD:
((`^`(MonodA,0) . `^`(MonodB,0)) . `^`(MonodC,0)) . `^`(MonodD,0);
qAB+qBC+qCD+qDA;
((MonodA . MonodB) . MonodC) . MonodD;
qAB+MonodB . qBC+(MonodB . MonodC) . qCD+((MonodB . MonodC) . MonodD) . qDA;
NV := Array([]):
NV0 := Array([]):
NVnumb := Array([]):
Monod := [MonodA, MonodB, MonodC, MonodD]:
for M in Monod do
    # tM ~ NullSpace:
    tM := (M-1)^2:
    tM := [(M-1)[1], (M-1)[2], (M-1)[3], (M-1)[4], RedVec(tM[1]), RedVec(tM[2]), RedVec(tM[3]), RedVec(tM[4])]:
    tL := LLL(Matrix([seq([seq(IdentityMatrix(8)[tx,ty],ty = 1 .. 8),
        `^`(10,20)*seq(tM[tx][tz],tz = 1 .. 4)],tx = 1 .. 8)]));
    nv := Array([]):
    for s from 1 to 8 do
        if not (tL[s,8+1]=0 and tL[s,8+2]=0 and tL[s,8+3]=0 and tL[s,8+4]=0) then
            Append(nv,Vector([seq(`^`(10,-20)*tL[s,tx+8],tx = 1 .. 4)]));
        end if:
    end do:
    nv0 := nv;
    Append(NVnumb,numelems(nv));
    Append(NV0,nv);
    nv := add(d(k)[temps]*nv[temps],temps = 1 .. numelems(nv));
    Append(NV,nv);
    print(nv0,"--->",nv);
end do:
NW := [wA, wB, wC, wD]:
A := Array([]):
nv := [seq(subs(k = j,NV[j]),j = 1 .. 4)]:
v := nv[4]+nv[3]+nv[2]+nv[1]:
tD := [seq(seq(d(k)[j] = 0,j = 1 .. NVnumb[k]),k = 1 .. 4)]:
tn := numelems(tD):
tM := Matrix([seq(seq([seq(IdentityMatrix(tn)[add(NVnumb[k],k = 1 .. tx-1)+ty,tz],tz = 1 .. tn),
    `^`(10,20)*seq(subs(tD,subs(d(tx)[ty] = 1,v[tz])),tz = 1 .. 4)],ty = 1 .. NVnumb[tx]),tx = 1 .. 4)]):
tL := LLL(tM,integer):
print(tL);
for s from 1 to tn do
    if tL[s,tn+1]=0 and tL[s,tn+2]=0 and tL[s,tn+3]=0 and tL[s,tn+4]=0 then
        v := [seq(seq(d(tx)[ty] = tL[s,ty+add(NVnumb[j],j = 1 .. tx-1)],ty = 1 .. NVnumb[tx]),tx = 1 .. 4)];
        # print;
        ta := subs(v,nv[1] . NW[1]+nv[2] . NW[2]+nv[3] . NW[3]+nv[4] . NW[4]):
        if ta<>0 then
            print(ta):
            Append(A,ta):
        end if:
    else
        break;
    end if:
end do:
numelems(A);
for s from 1 to numelems(A) do
    if abs(Re(A[s]))<10^(-30) then
        A[s] := I*Im(A[s]);
    end if:
    if abs(Im(A[s]))<10^(-30) then
        A[s] := Re(A[s]);
    end if:
end do:
Eres := Basis2MaxSearchCount(convert(A,list)):
evalf(Vector([Eres]),30);
# order of magnitude of error is about 10^-34
NV := Array([]):
NV0 := Array([]):
NVnumb := Array([]):
Monod := [MonodA, MonodB, MonodC, MonodD]:
for M in Monod do
    tM := (M-1)^(1+Rank((M-1)^2)):
    tL := LLL(Matrix([seq([seq(IdentityMatrix(4)[tx,ty],ty = 1 .. 4),
        `^`(10,20)*seq(tM[tx,tz],tz = 1 .. 4)],tx = 1 .. 4)])):
    nv := Array([]):
    for s from 1 to 4 do
        if tL[s,4+1]=0 and tL[s,4+2]=0 and tL[s,4+3]=0 and tL[s,4+4]=0 then
            Append(nv,Vector([seq(tL[s,tx],tx = 1 .. 4)]));
        end if:
    end do:
    nv0 := nv;
    Append(NVnumb,numelems(nv));
    Append(NV0,nv);
    nv := add(d(k)[temps]*nv[temps],temps = 1 .. numelems(nv));
    Append(NV,nv);
    print(nv0,"--->",nv);
end do:
NW := [wA, wB, wC, wD]:
A := Array([]):
nv := [seq(subs(k = j,NV[j]),j = 1 .. 4)]:
v := nv[4]+nv[3]+nv[2]+nv[1]:
tD := [seq(seq(d(k)[j] = 0,j = 1 .. NVnumb[k]),k = 1 .. 4)]:
tn := numelems(tD):
tM := Matrix([seq(seq([seq(IdentityMatrix(tn)[add(NVnumb[k],k = 1 .. tx-1)+ty,tz],tz = 1 .. tn),
    `^`(10,20)*seq(subs(tD,subs(d(tx)[ty] = 1,v[tz])),tz = 1 .. 4)],ty = 1 .. NVnumb[tx]),tx = 1 .. 4)]):
tL := LLL(tM,integer):
print(tL);
for s from 1 to tn do
    if tL[s,tn+1]=0 and tL[s,tn+2]=0 and tL[s,tn+3]=0 and tL[s,tn+4]=0 then
        v := [seq(seq(d(tx)[ty] = tL[s,ty+add(NVnumb[j],j = 1 .. tx-1)],ty = 1 .. NVnumb[tx]),tx = 1 .. 4)];
        # print;
        ta := subs(v,nv[1] . NW[1]+nv[2] . NW[2]+nv[3] . NW[3]+nv[4] . NW[4]):
        if ta<>0 then
            print(ta):
            Append(A,ta):
        end if:
    else
        break;
    end if:
end do:
numelems(A);
for s from 1 to numelems(A) do
    if abs(Re(A[s]))<10^(-30) then
        A[s] := I*Im(A[s]);
    end if:
    if abs(Im(A[s]))<10^(-30) then
        A[s] := Re(A[s]);
    end if:
end do:
Eres := Basis2MaxSearchCount(convert(A,list)):
evalf(Vector([Eres]),30);
# order of magnitude of error is about 10^-33.921
