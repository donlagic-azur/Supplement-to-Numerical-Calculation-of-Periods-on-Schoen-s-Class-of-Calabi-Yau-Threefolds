read "PeriodHeader.mpl":
with(PeriodHeader);
with(ArrayTools):
with(ListTools):
with(IntegerRelations):
with(LinearAlgebra):
interface(rtablesize = 70):
Digits := 100;
qAB := EnDim(Vector([-2.332521005560210786518919275617698668566-`*`(5.780791554177295764403353892979281055600,I),
    -6.997563016680632359556757826853096005698+`*`(0,I),
    -2.332521005560210786518919275617698668566+`*`(5.780791554177295764403353892979281055600,I)])):
qBC := EnDim(Vector([16.32764703892147550563243492932389067996-`*`(17.34237466253188729321006167893784316680,I),
    2.332521005560210786518919275617698668566-`*`(11.56158310835459152880670778595856211120,I),
    -2.332521005560210786518919275617698668566-`*`(5.780791554177295764403353892979281055600,I)])):
qCD := EnDim(Vector([2.332521005560210786518919275617698668566+`*`(34.68474932506377458642012335787568633360,I),
    6.997563016680632359556757826853096005698+`*`(5.780791554177295764403353892979281055600,I),
    2.332521005560210786518919275617698668566])):
qDA := EnDim(Vector([-16.32764703892147550563243492932389067996-`*`(11.56158310835459152880670778595856211120,I),
    -2.332521005560210786518919275617698668566+`*`(5.780791554177295764403353892979281055600,I),
    2.332521005560210786518919275617698668566])):
MonodA := EnDim(Matrix(3,3,[[4, 4, 1], [-2, -1, 0], [1, 0, 0]]));
MonodB := EnDim(Matrix(3,3,[[1, -4, 4], [2, -7, 6], [4, -12, 9]]));
MonodC := EnDim(Matrix(3,3,[[4, -36, 81], [2, -17, 36], [1, -8, 16]]));
MonodD := EnDim(Matrix(3,3,[[1, 16, 64], [0, 1, 8], [0, 0, 1]]));
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
# order of magnitude of error is about 10^-35
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
# order of magnitude of error is about 10^-34.523
