# Maple-readable numerical data extracted from tables/result blocks in:
# Complex numbers use Maple's imaginary unit I. Increase Digits before numerical work if needed.
restart:
Digits := 80:

NumericalTableData := table():

# Section 7.2: Type II final results.
NumericalTableData[TypeII] := table():


NumericalTableData[TypeII][z_to_minus_z_over_2z_plus_1] := table():
NumericalTableData[TypeII][z_to_minus_z_over_2z_plus_1][mobius] := "z -> -z/(2*z+1)":
NumericalTableData[TypeII][z_to_minus_z_over_2z_plus_1][modular_form] := "10/1":
NumericalTableData[TypeII][z_to_minus_z_over_2z_plus_1][IXhat_basis] := [
  10.1775023004434042782842363056 + 0*I,
  0 + 18.4818859696295144823284523403*I
]:
NumericalTableData[TypeII][z_to_minus_z_over_2z_plus_1][IX_basis] := [
  2.91210278431660900940154286414 + 0*I,
  10.1775023004434042782842363056 + 0*I,
  0 + 1.84818859696295144823284523403*I
]:
NumericalTableData[TypeII][z_to_minus_z_over_2z_plus_1][modular_periods_pi2_omega] := [
  1.52662534506651064174263544586 + 4.15842434316664075852390177654*I,
  0 - 1.38614144772221358617463392551*I
]:


NumericalTableData[TypeII][z_to_1_over_4z] := table():
NumericalTableData[TypeII][z_to_1_over_4z][mobius] := "z -> 1/(4*z)":
NumericalTableData[TypeII][z_to_1_over_4z][modular_form] := "21/2":
NumericalTableData[TypeII][z_to_1_over_4z][IXhat_basis] := [
  7.02101068370533129948622785532 + 56.5585057250058555647892632590*I,
  7.02101068370533129948622785532 - 56.5585057250058555647892632590*I
]:
NumericalTableData[TypeII][z_to_1_over_4z][IX_basis] := [
  7.02101068370533129948622785532 + 1.27530934398204965297269158496*I,
  7.02101068370533129948622785532 - 1.27530934398204965297269158496*I,
  0 + 5.85917529922885628954124316713*I
]:
NumericalTableData[TypeII][z_to_1_over_4z][modular_periods_pi2_omega] := [
  -0.752251144682714067802095841641 - 2.01994663303592341302818797353*I,
  2.75825419716995158194101808602 + 6.05983989910777023908456392067*I
]:


NumericalTableData[TypeII][z_to_minus_1_over_2_minus_z] := table():
NumericalTableData[TypeII][z_to_minus_1_over_2_minus_z][mobius] := "z -> -1/2 - z":
NumericalTableData[TypeII][z_to_minus_1_over_2_minus_z][modular_form] := "17/1":
NumericalTableData[TypeII][z_to_minus_1_over_2_minus_z][IXhat_basis] := [
  15.6841094746081896506771079054 + 0*I,
  0 + 17.0625437087419634029794240275*I
]:
NumericalTableData[TypeII][z_to_minus_1_over_2_minus_z][IX_basis] := [
  15.6841094746081896506771079054 + 0*I,
  3.25863063246525860751003087346 + 8.53127185437098170148971201373*I,
  3.25863063246525860751003087346 - 8.53127185437098170148971201373*I,
  0 + 2.27403622608502294651981388337*I
]:
NumericalTableData[TypeII][z_to_minus_1_over_2_minus_z][modular_periods_pi2_omega] := [
  1.38389201246542849858915657984 + 4.51655568760816678314161224223*I,
  -2.76778402493085699717831315967 - 10.5386299377523891606637618985*I
]:


NumericalTableData[TypeII][z_to_minus_1_over_4z_plus_2] := table():
NumericalTableData[TypeII][z_to_minus_1_over_4z_plus_2][mobius] := "z -> -1/(4*z+2)":
NumericalTableData[TypeII][z_to_minus_1_over_4z_plus_2][modular_form] := "73/1":
NumericalTableData[TypeII][z_to_minus_1_over_4z_plus_2][IXhat_basis] := [
  11.4118420461744051807307248949 + 0*I,
  0 + 4.21820528193890201592231653461*I
]:
NumericalTableData[TypeII][z_to_minus_1_over_4z_plus_2][IX_basis] := [
  1.15964224207708731147527428542 + 0*I,
  1.72900723259046311316331046608 + 0*I,
  5.70592102308720259036536244744 + 1.79835796276720374170952048120*I,
  5.70592102308720259036536244744 - 1.79835796276720374170952048120*I,
  0 + 4.21820528193890201592231653461*I
]:
NumericalTableData[TypeII][z_to_minus_1_over_4z_plus_2][modular_periods_pi2_omega] := [
  8.91061639221837116851577154289 + 0.520052705992467371826039094834*I,
  -2.81387886070053826374182256461 - 0.173350901997489123942012962848*I
]:


NumericalTableData[TypeII][z_to_minus_2z_plus_1_over_4z] := table():
NumericalTableData[TypeII][z_to_minus_2z_plus_1_over_4z][mobius] := "z -> -(2*z+1)/(4*z)":
NumericalTableData[TypeII][z_to_minus_2z_plus_1_over_4z][modular_form] := "73/1":
NumericalTableData[TypeII][z_to_minus_2z_plus_1_over_4z][IXhat_basis] := [
  45.6473681846976207229228995794 + 0*I,
  0 + 16.8728211277556080636892661384*I
]:
NumericalTableData[TypeII][z_to_minus_2z_plus_1_over_4z][IX_basis] := [
  4.63856896830834924590109714169 + 0*I,
  6.91602893036185245265324186434 + 0*I,
  22.8236840923488103614614497897 + 7.19343185106881496683808192482*I,
  22.8236840923488103614614497897 - 7.19343185106881496683808192482*I,
  0 + 16.8728211277556080636892661384*I
]:
NumericalTableData[TypeII][z_to_minus_2z_plus_1_over_4z][modular_periods_pi2_omega] := [
  8.91061639221837116851577154289 + 0.520052705992467371826039094834*I,
  -2.81387886070053826374182256461 - 0.173350901997489123942012962848*I
]:

