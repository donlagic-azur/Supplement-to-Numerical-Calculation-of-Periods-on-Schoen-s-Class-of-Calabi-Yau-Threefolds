This repository complements the paper
"Numerical Calculation of Periods on Schoen's Class of Calabi-Yau Threefolds"
by Azur Đonlagić.

It contains the final results obtained in the paper (following the notation presented in Section 7) as well as the intermediate results appearing in the explicit example covered in Section 4.
It also contains the following two Maple documents, whose Workbook files have also been exported in Maple-readable plain-text form to the corresponding folders.

1. The document "Partial Period Calculations" contains the first half of the calculation. The Workbook files correspond to the different fiber products of surfaces considered (two charts per fiber product, also multiple 3-forms in Type III) and their names follow these conventions: 

-- Type I:   Calc (abcd, abcd) corresponds to the self-fiber product of the Beauville surface with fibers of types I_a, I_b, I_c, I_d;
	     Calc (dcba, dcba) corresponds to the inverse chart of the same threefold

-- Type II:  Calc (10236, 01abc) corresponds to the permutation taking fibers of types I_2, I_3, I_6 to those of types I_a, I_b, I_c;
	     Calc (63201, cba10) corresponds to the inverse chart of the same threefold

-- Type III: Calc (X(N), X(N))_k corresponds to the k-th basis 3-form on self-fiber product of the module surface X(N);
	     Calc (X(N), X(N))_k inv corresponds to the inverse chart of the same threefold with the same 3-form

2. The document "Final Period Calculations" contains the second half of the calculation, obtaining the modules of periods from the partial periods between pairs of singular fibers. The Workbook file names follow the same convention, except that the two charts are not considered separately and that in the case of X(10) all 3-forms were extracted from the same computation done for a generic 3-form.

NOTE: These files have been computed in Maple 2023. Due to the inner workings of Maple, some calculations may not terminate properly in newer versions without tweaking certain parameters.