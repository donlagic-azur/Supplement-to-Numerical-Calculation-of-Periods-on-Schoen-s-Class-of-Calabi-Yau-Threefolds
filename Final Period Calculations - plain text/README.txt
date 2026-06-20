Period Calculations Extended - plain-text Maple input
=====================================================

The source Maple Workbook contains 24 calculation worksheets. Each worksheet
has been converted to a separate .mpl file so that the workbook's document
boundaries and independent calculations are preserved.

Reusable library:

    PeriodHeader.mpl

This file defines the PeriodHeader package and exports:

    EnDim
    RedVec
    Basis2MaxSearchCount
    LLLrepPSLQ
    QBasisFast
    BasisLLL

Each calculation file begins with:

    read "PeriodHeader.mpl":
    with(PeriodHeader);

Run calculation files with this directory as the current working directory,
for example:

    maple -q "Periods of (1236,1236).mpl"

Only executable Maple input was retained. Workbook formatting and stored
output cells are not part of the plain-text files.

The generated sources have been reformatted for human reading. Statements are
placed on separate lines, control-flow and procedure bodies are indented, and
large vectors and matrices are wrapped across continuation lines.
