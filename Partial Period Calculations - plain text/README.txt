Work in Progress - plain-text Maple input
=========================================

This directory was reconstructed from "Work in Progress.maple".

Contents
--------

- 43 files named "Calc (...).mpl": calculation worksheets with presentation
  output removed and Maple statements formatted for reading.
- ECHeader.mpl: elliptic-curve utilities recovered from the embedded package.
- DEHeader.mpl: differential-equation utilities recovered from the embedded
  package.
- PSHeader.mpl: project-specific utilities recovered from the embedded package.
- Headers.mpl: compatibility loader for all three reusable packages.

Run a calculation with this directory as Maple's current directory. Each
calculation explicitly reads the three package files before using them.

Validation
----------

The generated files were checked with Maple's mint parser. The reusable
packages were also loaded in Maple and their exported procedure names checked.
