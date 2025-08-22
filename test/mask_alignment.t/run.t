No args.

  $ mask_alignment 2> err
  [1]
  $ sed -E 's/(mask_alignment version).*/\1 REDACTED/' err
  mask_alignment version REDACTED
  
  usage: mask_alignment <aln.fa> <max_gap_percent> > <output.fa>
  
  example: mask_alignment silly.aln.fa 95 > silly.aln.masked.fa
           ^^^ That would make fasta file that has any column with >= 95% gaps removed.
  
  Note: To get rid of columns with all gaps, pass in 100.
  


Missing percent

  $ mask_alignment seqs.fa  2> err
  [1]
  $ sed -E 's/(mask_alignment version).*/\1 REDACTED/' err
  mask_alignment version REDACTED
  
  usage: mask_alignment <aln.fa> <max_gap_percent> > <output.fa>
  
  example: mask_alignment silly.aln.fa 95 > silly.aln.masked.fa
           ^^^ That would make fasta file that has any column with >= 95% gaps removed.
  
  Note: To get rid of columns with all gaps, pass in 100.
  

Bad percent.

  $ mask_alignment seqs.fa apple
  ERROR -- mask_ratio can't be coerced to Float: (Invalid_argument "Float.of_string apple")
  [1]

Percent too low.

  $ mask_alignment seqs.fa -1
  ERROR -- max_gap_percent must be in the range (0, 100], got -1.0
  [1]

Percent too high.

  $ mask_alignment seqs.fa 101
  ERROR -- max_gap_percent must be in the range (0, 100], got 101.0
  [1]

First column has 100% gaps.  The max_gap_percent says only columns
with greater than or equal to X percent will be masked.  So at 100,
just the first column is masked.

  $ mask_alignment seqs.fa 100 # drop first column
  >s1
  ABCDEF
  >s2
  -BCDEF
  >s3
  .-CDEF
  >s4
  --.DEF
  >s5
  .---EF
  >s6
  --.-.F

A column has 5/6 = 0.8333 gap ratio.

  $ mask_alignment seqs.fa 83.3 # drops the column
  >s1
  BCDEF
  >s2
  BCDEF
  >s3
  -CDEF
  >s4
  -.DEF
  >s5
  ---EF
  >s6
  -.-.F
  $ mask_alignment seqs.fa 83.4 # keeps the column
  >s1
  ABCDEF
  >s2
  -BCDEF
  >s3
  .-CDEF
  >s4
  --.DEF
  >s5
  .---EF
  >s6
  --.-.F

B column has 4/6 = 0.6666 gap ratio.

  $ mask_alignment seqs.fa 66.6 # drops the column
  >s1
  CDEF
  >s2
  CDEF
  >s3
  CDEF
  >s4
  .DEF
  >s5
  --EF
  >s6
  .-.F
  $ mask_alignment seqs.fa 66.7 # keeps the column
  >s1
  BCDEF
  >s2
  BCDEF
  >s3
  -CDEF
  >s4
  -.DEF
  >s5
  ---EF
  >s6
  -.-.F

C column has 3/6 = 0.5 gap ratio.

  $ mask_alignment seqs.fa 49 # drops the column
  >s1
  DEF
  >s2
  DEF
  >s3
  DEF
  >s4
  DEF
  >s5
  -EF
  >s6
  -.F
  $ mask_alignment seqs.fa 50 # drops the column
  >s1
  DEF
  >s2
  DEF
  >s3
  DEF
  >s4
  DEF
  >s5
  -EF
  >s6
  -.F
  $ mask_alignment seqs.fa 51 # keeps the column
  >s1
  CDEF
  >s2
  CDEF
  >s3
  CDEF
  >s4
  .DEF
  >s5
  --EF
  >s6
  .-.F

D column has 2/6 = 0.3333 gap ratio.

  $ mask_alignment seqs.fa 33 # drops the column
  >s1
  EF
  >s2
  EF
  >s3
  EF
  >s4
  EF
  >s5
  EF
  >s6
  .F
  $ mask_alignment seqs.fa 34 # keeps the column
  >s1
  DEF
  >s2
  DEF
  >s3
  DEF
  >s4
  DEF
  >s5
  -EF
  >s6
  -.F

E column has 1/6 = 0.1666 gap ratio.

  $ mask_alignment seqs.fa 16 # drops the column
  >s1
  F
  >s2
  F
  >s3
  F
  >s4
  F
  >s5
  F
  >s6
  F
  $ mask_alignment seqs.fa 17 # keeps the column
  >s1
  EF
  >s2
  EF
  >s3
  EF
  >s4
  EF
  >s5
  EF
  >s6
  .F

F column has 0/6 = 0.0 gap ratio.

  $ mask_alignment seqs.fa 0 # Zero is an error
  ERROR -- max_gap_percent must be in the range (0, 100], got 0.0
  [1]
  $ mask_alignment seqs.fa 1 # keeps the column
  >s1
  F
  >s2
  F
  >s3
  F
  >s4
  F
  >s5
  F
  >s6
  F
