**ExampleDb** - *Example database*

Description
--------------------

A small example database subset from Laserson and Vigneault et al, 2014.


Usage
--------------------
```
ExampleDb
```




Format
-------------------

A data.frame with the following columns:

+ `sequence_id`:                Sequence identifier
+ `sequence_alignment`:         IMGT-gapped observed sequence.
+ `germline_alignment`:         IMGT-gapped germline sequence.
+ `germline_alignment_d_mask`:  IMGT-gapped germline sequence with N, P and
D regions masked.
+ `v_call`:                     V region allele assignments.
+ `v_call_genotyped`:           TIgGER corrected V region allele assignment.
+ `d_call`:                     D region allele assignments.
+ `j_call`:                     J region allele assignments.
+ `junction`:                   Junction region sequence.
+ `junction_length`:            Length of the junction region in nucleotides.
+ `np1_length`:                 Number of nucleotides between V and D segments
+ `np2_length`:                 Number of nucleotides between D and J segments
+ `sample_id`:                  Sample identifier
+ `c_call`:                     C region assignment.
+ `duplicate_count`:            Copy number of the sequence
+ `locus`:                      Locus of the receptor



References
-------------------


1. Laserson U and Vigneault F, et al. High-resolution antibody dynamics of
vaccine-induced immune responses.
Proc Natl Acad Sci USA. 2014 111:4928-33.










