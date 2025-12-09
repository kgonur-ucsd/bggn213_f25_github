# Class 6: R Functions
Kavi (PID: A69046927)

- [Our first (silly) function](#our-first-silly-function)
- [A second function](#a-second-function)

All functions in R have at least 3 things:

- A **name**, we pick this and use it to call our function,
- Input **arguments** (there can be multiple)
- The **body** lines of R code that do the work

## Our first (silly) function

Write a function to add some numbers

``` r
add <- function(x, y=1) {
  x + y
}
```

Now we can call this function:

``` r
add(10, 100)
```

    [1] 110

## A second function

Write a function to generate random nucleotide sequences of a user
specified length:

The `sample()` function can be helpful here.

``` r
v <- sample(c("A","C","G","T"), size = 50, replace = TRUE)
```

I want the a 1 element leong characyer vector that looks like “CACAGC”,
not “C”,“A”,“C”,“A”,“G”,“C”.

``` r
paste(v,collapse="")
```

    [1] "TACTTTCCGAGCAGCTGTCCACATTAAAAGCCTGGAGGTCTGGGCTGTCT"

Combined, the code is:

``` r
generate_dna <- function(size=50) {
  v=sample(c("A","C","G","T"), size = size, replace = TRUE)
  paste(v,collapse="")
}
```

Test it:

``` r
generate_dna(60)
```

    [1] "CGCTGGGGACTATATGCACTCGCTCTTTACTTGTTTGGCAGCGGACGTTGGGGCCAGGAC"

``` r
fasta <- FALSE
if(TRUE) {
  cat("HELLO You!")
}
```

    HELLO You!

Add the ability to return a multi-element vector or a single elemet
fasta like vector.

``` r
fasta <- FALSE
generate_fasta <- function(size=50) {
  paste(sample(c("A","C","G","T"), size = size, replace = TRUE),collapse="") 
}
generate_me <- function(size=50) {
  sample(c("A","C","G","T"), size = size, replace = TRUE) 
}
if(fasta == TRUE) {
  generate_fasta()
} else{
  generate_me()
}
```

     [1] "A" "C" "C" "C" "G" "A" "C" "C" "G" "T" "A" "A" "C" "C" "G" "A" "A" "C" "C"
    [20] "A" "A" "C" "G" "T" "A" "T" "T" "T" "A" "G" "C" "G" "G" "C" "C" "T" "T" "C"
    [39] "A" "C" "C" "G" "T" "G" "T" "T" "G" "C" "C" "C"

Now to generate a protein sequence…

``` r
fasta <- FALSE
generate_fasta <- function(size=50) {
  paste(sample(c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"), size = size, replace = TRUE), collapse = "")
}
generate_me <- function(size=50) {
  sample(c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"), size = size, replace = TRUE)
}
if(fasta == TRUE) {
  generate_fasta()
} else {
  generate_me()
}
```

     [1] "F" "W" "Y" "H" "F" "R" "D" "P" "Q" "K" "G" "V" "P" "F" "Y" "E" "N" "L" "M"
    [20] "H" "V" "V" "N" "A" "V" "C" "H" "F" "M" "T" "E" "K" "T" "A" "V" "M" "V" "A"
    [39] "I" "C" "W" "K" "C" "K" "Q" "M" "H" "M" "H" "P"

Better way:

``` r
generate_protein <- function(size = 50, fasta = FALSE) {
  aa <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  seq <- sample(aa, size = size, replace = TRUE)
  if(fasta) {
    return(paste(seq, collapse = ""))
  } else {
    return(seq)
  }
}
```

``` r
generate_protein(60, fasta = FALSE)
```

     [1] "L" "W" "P" "K" "V" "E" "C" "K" "R" "Q" "C" "I" "D" "A" "G" "C" "S" "T" "K"
    [20] "L" "K" "E" "D" "P" "K" "S" "M" "N" "T" "G" "F" "I" "V" "Q" "G" "R" "I" "M"
    [39] "A" "H" "E" "W" "C" "R" "N" "I" "A" "I" "K" "V" "N" "S" "C" "D" "N" "K" "F"
    [58] "Y" "F" "F"

Use our new `generate_protein()` function to make random protein
sequences of length 6 to 12 (i.e. one length 6, one length 7, etc. up to
a length of 12.)

One way to do this is “brute force”.

``` r
generate_protein(6)
```

    [1] "H" "G" "R" "S" "L" "H"

``` r
generate_protein(7)
```

    [1] "Q" "P" "I" "Q" "G" "P" "W"

``` r
generate_protein(8)
```

    [1] "K" "D" "I" "F" "E" "I" "P" "W"

``` r
generate_protein(9)
```

    [1] "I" "M" "C" "M" "K" "E" "T" "K" "L"

Work smarter, not harder:

``` r
lengths <- 6:12
```

``` r
for(i in lengths) {
  cat(">",i,"\n",sep="")
  aa <- generate_protein(i)
  cat(paste(aa,collapse=""))
  cat("\n")
}
```

    >6
    ANYPRF
    >7
    PPGKTFA
    >8
    VDHAPIEW
    >9
    WVFKADQEH
    >10
    KCLNPHRLTA
    >11
    GLLWEASLQQN
    >12
    DYVICMGEAPTV

A third, and better, way to solve this is to use the `apply()` family of
functions, specifically the `sapply()` function in this case.

``` r
lengths <- 6:12
fasta_seqs <- sapply(lengths, function(i) {
  seq <- generate_protein(i)
  paste0(">", i, "\n", paste(seq, collapse = ""))
})
cat(paste(fasta_seqs, collapse = "\n"))
```

    >6
    CFTLKQ
    >7
    RADGLKI
    >8
    CEYYACEE
    >9
    TQAGVRQSG
    >10
    RYTLSERGPI
    >11
    CWLKNVIDMPD
    >12
    MMPQETCRLGSA
