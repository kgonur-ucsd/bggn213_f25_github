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

    [1] "ACCTTAACCGACTAACCTAATGTTGCAGGGACAGTAGCGAGCATCAAAAC"

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

    [1] "TCAATTCCTAGGCTTTTGTACTCCCCTACACTGCGACTTGCGCATAACCGGGACACCTCT"

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

     [1] "C" "C" "C" "A" "A" "C" "A" "G" "T" "C" "A" "T" "C" "C" "A" "T" "T" "A" "G"
    [20] "C" "T" "A" "G" "G" "A" "T" "G" "A" "G" "G" "C" "A" "C" "G" "C" "A" "T" "T"
    [39] "T" "T" "T" "T" "G" "G" "C" "T" "G" "G" "A" "T"

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

     [1] "W" "G" "D" "D" "L" "E" "Y" "S" "R" "V" "A" "L" "D" "Q" "P" "A" "M" "K" "E"
    [20] "H" "C" "M" "A" "F" "W" "D" "W" "Q" "T" "F" "H" "R" "A" "Y" "D" "H" "G" "D"
    [39] "G" "L" "W" "K" "H" "G" "Q" "F" "K" "I" "A" "N"

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

     [1] "P" "F" "N" "T" "Q" "T" "G" "N" "Q" "S" "L" "K" "Q" "A" "T" "A" "W" "P" "N"
    [20] "T" "H" "F" "K" "T" "F" "P" "V" "P" "W" "H" "C" "T" "P" "I" "W" "P" "P" "G"
    [39] "T" "P" "I" "F" "P" "F" "P" "I" "K" "Q" "F" "A" "T" "D" "Q" "H" "L" "A" "K"
    [58] "R" "A" "P"

Use our new `generate_protein()` function to make random protein
sequences of length 6 to 12 (i.e. one length 6, one length 7, etc. up to
a length of 12.)

One way to do this is “brute force”.

``` r
generate_protein(6)
```

    [1] "N" "W" "K" "C" "R" "D"

``` r
generate_protein(7)
```

    [1] "D" "P" "Q" "Y" "K" "P" "G"

``` r
generate_protein(8)
```

    [1] "E" "W" "E" "W" "G" "Y" "F" "K"

``` r
generate_protein(9)
```

    [1] "M" "Q" "Q" "L" "H" "V" "P" "M" "C"

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
    FSEIYQ
    >7
    WCLRMDK
    >8
    FWGAGGDA
    >9
    HYTRPCHRK
    >10
    LFHWTHEWTR
    >11
    PCMHSKEENLL
    >12
    WIECKDVQNWKS

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
    VVIAQH
    >7
    TKYRGYM
    >8
    TMGSARGK
    >9
    FHYRGDTEP
    >10
    SDNHQWEVRN
    >11
    MSWGMCFCLWL
    >12
    TYEAGECKKWKY
