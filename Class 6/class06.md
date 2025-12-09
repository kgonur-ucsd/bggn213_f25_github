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

    [1] "ACTGTGGGCCTGCTGCTAAGGATCTCCACCTCCGCTCCCTCGTGGAAATT"

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

    [1] "CTAGGTCACCGTTTATAGCCCCGAAGGGAAACATCACCAGTCTTCCCTTTGGGACAGCCG"

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

     [1] "T" "C" "G" "T" "A" "C" "C" "G" "G" "C" "G" "T" "T" "T" "C" "G" "T" "T" "T"
    [20] "G" "T" "T" "A" "T" "C" "A" "T" "T" "T" "A" "T" "C" "G" "C" "G" "A" "C" "T"
    [39] "A" "T" "A" "G" "G" "A" "T" "A" "T" "T" "G" "G"

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

     [1] "Y" "F" "T" "C" "Y" "N" "S" "N" "S" "C" "T" "C" "V" "E" "S" "V" "K" "G" "K"
    [20] "L" "S" "M" "C" "R" "R" "M" "V" "M" "C" "N" "N" "V" "N" "P" "M" "V" "Q" "F"
    [39] "V" "D" "F" "M" "V" "I" "P" "N" "C" "A" "M" "I"

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

     [1] "Y" "Y" "K" "P" "H" "D" "A" "C" "R" "S" "I" "H" "S" "F" "W" "S" "Q" "F" "V"
    [20] "S" "Q" "H" "P" "H" "G" "K" "I" "E" "I" "V" "I" "C" "K" "Q" "K" "H" "A" "E"
    [39] "D" "G" "Q" "M" "P" "Q" "P" "I" "K" "M" "W" "A" "C" "T" "T" "T" "W" "K" "F"
    [58] "V" "P" "S"

Use our new `generate_protein()` function to make random protein
sequences of length 6 to 12 (i.e. one length 6, one length 7, etc. up to
a length of 12.)

One way to do this is “brute force”.

``` r
generate_protein(6)
```

    [1] "L" "P" "L" "W" "V" "K"

``` r
generate_protein(7)
```

    [1] "I" "M" "G" "L" "V" "V" "V"

``` r
generate_protein(8)
```

    [1] "M" "P" "F" "Q" "L" "Y" "M" "G"

``` r
generate_protein(9)
```

    [1] "H" "H" "T" "Y" "H" "Q" "K" "F" "D"

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
    GSHGMG
    >7
    HLDGWYY
    >8
    AVYGSFVC
    >9
    IQETWKGLS
    >10
    PFHFGGWHQI
    >11
    CTLFQCCNGQW
    >12
    DFGCRQIMHEAV

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
    GVWLKL
    >7
    TWIGLRH
    >8
    RWSIAEPQ
    >9
    CPQWVPPGW
    >10
    QWCYHKPMII
    >11
    YRAYNDECYCP
    >12
    WLHTGCSGKRYT
