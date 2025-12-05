# Find a Gene Project
Kavi (PID: )

> 7.  Convert the alignment to FASTA if needed and use R/Bio3D to
>     calculate a sequence identity matrix.Generate a sequence identity
>     heatmap in R and include it in the report.

``` r
library(bio3d)

aln <- read.fasta("find-a-gene.fa")
idmat <- seqidentity(aln)

# widen plot margins so labels aren't cut off
par(mar = c(12, 12, 4, 4))  # bottom, left, top, right

heatmap(idmat,
        symm = TRUE,
        Rowv = NA,
        Colv = NA,
        scale = "none"
,
        cexRow = 0.9,
        cexCol = 0.9,
        margins = c(12,12),
        las = 2,         # rotate labels for readability
        main = "Sequence Identity Heatmap")
```

![](Q7-Find-a-Gene-Project_files/figure-commonmark/unnamed-chunk-1-1.png)

``` r
# restore defaults
par(mar = c(5,4,4,2))
```
