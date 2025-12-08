# Class 19
Kavi Gonur (PID: A69046927)

- [Background](#background)
- [The CMI-PB Project](#the-cmi-pb-project)
  - [Working with Dates](#working-with-dates)
- [Focus in IgG](#focus-in-igg)
  - [Differences between aP and wP?](#differences-between-ap-and-wp)
  - [Time course analysis](#time-course-analysis)
- [Time course of PT (Virulence Factor: Pertussis
  Toxin)](#time-course-of-pt-virulence-factor-pertussis-toxin)
- [System setup](#system-setup)

## Background

Pertussis (a.k.a Whooping Cough) is a highly infectious lung infection
caused by the bacteia *B. pertussis*.

The CDC tracks case numbers in the US and makes this data available
online:

> Q1. With the help of the R “addin” package datapasta assign the CDC
> pertussis case number data to a data frame called cdc and use ggplot
> to make a plot of cases numbers over time.

``` r
library(ggplot2)
```

    Warning: package 'ggplot2' was built under R version 4.5.2

``` r
cdcgraph <- ggplot(cdc) +
  aes(year, cases) +
  geom_point() +
  geom_line() +
  labs(x="Year",y="Number of cases") +
  theme_classic()

cdcgraph
```

![](class19_files/figure-commonmark/unnamed-chunk-2-1.png)

> Q2. Using the ggplot `geom_vline()` function add lines to your
> previous plot for the 1946 introduction of the wP vaccine and the 1996
> switch to aP vaccine (see example in the hint below). What do you
> notice?

``` r
cdcgraph +
  geom_vline(xintercept = 1946, linetype="dashed", color = "blue") +
  geom_vline(xintercept = 1996, linetype="dashed", color = "red") +
  geom_vline(xintercept = 2020, linetype="dashed", color = "grey") +
  geom_text(aes(x=1946, y=250000, label="wP"), 
            angle=90, vjust = -0.5, color="blue") +
  geom_text(aes(x=1996, y=250000, label="aP"), 
            angle=90, vjust = -0.5, color="red")
```

    Warning in geom_text(aes(x = 1946, y = 250000, label = "wP"), angle = 90, : All aesthetics have length 1, but the data has 102 rows.
    ℹ Please consider using `annotate()` or provide this layer with data containing
      a single row.

    Warning in geom_text(aes(x = 1996, y = 250000, label = "aP"), angle = 90, : All aesthetics have length 1, but the data has 102 rows.
    ℹ Please consider using `annotate()` or provide this layer with data containing
      a single row.

![](class19_files/figure-commonmark/unnamed-chunk-3-1.png)

> Q3. Describe what happened after the introduction of the aP vaccine?
> Do you have a possible explanation for the observed trend?

Maybe the vaccine wasn’t as effective as hoped, or there were changes in
vaccination rates, pathogen evolution, or reporting practices. Further
investigation would be needed to determine the exact cause.

## The CMI-PB Project

The CMI-PB project is a collaboration between researchers at UCSD and
the Scripps Institution of Oceanography to study the microbial
communities in the coastal waters of Southern California. The project
involves collecting water samples from various locations along the coast
and analyzing the microbial DNA using high-throughput sequencing
techniques.

They make their data aavailable via a JSON format running API. E=We can
read JSON format with the `read_json` function from the jsonlite R
package..

``` r
library(jsonlite)
subject <- read_json("https://www.cmi-pb.org/api/subject", simplifyVector = TRUE) 
head(subject, 3)
```

      subject_id infancy_vac biological_sex              ethnicity  race
    1          1          wP         Female Not Hispanic or Latino White
    2          2          wP         Female Not Hispanic or Latino White
    3          3          wP         Female                Unknown White
      year_of_birth date_of_boost      dataset
    1    1986-01-01    2016-09-12 2020_dataset
    2    1968-01-01    2019-01-28 2020_dataset
    3    1983-01-01    2016-10-10 2020_dataset

> Q4. How many aP and wP infancy vaccinated subjects are in the dataset?

``` r
table(subject$infancy_vac)
```


    aP wP 
    87 85 

``` r
subject$infancy_vac
```

      [1] "wP" "wP" "wP" "wP" "wP" "wP" "wP" "wP" "aP" "wP" "wP" "wP" "aP" "wP" "wP"
     [16] "wP" "wP" "aP" "wP" "wP" "wP" "wP" "wP" "wP" "wP" "wP" "aP" "wP" "aP" "wP"
     [31] "wP" "aP" "wP" "wP" "wP" "aP" "aP" "aP" "wP" "wP" "wP" "aP" "aP" "aP" "aP"
     [46] "aP" "aP" "aP" "aP" "aP" "aP" "aP" "aP" "aP" "aP" "aP" "aP" "aP" "aP" "aP"
     [61] "wP" "wP" "wP" "wP" "wP" "wP" "wP" "wP" "wP" "aP" "aP" "wP" "wP" "wP" "aP"
     [76] "aP" "wP" "wP" "wP" "wP" "wP" "aP" "aP" "aP" "aP" "aP" "aP" "aP" "aP" "aP"
     [91] "aP" "aP" "aP" "aP" "aP" "aP" "wP" "wP" "aP" "aP" "aP" "aP" "wP" "wP" "wP"
    [106] "aP" "aP" "wP" "wP" "aP" "wP" "aP" "aP" "wP" "aP" "aP" "aP" "aP" "aP" "wP"
    [121] "aP" "aP" "wP" "aP" "wP" "wP" "aP" "wP" "wP" "wP" "aP" "wP" "aP" "wP" "wP"
    [136] "wP" "aP" "aP" "wP" "aP" "wP" "aP" "aP" "aP" "aP" "wP" "aP" "wP" "wP" "wP"
    [151] "wP" "wP" "aP" "aP" "aP" "aP" "aP" "aP" "wP" "aP" "aP" "aP" "wP" "wP" "wP"
    [166] "aP" "aP" "wP" "aP" "wP" "wP" "wP"

> Q5. How many Male and Female subjects/patients are in the dataset?

``` r
table(subject$biological_sex)
```


    Female   Male 
       112     60 

> Q6. What is the breakdown of race and biological sex (e.g. number of
> Asian females, White males etc…)?

``` r
table(subject$race, subject$biological_sex)
```

                                               
                                                Female Male
      American Indian/Alaska Native                  0    1
      Asian                                         32   12
      Black or African American                      2    3
      More Than One Race                            15    4
      Native Hawaiian or Other Pacific Islander      1    1
      Unknown or Not Reported                       14    7
      White                                         48   32

Let’s read more tables

``` r
library(jsonlite)
specimen <- read_json("https://www.cmi-pb.org/api/v5_1/specimen", simplifyVector = TRUE)
ab_titer <- read_json("https://www.cmi-pb.org/api/v5_1/plasma_ab_titer",simplifyVector = TRUE)
```

### Working with Dates

> Q7. Using this approach determine (i) the average age of wP
> individuals, (ii) the average age of aP individuals; and (iii) are
> they significantly different?

``` r
library(lubridate)
```

    Warning: package 'lubridate' was built under R version 4.5.2


    Attaching package: 'lubridate'

    The following objects are masked from 'package:base':

        date, intersect, setdiff, union

``` r
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
subject$age <- today() - ymd(subject$year_of_birth)

# (i)
ap <- subject %>% filter(infancy_vac == "aP")
round(summary(time_length(ap$age, "years" )))
```

       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
         23      27      28      28      29      35 

``` r
 # (ii)
wp <- subject %>% filter(infancy_vac == "wP")
round(summary(time_length(wp$age, "years")))
```

       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
         23      33      35      37      40      58 

``` r
# (iii)
t.test(ap$age,wp$age)
```


        Welch Two Sample t-test

    data:  ap$age and wp$age
    t = -12.918 days, df = 104.03, p-value < 2.2e-16
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -3686.855 days -2705.535 days
    sample estimates:
    Time differences in days
    mean of x mean of y 
     10165.28  13361.47 

1)  28
2)  37
3)  yes, significantly different (p 2.2e-16)

> Q8. Determine the age of all individuals at time of boost?

``` r
subject$boost_age <- ymd(subject$date_of_boost) - ymd(subject$year_of_birth)
round(head(time_length(subject$boost_age,"years")))
```

    [1] 31 51 34 29 26 29

``` r
round(summary(time_length(subject$boost_age,"years")))
```

       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
         19      21      26      26      30      51 

> Q9a. With the help of a faceted boxplot or histogram, do you think
> these two groups are significantly different?

``` r
ggplot(subject) +
  aes(time_length(age, "year"),
      fill=as.factor(infancy_vac)) +
  geom_histogram(show.legend=FALSE) +
  facet_wrap(vars(infancy_vac), nrow=2) +
  xlab("Age in years")
```

    `stat_bin()` using `bins = 30`. Pick better value `binwidth`.

![](class19_files/figure-commonmark/unnamed-chunk-11-1.png)

I think they are significantly different!

Join (or link, or merge) using the

> Q9b. Complete the code to join specimen and subject tables to make a
> new merged data frame containing all specimen records along with their
> associated subject details: Q10. Now using the same procedure join
> meta with titer data so we can further analyze this data in terms of
> time of visit aP/wP, male/female etc.

``` r
library(dplyr)

meta <- inner_join(subject, specimen)
```

    Joining with `by = join_by(subject_id)`

``` r
head(meta)
```

      subject_id infancy_vac biological_sex              ethnicity  race
    1          1          wP         Female Not Hispanic or Latino White
    2          1          wP         Female Not Hispanic or Latino White
    3          1          wP         Female Not Hispanic or Latino White
    4          1          wP         Female Not Hispanic or Latino White
    5          1          wP         Female Not Hispanic or Latino White
    6          1          wP         Female Not Hispanic or Latino White
      year_of_birth date_of_boost      dataset        age  boost_age specimen_id
    1    1986-01-01    2016-09-12 2020_dataset 14586 days 11212 days           1
    2    1986-01-01    2016-09-12 2020_dataset 14586 days 11212 days           2
    3    1986-01-01    2016-09-12 2020_dataset 14586 days 11212 days           3
    4    1986-01-01    2016-09-12 2020_dataset 14586 days 11212 days           4
    5    1986-01-01    2016-09-12 2020_dataset 14586 days 11212 days           5
    6    1986-01-01    2016-09-12 2020_dataset 14586 days 11212 days           6
      actual_day_relative_to_boost planned_day_relative_to_boost specimen_type
    1                           -3                             0         Blood
    2                            1                             1         Blood
    3                            3                             3         Blood
    4                            7                             7         Blood
    5                           11                            14         Blood
    6                           32                            30         Blood
      visit
    1     1
    2     2
    3     3
    4     4
    5     5
    6     6

``` r
ab_data <- inner_join(meta,ab_titer)
```

    Joining with `by = join_by(specimen_id)`

``` r
head(ab_data)
```

      subject_id infancy_vac biological_sex              ethnicity  race
    1          1          wP         Female Not Hispanic or Latino White
    2          1          wP         Female Not Hispanic or Latino White
    3          1          wP         Female Not Hispanic or Latino White
    4          1          wP         Female Not Hispanic or Latino White
    5          1          wP         Female Not Hispanic or Latino White
    6          1          wP         Female Not Hispanic or Latino White
      year_of_birth date_of_boost      dataset        age  boost_age specimen_id
    1    1986-01-01    2016-09-12 2020_dataset 14586 days 11212 days           1
    2    1986-01-01    2016-09-12 2020_dataset 14586 days 11212 days           1
    3    1986-01-01    2016-09-12 2020_dataset 14586 days 11212 days           1
    4    1986-01-01    2016-09-12 2020_dataset 14586 days 11212 days           1
    5    1986-01-01    2016-09-12 2020_dataset 14586 days 11212 days           1
    6    1986-01-01    2016-09-12 2020_dataset 14586 days 11212 days           1
      actual_day_relative_to_boost planned_day_relative_to_boost specimen_type
    1                           -3                             0         Blood
    2                           -3                             0         Blood
    3                           -3                             0         Blood
    4                           -3                             0         Blood
    5                           -3                             0         Blood
    6                           -3                             0         Blood
      visit isotype is_antigen_specific antigen        MFI MFI_normalised  unit
    1     1     IgE               FALSE   Total 1110.21154       2.493425 UG/ML
    2     1     IgE               FALSE   Total 2708.91616       2.493425 IU/ML
    3     1     IgG                TRUE      PT   68.56614       3.736992 IU/ML
    4     1     IgG                TRUE     PRN  332.12718       2.602350 IU/ML
    5     1     IgG                TRUE     FHA 1887.12263      34.050956 IU/ML
    6     1     IgE                TRUE     ACT    0.10000       1.000000 IU/ML
      lower_limit_of_detection
    1                 2.096133
    2                29.170000
    3                 0.530000
    4                 6.205949
    5                 4.679535
    6                 2.816431

> Q11. How many specimens (i.e. entries in abdata) do we have for each
> isotype?

``` r
head(ab_data$isotype)
```

    [1] "IgE" "IgE" "IgG" "IgG" "IgG" "IgE"

> How many different antigens are there in the dataset?

``` r
unique(ab_data$antigen)
```

     [1] "Total"   "PT"      "PRN"     "FHA"     "ACT"     "LOS"     "FELD1"  
     [8] "BETV1"   "LOLP1"   "Measles" "PTM"     "FIM2/3"  "TT"      "DT"     
    [15] "OVA"     "PD1"    

``` r
ggplot(ab_data) +
  aes(MFI,antigen) +
  geom_boxplot() +
  theme_classic()
```

    Warning: Removed 1 row containing non-finite outside the scale range
    (`stat_boxplot()`).

![](class19_files/figure-commonmark/unnamed-chunk-16-1.png)

> Q12. What are the different \$dataset values in abdata and what do you
> notice about the number of rows for the most “recent” dataset?

``` r
table(ab_data$dataset)
```


    2020_dataset 2021_dataset 2022_dataset 2023_dataset 
           31520         8085         7301        15050 

There’s a lot more rows in the most recent dataset!

## Focus in IgG

IgG is crucial for long-term immunity and responding to bacterial and
viral infections

``` r
ab_data |>
  filter(isotype == "IgG") -> igg_data
```

Plot of antigen levels again but for IgG only

``` r
igg_dataplot <- ggplot(igg_data) +
  aes(x=MFI_normalised, y=antigen) +
  geom_boxplot() +
  theme_classic()
igg_dataplot
```

![](class19_files/figure-commonmark/unnamed-chunk-19-1.png)

### Differences between aP and wP?

We can color up by the `infancy_vac` values of “wP” or “aP”

``` r
igg_dataplot +
  aes(color=infancy_vac)
```

![](class19_files/figure-commonmark/unnamed-chunk-20-1.png)

We could “facet” by the “aP” vs “wP” column

``` r
igg_dataplot +
  aes(color=infancy_vac) +
  facet_wrap(~infancy_vac)
```

![](class19_files/figure-commonmark/unnamed-chunk-21-1.png)

### Time course analysis

> Q13. Complete the following code to make a summary boxplot of Ab titer
> levels (MFI) for all antigens:

We can use `visit` as a proxy for time here and facet our plots by this
value 1 to 8…

``` r
igg_data |>
  filter(visit %in% 1:8) |>
  ggplot() +
    aes(x=MFI_normalised,y=antigen,color=infancy_vac) +
    facet_wrap(~visit) +
    geom_boxplot() +
    theme_classic()
```

![](class19_files/figure-commonmark/unnamed-chunk-22-1.png)

> Q14. What antigens show differences in the level of IgG antibody
> titers recognizing them over time? Why these and not others?

> Q15. Filter to pull out only two specific antigens for analysis and
> create a boxplot for each. You can chose any you like. Below I picked
> a “control” antigen (“OVA”, that is not in our vaccines) and a clear
> antigen of interest (“PT”, Pertussis Toxin, one of the key virulence
> factors produced by the bacterium B. pertussis).

``` r
library(dplyr)
filter(igg_data, antigen=="PRN") %>%
  ggplot() +
  aes(MFI_normalised, col=infancy_vac) +
  geom_boxplot(show.legend = TRUE) +
  facet_wrap(vars(visit)) +
  theme_classic() +
  labs(title="PRN")
```

![](class19_files/figure-commonmark/unnamed-chunk-23-1.png)

``` r
filter(igg_data, antigen=="FIM2/3") %>%
  ggplot() +
  aes(MFI_normalised, col=infancy_vac) +
  geom_boxplot(show.legend = TRUE) +
  facet_wrap(vars(visit)) +
  theme_classic() +
  labs(title="FIM2/3")
```

![](class19_files/figure-commonmark/unnamed-chunk-23-2.png)

> Q16. What do you notice about these two antigens time courses and the
> PT data in particular?

Of the data presented in the example: PT levels overtime rise and exceed
OVA. In this dataset, FIM levels start out large but drop significantly
more than PRN (this could be becuse FIM just had a greater range of
values than PRN.)

> Q17. Do you see any clear difference in aP vs. wP responses?

PRN: Responses about the same FIM: For most part, aP \> wP.

## Time course of PT (Virulence Factor: Pertussis Toxin)

``` r
pt_2020 <- igg_data |>
  filter(antigen == "PT") |>
  filter(dataset == "2020_dataset")
pt_2021 <- igg_data |>
  filter(antigen == "PT") |>
  filter(dataset == "2021_dataset")

pt_2020
```

        subject_id infancy_vac biological_sex              ethnicity
    1            1          wP         Female Not Hispanic or Latino
    2            1          wP         Female Not Hispanic or Latino
    3            1          wP         Female Not Hispanic or Latino
    4            1          wP         Female Not Hispanic or Latino
    5            1          wP         Female Not Hispanic or Latino
    6            1          wP         Female Not Hispanic or Latino
    7            1          wP         Female Not Hispanic or Latino
    8            3          wP         Female                Unknown
    9            3          wP         Female                Unknown
    10           3          wP         Female                Unknown
    11           3          wP         Female                Unknown
    12           3          wP         Female                Unknown
    13           3          wP         Female                Unknown
    14           3          wP         Female                Unknown
    15           4          wP           Male Not Hispanic or Latino
    16           4          wP           Male Not Hispanic or Latino
    17           4          wP           Male Not Hispanic or Latino
    18           4          wP           Male Not Hispanic or Latino
    19           4          wP           Male Not Hispanic or Latino
    20           4          wP           Male Not Hispanic or Latino
    21           4          wP           Male Not Hispanic or Latino
    22           5          wP           Male Not Hispanic or Latino
    23           5          wP           Male Not Hispanic or Latino
    24           5          wP           Male Not Hispanic or Latino
    25           5          wP           Male Not Hispanic or Latino
    26           5          wP           Male Not Hispanic or Latino
    27           5          wP           Male Not Hispanic or Latino
    28           5          wP           Male Not Hispanic or Latino
    29           6          wP         Female Not Hispanic or Latino
    30           6          wP         Female Not Hispanic or Latino
    31           6          wP         Female Not Hispanic or Latino
    32           6          wP         Female Not Hispanic or Latino
    33           6          wP         Female Not Hispanic or Latino
    34           6          wP         Female Not Hispanic or Latino
    35           6          wP         Female Not Hispanic or Latino
    36           7          wP         Female     Hispanic or Latino
    37           7          wP         Female     Hispanic or Latino
    38           7          wP         Female     Hispanic or Latino
    39           7          wP         Female     Hispanic or Latino
    40           7          wP         Female     Hispanic or Latino
    41           7          wP         Female     Hispanic or Latino
    42           7          wP         Female     Hispanic or Latino
    43           9          aP           Male Not Hispanic or Latino
    44           9          aP           Male Not Hispanic or Latino
    45           9          aP           Male Not Hispanic or Latino
    46           9          aP           Male Not Hispanic or Latino
    47           9          aP           Male Not Hispanic or Latino
    48           9          aP           Male Not Hispanic or Latino
    49           9          aP           Male Not Hispanic or Latino
    50          10          wP         Female Not Hispanic or Latino
    51          10          wP         Female Not Hispanic or Latino
    52          10          wP         Female Not Hispanic or Latino
    53          10          wP         Female Not Hispanic or Latino
    54          10          wP         Female Not Hispanic or Latino
    55          10          wP         Female Not Hispanic or Latino
    56          11          wP         Female     Hispanic or Latino
    57          11          wP         Female     Hispanic or Latino
    58          11          wP         Female     Hispanic or Latino
    59          11          wP         Female     Hispanic or Latino
    60          11          wP         Female     Hispanic or Latino
    61          11          wP         Female     Hispanic or Latino
    62          11          wP         Female     Hispanic or Latino
    63          12          wP           Male Not Hispanic or Latino
    64          12          wP           Male Not Hispanic or Latino
    65          12          wP           Male Not Hispanic or Latino
    66          12          wP           Male Not Hispanic or Latino
    67          12          wP           Male Not Hispanic or Latino
    68          12          wP           Male Not Hispanic or Latino
    69          13          aP           Male Not Hispanic or Latino
    70          13          aP           Male Not Hispanic or Latino
    71          13          aP           Male Not Hispanic or Latino
    72          13          aP           Male Not Hispanic or Latino
    73          13          aP           Male Not Hispanic or Latino
    74          13          aP           Male Not Hispanic or Latino
    75          13          aP           Male Not Hispanic or Latino
    76          14          wP           Male Not Hispanic or Latino
    77          14          wP           Male Not Hispanic or Latino
    78          14          wP           Male Not Hispanic or Latino
    79          14          wP           Male Not Hispanic or Latino
    80          14          wP           Male Not Hispanic or Latino
    81          15          wP           Male Not Hispanic or Latino
    82          15          wP           Male Not Hispanic or Latino
    83          15          wP           Male Not Hispanic or Latino
    84          15          wP           Male Not Hispanic or Latino
    85          15          wP           Male Not Hispanic or Latino
    86          15          wP           Male Not Hispanic or Latino
    87          15          wP           Male Not Hispanic or Latino
    88          16          wP         Female     Hispanic or Latino
    89          16          wP         Female     Hispanic or Latino
    90          16          wP         Female     Hispanic or Latino
    91          16          wP         Female     Hispanic or Latino
    92          16          wP         Female     Hispanic or Latino
    93          16          wP         Female     Hispanic or Latino
    94          16          wP         Female     Hispanic or Latino
    95          17          wP         Female     Hispanic or Latino
    96          17          wP         Female     Hispanic or Latino
    97          17          wP         Female     Hispanic or Latino
    98          17          wP         Female     Hispanic or Latino
    99          17          wP         Female     Hispanic or Latino
    100         17          wP         Female     Hispanic or Latino
    101         17          wP         Female     Hispanic or Latino
    102         18          aP         Female     Hispanic or Latino
    103         18          aP         Female     Hispanic or Latino
    104         18          aP         Female     Hispanic or Latino
    105         18          aP         Female     Hispanic or Latino
    106         18          aP         Female     Hispanic or Latino
    107         18          aP         Female     Hispanic or Latino
    108         18          aP         Female     Hispanic or Latino
    109         19          wP           Male Not Hispanic or Latino
    110         19          wP           Male Not Hispanic or Latino
    111         19          wP           Male Not Hispanic or Latino
    112         19          wP           Male Not Hispanic or Latino
    113         19          wP           Male Not Hispanic or Latino
    114         19          wP           Male Not Hispanic or Latino
    115         19          wP           Male Not Hispanic or Latino
    116         20          wP         Female Not Hispanic or Latino
    117         20          wP         Female Not Hispanic or Latino
    118         20          wP         Female Not Hispanic or Latino
    119         20          wP         Female Not Hispanic or Latino
    120         20          wP         Female Not Hispanic or Latino
    121         20          wP         Female Not Hispanic or Latino
    122         20          wP         Female Not Hispanic or Latino
    123         21          wP           Male Not Hispanic or Latino
    124         21          wP           Male Not Hispanic or Latino
    125         21          wP           Male Not Hispanic or Latino
    126         21          wP           Male Not Hispanic or Latino
    127         21          wP           Male Not Hispanic or Latino
    128         21          wP           Male Not Hispanic or Latino
    129         21          wP           Male Not Hispanic or Latino
    130         22          wP         Female Not Hispanic or Latino
    131         22          wP         Female Not Hispanic or Latino
    132         22          wP         Female Not Hispanic or Latino
    133         22          wP         Female Not Hispanic or Latino
    134         22          wP         Female Not Hispanic or Latino
    135         22          wP         Female Not Hispanic or Latino
    136         22          wP         Female Not Hispanic or Latino
    137         23          wP         Female Not Hispanic or Latino
    138         23          wP         Female Not Hispanic or Latino
    139         23          wP         Female Not Hispanic or Latino
    140         23          wP         Female Not Hispanic or Latino
    141         23          wP         Female Not Hispanic or Latino
    142         23          wP         Female Not Hispanic or Latino
    143         23          wP         Female Not Hispanic or Latino
    144         24          wP         Female Not Hispanic or Latino
    145         24          wP         Female Not Hispanic or Latino
    146         24          wP         Female Not Hispanic or Latino
    147         24          wP         Female Not Hispanic or Latino
    148         24          wP         Female Not Hispanic or Latino
    149         24          wP         Female Not Hispanic or Latino
    150         24          wP         Female Not Hispanic or Latino
    151         25          wP         Female Not Hispanic or Latino
    152         25          wP         Female Not Hispanic or Latino
    153         25          wP         Female Not Hispanic or Latino
    154         25          wP         Female Not Hispanic or Latino
    155         25          wP         Female Not Hispanic or Latino
    156         25          wP         Female Not Hispanic or Latino
    157         25          wP         Female Not Hispanic or Latino
    158         26          wP         Female     Hispanic or Latino
    159         26          wP         Female     Hispanic or Latino
    160         26          wP         Female     Hispanic or Latino
    161         26          wP         Female     Hispanic or Latino
    162         26          wP         Female     Hispanic or Latino
    163         26          wP         Female     Hispanic or Latino
    164         26          wP         Female     Hispanic or Latino
    165         27          aP         Female Not Hispanic or Latino
    166         27          aP         Female Not Hispanic or Latino
    167         27          aP         Female Not Hispanic or Latino
    168         27          aP         Female Not Hispanic or Latino
    169         27          aP         Female Not Hispanic or Latino
    170         27          aP         Female Not Hispanic or Latino
    171         27          aP         Female Not Hispanic or Latino
    172         28          wP           Male                Unknown
    173         28          wP           Male                Unknown
    174         28          wP           Male                Unknown
    175         28          wP           Male                Unknown
    176         28          wP           Male                Unknown
    177         28          wP           Male                Unknown
    178         28          wP           Male                Unknown
    179         29          aP           Male     Hispanic or Latino
    180         29          aP           Male     Hispanic or Latino
    181         29          aP           Male     Hispanic or Latino
    182         29          aP           Male     Hispanic or Latino
    183         29          aP           Male     Hispanic or Latino
    184         29          aP           Male     Hispanic or Latino
    185         29          aP           Male     Hispanic or Latino
    186         30          wP         Female     Hispanic or Latino
    187         30          wP         Female     Hispanic or Latino
    188         30          wP         Female     Hispanic or Latino
    189         30          wP         Female     Hispanic or Latino
    190         30          wP         Female     Hispanic or Latino
    191         30          wP         Female     Hispanic or Latino
    192         30          wP         Female     Hispanic or Latino
    193         31          wP         Female Not Hispanic or Latino
    194         31          wP         Female Not Hispanic or Latino
    195         31          wP         Female Not Hispanic or Latino
    196         31          wP         Female Not Hispanic or Latino
    197         31          wP         Female Not Hispanic or Latino
    198         31          wP         Female Not Hispanic or Latino
    199         31          wP         Female Not Hispanic or Latino
    200         32          aP           Male Not Hispanic or Latino
    201         32          aP           Male Not Hispanic or Latino
    202         32          aP           Male Not Hispanic or Latino
    203         32          aP           Male Not Hispanic or Latino
    204         32          aP           Male Not Hispanic or Latino
    205         32          aP           Male Not Hispanic or Latino
    206         32          aP           Male Not Hispanic or Latino
    207         33          wP           Male     Hispanic or Latino
    208         33          wP           Male     Hispanic or Latino
    209         33          wP           Male     Hispanic or Latino
    210         33          wP           Male     Hispanic or Latino
    211         33          wP           Male     Hispanic or Latino
    212         33          wP           Male     Hispanic or Latino
    213         33          wP           Male     Hispanic or Latino
    214         34          wP         Female     Hispanic or Latino
    215         34          wP         Female     Hispanic or Latino
    216         34          wP         Female     Hispanic or Latino
    217         34          wP         Female     Hispanic or Latino
    218         34          wP         Female     Hispanic or Latino
    219         34          wP         Female     Hispanic or Latino
    220         34          wP         Female     Hispanic or Latino
    221         35          wP           Male                Unknown
    222         35          wP           Male                Unknown
    223         35          wP           Male                Unknown
    224         35          wP           Male                Unknown
    225         35          wP           Male                Unknown
    226         35          wP           Male                Unknown
    227         35          wP           Male                Unknown
    228         36          aP         Female     Hispanic or Latino
    229         36          aP         Female     Hispanic or Latino
    230         36          aP         Female     Hispanic or Latino
    231         36          aP         Female     Hispanic or Latino
    232         36          aP         Female     Hispanic or Latino
    233         36          aP         Female     Hispanic or Latino
    234         37          aP         Female Not Hispanic or Latino
    235         37          aP         Female Not Hispanic or Latino
    236         37          aP         Female Not Hispanic or Latino
    237         37          aP         Female Not Hispanic or Latino
    238         37          aP         Female Not Hispanic or Latino
    239         38          aP         Female Not Hispanic or Latino
    240         38          aP         Female Not Hispanic or Latino
    241         38          aP         Female Not Hispanic or Latino
    242         38          aP         Female Not Hispanic or Latino
    243         38          aP         Female Not Hispanic or Latino
    244         38          aP         Female Not Hispanic or Latino
    245         38          aP         Female Not Hispanic or Latino
    246         39          wP         Female Not Hispanic or Latino
    247         39          wP         Female Not Hispanic or Latino
    248         39          wP         Female Not Hispanic or Latino
    249         39          wP         Female Not Hispanic or Latino
    250         39          wP         Female Not Hispanic or Latino
    251         39          wP         Female Not Hispanic or Latino
    252         39          wP         Female Not Hispanic or Latino
    253         40          wP         Female Not Hispanic or Latino
    254         40          wP         Female Not Hispanic or Latino
    255         40          wP         Female Not Hispanic or Latino
    256         40          wP         Female Not Hispanic or Latino
    257         40          wP         Female Not Hispanic or Latino
    258         40          wP         Female Not Hispanic or Latino
    259         40          wP         Female Not Hispanic or Latino
    260         41          wP           Male Not Hispanic or Latino
    261         41          wP           Male Not Hispanic or Latino
    262         41          wP           Male Not Hispanic or Latino
    263         41          wP           Male Not Hispanic or Latino
    264         41          wP           Male Not Hispanic or Latino
    265         41          wP           Male Not Hispanic or Latino
    266         41          wP           Male Not Hispanic or Latino
    267         42          aP         Female Not Hispanic or Latino
    268         42          aP         Female Not Hispanic or Latino
    269         42          aP         Female Not Hispanic or Latino
    270         42          aP         Female Not Hispanic or Latino
    271         42          aP         Female Not Hispanic or Latino
    272         42          aP         Female Not Hispanic or Latino
    273         42          aP         Female Not Hispanic or Latino
    274         43          aP         Female Not Hispanic or Latino
    275         43          aP         Female Not Hispanic or Latino
    276         43          aP         Female Not Hispanic or Latino
    277         43          aP         Female Not Hispanic or Latino
    278         43          aP         Female Not Hispanic or Latino
    279         43          aP         Female Not Hispanic or Latino
    280         43          aP         Female Not Hispanic or Latino
    281         44          aP         Female     Hispanic or Latino
    282         44          aP         Female     Hispanic or Latino
    283         44          aP         Female     Hispanic or Latino
    284         44          aP         Female     Hispanic or Latino
    285         44          aP         Female     Hispanic or Latino
    286         44          aP         Female     Hispanic or Latino
    287         44          aP         Female     Hispanic or Latino
    288         45          aP         Female Not Hispanic or Latino
    289         45          aP         Female Not Hispanic or Latino
    290         45          aP         Female Not Hispanic or Latino
    291         45          aP         Female Not Hispanic or Latino
    292         45          aP         Female Not Hispanic or Latino
    293         45          aP         Female Not Hispanic or Latino
    294         46          aP         Female Not Hispanic or Latino
    295         46          aP         Female Not Hispanic or Latino
    296         46          aP         Female Not Hispanic or Latino
    297         46          aP         Female Not Hispanic or Latino
    298         46          aP         Female Not Hispanic or Latino
    299         47          aP         Female Not Hispanic or Latino
    300         47          aP         Female Not Hispanic or Latino
    301         47          aP         Female Not Hispanic or Latino
    302         47          aP         Female Not Hispanic or Latino
    303         47          aP         Female Not Hispanic or Latino
    304         47          aP         Female Not Hispanic or Latino
    305         47          aP         Female Not Hispanic or Latino
    306         48          aP         Female Not Hispanic or Latino
    307         48          aP         Female Not Hispanic or Latino
    308         48          aP         Female Not Hispanic or Latino
    309         48          aP         Female Not Hispanic or Latino
    310         48          aP         Female Not Hispanic or Latino
    311         48          aP         Female Not Hispanic or Latino
    312         48          aP         Female Not Hispanic or Latino
    313         49          aP         Female Not Hispanic or Latino
    314         49          aP         Female Not Hispanic or Latino
    315         49          aP         Female Not Hispanic or Latino
    316         49          aP         Female Not Hispanic or Latino
    317         49          aP         Female Not Hispanic or Latino
    318         49          aP         Female Not Hispanic or Latino
    319         49          aP         Female Not Hispanic or Latino
    320         50          aP         Female Not Hispanic or Latino
    321         50          aP         Female Not Hispanic or Latino
    322         50          aP         Female Not Hispanic or Latino
    323         50          aP         Female Not Hispanic or Latino
    324         50          aP         Female Not Hispanic or Latino
    325         50          aP         Female Not Hispanic or Latino
    326         50          aP         Female Not Hispanic or Latino
    327         51          aP           Male Not Hispanic or Latino
    328         51          aP           Male Not Hispanic or Latino
    329         51          aP           Male Not Hispanic or Latino
    330         51          aP           Male Not Hispanic or Latino
    331         51          aP           Male Not Hispanic or Latino
    332         52          aP           Male Not Hispanic or Latino
    333         52          aP           Male Not Hispanic or Latino
    334         52          aP           Male Not Hispanic or Latino
    335         52          aP           Male Not Hispanic or Latino
    336         52          aP           Male Not Hispanic or Latino
    337         52          aP           Male Not Hispanic or Latino
    338         52          aP           Male Not Hispanic or Latino
    339         53          aP         Female     Hispanic or Latino
    340         53          aP         Female     Hispanic or Latino
    341         53          aP         Female     Hispanic or Latino
    342         53          aP         Female     Hispanic or Latino
    343         53          aP         Female     Hispanic or Latino
    344         53          aP         Female     Hispanic or Latino
    345         53          aP         Female     Hispanic or Latino
    346         54          aP         Female Not Hispanic or Latino
    347         54          aP         Female Not Hispanic or Latino
    348         54          aP         Female Not Hispanic or Latino
    349         54          aP         Female Not Hispanic or Latino
    350         54          aP         Female Not Hispanic or Latino
    351         54          aP         Female Not Hispanic or Latino
    352         54          aP         Female Not Hispanic or Latino
    353         55          aP         Female Not Hispanic or Latino
    354         55          aP         Female Not Hispanic or Latino
    355         55          aP         Female Not Hispanic or Latino
    356         55          aP         Female Not Hispanic or Latino
    357         55          aP         Female Not Hispanic or Latino
    358         55          aP         Female Not Hispanic or Latino
    359         55          aP         Female Not Hispanic or Latino
    360         56          aP         Female Not Hispanic or Latino
    361         56          aP         Female Not Hispanic or Latino
    362         56          aP         Female Not Hispanic or Latino
    363         56          aP         Female Not Hispanic or Latino
    364         56          aP         Female Not Hispanic or Latino
    365         56          aP         Female Not Hispanic or Latino
    366         56          aP         Female Not Hispanic or Latino
    367         57          aP         Female Not Hispanic or Latino
    368         57          aP         Female Not Hispanic or Latino
    369         57          aP         Female Not Hispanic or Latino
    370         57          aP         Female Not Hispanic or Latino
    371         57          aP         Female Not Hispanic or Latino
    372         57          aP         Female Not Hispanic or Latino
    373         57          aP         Female Not Hispanic or Latino
    374         58          aP         Female     Hispanic or Latino
    375         58          aP         Female     Hispanic or Latino
    376         58          aP         Female     Hispanic or Latino
    377         58          aP         Female     Hispanic or Latino
    378         58          aP         Female     Hispanic or Latino
    379         58          aP         Female     Hispanic or Latino
    380         58          aP         Female     Hispanic or Latino
    381         59          aP         Female     Hispanic or Latino
    382         59          aP         Female     Hispanic or Latino
    383         59          aP         Female     Hispanic or Latino
    384         59          aP         Female     Hispanic or Latino
    385         59          aP         Female     Hispanic or Latino
    386         59          aP         Female     Hispanic or Latino
    387         59          aP         Female     Hispanic or Latino
    388         60          aP           Male     Hispanic or Latino
    389         60          aP           Male     Hispanic or Latino
    390         60          aP           Male     Hispanic or Latino
    391         60          aP           Male     Hispanic or Latino
    392         60          aP           Male     Hispanic or Latino
    393         60          aP           Male     Hispanic or Latino
    394         60          aP           Male     Hispanic or Latino
                                             race year_of_birth date_of_boost
    1                                       White    1986-01-01    2016-09-12
    2                                       White    1986-01-01    2016-09-12
    3                                       White    1986-01-01    2016-09-12
    4                                       White    1986-01-01    2016-09-12
    5                                       White    1986-01-01    2016-09-12
    6                                       White    1986-01-01    2016-09-12
    7                                       White    1986-01-01    2016-09-12
    8                                       White    1983-01-01    2016-10-10
    9                                       White    1983-01-01    2016-10-10
    10                                      White    1983-01-01    2016-10-10
    11                                      White    1983-01-01    2016-10-10
    12                                      White    1983-01-01    2016-10-10
    13                                      White    1983-01-01    2016-10-10
    14                                      White    1983-01-01    2016-10-10
    15                                      Asian    1988-01-01    2016-08-29
    16                                      Asian    1988-01-01    2016-08-29
    17                                      Asian    1988-01-01    2016-08-29
    18                                      Asian    1988-01-01    2016-08-29
    19                                      Asian    1988-01-01    2016-08-29
    20                                      Asian    1988-01-01    2016-08-29
    21                                      Asian    1988-01-01    2016-08-29
    22                                      Asian    1991-01-01    2016-08-29
    23                                      Asian    1991-01-01    2016-08-29
    24                                      Asian    1991-01-01    2016-08-29
    25                                      Asian    1991-01-01    2016-08-29
    26                                      Asian    1991-01-01    2016-08-29
    27                                      Asian    1991-01-01    2016-08-29
    28                                      Asian    1991-01-01    2016-08-29
    29                                      White    1988-01-01    2016-10-10
    30                                      White    1988-01-01    2016-10-10
    31                                      White    1988-01-01    2016-10-10
    32                                      White    1988-01-01    2016-10-10
    33                                      White    1988-01-01    2016-10-10
    34                                      White    1988-01-01    2016-10-10
    35                                      White    1988-01-01    2016-10-10
    36                         More Than One Race    1981-01-01    2016-11-07
    37                         More Than One Race    1981-01-01    2016-11-07
    38                         More Than One Race    1981-01-01    2016-11-07
    39                         More Than One Race    1981-01-01    2016-11-07
    40                         More Than One Race    1981-01-01    2016-11-07
    41                         More Than One Race    1981-01-01    2016-11-07
    42                         More Than One Race    1981-01-01    2016-11-07
    43                                      Asian    1996-01-01    2016-07-25
    44                                      Asian    1996-01-01    2016-07-25
    45                                      Asian    1996-01-01    2016-07-25
    46                                      Asian    1996-01-01    2016-07-25
    47                                      Asian    1996-01-01    2016-07-25
    48                                      Asian    1996-01-01    2016-07-25
    49                                      Asian    1996-01-01    2016-07-25
    50                                      Asian    1982-01-01    2016-07-25
    51                                      Asian    1982-01-01    2016-07-25
    52                                      Asian    1982-01-01    2016-07-25
    53                                      Asian    1982-01-01    2016-07-25
    54                                      Asian    1982-01-01    2016-07-25
    55                                      Asian    1982-01-01    2016-07-25
    56                    Unknown or Not Reported    1986-01-01    2016-08-29
    57                    Unknown or Not Reported    1986-01-01    2016-08-29
    58                    Unknown or Not Reported    1986-01-01    2016-08-29
    59                    Unknown or Not Reported    1986-01-01    2016-08-29
    60                    Unknown or Not Reported    1986-01-01    2016-08-29
    61                    Unknown or Not Reported    1986-01-01    2016-08-29
    62                    Unknown or Not Reported    1986-01-01    2016-08-29
    63                                      Asian    1982-01-01    2016-07-25
    64                                      Asian    1982-01-01    2016-07-25
    65                                      Asian    1982-01-01    2016-07-25
    66                                      Asian    1982-01-01    2016-07-25
    67                                      Asian    1982-01-01    2016-07-25
    68                                      Asian    1982-01-01    2016-07-25
    69                                      White    1997-01-01    2016-07-25
    70                                      White    1997-01-01    2016-07-25
    71                                      White    1997-01-01    2016-07-25
    72                                      White    1997-01-01    2016-07-25
    73                                      White    1997-01-01    2016-07-25
    74                                      White    1997-01-01    2016-07-25
    75                                      White    1997-01-01    2016-07-25
    76                                      White    1993-01-01    2016-08-15
    77                                      White    1993-01-01    2016-08-15
    78                                      White    1993-01-01    2016-08-15
    79                                      White    1993-01-01    2016-08-15
    80                                      White    1993-01-01    2016-08-15
    81                                      Asian    1989-01-01    2016-08-15
    82                                      Asian    1989-01-01    2016-08-15
    83                                      Asian    1989-01-01    2016-08-15
    84                                      Asian    1989-01-01    2016-08-15
    85                                      Asian    1989-01-01    2016-08-15
    86                                      Asian    1989-01-01    2016-08-15
    87                                      Asian    1989-01-01    2016-08-15
    88                    Unknown or Not Reported    1987-01-01    2016-07-25
    89                    Unknown or Not Reported    1987-01-01    2016-07-25
    90                    Unknown or Not Reported    1987-01-01    2016-07-25
    91                    Unknown or Not Reported    1987-01-01    2016-07-25
    92                    Unknown or Not Reported    1987-01-01    2016-07-25
    93                    Unknown or Not Reported    1987-01-01    2016-07-25
    94                    Unknown or Not Reported    1987-01-01    2016-07-25
    95                                      White    1980-01-01    2016-09-12
    96                                      White    1980-01-01    2016-09-12
    97                                      White    1980-01-01    2016-09-12
    98                                      White    1980-01-01    2016-09-12
    99                                      White    1980-01-01    2016-09-12
    100                                     White    1980-01-01    2016-09-12
    101                                     White    1980-01-01    2016-09-12
    102                   Unknown or Not Reported    1997-01-01    2016-08-29
    103                   Unknown or Not Reported    1997-01-01    2016-08-29
    104                   Unknown or Not Reported    1997-01-01    2016-08-29
    105                   Unknown or Not Reported    1997-01-01    2016-08-29
    106                   Unknown or Not Reported    1997-01-01    2016-08-29
    107                   Unknown or Not Reported    1997-01-01    2016-08-29
    108                   Unknown or Not Reported    1997-01-01    2016-08-29
    109                                     Asian    1994-01-01    2016-09-26
    110                                     Asian    1994-01-01    2016-09-26
    111                                     Asian    1994-01-01    2016-09-26
    112                                     Asian    1994-01-01    2016-09-26
    113                                     Asian    1994-01-01    2016-09-26
    114                                     Asian    1994-01-01    2016-09-26
    115                                     Asian    1994-01-01    2016-09-26
    116                                     White    1981-01-01    2016-08-29
    117                                     White    1981-01-01    2016-08-29
    118                                     White    1981-01-01    2016-08-29
    119                                     White    1981-01-01    2016-08-29
    120                                     White    1981-01-01    2016-08-29
    121                                     White    1981-01-01    2016-08-29
    122                                     White    1981-01-01    2016-08-29
    123                                     White    1983-01-01    2016-08-29
    124                                     White    1983-01-01    2016-08-29
    125                                     White    1983-01-01    2016-08-29
    126                                     White    1983-01-01    2016-08-29
    127                                     White    1983-01-01    2016-08-29
    128                                     White    1983-01-01    2016-08-29
    129                                     White    1983-01-01    2016-08-29
    130                                     White    1985-01-01    2016-08-29
    131                                     White    1985-01-01    2016-08-29
    132                                     White    1985-01-01    2016-08-29
    133                                     White    1985-01-01    2016-08-29
    134                                     White    1985-01-01    2016-08-29
    135                                     White    1985-01-01    2016-08-29
    136                                     White    1985-01-01    2016-08-29
    137                                     White    1991-01-01    2016-09-26
    138                                     White    1991-01-01    2016-09-26
    139                                     White    1991-01-01    2016-09-26
    140                                     White    1991-01-01    2016-09-26
    141                                     White    1991-01-01    2016-09-26
    142                                     White    1991-01-01    2016-09-26
    143                                     White    1991-01-01    2016-09-26
    144                                     Asian    1992-01-01    2016-09-13
    145                                     Asian    1992-01-01    2016-09-13
    146                                     Asian    1992-01-01    2016-09-13
    147                                     Asian    1992-01-01    2016-09-13
    148                                     Asian    1992-01-01    2016-09-13
    149                                     Asian    1992-01-01    2016-09-13
    150                                     Asian    1992-01-01    2016-09-13
    151                 Black or African American    1988-01-01    2016-09-13
    152                 Black or African American    1988-01-01    2016-09-13
    153                 Black or African American    1988-01-01    2016-09-13
    154                 Black or African American    1988-01-01    2016-09-13
    155                 Black or African American    1988-01-01    2016-09-13
    156                 Black or African American    1988-01-01    2016-09-13
    157                 Black or African American    1988-01-01    2016-09-13
    158                   Unknown or Not Reported    1983-01-01    2016-09-26
    159                   Unknown or Not Reported    1983-01-01    2016-09-26
    160                   Unknown or Not Reported    1983-01-01    2016-09-26
    161                   Unknown or Not Reported    1983-01-01    2016-09-26
    162                   Unknown or Not Reported    1983-01-01    2016-09-26
    163                   Unknown or Not Reported    1983-01-01    2016-09-26
    164                   Unknown or Not Reported    1983-01-01    2016-09-26
    165                                     Asian    1997-01-01    2016-09-26
    166                                     Asian    1997-01-01    2016-09-26
    167                                     Asian    1997-01-01    2016-09-26
    168                                     Asian    1997-01-01    2016-09-26
    169                                     Asian    1997-01-01    2016-09-26
    170                                     Asian    1997-01-01    2016-09-26
    171                                     Asian    1997-01-01    2016-09-26
    172                   Unknown or Not Reported    1982-01-01    2016-09-26
    173                   Unknown or Not Reported    1982-01-01    2016-09-26
    174                   Unknown or Not Reported    1982-01-01    2016-09-26
    175                   Unknown or Not Reported    1982-01-01    2016-09-26
    176                   Unknown or Not Reported    1982-01-01    2016-09-26
    177                   Unknown or Not Reported    1982-01-01    2016-09-26
    178                   Unknown or Not Reported    1982-01-01    2016-09-26
    179                                     White    1997-01-01    2016-09-26
    180                                     White    1997-01-01    2016-09-26
    181                                     White    1997-01-01    2016-09-26
    182                                     White    1997-01-01    2016-09-26
    183                                     White    1997-01-01    2016-09-26
    184                                     White    1997-01-01    2016-09-26
    185                                     White    1997-01-01    2016-09-26
    186                                     White    1988-01-01    2016-09-26
    187                                     White    1988-01-01    2016-09-26
    188                                     White    1988-01-01    2016-09-26
    189                                     White    1988-01-01    2016-09-26
    190                                     White    1988-01-01    2016-09-26
    191                                     White    1988-01-01    2016-09-26
    192                                     White    1988-01-01    2016-09-26
    193                                     Asian    1989-01-01    2016-09-26
    194                                     Asian    1989-01-01    2016-09-26
    195                                     Asian    1989-01-01    2016-09-26
    196                                     Asian    1989-01-01    2016-09-26
    197                                     Asian    1989-01-01    2016-09-26
    198                                     Asian    1989-01-01    2016-09-26
    199                                     Asian    1989-01-01    2016-09-26
    200 Native Hawaiian or Other Pacific Islander    1997-01-01    2016-10-24
    201 Native Hawaiian or Other Pacific Islander    1997-01-01    2016-10-24
    202 Native Hawaiian or Other Pacific Islander    1997-01-01    2016-10-24
    203 Native Hawaiian or Other Pacific Islander    1997-01-01    2016-10-24
    204 Native Hawaiian or Other Pacific Islander    1997-01-01    2016-10-24
    205 Native Hawaiian or Other Pacific Islander    1997-01-01    2016-10-24
    206 Native Hawaiian or Other Pacific Islander    1997-01-01    2016-10-24
    207                        More Than One Race    1990-01-01    2016-10-10
    208                        More Than One Race    1990-01-01    2016-10-10
    209                        More Than One Race    1990-01-01    2016-10-10
    210                        More Than One Race    1990-01-01    2016-10-10
    211                        More Than One Race    1990-01-01    2016-10-10
    212                        More Than One Race    1990-01-01    2016-10-10
    213                        More Than One Race    1990-01-01    2016-10-10
    214                   Unknown or Not Reported    1983-01-01    2016-10-24
    215                   Unknown or Not Reported    1983-01-01    2016-10-24
    216                   Unknown or Not Reported    1983-01-01    2016-10-24
    217                   Unknown or Not Reported    1983-01-01    2016-10-24
    218                   Unknown or Not Reported    1983-01-01    2016-10-24
    219                   Unknown or Not Reported    1983-01-01    2016-10-24
    220                   Unknown or Not Reported    1983-01-01    2016-10-24
    221                                     White    1991-01-01    2016-10-10
    222                                     White    1991-01-01    2016-10-10
    223                                     White    1991-01-01    2016-10-10
    224                                     White    1991-01-01    2016-10-10
    225                                     White    1991-01-01    2016-10-10
    226                                     White    1991-01-01    2016-10-10
    227                                     White    1991-01-01    2016-10-10
    228                                     White    1997-01-01    2016-10-24
    229                                     White    1997-01-01    2016-10-24
    230                                     White    1997-01-01    2016-10-24
    231                                     White    1997-01-01    2016-10-24
    232                                     White    1997-01-01    2016-10-24
    233                                     White    1997-01-01    2016-10-24
    234                        More Than One Race    1998-01-01    2016-11-07
    235                        More Than One Race    1998-01-01    2016-11-07
    236                        More Than One Race    1998-01-01    2016-11-07
    237                        More Than One Race    1998-01-01    2016-11-07
    238                        More Than One Race    1998-01-01    2016-11-07
    239                                     White    1997-01-01    2016-10-24
    240                                     White    1997-01-01    2016-10-24
    241                                     White    1997-01-01    2016-10-24
    242                                     White    1997-01-01    2016-10-24
    243                                     White    1997-01-01    2016-10-24
    244                                     White    1997-01-01    2016-10-24
    245                                     White    1997-01-01    2016-10-24
    246                                     White    1985-01-01    2016-10-24
    247                                     White    1985-01-01    2016-10-24
    248                                     White    1985-01-01    2016-10-24
    249                                     White    1985-01-01    2016-10-24
    250                                     White    1985-01-01    2016-10-24
    251                                     White    1985-01-01    2016-10-24
    252                                     White    1985-01-01    2016-10-24
    253                                     Asian    1994-01-01    2016-10-24
    254                                     Asian    1994-01-01    2016-10-24
    255                                     Asian    1994-01-01    2016-10-24
    256                                     Asian    1994-01-01    2016-10-24
    257                                     Asian    1994-01-01    2016-10-24
    258                                     Asian    1994-01-01    2016-10-24
    259                                     Asian    1994-01-01    2016-10-24
    260                                     White    1985-01-01    2016-11-07
    261                                     White    1985-01-01    2016-11-07
    262                                     White    1985-01-01    2016-11-07
    263                                     White    1985-01-01    2016-11-07
    264                                     White    1985-01-01    2016-11-07
    265                                     White    1985-01-01    2016-11-07
    266                                     White    1985-01-01    2016-11-07
    267                                     Asian    1997-01-01    2016-11-07
    268                                     Asian    1997-01-01    2016-11-07
    269                                     Asian    1997-01-01    2016-11-07
    270                                     Asian    1997-01-01    2016-11-07
    271                                     Asian    1997-01-01    2016-11-07
    272                                     Asian    1997-01-01    2016-11-07
    273                                     Asian    1997-01-01    2016-11-07
    274                        More Than One Race    1998-01-01    2016-11-07
    275                        More Than One Race    1998-01-01    2016-11-07
    276                        More Than One Race    1998-01-01    2016-11-07
    277                        More Than One Race    1998-01-01    2016-11-07
    278                        More Than One Race    1998-01-01    2016-11-07
    279                        More Than One Race    1998-01-01    2016-11-07
    280                        More Than One Race    1998-01-01    2016-11-07
    281                        More Than One Race    1998-01-01    2016-11-07
    282                        More Than One Race    1998-01-01    2016-11-07
    283                        More Than One Race    1998-01-01    2016-11-07
    284                        More Than One Race    1998-01-01    2016-11-07
    285                        More Than One Race    1998-01-01    2016-11-07
    286                        More Than One Race    1998-01-01    2016-11-07
    287                        More Than One Race    1998-01-01    2016-11-07
    288                                     Asian    1997-01-01    2016-11-28
    289                                     Asian    1997-01-01    2016-11-28
    290                                     Asian    1997-01-01    2016-11-28
    291                                     Asian    1997-01-01    2016-11-28
    292                                     Asian    1997-01-01    2016-11-28
    293                                     Asian    1997-01-01    2016-11-28
    294                   Unknown or Not Reported    1998-01-01    2016-11-07
    295                   Unknown or Not Reported    1998-01-01    2016-11-07
    296                   Unknown or Not Reported    1998-01-01    2016-11-07
    297                   Unknown or Not Reported    1998-01-01    2016-11-07
    298                   Unknown or Not Reported    1998-01-01    2016-11-07
    299                                     White    1996-01-01    2016-11-28
    300                                     White    1996-01-01    2016-11-28
    301                                     White    1996-01-01    2016-11-28
    302                                     White    1996-01-01    2016-11-28
    303                                     White    1996-01-01    2016-11-28
    304                                     White    1996-01-01    2016-11-28
    305                                     White    1996-01-01    2016-11-28
    306                                     White    1998-01-01    2017-01-17
    307                                     White    1998-01-01    2017-01-17
    308                                     White    1998-01-01    2017-01-17
    309                                     White    1998-01-01    2017-01-17
    310                                     White    1998-01-01    2017-01-17
    311                                     White    1998-01-01    2017-01-17
    312                                     White    1998-01-01    2017-01-17
    313                                     White    1997-01-01    2017-01-17
    314                                     White    1997-01-01    2017-01-17
    315                                     White    1997-01-01    2017-01-17
    316                                     White    1997-01-01    2017-01-17
    317                                     White    1997-01-01    2017-01-17
    318                                     White    1997-01-01    2017-01-17
    319                                     White    1997-01-01    2017-01-17
    320                                     Asian    1997-01-01    2016-11-28
    321                                     Asian    1997-01-01    2016-11-28
    322                                     Asian    1997-01-01    2016-11-28
    323                                     Asian    1997-01-01    2016-11-28
    324                                     Asian    1997-01-01    2016-11-28
    325                                     Asian    1997-01-01    2016-11-28
    326                                     Asian    1997-01-01    2016-11-28
    327                                     White    1997-01-01    2016-11-28
    328                                     White    1997-01-01    2016-11-28
    329                                     White    1997-01-01    2016-11-28
    330                                     White    1997-01-01    2016-11-28
    331                                     White    1997-01-01    2016-11-28
    332                        More Than One Race    1998-01-01    2017-01-03
    333                        More Than One Race    1998-01-01    2017-01-03
    334                        More Than One Race    1998-01-01    2017-01-03
    335                        More Than One Race    1998-01-01    2017-01-03
    336                        More Than One Race    1998-01-01    2017-01-03
    337                        More Than One Race    1998-01-01    2017-01-03
    338                        More Than One Race    1998-01-01    2017-01-03
    339                   Unknown or Not Reported    1998-01-01    2017-01-03
    340                   Unknown or Not Reported    1998-01-01    2017-01-03
    341                   Unknown or Not Reported    1998-01-01    2017-01-03
    342                   Unknown or Not Reported    1998-01-01    2017-01-03
    343                   Unknown or Not Reported    1998-01-01    2017-01-03
    344                   Unknown or Not Reported    1998-01-01    2017-01-03
    345                   Unknown or Not Reported    1998-01-01    2017-01-03
    346                                     Asian    1997-01-01    2017-01-17
    347                                     Asian    1997-01-01    2017-01-17
    348                                     Asian    1997-01-01    2017-01-17
    349                                     Asian    1997-01-01    2017-01-17
    350                                     Asian    1997-01-01    2017-01-17
    351                                     Asian    1997-01-01    2017-01-17
    352                                     Asian    1997-01-01    2017-01-17
    353                                     Asian    1997-01-01    2017-01-17
    354                                     Asian    1997-01-01    2017-01-17
    355                                     Asian    1997-01-01    2017-01-17
    356                                     Asian    1997-01-01    2017-01-17
    357                                     Asian    1997-01-01    2017-01-17
    358                                     Asian    1997-01-01    2017-01-17
    359                                     Asian    1997-01-01    2017-01-17
    360                                     Asian    1997-01-01    2017-01-30
    361                                     Asian    1997-01-01    2017-01-30
    362                                     Asian    1997-01-01    2017-01-30
    363                                     Asian    1997-01-01    2017-01-30
    364                                     Asian    1997-01-01    2017-01-30
    365                                     Asian    1997-01-01    2017-01-30
    366                                     Asian    1997-01-01    2017-01-30
    367                                     Asian    1996-01-01    2017-01-30
    368                                     Asian    1996-01-01    2017-01-30
    369                                     Asian    1996-01-01    2017-01-30
    370                                     Asian    1996-01-01    2017-01-30
    371                                     Asian    1996-01-01    2017-01-30
    372                                     Asian    1996-01-01    2017-01-30
    373                                     Asian    1996-01-01    2017-01-30
    374                   Unknown or Not Reported    1997-01-01    2017-01-30
    375                   Unknown or Not Reported    1997-01-01    2017-01-30
    376                   Unknown or Not Reported    1997-01-01    2017-01-30
    377                   Unknown or Not Reported    1997-01-01    2017-01-30
    378                   Unknown or Not Reported    1997-01-01    2017-01-30
    379                   Unknown or Not Reported    1997-01-01    2017-01-30
    380                   Unknown or Not Reported    1997-01-01    2017-01-30
    381                        More Than One Race    1997-01-01    2017-01-30
    382                        More Than One Race    1997-01-01    2017-01-30
    383                        More Than One Race    1997-01-01    2017-01-30
    384                        More Than One Race    1997-01-01    2017-01-30
    385                        More Than One Race    1997-01-01    2017-01-30
    386                        More Than One Race    1997-01-01    2017-01-30
    387                        More Than One Race    1997-01-01    2017-01-30
    388                                     White    1997-01-01    2017-01-30
    389                                     White    1997-01-01    2017-01-30
    390                                     White    1997-01-01    2017-01-30
    391                                     White    1997-01-01    2017-01-30
    392                                     White    1997-01-01    2017-01-30
    393                                     White    1997-01-01    2017-01-30
    394                                     White    1997-01-01    2017-01-30
             dataset        age  boost_age specimen_id actual_day_relative_to_boost
    1   2020_dataset 14586 days 11212 days           1                           -3
    2   2020_dataset 14586 days 11212 days           2                            1
    3   2020_dataset 14586 days 11212 days           3                            3
    4   2020_dataset 14586 days 11212 days           4                            7
    5   2020_dataset 14586 days 11212 days           5                           11
    6   2020_dataset 14586 days 11212 days           6                           32
    7   2020_dataset 14586 days 11212 days           7                          100
    8   2020_dataset 15682 days 12336 days          19                           -3
    9   2020_dataset 15682 days 12336 days          20                            1
    10  2020_dataset 15682 days 12336 days          21                            3
    11  2020_dataset 15682 days 12336 days          22                            7
    12  2020_dataset 15682 days 12336 days          23                           14
    13  2020_dataset 15682 days 12336 days          24                           30
    14  2020_dataset 15682 days 12336 days          25                           92
    15  2020_dataset 13856 days 10468 days          27                           -7
    16  2020_dataset 13856 days 10468 days          28                            1
    17  2020_dataset 13856 days 10468 days          29                            3
    18  2020_dataset 13856 days 10468 days          30                            8
    19  2020_dataset 13856 days 10468 days          31                           14
    20  2020_dataset 13856 days 10468 days          32                           32
    21  2020_dataset 13856 days 10468 days          33                          108
    22  2020_dataset 12760 days  9372 days          37                           -5
    23  2020_dataset 12760 days  9372 days          38                            1
    24  2020_dataset 12760 days  9372 days          39                            3
    25  2020_dataset 12760 days  9372 days          40                            8
    26  2020_dataset 12760 days  9372 days          41                           14
    27  2020_dataset 12760 days  9372 days          42                           30
    28  2020_dataset 12760 days  9372 days          43                           92
    29  2020_dataset 13856 days 10510 days          45                           -6
    30  2020_dataset 13856 days 10510 days          46                            1
    31  2020_dataset 13856 days 10510 days          47                            3
    32  2020_dataset 13856 days 10510 days          48                            7
    33  2020_dataset 13856 days 10510 days          49                           14
    34  2020_dataset 13856 days 10510 days          50                           31
    35  2020_dataset 13856 days 10510 days          51                           92
    36  2020_dataset 16412 days 13094 days          55                           -6
    37  2020_dataset 16412 days 13094 days          56                            3
    38  2020_dataset 16412 days 13094 days          57                            7
    39  2020_dataset 16412 days 13094 days          58                           14
    40  2020_dataset 16412 days 13094 days          59                           30
    41  2020_dataset 16412 days 13094 days          60                          100
    42  2020_dataset 16412 days 13094 days          61                          386
    43  2020_dataset 10934 days  7511 days          70                           -4
    44  2020_dataset 10934 days  7511 days          71                            1
    45  2020_dataset 10934 days  7511 days          72                            3
    46  2020_dataset 10934 days  7511 days          73                            7
    47  2020_dataset 10934 days  7511 days          74                           14
    48  2020_dataset 10934 days  7511 days          75                           30
    49  2020_dataset 10934 days  7511 days          76                          126
    50  2020_dataset 16047 days 12624 days          77                           -4
    51  2020_dataset 16047 days 12624 days          78                            1
    52  2020_dataset 16047 days 12624 days          79                            3
    53  2020_dataset 16047 days 12624 days          80                            7
    54  2020_dataset 16047 days 12624 days          81                           14
    55  2020_dataset 16047 days 12624 days          82                           31
    56  2020_dataset 14586 days 11198 days          87                          -12
    57  2020_dataset 14586 days 11198 days          88                            1
    58  2020_dataset 14586 days 11198 days          89                            3
    59  2020_dataset 14586 days 11198 days          90                            8
    60  2020_dataset 14586 days 11198 days          91                           14
    61  2020_dataset 14586 days 11198 days          92                           30
    62  2020_dataset 14586 days 11198 days          93                           91
    63  2020_dataset 16047 days 12624 days          96                           -4
    64  2020_dataset 16047 days 12624 days          97                            1
    65  2020_dataset 16047 days 12624 days          98                            3
    66  2020_dataset 16047 days 12624 days          99                            7
    67  2020_dataset 16047 days 12624 days         100                           14
    68  2020_dataset 16047 days 12624 days         101                           31
    69  2020_dataset 10568 days  7145 days         102                            0
    70  2020_dataset 10568 days  7145 days         103                            1
    71  2020_dataset 10568 days  7145 days         104                            3
    72  2020_dataset 10568 days  7145 days         105                            7
    73  2020_dataset 10568 days  7145 days         106                           14
    74  2020_dataset 10568 days  7145 days         107                           39
    75  2020_dataset 10568 days  7145 days         108                          126
    76  2020_dataset 12029 days  8627 days         109                           -5
    77  2020_dataset 12029 days  8627 days         110                            1
    78  2020_dataset 12029 days  8627 days         111                            3
    79  2020_dataset 12029 days  8627 days         112                            7
    80  2020_dataset 12029 days  8627 days         113                           14
    81  2020_dataset 13490 days 10088 days         114                            0
    82  2020_dataset 13490 days 10088 days         115                            1
    83  2020_dataset 13490 days 10088 days         116                            3
    84  2020_dataset 13490 days 10088 days         117                            7
    85  2020_dataset 13490 days 10088 days         118                           14
    86  2020_dataset 13490 days 10088 days         119                           31
    87  2020_dataset 13490 days 10088 days         120                           92
    88  2020_dataset 14221 days 10798 days         121                            0
    89  2020_dataset 14221 days 10798 days         122                            1
    90  2020_dataset 14221 days 10798 days         123                            3
    91  2020_dataset 14221 days 10798 days         124                            7
    92  2020_dataset 14221 days 10798 days         125                           14
    93  2020_dataset 14221 days 10798 days         126                           31
    94  2020_dataset 14221 days 10798 days         127                           92
    95  2020_dataset 16778 days 13404 days         131                          -40
    96  2020_dataset 16778 days 13404 days         132                            1
    97  2020_dataset 16778 days 13404 days         133                            3
    98  2020_dataset 16778 days 13404 days         134                            7
    99  2020_dataset 16778 days 13404 days         135                           14
    100 2020_dataset 16778 days 13404 days         136                           38
    101 2020_dataset 16778 days 13404 days         137                           91
    102 2020_dataset 10568 days  7180 days         138                          -12
    103 2020_dataset 10568 days  7180 days         139                            1
    104 2020_dataset 10568 days  7180 days         140                            3
    105 2020_dataset 10568 days  7180 days         141                            8
    106 2020_dataset 10568 days  7180 days         142                           14
    107 2020_dataset 10568 days  7180 days         143                           30
    108 2020_dataset 10568 days  7180 days         144                           93
    109 2020_dataset 11664 days  8304 days         146                          -34
    110 2020_dataset 11664 days  8304 days         147                            1
    111 2020_dataset 11664 days  8304 days         148                            3
    112 2020_dataset 11664 days  8304 days         149                            7
    113 2020_dataset 11664 days  8304 days         150                           14
    114 2020_dataset 11664 days  8304 days         151                           30
    115 2020_dataset 11664 days  8304 days         152                          120
    116 2020_dataset 16412 days 13024 days         153                            0
    117 2020_dataset 16412 days 13024 days         154                            1
    118 2020_dataset 16412 days 13024 days         155                            3
    119 2020_dataset 16412 days 13024 days         156                            8
    120 2020_dataset 16412 days 13024 days         157                           14
    121 2020_dataset 16412 days 13024 days         158                           37
    122 2020_dataset 16412 days 13024 days         159                           93
    123 2020_dataset 15682 days 12294 days         160                            0
    124 2020_dataset 15682 days 12294 days         161                            1
    125 2020_dataset 15682 days 12294 days         162                            3
    126 2020_dataset 15682 days 12294 days         163                            8
    127 2020_dataset 15682 days 12294 days         164                           14
    128 2020_dataset 15682 days 12294 days         165                           30
    129 2020_dataset 15682 days 12294 days         166                           93
    130 2020_dataset 14951 days 11563 days         167                            0
    131 2020_dataset 14951 days 11563 days         168                            1
    132 2020_dataset 14951 days 11563 days         169                            3
    133 2020_dataset 14951 days 11563 days         170                            8
    134 2020_dataset 14951 days 11563 days         171                           14
    135 2020_dataset 14951 days 11563 days         172                           30
    136 2020_dataset 14951 days 11563 days         173                           93
    137 2020_dataset 12760 days  9400 days         174                          -26
    138 2020_dataset 12760 days  9400 days         175                            1
    139 2020_dataset 12760 days  9400 days         176                            3
    140 2020_dataset 12760 days  9400 days         177                            7
    141 2020_dataset 12760 days  9400 days         178                           14
    142 2020_dataset 12760 days  9400 days         179                           37
    143 2020_dataset 12760 days  9400 days         180                          115
    144 2020_dataset 12395 days  9022 days         181                          -13
    145 2020_dataset 12395 days  9022 days         182                            0
    146 2020_dataset 12395 days  9022 days         183                            2
    147 2020_dataset 12395 days  9022 days         184                            6
    148 2020_dataset 12395 days  9022 days         185                           13
    149 2020_dataset 12395 days  9022 days         186                           29
    150 2020_dataset 12395 days  9022 days         187                           94
    151 2020_dataset 13856 days 10483 days         191                           -6
    152 2020_dataset 13856 days 10483 days         192                            0
    153 2020_dataset 13856 days 10483 days         193                            2
    154 2020_dataset 13856 days 10483 days         194                            6
    155 2020_dataset 13856 days 10483 days         195                           13
    156 2020_dataset 13856 days 10483 days         196                           55
    157 2020_dataset 13856 days 10483 days         197                          112
    158 2020_dataset 15682 days 12322 days         201                           -7
    159 2020_dataset 15682 days 12322 days         202                            1
    160 2020_dataset 15682 days 12322 days         203                            3
    161 2020_dataset 15682 days 12322 days         204                            7
    162 2020_dataset 15682 days 12322 days         205                           14
    163 2020_dataset 15682 days 12322 days         206                           30
    164 2020_dataset 15682 days 12322 days         207                          107
    165 2020_dataset 10568 days  7208 days         208                           -5
    166 2020_dataset 10568 days  7208 days         209                            1
    167 2020_dataset 10568 days  7208 days         210                            3
    168 2020_dataset 10568 days  7208 days         211                            7
    169 2020_dataset 10568 days  7208 days         212                           14
    170 2020_dataset 10568 days  7208 days         213                           30
    171 2020_dataset 10568 days  7208 days         214                          108
    172 2020_dataset 16047 days 12687 days         216                           -4
    173 2020_dataset 16047 days 12687 days         217                            1
    174 2020_dataset 16047 days 12687 days         218                            3
    175 2020_dataset 16047 days 12687 days         219                            7
    176 2020_dataset 16047 days 12687 days         220                           14
    177 2020_dataset 16047 days 12687 days         221                           36
    178 2020_dataset 16047 days 12687 days         222                          163
    179 2020_dataset 10568 days  7208 days         223                           -4
    180 2020_dataset 10568 days  7208 days         224                            1
    181 2020_dataset 10568 days  7208 days         225                            3
    182 2020_dataset 10568 days  7208 days         226                            7
    183 2020_dataset 10568 days  7208 days         227                           18
    184 2020_dataset 10568 days  7208 days         228                           37
    185 2020_dataset 10568 days  7208 days         229                           93
    186 2020_dataset 13856 days 10496 days         232                           -4
    187 2020_dataset 13856 days 10496 days         233                            1
    188 2020_dataset 13856 days 10496 days         234                            3
    189 2020_dataset 13856 days 10496 days         235                            7
    190 2020_dataset 13856 days 10496 days         236                           14
    191 2020_dataset 13856 days 10496 days         237                           32
    192 2020_dataset 13856 days 10496 days         238                          129
    193 2020_dataset 13490 days 10130 days         241                            0
    194 2020_dataset 13490 days 10130 days         242                            1
    195 2020_dataset 13490 days 10130 days         243                            3
    196 2020_dataset 13490 days 10130 days         244                            7
    197 2020_dataset 13490 days 10130 days         245                           14
    198 2020_dataset 13490 days 10130 days         246                           31
    199 2020_dataset 13490 days 10130 days         247                          428
    200 2020_dataset 10568 days  7236 days         248                          -19
    201 2020_dataset 10568 days  7236 days         249                            1
    202 2020_dataset 10568 days  7236 days         250                            3
    203 2020_dataset 10568 days  7236 days         251                            7
    204 2020_dataset 10568 days  7236 days         252                           16
    205 2020_dataset 10568 days  7236 days         253                           30
    206 2020_dataset 10568 days  7236 days         254                          112
    207 2020_dataset 13125 days  9779 days         255                           -6
    208 2020_dataset 13125 days  9779 days         256                            1
    209 2020_dataset 13125 days  9779 days         257                            3
    210 2020_dataset 13125 days  9779 days         258                            7
    211 2020_dataset 13125 days  9779 days         259                           15
    212 2020_dataset 13125 days  9779 days         260                           30
    213 2020_dataset 13125 days  9779 days         261                           92
    214 2020_dataset 15682 days 12350 days         266                          -18
    215 2020_dataset 15682 days 12350 days         267                            1
    216 2020_dataset 15682 days 12350 days         268                            3
    217 2020_dataset 15682 days 12350 days         269                           14
    218 2020_dataset 15682 days 12350 days         270                           30
    219 2020_dataset 15682 days 12350 days         271                           92
    220 2020_dataset 15682 days 12350 days         272                          402
    221 2020_dataset 12760 days  9414 days         274                           -4
    222 2020_dataset 12760 days  9414 days         275                            1
    223 2020_dataset 12760 days  9414 days         276                            3
    224 2020_dataset 12760 days  9414 days         277                            7
    225 2020_dataset 12760 days  9414 days         278                           14
    226 2020_dataset 12760 days  9414 days         279                           37
    227 2020_dataset 12760 days  9414 days         280                           94
    228 2020_dataset 10568 days  7236 days         281                           -6
    229 2020_dataset 10568 days  7236 days         282                            1
    230 2020_dataset 10568 days  7236 days         283                            3
    231 2020_dataset 10568 days  7236 days         284                            7
    232 2020_dataset 10568 days  7236 days         285                           14
    233 2020_dataset 10568 days  7236 days         286                           30
    234 2020_dataset 10203 days  6885 days         288                          -20
    235 2020_dataset 10203 days  6885 days         289                            1
    236 2020_dataset 10203 days  6885 days         290                            7
    237 2020_dataset 10203 days  6885 days         291                           30
    238 2020_dataset 10203 days  6885 days         292                           99
    239 2020_dataset 10568 days  7236 days         293                           -5
    240 2020_dataset 10568 days  7236 days         294                            1
    241 2020_dataset 10568 days  7236 days         295                            3
    242 2020_dataset 10568 days  7236 days         296                            7
    243 2020_dataset 10568 days  7236 days         297                           14
    244 2020_dataset 10568 days  7236 days         298                           36
    245 2020_dataset 10568 days  7236 days         299                          106
    246 2020_dataset 14951 days 11619 days         300                           -4
    247 2020_dataset 14951 days 11619 days         301                            1
    248 2020_dataset 14951 days 11619 days         302                            3
    249 2020_dataset 14951 days 11619 days         303                            7
    250 2020_dataset 14951 days 11619 days         304                           14
    251 2020_dataset 14951 days 11619 days         305                           30
    252 2020_dataset 14951 days 11619 days         306                           95
    253 2020_dataset 11664 days  8332 days         310                           -3
    254 2020_dataset 11664 days  8332 days         311                            1
    255 2020_dataset 11664 days  8332 days         312                            3
    256 2020_dataset 11664 days  8332 days         313                            7
    257 2020_dataset 11664 days  8332 days         314                           14
    258 2020_dataset 11664 days  8332 days         315                           35
    259 2020_dataset 11664 days  8332 days         316                           94
    260 2020_dataset 14951 days 11633 days         317                          -10
    261 2020_dataset 14951 days 11633 days         318                            1
    262 2020_dataset 14951 days 11633 days         319                            3
    263 2020_dataset 14951 days 11633 days         320                            7
    264 2020_dataset 14951 days 11633 days         321                           14
    265 2020_dataset 14951 days 11633 days         322                           30
    266 2020_dataset 14951 days 11633 days         323                           92
    267 2020_dataset 10568 days  7250 days         324                           -6
    268 2020_dataset 10568 days  7250 days         325                            1
    269 2020_dataset 10568 days  7250 days         326                            3
    270 2020_dataset 10568 days  7250 days         327                            7
    271 2020_dataset 10568 days  7250 days         328                           14
    272 2020_dataset 10568 days  7250 days         329                           30
    273 2020_dataset 10568 days  7250 days         330                          107
    274 2020_dataset 10203 days  6885 days         332                           -6
    275 2020_dataset 10203 days  6885 days         333                            1
    276 2020_dataset 10203 days  6885 days         334                            3
    277 2020_dataset 10203 days  6885 days         335                            7
    278 2020_dataset 10203 days  6885 days         336                           14
    279 2020_dataset 10203 days  6885 days         337                           32
    280 2020_dataset 10203 days  6885 days         338                          101
    281 2020_dataset 10203 days  6885 days         342                           -5
    282 2020_dataset 10203 days  6885 days         343                            1
    283 2020_dataset 10203 days  6885 days         344                            3
    284 2020_dataset 10203 days  6885 days         345                            7
    285 2020_dataset 10203 days  6885 days         346                           14
    286 2020_dataset 10203 days  6885 days         347                           30
    287 2020_dataset 10203 days  6885 days         348                          100
    288 2020_dataset 10568 days  7271 days         349                          -26
    289 2020_dataset 10568 days  7271 days         350                            1
    290 2020_dataset 10568 days  7271 days         351                            3
    291 2020_dataset 10568 days  7271 days         352                            7
    292 2020_dataset 10568 days  7271 days         353                           14
    293 2020_dataset 10568 days  7271 days         354                           99
    294 2020_dataset 10203 days  6885 days         355                           -4
    295 2020_dataset 10203 days  6885 days         356                            1
    296 2020_dataset 10203 days  6885 days         357                            3
    297 2020_dataset 10203 days  6885 days         358                            7
    298 2020_dataset 10203 days  6885 days         359                           14
    299 2020_dataset 10934 days  7637 days         360                          -13
    300 2020_dataset 10934 days  7637 days         361                            1
    301 2020_dataset 10934 days  7637 days         362                            3
    302 2020_dataset 10934 days  7637 days         363                            7
    303 2020_dataset 10934 days  7637 days         364                           14
    304 2020_dataset 10934 days  7637 days         365                           29
    305 2020_dataset 10934 days  7637 days         366                           94
    306 2020_dataset 10203 days  6956 days         369                          -63
    307 2020_dataset 10203 days  6956 days         370                            1
    308 2020_dataset 10203 days  6956 days         371                            7
    309 2020_dataset 10203 days  6956 days         372                            7
    310 2020_dataset 10203 days  6956 days         373                           14
    311 2020_dataset 10203 days  6956 days         374                           36
    312 2020_dataset 10203 days  6956 days         375                          105
    313 2020_dataset 10568 days  7321 days         376                          -56
    314 2020_dataset 10568 days  7321 days         377                            1
    315 2020_dataset 10568 days  7321 days         378                            7
    316 2020_dataset 10568 days  7321 days         379                            7
    317 2020_dataset 10568 days  7321 days         380                           14
    318 2020_dataset 10568 days  7321 days         381                           31
    319 2020_dataset 10568 days  7321 days         382                           90
    320 2020_dataset 10568 days  7271 days         385                           -6
    321 2020_dataset 10568 days  7271 days         386                            1
    322 2020_dataset 10568 days  7271 days         387                            3
    323 2020_dataset 10568 days  7271 days         388                            7
    324 2020_dataset 10568 days  7271 days         389                           15
    325 2020_dataset 10568 days  7271 days         390                           36
    326 2020_dataset 10568 days  7271 days         391                          116
    327 2020_dataset 10568 days  7271 days         392                           -6
    328 2020_dataset 10568 days  7271 days         393                            1
    329 2020_dataset 10568 days  7271 days         394                            3
    330 2020_dataset 10568 days  7271 days         395                            8
    331 2020_dataset 10568 days  7271 days         396                           14
    332 2020_dataset 10203 days  6942 days         397                          -34
    333 2020_dataset 10203 days  6942 days         398                            1
    334 2020_dataset 10203 days  6942 days         399                            3
    335 2020_dataset 10203 days  6942 days         400                            8
    336 2020_dataset 10203 days  6942 days         401                           14
    337 2020_dataset 10203 days  6942 days         402                           31
    338 2020_dataset 10203 days  6942 days         403                           90
    339 2020_dataset 10203 days  6942 days         405                          -28
    340 2020_dataset 10203 days  6942 days         406                            1
    341 2020_dataset 10203 days  6942 days         407                            3
    342 2020_dataset 10203 days  6942 days         408                            8
    343 2020_dataset 10203 days  6942 days         409                           14
    344 2020_dataset 10203 days  6942 days         410                           31
    345 2020_dataset 10203 days  6942 days         411                          100
    346 2020_dataset 10568 days  7321 days         412                          -36
    347 2020_dataset 10568 days  7321 days         413                            1
    348 2020_dataset 10568 days  7321 days         414                            7
    349 2020_dataset 10568 days  7321 days         415                            7
    350 2020_dataset 10568 days  7321 days         416                           14
    351 2020_dataset 10568 days  7321 days         417                           42
    352 2020_dataset 10568 days  7321 days         418                           94
    353 2020_dataset 10568 days  7321 days         419                           -8
    354 2020_dataset 10568 days  7321 days         420                            1
    355 2020_dataset 10568 days  7321 days         421                            7
    356 2020_dataset 10568 days  7321 days         422                            7
    357 2020_dataset 10568 days  7321 days         423                           14
    358 2020_dataset 10568 days  7321 days         424                           31
    359 2020_dataset 10568 days  7321 days         425                          107
    360 2020_dataset 10568 days  7334 days         427                          -18
    361 2020_dataset 10568 days  7334 days         428                            1
    362 2020_dataset 10568 days  7334 days         429                            3
    363 2020_dataset 10568 days  7334 days         430                            8
    364 2020_dataset 10568 days  7334 days         431                           14
    365 2020_dataset 10568 days  7334 days         432                           29
    366 2020_dataset 10568 days  7334 days         433                          116
    367 2020_dataset 10934 days  7700 days         434                           -6
    368 2020_dataset 10934 days  7700 days         435                            1
    369 2020_dataset 10934 days  7700 days         436                            4
    370 2020_dataset 10934 days  7700 days         437                            7
    371 2020_dataset 10934 days  7700 days         438                           14
    372 2020_dataset 10934 days  7700 days         439                           30
    373 2020_dataset 10934 days  7700 days         440                           95
    374 2020_dataset 10568 days  7334 days         441                           -5
    375 2020_dataset 10568 days  7334 days         442                            1
    376 2020_dataset 10568 days  7334 days         443                            3
    377 2020_dataset 10568 days  7334 days         444                            7
    378 2020_dataset 10568 days  7334 days         445                           14
    379 2020_dataset 10568 days  7334 days         446                           29
    380 2020_dataset 10568 days  7334 days         447                           92
    381 2020_dataset 10568 days  7334 days         450                           -5
    382 2020_dataset 10568 days  7334 days         451                            1
    383 2020_dataset 10568 days  7334 days         452                            3
    384 2020_dataset 10568 days  7334 days         453                            7
    385 2020_dataset 10568 days  7334 days         454                           14
    386 2020_dataset 10568 days  7334 days         455                           29
    387 2020_dataset 10568 days  7334 days         456                           92
    388 2020_dataset 10568 days  7334 days         458                           -4
    389 2020_dataset 10568 days  7334 days         459                            1
    390 2020_dataset 10568 days  7334 days         460                            3
    391 2020_dataset 10568 days  7334 days         461                            7
    392 2020_dataset 10568 days  7334 days         462                           14
    393 2020_dataset 10568 days  7334 days         463                           29
    394 2020_dataset 10568 days  7334 days         464                           98
        planned_day_relative_to_boost specimen_type visit isotype
    1                               0         Blood     1     IgG
    2                               1         Blood     2     IgG
    3                               3         Blood     3     IgG
    4                               7         Blood     4     IgG
    5                              14         Blood     5     IgG
    6                              30         Blood     6     IgG
    7                             120         Blood     7     IgG
    8                               0         Blood     1     IgG
    9                               1         Blood     2     IgG
    10                              3         Blood     3     IgG
    11                              7         Blood     4     IgG
    12                             14         Blood     5     IgG
    13                             30         Blood     6     IgG
    14                            120         Blood     7     IgG
    15                              0         Blood     1     IgG
    16                              1         Blood     2     IgG
    17                              3         Blood     3     IgG
    18                              7         Blood     4     IgG
    19                             14         Blood     5     IgG
    20                             30         Blood     6     IgG
    21                            120         Blood     7     IgG
    22                              0         Blood     1     IgG
    23                              1         Blood     2     IgG
    24                              3         Blood     3     IgG
    25                              7         Blood     4     IgG
    26                             14         Blood     5     IgG
    27                             30         Blood     6     IgG
    28                            120         Blood     7     IgG
    29                              0         Blood     1     IgG
    30                              1         Blood     2     IgG
    31                              3         Blood     3     IgG
    32                              7         Blood     4     IgG
    33                             14         Blood     5     IgG
    34                             30         Blood     6     IgG
    35                            120         Blood     7     IgG
    36                              0         Blood     1     IgG
    37                              1         Blood     2     IgG
    38                              3         Blood     3     IgG
    39                             14         Blood     4     IgG
    40                             30         Blood     5     IgG
    41                            120         Blood     6     IgG
    42                            386         Blood     7     IgG
    43                              0         Blood     1     IgG
    44                              1         Blood     2     IgG
    45                              3         Blood     3     IgG
    46                              7         Blood     4     IgG
    47                             14         Blood     5     IgG
    48                             30         Blood     6     IgG
    49                            120         Blood     7     IgG
    50                              0         Blood     1     IgG
    51                              1         Blood     2     IgG
    52                              3         Blood     3     IgG
    53                              7         Blood     4     IgG
    54                             14         Blood     5     IgG
    55                             30         Blood     6     IgG
    56                              0         Blood     1     IgG
    57                              1         Blood     2     IgG
    58                              3         Blood     3     IgG
    59                              7         Blood     4     IgG
    60                             14         Blood     5     IgG
    61                             30         Blood     6     IgG
    62                            120         Blood     7     IgG
    63                              0         Blood     1     IgG
    64                              1         Blood     2     IgG
    65                              3         Blood     3     IgG
    66                              7         Blood     4     IgG
    67                             14         Blood     5     IgG
    68                             30         Blood     6     IgG
    69                              0         Blood     1     IgG
    70                              1         Blood     2     IgG
    71                              3         Blood     3     IgG
    72                              7         Blood     4     IgG
    73                             14         Blood     5     IgG
    74                             30         Blood     6     IgG
    75                            120         Blood     7     IgG
    76                              0         Blood     1     IgG
    77                              1         Blood     2     IgG
    78                              3         Blood     3     IgG
    79                              7         Blood     4     IgG
    80                             14         Blood     5     IgG
    81                              0         Blood     1     IgG
    82                              1         Blood     2     IgG
    83                              3         Blood     3     IgG
    84                              7         Blood     4     IgG
    85                             14         Blood     5     IgG
    86                             30         Blood     6     IgG
    87                            120         Blood     7     IgG
    88                              0         Blood     1     IgG
    89                              1         Blood     2     IgG
    90                              3         Blood     3     IgG
    91                              7         Blood     4     IgG
    92                             14         Blood     5     IgG
    93                             30         Blood     6     IgG
    94                            120         Blood     7     IgG
    95                              0         Blood     1     IgG
    96                              1         Blood     2     IgG
    97                              3         Blood     3     IgG
    98                              7         Blood     4     IgG
    99                             14         Blood     5     IgG
    100                            30         Blood     6     IgG
    101                           120         Blood     7     IgG
    102                             0         Blood     1     IgG
    103                             1         Blood     2     IgG
    104                             3         Blood     3     IgG
    105                             7         Blood     4     IgG
    106                            14         Blood     5     IgG
    107                            30         Blood     6     IgG
    108                           120         Blood     7     IgG
    109                             0         Blood     1     IgG
    110                             1         Blood     2     IgG
    111                             3         Blood     3     IgG
    112                             7         Blood     4     IgG
    113                            14         Blood     5     IgG
    114                            30         Blood     6     IgG
    115                           120         Blood     7     IgG
    116                             0         Blood     1     IgG
    117                             1         Blood     2     IgG
    118                             3         Blood     3     IgG
    119                             7         Blood     4     IgG
    120                            14         Blood     5     IgG
    121                            30         Blood     6     IgG
    122                           120         Blood     7     IgG
    123                             0         Blood     1     IgG
    124                             1         Blood     2     IgG
    125                             3         Blood     3     IgG
    126                             7         Blood     4     IgG
    127                            14         Blood     5     IgG
    128                            30         Blood     6     IgG
    129                           120         Blood     7     IgG
    130                             0         Blood     1     IgG
    131                             1         Blood     2     IgG
    132                             3         Blood     3     IgG
    133                             7         Blood     4     IgG
    134                            14         Blood     5     IgG
    135                            30         Blood     6     IgG
    136                           120         Blood     7     IgG
    137                             0         Blood     1     IgG
    138                             1         Blood     2     IgG
    139                             3         Blood     3     IgG
    140                             7         Blood     4     IgG
    141                            14         Blood     5     IgG
    142                            30         Blood     6     IgG
    143                           120         Blood     7     IgG
    144                             0         Blood     1     IgG
    145                             1         Blood     2     IgG
    146                             3         Blood     3     IgG
    147                             7         Blood     4     IgG
    148                            14         Blood     5     IgG
    149                            30         Blood     6     IgG
    150                           120         Blood     7     IgG
    151                             0         Blood     1     IgG
    152                             1         Blood     2     IgG
    153                             3         Blood     3     IgG
    154                             7         Blood     4     IgG
    155                            14         Blood     5     IgG
    156                            30         Blood     6     IgG
    157                           120         Blood     7     IgG
    158                             0         Blood     1     IgG
    159                             1         Blood     2     IgG
    160                             3         Blood     3     IgG
    161                             7         Blood     4     IgG
    162                            14         Blood     5     IgG
    163                            30         Blood     6     IgG
    164                           120         Blood     7     IgG
    165                             0         Blood     1     IgG
    166                             1         Blood     2     IgG
    167                             3         Blood     3     IgG
    168                             7         Blood     4     IgG
    169                            14         Blood     5     IgG
    170                            30         Blood     6     IgG
    171                           120         Blood     7     IgG
    172                             0         Blood     1     IgG
    173                             1         Blood     2     IgG
    174                             3         Blood     3     IgG
    175                             7         Blood     4     IgG
    176                            14         Blood     5     IgG
    177                            30         Blood     6     IgG
    178                           120         Blood     7     IgG
    179                             0         Blood     1     IgG
    180                             1         Blood     2     IgG
    181                             3         Blood     3     IgG
    182                             7         Blood     4     IgG
    183                            14         Blood     5     IgG
    184                            30         Blood     6     IgG
    185                           120         Blood     7     IgG
    186                             0         Blood     1     IgG
    187                             1         Blood     2     IgG
    188                             3         Blood     3     IgG
    189                             7         Blood     4     IgG
    190                            14         Blood     5     IgG
    191                            30         Blood     6     IgG
    192                           120         Blood     7     IgG
    193                             0         Blood     1     IgG
    194                             1         Blood     2     IgG
    195                             3         Blood     3     IgG
    196                             7         Blood     4     IgG
    197                            14         Blood     5     IgG
    198                            30         Blood     6     IgG
    199                           428         Blood     8     IgG
    200                             0         Blood     1     IgG
    201                             1         Blood     2     IgG
    202                             3         Blood     3     IgG
    203                             7         Blood     4     IgG
    204                            14         Blood     5     IgG
    205                            30         Blood     6     IgG
    206                           120         Blood     7     IgG
    207                             0         Blood     1     IgG
    208                             1         Blood     2     IgG
    209                             3         Blood     3     IgG
    210                             7         Blood     4     IgG
    211                            14         Blood     5     IgG
    212                            30         Blood     6     IgG
    213                           120         Blood     7     IgG
    214                             0         Blood     1     IgG
    215                             1         Blood     2     IgG
    216                             3         Blood     3     IgG
    217                            14         Blood     4     IgG
    218                            30         Blood     5     IgG
    219                           120         Blood     6     IgG
    220                           402         Blood     7     IgG
    221                             0         Blood     1     IgG
    222                             1         Blood     2     IgG
    223                             3         Blood     3     IgG
    224                             7         Blood     4     IgG
    225                            14         Blood     5     IgG
    226                            30         Blood     6     IgG
    227                           120         Blood     7     IgG
    228                             0         Blood     1     IgG
    229                             1         Blood     2     IgG
    230                             3         Blood     3     IgG
    231                             7         Blood     4     IgG
    232                            14         Blood     5     IgG
    233                            30         Blood     6     IgG
    234                             0         Blood     1     IgG
    235                             1         Blood     2     IgG
    236                             7         Blood     3     IgG
    237                            30         Blood     4     IgG
    238                           120         Blood     5     IgG
    239                             0         Blood     1     IgG
    240                             1         Blood     2     IgG
    241                             3         Blood     3     IgG
    242                             7         Blood     4     IgG
    243                            14         Blood     5     IgG
    244                            30         Blood     6     IgG
    245                           120         Blood     7     IgG
    246                             0         Blood     1     IgG
    247                             1         Blood     2     IgG
    248                             3         Blood     3     IgG
    249                             7         Blood     4     IgG
    250                            14         Blood     5     IgG
    251                            30         Blood     6     IgG
    252                           120         Blood     7     IgG
    253                             0         Blood     1     IgG
    254                             1         Blood     2     IgG
    255                             3         Blood     3     IgG
    256                             7         Blood     4     IgG
    257                            14         Blood     5     IgG
    258                            30         Blood     6     IgG
    259                           120         Blood     7     IgG
    260                             0         Blood     1     IgG
    261                             1         Blood     2     IgG
    262                             3         Blood     3     IgG
    263                             7         Blood     4     IgG
    264                            14         Blood     5     IgG
    265                            30         Blood     6     IgG
    266                           120         Blood     7     IgG
    267                             0         Blood     1     IgG
    268                             1         Blood     2     IgG
    269                             3         Blood     3     IgG
    270                             7         Blood     4     IgG
    271                            14         Blood     5     IgG
    272                            30         Blood     6     IgG
    273                           120         Blood     7     IgG
    274                             0         Blood     1     IgG
    275                             1         Blood     2     IgG
    276                             3         Blood     3     IgG
    277                             7         Blood     4     IgG
    278                            14         Blood     5     IgG
    279                            30         Blood     6     IgG
    280                           120         Blood     7     IgG
    281                             0         Blood     1     IgG
    282                             1         Blood     2     IgG
    283                             3         Blood     3     IgG
    284                             7         Blood     4     IgG
    285                            14         Blood     5     IgG
    286                            30         Blood     6     IgG
    287                           120         Blood     7     IgG
    288                             0         Blood     1     IgG
    289                             1         Blood     2     IgG
    290                             3         Blood     3     IgG
    291                             7         Blood     4     IgG
    292                            14         Blood     5     IgG
    293                           120         Blood     6     IgG
    294                             0         Blood     1     IgG
    295                             1         Blood     2     IgG
    296                             3         Blood     3     IgG
    297                             7         Blood     4     IgG
    298                            14         Blood     5     IgG
    299                             0         Blood     1     IgG
    300                             1         Blood     2     IgG
    301                             3         Blood     3     IgG
    302                             7         Blood     4     IgG
    303                            14         Blood     5     IgG
    304                            30         Blood     6     IgG
    305                           120         Blood     7     IgG
    306                             0         Blood     1     IgG
    307                             1         Blood     2     IgG
    308                             3         Blood     3     IgG
    309                             7         Blood     4     IgG
    310                            14         Blood     5     IgG
    311                            30         Blood     6     IgG
    312                           120         Blood     7     IgG
    313                             0         Blood     1     IgG
    314                             1         Blood     2     IgG
    315                             3         Blood     3     IgG
    316                             7         Blood     4     IgG
    317                            14         Blood     5     IgG
    318                            30         Blood     6     IgG
    319                           120         Blood     7     IgG
    320                             0         Blood     1     IgG
    321                             1         Blood     2     IgG
    322                             3         Blood     3     IgG
    323                             7         Blood     4     IgG
    324                            14         Blood     5     IgG
    325                            30         Blood     6     IgG
    326                           120         Blood     7     IgG
    327                             0         Blood     1     IgG
    328                             1         Blood     2     IgG
    329                             3         Blood     3     IgG
    330                             7         Blood     4     IgG
    331                            14         Blood     5     IgG
    332                             0         Blood     1     IgG
    333                             1         Blood     2     IgG
    334                             3         Blood     3     IgG
    335                             7         Blood     4     IgG
    336                            14         Blood     5     IgG
    337                            30         Blood     6     IgG
    338                           120         Blood     7     IgG
    339                             0         Blood     1     IgG
    340                             1         Blood     2     IgG
    341                             3         Blood     3     IgG
    342                             7         Blood     4     IgG
    343                            14         Blood     5     IgG
    344                            30         Blood     6     IgG
    345                           120         Blood     7     IgG
    346                             0         Blood     1     IgG
    347                             1         Blood     2     IgG
    348                             3         Blood     3     IgG
    349                             7         Blood     4     IgG
    350                            14         Blood     5     IgG
    351                            30         Blood     6     IgG
    352                           120         Blood     7     IgG
    353                             0         Blood     1     IgG
    354                             1         Blood     2     IgG
    355                             3         Blood     3     IgG
    356                             7         Blood     4     IgG
    357                            14         Blood     5     IgG
    358                            30         Blood     6     IgG
    359                           120         Blood     7     IgG
    360                             0         Blood     1     IgG
    361                             1         Blood     2     IgG
    362                             3         Blood     3     IgG
    363                             7         Blood     4     IgG
    364                            14         Blood     5     IgG
    365                            30         Blood     6     IgG
    366                           120         Blood     7     IgG
    367                             0         Blood     1     IgG
    368                             1         Blood     2     IgG
    369                             3         Blood     3     IgG
    370                             7         Blood     4     IgG
    371                            14         Blood     5     IgG
    372                            30         Blood     6     IgG
    373                           120         Blood     7     IgG
    374                             0         Blood     1     IgG
    375                             1         Blood     2     IgG
    376                             3         Blood     3     IgG
    377                             7         Blood     4     IgG
    378                            14         Blood     5     IgG
    379                            30         Blood     6     IgG
    380                           120         Blood     7     IgG
    381                             0         Blood     1     IgG
    382                             1         Blood     2     IgG
    383                             3         Blood     3     IgG
    384                             7         Blood     4     IgG
    385                            14         Blood     5     IgG
    386                            30         Blood     6     IgG
    387                           120         Blood     7     IgG
    388                             0         Blood     1     IgG
    389                             1         Blood     2     IgG
    390                             3         Blood     3     IgG
    391                             7         Blood     4     IgG
    392                            14         Blood     5     IgG
    393                            30         Blood     6     IgG
    394                           120         Blood     7     IgG
        is_antigen_specific antigen          MFI MFI_normalised  unit
    1                  TRUE      PT 6.856614e+01     3.73699166 IU/ML
    2                  TRUE      PT 4.138442e+01     2.25553379 IU/ML
    3                  TRUE      PT 5.963762e+01     3.25036934 IU/ML
    4                  TRUE      PT 1.995177e+02    10.87411169 IU/ML
    5                  TRUE      PT 2.296037e+02    12.51386026 IU/ML
    6                  TRUE      PT 2.288146e+02    12.47085353 IU/ML
    7                  TRUE      PT 1.985562e+02    10.82171086 IU/ML
    8                  TRUE      PT 2.011607e+01     1.09636587 IU/ML
    9                  TRUE      PT 1.307630e+01     0.71268463 IU/ML
    10                 TRUE      PT 9.473398e+00     0.51631914 IU/ML
    11                 TRUE      PT 1.178683e+02     6.42405817 IU/ML
    12                 TRUE      PT 1.291980e+02     7.04154687 IU/ML
    13                 TRUE      PT 1.648531e+02     8.98482346 IU/ML
    14                 TRUE      PT 1.552430e+02     8.46105386 IU/ML
    15                 TRUE      PT 3.755222e+01     2.04667118 IU/ML
    16                 TRUE      PT 2.986034e+01     1.62744809 IU/ML
    17                 TRUE      PT 2.607730e+01     1.42126489 IU/ML
    18                 TRUE      PT 1.448853e+02     7.89654060 IU/ML
    19                 TRUE      PT 1.054266e+02     5.74595904 IU/ML
    20                 TRUE      PT 1.486644e+02     8.10250666 IU/ML
    21                 TRUE      PT 1.395544e+02     7.60599157 IU/ML
    22                 TRUE      PT 6.968565e+01     3.79800697 IU/ML
    23                 TRUE      PT 5.791494e+01     3.15648000 IU/ML
    24                 TRUE      PT 4.106182e+01     2.23795148 IU/ML
    25                 TRUE      PT 1.100317e+02     5.99694906 IU/ML
    26                 TRUE      PT 9.774326e+01     5.32720297 IU/ML
    27                 TRUE      PT 1.289116e+02     7.02594046 IU/ML
    28                 TRUE      PT 1.204283e+02     6.56358324 IU/ML
    29                 TRUE      PT 3.914130e+00     0.21332789 IU/ML
    30                 TRUE      PT 3.657555e+00     0.19934407 IU/ML
    31                 TRUE      PT 3.209682e+00     0.17493411 IU/ML
    32                 TRUE      PT 1.674964e+02     9.12888622 IU/ML
    33                 TRUE      PT 1.625000e+02     8.85657487 IU/ML
    34                 TRUE      PT 1.767945e+02     9.63565138 IU/ML
    35                 TRUE      PT 1.011998e+02     5.51558979 IU/ML
    36                 TRUE      PT 9.139656e+00     0.49812954 IU/ML
    37                 TRUE      PT 6.543368e+00     0.35662666 IU/ML
    38                 TRUE      PT 4.107241e+01     2.23852839 IU/ML
    39                 TRUE      PT 9.268763e+01     5.05166119 IU/ML
    40                 TRUE      PT 6.990745e+01     3.81009595 IU/ML
    41                 TRUE      PT 2.882964e+01     1.57127294 IU/ML
    42                 TRUE      PT 1.676128e+01     0.91352323 IU/ML
    43                 TRUE      PT 4.051410e+00     0.22080995 IU/ML
    44                 TRUE      PT 3.550618e+00     0.19351580 IU/ML
    45                 TRUE      PT 6.419504e+00     0.34987579 IU/ML
    46                 TRUE      PT 1.665922e+01     0.90796079 IU/ML
    47                 TRUE      PT 1.740970e+01     0.94886320 IU/ML
    48                 TRUE      PT 1.646815e+01     0.89754717 IU/ML
    49                 TRUE      PT 1.011236e+01     0.55114405 IU/ML
    50                 TRUE      PT 5.360000e-01     0.02921307 IU/ML
    51                 TRUE      PT 5.360000e-01     0.02921307 IU/ML
    52                 TRUE      PT 4.638202e+00     0.25279129 IU/ML
    53                 TRUE      PT 3.821589e+00     0.20828423 IU/ML
    54                 TRUE      PT 5.135951e+00     0.27991960 IU/ML
    55                 TRUE      PT 3.418229e+00     0.18630034 IU/ML
    56                 TRUE      PT 4.453207e+01     2.42708675 IU/ML
    57                 TRUE      PT 3.613811e+01     1.96959897 IU/ML
    58                 TRUE      PT 3.523752e+01     1.92051552 IU/ML
    59                 TRUE      PT 5.337719e+02    29.09163758 IU/ML
    60                 TRUE      PT 4.145139e+02    22.59183883 IU/ML
    61                 TRUE      PT 1.991238e+02    10.85264293 IU/ML
    62                 TRUE      PT 1.434999e+02     7.82103317 IU/ML
    63                 TRUE      PT 2.426124e+00     0.13222859 IU/ML
    64                 TRUE      PT 2.736522e+00     0.14914591 IU/ML
    65                 TRUE      PT 5.360000e-01     0.02921307 IU/ML
    66                 TRUE      PT 2.495638e+01     1.36017276 IU/ML
    67                 TRUE      PT 9.687427e+01     5.27984158 IU/ML
    68                 TRUE      PT 3.544540e+01     1.93184542 IU/ML
    69                 TRUE      PT 2.901240e+01     1.58123355 IU/ML
    70                 TRUE      PT 3.688192e+01     2.01013838 IU/ML
    71                 TRUE      PT 4.178865e+01     2.27756513 IU/ML
    72                 TRUE      PT 4.028783e+01     2.19576744 IU/ML
    73                 TRUE      PT 5.806193e+01     3.16449137 IU/ML
    74                 TRUE      PT 2.044619e+02    11.14358372 IU/ML
    75                 TRUE      PT 1.711338e+02     9.32713483 IU/ML
    76                 TRUE      PT 2.250771e+01     1.22671504 IU/ML
    77                 TRUE      PT 2.606545e+01     1.42061899 IU/ML
    78                 TRUE      PT 8.942424e-01     0.04873800 IU/ML
    79                 TRUE      PT 2.042297e+02    11.13092955 IU/ML
    80                 TRUE      PT 2.690013e+02    14.66110672 IU/ML
    81                 TRUE      PT 5.545653e+00     0.30224919 IU/ML
    82                 TRUE      PT 3.514701e+00     0.19155823 IU/ML
    83                 TRUE      PT 2.014169e+01     1.09776245 IU/ML
    84                 TRUE      PT 1.618016e+01     0.88185086 IU/ML
    85                 TRUE      PT 1.004895e+02     5.47687617 IU/ML
    86                 TRUE      PT 1.174317e+01     0.64002640 IU/ML
    87                 TRUE      PT 6.656937e+00     0.36281637 IU/ML
    88                 TRUE      PT 5.271166e+01     2.87289077 IU/ML
    89                 TRUE      PT 2.396038e+01     1.30588864 IU/ML
    90                 TRUE      PT 3.844803e+00     0.20954945 IU/ML
    91                 TRUE      PT 8.667540e+01     4.72398240 IU/ML
    92                 TRUE      PT 9.944042e+01     5.41970147 IU/ML
    93                 TRUE      PT 5.779780e+01     3.15009551 IU/ML
    94                 TRUE      PT 4.414412e+01     2.40594281 IU/ML
    95                 TRUE      PT 1.829349e+02     9.97031964 IU/ML
    96                 TRUE      PT 1.471305e+02     8.01890564 IU/ML
    97                 TRUE      PT 1.435148e+02     7.82184318 IU/ML
    98                 TRUE      PT 1.625990e+02     8.86197308 IU/ML
    99                 TRUE      PT 1.689009e+02     9.20543572 IU/ML
    100                TRUE      PT 1.298864e+02     7.07906798 IU/ML
    101                TRUE      PT 9.833932e+01     5.35968966 IU/ML
    102                TRUE      PT 4.074401e+01     2.22062993 IU/ML
    103                TRUE      PT 2.795453e+01     1.52357789 IU/ML
    104                TRUE      PT 2.132857e+01     1.16244994 IU/ML
    105                TRUE      PT 1.906822e+02    10.39256335 IU/ML
    106                TRUE      PT 1.829061e+02     9.96874747 IU/ML
    107                TRUE      PT 1.416498e+02     7.72019871 IU/ML
    108                TRUE      PT 1.116795e+02     6.08675823 IU/ML
    109                TRUE      PT 6.382116e+01     3.47838076 IU/ML
    110                TRUE      PT 6.693838e+01     3.64827565 IU/ML
    111                TRUE      PT 7.681929e+01     4.18680488 IU/ML
    112                TRUE      PT 1.688730e+02     9.20391366 IU/ML
    113                TRUE      PT 1.701837e+02     9.27535380 IU/ML
    114                TRUE      PT 8.618554e+01     4.69728443 IU/ML
    115                TRUE      PT 7.018579e+01     3.82526577 IU/ML
    116                TRUE      PT 8.912066e+01     4.85725419 IU/ML
    117                TRUE      PT 9.381605e+01     5.11316216 IU/ML
    118                TRUE      PT 2.483628e+01     1.35362715 IU/ML
    119                TRUE      PT 1.116587e+02     6.08562490 IU/ML
    120                TRUE      PT 1.332696e+02     7.26345933 IU/ML
    121                TRUE      PT 9.902579e+01     5.39710358 IU/ML
    122                TRUE      PT 6.589953e+01     3.59165598 IU/ML
    123                TRUE      PT 1.071943e+02     5.84230335 IU/ML
    124                TRUE      PT 9.663384e+00     0.52667375 IU/ML
    125                TRUE      PT 3.046665e+00     0.16604934 IU/ML
    126                TRUE      PT 2.561991e+01     1.39633657 IU/ML
    127                TRUE      PT 3.536885e+01     1.92767308 IU/ML
    128                TRUE      PT 2.291502e+01     1.24891414 IU/ML
    129                TRUE      PT 1.334159e+01     0.72714312 IU/ML
    130                TRUE      PT 6.387587e+01     3.48136253 IU/ML
    131                TRUE      PT 5.805951e+00     0.31643595 IU/ML
    132                TRUE      PT 1.966594e+00     0.10718333 IU/ML
    133                TRUE      PT 1.371755e+02     7.47634110 IU/ML
    134                TRUE      PT 8.998965e+01     4.90461599 IU/ML
    135                TRUE      PT 1.124450e+02     6.12847970 IU/ML
    136                TRUE      PT 7.609423e+01     4.14728768 IU/ML
    137                TRUE      PT 5.895038e+00     0.32129134 IU/ML
    138                TRUE      PT 6.333229e+00     0.34517365 IU/ML
    139                TRUE      PT 5.301253e+00     0.28892887 IU/ML
    140                TRUE      PT 7.992845e+01     4.35626048 IU/ML
    141                TRUE      PT 1.322196e+02     7.20623213 IU/ML
    142                TRUE      PT 1.287069e+02     7.01478560 IU/ML
    143                TRUE      PT 7.821953e+01     4.26312098 IU/ML
    144                TRUE      PT 3.774999e+00     0.20574500 IU/ML
    145                TRUE      PT 2.995838e+00     0.16327916 IU/ML
    146                TRUE      PT 2.867299e+00     0.15627354 IU/ML
    147                TRUE      PT 2.363176e+01     1.28797842 IU/ML
    148                TRUE      PT 1.600102e+01     0.87208731 IU/ML
    149                TRUE      PT 1.461359e+01     0.79647014 IU/ML
    150                TRUE      PT 1.257752e+01     0.68549982 IU/ML
    151                TRUE      PT 1.250248e+01     0.68141033 IU/ML
    152                TRUE      PT 5.855506e+00     0.31913680 IU/ML
    153                TRUE      PT 9.181613e+00     0.50041625 IU/ML
    154                TRUE      PT 2.309943e+01     1.25896521 IU/ML
    155                TRUE      PT 9.055533e+00     0.49354466 IU/ML
    156                TRUE      PT 1.030604e+01     0.56169992 IU/ML
    157                TRUE      PT 1.283892e+01     0.69974671 IU/ML
    158                TRUE      PT 1.201846e+01     0.65502996 IU/ML
    159                TRUE      PT 3.798321e+00     0.20701609 IU/ML
    160                TRUE      PT 1.145785e+01     0.62447577 IU/ML
    161                TRUE      PT 1.893884e+01     1.03220454 IU/ML
    162                TRUE      PT 4.574896e+01     2.49340964 IU/ML
    163                TRUE      PT 1.027561e+02     5.60041042 IU/ML
    164                TRUE      PT 3.381657e+01     1.84307076 IU/ML
    165                TRUE      PT 1.295780e+01     0.70622597 IU/ML
    166                TRUE      PT 3.798321e+00     0.20701609 IU/ML
    167                TRUE      PT 9.240235e+00     0.50361130 IU/ML
    168                TRUE      PT 5.360000e-01     0.02921307 IU/ML
    169                TRUE      PT 4.395969e+01     2.39589088 IU/ML
    170                TRUE      PT 4.423585e+01     2.41094249 IU/ML
    171                TRUE      PT 3.389009e+01     1.84707760 IU/ML
    172                TRUE      PT 2.492037e+01     1.35821005 IU/ML
    173                TRUE      PT 1.223861e+01     0.66702882 IU/ML
    174                TRUE      PT 2.153587e+01     1.17374816 IU/ML
    175                TRUE      PT 1.064620e+02     5.80239345 IU/ML
    176                TRUE      PT 1.164783e+02     6.34829906 IU/ML
    177                TRUE      PT 1.415592e+02     7.71525868 IU/ML
    178                TRUE      PT 8.791129e+01     4.79134103 IU/ML
    179                TRUE      PT 1.826877e+01     0.99568429 IU/ML
    180                TRUE      PT 1.437091e+01     0.78324334 IU/ML
    181                TRUE      PT 1.839418e+01     1.00251948 IU/ML
    182                TRUE      PT 3.514443e+01     1.91544177 IU/ML
    183                TRUE      PT 7.509077e+01     4.09259703 IU/ML
    184                TRUE      PT 7.084239e+01     3.86105190 IU/ML
    185                TRUE      PT 5.424624e+01     2.95652855 IU/ML
    186                TRUE      PT 1.023357e+01     0.55774987 IU/ML
    187                TRUE      PT 7.595422e+00     0.41396569 IU/ML
    188                TRUE      PT 3.704702e+00     0.20191367 IU/ML
    189                TRUE      PT 2.155469e+01     1.17477360 IU/ML
    190                TRUE      PT 6.269973e+01     3.41726062 IU/ML
    191                TRUE      PT 4.244389e+01     2.31327662 IU/ML
    192                TRUE      PT 1.543155e+01     0.84105042 IU/ML
    193                TRUE      PT 2.229580e+00     0.12151655 IU/ML
    194                TRUE      PT 2.970301e+00     0.16188733 IU/ML
    195                TRUE      PT 1.609199e+00     0.08770455 IU/ML
    196                TRUE      PT 1.825555e+01     0.99496392 IU/ML
    197                TRUE      PT 3.300308e+01     1.79873358 IU/ML
    198                TRUE      PT 2.724314e+01     1.48480542 IU/ML
    199                TRUE      PT 4.884167e+00     0.26619685 IU/ML
    200                TRUE      PT 3.019657e+01     1.64577368 IU/ML
    201                TRUE      PT 2.746527e+01     1.49691198 IU/ML
    202                TRUE      PT 2.919055e+01     1.59094337 IU/ML
    203                TRUE      PT 1.459129e+02     7.95254476 IU/ML
    204                TRUE      PT 2.727172e+02    14.86363460 IU/ML
    205                TRUE      PT 4.075844e+02    22.21416301 IU/ML
    206                TRUE      PT 1.678208e+02     9.14657064 IU/ML
    207                TRUE      PT 1.020937e+01     0.55643096 IU/ML
    208                TRUE      PT 1.151975e+01     0.62784956 IU/ML
    209                TRUE      PT 1.337095e+01     0.72874345 IU/ML
    210                TRUE      PT 3.134043e+01     1.70811604 IU/ML
    211                TRUE      PT 1.540787e+02     8.39759673 IU/ML
    212                TRUE      PT 1.262250e+02     6.87951327 IU/ML
    213                TRUE      PT 7.603800e+01     4.14422283 IU/ML
    214                TRUE      PT 2.169251e+01     1.18228510 IU/ML
    215                TRUE      PT 2.065205e+01     1.12557828 IU/ML
    216                TRUE      PT 2.410600e+01     1.31382533 IU/ML
    217                TRUE      PT 1.723144e+02     9.39148179 IU/ML
    218                TRUE      PT 1.407706e+02     7.67228073 IU/ML
    219                TRUE      PT 3.968596e+01     2.16296414 IU/ML
    220                TRUE      PT 2.751197e+01     1.49945725 IU/ML
    221                TRUE      PT 2.360127e+01     1.28631632 IU/ML
    222                TRUE      PT 3.265219e+01     1.77960983 IU/ML
    223                TRUE      PT 4.159569e+01     2.26704815 IU/ML
    224                TRUE      PT 7.345474e+01     4.00342993 IU/ML
    225                TRUE      PT 7.734049e+01     4.21521124 IU/ML
    226                TRUE      PT 6.973527e+01     3.80071177 IU/ML
    227                TRUE      PT 4.484448e+01     2.44411398 IU/ML
    228                TRUE      PT 1.940079e+01     1.05738162 IU/ML
    229                TRUE      PT 9.597456e+00     0.52308055 IU/ML
    230                TRUE      PT 1.087948e+01     0.59295350 IU/ML
    231                TRUE      PT 2.676206e+01     1.45858577 IU/ML
    232                TRUE      PT 4.140028e+01     2.25639777 IU/ML
    233                TRUE      PT 2.243961e+01     1.22300337 IU/ML
    234                TRUE      PT 7.505371e+00     0.40905771 IU/ML
    235                TRUE      PT 4.670559e+00     0.25455480 IU/ML
    236                TRUE      PT 2.580434e+01     1.40638789 IU/ML
    237                TRUE      PT 2.103225e+01     1.14629948 IU/ML
    238                TRUE      PT 2.126563e+01     1.15901953 IU/ML
    239                TRUE      PT 1.297263e+01     0.70703442 IU/ML
    240                TRUE      PT 9.933649e+00     0.54140374 IU/ML
    241                TRUE      PT 8.835470e+00     0.48155077 IU/ML
    242                TRUE      PT 9.024783e+01     4.91868692 IU/ML
    243                TRUE      PT 1.211634e+02     6.60364528 IU/ML
    244                TRUE      PT 5.811900e+01     3.16760167 IU/ML
    245                TRUE      PT 1.678167e+01     0.91463437 IU/ML
    246                TRUE      PT 7.433043e+00     0.40511568 IU/ML
    247                TRUE      PT 2.402022e+02    13.09150273 IU/ML
    248                TRUE      PT 4.231762e+00     0.23063950 IU/ML
    249                TRUE      PT 7.188113e+01     3.91766514 IU/ML
    250                TRUE      PT 1.580558e+02     8.61435564 IU/ML
    251                TRUE      PT 1.065845e+02     5.80906991 IU/ML
    252                TRUE      PT 4.548803e+01     2.47918867 IU/ML
    253                TRUE      PT 1.537562e+01     0.83800194 IU/ML
    254                TRUE      PT 1.791079e+01     0.97617404 IU/ML
    255                TRUE      PT 1.095040e+01     0.59681859 IU/ML
    256                TRUE      PT 9.736591e+01     5.30663653 IU/ML
    257                TRUE      PT 1.450357e+02     7.90473632 IU/ML
    258                TRUE      PT 1.085604e+02     5.91676096 IU/ML
    259                TRUE      PT 3.657272e+01     1.99328622 IU/ML
    260                TRUE      PT 2.599431e+01     1.41674185 IU/ML
    261                TRUE      PT 2.673853e+01     1.45730338 IU/ML
    262                TRUE      PT 2.973470e+01     1.62060080 IU/ML
    263                TRUE      PT 7.844682e+01     4.27550833 IU/ML
    264                TRUE      PT 1.985910e+02    10.82360502 IU/ML
    265                TRUE      PT 1.582781e+02     8.62647047 IU/ML
    266                TRUE      PT 1.062967e+02     5.79338157 IU/ML
    267                TRUE      PT 2.143544e+01     1.16827445 IU/ML
    268                TRUE      PT 1.991070e+01     1.08517287 IU/ML
    269                TRUE      PT 1.986568e+01     1.08271954 IU/ML
    270                TRUE      PT 1.462475e+02     7.97078041 IU/ML
    271                TRUE      PT 1.626313e+02     8.86373049 IU/ML
    272                TRUE      PT 1.405726e+02     7.66148966 IU/ML
    273                TRUE      PT 6.919271e+01     3.77114100 IU/ML
    274                TRUE      PT 1.842714e+01     1.00431571 IU/ML
    275                TRUE      PT 1.469892e+01     0.80112055 IU/ML
    276                TRUE      PT 1.418438e+01     0.77307707 IU/ML
    277                TRUE      PT 2.650881e+01     1.44478329 IU/ML
    278                TRUE      PT 4.347727e+01     2.36959803 IU/ML
    279                TRUE      PT 3.329199e+01     1.81448018 IU/ML
    280                TRUE      PT 2.161737e+01     1.17818975 IU/ML
    281                TRUE      PT 1.657742e+01     0.90350269 IU/ML
    282                TRUE      PT 1.471312e+01     0.80189472 IU/ML
    283                TRUE      PT 1.634493e+01     0.89083146 IU/ML
    284                TRUE      PT 2.308181e+02    12.58004938 IU/ML
    285                TRUE      PT 2.938809e+02    16.01709610 IU/ML
    286                TRUE      PT 1.799022e+02     9.80502883 IU/ML
    287                TRUE      PT 6.102118e+01     3.32577632 IU/ML
    288                TRUE      PT 1.891160e+00     0.10307198 IU/ML
    289                TRUE      PT 2.467442e+00     0.13448051 IU/ML
    290                TRUE      PT 5.360000e-01     0.02921307 IU/ML
    291                TRUE      PT 7.765449e+00     0.42323251 IU/ML
    292                TRUE      PT 3.849907e+01     2.09827627 IU/ML
    293                TRUE      PT 4.788488e+00     0.26098216 IU/ML
    294                TRUE      PT 1.805916e+02     9.84260193 IU/ML
    295                TRUE      PT 1.954764e+02    10.65385586 IU/ML
    296                TRUE      PT 1.731525e+02     9.43715476 IU/ML
    297                TRUE      PT 2.393218e+02    13.04351917 IU/ML
    298                TRUE      PT 2.012391e+02    10.96793431 IU/ML
    299                TRUE      PT 9.749316e+01     5.31357197 IU/ML
    300                TRUE      PT 9.446639e+01     5.14860711 IU/ML
    301                TRUE      PT 1.039093e+02     5.66326261 IU/ML
    302                TRUE      PT 2.251260e+02    12.26981730 IU/ML
    303                TRUE      PT 2.808927e+02    15.30921334 IU/ML
    304                TRUE      PT 2.154813e+02    11.74416023 IU/ML
    305                TRUE      PT 1.707925e+02     9.30853210 IU/ML
    306                TRUE      PT 5.360000e-01     0.02921307 IU/ML
    307                TRUE      PT 5.360000e-01     0.02921307 IU/ML
    308                TRUE      PT 5.360000e-01     0.02921307 IU/ML
    309                TRUE      PT 3.574487e+00     0.19481668 IU/ML
    310                TRUE      PT 5.360000e-01     0.02921307 IU/ML
    311                TRUE      PT 5.383180e+00     0.29339408 IU/ML
    312                TRUE      PT 1.736680e+00     0.09465256 IU/ML
    313                TRUE      PT 1.623514e+01     0.88484777 IU/ML
    314                TRUE      PT 1.359047e+01     0.74070760 IU/ML
    315                TRUE      PT 1.627634e+01     0.88709323 IU/ML
    316                TRUE      PT 1.023437e+02     5.57793682 IU/ML
    317                TRUE      PT 1.399996e+02     7.63025546 IU/ML
    318                TRUE      PT 1.158109e+02     6.31192769 IU/ML
    319                TRUE      PT 5.811900e+01     3.16760167 IU/ML
    320                TRUE      PT 7.738694e+00     0.42177432 IU/ML
    321                TRUE      PT 6.993442e+00     0.38115658 IU/ML
    322                TRUE      PT 3.478644e+00     0.18959303 IU/ML
    323                TRUE      PT 1.656378e+01     0.90275896 IU/ML
    324                TRUE      PT 5.288206e+01     2.88217808 IU/ML
    325                TRUE      PT 4.468982e+01     2.43568435 IU/ML
    326                TRUE      PT 2.099432e+01     1.14423265 IU/ML
    327                TRUE      PT 2.889154e+01     1.57464676 IU/ML
    328                TRUE      PT 2.841253e+01     1.54853942 IU/ML
    329                TRUE      PT 2.086143e+01     1.13698953 IU/ML
    330                TRUE      PT 3.855953e+01     2.10157135 IU/ML
    331                TRUE      PT 7.357354e+01     4.00990527 IU/ML
    332                TRUE      PT 5.363396e+01     2.92315807 IU/ML
    333                TRUE      PT 4.938502e+01     2.69158207 IU/ML
    334                TRUE      PT 5.360000e-01     0.02921307 IU/ML
    335                TRUE      PT 6.486963e+01     3.53552427 IU/ML
    336                TRUE      PT 1.099980e+02     5.99511003 IU/ML
    337                TRUE      PT 3.593304e+01     1.95842264 IU/ML
    338                TRUE      PT 2.182097e+01     1.18928641 IU/ML
    339                TRUE      PT 3.135649e+02    17.08991541 IU/ML
    340                TRUE      PT 1.626649e+02     8.86556036 IU/ML
    341                TRUE      PT 1.126574e+02     6.14005148 IU/ML
    342                TRUE      PT 1.235583e+02     6.73417266 IU/ML
    343                TRUE      PT 2.192793e+02    11.95116199 IU/ML
    344                TRUE      PT 1.401067e+02     7.63609789 IU/ML
    345                TRUE      PT 1.270450e+02     6.92420856 IU/ML
    346                TRUE      PT 6.455294e+01     3.51826438 IU/ML
    347                TRUE      PT 5.466709e+01     2.97946590 IU/ML
    348                TRUE      PT 3.128600e+01     1.70514930 IU/ML
    349                TRUE      PT 1.342113e+02     7.31478393 IU/ML
    350                TRUE      PT 1.723631e+02     9.39413524 IU/ML
    351                TRUE      PT 1.583413e+02     8.62991960 IU/ML
    352                TRUE      PT 1.506017e+02     8.20809541 IU/ML
    353                TRUE      PT 6.236853e+00     0.33992094 IU/ML
    354                TRUE      PT 8.250169e+00     0.44965073 IU/ML
    355                TRUE      PT 3.960090e+00     0.21583280 IU/ML
    356                TRUE      PT 4.716210e+01     2.57042853 IU/ML
    357                TRUE      PT 1.206523e+02     6.57578973 IU/ML
    358                TRUE      PT 9.477905e+01     5.16564782 IU/ML
    359                TRUE      PT 9.057606e+01     4.93657616 IU/ML
    360                TRUE      PT 1.520603e+02     8.28759185 IU/ML
    361                TRUE      PT 1.637493e+02     8.92466390 IU/ML
    362                TRUE      PT 1.690590e+02     9.21405607 IU/ML
    363                TRUE      PT 1.975188e+02    10.76516737 IU/ML
    364                TRUE      PT 2.025211e+02    11.03780589 IU/ML
    365                TRUE      PT 1.766280e+02     9.62657769 IU/ML
    366                TRUE      PT 1.500179e+02     8.17627776 IU/ML
    367                TRUE      PT 3.751622e+00     0.20447092 IU/ML
    368                TRUE      PT 1.427234e+00     0.07778712 IU/ML
    369                TRUE      PT 4.464164e+00     0.24330586 IU/ML
    370                TRUE      PT 6.028974e+01     3.28591154 IU/ML
    371                TRUE      PT 8.300331e+01     4.52384648 IU/ML
    372                TRUE      PT 7.641164e+01     4.16458701 IU/ML
    373                TRUE      PT 3.410543e+01     1.85881424 IU/ML
    374                TRUE      PT 1.598721e+01     0.87133475 IU/ML
    375                TRUE      PT 2.552458e+01     1.39114064 IU/ML
    376                TRUE      PT 3.062824e+01     1.66930048 IU/ML
    377                TRUE      PT 1.850249e+02    10.08422490 IU/ML
    378                TRUE      PT 2.233094e+02    12.17081132 IU/ML
    379                TRUE      PT 1.962044e+02    10.69353185 IU/ML
    380                TRUE      PT 8.720626e+01     4.75291554 IU/ML
    381                TRUE      PT 1.254752e+01     0.68386505 IU/ML
    382                TRUE      PT 2.148568e+01     1.17101229 IU/ML
    383                TRUE      PT 3.937680e+01     2.14611450 IU/ML
    384                TRUE      PT 1.058576e+02     5.76945289 IU/ML
    385                TRUE      PT 3.529587e+02    19.23695317 IU/ML
    386                TRUE      PT 1.533711e+04     0.02888606 IU/ML
    387                TRUE      PT 1.920847e+02    10.46900269 IU/ML
    388                TRUE      PT 5.360000e-01     0.02921307 IU/ML
    389                TRUE      PT 5.360000e-01     0.02921307 IU/ML
    390                TRUE      PT 5.360000e-01     0.02921307 IU/ML
    391                TRUE      PT 2.118372e+01     1.15455529 IU/ML
    392                TRUE      PT 3.155439e+01     1.71977739 IU/ML
    393                TRUE      PT 4.176455e+01     2.27625121 IU/ML
    394                TRUE      PT 1.048244e+01     0.57131410 IU/ML
        lower_limit_of_detection
    1                       0.53
    2                       0.53
    3                       0.53
    4                       0.53
    5                       0.53
    6                       0.53
    7                       0.53
    8                       0.53
    9                       0.53
    10                      0.53
    11                      0.53
    12                      0.53
    13                      0.53
    14                      0.53
    15                      0.53
    16                      0.53
    17                      0.53
    18                      0.53
    19                      0.53
    20                      0.53
    21                      0.53
    22                      0.53
    23                      0.53
    24                      0.53
    25                      0.53
    26                      0.53
    27                      0.53
    28                      0.53
    29                      0.53
    30                      0.53
    31                      0.53
    32                      0.53
    33                      0.53
    34                      0.53
    35                      0.53
    36                      0.53
    37                      0.53
    38                      0.53
    39                      0.53
    40                      0.53
    41                      0.53
    42                      0.53
    43                      0.53
    44                      0.53
    45                      0.53
    46                      0.53
    47                      0.53
    48                      0.53
    49                      0.53
    50                      0.53
    51                      0.53
    52                      0.53
    53                      0.53
    54                      0.53
    55                      0.53
    56                      0.53
    57                      0.53
    58                      0.53
    59                      0.53
    60                      0.53
    61                      0.53
    62                      0.53
    63                      0.53
    64                      0.53
    65                      0.53
    66                      0.53
    67                      0.53
    68                      0.53
    69                      0.53
    70                      0.53
    71                      0.53
    72                      0.53
    73                      0.53
    74                      0.53
    75                      0.53
    76                      0.53
    77                      0.53
    78                      0.53
    79                      0.53
    80                      0.53
    81                      0.53
    82                      0.53
    83                      0.53
    84                      0.53
    85                      0.53
    86                      0.53
    87                      0.53
    88                      0.53
    89                      0.53
    90                      0.53
    91                      0.53
    92                      0.53
    93                      0.53
    94                      0.53
    95                      0.53
    96                      0.53
    97                      0.53
    98                      0.53
    99                      0.53
    100                     0.53
    101                     0.53
    102                     0.53
    103                     0.53
    104                     0.53
    105                     0.53
    106                     0.53
    107                     0.53
    108                     0.53
    109                     0.53
    110                     0.53
    111                     0.53
    112                     0.53
    113                     0.53
    114                     0.53
    115                     0.53
    116                     0.53
    117                     0.53
    118                     0.53
    119                     0.53
    120                     0.53
    121                     0.53
    122                     0.53
    123                     0.53
    124                     0.53
    125                     0.53
    126                     0.53
    127                     0.53
    128                     0.53
    129                     0.53
    130                     0.53
    131                     0.53
    132                     0.53
    133                     0.53
    134                     0.53
    135                     0.53
    136                     0.53
    137                     0.53
    138                     0.53
    139                     0.53
    140                     0.53
    141                     0.53
    142                     0.53
    143                     0.53
    144                     0.53
    145                     0.53
    146                     0.53
    147                     0.53
    148                     0.53
    149                     0.53
    150                     0.53
    151                     0.53
    152                     0.53
    153                     0.53
    154                     0.53
    155                     0.53
    156                     0.53
    157                     0.53
    158                     0.53
    159                     0.53
    160                     0.53
    161                     0.53
    162                     0.53
    163                     0.53
    164                     0.53
    165                     0.53
    166                     0.53
    167                     0.53
    168                     0.53
    169                     0.53
    170                     0.53
    171                     0.53
    172                     0.53
    173                     0.53
    174                     0.53
    175                     0.53
    176                     0.53
    177                     0.53
    178                     0.53
    179                     0.53
    180                     0.53
    181                     0.53
    182                     0.53
    183                     0.53
    184                     0.53
    185                     0.53
    186                     0.53
    187                     0.53
    188                     0.53
    189                     0.53
    190                     0.53
    191                     0.53
    192                     0.53
    193                     0.53
    194                     0.53
    195                     0.53
    196                     0.53
    197                     0.53
    198                     0.53
    199                     0.53
    200                     0.53
    201                     0.53
    202                     0.53
    203                     0.53
    204                     0.53
    205                     0.53
    206                     0.53
    207                     0.53
    208                     0.53
    209                     0.53
    210                     0.53
    211                     0.53
    212                     0.53
    213                     0.53
    214                     0.53
    215                     0.53
    216                     0.53
    217                     0.53
    218                     0.53
    219                     0.53
    220                     0.53
    221                     0.53
    222                     0.53
    223                     0.53
    224                     0.53
    225                     0.53
    226                     0.53
    227                     0.53
    228                     0.53
    229                     0.53
    230                     0.53
    231                     0.53
    232                     0.53
    233                     0.53
    234                     0.53
    235                     0.53
    236                     0.53
    237                     0.53
    238                     0.53
    239                     0.53
    240                     0.53
    241                     0.53
    242                     0.53
    243                     0.53
    244                     0.53
    245                     0.53
    246                     0.53
    247                     0.53
    248                     0.53
    249                     0.53
    250                     0.53
    251                     0.53
    252                     0.53
    253                     0.53
    254                     0.53
    255                     0.53
    256                     0.53
    257                     0.53
    258                     0.53
    259                     0.53
    260                     0.53
    261                     0.53
    262                     0.53
    263                     0.53
    264                     0.53
    265                     0.53
    266                     0.53
    267                     0.53
    268                     0.53
    269                     0.53
    270                     0.53
    271                     0.53
    272                     0.53
    273                     0.53
    274                     0.53
    275                     0.53
    276                     0.53
    277                     0.53
    278                     0.53
    279                     0.53
    280                     0.53
    281                     0.53
    282                     0.53
    283                     0.53
    284                     0.53
    285                     0.53
    286                     0.53
    287                     0.53
    288                     0.53
    289                     0.53
    290                     0.53
    291                     0.53
    292                     0.53
    293                     0.53
    294                     0.53
    295                     0.53
    296                     0.53
    297                     0.53
    298                     0.53
    299                     0.53
    300                     0.53
    301                     0.53
    302                     0.53
    303                     0.53
    304                     0.53
    305                     0.53
    306                     0.53
    307                     0.53
    308                     0.53
    309                     0.53
    310                     0.53
    311                     0.53
    312                     0.53
    313                     0.53
    314                     0.53
    315                     0.53
    316                     0.53
    317                     0.53
    318                     0.53
    319                     0.53
    320                     0.53
    321                     0.53
    322                     0.53
    323                     0.53
    324                     0.53
    325                     0.53
    326                     0.53
    327                     0.53
    328                     0.53
    329                     0.53
    330                     0.53
    331                     0.53
    332                     0.53
    333                     0.53
    334                     0.53
    335                     0.53
    336                     0.53
    337                     0.53
    338                     0.53
    339                     0.53
    340                     0.53
    341                     0.53
    342                     0.53
    343                     0.53
    344                     0.53
    345                     0.53
    346                     0.53
    347                     0.53
    348                     0.53
    349                     0.53
    350                     0.53
    351                     0.53
    352                     0.53
    353                     0.53
    354                     0.53
    355                     0.53
    356                     0.53
    357                     0.53
    358                     0.53
    359                     0.53
    360                     0.53
    361                     0.53
    362                     0.53
    363                     0.53
    364                     0.53
    365                     0.53
    366                     0.53
    367                     0.53
    368                     0.53
    369                     0.53
    370                     0.53
    371                     0.53
    372                     0.53
    373                     0.53
    374                     0.53
    375                     0.53
    376                     0.53
    377                     0.53
    378                     0.53
    379                     0.53
    380                     0.53
    381                     0.53
    382                     0.53
    383                     0.53
    384                     0.53
    385                     0.53
    386                     0.53
    387                     0.53
    388                     0.53
    389                     0.53
    390                     0.53
    391                     0.53
    392                     0.53
    393                     0.53
    394                     0.53

``` r
pt_2021
```

        subject_id infancy_vac biological_sex              ethnicity
    1           61          wP         Female Not Hispanic or Latino
    2           61          wP         Female Not Hispanic or Latino
    3           61          wP         Female Not Hispanic or Latino
    4           61          wP         Female Not Hispanic or Latino
    5           61          wP         Female Not Hispanic or Latino
    6           61          wP         Female Not Hispanic or Latino
    7           61          wP         Female Not Hispanic or Latino
    8           62          wP         Female Not Hispanic or Latino
    9           62          wP         Female Not Hispanic or Latino
    10          62          wP         Female Not Hispanic or Latino
    11          62          wP         Female Not Hispanic or Latino
    12          62          wP         Female Not Hispanic or Latino
    13          62          wP         Female Not Hispanic or Latino
    14          62          wP         Female Not Hispanic or Latino
    15          63          wP         Female Not Hispanic or Latino
    16          63          wP         Female Not Hispanic or Latino
    17          63          wP         Female Not Hispanic or Latino
    18          63          wP         Female Not Hispanic or Latino
    19          63          wP         Female Not Hispanic or Latino
    20          63          wP         Female Not Hispanic or Latino
    21          63          wP         Female Not Hispanic or Latino
    22          64          wP           Male Not Hispanic or Latino
    23          64          wP           Male Not Hispanic or Latino
    24          64          wP           Male Not Hispanic or Latino
    25          64          wP           Male Not Hispanic or Latino
    26          64          wP           Male Not Hispanic or Latino
    27          64          wP           Male Not Hispanic or Latino
    28          64          wP           Male Not Hispanic or Latino
    29          65          wP           Male Not Hispanic or Latino
    30          65          wP           Male Not Hispanic or Latino
    31          65          wP           Male Not Hispanic or Latino
    32          65          wP           Male Not Hispanic or Latino
    33          65          wP           Male Not Hispanic or Latino
    34          65          wP           Male Not Hispanic or Latino
    35          65          wP           Male Not Hispanic or Latino
    36          66          wP         Female Not Hispanic or Latino
    37          66          wP         Female Not Hispanic or Latino
    38          66          wP         Female Not Hispanic or Latino
    39          66          wP         Female Not Hispanic or Latino
    40          66          wP         Female Not Hispanic or Latino
    41          66          wP         Female Not Hispanic or Latino
    42          66          wP         Female Not Hispanic or Latino
    43          67          wP         Female     Hispanic or Latino
    44          67          wP         Female     Hispanic or Latino
    45          67          wP         Female     Hispanic or Latino
    46          67          wP         Female     Hispanic or Latino
    47          67          wP         Female     Hispanic or Latino
    48          67          wP         Female     Hispanic or Latino
    49          67          wP         Female     Hispanic or Latino
    50          68          wP           Male     Hispanic or Latino
    51          68          wP           Male     Hispanic or Latino
    52          68          wP           Male     Hispanic or Latino
    53          68          wP           Male     Hispanic or Latino
    54          68          wP           Male     Hispanic or Latino
    55          68          wP           Male     Hispanic or Latino
    56          68          wP           Male     Hispanic or Latino
    57          69          wP         Female     Hispanic or Latino
    58          69          wP         Female     Hispanic or Latino
    59          69          wP         Female     Hispanic or Latino
    60          69          wP         Female     Hispanic or Latino
    61          69          wP         Female     Hispanic or Latino
    62          69          wP         Female     Hispanic or Latino
    63          69          wP         Female     Hispanic or Latino
    64          70          aP           Male Not Hispanic or Latino
    65          70          aP           Male Not Hispanic or Latino
    66          70          aP           Male Not Hispanic or Latino
    67          70          aP           Male Not Hispanic or Latino
    68          70          aP           Male Not Hispanic or Latino
    69          70          aP           Male Not Hispanic or Latino
    70          70          aP           Male Not Hispanic or Latino
    71          71          aP         Female Not Hispanic or Latino
    72          71          aP         Female Not Hispanic or Latino
    73          71          aP         Female Not Hispanic or Latino
    74          71          aP         Female Not Hispanic or Latino
    75          71          aP         Female Not Hispanic or Latino
    76          71          aP         Female Not Hispanic or Latino
    77          71          aP         Female Not Hispanic or Latino
    78          72          wP         Female Not Hispanic or Latino
    79          72          wP         Female Not Hispanic or Latino
    80          72          wP         Female Not Hispanic or Latino
    81          72          wP         Female Not Hispanic or Latino
    82          72          wP         Female Not Hispanic or Latino
    83          72          wP         Female Not Hispanic or Latino
    84          72          wP         Female Not Hispanic or Latino
    85          73          wP         Female Not Hispanic or Latino
    86          73          wP         Female Not Hispanic or Latino
    87          73          wP         Female Not Hispanic or Latino
    88          73          wP         Female Not Hispanic or Latino
    89          73          wP         Female Not Hispanic or Latino
    90          73          wP         Female Not Hispanic or Latino
    91          73          wP         Female Not Hispanic or Latino
    92          74          wP         Female Not Hispanic or Latino
    93          74          wP         Female Not Hispanic or Latino
    94          74          wP         Female Not Hispanic or Latino
    95          74          wP         Female Not Hispanic or Latino
    96          74          wP         Female Not Hispanic or Latino
    97          74          wP         Female Not Hispanic or Latino
    98          74          wP         Female Not Hispanic or Latino
    99          75          aP         Female Not Hispanic or Latino
    100         75          aP         Female Not Hispanic or Latino
    101         75          aP         Female Not Hispanic or Latino
    102         75          aP         Female Not Hispanic or Latino
    103         75          aP         Female Not Hispanic or Latino
    104         75          aP         Female Not Hispanic or Latino
    105         75          aP         Female Not Hispanic or Latino
    106         76          aP         Female Not Hispanic or Latino
    107         76          aP         Female Not Hispanic or Latino
    108         76          aP         Female Not Hispanic or Latino
    109         76          aP         Female Not Hispanic or Latino
    110         76          aP         Female Not Hispanic or Latino
    111         76          aP         Female Not Hispanic or Latino
    112         76          aP         Female Not Hispanic or Latino
    113         77          wP           Male Not Hispanic or Latino
    114         77          wP           Male Not Hispanic or Latino
    115         77          wP           Male Not Hispanic or Latino
    116         77          wP           Male Not Hispanic or Latino
    117         77          wP           Male Not Hispanic or Latino
    118         77          wP           Male Not Hispanic or Latino
    119         77          wP           Male Not Hispanic or Latino
    120         78          wP         Female Not Hispanic or Latino
    121         78          wP         Female Not Hispanic or Latino
    122         78          wP         Female Not Hispanic or Latino
    123         78          wP         Female Not Hispanic or Latino
    124         78          wP         Female Not Hispanic or Latino
    125         78          wP         Female Not Hispanic or Latino
    126         78          wP         Female Not Hispanic or Latino
    127         79          wP           Male Not Hispanic or Latino
    128         79          wP           Male Not Hispanic or Latino
    129         79          wP           Male Not Hispanic or Latino
    130         79          wP           Male Not Hispanic or Latino
    131         79          wP           Male Not Hispanic or Latino
    132         79          wP           Male Not Hispanic or Latino
    133         79          wP           Male Not Hispanic or Latino
    134         80          wP         Female Not Hispanic or Latino
    135         80          wP         Female Not Hispanic or Latino
    136         80          wP         Female Not Hispanic or Latino
    137         80          wP         Female Not Hispanic or Latino
    138         80          wP         Female Not Hispanic or Latino
    139         80          wP         Female Not Hispanic or Latino
    140         80          wP         Female Not Hispanic or Latino
    141         81          wP           Male Not Hispanic or Latino
    142         81          wP           Male Not Hispanic or Latino
    143         81          wP           Male Not Hispanic or Latino
    144         81          wP           Male Not Hispanic or Latino
    145         81          wP           Male Not Hispanic or Latino
    146         81          wP           Male Not Hispanic or Latino
    147         81          wP           Male Not Hispanic or Latino
    148         83          aP         Female Not Hispanic or Latino
    149         83          aP         Female Not Hispanic or Latino
    150         83          aP         Female Not Hispanic or Latino
    151         83          aP         Female Not Hispanic or Latino
    152         83          aP         Female Not Hispanic or Latino
    153         83          aP         Female Not Hispanic or Latino
    154         83          aP         Female Not Hispanic or Latino
    155         84          aP         Female Not Hispanic or Latino
    156         84          aP         Female Not Hispanic or Latino
    157         84          aP         Female Not Hispanic or Latino
    158         84          aP         Female Not Hispanic or Latino
    159         84          aP         Female Not Hispanic or Latino
    160         84          aP         Female Not Hispanic or Latino
    161         84          aP         Female Not Hispanic or Latino
    162         85          aP         Female     Hispanic or Latino
    163         85          aP         Female     Hispanic or Latino
    164         85          aP         Female     Hispanic or Latino
    165         85          aP         Female     Hispanic or Latino
    166         85          aP         Female     Hispanic or Latino
    167         85          aP         Female     Hispanic or Latino
    168         85          aP         Female     Hispanic or Latino
    169         86          aP         Female Not Hispanic or Latino
    170         86          aP         Female Not Hispanic or Latino
    171         86          aP         Female Not Hispanic or Latino
    172         86          aP         Female Not Hispanic or Latino
    173         86          aP         Female Not Hispanic or Latino
    174         86          aP         Female Not Hispanic or Latino
    175         86          aP         Female Not Hispanic or Latino
    176         89          aP         Female Not Hispanic or Latino
    177         89          aP         Female Not Hispanic or Latino
    178         89          aP         Female Not Hispanic or Latino
    179         89          aP         Female Not Hispanic or Latino
    180         89          aP         Female Not Hispanic or Latino
    181         89          aP         Female Not Hispanic or Latino
    182         89          aP         Female Not Hispanic or Latino
    183         90          aP         Female Not Hispanic or Latino
    184         90          aP         Female Not Hispanic or Latino
    185         90          aP         Female Not Hispanic or Latino
    186         90          aP         Female Not Hispanic or Latino
    187         90          aP         Female Not Hispanic or Latino
    188         90          aP         Female Not Hispanic or Latino
    189         90          aP         Female Not Hispanic or Latino
    190         91          aP           Male                Unknown
    191         91          aP           Male                Unknown
    192         91          aP           Male                Unknown
    193         91          aP           Male                Unknown
    194         91          aP           Male                Unknown
    195         91          aP           Male                Unknown
    196         91          aP           Male                Unknown
    197         92          aP         Female     Hispanic or Latino
    198         92          aP         Female     Hispanic or Latino
    199         92          aP         Female     Hispanic or Latino
    200         92          aP         Female     Hispanic or Latino
    201         92          aP         Female     Hispanic or Latino
    202         92          aP         Female     Hispanic or Latino
    203         92          aP         Female     Hispanic or Latino
    204         93          aP         Female Not Hispanic or Latino
    205         93          aP         Female Not Hispanic or Latino
    206         93          aP         Female Not Hispanic or Latino
    207         93          aP         Female Not Hispanic or Latino
    208         93          aP         Female Not Hispanic or Latino
    209         93          aP         Female Not Hispanic or Latino
    210         93          aP         Female Not Hispanic or Latino
    211         94          aP           Male Not Hispanic or Latino
    212         94          aP           Male Not Hispanic or Latino
    213         94          aP           Male Not Hispanic or Latino
    214         94          aP           Male Not Hispanic or Latino
    215         94          aP           Male Not Hispanic or Latino
    216         94          aP           Male Not Hispanic or Latino
    217         94          aP           Male Not Hispanic or Latino
    218         95          aP         Female     Hispanic or Latino
    219         95          aP         Female     Hispanic or Latino
    220         95          aP         Female     Hispanic or Latino
    221         95          aP         Female     Hispanic or Latino
    222         95          aP         Female     Hispanic or Latino
    223         95          aP         Female     Hispanic or Latino
    224         95          aP         Female     Hispanic or Latino
    225         96          aP           Male     Hispanic or Latino
    226         96          aP           Male     Hispanic or Latino
    227         96          aP           Male     Hispanic or Latino
    228         96          aP           Male     Hispanic or Latino
    229         96          aP           Male     Hispanic or Latino
    230         96          aP           Male     Hispanic or Latino
    231         96          aP           Male     Hispanic or Latino
                                             race year_of_birth date_of_boost
    1                     Unknown or Not Reported    1987-01-01    2019-04-08
    2                     Unknown or Not Reported    1987-01-01    2019-04-08
    3                     Unknown or Not Reported    1987-01-01    2019-04-08
    4                     Unknown or Not Reported    1987-01-01    2019-04-08
    5                     Unknown or Not Reported    1987-01-01    2019-04-08
    6                     Unknown or Not Reported    1987-01-01    2019-04-08
    7                     Unknown or Not Reported    1987-01-01    2019-04-08
    8                                       Asian    1993-01-01    2018-11-26
    9                                       Asian    1993-01-01    2018-11-26
    10                                      Asian    1993-01-01    2018-11-26
    11                                      Asian    1993-01-01    2018-11-26
    12                                      Asian    1993-01-01    2018-11-26
    13                                      Asian    1993-01-01    2018-11-26
    14                                      Asian    1993-01-01    2018-11-26
    15                                      White    1995-01-01    2018-11-26
    16                                      White    1995-01-01    2018-11-26
    17                                      White    1995-01-01    2018-11-26
    18                                      White    1995-01-01    2018-11-26
    19                                      White    1995-01-01    2018-11-26
    20                                      White    1995-01-01    2018-11-26
    21                                      White    1995-01-01    2018-11-26
    22                                      Asian    1993-01-01    2018-11-26
    23                                      Asian    1993-01-01    2018-11-26
    24                                      Asian    1993-01-01    2018-11-26
    25                                      Asian    1993-01-01    2018-11-26
    26                                      Asian    1993-01-01    2018-11-26
    27                                      Asian    1993-01-01    2018-11-26
    28                                      Asian    1993-01-01    2018-11-26
    29                                      White    1990-01-01    2018-12-03
    30                                      White    1990-01-01    2018-12-03
    31                                      White    1990-01-01    2018-12-03
    32                                      White    1990-01-01    2018-12-03
    33                                      White    1990-01-01    2018-12-03
    34                                      White    1990-01-01    2018-12-03
    35                                      White    1990-01-01    2018-12-03
    36                  Black or African American    1976-01-01    2018-12-03
    37                  Black or African American    1976-01-01    2018-12-03
    38                  Black or African American    1976-01-01    2018-12-03
    39                  Black or African American    1976-01-01    2018-12-03
    40                  Black or African American    1976-01-01    2018-12-03
    41                  Black or African American    1976-01-01    2018-12-03
    42                  Black or African American    1976-01-01    2018-12-03
    43                                      White    1972-01-01    2019-01-28
    44                                      White    1972-01-01    2019-01-28
    45                                      White    1972-01-01    2019-01-28
    46                                      White    1972-01-01    2019-01-28
    47                                      White    1972-01-01    2019-01-28
    48                                      White    1972-01-01    2019-01-28
    49                                      White    1972-01-01    2019-01-28
    50                                      White    1972-01-01    2019-01-28
    51                                      White    1972-01-01    2019-01-28
    52                                      White    1972-01-01    2019-01-28
    53                                      White    1972-01-01    2019-01-28
    54                                      White    1972-01-01    2019-01-28
    55                                      White    1972-01-01    2019-01-28
    56                                      White    1972-01-01    2019-01-28
    57                                      White    1990-01-01    2019-01-28
    58                                      White    1990-01-01    2019-01-28
    59                                      White    1990-01-01    2019-01-28
    60                                      White    1990-01-01    2019-01-28
    61                                      White    1990-01-01    2019-01-28
    62                                      White    1990-01-01    2019-01-28
    63                                      White    1990-01-01    2019-01-28
    64              American Indian/Alaska Native    1998-01-01    2019-01-28
    65              American Indian/Alaska Native    1998-01-01    2019-01-28
    66              American Indian/Alaska Native    1998-01-01    2019-01-28
    67              American Indian/Alaska Native    1998-01-01    2019-01-28
    68              American Indian/Alaska Native    1998-01-01    2019-01-28
    69              American Indian/Alaska Native    1998-01-01    2019-01-28
    70              American Indian/Alaska Native    1998-01-01    2019-01-28
    71                                      White    1998-01-01    2019-01-28
    72                                      White    1998-01-01    2019-01-28
    73                                      White    1998-01-01    2019-01-28
    74                                      White    1998-01-01    2019-01-28
    75                                      White    1998-01-01    2019-01-28
    76                                      White    1998-01-01    2019-01-28
    77                                      White    1998-01-01    2019-01-28
    78                                      White    1991-01-01    2019-02-25
    79                                      White    1991-01-01    2019-02-25
    80                                      White    1991-01-01    2019-02-25
    81                                      White    1991-01-01    2019-02-25
    82                                      White    1991-01-01    2019-02-25
    83                                      White    1991-01-01    2019-02-25
    84                                      White    1991-01-01    2019-02-25
    85                                      White    1995-01-01    2019-02-25
    86                                      White    1995-01-01    2019-02-25
    87                                      White    1995-01-01    2019-02-25
    88                                      White    1995-01-01    2019-02-25
    89                                      White    1995-01-01    2019-02-25
    90                                      White    1995-01-01    2019-02-25
    91                                      White    1995-01-01    2019-02-25
    92                                      White    1995-01-01    2019-02-25
    93                                      White    1995-01-01    2019-02-25
    94                                      White    1995-01-01    2019-02-25
    95                                      White    1995-01-01    2019-02-25
    96                                      White    1995-01-01    2019-02-25
    97                                      White    1995-01-01    2019-02-25
    98                                      White    1995-01-01    2019-02-25
    99  Native Hawaiian or Other Pacific Islander    1998-01-01    2019-02-25
    100 Native Hawaiian or Other Pacific Islander    1998-01-01    2019-02-25
    101 Native Hawaiian or Other Pacific Islander    1998-01-01    2019-02-25
    102 Native Hawaiian or Other Pacific Islander    1998-01-01    2019-02-25
    103 Native Hawaiian or Other Pacific Islander    1998-01-01    2019-02-25
    104 Native Hawaiian or Other Pacific Islander    1998-01-01    2019-02-25
    105 Native Hawaiian or Other Pacific Islander    1998-01-01    2019-02-25
    106                                     Asian    1998-01-01    2019-02-25
    107                                     Asian    1998-01-01    2019-02-25
    108                                     Asian    1998-01-01    2019-02-25
    109                                     Asian    1998-01-01    2019-02-25
    110                                     Asian    1998-01-01    2019-02-25
    111                                     Asian    1998-01-01    2019-02-25
    112                                     Asian    1998-01-01    2019-02-25
    113                                     White    1988-01-01    2019-03-18
    114                                     White    1988-01-01    2019-03-18
    115                                     White    1988-01-01    2019-03-18
    116                                     White    1988-01-01    2019-03-18
    117                                     White    1988-01-01    2019-03-18
    118                                     White    1988-01-01    2019-03-18
    119                                     White    1988-01-01    2019-03-18
    120                                     White    1993-01-01    2019-03-18
    121                                     White    1993-01-01    2019-03-18
    122                                     White    1993-01-01    2019-03-18
    123                                     White    1993-01-01    2019-03-18
    124                                     White    1993-01-01    2019-03-18
    125                                     White    1993-01-01    2019-03-18
    126                                     White    1993-01-01    2019-03-18
    127                                     White    1987-01-01    2019-03-18
    128                                     White    1987-01-01    2019-03-18
    129                                     White    1987-01-01    2019-03-18
    130                                     White    1987-01-01    2019-03-18
    131                                     White    1987-01-01    2019-03-18
    132                                     White    1987-01-01    2019-03-18
    133                                     White    1987-01-01    2019-03-18
    134                                     Asian    1992-01-01    2019-03-18
    135                                     Asian    1992-01-01    2019-03-18
    136                                     Asian    1992-01-01    2019-03-18
    137                                     Asian    1992-01-01    2019-03-18
    138                                     Asian    1992-01-01    2019-03-18
    139                                     Asian    1992-01-01    2019-03-18
    140                                     Asian    1992-01-01    2019-03-18
    141                                     White    1993-01-01    2019-03-18
    142                                     White    1993-01-01    2019-03-18
    143                                     White    1993-01-01    2019-03-18
    144                                     White    1993-01-01    2019-03-18
    145                                     White    1993-01-01    2019-03-18
    146                                     White    1993-01-01    2019-03-18
    147                                     White    1993-01-01    2019-03-18
    148                                     White    1999-01-01    2019-04-08
    149                                     White    1999-01-01    2019-04-08
    150                                     White    1999-01-01    2019-04-08
    151                                     White    1999-01-01    2019-04-08
    152                                     White    1999-01-01    2019-04-08
    153                                     White    1999-01-01    2019-04-08
    154                                     White    1999-01-01    2019-04-08
    155                        More Than One Race    1997-01-01    2019-04-08
    156                        More Than One Race    1997-01-01    2019-04-08
    157                        More Than One Race    1997-01-01    2019-04-08
    158                        More Than One Race    1997-01-01    2019-04-08
    159                        More Than One Race    1997-01-01    2019-04-08
    160                        More Than One Race    1997-01-01    2019-04-08
    161                        More Than One Race    1997-01-01    2019-04-08
    162                                     White    2000-01-01    2019-04-29
    163                                     White    2000-01-01    2019-04-29
    164                                     White    2000-01-01    2019-04-29
    165                                     White    2000-01-01    2019-04-29
    166                                     White    2000-01-01    2019-04-29
    167                                     White    2000-01-01    2019-04-29
    168                                     White    2000-01-01    2019-04-29
    169                                     Asian    1998-01-01    2019-04-29
    170                                     Asian    1998-01-01    2019-04-29
    171                                     Asian    1998-01-01    2019-04-29
    172                                     Asian    1998-01-01    2019-04-29
    173                                     Asian    1998-01-01    2019-04-29
    174                                     Asian    1998-01-01    2019-04-29
    175                                     Asian    1998-01-01    2019-04-29
    176                                     Asian    1997-01-01    2019-06-03
    177                                     Asian    1997-01-01    2019-06-03
    178                                     Asian    1997-01-01    2019-06-03
    179                                     Asian    1997-01-01    2019-06-03
    180                                     Asian    1997-01-01    2019-06-03
    181                                     Asian    1997-01-01    2019-06-03
    182                                     Asian    1997-01-01    2019-06-03
    183                                     Asian    1999-01-01    2019-06-03
    184                                     Asian    1999-01-01    2019-06-03
    185                                     Asian    1999-01-01    2019-06-03
    186                                     Asian    1999-01-01    2019-06-03
    187                                     Asian    1999-01-01    2019-06-03
    188                                     Asian    1999-01-01    2019-06-03
    189                                     Asian    1999-01-01    2019-06-03
    190                   Unknown or Not Reported    1998-01-01    2019-06-03
    191                   Unknown or Not Reported    1998-01-01    2019-06-03
    192                   Unknown or Not Reported    1998-01-01    2019-06-03
    193                   Unknown or Not Reported    1998-01-01    2019-06-03
    194                   Unknown or Not Reported    1998-01-01    2019-06-03
    195                   Unknown or Not Reported    1998-01-01    2019-06-03
    196                   Unknown or Not Reported    1998-01-01    2019-06-03
    197                                     White    2000-01-01    2019-06-24
    198                                     White    2000-01-01    2019-06-24
    199                                     White    2000-01-01    2019-06-24
    200                                     White    2000-01-01    2019-06-24
    201                                     White    2000-01-01    2019-06-24
    202                                     White    2000-01-01    2019-06-24
    203                                     White    2000-01-01    2019-06-24
    204                        More Than One Race    1996-01-01    2019-06-24
    205                        More Than One Race    1996-01-01    2019-06-24
    206                        More Than One Race    1996-01-01    2019-06-24
    207                        More Than One Race    1996-01-01    2019-06-24
    208                        More Than One Race    1996-01-01    2019-06-24
    209                        More Than One Race    1996-01-01    2019-06-24
    210                        More Than One Race    1996-01-01    2019-06-24
    211                   Unknown or Not Reported    1999-01-01    2019-06-24
    212                   Unknown or Not Reported    1999-01-01    2019-06-24
    213                   Unknown or Not Reported    1999-01-01    2019-06-24
    214                   Unknown or Not Reported    1999-01-01    2019-06-24
    215                   Unknown or Not Reported    1999-01-01    2019-06-24
    216                   Unknown or Not Reported    1999-01-01    2019-06-24
    217                   Unknown or Not Reported    1999-01-01    2019-06-24
    218                   Unknown or Not Reported    1998-01-01    2019-06-24
    219                   Unknown or Not Reported    1998-01-01    2019-06-24
    220                   Unknown or Not Reported    1998-01-01    2019-06-24
    221                   Unknown or Not Reported    1998-01-01    2019-06-24
    222                   Unknown or Not Reported    1998-01-01    2019-06-24
    223                   Unknown or Not Reported    1998-01-01    2019-06-24
    224                   Unknown or Not Reported    1998-01-01    2019-06-24
    225                   Unknown or Not Reported    2000-01-01    2019-06-24
    226                   Unknown or Not Reported    2000-01-01    2019-06-24
    227                   Unknown or Not Reported    2000-01-01    2019-06-24
    228                   Unknown or Not Reported    2000-01-01    2019-06-24
    229                   Unknown or Not Reported    2000-01-01    2019-06-24
    230                   Unknown or Not Reported    2000-01-01    2019-06-24
    231                   Unknown or Not Reported    2000-01-01    2019-06-24
             dataset        age  boost_age specimen_id actual_day_relative_to_boost
    1   2021_dataset 14221 days 11785 days         468                           -4
    2   2021_dataset 14221 days 11785 days         469                            1
    3   2021_dataset 14221 days 11785 days         470                            3
    4   2021_dataset 14221 days 11785 days         471                            7
    5   2021_dataset 14221 days 11785 days         472                           14
    6   2021_dataset 14221 days 11785 days         473                           30
    7   2021_dataset 14221 days 11785 days         474                           91
    8   2021_dataset 12029 days  9460 days         475                            0
    9   2021_dataset 12029 days  9460 days         476                            1
    10  2021_dataset 12029 days  9460 days         477                            3
    11  2021_dataset 12029 days  9460 days         478                            7
    12  2021_dataset 12029 days  9460 days         479                           14
    13  2021_dataset 12029 days  9460 days         480                           30
    14  2021_dataset 12029 days  9460 days         481                          101
    15  2021_dataset 11299 days  8730 days         483                            0
    16  2021_dataset 11299 days  8730 days         484                            1
    17  2021_dataset 11299 days  8730 days         485                            3
    18  2021_dataset 11299 days  8730 days         486                            7
    19  2021_dataset 11299 days  8730 days         487                           14
    20  2021_dataset 11299 days  8730 days         488                           38
    21  2021_dataset 11299 days  8730 days         489                          121
    22  2021_dataset 12029 days  9460 days         490                            0
    23  2021_dataset 12029 days  9460 days         491                            1
    24  2021_dataset 12029 days  9460 days         492                            3
    25  2021_dataset 12029 days  9460 days         493                            7
    26  2021_dataset 12029 days  9460 days         494                           14
    27  2021_dataset 12029 days  9460 days         495                           30
    28  2021_dataset 12029 days  9460 days         496                          101
    29  2021_dataset 13125 days 10563 days         498                            0
    30  2021_dataset 13125 days 10563 days         499                            1
    31  2021_dataset 13125 days 10563 days         500                            3
    32  2021_dataset 13125 days 10563 days         501                            7
    33  2021_dataset 13125 days 10563 days         502                           14
    34  2021_dataset 13125 days 10563 days         503                           37
    35  2021_dataset 13125 days 10563 days         504                           98
    36  2021_dataset 18239 days 15677 days         506                            0
    37  2021_dataset 18239 days 15677 days         507                            1
    38  2021_dataset 18239 days 15677 days         508                            3
    39  2021_dataset 18239 days 15677 days         509                            7
    40  2021_dataset 18239 days 15677 days         510                           14
    41  2021_dataset 18239 days 15677 days         511                           31
    42  2021_dataset 18239 days 15677 days         512                          101
    43  2021_dataset 19700 days 17194 days         513                            0
    44  2021_dataset 19700 days 17194 days         514                            1
    45  2021_dataset 19700 days 17194 days         515                            3
    46  2021_dataset 19700 days 17194 days         516                            7
    47  2021_dataset 19700 days 17194 days         517                           14
    48  2021_dataset 19700 days 17194 days         518                           30
    49  2021_dataset 19700 days 17194 days         519                           93
    50  2021_dataset 19700 days 17194 days         521                            0
    51  2021_dataset 19700 days 17194 days         522                            1
    52  2021_dataset 19700 days 17194 days         523                            3
    53  2021_dataset 19700 days 17194 days         524                            7
    54  2021_dataset 19700 days 17194 days         525                           14
    55  2021_dataset 19700 days 17194 days         526                           30
    56  2021_dataset 19700 days 17194 days         527                           93
    57  2021_dataset 13125 days 10619 days         529                            0
    58  2021_dataset 13125 days 10619 days         530                            1
    59  2021_dataset 13125 days 10619 days         531                            3
    60  2021_dataset 13125 days 10619 days         532                            7
    61  2021_dataset 13125 days 10619 days         533                           14
    62  2021_dataset 13125 days 10619 days         534                           32
    63  2021_dataset 13125 days 10619 days         535                           91
    64  2021_dataset 10203 days  7697 days         537                            0
    65  2021_dataset 10203 days  7697 days         538                            1
    66  2021_dataset 10203 days  7697 days         539                            3
    67  2021_dataset 10203 days  7697 days         540                            7
    68  2021_dataset 10203 days  7697 days         541                           14
    69  2021_dataset 10203 days  7697 days         542                           32
    70  2021_dataset 10203 days  7697 days         543                           93
    71  2021_dataset 10203 days  7697 days         546                            0
    72  2021_dataset 10203 days  7697 days         547                            1
    73  2021_dataset 10203 days  7697 days         548                            3
    74  2021_dataset 10203 days  7697 days         549                            7
    75  2021_dataset 10203 days  7697 days         550                           14
    76  2021_dataset 10203 days  7697 days         551                           37
    77  2021_dataset 10203 days  7697 days         552                          108
    78  2021_dataset 12760 days 10282 days         554                            0
    79  2021_dataset 12760 days 10282 days         555                            1
    80  2021_dataset 12760 days 10282 days         556                            3
    81  2021_dataset 12760 days 10282 days         557                            7
    82  2021_dataset 12760 days 10282 days         558                           14
    83  2021_dataset 12760 days 10282 days         559                           29
    84  2021_dataset 12760 days 10282 days         560                           94
    85  2021_dataset 11299 days  8821 days         562                            0
    86  2021_dataset 11299 days  8821 days         563                            1
    87  2021_dataset 11299 days  8821 days         564                            3
    88  2021_dataset 11299 days  8821 days         565                            7
    89  2021_dataset 11299 days  8821 days         566                           14
    90  2021_dataset 11299 days  8821 days         567                           37
    91  2021_dataset 11299 days  8821 days         568                           98
    92  2021_dataset 11299 days  8821 days         569                            0
    93  2021_dataset 11299 days  8821 days         570                            1
    94  2021_dataset 11299 days  8821 days         571                            3
    95  2021_dataset 11299 days  8821 days         572                            7
    96  2021_dataset 11299 days  8821 days         573                           14
    97  2021_dataset 11299 days  8821 days         574                           29
    98  2021_dataset 11299 days  8821 days         575                           94
    99  2021_dataset 10203 days  7725 days         577                            0
    100 2021_dataset 10203 days  7725 days         578                            1
    101 2021_dataset 10203 days  7725 days         579                            3
    102 2021_dataset 10203 days  7725 days         580                            7
    103 2021_dataset 10203 days  7725 days         581                           14
    104 2021_dataset 10203 days  7725 days         582                           29
    105 2021_dataset 10203 days  7725 days         583                           94
    106 2021_dataset 10203 days  7725 days         585                            0
    107 2021_dataset 10203 days  7725 days         586                            1
    108 2021_dataset 10203 days  7725 days         587                            3
    109 2021_dataset 10203 days  7725 days         588                            7
    110 2021_dataset 10203 days  7725 days         589                           14
    111 2021_dataset 10203 days  7725 days         590                           30
    112 2021_dataset 10203 days  7725 days         591                           93
    113 2021_dataset 13856 days 11399 days         593                            0
    114 2021_dataset 13856 days 11399 days         594                            1
    115 2021_dataset 13856 days 11399 days         595                            3
    116 2021_dataset 13856 days 11399 days         596                            7
    117 2021_dataset 13856 days 11399 days         597                           14
    118 2021_dataset 13856 days 11399 days         598                           31
    119 2021_dataset 13856 days 11399 days         599                           94
    120 2021_dataset 12029 days  9572 days         601                            0
    121 2021_dataset 12029 days  9572 days         602                            1
    122 2021_dataset 12029 days  9572 days         603                            3
    123 2021_dataset 12029 days  9572 days         604                            7
    124 2021_dataset 12029 days  9572 days         605                           14
    125 2021_dataset 12029 days  9572 days         606                           31
    126 2021_dataset 12029 days  9572 days         607                           92
    127 2021_dataset 14221 days 11764 days         608                            0
    128 2021_dataset 14221 days 11764 days         609                            1
    129 2021_dataset 14221 days 11764 days         610                            3
    130 2021_dataset 14221 days 11764 days         611                            7
    131 2021_dataset 14221 days 11764 days         612                           14
    132 2021_dataset 14221 days 11764 days         613                           31
    133 2021_dataset 14221 days 11764 days         614                           92
    134 2021_dataset 12395 days  9938 days         616                            0
    135 2021_dataset 12395 days  9938 days         617                            1
    136 2021_dataset 12395 days  9938 days         618                            3
    137 2021_dataset 12395 days  9938 days         619                            7
    138 2021_dataset 12395 days  9938 days         620                           14
    139 2021_dataset 12395 days  9938 days         621                           31
    140 2021_dataset 12395 days  9938 days         622                           92
    141 2021_dataset 12029 days  9572 days         623                            0
    142 2021_dataset 12029 days  9572 days         624                            1
    143 2021_dataset 12029 days  9572 days         625                            3
    144 2021_dataset 12029 days  9572 days         626                            7
    145 2021_dataset 12029 days  9572 days         627                           14
    146 2021_dataset 12029 days  9572 days         628                           36
    147 2021_dataset 12029 days  9572 days         629                          100
    148 2021_dataset  9838 days  7402 days         636                            0
    149 2021_dataset  9838 days  7402 days         637                            1
    150 2021_dataset  9838 days  7402 days         638                            3
    151 2021_dataset  9838 days  7402 days         639                            7
    152 2021_dataset  9838 days  7402 days         640                           14
    153 2021_dataset  9838 days  7402 days         641                           30
    154 2021_dataset  9838 days  7402 days         642                           99
    155 2021_dataset 10568 days  8132 days         643                            0
    156 2021_dataset 10568 days  8132 days         644                            1
    157 2021_dataset 10568 days  8132 days         645                            3
    158 2021_dataset 10568 days  8132 days         646                            7
    159 2021_dataset 10568 days  8132 days         647                           15
    160 2021_dataset 10568 days  8132 days         648                           28
    161 2021_dataset 10568 days  8132 days         649                           95
    162 2021_dataset  9473 days  7058 days         650                            0
    163 2021_dataset  9473 days  7058 days         651                            1
    164 2021_dataset  9473 days  7058 days         652                            3
    165 2021_dataset  9473 days  7058 days         653                            7
    166 2021_dataset  9473 days  7058 days         654                           14
    167 2021_dataset  9473 days  7058 days         655                           30
    168 2021_dataset  9473 days  7058 days         656                          150
    169 2021_dataset 10203 days  7788 days         657                            0
    170 2021_dataset 10203 days  7788 days         658                            1
    171 2021_dataset 10203 days  7788 days         659                            3
    172 2021_dataset 10203 days  7788 days         660                            7
    173 2021_dataset 10203 days  7788 days         661                           14
    174 2021_dataset 10203 days  7788 days         662                           30
    175 2021_dataset 10203 days  7788 days         663                          102
    176 2021_dataset 10568 days  8188 days         674                            0
    177 2021_dataset 10568 days  8188 days         675                            1
    178 2021_dataset 10568 days  8188 days         676                            3
    179 2021_dataset 10568 days  8188 days         677                            7
    180 2021_dataset 10568 days  8188 days         678                           14
    181 2021_dataset 10568 days  8188 days         679                           29
    182 2021_dataset 10568 days  8188 days         680                          112
    183 2021_dataset  9838 days  7458 days         681                            0
    184 2021_dataset  9838 days  7458 days         682                            1
    185 2021_dataset  9838 days  7458 days         683                            3
    186 2021_dataset  9838 days  7458 days         684                            7
    187 2021_dataset  9838 days  7458 days         685                           14
    188 2021_dataset  9838 days  7458 days         686                           29
    189 2021_dataset  9838 days  7458 days         687                          112
    190 2021_dataset 10203 days  7823 days         688                            0
    191 2021_dataset 10203 days  7823 days         689                            1
    192 2021_dataset 10203 days  7823 days         690                            3
    193 2021_dataset 10203 days  7823 days         691                            7
    194 2021_dataset 10203 days  7823 days         692                           14
    195 2021_dataset 10203 days  7823 days         693                           37
    196 2021_dataset 10203 days  7823 days         694                          107
    197 2021_dataset  9473 days  7114 days         695                            0
    198 2021_dataset  9473 days  7114 days         696                            1
    199 2021_dataset  9473 days  7114 days         697                            3
    200 2021_dataset  9473 days  7114 days         698                            7
    201 2021_dataset  9473 days  7114 days         699                           14
    202 2021_dataset  9473 days  7114 days         700                           30
    203 2021_dataset  9473 days  7114 days         701                           93
    204 2021_dataset 10934 days  8575 days         702                            0
    205 2021_dataset 10934 days  8575 days         703                            1
    206 2021_dataset 10934 days  8575 days         704                            3
    207 2021_dataset 10934 days  8575 days         705                            7
    208 2021_dataset 10934 days  8575 days         706                           14
    209 2021_dataset 10934 days  8575 days         707                           30
    210 2021_dataset 10934 days  8575 days         708                           95
    211 2021_dataset  9838 days  7479 days         709                            0
    212 2021_dataset  9838 days  7479 days         710                            1
    213 2021_dataset  9838 days  7479 days         711                            3
    214 2021_dataset  9838 days  7479 days         712                            7
    215 2021_dataset  9838 days  7479 days         713                           14
    216 2021_dataset  9838 days  7479 days         714                           30
    217 2021_dataset  9838 days  7479 days         715                           93
    218 2021_dataset 10203 days  7844 days         716                            0
    219 2021_dataset 10203 days  7844 days         717                            1
    220 2021_dataset 10203 days  7844 days         718                            3
    221 2021_dataset 10203 days  7844 days         719                            7
    222 2021_dataset 10203 days  7844 days         720                           14
    223 2021_dataset 10203 days  7844 days         721                           30
    224 2021_dataset 10203 days  7844 days         722                           93
    225 2021_dataset  9473 days  7114 days         723                            0
    226 2021_dataset  9473 days  7114 days         724                            1
    227 2021_dataset  9473 days  7114 days         725                            3
    228 2021_dataset  9473 days  7114 days         726                            7
    229 2021_dataset  9473 days  7114 days         727                           14
    230 2021_dataset  9473 days  7114 days         728                           30
    231 2021_dataset  9473 days  7114 days         729                           93
        planned_day_relative_to_boost specimen_type visit isotype
    1                               0         Blood     1     IgG
    2                               1         Blood     2     IgG
    3                               3         Blood     3     IgG
    4                               7         Blood     4     IgG
    5                              14         Blood     5     IgG
    6                              30         Blood     6     IgG
    7                             120         Blood     7     IgG
    8                               0         Blood     1     IgG
    9                               1         Blood     2     IgG
    10                              3         Blood     3     IgG
    11                              7         Blood     4     IgG
    12                             14         Blood     5     IgG
    13                             30         Blood     6     IgG
    14                            120         Blood     7     IgG
    15                              0         Blood     1     IgG
    16                              1         Blood     2     IgG
    17                              3         Blood     3     IgG
    18                              7         Blood     4     IgG
    19                             14         Blood     5     IgG
    20                             30         Blood     6     IgG
    21                            120         Blood     7     IgG
    22                              0         Blood     1     IgG
    23                              1         Blood     2     IgG
    24                              3         Blood     3     IgG
    25                              7         Blood     4     IgG
    26                             14         Blood     5     IgG
    27                             30         Blood     6     IgG
    28                            120         Blood     7     IgG
    29                              0         Blood     1     IgG
    30                              1         Blood     2     IgG
    31                              3         Blood     3     IgG
    32                              7         Blood     4     IgG
    33                             14         Blood     5     IgG
    34                             30         Blood     6     IgG
    35                            120         Blood     7     IgG
    36                              0         Blood     1     IgG
    37                              1         Blood     2     IgG
    38                              3         Blood     3     IgG
    39                              7         Blood     4     IgG
    40                             14         Blood     5     IgG
    41                             30         Blood     6     IgG
    42                            120         Blood     7     IgG
    43                              0         Blood     1     IgG
    44                              1         Blood     2     IgG
    45                              3         Blood     3     IgG
    46                              7         Blood     4     IgG
    47                             14         Blood     5     IgG
    48                             30         Blood     6     IgG
    49                            120         Blood     7     IgG
    50                              0         Blood     1     IgG
    51                              1         Blood     2     IgG
    52                              3         Blood     3     IgG
    53                              7         Blood     4     IgG
    54                             14         Blood     5     IgG
    55                             30         Blood     6     IgG
    56                            120         Blood     7     IgG
    57                              0         Blood     1     IgG
    58                              1         Blood     2     IgG
    59                              3         Blood     3     IgG
    60                              7         Blood     4     IgG
    61                             14         Blood     5     IgG
    62                             30         Blood     6     IgG
    63                            120         Blood     7     IgG
    64                              0         Blood     1     IgG
    65                              1         Blood     2     IgG
    66                              3         Blood     3     IgG
    67                              7         Blood     4     IgG
    68                             14         Blood     5     IgG
    69                             30         Blood     6     IgG
    70                            120         Blood     7     IgG
    71                              0         Blood     1     IgG
    72                              1         Blood     2     IgG
    73                              3         Blood     3     IgG
    74                              7         Blood     4     IgG
    75                             14         Blood     5     IgG
    76                             30         Blood     6     IgG
    77                            120         Blood     7     IgG
    78                              0         Blood     1     IgG
    79                              1         Blood     2     IgG
    80                              3         Blood     3     IgG
    81                              7         Blood     4     IgG
    82                             14         Blood     5     IgG
    83                             30         Blood     6     IgG
    84                            120         Blood     7     IgG
    85                              0         Blood     1     IgG
    86                              1         Blood     2     IgG
    87                              3         Blood     3     IgG
    88                              7         Blood     4     IgG
    89                             14         Blood     5     IgG
    90                             30         Blood     6     IgG
    91                            120         Blood     7     IgG
    92                              0         Blood     1     IgG
    93                              1         Blood     2     IgG
    94                              3         Blood     3     IgG
    95                              7         Blood     4     IgG
    96                             14         Blood     5     IgG
    97                             30         Blood     6     IgG
    98                            120         Blood     7     IgG
    99                              0         Blood     1     IgG
    100                             1         Blood     2     IgG
    101                             3         Blood     3     IgG
    102                             7         Blood     4     IgG
    103                            14         Blood     5     IgG
    104                            30         Blood     6     IgG
    105                           120         Blood     7     IgG
    106                             0         Blood     1     IgG
    107                             1         Blood     2     IgG
    108                             3         Blood     3     IgG
    109                             7         Blood     4     IgG
    110                            14         Blood     5     IgG
    111                            30         Blood     6     IgG
    112                           120         Blood     7     IgG
    113                             0         Blood     1     IgG
    114                             1         Blood     2     IgG
    115                             3         Blood     3     IgG
    116                             7         Blood     4     IgG
    117                            14         Blood     5     IgG
    118                            30         Blood     6     IgG
    119                           120         Blood     7     IgG
    120                             0         Blood     1     IgG
    121                             1         Blood     2     IgG
    122                             3         Blood     3     IgG
    123                             7         Blood     4     IgG
    124                            14         Blood     5     IgG
    125                            30         Blood     6     IgG
    126                           120         Blood     7     IgG
    127                             0         Blood     1     IgG
    128                             1         Blood     2     IgG
    129                             3         Blood     3     IgG
    130                             7         Blood     4     IgG
    131                            14         Blood     5     IgG
    132                            30         Blood     6     IgG
    133                           120         Blood     7     IgG
    134                             0         Blood     1     IgG
    135                             1         Blood     2     IgG
    136                             3         Blood     3     IgG
    137                             7         Blood     4     IgG
    138                            14         Blood     5     IgG
    139                            30         Blood     6     IgG
    140                           120         Blood     7     IgG
    141                             0         Blood     1     IgG
    142                             1         Blood     2     IgG
    143                             3         Blood     3     IgG
    144                             7         Blood     4     IgG
    145                            14         Blood     5     IgG
    146                            30         Blood     6     IgG
    147                           120         Blood     7     IgG
    148                             0         Blood     1     IgG
    149                             1         Blood     2     IgG
    150                             3         Blood     3     IgG
    151                             7         Blood     4     IgG
    152                            14         Blood     5     IgG
    153                            30         Blood     6     IgG
    154                           120         Blood     7     IgG
    155                             0         Blood     1     IgG
    156                             1         Blood     2     IgG
    157                             3         Blood     3     IgG
    158                             7         Blood     4     IgG
    159                            14         Blood     5     IgG
    160                            30         Blood     6     IgG
    161                           120         Blood     7     IgG
    162                             0         Blood     1     IgG
    163                             1         Blood     2     IgG
    164                             3         Blood     3     IgG
    165                             7         Blood     4     IgG
    166                            14         Blood     5     IgG
    167                            30         Blood     6     IgG
    168                           120         Blood     7     IgG
    169                             0         Blood     1     IgG
    170                             1         Blood     2     IgG
    171                             3         Blood     3     IgG
    172                             7         Blood     4     IgG
    173                            14         Blood     5     IgG
    174                            30         Blood     6     IgG
    175                           120         Blood     7     IgG
    176                             0         Blood     1     IgG
    177                             1         Blood     2     IgG
    178                             3         Blood     3     IgG
    179                             7         Blood     4     IgG
    180                            14         Blood     5     IgG
    181                            30         Blood     6     IgG
    182                           120         Blood     7     IgG
    183                             0         Blood     1     IgG
    184                             1         Blood     2     IgG
    185                             3         Blood     3     IgG
    186                             7         Blood     4     IgG
    187                            14         Blood     5     IgG
    188                            30         Blood     6     IgG
    189                           120         Blood     7     IgG
    190                             0         Blood     1     IgG
    191                             1         Blood     2     IgG
    192                             3         Blood     3     IgG
    193                             7         Blood     4     IgG
    194                            14         Blood     5     IgG
    195                            30         Blood     6     IgG
    196                           120         Blood     7     IgG
    197                             0         Blood     1     IgG
    198                             1         Blood     2     IgG
    199                             3         Blood     3     IgG
    200                             7         Blood     4     IgG
    201                            14         Blood     5     IgG
    202                            30         Blood     6     IgG
    203                           120         Blood     7     IgG
    204                             0         Blood     1     IgG
    205                             1         Blood     2     IgG
    206                             3         Blood     3     IgG
    207                             7         Blood     4     IgG
    208                            14         Blood     5     IgG
    209                            30         Blood     6     IgG
    210                           120         Blood     7     IgG
    211                             0         Blood     1     IgG
    212                             1         Blood     2     IgG
    213                             3         Blood     3     IgG
    214                             7         Blood     4     IgG
    215                            14         Blood     5     IgG
    216                            30         Blood     6     IgG
    217                           120         Blood     7     IgG
    218                             0         Blood     1     IgG
    219                             1         Blood     2     IgG
    220                             3         Blood     3     IgG
    221                             7         Blood     4     IgG
    222                            14         Blood     5     IgG
    223                            30         Blood     6     IgG
    224                           120         Blood     7     IgG
    225                             0         Blood     1     IgG
    226                             1         Blood     2     IgG
    227                             3         Blood     3     IgG
    228                             7         Blood     4     IgG
    229                            14         Blood     5     IgG
    230                            30         Blood     6     IgG
    231                           120         Blood     7     IgG
        is_antigen_specific antigen         MFI MFI_normalised unit
    1                 FALSE      PT  112.750000     1.00000000  MFI
    2                 FALSE      PT  111.250000     0.98669623  MFI
    3                 FALSE      PT  125.500000     1.11308204  MFI
    4                 FALSE      PT  224.250000     1.98891353  MFI
    5                 FALSE      PT  304.000000     2.69623060  MFI
    6                 FALSE      PT  274.000000     2.43015521  MFI
    7                 FALSE      PT  171.750000     1.52328160  MFI
    8                 FALSE      PT  124.951277     1.10821532  MFI
    9                 FALSE      PT  275.201277     2.44080956  MFI
    10                FALSE      PT  263.201277     2.33437940  MFI
    11                FALSE      PT  852.201277     7.55832619  MFI
    12                FALSE      PT 1548.451277    13.73349248  MFI
    13                FALSE      PT 1194.451277    10.59380290  MFI
    14                FALSE      PT  679.951277     6.03061000  MFI
    15                FALSE      PT  278.701277     2.47185168  MFI
    16                FALSE      PT  441.451277     3.91531066  MFI
    17                FALSE      PT  498.451277     4.42085390  MFI
    18                FALSE      PT 1526.951277    13.54280512  MFI
    19                FALSE      PT 1916.701277    16.99956787  MFI
    20                FALSE      PT 1150.701277    10.20577630  MFI
    21                FALSE      PT  748.951277     6.64258339  MFI
    22                FALSE      PT  353.201277     3.13260556  MFI
    23                FALSE      PT  506.701277     4.49402463  MFI
    24                FALSE      PT  476.701277     4.22794925  MFI
    25                FALSE      PT  443.201277     3.93083173  MFI
    26                FALSE      PT  558.451277     4.95300468  MFI
    27                FALSE      PT  458.451277     4.06608672  MFI
    28                FALSE      PT  408.451277     3.62262774  MFI
    29                FALSE      PT   83.750000     0.74279379  MFI
    30                FALSE      PT   83.750000     0.74279379  MFI
    31                FALSE      PT  102.000000     0.90465632  MFI
    32                FALSE      PT  504.500000     4.47450111  MFI
    33                FALSE      PT 1837.750000    16.29933481  MFI
    34                FALSE      PT 1243.250000    11.02660754  MFI
    35                FALSE      PT  780.750000     6.92461197  MFI
    36                FALSE      PT  169.750000     1.50554324  MFI
    37                FALSE      PT  147.750000     1.31042129  MFI
    38                FALSE      PT  149.750000     1.32815965  MFI
    39                FALSE      PT 1849.000000    16.39911308  MFI
    40                FALSE      PT 2059.500000    18.26607539  MFI
    41                FALSE      PT 1194.750000    10.59645233  MFI
    42                FALSE      PT  417.500000     3.70288248  MFI
    43                FALSE      PT   71.200062     0.63148614  MFI
    44                FALSE      PT   66.450062     0.58935754  MFI
    45                FALSE      PT   64.450062     0.57161918  MFI
    46                FALSE      PT  369.450062     3.27671896  MFI
    47                FALSE      PT 2128.450062    18.87760587  MFI
    48                FALSE      PT 1606.200062    14.24567683  MFI
    49                FALSE      PT  821.450062     7.28558814  MFI
    50                FALSE      PT  121.447441     1.07713917  MFI
    51                FALSE      PT  120.197441     1.06605269  MFI
    52                FALSE      PT  124.197441     1.10152941  MFI
    53                FALSE      PT  939.697441     8.33434538  MFI
    54                FALSE      PT 1565.447441    13.88423451  MFI
    55                FALSE      PT 1086.947441     9.64033207  MFI
    56                FALSE      PT  535.947441     4.75341411  MFI
    57                FALSE      PT  124.450062     1.10376995  MFI
    58                FALSE      PT  125.950062     1.11707372  MFI
    59                FALSE      PT  121.950062     1.08159701  MFI
    60                FALSE      PT  719.450062     6.38093182  MFI
    61                FALSE      PT 1150.200062    10.20133093  MFI
    62                FALSE      PT  782.700062     6.94190743  MFI
    63                FALSE      PT  395.200062     3.50510033  MFI
    64                FALSE      PT   43.451277     0.38537718  MFI
    65                FALSE      PT   64.701277     0.57384725  MFI
    66                FALSE      PT   56.701277     0.50289381  MFI
    67                FALSE      PT   50.701277     0.44967874  MFI
    68                FALSE      PT  201.201277     1.78449027  MFI
    69                FALSE      PT  121.451277     1.07717319  MFI
    70                FALSE      PT  156.701277     1.38981177  MFI
    71                FALSE      PT  732.701277     6.49845922  MFI
    72                FALSE      PT  959.451277     8.50954570  MFI
    73                FALSE      PT  740.201277     6.56497807  MFI
    74                FALSE      PT  723.201277     6.41420202  MFI
    75                FALSE      PT  916.951277     8.13260556  MFI
    76                FALSE      PT  637.201277     5.65145257  MFI
    77                FALSE      PT  631.201277     5.59823749  MFI
    78                FALSE      PT  366.750000     3.25277162  MFI
    79                FALSE      PT  285.500000     2.53215078  MFI
    80                FALSE      PT  359.500000     3.18847007  MFI
    81                FALSE      PT  963.750000     8.54767184  MFI
    82                FALSE      PT 1736.250000    15.39911308  MFI
    83                FALSE      PT 1120.750000     9.94013304  MFI
    84                FALSE      PT  610.000000     5.41019956  MFI
    85                FALSE      PT  303.197441     2.68911256  MFI
    86                FALSE      PT  325.197441     2.88423451  MFI
    87                FALSE      PT  344.947441     3.05940081  MFI
    88                FALSE      PT  981.947441     8.70906821  MFI
    89                FALSE      PT 1349.697441    11.97070901  MFI
    90                FALSE      PT  857.947441     7.60928994  MFI
    91                FALSE      PT  640.947441     5.68467797  MFI
    92                FALSE      PT   48.947441     0.43412365  MFI
    93                FALSE      PT   53.947441     0.47846954  MFI
    94                FALSE      PT   58.947441     0.52281544  MFI
    95                FALSE      PT  250.947441     2.22569793  MFI
    96                FALSE      PT  385.697441     3.42081988  MFI
    97                FALSE      PT  227.947441     2.02170679  MFI
    98                FALSE      PT  141.447441     1.25452276  MFI
    99                FALSE      PT   14.750000     0.13082040  MFI
    100               FALSE      PT    9.750000     0.08647450  MFI
    101               FALSE      PT   11.750000     0.10421286  MFI
    102               FALSE      PT   52.250000     0.46341463  MFI
    103               FALSE      PT  136.250000     1.20842572  MFI
    104               FALSE      PT  120.500000     1.06873614  MFI
    105               FALSE      PT   65.750000     0.58314856  MFI
    106               FALSE      PT 1195.500000    10.60310421  MFI
    107               FALSE      PT 1105.750000     9.80709534  MFI
    108               FALSE      PT 1124.250000     9.97117517  MFI
    109               FALSE      PT 1129.750000    10.01995565  MFI
    110               FALSE      PT 1247.000000    11.05986696  MFI
    111               FALSE      PT 1028.750000     9.12416851  MFI
    112               FALSE      PT  582.750000     5.16851441  MFI
    113               FALSE      PT   69.250000     0.61419069  MFI
    114               FALSE      PT   91.250000     0.80931264  MFI
    115               FALSE      PT   62.250000     0.55210643  MFI
    116               FALSE      PT  224.000000     1.98669623  MFI
    117               FALSE      PT  506.500000     4.49223947  MFI
    118               FALSE      PT  307.250000     2.72505543  MFI
    119               FALSE      PT  248.750000     2.20620843  MFI
    120               FALSE      PT  112.000000     0.99334812  MFI
    121               FALSE      PT  122.000000     1.08203991  MFI
    122               FALSE      PT  106.250000     0.94235033  MFI
    123               FALSE      PT  495.000000     4.39024390  MFI
    124               FALSE      PT  831.500000     7.37472284  MFI
    125               FALSE      PT  550.500000     4.88248337  MFI
    126               FALSE      PT  371.000000     3.29046563  MFI
    127               FALSE      PT  451.250000     4.00221729  MFI
    128               FALSE      PT  390.750000     3.46563193  MFI
    129               FALSE      PT  446.250000     3.95787140  MFI
    130               FALSE      PT  489.750000     4.34368071  MFI
    131               FALSE      PT 1019.250000     9.03991131  MFI
    132               FALSE      PT 1023.250000     9.07538803  MFI
    133               FALSE      PT  874.250000     7.75388027  MFI
    134               FALSE      PT   43.250000     0.38359202  MFI
    135               FALSE      PT   46.750000     0.41463415  MFI
    136               FALSE      PT   41.250000     0.36585366  MFI
    137               FALSE      PT  133.500000     1.18403548  MFI
    138               FALSE      PT  248.500000     2.20399113  MFI
    139               FALSE      PT  291.500000     2.58536585  MFI
    140               FALSE      PT  208.250000     1.84700665  MFI
    141               FALSE      PT  284.250000     2.52106430  MFI
    142               FALSE      PT  291.250000     2.58314856  MFI
    143               FALSE      PT  299.500000     2.65631929  MFI
    144               FALSE      PT  331.500000     2.94013304  MFI
    145               FALSE      PT  516.750000     4.58314856  MFI
    146               FALSE      PT  477.750000     4.23725055  MFI
    147               FALSE      PT  414.000000     3.67184035  MFI
    148               FALSE      PT  157.950062     1.40088747  MFI
    149               FALSE      PT  147.200062     1.30554379  MFI
    150               FALSE      PT  152.700062     1.35432428  MFI
    151               FALSE      PT  331.200062     2.93747284  MFI
    152               FALSE      PT  436.450062     3.87095399  MFI
    153               FALSE      PT  381.950062     3.38758370  MFI
    154               FALSE      PT  213.200062     1.89090964  MFI
    155               FALSE      PT   24.201277     0.21464548  MFI
    156               FALSE      PT   59.701277     0.52950135  MFI
    157               FALSE      PT   56.701277     0.50289381  MFI
    158               FALSE      PT  152.701277     1.35433505  MFI
    159               FALSE      PT  474.951277     4.21242818  MFI
    160               FALSE      PT  289.201277     2.56497807  MFI
    161               FALSE      PT  217.701277     1.93083173  MFI
    162               FALSE      PT   29.250000     0.25942350  MFI
    163               FALSE      PT   32.750000     0.29046563  MFI
    164               FALSE      PT   32.250000     0.28603104  MFI
    165               FALSE      PT  107.750000     0.95565410  MFI
    166               FALSE      PT  168.750000     1.49667406  MFI
    167               FALSE      PT  492.750000     4.37028825  MFI
    168               FALSE      PT   42.250000     0.37472284  MFI
    169               FALSE      PT   26.450062     0.23459035  MFI
    170               FALSE      PT   23.450062     0.20798281  MFI
    171               FALSE      PT  373.700062     3.31441297  MFI
    172               FALSE      PT  471.700062     4.18359257  MFI
    173               FALSE      PT  615.200062     5.45631984  MFI
    174               FALSE      PT  337.200062     2.99068791  MFI
    175               FALSE      PT  116.200062     1.03059922  MFI
    176               FALSE      PT    5.447441     0.04831433  MFI
    177               FALSE      PT    6.447441     0.05718351  MFI
    178               FALSE      PT    5.197441     0.04609704  MFI
    179               FALSE      PT   67.197441     0.59598617  MFI
    180               FALSE      PT  133.947441     1.18800391  MFI
    181               FALSE      PT   82.947441     0.73567575  MFI
    182               FALSE      PT   50.947441     0.45186201  MFI
    183               FALSE      PT   29.950062     0.26563248  MFI
    184               FALSE      PT   27.950062     0.24789412  MFI
    185               FALSE      PT   30.450062     0.27006707  MFI
    186               FALSE      PT   76.950062     0.68248392  MFI
    187               FALSE      PT  221.950062     1.96851497  MFI
    188               FALSE      PT  193.700062     1.71796064  MFI
    189               FALSE      PT   67.950062     0.60266131  MFI
    190               FALSE      PT  101.697441     0.90197287  MFI
    191               FALSE      PT   83.947441     0.74454493  MFI
    192               FALSE      PT  106.697441     0.94631877  MFI
    193               FALSE      PT  270.447441     2.39864693  MFI
    194               FALSE      PT  487.447441     4.32325890  MFI
    195               FALSE      PT  385.947441     3.42303717  MFI
    196               FALSE      PT  252.697441     2.24121899  MFI
    197               FALSE      PT   20.250000     0.17960089  MFI
    198               FALSE      PT   22.250000     0.19733925  MFI
    199               FALSE      PT   25.250000     0.22394678  MFI
    200               FALSE      PT  180.250000     1.59866962  MFI
    201               FALSE      PT  377.750000     3.35033259  MFI
    202               FALSE      PT  323.500000     2.86917960  MFI
    203               FALSE      PT  182.250000     1.61640798  MFI
    204               FALSE      PT  223.750000     1.98447894  MFI
    205               FALSE      PT  215.750000     1.91352550  MFI
    206               FALSE      PT  220.750000     1.95787140  MFI
    207               FALSE      PT  586.250000     5.19955654  MFI
    208               FALSE      PT  645.250000     5.72283814  MFI
    209               FALSE      PT  596.250000     5.28824834  MFI
    210               FALSE      PT  396.250000     3.51441242  MFI
    211               FALSE      PT  599.000000     5.31263858  MFI
    212               FALSE      PT  557.250000     4.94235033  MFI
    213               FALSE      PT  519.250000     4.60532151  MFI
    214               FALSE      PT  709.750000     6.29490022  MFI
    215               FALSE      PT  993.750000     8.81374723  MFI
    216               FALSE      PT  877.000000     7.77827051  MFI
    217               FALSE      PT  639.250000     5.66962306  MFI
    218               FALSE      PT  441.697441     3.91749393  MFI
    219               FALSE      PT  494.197441     4.38312586  MFI
    220               FALSE      PT  484.947441     4.30108595  MFI
    221               FALSE      PT  661.947441     5.87093074  MFI
    222               FALSE      PT 1078.697441     9.56716134  MFI
    223               FALSE      PT  856.697441     7.59820347  MFI
    224               FALSE      PT  651.947441     5.78223895  MFI
    225               FALSE      PT   78.750000     0.69844789  MFI
    226               FALSE      PT   87.750000     0.77827051  MFI
    227               FALSE      PT   78.750000     0.69844789  MFI
    228               FALSE      PT  140.250000     1.24390244  MFI
    229               FALSE      PT  386.000000     3.42350333  MFI
    230               FALSE      PT  331.000000     2.93569845  MFI
    231               FALSE      PT  304.750000     2.70288248  MFI
        lower_limit_of_detection
    1                   5.197441
    2                   5.197441
    3                   5.197441
    4                   5.197441
    5                   5.197441
    6                   5.197441
    7                   5.197441
    8                   5.197441
    9                   5.197441
    10                  5.197441
    11                  5.197441
    12                  5.197441
    13                  5.197441
    14                  5.197441
    15                  5.197441
    16                  5.197441
    17                  5.197441
    18                  5.197441
    19                  5.197441
    20                  5.197441
    21                  5.197441
    22                  5.197441
    23                  5.197441
    24                  5.197441
    25                  5.197441
    26                  5.197441
    27                  5.197441
    28                  5.197441
    29                  5.197441
    30                  5.197441
    31                  5.197441
    32                  5.197441
    33                  5.197441
    34                  5.197441
    35                  5.197441
    36                  5.197441
    37                  5.197441
    38                  5.197441
    39                  5.197441
    40                  5.197441
    41                  5.197441
    42                  5.197441
    43                  5.197441
    44                  5.197441
    45                  5.197441
    46                  5.197441
    47                  5.197441
    48                  5.197441
    49                  5.197441
    50                  5.197441
    51                  5.197441
    52                  5.197441
    53                  5.197441
    54                  5.197441
    55                  5.197441
    56                  5.197441
    57                  5.197441
    58                  5.197441
    59                  5.197441
    60                  5.197441
    61                  5.197441
    62                  5.197441
    63                  5.197441
    64                  5.197441
    65                  5.197441
    66                  5.197441
    67                  5.197441
    68                  5.197441
    69                  5.197441
    70                  5.197441
    71                  5.197441
    72                  5.197441
    73                  5.197441
    74                  5.197441
    75                  5.197441
    76                  5.197441
    77                  5.197441
    78                  5.197441
    79                  5.197441
    80                  5.197441
    81                  5.197441
    82                  5.197441
    83                  5.197441
    84                  5.197441
    85                  5.197441
    86                  5.197441
    87                  5.197441
    88                  5.197441
    89                  5.197441
    90                  5.197441
    91                  5.197441
    92                  5.197441
    93                  5.197441
    94                  5.197441
    95                  5.197441
    96                  5.197441
    97                  5.197441
    98                  5.197441
    99                  5.197441
    100                 5.197441
    101                 5.197441
    102                 5.197441
    103                 5.197441
    104                 5.197441
    105                 5.197441
    106                 5.197441
    107                 5.197441
    108                 5.197441
    109                 5.197441
    110                 5.197441
    111                 5.197441
    112                 5.197441
    113                 5.197441
    114                 5.197441
    115                 5.197441
    116                 5.197441
    117                 5.197441
    118                 5.197441
    119                 5.197441
    120                 5.197441
    121                 5.197441
    122                 5.197441
    123                 5.197441
    124                 5.197441
    125                 5.197441
    126                 5.197441
    127                 5.197441
    128                 5.197441
    129                 5.197441
    130                 5.197441
    131                 5.197441
    132                 5.197441
    133                 5.197441
    134                 5.197441
    135                 5.197441
    136                 5.197441
    137                 5.197441
    138                 5.197441
    139                 5.197441
    140                 5.197441
    141                 5.197441
    142                 5.197441
    143                 5.197441
    144                 5.197441
    145                 5.197441
    146                 5.197441
    147                 5.197441
    148                 5.197441
    149                 5.197441
    150                 5.197441
    151                 5.197441
    152                 5.197441
    153                 5.197441
    154                 5.197441
    155                 5.197441
    156                 5.197441
    157                 5.197441
    158                 5.197441
    159                 5.197441
    160                 5.197441
    161                 5.197441
    162                 5.197441
    163                 5.197441
    164                 5.197441
    165                 5.197441
    166                 5.197441
    167                 5.197441
    168                 5.197441
    169                 5.197441
    170                 5.197441
    171                 5.197441
    172                 5.197441
    173                 5.197441
    174                 5.197441
    175                 5.197441
    176                 5.197441
    177                 5.197441
    178                 5.197441
    179                 5.197441
    180                 5.197441
    181                 5.197441
    182                 5.197441
    183                 5.197441
    184                 5.197441
    185                 5.197441
    186                 5.197441
    187                 5.197441
    188                 5.197441
    189                 5.197441
    190                 5.197441
    191                 5.197441
    192                 5.197441
    193                 5.197441
    194                 5.197441
    195                 5.197441
    196                 5.197441
    197                 5.197441
    198                 5.197441
    199                 5.197441
    200                 5.197441
    201                 5.197441
    202                 5.197441
    203                 5.197441
    204                 5.197441
    205                 5.197441
    206                 5.197441
    207                 5.197441
    208                 5.197441
    209                 5.197441
    210                 5.197441
    211                 5.197441
    212                 5.197441
    213                 5.197441
    214                 5.197441
    215                 5.197441
    216                 5.197441
    217                 5.197441
    218                 5.197441
    219                 5.197441
    220                 5.197441
    221                 5.197441
    222                 5.197441
    223                 5.197441
    224                 5.197441
    225                 5.197441
    226                 5.197441
    227                 5.197441
    228                 5.197441
    229                 5.197441
    230                 5.197441
    231                 5.197441

``` r
pt_2020 |>
  ggplot() +
    aes(planned_day_relative_to_boost,
        MFI_normalised, 
        color=infancy_vac, 
        group=subject_id) +
    geom_point() +
    geom_line() +
    theme_classic() +
    geom_vline(xintercept=0, col="darkgreen") +
    geom_vline(xintercept =14, col="darkgreen") +
    labs(title="2020 dataset IgG PT",
       subtitle = "Dashed lines indicate day 0 (pre-boost) and 14 (apparent peak levels)",x="Days post booster", y="Normalized MFI for PT IgG")
```

![](class19_files/figure-commonmark/unnamed-chunk-25-1.png)

``` r
pt_2020 |>
  ggplot() +
    aes(planned_day_relative_to_boost,
        MFI_normalised, 
        color=infancy_vac, 
        group=subject_id) +
    geom_point() +
    geom_line() +
    theme_classic() +
    geom_vline(xintercept=0, col="darkgreen") +
    geom_vline(xintercept =14, col="darkgreen") +
    scale_x_continuous(limits = c(0, 125)) +
    labs(title="2020 dataset IgG PT",
       subtitle = "Dashed lines indicate day 0 (pre-boost) and 14 (apparent peak levels), \n scale cropped at 125 days for better comparison with 2021",x="Days post booster", y="Normalized MFI for PT IgG")
```

    Warning: Removed 3 rows containing missing values or values outside the scale range
    (`geom_point()`).

    Warning: Removed 3 rows containing missing values or values outside the scale range
    (`geom_line()`).

![](class19_files/figure-commonmark/unnamed-chunk-25-2.png)

``` r
pt_2021 |>
  ggplot() +
    aes(planned_day_relative_to_boost,
        MFI_normalised, 
        color=infancy_vac, 
        group=subject_id) +
    geom_point() +
    geom_line() +
    theme_classic() +
    geom_vline(xintercept=0, col="darkgreen") +
    geom_vline(xintercept =14, col="darkgreen") +
    labs(title="2021 dataset IgG PT",
       subtitle = "Dashed lines indicate day 0 (pre-boost) and 14 (apparent peak levels)",x="Days post booster", y="Normalized MFI for PT IgG")
```

![](class19_files/figure-commonmark/unnamed-chunk-25-3.png)

> Q18. Does this trend look similar for the 2020 dataset?

Hard to tell - the 2020 trend looks quite smushed bnecause of twice the
number of days being recorded. I would guess probably.

Actually, no, MFIs were HIGHER in 2020 at 14 days.

## System setup

``` r
sessionInfo()
```

    R version 4.5.1 (2025-06-13 ucrt)
    Platform: x86_64-w64-mingw32/x64
    Running under: Windows 11 x64 (build 26200)

    Matrix products: default
      LAPACK version 3.12.1

    locale:
    [1] LC_COLLATE=English_United States.utf8 
    [2] LC_CTYPE=English_United States.utf8   
    [3] LC_MONETARY=English_United States.utf8
    [4] LC_NUMERIC=C                          
    [5] LC_TIME=English_United States.utf8    

    time zone: America/Los_Angeles
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
    [1] dplyr_1.1.4     lubridate_1.9.4 jsonlite_2.0.0  ggplot2_4.0.1  

    loaded via a namespace (and not attached):
     [1] vctrs_0.6.5        cli_3.6.5          knitr_1.50         rlang_1.1.6       
     [5] xfun_0.54          generics_0.1.4     S7_0.2.1           labeling_0.4.3    
     [9] glue_1.8.0         htmltools_0.5.8.1  scales_1.4.0       rmarkdown_2.30    
    [13] grid_4.5.1         evaluate_1.0.5     tibble_3.3.0       fastmap_1.2.0     
    [17] yaml_2.3.10        lifecycle_1.0.4    compiler_4.5.1     RColorBrewer_1.1-3
    [21] timechange_0.3.0   pkgconfig_2.0.3    rstudioapi_0.17.1  farver_2.1.2      
    [25] digest_0.6.39      R6_2.6.1           tidyselect_1.2.1   pillar_1.11.1     
    [29] magrittr_2.0.4     withr_3.0.2        tools_4.5.1        gtable_0.3.6      
