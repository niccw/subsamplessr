# subsamplessr

#### Description

This is a simple function subsample a portion of ssr data and calculate the mean and sd of the subsampled set

#### Usage

```R
subsample_countmean(filepath, pct = 0.1, times = 20,
  outputname = "sample")
```

#### Arguments

`pct`	Percentage of data to subsample
`times`	Number of run
`outputname` Prefix of output files

#### Output

This function output 4 csv files:

1) full_locustable: A locustable of the full dataset

2) subsample_locustable: locustable of all subsample dataset

3) subsample_lis

4)subsample_locusallelmeansd: allel mean and sd of all subsample dataset



- Genind object of the whole dataset
- List object of all subsample data 
- Matrix descriping meand and sd of the subsample dataset