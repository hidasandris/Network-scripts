# Network-scripts
Scripts for network science.

## TI, WI calculator for R
> Input data frame should be without header, the columns should be ordered: node name of species A, node name of species B, weight. Symmetrization_method is optional, can be either 'sum' (default), 'average', 'max', 'min' or 'difference.' Multilink edges will be summed.
>
> If input data does not contain weights (3rd column), then only TI will be calculated.

    calculate_TI_WI <- function(input_data, steps, symmetrization_method="sum") {
        ...
    }

Input data should be an edge list without header:

| node A | node B | weight |
| ------ | ------ | ------ |
| A      | B      | 1.12   |
| B      | C      | 0.12   |
| A      | D      | 2.23   |

### Use the script like this:
1. Import the function calculate_TI_WI.R.
2. Import the dataset to a data frame you want to use.
3. Calculate the TI and WI indeces with 3 steps.
```
source("calculate_TI_WI.R")
df = read.table("abcd.txt")
result = calculate_TI_WI(df, 3)
```
