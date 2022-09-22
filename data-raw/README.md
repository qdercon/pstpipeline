# Using the raw data outside of pstpipeline

The task data is stored here as an R data object (```all_res_split.RDS```). The raw experiment outputs for each partcipant are too large to store on GitHub; this data is outputted from the ```import_multiple``` parsing function. The other two ```.csv``` files are from an earlier study by [Gillan *et al.* (2016)](https://elifesciences.org/articles/11305), the use of which is explained further in [this notebook](https://github.com/qdercon/pstpipeline/blob/main/notebooks/data_cleaning_factor_derivation.ipynb).

The data object is structured as a nested list.  As analyses were done separately in distancing participants, it is first separated by group (```non-distanced``` and ```distanced```). Within each group, there are then four data frames: ```ppt_info``` (demographic and exclusion info); ```training``` (training phase data); ```test``` (test phase data); and ```gillan_questions``` (psychiatric questionnaire question answers).

Note that with the exception of questionnaire questions, the data includes all 995 participants, so you may want to filter out excluded participants. The ```ppt_info``` data frame includes a column called ```exclusion``` which is ```TRUE``` for excluded participants and ```FALSE``` for non-excluded participants (see the [preprint](https://psyarxiv.com/jmnek) for more information on exclusion criteria).

If you use the data for any published work, please consider citing the preprint:
> Q. Dercon<sup>†</sup>, S. Z. Mehrhof<sup>†</sup>, T. R. Sandhu, C. Hitchcock, R. P. Lawson, D. A. Pizzagalli, T. Dalgleish, C. L. Nord. A Core Component of Psychological Therapy Causes Adaptive Changes in Computational Learning Mechanisms. *PsyArXiv* (2022). https://psyarxiv.com/jmnek.

## In R

The ```.RDS``` file can be easily loaded into R using the following code:

```
all_res_split <- readRDS("data-raw/all_res_split.RDS")
```
Then, the various data frames can be accessed using ```$``` or ```[[ ]]``` syntax. For example, to access the ```training``` data frame for the ```non-distanced``` group, you would call:

```
training_data <- all_res_split$non_distanced$training
```
If required, an easy way to remove the distancing group nesting would be as follows (the pipe requires ```R >= 4.1.0```):

```
all_res <- purrr::transpose(all_res_split) |> purrr::map(dplyr::bind_rows)
```

## In Python

It is probably easiest to first import the data into R as detailed above and then save the data you require as ```.csv``` files. However, should you wish to do this in Python, you can use the ```rpy2``` package to convert the RDS into ```pandas``` data frames. This can be run in Google Colab to avoid local installs.

```
import pandas as pd
import rpy2.robjects as ro
import urllib

from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

## download data from GitHub
url = "https://github.com/qdercon/pstpipeline/blob/main/data-raw/all_res_split.RDS?raw=true"
urllib.request.urlretrieve(url, "all_res_split.RDS")

## load RDS file
readRDS = ro.r['readRDS']
r_df = readRDS('all_res_split.RDS')

## convert to rpy2 ListVector that pandas can read
with localconverter(ro.default_converter + pandas2ri.converter):
  non_distanced = ro.conversion.rpy2py(r_df[0])
  distanced = ro.conversion.rpy2py(r_df[0])

## convert training data to pandas data frames
with localconverter(ro.default_converter + pandas2ri.converter):
  non_dis_train = ro.conversion.rpy2py(non_distanced.rx('training')[0])
  dis_train = ro.conversion.rpy2py(distanced.rx('training')[0])

non_dis_train.head()
```
