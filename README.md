# Contribution of systemic and somatic factors to clinical response and resistance in urothelial cancer: An exploratory multi-omic analysis

## Running the Code

To get started, install the requirements:

```
pip install -r requirements.txt
```

Create an `ENV.sh` modeled after `ENV_TEMPLATE.sh`, pointing to the data that you have available, and then call `run.sh` to `source` your `ENV.sh` in the context of a new Jupyter notebook.

### rpy2

Certain notebooks in this repo require `rpy2`. `rpy2` will be installed via the `requirements.txt`, but depending on your environment may require additional setup.

In order to execute the notebooks, you will need to have:

1. a functioning R install, preferably of a recent version of R. 
2. certain R libraries commonly used. 

Specifically:
```
# install R, if you haven't already
sudo apt-get install r-base r-base-dev

# create & set personal rlib directory, if not already done
# easiest way to do this is to open an interactive R console, and run 
> options(repos = 'https://cran.rstudio.com')
> install.packages('ggplot2')

## install certain R packages
R -e "options(repos = 'https://cran.rstudio.com'); install.packages(c('dplyr','survival','tidyr','ggplot2'));"
```
