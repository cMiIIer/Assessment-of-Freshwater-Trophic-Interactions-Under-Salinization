# Project Objective

The objective of this study is to quantify the trophic interactions of a
model freshwater community under the effects of salinzation. There are
two main objectives that were the goals of the  data analysis phase.

1.  Quantify the complex relationships of the system, looking how
    physiochemical variables, primarily salt content, effects biomass
    accumulation and transfer across our trophic levels.

2.  Quantify success of our priamry species of interest, *Hyla cinerea* at
    our various salinity treatments (i.e. survial and metamorphosis)

# Our Findings

To learn avout the project, we have attatched the Thesis manuscript. This is 
where you can find the specifics of our methods and see what we found 
with our analysis. 

If you are simply looking for our code structures, you need to read it as 
we provide them in the repository.

   
# Replication

If you are looking to replicated our methods or results, either for your
own research or for validation of our work, we provide the code and data
in this repository. None of the analystical methods we employed were
based on random processes (e.g. simulated data or random walk estimations)
nor permutational methods, so results should be primarily consistent. In the
case of our MCMC alrogrithms, where posterior draws can be inconsistent, we 
employed a seed.

## Code Base Structure and Dependencies

The data and code are stored within the eponomyously named code structures. For 
objective one of the study, please refer to files denoted with "SEM" and 
objective two "TTE". 

We employ mtwo programming languages, R and Stan, we work with both 
using the RStudio IDE. For replication, you need work with the Stan
code directly, just use the R files as it is how we interact with our
Stan files. We use a suite of packages in R, most notably ones
to interface with the Stan ecosystem, so you will need the following packages
installed in your R enviornment to succesfully run our code. All packages 
are loaded within the R scripts, but if you want to install them before hand
if you do not have them, we provide our package library below.

### Requried Packages

    install.packages(
      c(
        "tidyverse", 
        "rstan", 
        "loo", 
        "glmmTMB"
        )
      )

## Data Strcuture

The data were cleaned and formatted beforehand, so you need not run make any changes
to the data except what has been done in code base. All the code has been formattted
to run sequentially, so there is no need to run things out of order. 

## Working Directory

You can download the entire repository and use it as a project, as we included an R project file to help
with integrating this into your own R workspace.


# Acknowledgements

I would like to thank my Advisor, Dr. Allison Welch, and two
undergraduate coauthors: Maya Mylott and Alex Barron. This would also not 
be possible without my advisory commitee Dr. A. Challen Hyman, Dr. Tony Harold, Dr. Norm Levine and
Dr. Dan McGlinn. Further support and assistance was provided by Greg Townsley, Pete Meier, Dr. Bob 
Podoslky and Dr. Freya Rowland.
