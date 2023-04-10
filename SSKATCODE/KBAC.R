if (!require("devtools", character.only=TRUE, quietly=TRUE)) {
  install.packages("devtools")
}
remotes::install_github("AssotesteR")
library("devtools")

download.file('https://github.com/gaow/kbac', 
              f <- tempfile())
unzip(f, exdir=tempdir())
file.copy(file.path(tempdir(), '.RData'), 'bivpois.RData')
# the above copies the .RData file to a file called bivpois.RData in your current 
# working directory.
load('bivpois.RData')

install.packages("KBAC", dependencies = TRUE)


install.packages("mvtnorm")

install.packages("AssotesteR")
install.packages("https://cran.r-project.org/package=AssotesteR")
