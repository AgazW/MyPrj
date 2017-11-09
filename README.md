 
To get started with the tool remotely (not as web application), install R base from the link https://cran.r-project.org/bin/windows/base/. Install RStudio compatible to your machine from https://www.rstudio.com/products/rstudio/ download/. Once you are done, start RStudio and enter the following commands to install the needed packages.

    install.packages("shiny")

    install.packages("markdown")

    install.packages("plyr")

    install.packages("data.table")

    install.packages("dplyr")

## Affy package
    source("http://bioconductor.org/biocLite.R")
    biocLite("affy")

## simpleaffy package
    source("http://bioconductor.org/biocLite.R")
    biocLite("simpleaffy")

## lumi package
    source("http://bioconductor.org/biocLite.R")
    biocLite("lumi")
 
Once you are done with all the installations put the tool folder in the current working directory. You can also set your working directory at the top of the Rstudio session -> Set Working Directory -> Choose Directory. After setting you can check the current working directory by getwd(). To run the program, load shiny package by typing  library(shiny) in your R consol and then enter runApp('shinyMDE') i,e the name of tool folder into R, open in browser and work.

Check [wiki](https://github.com/AgazW/shinyMDE/wiki) for more help.
