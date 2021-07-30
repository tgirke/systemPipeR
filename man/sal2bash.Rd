\name{sal2bash}
\alias{sal2bash}
\title{
Translate SYSargsList back to a bash workflow
}
\description{
    This function takes a SYSargsList object and translate it to an 
    executable bash script, so one can run the workflow without loading SPR or 
    using an R console.
}
\usage{
sal2bash(sal, out_dir = ".", bash_path = "/bin/bash", stop_on_error = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sal}{
    \code{SYSargsList} object.
}
  \item{out_dir}{
    string, a relative or absolute path to a directory. If the 
    directory does not exist, this function will first try to create it. 
}
  \item{bash_path}{
    string, the path to the bash executable program
}
  \item{stop_on_error}{
  bool, should the bash script stop if any error happens in running
}
}
\details{
    ## out files 
    
    1. The main executable bash file will be created to the root of `out_dir`
    
    2. All R steps will be stored as R scripts and along with other supporting 
    files inside a folder called `spr_bash` under `out_dir`
    
    3. Not all R steps will have an individual file. This function will "collapse"
    adjacent R steps into one file as much as possible. Namely, if there is no 
    sysArgs steps in between, R steps will be merged into one file, otherwise 
    they will be in different files. 
    
    ## R steps 
    
    Similarly as running the workflow in R console, all R steps will share the same 
    environment variables and loaded packages. This is done by loading and saving 
    the R environment into a file `spr_wf.RData` before and after the R script execution.
    Therefore, it will be good to keep all R steps bundle together as much as possible 
    to avoid the package and environment loading/saving overhead time. 
    
    Initially, this environment only contains the SYSargsList object that was used to 
    create the bash script. Note: the SYSargsList object name will be the same 
    as what you pass to `sal2bash`. If you have `sal2bash(my_sal)`, then in the 
    `spr_wf.RData` there will be an object called `my_sal` saved there. It is important 
    to keep using the same name for the SYSargsList object managing the workflow 
    and in workflow running. For example, if you have an R step that requires to 
    query a column from the outfiles of the SAL, `getColumn(sal, "FileName1")`, but you 
    pass `my_sal` to `sal2bash(my_sal)`, this will cause this R step cannot find 
    the SAL object when run the workflow from bash.  
    
    ## Execution
    
    This way of execution is not able to handle complex dependency graphs. The original 
    step dependencies from SAL object will be ignored, so all steps will be executed 
    in a linear manner. It is recommended to adjust the workflow order before using 
    this function. 
}
\value{
    no return
}
\author{
    Le Zhang and Daniela Cassol
}
\examples{
file_path <- system.file("extdata/spr_simple_wf.Rmd", package="systemPipeR")
sal <- SPRproject(overwrite = TRUE)
sal <- importWF(sal, file_path)
sal2bash(sal)
}
\keyword{ utilities }

