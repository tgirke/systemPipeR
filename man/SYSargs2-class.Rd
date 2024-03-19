\name{SYSargs2-class}
\Rdversion{1.1}
\docType{class}
% Class
\alias{SYSargs2-class}
% Subsetting and Quick access
\alias{[,SYSargs2,ANY,ANY,ANY-method}
\alias{[[,SYSargs2,ANY,missing-method}
\alias{[[<-,SYSargs2,ANY,ANY,ANY-method}
\alias{[[<-,SYSargs2,ANY,ANY-method}
\alias{$,SYSargs2-method}
\alias{length,SYSargs2-method}
\alias{names,SYSargs2-method}
\alias{show,SYSargs2-method}
% Coercion
\alias{coerce,list,SYSargs2-method}
\alias{coerce,SYSargs2,list-method}
\alias{sysargs2}
\alias{sysargs2,SYSargs2-method}
%Accessors
\alias{targets}
\alias{targets,SYSargs2-method}
\alias{targetsheader,SYSargs2-method}
\alias{modules}
\alias{modules,SYSargs2-method}
\alias{infile1,SYSargs2-method}
\alias{infile2,SYSargs2-method}
\alias{wf}
\alias{wf,SYSargs2-method}
\alias{clt}
\alias{clt,SYSargs2-method}
%\alias{yamlinput}
\alias{yamlinput,SYSargs2-method}
%\alias{cmdlist}
\alias{cmdlist,SYSargs2-method}
\alias{input}
\alias{input,SYSargs2-method}
\alias{output}
\alias{output,SYSargs2-method}
\alias{files}
\alias{files,SYSargs2-method}
\alias{inputvars}
\alias{inputvars,SYSargs2-method}
\alias{cmdToCwl}
\alias{cmdToCwl,SYSargs2-method}
\alias{status}
\alias{status,SYSargs2-method}
\alias{internal_outfiles}
\alias{internal_outfiles,SYSargs2-method}
%Methods
%\alias{baseCommand}
\alias{baseCommand,SYSargs2-method}
\alias{SampleName,SYSargs2-method}
% Replacement
%\alias{yamlinput<-}
\alias{yamlinput<-,SYSargs2-method}
\alias{yamlinput<-,SYSargs2,ANY,ANY-method}
\alias{cmdToCwl<-}
\alias{cmdToCwl<-,SYSargs2-method}

\title{Class \code{"SYSargs2"}}
\description{
\code{SYSargs2 class} stores all the information and instructions needed for 
processing a set of input files with a specific command-line or a series of 
command-line within a workflow. The \code{SYSargs2} S4 class object is created 
from the \code{loadWF} and \code{renderWF} function, which populates all the 
command-line for each sample in each step of the particular workflow. Each sample 
level input/outfile operation uses its own \code{SYSargs2} instance. The output 
of \code{SYSargs2} define all the expected output files for each step in the 
workflow, which usually it is the sample input for the next step in an \code{SYSargs2}
instance. By chaining several \code{SYSargs2} steps together one can construct 
complex workflows involving many sample-level input/output file operations with 
any combination of command-line or R-based software. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SYSargs2", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{targets}:}{Object of class \code{"list"} storing data from each sample from \code{targets} file }
    \item{\code{targetsheader}:}{Object of class \code{"list"} storing header/comment lines of targets file }
    \item{\code{modules}:}{Object of class \code{"list"} storing software versions from module system }
    \item{\code{wf}:}{Object of class \code{"list"} storing data from Workflow \code{CWL parameters} file }
    \item{\code{clt}:}{Object of class \code{"list"} storing data from each \code{CommandLineTool}
                        substep in the Workflow or the single \code{CommandLineTool} \code{CWL parameters} file }
    \item{\code{yamlinput}:}{Object of class \code{"list"} storing data from \code{input (*.yml)} file }
    \item{\code{cmdlist}:}{Object of class \code{"list"} storing all command-line used to execute external software }
    \item{\code{input}:}{Object of class \code{"list"} storing data from each target defined in \code{inputvars} }
    \item{\code{output}:}{Object of class \code{"list"} paths to final outputs files }
    \item{\code{files}:}{Object of class \code{"list"} paths to \code{input} and \code{CWL parameters} files }
    \item{\code{inputvars}:}{Object of class \code{"list"} storing data from each \code{inputvars} }
    \item{\code{cmdToCwl}:}{Object of class \code{"list"} storing data from each \code{cmdToCwl} }
    \item{\code{status}:}{Object of class \code{"list"} storing data from each \code{status} }
    \item{\code{internal_outfiles}:}{Object of class \code{"list"} storing raw data from each \code{output} }
  }
}
\section{Methods}{
  \describe{
    \item{[}{ Subsetting of class with bracket operator. }
    \item{[[}{ Subsetting of class with bracket operator. }
    \item{[[<-}{ Replacement method for \code{"SYSargs2"} class. }
    \item{$}{ Extracting slots elements by name. }
    \item{length}{ Extracts number of samples. }
    \item{names}{ Extracts slot names. }
    \item{show}{ Summary view of \code{SYSargs2} objects. }

    % Coerce
    \item{coerce}{ \code{signature(from = "list", to = "SYSargs2")}: \code{as(list, "SYSargs2")} }
    \item{coerce}{ \code{signature(from = "SYSargs2", to = "list")}\code{as(SYSargs2, "list")} }
    \item{coerce}{ \code{signature(from = "SYSargs2", to = "DataFrame")}: \code{as(x, "DataFrame")}; for targets slot. }
    \item{sysargs2}{ Coerce back to list \code{as(SYSargs2, "list")} }
    
    % Accessors
    \item{targets}{ Extract data from \code{targets} slot. }
    \item{targetsheader}{ Extracts data from \code{targetsheader} slot. }
    \item{modules}{ Extracts data from \code{modules} slot. }
    \item{wf}{ Extracts data from \code{wf} slot. }
    \item{clt}{ Extracts data from \code{clt} slot. }
    \item{yamlinput}{ Extracts data from \code{yamlinput} slot. }
    \item{cmdlist}{ Extracts data from \code{cmdlist} slot. }
    \item{input}{ Extracts data from \code{input} slot. }
    \item{output}{ Extracts data from \code{cmdlist} slot. }
    \item{files}{ Extracts data from \code{files} slot. }
    \item{inputvars}{ Extracts data from \code{inputvars} slot. }
    \item{cmdToCwl}{ Extracts data from \code{cmdToCwl} slot. }
    \item{status}{ Extracts data from \code{status} slot. }
    
    % Methods
    \item{infile1}{extracting paths to first FASTQ file. }
    \item{infile2}{extracting paths to second FASTQ file if data is PE. }
    \item{baseCommand}{ Extracts baseCommand from command-line used to execute external software. }
    \item{SampleName}{ Extracts all samples names. }

        % Replacement
    \item{yamlinput<-}{ Replacement method for \code{yamlinput} slot input. }
	 }
}
\author{
Daniela Cassol and Thomas Girke
}
\seealso{
\code{loadWF} and \code{renderWF} and \code{runCommandline} and \code{clusterRun}
}
\examples{
showClass("SYSargs2")

## Construct SYSargs2 object from CWl param, CWL input, and targets files 
targetspath <- system.file("extdata/cwl/example/targets_example.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl", package="systemPipeR")
WF <- loadWorkflow(targets=targetspath, 
                  wf_file="example/workflow_example.cwl", 
                  input_file="example/example.yml", 
                  dir_path=dir_path)
WF <- renderWF(WF, inputvars=c(Message = "_STRING_", SampleName = "_SAMPLE_"))
WF

## Methods
names(WF)
length(WF)
baseCommand(WF)
SampleName(WF)

## Accessors
targets(WF)
targetsheader(WF)
modules(WF)
yamlinput(WF)
cmdlist(WF)
input(WF)
output(WF)
files(WF)
inputvars(WF)
cmdToCwl(WF)
status(WF)

## The subsetting operator '[' allows to select specific command-line/sample
WF2 <- WF[1:2]

\dontrun{
## Execute SYSargs2 on single machine
WF2 <- runCommandline(WF2)
}
\dontrun{
## Execute SYSargs2 on multiple machines of a compute cluster. The following
## example uses the conf and template files for the Slurm scheduler. Please 
## read the instructions on how to obtain the corresponding files for other schedulers.  
file.copy(system.file("extdata", ".batchtools.conf.R", package="systemPipeR"), ".")
file.copy(system.file("extdata", "batchtools.slurm.tmpl", package="systemPipeR"), ".")
resources <- list(walltime=120, ntasks=1, ncpus=4, memory=1024) 
reg <- clusterRun(WF, FUN = runCommandline, conffile=".batchtools.conf.R", 
                  template="batchtools.slurm.tmpl", Njobs=2, runid="01", resourceList=resources)

## Monitor progress of submitted jobs
getStatus(reg=reg)

}
}
\keyword{ classes }
