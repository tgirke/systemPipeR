\name{SYSargs2-class}
\Rdversion{1.1}
\docType{class}
% Class
\alias{SYSargs2-class}
% Subsetting and Quick access
\alias{[,SYSargs2,ANY,ANY,ANY-method}
\alias{$,SYSargs2-method}

% Coercion
\alias{coerce,list,SYSargs2-method}
\alias{coerce,SYSargs2,list-method}
% Replacement
\alias{[[,SYSargs2,ANY,missing-method}
\alias{[[<-,SYSargs2,ANY,ANY,ANY-method}

%Accessors
\alias{sysargs2}
\alias{sysargs2,SYSargs2-method}
\alias{length,SYSargs2-method}
\alias{names,SYSargs2-method}
\alias{show,SYSargs2-method}
\alias{modules,SYSargs2-method}
\alias{targets}
\alias{targets,SYSargs2-method}
\alias{targetsheader,SYSargs2-method}
\alias{infile1,SYSargs2-method}
\alias{infile2,SYSargs2-method}
\alias{wf}
\alias{wf,SYSargs2-method}
\alias{clt}
\alias{clt,SYSargs2-method}
\alias{yamlinput}
\alias{yamlinput,SYSargs2-method}
\alias{yamlinput<-}
\alias{yamlinput<-,SYSargs2-method}
\alias{yamlinput<-,SYSargs2,ANY,ANY-method}
\alias{cmdlist}
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
\alias{baseCommand}
\alias{baseCommand,SYSargs2-method}

\alias{SampleName,SYSargs2-method}


\title{Class \code{"SYSargs2"}}
\description{
\code{SYSargs2 class} stores all the information and instructions needed for processing a set of input files with a specific command-line or a series of command-line within a workflow. The \code{SYSargs2} S4 class object is created from the \code{loadWorkflow} and \code{renderWF} function, which populates all the command-line for each sample in each step of the particular workflow. Each sample level input/outfile operation uses its own \code{SYSargs2} instance. The output of \code{SYSargs2} define all the expected output files for each step in the workflow, which usually it is the sample input for the next step in an \code{SYSargs2} instance. Between different instances, this connectivity is established by writing the subsetting output with the \code{writeTargetsout} function to a new targets file that serves as input to the next \code{loadWorkflow} and \code{renderWF} call. By chaining several \code{SYSargs2} steps together one can construct complex workflows involving many sample-level input/output file operations with any combination of command-line or R-based software. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SYSargs2", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{targets}:}{Object of class \code{"list"} storing data from each sample from \code{targets} file }
    \item{\code{targetsheader}:}{Object of class \code{"list"} storing header/comment lines of targets file }
    \item{\code{modules}:}{Object of class \code{"list"} storing software versions from module system }
    \item{\code{wf}:}{Object of class \code{"list"} storing data from Workflow \code{CWL param} file }
    \item{\code{clt}:}{Object of class \code{"list"} storing data from each CommandLineTool step in the Workflow or the single CommandLineTool \code{CWL param} file }
    \item{\code{yamlinput}:}{Object of class \code{"list"} storing data from \code{input} file }
    \item{\code{cmdlist}:}{Object of class \code{"list"} storing all command-line used to execute external software }
    \item{\code{input}:}{Object of class \code{"list"} storing data from each target defined in \code{inputvars} }
    \item{\code{output}:}{Object of class \code{"list"} paths to final outputs files }
    \item{\code{files}:}{Object of class \code{"list"} paths to \code{input} and \code{CWL param} files }
    \item{\code{inputvars}:}{Object of class \code{"list"} storing data from each \code{inputvars} }
    \item{\code{cmdToCwl}:}{Object of class \code{"list"} storing data from each \code{cmdToCwl} }
    \item{\code{status}:}{Object of class \code{"list"} storing data from each \code{status} }
    \item{\code{internal_outfiles}:}{Object of class \code{"list"} storing raw data from each \code{output} }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "SYSargs2")}: subsetting of class with bracket operator }
    \item{[[}{\code{signature(x = "SYSargs2", i = "ANY", j = "missing")}: subsetting of class with bracket operator }
    \item{[[<-}{\code{signature(x = "SYSargs2")}: replacement method for SYSargs2 class }
    \item{$}{\code{signature(x = "SYSargs2")}: extracting slots elements by name }
    \item{clt}{\code{signature(x = "SYSargs2")}: extracts data from \code{clt} slot }
    \item{cmdlist}{\code{signature(x = "SYSargs2")}: extracts data from \code{cmdlist} slot }
    \item{coerce}{\code{signature(from = "list", to = "SYSargs2")}: \code{as(list, "SYSargs2")} }
    \item{files}{\code{signature(x = "SYSargs2")}: extracts data from \code{files} slot }
    \item{infile1}{\code{signature(x = "SYSargs2")}: extracting paths to first FASTQ file }
    \item{infile2}{\code{signature(x = "SYSargs2")}: extracting paths to second FASTQ file if data is PE }
    \item{input}{\code{signature(x = "SYSargs2")}: extracts data from \code{input} slot }
    \item{inputvars}{\code{signature(x = "SYSargs2")}: extracts data from \code{inputvars} slot }
    \item{length}{\code{signature(x = "SYSargs2")}: extracts number of samples }
    \item{modules}{\code{signature(x = "SYSargs2")}: extracts data from \code{modules} slot }
    \item{names}{\code{signature(x = "SYSargs2")}: extracts slot names }
    \item{output}{\code{signature(x = "SYSargs2")}: extracts data from \code{cmdlist} slot }
    \item{show}{\code{signature(object = "SYSargs2")}: summary view of \code{SYSargs2} objects }
    \item{sysargs2}{\code{signature(x = "SYSargs2")}: Coerce back to list \code{as(SYSargs2, "list")} }
    \item{targets}{\code{signature(x = "SYSargs2")}: extract data from \code{targets} slot }
    \item{targetsheader}{\code{signature(x = "SYSargs2")}: extracts data from \code{targetsheader} slot }
    \item{wf}{\code{signature(x = "SYSargs2")}:  extracts data from \code{wf} slot }
    \item{yamlinput}{\code{signature(x = "SYSargs2")}:  extracts data from \code{yamlinput} slot }
    \item{cmdToCwl}{\code{signature(x = "SYSargs2")}:  extracts data from \code{cmdToCwl} slot }
    \item{status}{\code{signature(x = "SYSargs2")}:  extracts data from \code{status} slot }
    \item{yamlinput<-}{\code{signature(x = "SYSargs2")}: replacement method for \code{yamlinput} slot }
	 }
}
\author{
Daniela Cassol and Thomas Girke
}
\seealso{
\code{loadWorkflow} and \code{renderWF} and \code{runCommandline} and \code{clusterRun}
}
\examples{
showClass("SYSargs2")

## Construct SYSargs2 object from CWl param, CWL input, and targets files 
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl", package="systemPipeR")
WF <- loadWorkflow(targets=targets, wf_file="hisat2/hisat2-mapping-se.cwl", 
                  input_file="hisat2/hisat2-mapping-se.yml", dir_path=dir_path)
WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
WF
names(WF); modules(WF); targets(WF)[1]; cmdlist(WF)[1:2]; output(WF)

## The subsetting operator '[' allows to select specific command-line
cmdlist(WF)[1:2]

\dontrun{
## Execute SYSargs2 on single machine
WF <- runCommandline(args=WF)

## Execute SYSargs2 on multiple machines of a compute cluster. The following
## example uses the conf and template files for the Slurm scheduler. Please 
## read the instructions on how to obtain the corresponding files for other schedulers.  
file.copy(system.file("extdata", ".batchtools.conf.R", package="systemPipeR"), ".")
file.copy(system.file("extdata", "batchtools.slurm.tmpl", package="systemPipeR"), ".")
resources <- list(walltime=120, ntasks=1, ncpus=4, memory=1024) 
reg <- clusterRun(args, FUN = runCommandline, conffile=".batchtools.conf.R", template="batchtools.slurm.tmpl", Njobs=18, runid="01", resourceList=resources)

## Monitor progress of submitted jobs
getStatus(reg=reg)

## Updates the path in the object \code{output(WF)}
WF <- output_update(WF, dir=FALSE, replace=TRUE, extension=c(".sam", ".bam"))

## Alignment stats
read_statsDF <- alignStats(WF) 
read_statsDF <- cbind(read_statsDF[targets$FileName,], targets)
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
}
}
\keyword{ classes }
