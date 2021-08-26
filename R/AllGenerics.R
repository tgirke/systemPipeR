#####################################
## Generic Definitions for SYSargs ##
#####################################
## Methods to return SYSargs components
setGeneric(name = "targetsin", def = function(x) standardGeneric("targetsin"))
setGeneric(name = "targetsout", def = function(x) standardGeneric("targetsout"))
setGeneric(name = "targetsheader", def = function(x, ...) standardGeneric("targetsheader"))
setGeneric(name = "modules", def = function(x) standardGeneric("modules"))
setGeneric(name = "software", def = function(x) standardGeneric("software"))
setGeneric(name = "cores", def = function(x) standardGeneric("cores"))
setGeneric(name = "other", def = function(x) standardGeneric("other"))
setGeneric(name = "reference", def = function(x) standardGeneric("reference"))
setGeneric(name = "results", def = function(x) standardGeneric("results"))
setGeneric(name = "infile1", def = function(x, ...) standardGeneric("infile1"))
setGeneric(name = "infile2", def = function(x, ...) standardGeneric("infile2"))
setGeneric(name = "outfile1", def = function(x) standardGeneric("outfile1"))
setGeneric(name = "SampleName", def = function(x, ...) standardGeneric("SampleName"))
setGeneric(name = "sysargs", def = function(x) standardGeneric("sysargs"))
setGeneric(name = "outpaths", def = function(x) standardGeneric("outpaths"))

######################################
## Generic Definitions for SYSargs2 ##
######################################
## Methods to return SYSargs2 components
setGeneric(name = "targets", def = function(x) standardGeneric("targets"))
setGeneric(name = "wf", def = function(x) standardGeneric("wf"))
setGeneric(name = "clt", def = function(x) standardGeneric("clt"))
setGeneric(name = "yamlinput", def = function(x, ...) standardGeneric("yamlinput"))
setGeneric(name = "cmdlist", def = function(x, ...) standardGeneric("cmdlist"))
setGeneric(name = "input", def = function(x) standardGeneric("input"))
setGeneric(name = "output", def = function(x) standardGeneric("output"))
setGeneric(name = "files", def = function(x) standardGeneric("files"))
setGeneric(name = "inputvars", def = function(x) standardGeneric("inputvars"))
setGeneric(name = "cmdToCwl", def = function(x) standardGeneric("cmdToCwl"))
setGeneric(name = "status", def = function(x) standardGeneric("status"))
## Coerce
setGeneric(name = "sysargs2", def = function(x) standardGeneric("sysargs2"))
## Accessors
setGeneric(name = "baseCommand", def = function(x, ...) standardGeneric("baseCommand"))
## Replacement methods
setGeneric(name = "yamlinput<-", def = function(x, paramName, ..., value) standardGeneric("yamlinput<-"))
setGeneric(name = "cmdToCwl<-", def = function(x, ..., value) standardGeneric("cmdToCwl<-"))

#########################################
## Generic Definitions for SYSargsList ##
#########################################
## Methods to return SYSargsList components
setGeneric(name = "stepsWF", def = function(x) standardGeneric("stepsWF"))
setGeneric(name = "statusWF", def = function(x) standardGeneric("statusWF"))
setGeneric(name = "targetsWF", def = function(x) standardGeneric("targetsWF"))
setGeneric(name = "outfiles", def = function(x) standardGeneric("outfiles"))
setGeneric(name = "SEobj", def = function(x) standardGeneric("SEobj"))
setGeneric(name = "dependency", def = function(x) standardGeneric("dependency"))
setGeneric(name = "projectInfo", def = function(x) standardGeneric("projectInfo"))
setGeneric(name = "runInfo", def = function(x) standardGeneric("runInfo"))
## Accessors
setGeneric(name = "stepName", def = function(x) standardGeneric("stepName"))
setGeneric(name = "getColumn", def = function(x, step, position = c("outfiles", "targetsWF"), column = 1, names = SampleName(x, step)) standardGeneric("getColumn"))
# setGeneric(name = "updateColumn", def = function(x, step, df, position = c("outfiles", "targetsWF")) standardGeneric("updateColumn"))
setGeneric(name = "viewEnvir", def = function(x) standardGeneric("viewEnvir"))
setGeneric(name = "copyEnvir", def = function(x, list = character(), new.env = globalenv(), silent=FALSE) standardGeneric("copyEnvir"))
## Coerce back to list: as(SYSargsList, "list")
setGeneric(name = "sysargslist", def = function(x) standardGeneric("sysargslist"))
## Replacement methods
setGeneric(name = "updateColumn<-", def = function(x, step, position = c("outfiles", "targetsWF"), value) standardGeneric("updateColumn<-"))
setGeneric(name = "appendStep<-", def = function(x, after = length(x), ..., value) standardGeneric("appendStep<-"))
setGeneric(name = "replaceStep<-", def = function(x, step, step_name = "default", value) standardGeneric("replaceStep<-"))
setGeneric(name = "stepsWF<-", def = function(x, step, ..., value) standardGeneric("stepsWF<-"))
setGeneric(name = "renameStep<-", def = function(x, step, ..., value) standardGeneric("renameStep<-"))
setGeneric(name = "statusWF<-", def = function(x, step, ..., value) standardGeneric("statusWF<-"))
setGeneric(name = "dependency<-", def = function(x, step, ..., value) standardGeneric("dependency<-"))
setGeneric(name = "updateStatus<-", def = function(x, step, ..., value) standardGeneric("updateStatus<-"))

######################################
## Generic Definitions for LineWise ##
######################################
## Methods to return LineWise components
setGeneric(name = "codeLine", def = function(x, ...) standardGeneric("codeLine"))
setGeneric(name = "rmdPath", def = function(x) standardGeneric("rmdPath"))
setGeneric(name = "codeChunkStart", def = function(x) standardGeneric("codeChunkStart"))
## Coerce back to list: as(LineWise, "list")
setGeneric(name = "linewise", def = function(x) standardGeneric("linewise"))
## Replacement methods
setGeneric(name = "replaceCodeLine<-", def = function(x, line, ..., value) standardGeneric("replaceCodeLine<-"))
setGeneric(name = "appendCodeLine<-", def = function(x, after = length(x), ..., value) standardGeneric("appendCodeLine<-"))
