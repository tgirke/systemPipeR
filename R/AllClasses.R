#######################################
## Class Definitions for SYSargsList ##
#######################################
## Define SYSargsList class
setClass("SYSargsList", slots = c(
    stepsWF = "list",
    statusWF = "list",
    targetsWF = "list",
    outfiles = "list",
    SEobj = "list",
    dependency = "list",
    targets_connection = "list",
    projectInfo = "list",
    runInfo = "list"
))

####################################
## Class Definitions for LINEWISE ##
####################################
## Define LineWise class
setClass("LineWise", slots = c(
    codeLine = "expression",
    codeChunkStart = "integer",
    #rmdPath = "character",
    stepName = "character",
    dependency = "list",
    status = "list", 
    files = "list", 
    runInfo = "list"
))

####################################
## Class Definitions for SYSargs2 ##
####################################
## Define SYSargs2 class
setClass("SYSargs2", slots = c(
    targets = "list",
    targetsheader = "list",
    modules = "list",
    wf = "list",
    clt = "list",
    yamlinput = "list",
    cmdlist = "list",
    input = "list",
    output = "list",
    files = "list",
    inputvars = "list",
    cmdToCwl = "list",
    status = "list",
    internal_outfiles = "list"
))

####################################
## Class Definitions for SYSargs ##
###################################
## Define SYSargs class
setClass("SYSargs", representation(
    targetsin = "data.frame",
    targetsout = "data.frame",
    targetsheader = "character",
    modules = "character",
    software = "character",
    cores = "numeric",
    other = "character",
    reference = "character",
    results = "character",
    infile1 = "character",
    infile2 = "character",
    outfile1 = "character",
    sysargs = "character",
    outpaths = "character"
))

