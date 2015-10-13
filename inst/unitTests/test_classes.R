testDir <- system.file("unitTests", package="systemPipeR")
test.import <- function(){
	## Test import for generating SYSargs instance
	args <- systemArgs(sysma=file.path(testDir, "tophat.param"), mytargets=file.path(testDir, "targets.txt"))	
	checkTrue(class(args)=="SYSargs")
}
