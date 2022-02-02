## Required libraries
library('getopt')

## Specify parameters
spec <- matrix(c(
	'spefile', 's', 2, 'character', 'SPE file name',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## Rename from spe_targeted to spe to simplify the code
if(opt$spefile == "spe_targeted_postqc.Rdata") {
    spe <- spe_targeted
}



## Save new SPE objects
if(opt$spefile == "spe_targeted_postqc.Rdata") {
    spe_targeted <- spe
    ## First time switching the order of the keywords: now targeted is at the
    ## end, which will make it easier to sort the spe files.
    save(spe_targeted, file = file.path(dir_rdata, "spe_harmony_targeted.Rdata"))
} else {
    ## First time using "wholegenome" in the spe name, to clearly identify it
    save(spe, file.path(dir_rdata, "spe_harmony_wholegenome.Rdata"))
}
