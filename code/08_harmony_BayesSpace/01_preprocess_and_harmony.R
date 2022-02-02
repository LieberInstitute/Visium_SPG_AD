## Required libraries
library('getopt')
library('devtools')

## Specify parameters
spec <- matrix(c(
	'outdir', 'o', 2, 'character', 'Path where the output of infer_experiment.py was saved. Defaults to HISAT2_out/infer_strandess',
    'pattern', 'p', 2, 'character', 'Name of the pattern file. Defaults to inferred_strandness_pattern.txt',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}
