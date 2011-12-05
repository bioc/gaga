groups2int <- function(groups,patterns) {
#check that names in groups and patterns match
if (is.null(colnames(patterns))) stop('You must specify colnames(patterns)')
if (sum(unique.default(groups)[order(unique.default(groups))]==colnames(patterns)[order(colnames(patterns))])<ncol(patterns)) stop('Group names in colnames(patterns) do no match group names indicated in groups')
#convert groups to integer vector
groupsr <- integer(length(groups))
for (i in 1:ncol(patterns)) { groupsr[groups==colnames(patterns)[i]] <- i-1 }
groupsr <- as.integer(groupsr)
return(groupsr)
}
