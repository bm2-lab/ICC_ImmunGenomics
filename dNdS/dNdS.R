args <- commandArgs()

input_file=args[6]

out_file=args[7]

library("dndscv")

mutations=read.table(input_file,header=T,sep='\t')

dndsout = dndscv(mutations)

sel_cv = dndsout$sel_cv


write.table(dndsout$globaldnds,file=out_file,sep="\t",col.names = TRUE,quote = FALSE,row.names =FALSE)