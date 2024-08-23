library(optparse)
# Define options
option_list <- list(
make_option(c("-c", "--crm"), type = "character", default = "Y", help = "Chromosome"),
make_option(c("-p", "--positionlist"), type = "character", default = "None", help = "Comma separated list of chromosome positions"),
make_option(c("-m", "--datadir"), type = "character", default = "None", help = "Directory with hypergeometric data"),
make_option(c("-o", "--outputdir"), type = "character", default = "/tmp/plotX", help = "Directory for the result plot"),
make_option(c("-t", "--type"), type = "character", default = "median", help = "DataType - median or hypergeom")
)
# Create the parser object
opt_parser <- OptionParser(option_list = option_list)
# Parse the arguments
opt <- parse_args(opt_parser)

# Access the arguments
#cat("Chromosome:", opt$crm, "\n")
#cat("datadir:", opt$datadir, "\n")
#cat("outputDir:", opt$outputdir, "\n")
#cat("positionlist:", opt$positionlist, "\n")

datadir <- opt$datadir
crm <- opt$crm
outdir <- opt$outputdir
positions <- strsplit(opt$positionlist,",")[[1]]
dataType <- opt$type

stopifnot(positions[1] >0)

## set working directory
setwd(outdir)
# import necessary libraries
library(ggplot2)
library(data.table)
#library(ggrepel)
library("scales")
## iterate over all chromosome file (crm is chromosome)
if (dataType == "median") {
        datPath=paste(datadir,crm,'_median.txt.gz', sep='') #D
} else if (dataType == "hypergeom") {
        datPath=paste(datadir,crm,'/hypergeometric.txt.gz', sep='') # ED
} else {
        stop ("No valid dataType (median or hypergeom) specified")
}

print(datPath)
print(dataType)
## Read file of the chromosome
print(paste('working on chromosome: ',crm, sep='') )
dt = fread(paste(datPath,sep=''))
print(paste('number of rows: ',nrow(dt), sep='')) # print number of rows

## depending on which file you are reading, the format of the file will differ. There are two types of file
# 1- median file: it has only the chromosome position (BP) and median depth over samples for that base.
if (dataType=='median'){
        colnames(dt)=c('BP','D') # set column names as it does not have any
        dt=subset(dt, D>0) # remove the zeros depths as they will not be shown in the figure anyway
        max_pos=max(dt$BP) # max position to be set as x_lim
        print(paste('max_pos=',max_pos))
}

# 2- hypergeometric file: it has only the chromosome position (BP) and median depth over samples for that base,
# also it has the results of the hypergeometric probability as follows: expected depth if we sequence wit a lower depth (ED)',
# probabilty of that depth ('PED),probabilty of calling a variant withing the expected depth ('Pck')
if (dataType=='hypergeom'){
        colnames(dt)=c('CHR','BP','D','ED','PED','Pck') # set column names as it does not have any
        max_pos=max(dt$BP) # max position to be set as x_lim
        dt=subset(dt, D>0) # remove the zeros depths as they will not be shown in the figure anyway
        dt=dt[,c('BP','D','ED','PED','Pck')] # pick only the position, depth and expected depth to vistualize
}
####

print(paste('Number of reads : ',sum(dt$D), sep='')) # print number of rows
print(paste('Number of distinct reads : ',sum(unique(dt$D)), sep='')) # print number of rows

if (dataType=='hypergeom'){
        print(paste('Number of reads : ',sum(dt$D), sep='')) # print number of rows
        print(paste('Number of distinct reads : ',sum(unique(dt$D)), sep='')) # print number of rows
        print(paste('Number of predicted reads : ',sum(dt$ED), sep='')) # print number of rows
        print(paste('Number of predicted distinct reads : ',sum(unique(dt$ED)), sep='')) # print number of rows
}

if (dataType=='median'){
        spots=subset(dt,BP %in% positions)
        print (spots)
#        if (length(spots$D) == 0){
#                stop (paste ('No D value for BP:'))
#        }
        spots$info=paste(spots$BP,'(',spots$D,')',sep='')

        ##log10
        dt$D=log10(dt$D+1)
        spots$D=log10(spots$D+1)
        #print (paste('spotsdD = ',spots$D))

        png(paste('manhattenPlot_chrom_4Figs_1_',crm,'.png',sep='')) # save the resuls in a pdf
        p <-ggplot(dt, aes(BP,D, color = D, alpha = 0.4),label=info)+
                geom_point(shape = 16, size = 1, show.legend = FALSE)+geom_point(data=spots, aes(BP,D, color = D, alpha = 0.4),color='red',alpha =1, size=2) +
                theme_minimal() + scale_color_gradient(low = "#32aeff", high = "#f2aeff") +
#                scale_alpha(range = c(.25, .6))+geom_text_repel(data=spots,aes(BP,D, color = Pty,label=info),hjust=0, vjust=0,color='black',alpha =2.7,  size=30,angle = 0,max.overlaps=100)+
                scale_alpha(range = c(.25, .6))+
                        scale_x_continuous(name="Base position", labels = comma,limits = c(0, max_pos))+  scale_y_continuous(name="Log10 (depth)")
                print(p +ggtitle(paste('Chromosome: ',crm, sep=''))+theme(plot.title = element_text(hjust = 0.5))+ geom_hline(yintercept=log10(20), linetype="dashed", color = "green")+ geom_hline(yintercept=log10(5), linetype="dashed", color = "green"))
                dev.off() # save the pdf file
}

if (dataType=='hypergeom'){
        spots=subset(dt,BP %in% positions)
        dt=subset(dt, ED>0) # remove the zero rows
        dt$ED=log10(dt$ED+1)
        dt$D=log10(dt$D+1)
        spots$D=log10(spots$D+1)
        spots$ED=log10(spots$ED+1)

#        if (length(spots$ED) == 0){
#                stop (paste ('No ED value for BP:', position))
#       }
        spots$info=paste(spots$BP,'(',spots$ED,')',sep='')

        ## The user can specify postions on that chromosome to plot
# median only plot
        png(paste('manhattenPlot_median_',crm,'.png',sep='')) # save the resuls in a pdf
        p <-ggplot(dt, aes(BP,D, color = D, alpha = 0.4),label=info)+
                geom_point(shape = 16, size = 1, show.legend = FALSE)+geom_point(data=spots, aes(BP,D, color = D, alpha = 0.4),color='red',alpha =1, size=2) +
                theme_minimal() + scale_color_gradient(low = "#32aeff", high = "#f2aeff") +
#                scale_alpha(range = c(.25, .6))+geom_text_repel(data=spots,aes(BP,D, color = Pty,label=info),hjust=0, vjust=0,color='black',alpha =2.7,  size=30,angle = 0,max.overlaps=100)+
                scale_alpha(range = c(.25, .6))+
                        scale_x_continuous(name="Base position", labels = comma,limits = c(0, max_pos))+  scale_y_continuous(name="Log10 (depth)")
                print(p +ggtitle(paste('Chromosome: ',crm, sep=''))+theme(plot.title = element_text(hjust = 0.5))+ geom_hline(yintercept=log10(20), linetype="dashed", color = "green")+ geom_hline(yintercept=log10(5), linetype="dashed", color = "green"))
                dev.off() # save the pdf file

# hypergeom only plot
        png(paste('manhattenPlot_hypergeom_',crm,'.png',sep='')) # save the resuls in a pdf
        # plot
        p <-ggplot(dt, aes(BP,ED, color = ED, alpha = 0.4),label=BP)+
                geom_point(shape = 16, size = 1, show.legend = FALSE)+geom_point(data=spots, aes(BP,ED, color = ED, alpha =0.4),color='red',alpha =1, size=2)+
                theme_minimal() + scale_color_gradient(low = "#32aeff", high = "#f2aeff") +
                scale_alpha(range = c(.25, .6))+
                        scale_x_continuous(name="Base position", labels = comma,limits = c(0, max_pos))+  scale_y_continuous(name="log10(Expected Depth)")
                print(p +ggtitle(paste('Chromosome: ',crm, sep=''))+theme(plot.title = element_text(hjust = 0.5))+ geom_hline(yintercept=log10(20), linetype="dashed", color = "green")+ geom_hline(yintercept=log10(5), linetype="dashed", color = "green"))

                dev.off() # save the pdf file

}
