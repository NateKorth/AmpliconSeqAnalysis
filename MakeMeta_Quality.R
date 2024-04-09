getwd()
#should be AmpliconTest

library("tidyverse")
library("dada2")
library("gridExtra")
no_of_cores=16

#Generate metadata if it doesn't exist
files<-list.files("./RawData")
R1<-files[grep("R1",files)]
R2<-files[grep("R2",files)]

AmpliconTest<-as.data.frame(cbind(R1,R2))

#Following steps are different by project, here we'll try our best to extract information from the file names:
AmpliconTest$ID <- sub("_S.*", "", AmpliconTest$R1)
AmpliconTest$ID2 <- sub("-.*", "", AmpliconTest$R1)
AmpliconTest$Day <- sub(".*-", "", AmpliconTest$R1)
AmpliconTest$Day <- sub("_S.*", "", AmpliconTest$Day)
AmpliconTest$Subject <- substr(AmpliconTest$ID2, 1, 3) # Extracts the first 3 characters
AmpliconTest$Treatment <- substr(AmpliconTest$ID2, 4, 5)
  
#Write a file for any additional customization in excel an usage in downstream analysis:  
#write.table(AmpliconTest,"AmpliconTest.csv",sep=",",row.names = FALSE,)

#this file should contain the names of forward and reverse reads:
metadata_df<-read.csv("AmpliconTest.csv")

#Add relative file path to the reads for the next step: 

#forward reads
fnFs<-paste("./RawData/",metadata_df$R1,sep="")
#reverse reads
fnRs<-paste("./RawData/",metadata_df$R2,sep="")


fPlot <- dada2::plotQualityProfile(fnFs, aggregate = TRUE) +
    ggtitle("Forward") +
        geom_hline(yintercept =  30, colour = "blue")
rPlot <- dada2::plotQualityProfile(fnRs, aggregate = TRUE) +
        ggtitle("Reverse") +
        geom_hline(yintercept =  30, colour = "blue")

g<-arrangeGrob(fPlot, rPlot, nrow = 1)
ggsave("AveQuality.pdf",plot=g)

#continue to ASVpick.R and adjust trimming values based 
