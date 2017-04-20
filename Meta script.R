# Kinship table

library(kinship2)
library(MCMCglmm)
library(pedigreemm)

d <- as.matrix(read.csv("g:/projects/meta-analysis/kinship.csv", header = TRUE))
d <- as.data.frame(d)
head(d)
str(d)


ped <- pedigree(d$Sire, d$Dam, d$Subject)


#transforming Raw scores into single column

d <- as.matrix(read.table("g:/projects/meta-analysis/collapsable.txt", header = FALSE), header = FALSE, sep="\t")
d<-matrix(d, ncol=1)
write.table(d, "g:/projects/meta-analysis/OneColumn.txt", sep="\t")
str(d)

#transforming Adjusted scores into single column
d <- as.matrix(read.table("g:/projects/meta-analysis/Adjcollapsable.txt", header = FALSE), header = FALSE, sep="\t")
d<-matrix(d, ncol=1)
write.table(d, "g:/projects/meta-analysis/OneAdjustedColumn.txt", sep="\t")

str(d)

#data set filling context + task columns

d <- as.matrix(read.table("g:/projects/meta-analysis/MetaData.txt", header = TRUE),  sep="\t", na.strings=c("#VALUE!"))
d <- as.data.frame(d)
head(d)

Tasks<- cbind("Diet","Tool", "Tool", "Tool","Tool", "Token", "Tool","Tool",
              "Tool",  NA,"Twoact", "Twoact", "Twoact","Tool","Tool","Tool",
               "Tool","Box", "Box", "Box","Box","Box","Box","Box","Box","Diet",
               "Diet", "Tool", "Box", "Box","Box","Box","Box",NA, "Box","Box")
str(Tasks)
d$Task <- rep(Tasks, each = 186)

Contexts<- cbind("Video", "Video", "Video","Video","Video","Group", "Ghost","Ghost", "Group",
                 "Video", "Dyad", "Ghost", "Ghost", "Group", "Group", "Group", "Group", "Group",
                 "Group","Group","Group","Group","Group", "Dyad", "Group","Group","Group","Group",
                 "Group","Group","Group","Group","Group", NA, "Group","Group")
str(Contexts)
d$Context <- rep(Contexts, each = 186)
write.table(d, "g:/projects/meta-analysis/Task_Context_output.txt", sep="\t")
