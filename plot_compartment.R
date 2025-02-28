######### Compartment plots ##########
##################################################################
# Developer: Milad Mokhtaridoost, milad.mokhtaridoost@sickkids.ca
##################################################################

library(ggplot2)
library(tidyr)

data <- read.csv("HiC_Gillespie_compartments1000000.csv")
#data <- final_compartment ## final combined data

mydata <- data[which(is.na(data$E1) == FALSE),]

chrs <- c(1:22, "X", "Y")
for(chr in chrs){
mydata1 <- mydata[which(mydata$chrom == sprintf("chr%s",chr)),]

for(i in 1:nrow(mydata1)){
if(mydata1$E1[i] > 0){
  mydata1$compartment[i] <- "A" } else {
    mydata1$compartment[i] <- "B" 
  }
}

mydata1$compartment <- as.factor(mydata1$compartment)
mydata1$start <- mydata1$start/1000000

p <- mydata1 %>% 
  ggplot(., aes(x = start, y = E1)) +
  geom_line() +
  geom_point(aes(color = compartment)) + 
  geom_hline(yintercept=0, color = "red")+ 
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(mydata1$start))) 

ggsave(sprintf("average_compartment_1MB_chr%s.pdf", chr), p, width = 14, height = 6, dpi = 300)
}
