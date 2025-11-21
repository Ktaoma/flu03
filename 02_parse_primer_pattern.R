library(data.table)
library(dplyr)


args <- commandArgs(trailingOnly = TRUE)

F_ <- fread(paste0(args[1],"/forward.tsv.gz"))
R_ <- fread(paste0(args[1],"/reverse.tsv.gz"))

F_seq <- table(F_$seqID,F_$strand) |> as.data.frame() |> reshape2::dcast(Var1~Var2,value.var = "Freq") 
head(F_seq)
colnames(F_seq) <- c("seq","F-","F+")

R_seq <- table(R_$seqID,R_$strand) |> as.data.frame() |> reshape2::dcast(Var1~Var2,value.var = "Freq") 
head(R_seq)
colnames(R_seq) <- c("seq","R-","R+")

dat <- left_join(F_seq,R_seq,by = "seq") 
dat[is.na(dat)] <- 0
dat_final <- paste0(dat$`F+`,"_",dat$`R-`,"_",dat$`R+`,"_",dat$`F-`) |> table() |> as.data.frame() 


if(args[2] == 'both'){

    set1 <- dat[dat$`F+` == 1 & dat$`R-` == 1 & dat$`R+` == 0 & dat$`F-` == 0,1] %>% as.data.frame()
    set2 <- dat[dat$`F+` == 0 & dat$`R-` == 0 & dat$`R+` == 1 & dat$`F-` == 1,1] %>% as.data.frame()
    rbind(set1,set2) |> write.table(paste0(args[1],"/QC_id.tsv"),quote = FALSE, row.names = FALSE, col.names = FALSE)

} else if (args[2] == 'both_F'){
   
    set1 <- dat[dat$`F+` == 1 & dat$`R-` == 1 & dat$`R+` == 0 & dat$`F-` == 0,1] %>% as.data.frame()
    set2 <- dat[dat$`F+` == 0 & dat$`R-` == 0 & dat$`R+` == 1 & dat$`F-` == 1,1] %>% as.data.frame()
    set3 <- dat[dat$`F+` == 1 & dat$`R-` == 0 & dat$`R+` == 0 & dat$`F-` == 0,1] %>% as.data.frame()
    rbind(set1,set2,set3) |> write.table(paste0(args[1],"/QC_id.tsv"),quote = FALSE, row.names = FALSE, col.names = FALSE)

} else {

    set1 <- dat[dat$`F+` == 1 & dat$`R-` == 1 & dat$`R+` == 0 & dat$`F-` == 0,1] %>% as.data.frame()
    set2 <- dat[dat$`F+` == 0 & dat$`R-` == 0 & dat$`R+` == 1 & dat$`F-` == 1,1] %>% as.data.frame()
    set3 <- dat[dat$`F+` == 1 & dat$`R-` == 0 & dat$`R+` == 0 & dat$`F-` == 0,1] %>% as.data.frame()
    set4 <- dat[dat$`F+` == 0 & dat$`R-` == 1 & dat$`R+` == 0 & dat$`F-` == 0,1] %>% as.data.frame()
    rbind(set1,set2,set3) |> write.table(paste0(args[1],"/QC_id.tsv"),quote = FALSE, row.names = FALSE, col.names = FALSE)

}

