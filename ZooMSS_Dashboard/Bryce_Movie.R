N <- data$model$N

dt <- data$model$param$dt*data$model$param$isave

df <- data.frame()
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
for (class in 1:dim(N)[2]){
  y_height = max(N[,class,])
  for (i in 1:dim(N)[1]){
    df <- data.frame(size_class=1:dim(N)[3], Ns = N[i,class,])
    g <- ggplot(df, aes(size_class, Ns))
    g + geom_col(fill=cbbPalette[1+class%%8], color=cbbPalette[1+class%%8]) + scale_y_continuous(limit = c(0, y_height*1.2)) + ggtitle(paste(c("Class", class)), subtitle=paste("Year:",dt*i))
    ff_format = toString(i)
    ff_format = paste(c(replicate(4-nchar(ff_format), "0"), ff_format), collapse="")
    ggsave(paste("images/class", toString(class), "save", ff_format, ".png", collapse="", sep=""))
    print((class*dim(N)[1] + i)/(dim(N)[2]*dim(N)[1]))
  }
}