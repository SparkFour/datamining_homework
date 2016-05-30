#lib
library(car)
library(DMwR)
library(lattice)

#change work dir
script.dir <- dirname(sys.frame(1)$ofile)
setwd(script.dir)

##################### const variable #####################
original_filename <- "Analysis.txt"
col_names <- c("season", "size", "speed", "mxPH", "mnO2", "Cl", "NO3", "NH4", "oPO4", "PO4", "Chla", "a1", "a2", "a3", "a4", "a5", "a6", "a7")
na_strings <- c("XXXXXXX", "NA")
strategy_miss_attr <- c("剔除", "高频属性代替", "属性之间的相关性", "对象之间相似性")
name_of_strategy <- c("delete", "most", "correlation", "similarity")


#获取数据集的摘要和属性值可视化信息
get_summary_visualinfo <- function(filename, path="."){
  #check file
  if (file.exists(filename)==FALSE) {
    msg <- paste("Data file doesn't exist: ", filename, sep = "")
    stop(msg)
  }
  if (file.exists(path)==FALSE){
    dir.create(path)
  }
  #read
  algae <- read.table(filename, col.names = col_names, na.strings = na_strings)
  
  #summary
  data_summary <- summary(algae)
  #print(data_summary)
  summary_filename <- paste(path, "/", "summary.csv", sep="")
  write.csv(data_summary, file = summary_filename)
  jpeg(paste(path, "/summary.jpg", sep=""))
  plot(algae)
  dev.off()
  #visual
  #numeric attribute
  p_names <- col_names[4:18]
  p_data <- algae[, 4:18]
  for (i in 1:1){
    value <- p_data[,i]
    #hist
    fig_hist_filename <- paste(path, "/hist_", p_names[i], ".jpg", sep="")
    jpeg(fig_hist_filename)
    hist(value, probability = T, main = paste("histogram of ", p_names[i], sep=""))
    rug(value, side = 1)
    dev.off()
    #qqPlot
    fig_qq_filename <- paste(path, "/QQ_", p_names[i], ".jpg", sep="")
    jpeg(fig_qq_filename)
    qqPlot(value)
    dev.off()
    #boxplot
    fig_box_filename <- paste(path, "/box_", p_names[i], ".jpg", sep="")
    jpeg(fig_box_filename)
    boxplot(value)
    rug(value, side=4)
    abline(h=mean(value, na.rm=T), lty=2)
    dev.off()
  }
  #non-numeric attribute
  for (i in 1:3){
    #bwplot
    fig_bw_filename <- paste(path, "/bw_", p_names[i], ".jpg", sep="")
    jpeg(fig_bw_filename)
    bwplot(algae$size~a1, data=algae, xlab = "xx", ylab = "xx")
    dev.off()
    print(fig_bw_filename)
  }
  return()
}

#取得高频元素(数字或字符串),可设置去掉缺失值
get_most_appear_value <- function(x, na_rm=F, na_string="NA") {
  #
  if (na_rm==T){
    x <- x[x!=na_string]
  }
  #
  if (mode(x)=='numeric') {
    return(as.numeric(names(table(x)))[which.max(table(x))])
  }
  else {
    return(names(table(x))[which.max(table(x))])
  }
}


#四种方案，处理缺失属性
solve_miss_attr <- function(filename){
  #check
  if (file.exists(filename)==FALSE) {
    msg <- paste("Data file doesn't exist: ", filename, sep = "")
    stop(msg)
  }
  #read
  algae <- read.table(filename, col.names = col_names, na.strings = na_strings)
  #strategy 1
  algae_1 <- na.omit(algae)
  write.table(algae_1, paste(name_of_strategy[1], "_data.txt", sep=""), quote = F, col.names = F, row.names = F)
  #strategy 2
  algae_2 <- algae
  for (i in 1:18) {
    algae_2[is.na(algae_2[i]), i] <- get_most_appear_value(algae_2[i], T, na_strings)
  }
  write.table(algae_2, paste(name_of_strategy[2], "_data.txt", sep=""), quote = F, col.names = F, row.names = F)
  #strategy 3
  #通过如下指令可视化，发现变量之间的相关性,发现PO4与oPO4具有相关性
  symnum(cor(algae[,4:18], use = "complete.obs"))
  #去除缺失太多属性的值
  algae_3 <- algae[-manyNAs(algae, 0.2), ]
  #建立线性模型，可输出查看，PO4 = 1.29 * oPO4 + 42.9
  xx = lm(PO4~oPO4, data=algae)
  #print(xx)
  #
  algae_3[is.na(algae_3["PO4"]), "PO4"] <- 1.29 * algae_3[is.na(algae_3["PO4"]), "oPO4"] + 42.9
  write.table(algae_3, paste(name_of_strategy[3], "_data.txt", sep=""), quote = F, col.names = F, row.names = F)
  #strategy 4
  algae_4 <- algae[-manyNAs(algae, 0.2), ]
  #使用最近的k个样本
  algae_4 <- knnImputation(algae_4, k=1)
  write.table(algae_4, paste(name_of_strategy[4], "_data.txt", sep=""), quote = F, col.names = F, row.names = F, sep = "\t")
  return()
}


################# start ########################


get_summary_visualinfo(original_filename, path = "original")
solve_miss_attr(original_filename)
for (i in 1:4) {
  get_summary_visualinfo(paste(name_of_strategy[i], "_data.txt", sep=""), path=name_of_strategy[i])
}

#条件盒图上面的代码保存不了图片，原因还未查清，使用下面代码在控制台中执行
setwd('~/Study/研究生课程/数据挖掘/作业/1/')
data <- read.table("Analysis.txt", col.names=c("season", "size", "speed", "mxPH", "mnO2", "Cl", "NO3", "NH4", "oPO4", "PO4", "Chla", "a1", "a2", "a3", "a4", "a5", "a6", "a7"), na.strings=c("XXXXXXX", "NA"))
bwplot(data$season~a1, data=data, xlab='a1', ylab='season', add=T)
bwplot(data$season~a2, data=data, xlab='a2', ylab='season', add=T)
bwplot(data$season~a3, data=data, xlab='a3', ylab='season', add=T)
bwplot(data$season~a4, data=data, xlab='a4', ylab='season', add=T)
bwplot(data$season~a5, data=data, xlab='a5', ylab='season', add=T)
bwplot(data$season~a6, data=data, xlab='a6', ylab='season', add=T)
bwplot(data$season~a7, data=data, xlab='a7', ylab='season', add=T)

bwplot(data$size~a1, data=data, xlab='a1', ylab='size', add=T)
bwplot(data$size~a2, data=data, xlab='a2', ylab='size', add=T)
bwplot(data$size~a3, data=data, xlab='a3', ylab='size', add=T)
bwplot(data$size~a4, data=data, xlab='a4', ylab='size', add=T)
bwplot(data$size~a5, data=data, xlab='a5', ylab='size', add=T)
bwplot(data$size~a6, data=data, xlab='a6', ylab='size', add=T)
bwplot(data$size~a7, data=data, xlab='a7', ylab='size', add=T)

bwplot(data$speed~a1, data=data, xlab='a1', ylab='speed', add=T)
bwplot(data$speed~a2, data=data, xlab='a2', ylab='speed', add=T)
bwplot(data$speed~a3, data=data, xlab='a3', ylab='speed', add=T)
bwplot(data$speed~a4, data=data, xlab='a4', ylab='speed', add=T)
bwplot(data$speed~a5, data=data, xlab='a5', ylab='speed', add=T)
bwplot(data$speed~a6, data=data, xlab='a6', ylab='speed', add=T)
bwplot(data$speed~a7, data=data, xlab='a7', ylab='speed', add=T)






