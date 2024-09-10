library("XenofilteR")
bp.param <- SnowParam(workers = 20, type = "SOCK")
sample.list <- read.table(file ='./sample.list',
                           header = F, sep = ",")
XenofilteR(sample.list, destination.folder = "./", bp.param = bp.param)
