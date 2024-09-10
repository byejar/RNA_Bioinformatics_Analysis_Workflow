library("XenofilteR")
bp.param <- SnowParam(workers = 60, type = "SOCK")
m_sample.list <- read.table(file ='m_sample.list',
                           header = F, sep = ",")
XenofilteR(m_sample.list, destination.folder = "/home/wus/2023_3_19-009-cxl_RNA/TPL202303714+TPL202303545/pipline_2023-3-23_RNA-analysis/5_XenofilteR/mouse", bp.param = bp.param)
