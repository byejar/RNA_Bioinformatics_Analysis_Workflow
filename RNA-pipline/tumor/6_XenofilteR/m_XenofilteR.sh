library("XenofilteR")
bp.param <- SnowParam(workers = 20, type = "SOCK")
m_sample.list <- read.table(file ='/home/wus/scRNA_program_blackmouseRNA/Nude/5_XenofilteR/m_sample.list',
                           header = F, sep = ",")
XenofilteR(m_sample.list, destination.folder = "/home/wus/scRNA_program_blackmouseRNA/Nude/5_XenofilteR/mouse", bp.param = bp.param)
