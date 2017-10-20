

data_folder <- 'data_raw/'
data_set = "GSM2396857_dc_0hr"

counts.dt <- fread(paste0('gzip -dc ', data_folder, data_set,'.mtx.txt.gz'))
setnames(counts.dt,c("gene_id","cell_id","count"))
n_genes <- counts.dt[1,gene_id]
n_cells <- counts.dt[1,cell_id]
n_not_no_reads <- counts.dt[1,count]
counts.dt <- counts.dt[-1]

# Load ID-cell and ID-gene mapping
cellnames.dt <- fread(paste0('gzip -dc ', data_folder, data_set, '_cellnames.csv.gz'), header=TRUE)
setnames(cellnames.dt, c("cell_id", "cell"))
cellnames.dt[,cell_id:=cell_id+1] # changing to one_based

cellnames.dt[,batch:=stringr::str_match(cell,".*_(.*_.*)")[,2]]
cellnames.dt[,cell_bc:=stringr::str_match(cell,"(.*)_.*_.*")[,2]]

counts.dt <- counts.dt[cellnames.dt, on=.(cell_id)]


ddt = counts.dt[, .(total_counts = sum(count)), by=.(cell_id, batch)]
ddt

figure("batch effect on total_counts",
ggplot(ddt, aes(x=total_counts,color=batch))+geom_density()+scale_x_log10()+
  geom_jitter_normal(aes(y=-as.integer(as.factor(batch))/20.0),height=0.01,width=0,size=1)
)


dt = counts.dt[, .(mean_counts = mean(count)), by=.(gene_id, batch)]
dt = dcast(dt, gene_id~batch, value.var="mean_counts",fill=0)
ggplot(dt, aes(y=dc0h_G9, x=dc0h_E8))+geom_point()+scale_x_log10()+scale_y_log10()+geom_abline(slope=dt[,sum(dc0h_G9)/sum(dc0h_E8)])+geom_smooth()
