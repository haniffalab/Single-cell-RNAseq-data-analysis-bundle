args = commandArgs(trailingOnly=T)
args = paste(args, collapse = "")
args = unlist(strsplit(args, ";"))

arguments.list = "
seurat.addr.arg = args[1]
set.ident.arg   = args[2]
"

expected_arguments = unlist(strsplit(arguments.list, "\n"))
expected_arguments = expected_arguments[!(expected_arguments == "")]

if(length(args) != length(expected_arguments)){
  error.msg = sprintf('This pipeline requires %s parameters', as.character(length(expected_arguments)))
  expected_arguments = paste(unlist(lapply(strsplit(expected_arguments, ".arg"), "[", 1)), collapse = "\n")
  stop(sprintf('This pipeline requires %s parameters: '))
}

eval(parse(text = arguments.list))

for(n in 1:length(expected_arguments)){
  argument = expected_arguments[n]
  argument = gsub(pattern=" ", replacement="", x=argument)
  argument.name = unlist(strsplit(argument, "="))[1]
  variable.name = gsub(pattern=".arg", replacement="", argument.name)
  argument.content = eval(parse(text = argument.name))
  eval(parse(text = argument.content))
  if (!exists(variable.name)){
    stop(sprintf("Argument %s not passed. Stopping ... ", variable.name))
  }
}

# create required folders for output and work material
output_folder = gsub(pattern="^\\d+_", replacement="", x=basename(getwd()))
output_folder = paste(output_folder, seurat.addr, sep = "_")
c.time = Sys.time()
c.time = gsub(pattern=" BST", replacement="", x=c.time)
c.time = gsub(pattern=":", replacement="", x=c.time)
c.time = gsub(pattern=" ", replacement="", x=c.time)
c.time = gsub(pattern="-", replacement="", x=c.time)
c.time = substr(x=c.time, start=3, stop=nchar(c.time))
output_folder = paste(output_folder, c.time, sep = "_")
output_folder = file.path("../../output", output_folder)
dir.create(output_folder)

passed.threshold         = .5
competition.threshold    = .8
expansion                = 4

library(Seurat)
library(plyr)
library(dplyr)

source("../../tools/bunddle_utils.R")

seurat.addr = file.path("../../data", seurat.addr)

gene_mean_expression = function(seurat.obj){
  no.genes = nrow(seurat.obj@data)
  start_index = 1
  while (start_index < no.genes){
    end_index = start_index + 999
    end_index = min(end_index, no.genes)
    expression.data_ = data.matrix(seurat.obj@data[start_index:end_index, ])
    expression.data_ = t(expression.data_)
    expression.data_ = as.data.frame(expression.data_)
    expression.data_ = cbind(data.frame(CellLabels = as.vector(seurat.obj@ident)), expression.data_)
    expression.data_ = aggregate(expression.data_[2:dim(expression.data_)[2]], list(expression.data_$CellLabels), mean)
    expression.data_ = cbind(data.frame(CellType = expression.data_$Group.1), expression.data_[, 2:dim(expression.data_)[2]])
    rownames(expression.data_) = expression.data_$CellType
    expression.data_ = expression.data_[, 2:ncol(expression.data_)]
    print(start_index)
    if (start_index == 1){
      expression.data = expression.data_
    }else{
      expression.data = cbind(expression.data, expression.data_)
    }
    start_index = start_index + 1000
  }
  expression.data
}

print("Loading data ... ")
seurat.obj = readRDS(seurat.addr)
seurat.obj = SetAllIdent(seurat.obj, id=set.ident)
print("Data loaded.")

# get aggregated expression matrix
expression.data = gene_mean_expression(seurat.obj)

gene.names = colnames(expression.data)
cell.types = rownames(expression.data)
to.remove = grep(pattern="^MT-", gene.names)
to.remove = c(to.remove, grep(pattern="^RPS", gene.names))
gene.names = gene.names[-to.remove] # remove mito and ribo genes
expression.data = expression.data[, gene.names]
n.cell.types = length(cell.types)
super.markers = rep("", length(cell.types))
names(super.markers) = cell.types

gene.max.expression  = apply(expression.data, 2, max)
expression.data.norm = t(t(expression.data) / gene.max.expression)

for(j in seq_along(gene.names)){
  gene.name = gene.names[j]
  super.vector = expression.data.norm[, gene.name]
  super.order  = order(super.vector, decreasing=T)
  cell.type    = names(super.vector[super.order[1]])
  expr.val     = as.vector(expression.data[cell.type, gene.name] )
  competition = all(super.vector[super.order[expansion:length(super.order)]] < competition.threshold)
  if(competition & expr.val > passed.threshold){
    super.markers[cell.type] = paste(super.markers[cell.type], gene.name, sep = ", ")
  }
}

# order markers
for (i in seq_along(super.markers)){
  ms        = super.markers[i]
  cell.type = names(ms)[1]
  ms = unlist(strsplit(as.character(ms), split=", "))
  ms = ms[ms != ""]
  expression.row = expression.data[cell.type, ms]
  ms = ms[rev(order(expression.row))]
  ms = paste(ms, collapse = ", ")
  super.markers[i] = ms
}

# make classification markers
classification.markers = super.markers

seurat.obj = FindVariableGenes(object = seurat.obj, mean.function = ExpMean, 
                               dispersion.function = LogVMR, x.low.cutoff = .0125, 
                               x.high.cutoff = 9, y.cutoff = .2, do.plot=F)
seurat.obj@var.genes = seurat.obj@var.genes[-c(grep("^MT-", seurat.obj@var.genes), grep("^RPS", seurat.obj@var.genes))]

for (k in seq_along(cell.types)){
  silent.variable = tryCatch({
    cell.type = cell.types[k]
    eval(parse(text = sprintf("cells.markers = super.markers[['%s']]", cell.type)))
    cells.markers = unlist(strsplit(cells.markers, split=", "))
    cor.expression.data = expression.data[, cells.markers]
    compare.to = cell.types[cell.types != cell.type]
    cor.dist = cor(t(cor.expression.data), method="spearman")[cell.type, compare.to]
    compare.to = names(cor.dist)[cor.dist > .3]
    if(length(compare.to) > .45){
      cor.dist = cor.dist[compare.to]
      cor.dist = cor.dist[order(cor.dist, decreasing=T)]
      compare.to = names(cor.dist)[1:3]
    }
    if (length(compare.to) > 0 & !any(is.na(compare.to))){
      print(sprintf("Comparing %s to: %s", cell.type, paste(compare.to, collapse = ", ")))
      DEGs = list()
      for (m in seq_along(compare.to)){
        DEG = FindMarkers(object=seurat.obj, ident.2=compare.to[m], ident.1=cell.type,max.cells.per.ident=500, only.pos=T, min.pct=.5, genes.use=seurat.obj@var.genes, )
        DEG$cluster = compare.to[m]
        DEG$gene = rownames(DEG)
        DEGs[[m]] = DEG
      }
      DEGs = Reduce(f=rbind, x=DEGs)
      additional.markers = (DEGs %>% group_by(cluster) %>% top_n(5, avg_logFC))$gene
      new.markers = unique(c(cells.markers, additional.markers))
      new.markers = paste(new.markers, collapse = ", ")
      eval(parse(text = sprintf("classification.markers[['%s']] = new.markers", cell.type)))
    }
  },
  warning = function(warn_cond){print("")},
  error   = function(err_cond){print("")})
}

# order markers
for (i in seq_along(classification.markers)){
  ms        = classification.markers[i]
  cell.type = names(ms)[1]
  ms = unlist(strsplit(as.character(ms), split=", "))
  ms = ms[ms != ""]
  expression.row = expression.data[cell.type, ms]
  ms = ms[rev(order(expression.row))]
  ms = paste(ms, collapse = ", ")
  classification.markers[i] = ms
}

saveRDS(classification.markers, file.path(output_folder,"signatures.RDS"))

sig.df = data.frame(CellTypes = names(classification.markers), Signatures = classification.markers)
write.csv(sig.df, file.path(output_folder,"Signatures.csv"), row.names = F)

print("Ended beautifully ... ")

