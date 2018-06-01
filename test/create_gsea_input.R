anno <- read.delim(file.choose(), header=T, skip=1, row.names=1)
u_anno <- as.character(unique(anno$subclass))
for (row in u_anno){
  temp_subset <- subset(anno, anno$subclass==row)
  print(row)
  print(nrow(temp_subset))
}

# create subset to get all genes matching annotation subclass
# 