library(crumblr)
library(variancePartition)

for (region in c("C")){
  path_counts = list.files(path='/home/MinaRyten/Melissa/teamWood/Abundance_final_data/processed_data/counts', pattern=paste0('*.csv'), all.files=FALSE, full.names=TRUE)
  for (n in 1:length(path_counts)) {
    print(n)
    print(path_counts[n])
    counts = read.csv(path_counts[n], row.names = 'sample_id')
    counts[is.na(counts)] <- 0
    split = strsplit(path_counts[n], 'counts_')
    covariates = read.csv(paste0('/home/MinaRyten/Melissa/teamWood/Abundance_final_data/processed_data/covariates/covariates_', split[[1]][2]), row.names = 'sample_id')
    cobj = crumblr(counts)
    form = ~ group + (1|sex) + (1|batch) + scale(total_deduplicated_percentage)
    print(split[[1]][2])
    print(form)
    fit = dream(cobj, form, covariates)
    print(attr(fit, 'errors'))
    fit = eBayes(fit)
    res = topTable(fit, coef="groupPD", number=Inf, sort.by="none")
    write.csv(res, file = paste0('/home/MinaRyten/Melissa/teamWood/Abundance_final_data/processed_data/results/results_crumblr_', split[[1]][2]))
  }
}

warnings()