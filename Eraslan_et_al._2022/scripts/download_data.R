classic_df <- "https://storage.googleapis.com/adult-gtex/single-cell/v9/snrna-seq-data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad"
immune_df <- "https://storage.googleapis.com/adult-gtex/single-cell/v9/snrna-seq-data/GTEx_8_tissues_snRNAseq_immune_atlas_071421.public_obs.h5ad"

urls <- c(classic_df, immune_df)

for (url in urls) {
    download.file(url, destfile = basename(url), method="curl", extra="-k")
    }

  