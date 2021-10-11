library(readr)
library(plyr)

# motif -------------------------------------------------------------------

dirs = list.files(path = "matrices/3mermotif/", full.names = TRUE, recursive = FALSE)

final = data.frame()
final_WXS = data.frame()
final_WGS = data.frame()

for (dir in dirs) {
  print(dir)
  a = read_tsv(dir)
  if (grepl("WGS", dir))
    final_WGS = rbind.fill(final_WGS, a)
  else if (grepl("WXS", dir))
    final_WXS = rbind.fill(final_WXS, a)
  else
    final = rbind.fill(final, a)
  
  # print(paste(nrow(final), "   WXS:", nrow(final_WXS), "   WGS:", nrow(final_WGS)))
}

final[is.na(final)] <- 0
final_WGS[is.na(final_WGS)] <- 0
final_WXS[is.na(final_WXS)] <- 0

write_tsv(final, file.path("matrices/", "3mermotif.tsv.gz"))
write_tsv(final_WGS, file.path("matrices/", "3mermotif_WGS.tsv.gz"))
write_tsv(final_WXS, file.path("matrices/", "3mermotif_WXS.tsv.gz"))




# basic -------------------------------------------------------------------

dirs = list.files(path = "matrices/allgenemotifs/", full.names = TRUE, recursive = FALSE)

final = data.frame()

for (dir in dirs) {
  print(dir)
  a = read_tsv(dir)
  final = rbind.fill(final, a)
}

final[is.na(final)] <- 0

write_tsv(final, file.path("matrices/", "allgenemotifs.tsv.gz"))

