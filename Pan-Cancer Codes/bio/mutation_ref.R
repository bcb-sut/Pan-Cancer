

mutation_ref <- function() {
  ref_df = data.frame(matrix(ncol = 2))
  colnames(ref_df) = c('motif_3mer', 'motif_class')
  mutation_base <- c('CA', 'CG', 'CT', 'TA', 'TC', 'TG')
  for (i in 1:6) {
    mutation <- mutation_base[i]
    for (left_base in c('A', 'C', 'G', 'T')){
      for (right_base in c('A', 'C', 'G', 'T')){
        ref_df = rbind(ref_df, c(motif_3mer = paste0(mutation, '-', left_base, '.', right_base), motif_class = i))
      }
    }
  }
  ref_df = ref_df %>% drop_na()
}
