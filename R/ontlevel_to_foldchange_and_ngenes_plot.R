#' Ontology level to fold change & ontology level to n genes facet plot
#'
#' This creates two plots faceted which show the relationship between ontology
#' level and fold change in specific expression and number of associated genes
#' in the rare disease EWCE results and HPO. The phenotype name column is called
#' list.
#' @param all_results_merged HPO EWCE results <data.frame>
#' @param hpo The HPO data object from ontologyIndex
#' @param phenotype_to_genes phenotype gene list annotations from hpo
#' phenotype_to_genes.txt
#' @import dplyr
#' @phenotype_to_genes The gene list annotations from HPO phenotype_to_genes.txt
#' <data.frame>
#' @export
ontlevel_to_foldchange_and_ngenes_plot <- function(all_results_merged,
                                               hpo,
                                               phenotype_to_genes){
  ont_level_df = data.frame()
  for (p in unique(all_results_merged$list)) {
    #hpo_id = paste(hpo$id[match(p,hpo$name)])
    hpo_id = get_hpo_termID(p,phenotype_to_genes)
    if (hpo_id != "NA") {
      mean_fold = mean(all_results_merged[all_results_merged$list == p & all_results_merged$q < 0.05, ]$fold_change, na.rm=TRUE)
      n_genes = length(unique(phenotype_to_genes$Gene[phenotype_to_genes$Phenotype==p]))
      #hpo_id = paste(hpo$id[match(p,hpo$name)])
      hpo_id = get_hpo_termID(p,phenotype_to_genes)
      ont_level = get_ont_level(hpo, hpo_id)
      ont_level_df = rbind(ont_level_df,
                           data.frame("phenotype"=p,
                                      "hpo_id"=hpo_id,
                                      "ont_level"=ont_level,
                                      "n_genes"=n_genes,
                                      "mean_fold"=mean_fold))
    }
  }
  # plot
  ont_level_df2 = ont_level_df[ont_level_df$mean_fold != "NaN",]
  ont_level_df2 = ont_level_df2 %>%
    group_by(ont_level) %>%
    summarize(mean_fold_change = mean(mean_fold, na.rm=TRUE))

  lvlplt1 = ggplot(ont_level_df2,aes(x=ont_level,y=mean_fold_change)) +
    geom_point() +
    geom_line() +
    cowplot::theme_cowplot() +
    labs(y = "Mean fold change", x = "Ontology level")

  ont_level_df2 = ont_level_df[ont_level_df$mean_fold != "NaN",]
  ont_level_df2 = ont_level_df2 %>%
    group_by(ont_level) %>%
    summarize(mean_n_genes = mean(n_genes, na.rm=TRUE))

  lvlplt2 = ggplot(ont_level_df2,aes(x=ont_level,y=mean_n_genes)) +
    geom_point() +
    geom_line() +
    cowplot::theme_cowplot() +
    scale_x_continuous() +
    labs(y = "Mean number of genes", x = "Ontology level")

  return(cowplot::plot_grid(plotlist=list(lvlplt1,lvlplt2),align = "hv",nrow = 1,labels = c("A","B")))
}

