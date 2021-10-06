#' Proportion of expected enrichments in a cell type by significance plot
#'
#' This allows you to select a group of phenotypes by from a particular HPO
#' branch (or selection of phenotypes) and a cell type (or collection of cells)
#' that you expected is often associated with that phenotype.
#' The plot produced shows the proportion of enrichments in that cell type that
#' are from the expected branch at different levels of significance.
#' Typically as you go towards more significant
#' phenotypes, you would expecte a higher number of them to be from the
#' expected branch. E.g. we may see many "Abnormality of the nervous system"
#' sub-phenotypes significantly enriched in excitatory neurons.
#' @param all_results_merged Rare disease EWCE results <data.frame>
#' @param target_cells vector of target cell names <vector<string>>
#' @param cell_type_description a description of the target cells
#' (e.g. "Cells of the nervous system") <string>
#' @param HPO_Ids vector of HPO Ids of expected phenotypes <vector<string>>
#' @param phenotype_description a description of the expected phenotypes
#' (e.g. "Nervous system abnormalities") <string>
#' @param wes_color_palette a color palette from wesanderson package <string>
#' @param n_colors The number of colors to select from the palette
#' (see docs for wesanderson package) <int>
#' @param color_expected_phenotypes The index of the color from the selected
#' palette that you want to represent expected phenotypes <int>
#' @param color_other_phenotypes The index of the color from selected palette
#' that you want to represent other phenotypes <int>
#' @param hpo The HPO data object from ontologyIndex package
#' @import wesanderson
#' @import ggplot2
#' @examples \dontrun{
#' # Single plot using multiple related cell types,
#' # note "[[1]]" is to select the plot.
#' proportion_of_expected_enrichments_plot(all_results_merged,
#' target_cells=c("Excitatory neurons","Limbic system neurons"),
#' cell_type_description = "Neuronal cells",
#' HPO_Ids = ontologyIndex::get_descendants(
#'   hpo,hpo$id[match("Abnormality of the nervous system",hpo$name)]),
#' phenotype_description = "Nervous system phenotypes",
#' wes_color_palette="Darjeeling1",
#' n_colors = 4,
#' color_expected_phenotypes = 1,color_other_phenotypes = 4)[[1]]
#'
#' # Plot multiple facted plots and run pearsons correlation on pvalue~n expected phenos
#' plot_branches = c("Abnormality of the nervous system",
#'                   "Abnormality of the cardiovascular system",
#'                   "Abnormality of the immune system")
#' expected_cells = c("Excitatory neurons",
#'                    "Cardiomyocytes",
#'                    "Antigen presenting cells")
#' correlation_results = data.frame()
#' proportion_plots = list()
#' for (i in seq(length(plot_branches))) {
#'   cur_plot <- proportion_of_expected_enrichments_plot(all_results_merged,
#'                                                       target_cells=c(expected_cells[i]),
#'                                                       cell_type_description = expected_cells[i],
#'                                                       HPO_Ids = ontologyIndex::get_descendants(
#'                                                       hpo,hpo$id[match(plot_branches[i],hpo$name)]),
#'                                                       phenotype_description = plot_branches[i],
#'                                                       wes_color_palette="Darjeeling1",
#'                                                       n_colors = 4,
#'                                                       color_expected_phenotypes = i,
#'                                                       color_other_phenotypes = 4)
#'   proportion_plots[[i]] = cur_plot[[1]]
#'   correlation_results = rbind(correlation_results,cur_plot[[2]])
#' }
#' print(cowplot::plot_grid(
#' plotlist=proportion_plots, align = "h",nrow = 1,labels=c("A","B","C")))
#' pearsons_expected = cor.test(
#' correlation_results$minus_log_p[correlation_results$cur_branch != "Other"],
#' correlation_results$prop[correlation_results$cur_branch != "Other"],method="pearson")
#' }
#' @export
proportion_of_expected_enrichments_plot <- function(all_results_merged,
                                                    hpo,
                                                    target_cells=c("Excitatory neurons","Limbic system neurons"),
                                                    cell_type_description = "Neuronal cells",
                                                    HPO_Ids = ontologyIndex::get_descendants(hpo,hpo$id[match("Abnormality of the nervous system",hpo$name)]),
                                                    phenotype_description = "Nervous system phenotypes",
                                                    wes_color_palette="Darjeeling1",
                                                    n_colors = 4,
                                                    color_expected_phenotypes = 1,
                                                    color_other_phenotypes = 4
                                                    ) {
  cur_descendants_names = c()
  for(d in HPO_Ids) {
    cur_descendants_names = append(cur_descendants_names, paste(hpo$name[d]))
  }
  all_results_merged$cur_branch = "Other"
  all_results_merged[all_results_merged$list %in% cur_descendants_names,]$cur_branch = phenotype_description
  all_results_merged$cur_branch = factor(all_results_merged$cur_branch)
  all_results_merged$minus_log_p = round(-log10(all_results_merged$p+0.000001),digits=0)
  proportional_results = data.frame()
  for (logp in unique(all_results_merged$minus_log_p)) {
    cur = all_results_merged[all_results_merged$minus_log_p == logp, ]
    prop_branch = 100 * (length(cur$cur_branch[cur$cur_branch==phenotype_description & cur$CellType %in% target_cells])/length(cur$cur_branch[cur$CellType %in% target_cells]))
    prop_other = 100 * (length(cur$cur_branch[cur$cur_branch == "Other" & cur$CellType %in% target_cells])/length(cur$cur_branch[cur$CellType %in% target_cells]))
    proportional_results = rbind(proportional_results,
                                 data.frame("cur_branch" = phenotype_description,
                                            "prop" = prop_branch,
                                            "minus_log_p" = logp))
    proportional_results = rbind(proportional_results,
                                 data.frame("cur_branch" = "Other",
                                            "prop" = prop_other,
                                            "minus_log_p" = logp))
  }
  color_pal = wesanderson::wes_palette(wes_color_palette,n_colors)
  prop_plot <- ggplot(proportional_results, aes(x= minus_log_p, y = prop, color = cur_branch)) +
    geom_point(size = 3) +
    geom_line(mapping = aes(linetype= cur_branch),size=0.6) +
    labs(color="", linetype="") +
    ylab(paste0('% Enrichments in "', cell_type_description,'"')) +
    xlab("-log10(p)") +
    scale_y_continuous(limits = c(0,100)) +
    cowplot::theme_cowplot() +
    scale_color_manual(values = c(color_pal[color_expected_phenotypes],color_pal[color_other_phenotypes]))  +
    theme(legend.position = "bottom", legend.text = element_text(size = 6.5), title = element_text(size = 9))
  return(list(prop_plot,proportional_results))
}
