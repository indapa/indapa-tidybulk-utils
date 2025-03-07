library(tidyverse)
library(tidybulk)
library(stringr)
library(ggrepel)





# Set theme
custom_theme <-
  list(
    
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_line(size = 0.1),
        text = element_text(size = 12),
        legend.position = "bottom",
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
      )
  )

make_tidybulk <- function(df, sample_col , transcript_col, abundance_col) {

  #' returns a tidybulk object from dataframe (tibble)
  #' @param df Dataframe/tibble with at least 3 columns for sampleID, gene symbol, and transcript abundance
  #' @param sample_col The name of the sample column
  #' @param transcript_col The name of the gene/transcript column
  #' @param abundance_col The name of the gene/transcript abundance column (raw counts)
  #' @export

  sample_col <- enquo(sample_col)
  transcript_col <- enquo(transcript_col)
  abundance_col <- enquo(abundance_col)

  return (tidybulk(df, .sample = !!sample_col, .transcript = !!transcript_col, .abundance = !!abundance_col))

}

make_abundance_density_plot <- function(tb_obj, raw_count_col, scaled_count_col, sampleID_col) {

  #' returns a ggplot2 showing density plots of raw and scaled abundance counts from tidybulk object
  #' @param tb_obj tidybulk object with raw and scaled counts
  #' @param raw_count_col column name contaning raw counts
  #' @param scaled_count_col column name containing scaled counts
  #' @param sampleID_col column name name listing containing sampleIDs of analysis
  #' @export

  raw_count_col <- enquo(raw_count_col)
  scaled_count_col <- enquo(scaled_count_col)
  sampleID_col <- enquo(sampleID_col)

  density_plot <- tb_obj  %>%
    pivot_longer(cols = c(!!raw_count_col, !!scaled_count_col), names_to = "source", values_to = "abundance") %>%
    # Plotting
    ggplot(aes(x = abundance + 1, color = !!sampleID_col)) +
    geom_density() +
    facet_wrap(~source) +
    scale_x_log10() +
    custom_theme

  return (density_plot)
}

make_pca_plot <- function(tb_obj, color_factor, shape_factor, label_col, dim1, dim2) {

  #' returns PCA plot of tidybulk abundance scaled object
  #' @param tb_obj tidybulk object with raw and scaled counts
  #' @param color_factor for PCA plot data points
  #' @param shape_factor for PC plot data points
  #' @param label_col for labelling points on data points
  #' @param dim1 PC for x-axis
  #' @param dim2 PC for y-axis
  #' @export

  color_factor <- enquo(color_factor)
  shape_factor <- enquo(shape_factor)
  label_col <- enquo(label_col)
  dim1 <- enquo(dim1)
  dim2 <- enquo(dim2)

  pca_plot <- tb_obj %>%
    reduce_dimensions(method = "PCA") %>%
    pivot_sample() %>%
    ggplot(aes(x = !!dim1, y = !!dim2, colour = !!color_factor, shape = !!shape_factor)) +
    geom_point() +
    geom_text_repel(aes(label = !!label_col), show.legend = FALSE) +
    custom_theme

  return(pca_plot)
}

ciber_deconvolve <- function( tb_obj) {

  #' Perform cibersort deconvolution on tidybulk object
  #' @param  tb_object tidybulk object
  #' returns dataframe with cibersort deconvolution results
  #' @export

  ciber_res <- tb_obj  %>%
    deconvolve_cellularity(prefix="ciber_",cores=1) %>%
    pivot_sample() %>%
    pivot_longer(contains("ciber_"),
                 names_to= "cell_type",
                 values_to="proportion") %>%
    mutate(cell_type = str_replace(cell_type, "ciber_", ""))

  return(ciber_res)
}

make_ciber_plot <- function(ciber_res, category_factor, fill_factor) {

  #' makes stacked bar plot showing cell type proportions of sampleIDs, faceted by subjectID
  #' @param ciber_res is datframe of cibersort results returned from ciber_deconvolve
  #' @param category_factor for plotting cell proportions  (e.g Treated)
  #' @param fill factor to fill the stacked proportions
  #' returns ggplot2 object
  #' @export

  category_factor <- enquo(category_factor)

  fill_factor <- enquo(fill_factor)



  ciber_cell_prop_stacked_bar <- ciber_res %>%

    ggplot(aes(x = !!category_factor,
               y = proportion,
               fill = !!fill_factor)) +
    geom_bar(stat = "identity") +
    custom_theme +
    facet_grid( ~ subjectID) +
    theme(strip.text.x = element_text(size = 10))

  return (ciber_cell_prop_stacked_bar)

}

make_ciber_volcano_table <- function(ciber_res_paired_df, factor_one, factor_two, wilcox_df) {

  #' return dataframe of cell proportions log2FC of cibersort plot of cell type proportions.
  #' @param ciber_res_paired_df is dataframe with cell proportions of two conditions
  #' @param factor_one is the name of the column in ciber_res_paired representing condition one
  #' @param factor_two is the name of the column in ciber_res_paired representing condition one. Will calculate foldchange of factor_two/factor_one
  #' @param wilcox_df is dataframe of listing statistic, p.value, method, alternative (hypothesis) of cell_type proportions from cibersort results
  #' @export

  factor_one <- enquo(factor_one)
  factor_two <- enquo(factor_two)

  volcano_df <- ciber_res_paired_df %>%
    mutate(psuedo_factor_one=!!factor_one+.0001, psuedo_factor_two=!!factor_two+.0001) %>%
    mutate(logFC=log2(psuedo_factor_two/psuedo_factor_one)) %>%
    group_by(cell_type) %>%
    summarise(cell_type_logFC= mean(logFC)) %>%
    merge(wilcox_df)

  return (volcano_df)
}

make_ciber_paired_box_point_plot <- function(ciber_res_paired_df, cell_type_list ) {
  #' return ggplot2 paired plot of cibersort cell fractions Pre_treatment Post_treatment for paired samples
  #' @param  ciber_res_paired_df dataframe retured from
  #' @param cell_type_list list of cell types to analyze
  #' @export

  ciber_df <- ciber_res_paired_df %>%
    pivot_longer(c(Pre_treatment, Post_treatment), names_to = "Treatment", values_to="value") %>%
    filter(cell_type %in% cell_type_list) %>%
    mutate(value=as.numeric(value),
           Treatment=as.factor(Treatment))

  ciber_df$Treatment=fct_relevel(ciber_df$Treatment,c("Pre_treatment", "Post_treatment"))

  p <- ciber_df %>%
    ggplot(aes(Treatment, value, fill=Treatment)) +
    geom_boxplot() +
    geom_line(aes(group = subjectID), linetype=2, size=.8) +
    geom_point(aes(fill=Treatment,group=Treatment),size=1.5,shape=21) +
    facet_wrap(~cell_type , scales = "free") +
    theme(strip.text.x = element_text(size = 10))  +
    theme(legend.position="none")

  return(p)
}


make_ciber_volcano_plot <- function(ciber_res_paired_df, factor_one, factor_two, wilcox_df) {

  #' make a volcano plot of cell type proportions.
  #' @param ciber_res_paired_df is dataframe with cell proportions of two conditions
  #' @param factor_one is the name of the column in ciber_res_paired representing condition one
  #' @param factor_two is the name of the column in ciber_res_paired representing condition one. Will calculate foldchange of factor_two/factor_one
  #' @param wilcox_df is dataframe of listing statistic, p.value, method, alternative (hypothesis) of cell_type proportions from cibersort results
  #' @export

  factor_one <- enquo(factor_one)
  factor_two <- enquo(factor_two)

  volcano_plot <- ciber_res_paired_df %>%
    mutate(psuedo_factor_one=!!factor_one+.0001, psuedo_factor_two=!!factor_two+.0001) %>%
    mutate(logFC=log2(psuedo_factor_two/psuedo_factor_one)) %>%
    group_by(cell_type) %>%
    summarise(cell_type_logFC= mean(logFC)) %>%
    merge(wilcox_df) %>%
    ggplot(aes(cell_type_logFC, p.value, label=cell_type)) +

    geom_point(aes(size=3)) +
    xlab("log2 FC in cell type proportion") +
    geom_text_repel(size = 4, max.overlaps = Inf) +
    ggtitle("Volcano plot Cell type proportions", subtitle="Paired analysis") +

    scale_y_continuous(trans = "log10_reverse") +
    geom_vline(xintercept = 1, linetype="dotted", color = "black", size=.8)  +
    geom_vline(xintercept = -1, linetype="dotted", color = "black", size=.8) +
    theme(
      axis.text.x = element_text(size = 16, vjust = 0.65),
      axis.text.y = element_text(size = 16),
      axis.title = element_text(size = 16, color = 'black'),
      strip.text.x = element_text(size = 25),
      strip.background = element_blank(),
      legend.position="none")

  return (volcano_plot)
}

make_tcon_treg_proportions_table <- function(ciber_res, category_factor, subject_column, sample_column) {
  #' return dataframe Tcon/Treg cell type proportions
  #' @param ciber_res dataframe with tidybulk cibersort results
  #' @param category_factor column representing with category levels you want Treg/Tcon cell proportions for
  #' @param subject_column column representing with subjectID
  #' @param sample_column column representing with sample name
  #' @export

  category_factor <- enquo(category_factor)
  subject_column <- enquo(subject_column)
  sample_column <- enquo(sample_column)

  #get the Tcells and then classify as Treg or Tcon
  tcon_treg_res <- ciber_res %>%
    filter(str_detect(cell_type, "T cells")) %>%
    filter(cell_type != "T cells follicular helper") %>%
    filter(cell_type != "T cells gamma delta") %>%
    mutate(Tcell_type = ifelse(cell_type == 'T cells regulatory (Tregs)', "Treg", 'Tcon')) %>%
    group_by(!!sample_column, !!subject_column, !!category_factor, Tcell_type) %>% summarise(Tcell_type_prop=sum(proportion)) %>%
    pivot_wider(names_from = Tcell_type, values_from = Tcell_type_prop)

  return (tcon_treg_res)

}

make_tcon_treg_proportions_plot <- function(tcon_treg_table, shape_factor, color_factor, subtitle_string ='' ) {

  #' make scatter plot of Tcon vs Treg cell type proproportions from cibersort results
  #' @param tcon_treg_table from make_tcon_treg_proportions_table
  #' @param shape_factor name of column that you want to map shape of points to in plot
  #' @param color_factor name of column you want to map color of points to in plot
  #' @export


  shape_factor <- enquo(shape_factor)
  color_factor <- enquo(color_factor)

  tcon_treg_plot <-  ggplot(tcon_treg_table, aes(Tcon, Treg, shape=!!shape_factor, colour=!!color_factor)) +
    geom_point() +
    ggtitle("Tcon vs Treg cell type proportions", subtitle = subtitle_string)

  return(tcon_treg_plot)

}

get_tidybulk_topgenes_pval <- function(counts_de, n_top=10 ) {

  #' return n top differentially expressed gene symbols from tidybulk differential gene expression results 
  #' @param counts_de tidybulk differential expression result
  #' @param n_top return the n_top gene symbols
  #' @export
  
    res <- counts_de %>%
    arrange(PValue) %>%
    head(n_top) %>%
    pull(gene)
    
    return (res)
}

make_tidybulk_volcano_plot_top_genes <- function(counts_de, gene_symbols, FC_lim = 1,title_string="Tidybulk volcano plot", subtitle_string="" ) {
  
  #' return ggplot volcano plot from tidybulk differential expression result, labeling gene_symbols  
  #' @param counts_de tidybulk differential expression result
  #' @param gene_symbols list of gene symbols to label in volcano plot
  #' @param FDR_lim default is 0.05
  #' @param FC_lim fold change default is 2
  #' @export
  
  volcano_plot <- counts_de %>%
     mutate(gene = ifelse(gene %in% gene_symbols, as.character(gene), "")) %>%
    mutate(topgene = gene %in% gene_symbols) %>%
    ggplot(aes(x = logFC, y = PValue, label = gene)) +
    geom_point(aes(color = topgene, size=topgene, alpha = topgene)) +
    geom_text_repel(size =7, max.overlaps = Inf, box.padding = 0.5) + 
    
    # Custom scales
    custom_theme +
    scale_y_continuous(trans = "log10_reverse") +
    scale_color_manual(values = c("black", "#e11f28")) +
    scale_size_discrete(range = c(0, 2)) + 
    geom_vline(xintercept = FC_lim, linetype="dotted", color = "black", size=.8)  + 
    geom_vline(xintercept = -FC_lim, linetype="dotted", color = "black", size=.8) + 
    xlab("log2FoldChange") + 
    ylab("p.value") +
    ggtitle("Volcano plot", subtitle = subtitle_string )  + 
    theme(
      axis.text.x = element_text(size = 20, vjust = 0.65),
      axis.text.y = element_text(size = 20), 
      axis.title = element_text(size = 20, color = 'black'),
      strip.text.x = element_text(size = 30),
      strip.background = element_blank(), 
      legend.position="none")
    
  
    return (volcano_plot)
  
}

make_tidybulk_volcano_plot <- function(counts_de,pval_limit= 0.05, FC_lim=1, title_string="Tidybulk volcano plot", subtitle_string="") {
  
  #' return ggplot volcano plot from tidybulk differential expression result  
  #' @param counts_de tidybulk differential expression result
  #' @param top_gene_symbols list of gene symbols to label in volcano plot
  #' @param FDR_lim default is 0.05
  #' @param FC_lim fold change default is 2
  #' @export
  #' 
  #' 
  
  
  counts_de$significant <- ifelse(counts_de$PValue < pval_limit & abs(counts_de$logFC) >= FC_lim, TRUE, FALSE)
  top_50_sig <- counts_de %>% filter(significant == TRUE) %>%  arrange(PValue) %>% head(50)
  volcano_plot <- counts_de %>%
   
    mutate(gene = ifelse(significant == TRUE, as.character(gene), "")) %>%
    
    # Plot
    ggplot(aes(x = logFC, y = PValue, colour= significant)) +
    geom_point(aes(color = significant, size = significant, alpha = significant)) +
    geom_text_repel(
     data = top_50_sig,
      aes(label = gene),
      size = 4,
     max.overlaps = Inf,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")) + 
    geom_vline(xintercept =FC_lim , linetype="dotted", color = "black", size=.8)  + 
    geom_vline(xintercept = -1 *FC_lim, linetype="dotted", color = "black", size=.8) + 
    ggtitle(label = title_string, subtitle = subtitle_string) + 
    
    # Custom scales
    custom_theme +
    scale_y_continuous(trans = "log10_reverse") +
    scale_color_manual(values = c("black", "#e11f28")) +
    scale_size_discrete(range = c(0, 2)) + 
    
    theme(
      axis.text.x = element_text(size = 20, vjust = 0.65),
      axis.text.y = element_text(size = 20), 
      axis.title = element_text(size = 20, color = 'black'),
      strip.text.x = element_text(size = 30),
      strip.background = element_blank(), 
      legend.position="none")
  
    return (volcano_plot)
  
}