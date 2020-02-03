# devtools::install_github("omarwagih/ggseqlogo")
require(ggplot2)
require(ggseqlogo)
library(readr) # for data parsing / loading in data files
library(tidyverse) # for organizing data structures

main <- function()
{
  ###----------- create sequence logos for promoter data -----------###
  EPDnew_promoters_all_df <- read_csv("../../data/promoters_for_ML/EPDnew_promoters_for_ML.csv")
  RegulonDB_promoters_all_df <- read_csv("../../data/promoters_for_ML/RegulonDB_promoters_for_ML.csv")
  SEQ_LOGO_REGDB_DIR <- '../../figs/sequence_logos/RegulonDB/'
  SEQ_LOGO_EPDNEW_DIR <- '../../figs/sequence_logos/EPDnew/'
  
  # specify color-blind friendly colors for sequence logos
  color_scheme <- make_col_scheme(chars=c('A', 'T', 'C', 'G'),
                        cols=c('#E5A224', '#2C9F74', '#0277B4', '#CC7DAA'))
  
  # specify new indices with respect to TSS
  orig_idx <- seq(from = 1, to = 41, by = 5)
  new_idx <- seq(from = -40, to = 0, by = 5)
  
  ###----------- create sequence logos for RegulonDB -----------###
  # RegulonDB sigma factors
  sigma_factors_Ecoli = c('unknown', 
                          'Sigma70',
                          'Sigma54',
                          'Sigma38',
                          'Sigma32',
                          'Sigma28',
                          'Sigma24')
                          #'Sigma19')
  
  # create sequence logo for each sigma factor
  for (sigma_factor in sigma_factors_Ecoli) {
    
    # grab positive and negative permuted promoters
    RegulonDB_promoters_positive_idx <- which(
                                              RegulonDB_promoters_all_df$`class` == 1 & 
                                              RegulonDB_promoters_all_df$`database/source` == 'RegulonDB' &
                                              apply(RegulonDB_promoters_all_df, 1, function(x) any(grepl(sigma_factor, x)))
                                              )
    RegulonDB_promoters_negative_permuted_idx <- which(
                                                      RegulonDB_promoters_all_df$`class` == 0 & 
                                                      RegulonDB_promoters_all_df$`database/source` == 'RegulonDB (permuted)' &
                                                      apply(RegulonDB_promoters_all_df, 1, function(x) any(grepl(sigma_factor, x)))
                                                      )
    # grab nucleotide sequences
    RegulonDB_promoters_positive_df <- RegulonDB_promoters_all_df[RegulonDB_promoters_positive_idx,]
    RegulonDB_promoters_neg_permuted_df <- RegulonDB_promoters_all_df[RegulonDB_promoters_negative_permuted_idx,]
    RegulonDB_promoters_pos_seq <- RegulonDB_promoters_positive_df$`DNA sequence`
    RegulonDB_promoters_neg_perm_seq <- RegulonDB_promoters_neg_permuted_df$`DNA sequence`
  
    # get segments from -40 to +2
    RegulonDB_promoters_pos_seq <- sapply(RegulonDB_promoters_pos_seq,substr, start=20, stop=62)
    RegulonDB_promoters_neg_perm_seq <- sapply(RegulonDB_promoters_neg_perm_seq,substr, start=20, stop=62)
  
    # plot sequence logos
    pos_logo <- ggseqlogo( RegulonDB_promoters_pos_seq, font = 'roboto_bold',col_scheme=color_scheme )  + scale_x_continuous(breaks = orig_idx,labels = new_idx )+ ylim(0, 2)+
      theme_classic(base_size = 24)
    neg_perm_logo <- ggseqlogo( RegulonDB_promoters_neg_perm_seq, font = 'roboto_bold',col_scheme=color_scheme) +labs(x='Nucleotide position (with respect to TSS)') +scale_x_continuous(breaks = orig_idx,labels = new_idx )+ ylim(0, 2)+
      theme_classic(base_size = 24)
    
    final_plot <- gridExtra::grid.arrange(pos_logo, neg_perm_logo) 
    plot_name <- paste('RegulonDB', sigma_factor,sep='_')
    
    # save plots
    ggsave(paste(SEQ_LOGO_REGDB_DIR, plot_name,'.pdf', sep=''), plot = final_plot, 
           device = "pdf",
           height = 8,
           width = 10,
           scale = 1,
           dpi = 400)
    
    ggsave(paste(SEQ_LOGO_REGDB_DIR, plot_name,'.png', sep=''), plot = final_plot, 
           device = "png",
           height = 8,
           width = 10,
           scale = 1,
           dpi = 400)
  }
  
  ###----------- create sequence logos for EPDnew -----------###
  # binding motifs in EPDnew data set 
  motifs_EPDnew_withNA <- unique(EPDnew_promoters_all_df$`sigma factor/motif`)
  motifs_EPDnew <- motifs_EPDnew_withNA[!is.na(motifs_EPDnew_withNA)]
  
  # organisms in EPDnew data set
  organism_EPDnew_withNA <- unique(EPDnew_promoters_all_df$`organism`)
  organism_EPDnew <- organism_EPDnew_withNA[!is.na(organism_EPDnew_withNA)]
  
  # create sequence logo for every pair of motif and organism
  for (motif in motifs_EPDnew) {
    for (organism in organism_EPDnew){
      
      # p. falciparum does not have TATA and nonTATA motifs. Similarly, only p. falciparum is the only organism with 
      # Inr-pfa and nonInr-pfa
      if ((motif == 'Inr-pfa' | motif == 'nonInr-pfa') & (organism != 'p_falciparum')) next
      if ((motif == 'TATA' | motif == 'nonTATA') & (organism == 'p_falciparum')) next
      
      EPDnew_promoters_positive_idx <- which(
                                            EPDnew_promoters_all_df$`class` == 1 & 
                                            EPDnew_promoters_all_df$`database/source` == 'EPDnew' &
                                            EPDnew_promoters_all_df$`sigma factor/motif` == motif &
                                            EPDnew_promoters_all_df$organism == organism
                                            )
      
      EPDnew_promoters_negative_permuted_idx <- which(
                                            EPDnew_promoters_all_df$`class` == 0 & 
                                            EPDnew_promoters_all_df$`database/source` == 'EPDnew (permuted)' &
                                            EPDnew_promoters_all_df$`sigma factor/motif` == motif &
                                            EPDnew_promoters_all_df$organism == organism
                                            )    
      
      # break loop if no promoters were found
      if (length(EPDnew_promoters_positive_idx) == 0) next
      
      
      EPDnew_promoters_positive_df <- EPDnew_promoters_all_df[EPDnew_promoters_positive_idx,]
      EPDnew_promoters_neg_permuted_df <- EPDnew_promoters_all_df[EPDnew_promoters_negative_permuted_idx,]
      EPDnew_promoters_pos_seq <- EPDnew_promoters_positive_df$`DNA sequence`
      EPDnew_promoters_neg_perm_seq <- EPDnew_promoters_neg_permuted_df$`DNA sequence`
  
      EPDnew_promoters_pos_seq <- sapply(EPDnew_promoters_pos_seq, substr, start=460, stop=502)
      EPDnew_promoters_neg_perm_seq <- sapply(EPDnew_promoters_neg_perm_seq, substr, start=460, stop=502)
  
      # plot sequence logos
      pos_logo <- ggseqlogo( EPDnew_promoters_pos_seq, font = 'roboto_bold',col_scheme=color_scheme )  + scale_x_continuous(breaks = orig_idx,labels = new_idx )+ ylim(0, 2)+
        theme_classic(base_size = 24)
      neg_perm_logo <- ggseqlogo( EPDnew_promoters_neg_perm_seq, font = 'roboto_bold',col_scheme=color_scheme ) +labs(x='Nucleotide position (with respect to TSS)')+ scale_x_continuous(breaks = orig_idx,labels = new_idx )+ ylim(0, 2)+
        theme_classic(base_size = 24)
  
      final_plot <- gridExtra::grid.arrange(pos_logo, neg_perm_logo) 
      plot_name <- paste('EPDnew', organism, motif,sep='_')
      
      ggsave(paste(SEQ_LOGO_EPDNEW_DIR, plot_name,'.pdf', sep=''), plot = final_plot, 
             device = "pdf",
             height = 8,
             width = 10,
             scale = 1,
             dpi = 400)
      
      ggsave(paste(SEQ_LOGO_EPDNEW_DIR, plot_name,'.png', sep=''), plot = final_plot, 
             device = "png",
             height = 8,
             width = 10,
             scale = 1,
             dpi = 400)
      }
  }
}
###----------- sequence logos for random promoters -----------###

# pretty much just empty, so not interesting to plot
# RegulonDB_promoters_negative_random_idx <- which(
#                                                   RegulonDB_promoters_all_df$`class` == 0 &
#                                                   RegulonDB_promoters_all_df$`database/source` == 'random'
#                                                   )
# 
# RegulonDB_promoters_neg_random_df <- RegulonDB_promoters_all_df[RegulonDB_promoters_negative_random_idx,]
# 
# RegulonDB_promoters_neg_random_seq <- RegulonDB_promoters_neg_random_df$`DNA sequence`
# 
# RegulonDB_promoters_neg_random_seq <- sapply(RegulonDB_promoters_neg_random_seq,substr, start=20, stop=62)
# neg_rand_logo <- ggseqlogo( RegulonDB_promoters_neg_random_seq, font = 'roboto_bold',col_scheme=color_scheme )+labs(x='Nucleotide position (with respect to TSS)')+ scale_x_continuous(breaks = orig_idx,labels = new_idx )+ ylim(0, 2)+
#   theme_classic(base_size = 24)
# neg_rand_logo
main()
