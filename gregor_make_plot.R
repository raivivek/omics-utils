#! /usr/bin/env Rscript
#
# The Parker Lab, 2018
# theparkerlab.org
#

suppressMessages(library("argparse"))
suppressMessages(library("tidyverse"))

parser = ArgumentParser()

parser$add_argument("input", help = '[Required] Path to the GREGOR output directory')
parser$add_argument("--plot", action = 'store',
                    help = '[Optional] Output plots to given directory, if supplied')

args = parser$parse_args()


get_stdev <- function(name, root) {
  nf_path <- file.path(root, name, 'neighborhoodFile.txt')
  nf <- read_tsv(nf_path, skip = 1, col_names = c('x', 'y'), col_types = cols(
    x = col_integer(),
    y = col_double()
  ))
  sum(nf[1]*nf[2]*(1 - nf[2]))^0.5
}


# prepare summary
get_summary <- function(out_dir) {
  summary_path = file.path(out_dir, 'StatisticSummaryFile.txt')
  df = read_tsv(summary_path)
  mutate(df, Stdev = sapply(Bed_File, get_stdev, root=out_dir),
         Zscore = (InBed_Index_SNP - ExpectNum_of_InBed_SNP)/Stdev,
         log2_foldChange = log2(InBed_Index_SNP/ExpectNum_of_InBed_SNP))
}

#summary <- get_summary(args$input)
out_dir <- '/lab/work/vivekrai/liver-project/work/enrichment/random/output/'
summary <- get_summary(out_dir)
summary$trait <- 'random'
summary1 <- get_summary('/lab/work/vivekrai/liver-project/work/enrichment/arthritis/output/')
summary1$trait <- 'arthritis'
summary2 <- get_summary('/lab/work/vivekrai/liver-project/work/enrichment/diamante/output/')
summary2$trait <- 't2d'
summary3 <- get_summary('/lab/work/vivekrai/liver-project/work/enrichment/af/output/')
summary3$trait <- 'atrial_fib'

summary <- rbind(summary, summary1, summary2, summary3)


print(as.data.frame(summary), row.names = F)
theme_ms <- function(base_size=14, base_family="Helvetica") {
  (theme_classic(base_size = base_size, base_family = base_family) +
     theme(text=element_text(color="black"),
           plot.title = element_text(family = base_family, face="bold", size = rel(1.4)),
           plot.background = element_rect(fill = "transparent", color = NA),
           axis.title=element_text(family = base_family, size = rel(1.2)),
           axis.text=element_text(color = "black"),
           legend.title=element_text(size = rel(1.2), face="italic"),
           legend.text=element_text(size = rel(1.1)),
           legend.background=element_rect(fill="transparent"),
           legend.key = element_blank(),
           legend.key.size = unit(1.2, 'lines'),
           panel.border=element_blank(),
           #panel.grid = element_line(color = 'grey96'),
           panel.grid = element_blank(),
           panel.grid.minor = element_blank(),
           panel.spacing = unit(1, "lines"),
           strip.background = element_blank(),
           strip.text = element_text(face = "bold", hjust = 0),
           #axis.line = element_line(color = 'gray')
     )
   )
}
theme_set(theme_ms())

summary <- summary %>% mutate(N = as.integer(str_match(Bed_File, "\\d+"))) %>% arrange(trait, N)
x <- ggplot(summary) + geom_point(aes(N, 2^log2_foldChange)) + ylim(1, 1.75) +
  geom_smooth(aes(N, 2^log2_foldChange), se=F) +
  #geom_text_repel(aes(N, 2^log2_foldChange, label = InBed_Index_SNP), size=3) + 
  scale_x_continuous(breaks=c(1:12)) + facet_wrap(~trait, ncol=4) + labs( y = "fold enrichment") + guides(size = F)
y <- ggplot(summary, aes(N, -log10(PValue))) + geom_point() + geom_smooth(se = F) + geom_hline(yintercept = 1.3, linetype = 'dashed') + scale_x_continuous(breaks=c(1:12)) + facet_wrap(~trait, ncol=4) + labs( y = "log10 P-Value") + theme(strip.text = element_blank())
n <- ggplot(summary, aes(N, InBed_Index_SNP)) + geom_bar(stat = 'identity', color='grey', fill=NA, width=1) + facet_wrap(~trait, ncol=4) + scale_x_continuous(breaks=c(1:12))
n
z <- grid.arrange(x, y, n)

ggplot(summary) + geom_point(aes(N, InBed_Index_SNP, color = trait, group = trait)) +
  geom_errorbar(aes(x = N, ymin = (InBed_Index_SNP - 1.96 * Stdev), ymax = (InBed_Index_SNP + 1.96 * Stdev)), width=.2) +
  # geom_text(aes(N, -log10(PValue), label = InBed_Index_SNP), vjust=-1.5) +
  theme_ms() + theme(legend.title = element_blank()) + scale_x_continuous(breaks = c(1:12))
library(gridExtra)
library(ggrepel)

ggsave(plot = z, file.path('tmp', 'enrichment-stj-gwas-2.pdf'), width = 10, height = 8, dpi=300)
ggsave(plot = z, file.path('tmp', 'enrichment-stj-gwas-2.png'), width = 10, height = 8, dpi=300)
