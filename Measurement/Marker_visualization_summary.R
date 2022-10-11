library(tidyverse)
library(Matrix)
library(affy)
library(stringr)



# 12.10 & 12.23
B.m.selected <- read.csv('D:/SRTP2021/12/df_genes_final.csv')
B.m.selected <- B.m.selected[,-1]

min_max_nor <- function(x){
  (x-min(x))/(max(x)-min(x))
}

B.m.selected.nor <- as.data.frame(apply(B.m.selected, 2, FUN=min_max_nor))
B.m.selected.nor$target <- B.m.selected$target

  # normalize with internal ctrl
normalize_internal <- function(x){
  (x/B.m.selected.nor$Gapdh)
}
B.m.selected.nor.in <- as.data.frame(apply(B.m.selected.nor, 2, FUN=normalize_internal))
B.m.selected.nor.in <- B.m.selected.nor.in[,-16]
B.m.selected.nor.in <- gather(B.m.selected.nor.in, key='gene', value='expr', c(-'target'))
B.m.selected.nor.in$target <- B.m.selected$target

B.m.selected.nor <- gather(B.m.selected.nor, key='gene', value='expr', c(-'target'))

  # normalized with other markers
ggplot(B.m.selected.nor) +
  geom_boxplot(aes(x=target, group=target, y=expr),alpha=0.1) +
  facet_wrap(~gene, ncol=4) +
  xlab('Age') +
  ylab('Normalized Expression') +
  ggtitle('B Marrow Expression') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(B.m.selected) +
  geom_boxplot(aes(x=target,group=target, y = Cdkn1a))

  # unnormalized with internal controls
B.m.selected.gather <- gather(B.m.selected, key='gene', value='expr', c('Cd52',
                                                                        'Actb','Ly6e'))
ggplot(B.m.selected.gather) +
  geom_boxplot(aes(x=target, group=target, y=expr),alpha=0.1) +
  facet_wrap(~gene) +
  xlab('Age') +
  ylab('Raw Counts') +
  ggtitle('B Marrow Raw Counts with Internal Controls')


# B lung
B.lu.selected <- read.csv('D:/SRTP2021/12/df_lung_genes_final.csv')
B.lu.selected.nor <- as.data.frame(apply(B.lu.selected, 2, FUN=min_max_nor))
B.lu.selected.nor$target <- B.lu.selected$target
B.lu.selected.nor <- gather(B.lu.selected.nor, key='gene', value='expr', c(-'target'))

B.lu.selected.nor %>%
  mutate(across(gene, factor, levels=c('Cyba','Ifitm2','Ly6e',
                                       'Plaur','Gapdh','Trp53','Cdkn1a','Cdkn2a'))) %>%
ggplot(aes(color=target)) +
  geom_boxplot(aes(x=target, group=target, y=expr),alpha=0.1) +
  facet_wrap(~gene, ncol=4) +
  scale_color_continuous(low="blue",high = "red") +
  xlab('Age') +
  ylab('Normalized Expression') +
  ggtitle('B Lung Expression') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
    
  # log2 transform
ggplot(B.lu.selected.nor) +
  geom_boxplot(aes(x=target, group=target, y=expr),alpha=0.1) +
  # scale_y_continuous(trans = 'log2') +
  facet_wrap(~gene) +
  xlab('Age') +
  ylab('Normalized Expression') +
  ggtitle('B Lung Expression')

  # unnormalized with internal controls
B.lu.selected.gather <- gather(B.lu.selected, key='gene', value='expr', c('Cd52',
                                                                        'Gapdh','Ly6e'))
ggplot(B.lu.selected.gather) +
  geom_boxplot(aes(x=target, group=target, y=expr),alpha=0.1) +
  facet_wrap(~gene) +
  xlab('Age') +
  ylab('Raw Counts') +
  ggtitle('B Lung Raw Counts with Internal Controls')


# B liver
B.li.selected <- read.csv('D:/SRTP2021/12/df_liver_genes_final.csv')
B.li.selected.nor <- as.data.frame(apply(B.li.selected, 2, FUN=min_max_nor))
B.li.selected.nor$target <- B.li.selected$target
B.li.selected.nor <- gather(B.li.selected.nor, key='gene', value='expr', c(-'target'))

B.li.selected.nor %>%
  mutate(across(gene, factor, levels=c('Cyba','Ifitm2','Ly6e',
                                       'Plaur','Gapdh','Trp53','Cdkn1a','Cdkn2a'))) %>%
ggplot(aes(color=target)) +
  geom_boxplot(aes(x=target, group=target, y=expr),alpha=0.1) +
  facet_wrap(~gene, ncol=4) +
  scale_color_continuous(low="blue",high = "red") +
  xlab('Age') +
  ylab('Normalized Expression') +
  ggtitle('B Liver Expression') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

  # unnormalized with internal controls
B.li.selected.gather <- gather(B.li.selected, key='gene', value='expr', c('Cd52',
                                                                          'Actb','Ly6e'))
ggplot(B.li.selected.gather) +
  geom_boxplot(aes(x=target, group=target, y=expr),alpha=0.1) +
  facet_wrap(~gene) +
  xlab('Age') +
  ylab('Raw Counts') +
  ggtitle('B Liver Raw Counts with Internal Controls')

# B spleen
B.s.selected <- read.csv('D:/SRTP2021/12/spleen/subset_more.csv')
B.s.selected <- B.s.selected[,-1]
B.s.selected.nor <- as.data.frame(apply(B.s.selected, 2, FUN=min_max_nor))
B.s.selected.nor$target <- B.s.selected$target
B.s.selected.nor <- gather(B.s.selected.nor, key='gene', value='expr', c(-'target'))

ggplot(B.s.selected.nor) +
  geom_boxplot(aes(x=target, group=target, y=expr),alpha=0.1) +
  facet_wrap(~gene) +
  xlab('Age') +
  ylab('Normalized Expression') +
  ggtitle('B Spleen Expression') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# B heart
markers <- 
df_selected_final <- read.csv("D:/SRTP2021/12/heart/df_heart.csv")
  # NO VSIR
df_selected_final.fil <- df_selected_final %>%
  subset(select=c("Cdkn1a", "Ly6e", "Trp53", "Cyba", "Ifitm2", "Gapdh", 
                  "Icam1","Anxa3","Tspan8","Fxyd5","target"))
B.h.selected.nor <- as.data.frame(apply(df_selected_final.fil[,-11], 
                                        2, FUN=min_max_nor))
B.h.selected.nor$target <- df_selected_final.fil$target
B.h.selected.nor <- gather(B.h.selected.nor, key='gene', value='expr', c(-'target'))
# write.csv(B.h.selected.nor, 'D:/SRTP2021/12/heart/11_markers.csv')

ggplot(B.h.selected.nor) +
  geom_boxplot(aes(x=target, group=target, y=expr),alpha=0.1) +
  facet_wrap(~gene) +
  xlab('Age') +
  ylab('Normalized Expression') +
  ggtitle('B Heart Expression') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



# E
E <- read.csv("D:\\SRTP2021\\14_human_proteome\\merged_data.csv")
E <- E[,-1]
rownames(E) <- E[,1]
E <- E[,-1]
E <- as.data.frame(t(E))
E <- E[,c('CD52','SLC14A1','HCST','CD69','S100A6','PLAUR',
           'CD19','CD22','ACTB','GAPDH')] # A: 111 129 Y: 161 141
E$age <- c(rep(59.2,111), rep(62.2,129), rep(30.6,161), rep(21.2,141))

E_norm <- as.data.frame(apply(E, 2, FUN=min_max_nor))
E_norm$age <- E$age

E_gather <- gather(E_norm,key='gene',value='expr',c(-11))

ggplot(E_gather) +
  geom_boxplot(aes(x=age,group=age, y=expr),alpha=0.1) +
  facet_wrap(~gene) +
  xlab('Age') +
  ylab('Normalized Expression') +
  ggtitle('E Expression')

  # unnormalized with internal controls
E_gather_unnor <- gather(E,key='gene',value='expr',c(-11))

ggplot(E_gather_unnor) +
  geom_boxplot(aes(x=age, group=age, y=expr),alpha=0.1) +
  facet_wrap(~gene) +
  xlab('Age') +
  ylab('Raw Counts') +
  ggtitle('E Raw Counts with Internal Controls')



# FX
F_ <- read.csv('D:\\SRTP2021\\15_human_rnaseq\\four_sampleX.csv')
F_ <- na.omit(F_)
rownames(F_) <- F_$name
F_ <- F_[,colnames(F_)!='name']

F_ <- as.data.frame(t(F_))
F_ <- apply(F_, 2, as.integer)
F_norm <- as.data.frame(apply(F_, 2, FUN=min_max_nor))
F_norm$aged <- c(rep('24',234), rep('28',64), rep('67',94),rep('84',67))

F_norm_gathered <- gather(F_norm, key='gene', value='expr', c(-'aged'))

ggplot(F_norm_gathered) +
  geom_boxplot(aes(x=aged, y=expr),alpha=0.1) +
  facet_wrap(~gene) +
  xlab('Age') +
  ylab('Normalized Expression') +
  ggtitle('F Expression') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

  # unnormalized with internal controls
F_$aged <- c(rep('Young',234), rep('Aged',67))
F_gather_unnor <- gather(F_ ,key='gene',value='expr',c(-'aged'))

ggplot(F_gather_unnor) +
  geom_boxplot(aes(x=aged, group=aged, y=expr),alpha=0.1) +
  facet_wrap(~gene) +
  xlab('Age') +
  ylab('Raw Counts') +
  ggtitle('F Raw Counts with Internal Controls')





# C
C <- read.csv('D:/SRTP2021/7/data.csv')
C$aged <- c(rep(0, 78), rep(1, 71))
C <- C[,-1]
C.sel <- C[,c('aged', 'Icam2', 'Rpsa')]
C.sel.nor <- as.data.frame(apply(C.sel, 2, FUN=min_max_nor))
C.sel.nor$aged <- c(rep('young', 78), rep('aged', 71))
C.sel.nor.ga <- gather(C.sel.nor, key='gene', value='expr', c(-'aged'))

ggplot(C.sel.nor.ga) +
  geom_boxplot(aes(x=aged, y=expr),alpha=0.7) +
  facet_wrap(~gene) +
  xlab('Age') +
  ylab('Normalized Expression') +
  ggtitle('C Expression (n=149)')



# G
G <- read.csv('D:/SRTP2021/16_rat/G_plot.csv')
G <- as.data.frame(t(G))
colnames(G) <- G[1417,]
G <- G[-1417,]
G <- as.data.frame(apply(G, 2, as.integer))
G.nor <- as.data.frame(apply(G, 2, FUN=min_max_nor))
G.nor$aged <- c(rep('aged',727), rep('young',689))

G.nor.ga <- gather(G.nor, key='gene', value='expr', c(-'aged'))

ggplot(G.nor.ga) +
  geom_boxplot(aes(x=aged, y=expr),alpha=0.7) +
  facet_wrap(~gene) +
  xlab('Age') +
  ylab('Normalized Expression') +
  ggtitle('G Expression') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# COVID
covid <- read.csv('D:/SRTP2021/17_covid/df_plot.csv')
covid_filt <- covid[,c('PLAUR', 'CYBA', 'IFITM2', 'target', 'CDKN2A', 'CDKN1A', 'GAPDH', 'TP53', 'LY6E')]
covid.nor <- as.data.frame(apply(covid_filt, 2, FUN=min_max_nor))
covid.nor$target <- covid$target
covid.gath <- gather(covid.nor, key='gene', value='expr', c(-'target'))

covid.gath %>%
  mutate(across(gene, factor, levels=c('CYBA','IFITM2','LY6E',
                                       'PLAUR','GAPDH','TP53','CDKN1A','CDKN2A'))) %>%
ggplot() +
  geom_boxplot(aes(group=target, x=target, y=expr, color=target),alpha=0.7) +
  facet_wrap(~gene, ncol=4) +
  scale_color_continuous(low="blue",high = "red") +
  xlab('Age') +
  ylab('Normalized Expression') +
  ggtitle('COVID Expression') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Volcano Plot
library(EnhancedVolcano)
res <- read.csv('D:\\SRTP2021\\16_rat\\res.csv')
res_filter <- res[res$pvalue<1e-30,]
EnhancedVolcano(res_filter, lab=res_filter$X, x='log2FoldChange', y='pvalue')


# G heart
plot_plotdf <- function(filename, organ) {
  G_heart <- read.csv(paste('D:/SRTP2021/22_rat_revisited/',filename,'.csv', sep = ''))
  G_heart <- G_heart[,-1]
  g.heart <- as.data.frame(apply(G_heart[,1:11], 2, FUN=min_max_nor))
  g.heart$target <- G_heart$target
  g.gath <- gather(g.heart, key='gene', value='expr', c(-'target'))
  
  ggplot(g.gath) +
    geom_boxplot(aes(group=target, x=target, y=expr, color=target),alpha=0.7) +
    facet_wrap(~gene, ncol=6) +
    xlab('Age') +
    ylab('Normalized Expression') +
    ggtitle(paste('G ',organ,' Expression',sep = '')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
}
plot_plotdf('heart_11_markers', 'Aorta')

# G liver
plot_plotdf("G_liver_8", 'Liver')

# Monkey CA & AA
plot_plotdf_monkey <- function(filename, organ) {
  G_heart <- read.csv(paste('D:/SRTP2021/23_monkey/',filename,'.csv', sep = ''))
  G_heart <- G_heart[,-1]
  G_heart$target <- str_sub(G_heart$target, 1, 1)
  g.heart <- as.data.frame(apply(G_heart[,1:10], 2, FUN=min_max_nor))
  g.heart$target <- G_heart$target
  g.gath <- gather(g.heart, key='gene', value='expr', c(-'target'))
  
  ggplot(g.gath) +
    geom_boxplot(aes(group=target, x=target, y=expr, color=target),alpha=0.7) +
    facet_wrap(~gene, ncol=5) +
    xlab('Age') +
    ylab('Normalized Expression') +
    ggtitle(paste('Monkey ',organ,' Expression',sep = '')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
}
plot_plotdf_monkey("AA_11_markers", "Aortic Artery")
plot_plotdf_monkey("CA_11_markers", "Coronary Artery")

