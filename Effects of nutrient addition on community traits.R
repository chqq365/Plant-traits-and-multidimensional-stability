
## use sites bnch.us, bogong.au, burrawan.au, cbgb.us, cowi.ca, hopl.us, look.us, mcla.us, sgs.us, sier.us

rm(list=ls())
##### Load packages
library(tidyverse);library(xtable);library(factoextra);library(FactoMineR)
library(ggplot2);library(cowplot);library(GGally); library(svglite)
library(nlme);library(piecewiseSEM); library(FD)

# set up the work directory 
dir.data<-"C:/Users/chqq3/work/NutNet data/"
dir.graphs<-"C:/Users/chqq3/work/traits and stability/graphs/"
setwd(dir.graphs)

## environment needed
pd<-position_dodge(width=0.5)
mt<-theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
###################################################################################################
#############################leaf trait and cover data from NutNet ################################
###################################################################################################

avg.traits<-read.csv(file="five leaf traits from 10 NutNet sites.csv")%>%mutate(X=NULL)

## calculate an integrated leaf trait (slow-fast trait)
## check for NAs
w.nutnet2<-avg.traits%>% pivot_wider(names_from ="TraitName1", values_from = "avg.value" )%>%
  filter(!if_any(c("Leaf C", "Leaf N", "Leaf P", "Leaf K", "SLA"), is.na))
length(unique(w.nutnet2$site_code))## summ.za does not have any missing data 

d7<-read.csv(file="cover and biomass data from 10 NutNet sites.csv")
d.sites.select<-d7
(n.spp<-length(unique(d.sites.select$standard_taxon)))## 
spp.with.field.traits<-w.nutnet2%>%ungroup()%>%select(standard_taxon)%>%distinct()# species with all leaf traits

# get fast-slow traits
pca.r<-PCA(w.nutnet2[,c("Leaf C", "Leaf N", "Leaf P", "Leaf K", "SLA")], scale=T, graph=F)
pca.r.without.K<-PCA(w.nutnet2[,c("Leaf C", "Leaf N", "Leaf P", "SLA")], scale=T, graph=F)
summary(pca.r)##
(loadings <- as.data.frame(pca.r$ind$coord))
(loadings_without_K <- as.data.frame(pca.r.without.K$ind$coord))
## add this fast and slow species 
w.nutnet2$slow_fast<-loadings$Dim.1
w.nutnet2$slow_fast_without_K<-loadings_without_K$Dim.1

## Make a biplot of individuals and variables, show species in control and NPK by different color 
w.nutnet2[,1:3]<-lapply(w.nutnet2[,1:3], factor)
(pp.pca.field.traits<-fviz_pca_biplot(pca.r, label = "var", alpha.ind =0.5, pointsize=2, habillage = w.nutnet2$treatment, 
                                      title = paste0("PCA based on ", nrow(spp.with.field.traits), " species having all traits (in total ",n.spp, " species occurred)"))+ theme_bw(base_size=18)+mt)
ggsave(pp.pca.field.traits, height=6.64, width=13.3, file="pac for five leaf traits measured at nutnet sites.png")

# add the slow-fast trait to other single trait values 
avg.traits1<-w.nutnet2%>%select(site_code, treatment, standard_taxon, slow_fast, slow_fast_without_K)%>%
  pivot_longer(names_to="TraitName1", values_to = "avg.value", cols = c("slow_fast", "slow_fast_without_K"))%>%
  select(site_code, treatment, standard_taxon, TraitName1, avg.value)%>%
  bind_rows(avg.traits)%>%mutate(trt=ifelse(treatment=="Ambient", "Control", "NPK"), TraitName=TraitName1, TraitName1=NULL)
length(unique(avg.traits1$site_code))
# check percent species with traits
field.trait.inf1<-avg.traits1%>%group_by(TraitName, trt)%>%summarise(percent.species.with.traits=round((length(avg.value)/n.spp)*100), 2)%>%
  pivot_wider(names_from = "trt", values_from ="percent.species.with.traits" )%>%mutate(percent.species.with.traits.control.npk=paste(Control, NPK, sep="/"))%>%
  select(TraitName, percent.species.with.traits.control.npk)

###########################################################################################
#######################calculate functional diversity######################################
###########################################################################################
# add species with traits to the cover data 
# focus on cover data from year_trt 1, 
d8.l<-d7%>%filter(max_cover>0)%>%select(site_code, trt, standard_taxon, block, plot,  trt, year_trt , year, richness, max_cover)%>%unique()%>%
  group_by(site_code, block, plot,  trt, year_trt)%>%mutate(total_cover=sum(max_cover))%>%
  filter(year_trt>0)%>%distinct()%>%
  mutate(plot.id=paste(site_code, block, plot,  trt, year_trt, sep="_"))%>%
  merge(avg.traits1, by=c("site_code", "trt", "standard_taxon"))
length(unique(avg.traits1$site_code))
length(unique(d8.l$standard_taxon))

hist(d8.l$total_cover)
check.na.t<-d8.l[is.na(d8.l$avg.value),]
check.na.c<-d8.l[is.na(d8.l$max_cover),]

fun<-c()
for(tr in unique(d8.l$TraitName)){
  for (pl in unique(d8.l$plot.id)){
    # pl<-"mtca.au_1_2_NPK_5" ; tr<-"leaf_pct_N"
    temp.data<-d8.l%>%filter(plot.id==pl & TraitName==tr)
    sp.tr<-temp.data$avg.value
    names(sp.tr)<-temp.data$standard_taxon
    abund1<-temp.data$max_cover
    names(abund1)<-temp.data$standard_taxon
    total_cover_traits<-sum(abund1) # total cover with traits
    # 1 species is not possible to calculate fd, but possible to calculate cwm, two species with same trait value is not possible either
    if(length(unique(sp.tr))<2) 
      next
    (fd.temp <- dbFD(sp.tr, abund1, calc.FRic = F, calc.FDiv = F, w.abun=T))##  
    temp.traits<-data.frame(plot.id=pl, f.dispersion=fd.temp$FDis, cwm1=fd.temp$CWM[,1],
                            f.richness=fd.temp$sing.sp, total_cover_traits=total_cover_traits, trait.name=tr)
    fun<-rbind(fun, temp.traits)
  }
}
colnames(fun)
fun1<-d8.l%>%select(plot.id, site_code, block, plot, trt, year_trt, richness, total_cover)%>%distinct()%>%merge(fun, by=c("plot.id"))
# check the traits used 
table(fun1$trait.name)
write.csv(fun1, file="fd and cwm using traits measured at NutNet sites.csv")

###########################################################################################
#################calculate functional diversity for shared species ########################
###########################################################################################
colnames(avg.traits1)
avg.traits1.share<-avg.traits1%>% ungroup()%>% select(site_code, trt, standard_taxon, TraitName, avg.value)%>%
  pivot_wider(names_from = trt, values_from=avg.value)%>%
  filter(!is.na(Control))%>%filter(!is.na(NPK))%>%
  pivot_longer(names_to = "trt", values_to = "avg.value", cols = c("Control", "NPK"))

d8.l.share<-d8.l%>% group_by(site_code,   block, trt,   year_trt , standard_taxon)%>%summarise(max_cover= mean(max_cover))%>%
  pivot_wider(names_from = trt, values_from=max_cover)%>%
  #  filter(!is.na(Control))%>%filter(!is.na(NPK))%>%
  filter(Control > 0)%>%filter(NPK >0)%>%
  pivot_longer(names_to = "trt", values_to = "max_cover", cols = c("Control", "NPK")) %>%
  group_by(site_code, block,  trt, year_trt)%>%mutate(total_cover=sum(max_cover))%>%
  mutate(plot.id=paste(site_code, block,  trt, year_trt, sep="_"))%>%mutate(block.id=paste0(site_code, block, year_trt))%>%
  merge(avg.traits1.share, by=c("site_code", "trt", "standard_taxon"))
length(unique(d8.l.share$standard_taxon))

fun.share<-c()
for(tr in unique(d8.l.share$TraitName)){
  for (pl in unique(d8.l.share$plot.id)){
    # pl<-"bnch.us_1_Control_13" ; tr<-"slow_fast"
    temp.data<-d8.l.share%>%filter(plot.id==pl & TraitName==tr)
    sp.tr<-temp.data$avg.value
    names(sp.tr)<-temp.data$standard_taxon
    abund1<-temp.data$max_cover
    names(abund1)<-temp.data$standard_taxon
    total_cover_traits<-sum(abund1) # total cover with traits
    # 1 species is not possible to calculate fd, but possible to calculate cwm, two species with same trait value is not possible either
    if(length(unique(sp.tr))<2) 
      next
    (fd.temp <- dbFD(sp.tr, abund1, calc.FRic = F, calc.FDiv = F, w.abun=T))##  
    temp.traits<-data.frame(plot.id=pl, f.dispersion=fd.temp$FDis, cwm1=fd.temp$CWM[,1],
                            f.richness=fd.temp$sing.sp, total_cover_traits=total_cover_traits, trait.name=tr)
    fun.share<-rbind(fun.share, temp.traits)
  }
}
fun.share2<-d8.l.share%>%select(plot.id, site_code, block, trt, year_trt, total_cover)%>%distinct()%>%merge(fun.share, by=c("plot.id"))
write.csv(fun.share2, file="fd and cwm using traits measured at NutNet sites for shared species.csv")
fun.share2%>%filter(is.na(cwm1))

# sort the percent cover for species with traits 
colnames(fun1)
percent.coverage<-fun1%>%merge(fun.share2%>%select(trait.name, site_code, block,  trt, year_trt, total_cover_traits)%>%distinct(), by=c("trait.name", "site_code", "block",  "trt", "year_trt"))%>%
  mutate(percent.cover.with.traits0=total_cover_traits.x /total_cover , percent.cover.share.with.traits0=total_cover_traits.y /total_cover )%>%
  group_by(trait.name, site_code, trt)%>%summarise(percent.cover.with.traits=mean(percent.cover.with.traits0), percent.cover.share.with.traits=mean(percent.cover.share.with.traits0))%>%
  mutate(percent.cover.with.traits1=round(percent.cover.with.traits, 2), percent.cover.share.with.traits1=round(percent.cover.share.with.traits, 2))%>%
  select(trait.name, site_code, trt, percent.cover.with.traits1, percent.cover.share.with.traits1)%>%
  pivot_wider(names_from = "trt", values_from = c("percent.cover.with.traits1", "percent.cover.share.with.traits1")) %>%
  mutate(percent.cover.with.traits=paste(percent.cover.with.traits1_Control, percent.cover.with.traits1_NPK, sep="/"), 
         percent.cover.with.traits.shared.species=paste(percent.cover.share.with.traits1_Control, percent.cover.share.with.traits1_NPK, sep="/"))

## using the cover from the NPK and traits from ambient to calculate community traits induced by 
## species turnover and intraspecific trait variation 
fun.share.npk<-c()
for(tr in unique(d8.l.share$TraitName)){
  for (pl in unique(d8.l.share$block.id)){
    # pl<-"bnch.us113" ; tr<-"slow_fast"
    temp.data<-d8.l.share%>%filter(block.id==pl & TraitName==tr & trt=="Control")
    sp.tr<-temp.data$avg.value
    names(sp.tr)<-temp.data$standard_taxon
    temp.data1<-d8.l.share%>%filter(block.id==pl & TraitName==tr & trt=="NPK")
    abund1<-temp.data1$max_cover
    names(abund1)<-temp.data$standard_taxon
    # 1 species is not possible to calculate fd, but possible to calculate cwm, two species with same trait value is not possible either
    if( length(unique(sp.tr))<2) 
      next
    (fd.temp <- dbFD(sp.tr, abund1, calc.FRic = F, calc.FDiv = F, w.abun=T))##  
    temp.traits<-data.frame(block.id=pl, f.dispersion=fd.temp$FDis, cwm1=fd.temp$CWM[,1],
                            f.richness=fd.temp$sing.sp, trait.name=tr)
    fun.share.npk<-rbind(fun.share.npk, temp.traits)
  }
}
colnames(fun.share.npk)
fun.share.npk2<-d8.l.share%>%select("site_code", "block", "trt", "year_trt", "block.id")%>%filter(trt=="NPK")%>%distinct()%>%
  merge(fun.share.npk, by=c("block.id"))%>%mutate(block.id=NULL)
fun.share.npk2%>%filter(is.na(cwm1))

fun.share2.w<-fun.share2%>%select("trait.name", "site_code", "block",   "year_trt",  "trt", "cwm1")%>%distinct()%>%
  pivot_wider(names_from = trt, values_from=cwm1)%>%
  merge(y=fun.share.npk2[,c("trait.name", "site_code", "block",   "year_trt",  "cwm1")], by=c("trait.name", "site_code", "block",   "year_trt"))%>%
  mutate(cwm.control=Control, cwm.npk=NPK, cwm.npk1=cwm1, Control=NULL, NPK=NULL, cwm1=NULL)%>%
  mutate(total=cwm.npk-cwm.control, turn=cwm.npk1-cwm.control, intra=cwm.npk-cwm.npk1)

###################################################################################################
#############################individual trait from global databases ###############################
###################################################################################################
dir.data.global.traits<-"C:/Users/chqq3/work/plant traits compiled for NutNet/"
global.traits.individuals<-read.csv(paste0(dir.data.global.traits, "combining traits at individual level from TRY, BIEN, Aus, and NutNet.csv"))
global.traits.individuals.exclude.nutnet<-global.traits.individuals%>%filter(database!="NutNet")

############### check trait records for species occur at 10 nutnet sites that with on-site trait measurements ###########
all.traits.for.nutnet<-global.traits.individuals.exclude.nutnet%>%filter(standard_taxon%in%d.sites.select$standard_taxon)%>%
  # filter(!is.na(continent))%>% # reduce quite records, maybe no need to filter out data without geolocation
  arrange(standard_taxon, TraitName1)
colnames(all.traits.for.nutnet)

# average over databases 
all.traits.for.nutnet.avg<-all.traits.for.nutnet%>%group_by(standard_taxon, TraitName1)%>%summarise(avg.value=mean(StdValue_avg_across_contributors))
table(all.traits.for.nutnet.avg$TraitName1) # very few records for leaf K
# maybe focus on leaf traits that were also measured in the field of NutNet sites 
# also add LDMC 
sp.level.trait.nutnet<-all.traits.for.nutnet.avg%>%filter(TraitName1 %in% c("Leaf C", "Leaf K",  "Leaf N",  "Leaf P", "SLA", "LDMC"))
length(unique(sp.level.trait.nutnet$standard_taxon))
check.overlap.species<-sp.level.trait.nutnet%>%filter(standard_taxon %in% avg.traits1$standard_taxon)
length(unique(check.overlap.species$standard_taxon))/length(unique(avg.traits1$standard_taxon))
# check how many species have all these five leaf traits
check.species.with.all.traits<-sp.level.trait.nutnet%>% pivot_wider(names_from ="TraitName1", values_from = "avg.value" )%>%
  filter(!if_any(c("Leaf C", "Leaf N", "Leaf P", "Leaf K", "SLA", "LDMC"), is.na))

# focus on leaf traits that were measured in the field of NutNet sites for comparison
all.traits1a<-global.traits.individuals.exclude.nutnet%>%filter(TraitName1 %in% c("Leaf C", "Leaf K",  "Leaf N",  "Leaf P", "SLA", "LDMC"))
all.traits2<-all.traits1a%>%mutate(genus=gsub( " .*$", "", standard_taxon)) # add genus data
unique(all.traits2$TraitName1)


# fill NAs for each trait
all.traits.fill<-c()
for(tr in unique(all.traits2$TraitName1)){
  # tr<-"Leaf N"
  sp.trait<-sp.level.trait.nutnet%>%filter(TraitName1==tr)
  fill.trait<-d.sites.select%>%select(standard_taxon)%>%distinct()%>%filter(!standard_taxon%in%sp.trait$standard_taxon)%>%
    mutate(genus=gsub( " .*$", "", standard_taxon))
  for(g in unique(fill.trait$genus)){
    # g<-"POA"
    # first summarize mean trait value for each species, then...
    original.traits<-all.traits2%>%filter(TraitName1==tr & genus==g)%>%group_by(TraitName1, genus, standard_taxon)%>%
      dplyr::summarise(trait_value_avg_species=mean(StdValue_avg_across_contributors))%>%
      group_by(genus)%>%dplyr::summarise(trait_value_avg=mean(trait_value_avg_species))
    fill.trait1<-fill.trait%>%filter(genus==g)%>%mutate(avg.value=ifelse(nrow(original.traits)==0, NA,  original.traits$trait_value_avg), data.type="average.values.at.genus.level")
    fill.trait1$TraitName1<-tr
    all.traits.fill<-rbind(all.traits.fill, fill.trait1)
  }
}
all.traits.fill1<-all.traits.fill%>%filter(!is.na(avg.value))%>%select(standard_taxon, TraitName1, avg.value, data.type)

# add this data frame to the one with trait data at species level
sp.and.genus.level.traits.nutnet<-sp.level.trait.nutnet%>%mutate(data.type="traits.at.species.level")%>%bind_rows(all.traits.fill1)
## check records for each species 
table(sp.and.genus.level.traits.nutnet$TraitName1, sp.and.genus.level.traits.nutnet$data.type)

###################################################################################################
######################check trait values, distributions, and correlations##########################
###################################################################################################
sp.nutnet1<-sp.and.genus.level.traits.nutnet%>%filter(!is.na(TraitName1))
colnames(sp.nutnet1)
table(sp.nutnet1$TraitName1)
# check data distribution 
ggplot(sp.nutnet1)+geom_histogram(aes(x=avg.value))+facet_wrap(~TraitName1, scale="free_x")

# get traits for these species, delete extreme values by using 5-95% quantile
sp.nutnet2<-sp.nutnet1%>%group_by(TraitName1)%>%filter(avg.value < quantile(avg.value, 0.95))%>%
  filter(avg.value > quantile(avg.value, 0.05))%>%mutate(extreme.values="No")
table(sp.nutnet2$TraitName1)
length(unique(sp.nutnet2$standard_taxon))

# check trait value distribution including and excluding extreme values 
sp.nutnet_12<-sp.nutnet1%>%mutate(extreme.values="Yes")%>%bind_rows(sp.nutnet2)
sp.nutnet_12$extreme.values<-factor(sp.nutnet_12$extreme.values, levels = c("Yes", "No"))
(pp.distribution<-ggplot(sp.nutnet_12)+theme_bw(base_size = 20)+mt+
    geom_density(aes(x=avg.value, color=extreme.values, fill=extreme.values))+
    facet_wrap(~TraitName1, scales="free"))
ggsave(pp.distribution, width=13.3,height=6.64, file="distribution of compiled trait values including and excluding extreme values for 10 nutnet sites.png")

# change data to wide version, focus on traits with high coverage for species 
sp.nutnet.w<-sp.nutnet2%>%mutate(data.type=NULL, extreme.values=NULL)%>%filter(TraitName1 %in%c("LDMC", "SLA", "Leaf N", "Leaf C", "Leaf P", "Leaf K"))%>%
  pivot_wider(names_from = TraitName1, values_from = avg.value)
colnames(sp.nutnet.w)<-c("standard_taxon", "LDMC", "Leaf C", "Leaf K",  "Leaf N",  "Leaf P", "SLA")

(pp1<-ggpairs(sp.nutnet.w[,c("LDMC", "SLA", "Leaf C", "Leaf N", "Leaf P", "Leaf K")])+theme_bw(base_size = 18)+
    labs(title=paste0("pairwise correlation based on ", nrow(sp.nutnet.w), " NutNet species (in total ", n.spp, " species occurred at ", length(unique(d.sites.select$site_code)), " nutnet sites)")))
ggsave(pp1, width=13.3,height=6.64, file="correlation among compiled traits at 10 nutnet sites.png")

## check PCA plot using FIVE leaf traits 
sp.nutnet.w.na.omit<-sp.nutnet.w %>% filter_at(vars("Leaf C", "Leaf N", "Leaf P","Leaf K", "SLA"),all_vars(!is.na(.)))
pca.r<-PCA(sp.nutnet.w.na.omit[,c("SLA", "Leaf C", "Leaf N", "Leaf P", "Leaf K")], scale=T, graph=F)
summary(pca.r)## 
# visualize 
(pp2<-fviz_pca_biplot(pca.r, label = "var", alpha.ind =0.5, pointsize=2, title = paste0("PCA based on ", nrow(sp.nutnet.w.na.omit), " species having all traits (in total ", n.spp, " species occurred)"))+theme_bw(base_size=18)+mt)
ggsave(pp2, width=13.3, height=6.64, file="pca among compiled traits at 10 nutnet sites.png")

# also include other scenarios (delete leaf K; including it will dropped number of species strongly)
sp.nutnet.w.na.omit.without.K<-sp.nutnet.w %>% filter_at(vars("Leaf C", "Leaf N", "Leaf P", "SLA"),all_vars(!is.na(.)))
pca.r.without.K<-PCA(sp.nutnet.w.na.omit.without.K[,c("Leaf C", "Leaf N", "Leaf P", "SLA")], scale=T, graph=F)
summary(pca.r.without.K)## 
(pp2<-fviz_pca_biplot(pca.r.without.K, label = "var", alpha.ind =0.5, pointsize=2, title = paste0("PCA based on ", nrow(sp.nutnet.w.na.omit.without.K), " species having all traits (in total ", n.spp, " species occurred)"))+theme_bw(base_size=18)+mt)

# also include other scenarios (include LDMC)
sp.nutnet.w.na.omit.with.LDMC<-sp.nutnet.w %>% filter_at(vars("Leaf C", "Leaf N", "Leaf P", "SLA", "LDMC"),all_vars(!is.na(.)))
pca.r.with.LDMC<-PCA(sp.nutnet.w.na.omit.with.LDMC[,c("Leaf C", "Leaf N", "Leaf P", "SLA", "LDMC")], scale=T, graph=F)
summary(pca.r.with.LDMC)## 
(pp2<-fviz_pca_biplot(pca.r.with.LDMC, label = "var", alpha.ind =0.5, pointsize=2, title = paste0("PCA based on ", nrow(sp.nutnet.w.na.omit.with.LDMC), " species having all traits (in total ", n.spp, " species occurred)"))+theme_bw(base_size=18)+mt)
ggsave(pp2, width=11,height=6.64, file="pca among compiled traits at 10 nutnet sites.png")

# including all these traits reduces number of species strongly. 
(loadings <- as.data.frame(pca.r$ind$coord))
(loadings_without_K <- as.data.frame(pca.r.without.K$ind$coord))
(loadings_with_LDMC <- as.data.frame(pca.r.with.LDMC$ind$coord))
# better do not include leaf K; this trait has few records from global databases
# sp.nutnet.w.na.omit.without.K$slow_fast<-loadings$Dim.1
sp.nutnet.w.na.omit.without.K$slow_fast_without_K<-loadings_without_K$Dim.1
sp.nutnet.w.na.omit.with.LDMC$slow_fast_with_LDMC<-loadings_with_LDMC$Dim.1

# add integrated traits from pac to other individual traits 
sp.nutnet3<-sp.nutnet.w.na.omit.without.K%>%select(standard_taxon, slow_fast_without_K)%>%
  mutate(data.type="pca", TraitName1="slow_fast", avg.value=slow_fast_without_K)%>%
  select(standard_taxon, TraitName1, avg.value, data.type)%>%
  bind_rows(sp.nutnet.w.na.omit.with.LDMC%>%select(standard_taxon, slow_fast_with_LDMC)%>%mutate(data.type="pca", TraitName1="slow_fast_with_LDMC", avg.value=slow_fast_with_LDMC)%>%select(standard_taxon, TraitName1, avg.value, data.type) )%>%
  bind_rows(sp.nutnet2%>%mutate(extreme.values=NULL)%>%filter(TraitName1 %in%c("LDMC", "SLA", "Leaf N", "Leaf C", "Leaf P")))
colnames(sp.nutnet3); colnames(sp.level.trait.nutnet)

sp.nutnet4<-sp.nutnet3%>%mutate(data.type1="species.and.genus")%>%
  bind_rows(sp.nutnet2%>%mutate(extreme.values=NULL)%>%filter(data.type=="traits.at.species.level")%>%filter(TraitName1 %in%c("LDMC", "SLA", "Leaf N", "Leaf C", "Leaf P"))%>%mutate(data.type1="species only"))

# sort the trait information in a table 
trait.inf1<-sp.nutnet4%>%group_by(data.type1, TraitName1)%>%summarise(percent.species.with.traits=length(avg.value)/n.spp)
trait.inf2<-sp.nutnet4%>%group_by(data.type1, TraitName1, data.type)%>%summarise(percent.species.with.traits=length(avg.value)/n.spp)%>%
  pivot_wider(names_from = data.type, values_from = percent.species.with.traits)
trait.inf4<-trait.inf1%>%merge(y=trait.inf2, by=c("data.type1", "TraitName1"))%>%arrange(data.type1, -percent.species.with.traits)

trait.inf4[,3:6]<-lapply(trait.inf4[,3:6], function(x){x*100} )
trait.inf4[,3:6]<-round(trait.inf4[,3:6], 0)
colnames(trait.inf4)<-c("Data type", "TraitName", "percent species with traits",   "percent species with traits at species level",                               
                        "percent species with traits at genue level",      "percent species with traits for LES" )              
write.csv(trait.inf4, file="summary of compiled trait information at 10 nutnet sites.csv")

###################################################################################################
#############################cwm and fd  at the plot level over time###############################
###################################################################################################
# add species with traits to the cover data 
# focus on cover data from year_trt 1, 
colnames(d.sites.select);colnames(sp.nutnet3)
d.sites.select.l<-d.sites.select%>%filter(max_cover>0)%>%select(site_code, trt, standard_taxon, block, plot,trt, year_trt , year,  max_cover)%>%unique()%>%
  group_by(site_code, block,plot, trt, year_trt)%>%mutate(total_cover=sum(max_cover))%>%
  filter(year_trt>0 & trt %in% c("Control", "NPK"))%>%
  mutate(plot.id=paste0(site_code, block, plot, trt, year_trt))%>%
  merge(sp.nutnet4, by=c("standard_taxon"))
table(d.sites.select.l$TraitName1)

hist(d.sites.select.l$total_cover)
check.na.t<-d.sites.select.l[is.na(d.sites.select.l$avg.value),]
check.na.c<-d.sites.select.l[is.na(d.sites.select.l$max_cover),]

fun.global<-c()
for(dt in unique(d.sites.select.l$data.type1)){
  for(tr in unique(d.sites.select.l$TraitName1)){
    for (pl in unique(d.sites.select.l$plot.id)){
      # pl<-"azi.cn11CControl1"; tr<-"Leaf N"
      temp.data<-d.sites.select.l%>%filter(data.type1==dt & plot.id==pl & TraitName1==tr)
      sp.tr<-temp.data$avg.value
      names(sp.tr)<-temp.data$standard_taxon
      abund1<-temp.data$max_cover
      names(abund1)<-temp.data$standard_taxon
      total_cover_traits<-sum(abund1) # total cover with traits
      # 1 species is not possible to calculate fd, but possible to calculate cwm, two species with same trait value is not possible either
      if(length(unique(sp.tr))<2) 
        next
      (fd.temp <- dbFD(sp.tr, abund1, calc.FRic = F, calc.FDiv = F, w.abun=T))##  
      temp.traits<-data.frame(plot.id=pl, f.dispersion=fd.temp$FDis, cwm1=fd.temp$CWM[,1],
                              f.richness=fd.temp$sing.sp, total_cover_traits=total_cover_traits, trait.name=tr, data.type1=dt)
      fun.global<-rbind(fun.global, temp.traits)
    }
  }
}
colnames(fun.global)
fun.global1<-d.sites.select.l%>%select(site_code, block, plot, trt, year_trt, plot.id, total_cover)%>%distinct()%>% merge(fun.global, by=c("plot.id"))
# check the traits used 
table(fun.global1$trait.name)
write.csv(fun.global1, file="fd and cwm based on compiled traits.csv")

# sort the percent cover for species with traits 
percent.coverage.global<-fun.global1%>%
  mutate(percent.cover.with.traits0=total_cover_traits /total_cover)%>%
  group_by(data.type1, trait.name, site_code, trt)%>%summarise(percent.cover.with.traits=mean(percent.cover.with.traits0))%>%
  mutate(percent.cover.with.traits1=round(percent.cover.with.traits, 2))%>%
  select(data.type1, trait.name, site_code, trt, percent.cover.with.traits1)%>%
  pivot_wider(names_from = "trt", values_from = c("percent.cover.with.traits1")) %>%
  mutate(percent.cover.with.traits=paste(Control, NPK, sep="/"))
write.csv(percent.coverage.global, file="summary of percent species covered by compiled traits at 10 nutnet sites.csv")

###########################################################################################
#################treatment effects on community traits and richness########################
#######################calculate cwm and fd based on shared species########################
colnames(fun.share2.w)
unique(fun.share2.w$trait.name)
fun.share2.w1.l<-fun.share2.w%>%
  group_by(trait.name, site_code, block)%>%dplyr::summarise_at(c("cwm.control", "cwm.npk", "cwm.npk1"), list(mean))%>%
  select(trait.name, site_code, block,cwm.control, cwm.npk, cwm.npk1)%>%
  pivot_longer(cols = cwm.control:cwm.npk1)%>%
  mutate(trt=gsub("cwm.", "", name))%>% 
  mutate(trt1=ifelse(trt=="npk", "npk2", trt))%>%
  mutate(total=ifelse((name=="cwm.npk"|name=="cwm.control"), 1, 0), turn=ifelse((name=="cwm.npk1"|name=="cwm.control"), 1, 0), intra=ifelse((name=="cwm.npk1"|name=="cwm.npk"), 1, 0))
colnames(fun.share2.w1.l)

tt<-data.frame()
for (tr in unique(fun.share2.w1.l$trait.name)){ 
  for (i in c("total", "intra", "turn")){
    # i<-"total"; tr<- "Leaf N" 
    d.e<-fun.share2.w1.l%>%filter(trait.name==tr)%>%select(trait.name, site_code, block,   name,   value,    trt1, i)
    ## rename
    colnames(d.e)[7]<-"variable"
    d.e1<-d.e%>%filter(variable==1)
    mod<-lme(value ~ trt1, random= ~1|site_code/block, data=d.e1, na.action = na.omit)
    t<-xtable(summary(mod)$tTable)
    rs<-rsquared(mod)
    t$r2m<-rs$Marginal
    t$r2c<-rs$Conditional
    t$terms<-rownames(t)
    t$variable<-i
    t$trait.name<-tr
    tt<-rbind(tt, t)
  }
}
tt[,1:7]<-lapply(tt[,1:7], function(x) round(x, 3))

# functional diversity using shared species 
colnames(fun.share2)
fun.share3<-fun.share2%>%
  group_by(trait.name, site_code, trt, block)%>%dplyr::summarise_at(c("cwm1", "f.dispersion"), list(mean))%>%
  ungroup()%>%mutate(cwm1=NULL, ind="f.dispersion.share", values=f.dispersion)%>%select(trait.name, site_code, trt, block, ind, values)

# calculate cwm and fd using all species with trait data
fun2.l<-fun1%>%
  group_by(trait.name, site_code, trt, block)%>%dplyr::summarise_at(c("cwm1", "f.dispersion"), list(mean))%>%
  pivot_longer(cols=c("cwm1", "f.dispersion"), names_to="ind", values_to="values")%>%arrange(ind)

# calculate cwm and fd using trait data from global databases 
fun.global2.l<-fun.global1%>%
  group_by(data.type1, trait.name, site_code, trt, block)%>%dplyr::summarise_at(c("cwm1", "f.dispersion"), list(mean))%>%
  pivot_longer(cols=c("cwm1", "f.dispersion"), names_to="ind", values_to="values")%>%arrange(data.type1, ind)%>%mutate(ind=paste0("g_", ind, "_", data.type1))%>%ungroup()%>%select(- data.type1)

## species richness 
## select sites and log transform the diversity data because richness is not on the same scale with functional traits
rich<-d.sites.select%>%select(site_code, trt, block, year, richness)%>%distinct()%>%
  group_by(site_code, trt, block)%>%summarise(diversity=mean(richness))%>%
  mutate(values=log(diversity), ind="richness", trait.name="")%>%select(trait.name, site_code, trt, block, ind, values)
# all them up
colnames(fun2.l)
colnames(fun.global2.l)
div.facet<-fun.share3%>%bind_rows(fun2.l)%>%bind_rows(fun.global2.l)%>%bind_rows(rich)%>%mutate(ind1=paste0(trait.name, ind))

## perhaps run the models in one go 
tt.a<-data.frame()
for (i in unique(div.facet$ind1)){
  # i<-"cwm1" ; tr<- "Leaf N" 
  d.e<-div.facet%>%filter(ind1==i)
  mod<-lme(values ~ trt, random= ~1|site_code/block, data=d.e, na.action = na.omit)
  t<-xtable(summary(mod)$tTable)
  rs<-rsquared(mod)
  t$r2m<-rs$Marginal
  t$r2c<-rs$Conditional
  t$terms<-rownames(t)
  t$ind1<-i
  tt.a<-rbind(tt.a, t)
}
tt.a[,1:7]<-lapply(tt.a[,1:7], function(x) round(x, 3))
tt.a1<-div.facet%>%select(ind1, trait.name, ind)%>%distinct()%>%merge(tt.a, by=c("ind1"))%>%
  select(Value, Std.Error, DF, 't-value', 'p-value', r2m, r2c, terms, ind, trait.name)%>%dplyr::rename(variable=ind)

####################################### add all model output ######################################
colnames(tt.a1); colnames(tt)
unique(tt$trait.name)
## delete intercept 
tt2<-tt%>%bind_rows(tt.a1)%>%filter(terms!="(Intercept)")%>%mutate(lower=Value-1.96*Std.Error, upper=Value+1.96*Std.Error)%>%
  mutate(variable1=case_when(variable %in% c("total", "cwm1", "g_cwm1_species only", "g_cwm1_species.and.genus" )~"CWM ", 
                             variable %in% c("intra")~"CWM induced by ITS",
                             variable %in% c("turn")~"CWM induced by replacement",
                             variable %in% c("f.dispersion", "f.dispersion.share", "g_f.dispersion_species only", "g_f.dispersion_species.and.genus" )~"FD",
                             variable %in% c("richness")~"Log (richness)"))%>%
  mutate(species=case_when(variable %in% c("total", "intra", "turn", "f.dispersion.share")~"Shared species with traits", 
                           variable %in% c("cwm1", "f.dispersion")~"Species with traits", 
                           variable %in% c("g_cwm1_species only", "g_f.dispersion_species only" )~"Species with traits extracted (species only)",
                           variable %in% c("g_cwm1_species.and.genus", "g_f.dispersion_species.and.genus" )~"Species with traits extracted (species and genus)",
                           variable %in% c("richness")~"All species"))%>%
  mutate(trait.name1=ifelse((trait.name %in% c("slow_fast")), "LES", trait.name))%>%
  mutate(trait.name2=ifelse((trait.name1 %in% c("slow_fast_without_K")), "LES without K", trait.name1))%>%
  mutate(trait.name3=ifelse((trait.name2 %in% c("slow_fast_with_LDMC")), "LES with LDMC", trait.name2))

## relevel
tt2$species<-factor(tt2$species, levels=c("All species", "Species with traits", "Shared species with traits", "Species with traits extracted (species only)", "Species with traits extracted (species and genus)"))
tt2$trait.name3<-factor(tt2$trait.name3, levels=c("LES", "LES without K", "LES with LDMC", "Leaf C", "Leaf K", "Leaf N", "Leaf P", "SLA", ""))
## better plot all variables together in one figure

## adding data points from individual sites
# for shared species with traits
div.facet.1<-div.facet%>%
  pivot_wider(names_from=trt, values_from=values)%>%mutate(dif=NPK-Control)%>% select(trait.name, site_code, block, ind, dif)
# add all together 
raw.div<-fun.share2.w%>%group_by(trait.name, site_code, block)%>%summarise_at(c("total", "intra", "turn"), list(mean))%>%
  pivot_longer(names_to="ind", values_to="dif", cols = total:turn)%>%select(trait.name, site_code, block, ind, dif)%>%arrange(ind)%>%
  bind_rows(div.facet.1)%>%ungroup()%>%
  mutate(variable1=case_when(ind %in% c("total", "cwm1", "g_cwm1_species only", "g_cwm1_species.and.genus" )~"CWM ", 
                             ind %in% c("intra")~"CWM induced by ITS",
                             ind %in% c("turn")~"CWM induced by replacement",
                             ind %in% c("f.dispersion", "f.dispersion.share", "g_f.dispersion_species only", "g_f.dispersion_species.and.genus")~"FD",
                             ind %in% c("richness")~"Log (richness)"))%>%
  mutate(species=case_when(ind %in% c("total", "intra", "turn", "f.dispersion.share")~"Shared species with traits", 
                           ind %in% c("cwm1", "f.dispersion")~"Species with traits", 
                           ind %in% c("g_cwm1_species only", "g_f.dispersion_species only" )~"Species with traits extracted (species only)",
                           ind %in% c("g_cwm1_species.and.genus", "g_f.dispersion_species.and.genus" )~"Species with traits extracted (species and genus)",
                           ind %in% c("richness")~"All species"))%>%
  mutate(trait.name1=ifelse((trait.name %in% c("slow_fast")), "LES", trait.name))%>%
  mutate(trait.name2=ifelse((trait.name1 %in% c("slow_fast_without_K")), "LES without K", trait.name1))%>%
  mutate(trait.name3=ifelse((trait.name2 %in% c("slow_fast_with_LDMC")), "LES with LDMC", trait.name2))

## relevel
raw.div$species<-factor(raw.div$species, levels=c("All species", "Species with traits", "Shared species with traits", "Species with traits extracted (species only)", "Species with traits extracted (species and genus)"))
raw.div$trait.name3<-factor(raw.div$trait.name3, levels=c("LES", "LES without K", "LES with LDMC", "Leaf C", "Leaf K", "Leaf N", "Leaf P", "SLA", ""))
unique(raw.div$variable1)
unique(tt2$variable1)
trait.interest<-c("", "slow_fast", "Leaf C",  "Leaf N",  "Leaf P",  "Leaf K", "SLA")
tt2.measured<-tt2%>%filter(species %in% c("Species with traits", "Shared species with traits"))%>%filter(trait.name %in% trait.interest)%>%
  filter(!(variable1=="FD" & species=="Shared species with traits"))%>%
  mutate(diversity.type=ifelse(variable1=="FD", "trait.diversity", "trait.mean"))

raw.div.measured<-raw.div%>%filter(species %in% c("Species with traits", "Shared species with traits"))%>%filter(trait.name %in% trait.interest)%>%
  filter(!(variable1=="FD" & species=="Shared species with traits"))%>%
  mutate(diversity.type=ifelse(variable1=="FD", "trait.diversity", "trait.mean"))
tt2.measured$trait.name3<-factor(tt2.measured$trait.name3, levels = c("LES", "Leaf C",  "Leaf N",  "Leaf P",  "Leaf K", "SLA"))
raw.div.measured$trait.name3<-factor(raw.div.measured$trait.name3, levels = c("LES", "Leaf C",  "Leaf N",  "Leaf P",  "Leaf K", "SLA"))

# adjust decimal in y axis  
scaleFUN <- function(x) sprintf("%.1f", x)
raw.div.measured.fd<- raw.div.measured%>%filter(diversity.type %in% c("trait.diversity"))
(pp.div1<-tt2.measured %>%filter(diversity.type %in% c("trait.diversity"))%>%
    ggplot(aes(variable1, Value))+theme_bw(base_size = 12)+mt+
    geom_violin(raw.div.measured.fd, mapping = aes(variable1, dif), color="light grey")+
    geom_point(raw.div.measured.fd, mapping = aes(variable1, dif), size=0.5, alpha=0.15)+
    geom_point(position=pd, size=2, pch=19)+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, position=pd) +
    facet_grid(trait.name3 ~ species, switch = "y", space = "free_x")+
    geom_hline(yintercept = 0, linetype ="dotted") +
    scale_y_continuous(labels=scaleFUN)+ 
    theme(legend.position = "none", axis.text.x = element_text(angle=8, vjust = 1.05, hjust=1), strip.placement = "outside", strip.background.y = element_blank())+
    labs(x=NULL, y=NULL))
raw.div.measured.trait.mean<- raw.div.measured%>%filter(diversity.type %in% c("trait.mean"))
(pp.div2<-tt2.measured %>%filter(diversity.type %in% c("trait.mean"))%>%
    ggplot(aes(variable1, Value))+theme_bw(base_size = 12)+mt+
    geom_violin(raw.div.measured.trait.mean, mapping = aes(variable1, dif), color="light grey")+
    geom_point(raw.div.measured.trait.mean, mapping = aes(variable1, dif), size=0.5, alpha=0.15)+
    geom_point(position=pd, size=2, pch=19)+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, position=pd) +
    facet_grid(trait.name3 ~ species, scales="free", switch = "y", space = "free_x")+
    geom_hline(yintercept = 0, linetype ="dotted") +
    scale_y_continuous(labels=scaleFUN)+ 
    theme(legend.position = "none", axis.text.x = element_text(angle=8, vjust = 1.05, hjust=1), strip.placement = "outside", strip.background.y = element_blank())+
    labs(x=NULL, y=NULL))
ggsave(pp.div1, width = 5, height = 14.5, unit="cm",  file="treatment effects on trait diversity.svg", device = "svg")
ggsave(pp.div2, width = 13, height = 15,  unit="cm",  file="treatment effects on trait mean.svg", device = "svg")


## save the model output 
ttt<-tt%>%bind_rows(tt.a1)%>%
  mutate(variable1=case_when(variable %in% c("total", "cwm1", "g_cwm1_species only", "g_cwm1_species.and.genus" )~"CWM ", 
                             variable %in% c("intra")~"CWM induced by ITS",
                             variable %in% c("turn")~"CWM induced by replacement",
                             variable %in% c("f.dispersion", "f.dispersion.share", "g_f.dispersion_species only", "g_f.dispersion_species.and.genus" )~"FD",
                             variable %in% c("richness")~"Log (richness)"))%>%
  mutate(species=case_when(variable %in% c("total", "intra", "turn", "f.dispersion.share")~"Shared species with traits", 
                           variable %in% c("cwm1", "f.dispersion")~"Species with traits", 
                           variable %in% c("g_cwm1_species only", "g_f.dispersion_species only" )~"Species with traits extracted (species only)",
                           variable %in% c("g_cwm1_species.and.genus", "g_f.dispersion_species.and.genus" )~"Species with traits extracted (species and genus)",
                           variable %in% c("richness")~"All species"))%>%
  mutate(trait.name1=ifelse((trait.name %in% c("slow_fast")), "LES", trait.name))%>%
  mutate(trait.name2=ifelse((trait.name1 %in% c("slow_fast_without_K")), "LES without K", trait.name1))%>%
  mutate(trait.name3=ifelse((trait.name2 %in% c("slow_fast_with_LDMC")), "LES with LDMC", trait.name2))%>%
  mutate(terms1=ifelse(terms %in% c("(Intercept)"), "(Intercept)", "trtNPK"))%>%
  filter(!(variable1=="FD" & species=="Shared species with traits"))%>%
  arrange(species, trait.name3)%>%
  select("species", "trait.name3", "variable1", "terms1",  "Value",     "Std.Error", "DF",        "t-value",   "p-value",   "r2m",       "r2c")

col<-str_to_title(colnames(ttt))
colnames(ttt)<-col
ttt[,5:11]<-lapply(ttt[,5:11], function(x) round(x, 2))
ttt.sig<-ttt%>%filter(`P-Value` <= 0.1)%>%filter(Terms1!="(Intercept)")%>%arrange(Variable1, Species, Trait.name3)
write.csv(ttt, file="treatment effects on cwm and fd.csv")

## save the average data of CWM and its partition for each site
unique(raw.div$ind)
unique(raw.div$species)
unique(raw.div$trait.name3)

cwm.and.its.partition<-raw.div%>%  filter(ind%in% c("cwm1", "total",  "turn", "intra"))%>%
  filter(species %in% c("Shared species with traits", "Species with traits"))%>%
  group_by(trait.name1, site_code, ind)%>%summarise(dif.avg=mean(dif))%>%
  select(trait.name1, site_code, ind, dif.avg)%>%ungroup()%>%pivot_wider(names_from = ind, values_from = dif.avg)%>%
  select(trait.name1, site_code, cwm1, total, turn, intra)%>%filter(trait.name1 %in% c("LES"))
cwm.and.its.partition[,3:6]<-round(cwm.and.its.partition[,3:6], 2)
colMeans(cwm.and.its.partition[,3:6])
with(cwm.and.its.partition, t.test(total, cwm1, paired = T)) # marginal significant 

fd.all.and.share<- raw.div%>%  filter(ind%in% c("f.dispersion",  "f.dispersion.share"))%>%
  filter(species %in% c("Shared species with traits", "Species with traits"))%>%
  group_by(trait.name1, site_code, ind)%>%summarise(dif.avg=mean(dif))%>%
  select(trait.name1, site_code, ind, dif.avg)%>%ungroup()%>%pivot_wider(names_from = ind, values_from = dif.avg)%>%
  filter(trait.name1 %in% c("LES"))
colMeans(fd.all.and.share[,3:4])
with(fd.all.and.share, t.test(f.dispersion.share, f.dispersion, paired = T)) # no difference 

###########################################################################################
######################### sort out traits data in a table #################################
###########################################################################################
## check how many species have trait data for these 10 sites in ambient and nutrient addition
ns4<-avg.traits%>%group_by(site_code, treatment)%>%dplyr::summarise(N=length(unique(standard_taxon)))%>%pivot_wider(names_from = treatment, values_from=N)%>%
  mutate(richness.ambient.nutrient=paste(Ambient, `Nutrient addition`, sep="/"))
## check how many species measured for traits at each site 
ns6<-avg.traits%>%group_by(site_code)%>%dplyr::summarise(N=length(unique(standard_taxon)))
range(ns6$N)## 
ns7<-ns4%>%select(site_code, richness.ambient.nutrient)%>%bind_cols(ns6[,2])%>%dplyr::rename(total.species.with.traits=N)## many species found in both treatments
## add years in which traits were measured 
## add proportion data 
site.infor1<-ns7%>%
  merge(y=percent.coverage%>%filter(trait.name %in% c("slow_fast"))%>%ungroup()%>%select("site_code", "percent.cover.with.traits", "percent.cover.with.traits.shared.species"), by=c("site_code"))
# add CWM and its partition at each site 
site.infor2<-site.infor1%>%merge(y=cwm.and.its.partition, by=c("site_code"))
site.infor3<-site.infor2%>%
  select("site_code",  "richness.ambient.nutrient",  "percent.cover.with.traits" ,  "percent.cover.with.traits.shared.species" , "cwm1",   "total",   "turn",  "intra")            
write.csv(site.infor3, file="trait information at each site.csv")

# the end 