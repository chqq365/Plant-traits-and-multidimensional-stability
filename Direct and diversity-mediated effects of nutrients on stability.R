## use sites bnch.us, bogong.au, burrawan.au, cbgb.us, cowi.ca, hopl.us, look.us, mcla.us, sgs.us, sier.us

rm(list=ls())
## open the libraries
library(tidyverse);library(xtable)
library(ggplot2);library(cowplot);library(scales);library(ggthemes)
library(nlme);library(piecewiseSEM);library(car);library(AICcmodavg)

## environment needed
pd<-position_dodge(width=0.5)
mt<-theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
dir.data<-"C:/Users/chqq3/work/NutNet data/"
dir.graphs<-"C:/Users/chqq3/work/traits and stability/graphs/"
setwd(dir.graphs)

# select sites 
d7<-read.csv(file="cover and biomass data from 10 NutNet sites.csv")
d8<-d7%>%filter(site_code %in% c("bnch.us", "bogong.au", "burrawan.au", "cbgb.us", "cowi.ca", "hopl.us", "look.us", "mcla.us", "sgs.us", "sier.us"))%>%distinct()
unique(d8$site_code)
range(d8$year)
## find the pre-treatment and end of the experiments
## note that we have data from look.us from year_trt 2!!!!
start.end.year<-d8%>%group_by(site_code) %>% summarise_at("year", list(max, min))%>%
  mutate(min=ifelse(site_code=="look.us", fn2-2, fn2-1))%>%select(-fn2)
colnames(start.end.year)<-c("site_code", "site_lyear","site_fyear" )# site_lyear already include the pre-treatment year

# growing season for each site using data from Sid 
# the harvest month refer to the end of growing season in Sid's data 
gs<-read.csv(file="growing seasons at 10 NutNet sites.csv")
gs4<-gs%>%dplyr::select(site_code, gs_start_month, harvest_month)%>%
  filter(site_code %in% d8$site_code)%>%mutate(harvest_month1=harvest_month)%>%
  merge(y=start.end.year, by=c("site_code"))%>%select(site_code, site_fyear, site_lyear, gs_start_month, harvest_month)
#################################################################################################
###################################### water balance data #######################################
#################################################################################################
## get the data of precipitation and potential evaporation
prec <- read.csv(paste0(dir.data, "CRU-monthly-pre-1901-2022.csv"), header = T) %>%dplyr::rename(prec = value, year=years, month=months)
table(prec$year)
pet<-read.csv(paste0(dir.data, "CRU-monthly-pet-1901-2022.csv"), header = T, sep=",") %>%dplyr::rename(pet = value, year=years, month=months)
### note that potential evaporation was calculated as average per day. 
month.df <- data.frame(month = c("Jan","Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),
                       month1 = seq(1, 12, 1),
                       days = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))

dd0 <- prec %>% merge(y = pet, by=c("site_code", "year", "month", "date")) %>%
  merge(y = month.df, by=c("month"))  %>%
  mutate(pet1 = pet * days) %>%
  arrange(site_code, year, month1)

## add the growing season data 
dd<-merge(dd0, gs4[,c("site_code", "gs_start_month", "harvest_month")], by="site_code")
unique(dd$site_code)
##
# note, at some sites, growing season start in Oct. or Nov. and peak
# biomass season was in the next year
# do it separately for these two situations
dd1 <- dd %>% mutate(dif = harvest_month - gs_start_month, one.v.two = ifelse(dif < 0, "two", "one"))
check<-dd1%>%select(site_code, one.v.two)%>%distinct()
table(check$one.v.two)
range(dd1$year)
## summarize growing season precipitation and evapotranspiration for sites with growing season within one year
colnames(dd1)
s_dd2.within <- dd1 %>%
  filter(one.v.two == "one") %>%filter(month1>=gs_start_month & month1 <= harvest_month)%>%
  group_by(site_code, year) %>% 
  summarise_at(c("prec", "pet1"), list(sum))
unique(s_dd2.within$site_code) ##
range(s_dd2.within$year)

## summarize growing season precipitation and evapotranspiration for sites with growing season across two years
## months in the late season of the previous year should be counted to the next year in which peak biomass is recorded!!!! 
## in this case, data in 1901 would be incomplete for the peak biomass in that year, because the starting 
## growing months in the previous year (1900) were not available. Thus, delete data in 1901.
## summarize growing season precipitation and evapotranspiration during growing season 
s_dd2.across <- dd1 %>%
  filter(one.v.two == "two")%>%filter(month1>=gs_start_month | month1 <= harvest_month)%>%
  mutate(year1 = ifelse(month1 >= gs_start_month, year+1, year)) %>%
  filter(year1 != 1901) %>% 
  group_by(site_code, year1) %>% 
  summarise_at(c("prec", "pet1"), list(sum))
unique(s_dd2.across$site_code) ##
range(s_dd2.across$year1)

## for sites where growing season is within a year, whole growing season was recorded in 1901,  
## to be consistent, delete data in 1901 
## merge all the sites 
# Compute water balance (BAL) and standardize it 
s_dd3 <- s_dd2.within %>% filter(year!=1901) %>% dplyr::rename(year1=year) %>% bind_rows(s_dd2.across) %>%
  mutate(dif.re0 = prec - pet1) %>% group_by(site_code) %>%
  mutate(dif.re = scale(dif.re0))
length(unique(s_dd3$site_code))
## check the means 
che<-s_dd3%>%group_by(site_code)%>%summarise_at("dif.re", list(mean, sd)) 
range(che$fn1);range(che$fn2)
####calculate water balance from 2005 to 2020 for each site
water.balance<-s_dd3%>%filter(year1>=2005 & year1<=2022)%>%group_by(site_code)%>%dplyr::summarise(avg.wb=mean(dif.re))

#################################################################################################
##########################moderate and extreme dry and wet climate events #######################
#################################################################################################
## add the start and end of the experiments 
s_dd5<-s_dd3%>%merge(y=gs4[,c("site_code", "site_fyear", "site_lyear")], by=c("site_code"))%>%filter(year1>=site_fyear & year1<=site_lyear)
### define climate events based on 0.67 
s_dd5$climate1<-"Normal"
s_dd5$climate1[s_dd5$dif.re> 1.28]<-"Extreme wet"
s_dd5$climate1[s_dd5$dif.re< -1.28]<-"Extreme dry"
s_dd5$climate1[s_dd5$dif.re> 0.67 & s_dd5$dif.re<1.28]<-"Moderate wet"
s_dd5$climate1[s_dd5$dif.re< -0.67 & s_dd5$dif.re>-1.28]<-"Moderate dry"
## check whether all sites have extreme and normal years 
s_dd5_n<-s_dd5%>%filter(climate1=="Normal" & year1>site_fyear)%>%select(site_code)%>%distinct()# site_lyear already include the pre-treatment year
s_dd5_e<-s_dd5%>%filter(climate1 %in% c("Extreme wet", "Extreme dry"))%>%filter(year1>site_fyear)%>%select(site_code)%>%distinct()
## find sites without normal years or extreme years 
ss1<-unique(s_dd5_n$site_code)
ss2<-unique(s_dd5_e$site_code)
ss3<-levels(factor(dd1$site_code))
(sss1 <- ss3[which(!ss3 %in% ss1)])## 
(sss2 <- ss3[which(!ss3 %in% ss2)])## 
## bind these two data sets
(sss3<-c(sss1, sss2))
## delete sites without extreme or normal events, also add year_trt
s_dd6<-s_dd5%>%filter(!site_code%in%sss3)%>%mutate(year_trt=year1-site_fyear)%>%dplyr::rename(spei=dif.re, year=year1, climate=climate1)
colnames(s_dd6)
# get sites and years with biomass data and indicate missing biomass data 
intact.bio<-d8%>%select(site_code, year_trt)%>%distinct()%>%mutate(bio.data="Intact")
s_dd7<-s_dd6%>%merge(y=intact.bio, by=c("site_code", "year_trt"), all.x=T)%>%
  mutate(bio.data1=ifelse(is.na(bio.data), "Missing", bio.data))%>%
  mutate(bio.data2=ifelse(year_trt==0, "Intact", bio.data1))%>%arrange(site_code)
# relevel climate events
s_dd7$climate<-factor(s_dd7$climate, levels=c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))
(pp<-ggplot(s_dd7, aes(year_trt, spei, color=climate, shape=bio.data2))+geom_point(size=4)+theme_bw(base_size = 15)+
    facet_wrap(~site_code, scale="free_x", ncol=5)+
    scale_x_continuous(expand=c(0,0), limits = c(-0.5,15.5), breaks = seq(0,16,2))+
    scale_color_manual(values = c("#E69F00", "grey88", "#000000", "grey88", "#0072B2"))+
    scale_shape_manual(values=c(19, 17))+
    geom_vline(xintercept = 0, linetype="dashed")+
    labs(x=NULL, y="SPEI", shape="Biomass data", color="Climate extremes"))
ggsave(pp, width=13.3, height=6.64, dpi=600, file="moderate and extreme climate events.png")
table(s_dd7$climate)

###########################################################################################
#####biomass under dry and wet events in each year at each site using cutoff of 1.28sd
###########################################################################################
d.select<-s_dd7%>%select("site_code", "year", "spei", "climate")%>%distinct()%>%
  merge(y=d7%>%select("site_code", "year", "year_trt", "block", "trt", "live_mass")%>%distinct(), by=c("site_code", "year"))
# calculate the mean biomass 
d.select_normal_avg<-d.select%>%filter(climate=="Normal")%>%group_by(trt)%>%dplyr::summarise(avg_bio=mean(live_mass))
d.select_normal_avg$avg_bio
d.select_site_avg<-d.select%>%filter(climate=="Normal")%>%group_by(site_code, trt)%>% dplyr::summarise(avg_bio=mean(live_mass))

## plot biomass over time at each site
(pp.each.site<-ggplot(d.select, aes(year_trt, live_mass, color=climate, shape=trt))+theme_bw(base_size = 16)+mt+
    facet_wrap(~site_code, scales="free_y", nrow=2)+
    stat_summary(fun=mean, geom="point", size=2, position=pd, alpha=0.6)+
    stat_summary(fun.data =mean_cl_boot, geom="errorbar", width=0.1, position=pd)+
    geom_hline(d.select_site_avg, mapping=aes(yintercept =avg_bio, linetype=trt), alpha=0.6)+
    scale_color_manual(values = c("#E69F00", "grey88", "#000000", "grey88", "#0072B2"))+
    scale_x_continuous(expand=c(0,0), limits = c(0,16), breaks = seq(0,16,2))+
    theme(axis.text.y = element_text(angle=45), panel.spacing.x=unit(-0.2, "lines"))+
    labs(x=NULL, y="Aboveground biomass", color="Climate\n extremes", shape="Treatments", linetype="Treatments"))
ggsave(pp.each.site, width=13.3, height=6.64, dpi=600, file="above ground biomass at each site using cutoff of 1.28.png")

###########################################################################################
# biomass under different climate events under and one year after climate events using cutoff of 1.28sd
###########################################################################################
## select biomass During 
d.select.u1<-d.select%>% filter(grepl("Extreme", climate))%>%select("site_code", "trt", "block", "year_trt", "climate", "live_mass")%>%distinct()%>%
  mutate(events="During")
## summarize biomass one year after 
## add biomass During 
## calculate biomass deviation from the normal level for each site each treatment 

d.select_block_avg<-d.select%>%filter(climate=="Normal")%>%select("site_code", "trt", "block", "year_trt", "climate", "live_mass")%>%distinct()%>%
  group_by(site_code, block, trt)%>% dplyr::summarise(avg_bio=mean(live_mass))
# add average biomass During different treatments different blocks 
# calculate average across sites 
d.select_block_sum<-d.select_block_avg%>%group_by(trt)%>%summarise(avg.mass=mean(avg_bio))
(d.select_block_sum$avg.mass[d.select_block_sum$trt=="NPK"]-d.select_block_sum$avg.mass[d.select_block_sum$trt=="Control"])/d.select_block_sum$avg.mass[d.select_block_sum$trt=="Control"]

# when two climate extremes happen consecutively, recovery only calculate for the later year
# climate events happen in the last year of a site cannot calculate recovery because no data the year after
# sites with only one normal growing season cannot calculate recovery, because |y(e+1)-y(n)| cannot be calculated 
max.year<-d.select%>%select(site_code, year_trt)%>%group_by(site_code)%>%summarise(n.max=max(year_trt), length.year=length(unique(year_trt)))
years.recovery<-d.select%>% filter(grepl("Extreme", climate))%>%select(site_code, year_trt, climate)%>%distinct()%>%arrange(site_code, year_trt)%>%
  group_by(site_code)%>% mutate(consecutive = lead(year_trt)-year_trt)%>%
  mutate(consecutive1=ifelse(is.na(consecutive), -999, consecutive))%>%
  filter(consecutive1!=1)%>%merge(max.year, by=c("site_code"))%>%filter(year_trt!=n.max)
length(unique(years.recovery$site_code))
d.select.au<-years.recovery%>%select(site_code, year_trt, climate)%>%distinct()%>%
  arrange(site_code, year_trt)%>%
  mutate(year_trt1=year_trt+1)%>%
  mutate(year_trt=year_trt1, year_trt1=NULL)%>%
  merge(y=d.select[,c("site_code", "trt", "year_trt", "block", "live_mass")], by=c("site_code",  "year_trt"))%>%distinct()%>%
  #mutate(events=ifelse(climate=="Dry", "One year after dry",  "One year after wet"))%>%
  mutate(events="One year after")%>%
  bind_rows(d.select.u1)%>%
  merge(d.select_block_avg, by=c("site_code", "block", "trt"))%>%mutate(dev_mass=abs(live_mass-avg_bio))%>%
  mutate(avg_bio=NULL)%>%
  pivot_longer(cols=c("live_mass", "dev_mass"))%>%arrange(name, trt, events)%>%
  mutate(reference.line=ifelse(name=="dev_mass", NA, ifelse((name=="live_mass" & trt=="Control"), d.select_block_sum$avg.mass[d.select_block_sum$trt=="Control"], d.select_block_sum$avg.mass[d.select_block_sum$trt=="NPK"])))

# relevel the variables 
d.select.au$name1<-ifelse(d.select.au$name=="dev_mass", "Deviation in biomass", "Aboveground biomass")
d.select.au$name1<-factor(d.select.au$name1, levels=c("Aboveground biomass", "Deviation in biomass"))
unique(d.select.au$events)
d.select.au$events<-factor(d.select.au$events, levels = c("During",   "One year after"))

## biomass in each climate event average across sites During and one year after
(pp.u.o<-d.select.au%>%filter(!events %in% c("During the last"))%>%mutate(climate1=ifelse(climate=="Extreme dry", "Dry", "Wet"))%>%
    ggplot(aes(events, value, color=trt))+theme_bw(base_size =22)+mt+
    scale_color_manual(values = c("grey", "black"))+
    facet_grid(name1~climate1, scales="free",  switch = "y", space="free_x")+
    geom_violin(mapping = aes(events, value, color=trt), position=pd)+
    geom_point(mapping = aes(events, value, color=trt), position=pd, alpha=0.15)+
    stat_summary(fun=mean, geom="point", size=5, position=pd)+
    stat_summary(fun.data =mean_cl_boot, geom="errorbar", width=0.1, position=pd)+
    geom_hline(aes(yintercept=reference.line, color=trt), linetype="dashed", linewidth=1, alpha=0.6)+
    theme(axis.text.x = element_text(angle=8, vjust = 1.05, hjust=1), strip.placement = "outside", strip.background.y = element_blank())+
    labs(x=NULL, y=NULL, color=NULL))
ggsave(pp.u.o, dpi=600, width=13.3, height=6.64, file="aboveground biomass during and one year after.png")

# check percentage change in biomass difference and deviation across all sites
d.select.au_site<-d.select.au%>%group_by(trt, events, climate, name)%>%summarise(avg=mean(value))%>%
  group_by(name, events, climate)%>%mutate(proposition.change=(lead(avg)-avg)/avg)%>%arrange(name, events, climate)

d.select.au_site[,5:6]<-round(d.select.au_site[,5:6], 2)
names(d.select.au_site)<-c("Trt", "Event time",  "Climate",  "variable", "true values", "proportion change")
write.csv(d.select.au_site, file="summarizing biomass different and deviation over sites.csv")

###########################################################################################
########### formal tests for normal levels and deviation from the normal levels ###########
###########################################################################################
# do a formal test to see whether nutrient addition increased normal-level community aspects  
colnames(d.select_block_avg)
normal.levels<-d.select_block_avg%>% mutate(avg.value1=log(avg_bio))

mod<-lme(avg.value1 ~ trt, random = ~1|site_code/block, data=normal.levels)
# plot(mod)
t<-xtable(summary(mod)$tTable)
rs<-rsquared(mod)
get.sd.sites<-as.numeric(VarCorr(mod)[,"StdDev"][2])
get.sd.block<-as.numeric(VarCorr(mod)[,"StdDev"][4])
eff.normal1<-t%>%mutate(r2m=rs$Marginal, r2c=rs$Conditional, get.sd.sites=get.sd.sites, get.sd.block=get.sd.block,  terms=rownames(.))%>%mutate(across(1:9, \(x) round(x, 2)))%>%dplyr::rename(p='p-value')

# do a formal test to see whether nutrient addition increased deviation during and after a dry and wet growing season relative to normal levels 
colnames(d.select.au)
# add 0.001 because some deviations were 0, this will lead to "Inf" after log transformation 
d.select.au1<-d.select.au%>%select(site_code, year_trt,  climate, trt,  block,  value,   events, climate, name)%>%
  mutate(variable.id=paste(name, events, climate))%>%mutate(value1=log(value))

eff.dif.dev<-c()
for(va in unique(d.select.au1$variable.id)){
  # va<-"composition dev_composition One year after Wet"
  d.select.au1.temp<-d.select.au1%>%filter(variable.id==va)
  mod<-lme(value1 ~ trt, random = ~1|site_code/block/year_trt, data=d.select.au1.temp)
  # plot(mod)
  t<-xtable(summary(mod)$tTable)
  rs<-rsquared(mod)
  get.sd.sites<-as.numeric(VarCorr(mod)[,"StdDev"][2])
  get.sd.block<-as.numeric(VarCorr(mod)[,"StdDev"][4])
  get.sd.year<-as.numeric(VarCorr(mod)[,"StdDev"][6])
  
  t1<-t%>%mutate(r2m=rs$Marginal, r2c=rs$Conditional, get.sd.sites=get.sd.sites, get.sd.block=get.sd.block, get.sd.year=get.sd.year,  terms=rownames(.),  variable.id=va)
  
  eff.dif.dev<-rbind(eff.dif.dev, t1) 
}

eff.dif.dev1<-eff.dif.dev%>%mutate(across(1:10, \(x) round(x, 2)))%>%dplyr::rename(p='p-value')%>%
  merge(d.select.au1%>%select(name, events, climate, variable.id)%>%distinct(), by=c("variable.id"))%>%
  select(name, events, climate, terms, Value, Std.Error, DF, 't-value', p, r2m,  r2c, get.sd.year, get.sd.block, get.sd.sites)
eff.dif.dev.sig<-eff.dif.dev1 %>% filter(p<=0.05 & terms!="(Intercept)")

# add the normal and extreme test together 
colnames(eff.normal1)
colnames(eff.dif.dev1)
eff.normal.dev<-eff.dif.dev1%>%mutate(events1=gsub(" extreme", "", events), events2= paste(events, climate))%>%
  select(name, events2, terms, Value, Std.Error, DF, 't-value', p, r2m,  r2c, get.sd.year, get.sd.block, get.sd.sites)%>%
  rbind(eff.normal1%>%mutate(name="average", events2="During normal", get.sd.year="")%>%select(name, events2, terms, Value, Std.Error, DF, 't-value', p, r2m,  r2c, get.sd.year, get.sd.block, get.sd.sites))
eff.normal.dev %>% filter(p<=0.05 & terms!="(Intercept)")

write.csv(eff.normal.dev, file="full table of treatment effects on normal levels and deviation of aboveground biomass.csv")

###########################################################################################
################## calculate resistance and recovery 
###########################################################################################
pre<-d.select%>%mutate(plot.id=paste(site_code, block, trt, sep="_"))
## calculate mean biomass in normal years 
pre_norm<-pre%>%filter(climate=="Normal")%>% group_by(site_code, trt, block) %>% dplyr::summarise(avg.biomass=mean(live_mass))
pre_r1<-pre%>%merge(pre_norm, by=c("site_code", "trt", "block"))%>%mutate(dev.mass=abs(live_mass - avg.biomass))%>%
  mutate(resistance=avg.biomass/abs(live_mass - avg.biomass))
## get resistance data by deleting normal climate events 
resis1<-pre_r1%>%filter(climate!="Normal")%>%mutate(stability.facets="resistance")%>%dplyr::rename(values=resistance)%>%
  select(site_code, trt, block, plot.id, year, year_trt, climate, spei, stability.facets, values)

# calculate recovery 
recov<-c()
for(pl in unique(pre_r1$plot.id)){
  # pl<-"bnch.us_1_Control"  
  t.data2<-pre_r1%>%filter(plot.id==pl)%>%arrange(year_trt)
  
  t.recovery<-c()
  for(i in 1:nrow(t.data2)){
    avg.biomass<-t.data2$avg.biomass[1]
    t.recovery[i]<-abs(t.data2$live_mass[i] - avg.biomass)/abs(t.data2$live_mass[i+1]-avg.biomass)
  }
  # make sure that recovery is always calculated as first year compared with the second year
  t.data3<-t.data2%>%mutate(consecutive=lead(year_trt)- year_trt, recovery=ifelse(consecutive==1, t.recovery, NA))
  recov<-rbind(recov, t.data3)
}  

recov1<-recov%>%filter(climate!="Normal")%>%mutate(stability.facets="recovery")%>%dplyr::rename(values=recovery)%>%
  select(site_code, trt, block, plot.id, year, year_trt, climate, spei, stability.facets, values)
# add resistance and recovery together 
resis.recov1<-resis1%>%bind_rows(recov1)

############################################################################################
######################calculate invariability during experimental years ####################
############################################################################################
# delete years at some sites that were not included for resistance and recovery 
stb<- pre%>%filter(site_code %in% d.select$site_code)%>%
  group_by(site_code, block, trt, plot.id) %>% 
  summarise_at(c("live_mass"), list(N=~length(.), avg=~mean(.), sd=~sd(.), invariability=~mean(.)/sd(.)))
range(stb$N)
####calculate detrended stability
## a linear term may be enough for most sites, a quadratic term may be ok as well 
sd.d<-c()
for(pl in unique(pre$plot.id)){
  # pl<-"msla_2.us_2_Control"
  data<-pre %>%filter(plot.id==pl)
  if(nrow(data)<3) next
  mod<-lm(live_mass~year_trt, data=data, na.action = na.omit)
  sd.lm<-sd(residuals(mod))
  mod1<-lm(live_mass~year_trt+I(year_trt^2), data=data, na.action = na.omit)
  sd.q<-sd(residuals(mod1))
  temp<-data.frame(plot.id=pl, sd.lm=sd.lm, sd.q=sd.q)
  sd.d<-rbind(sd.d, temp)
}
## add detrended sd to stb data
stb2<-stb%>%merge(y=sd.d, by=c("plot.id"))%>% mutate(invariability.d=avg/sd.lm)%>%
  dplyr::select(site_code, trt, block, invariability, invariability.d)%>%
  pivot_longer(cols = c("invariability", "invariability.d"), names_to ="stability.facets1" , values_to ="values")

############################################################################################
############################ add all stability facets together  ############################
############################################################################################

# calculate mean resistance and recovery based on different cutoffs but keep normal years the same 
data.stability.facets<-c()
for(cut in c(1.28, 1.5, 2)){
  # cut<-1.28
  # when two same climate extremes happen consecutively, recovery only calculate for the later year 
  # the later year can be either a less extreme year or normal year
  resis.recov2<-resis.recov1%>%mutate(spei1=as.numeric(spei), climate1=case_when((spei1>=cut)~"Wet",
                                                                                 (spei1<=- cut)~"Dry",
                                                                                 TRUE~"between.normal.and.extreme"))%>%
    mutate(climate.num=case_when((spei1>=cut)~1,
                                 (spei1<=- cut)~ -1,
                                 TRUE~0))%>%
    arrange(stability.facets, plot.id, year_trt)%>%group_by(stability.facets, plot.id)%>%
    mutate(consecutive = lead(year_trt)-year_trt, climate.type=lead(climate.num)-climate.num)%>%
    mutate(consecutive1=ifelse(is.na(consecutive), -999, consecutive))%>%
    mutate(values1=ifelse((stability.facets=="recovery" & consecutive1==1 & climate.type==0), NA, values))%>%
    filter(climate1%in% c("Wet", "Dry"))
  # check the data 
  check.recovery<-resis.recov2%>%filter(!is.na(values))%>%filter(is.na(values1))%>%ungroup()%>%
    select(site_code, year_trt)%>%distinct()
  # average resistance during dry and wet and recovery from dry and wet for each site 
  s_resis.recov2<-resis.recov2%>%mutate(stability.facets1=paste(stability.facets, climate1, sep="_"))%>%
    filter(!values1%in%c("Inf", "-Inf"))%>%filter(!is.na(values1))%>%
    group_by(site_code, trt, block, stability.facets1)%>%
    summarise_at("values1", list(mean))%>%dplyr::rename(values=values1)
  # add invariability 
  all.stab<-stb2%>%filter(site_code%in%s_resis.recov2$site_code)%>%rbind(s_resis.recov2)             
  all.stab$cutoff<-cut
  data.stability.facets<-rbind(data.stability.facets, all.stab)  
}

#log transform stability facets
data.stability.facets.1<-data.stability.facets%>%mutate(log.value=log(values))%>%
  mutate(stability.facets.log=case_when(stability.facets1=="invariability"~"log.invariability",
                                        stability.facets1=="invariability.d"~"log.invariability.d",
                                        stability.facets1=="resistance_Dry"~"log.resistance.d",
                                        stability.facets1=="resistance_Wet"~"log.resistance.w",
                                        stability.facets1=="recovery_Dry"~"log.recovery.d",
                                        TRUE~"log.recovery.w"))
# check distribution of the data  
(p1<-data.stability.facets.1%>%ggplot(aes(x=log.value, color=trt))+geom_density()+facet_grid(cutoff~stability.facets1))
#############################################################################################
###################### effects of nutrient addition on stability facets #####################
#############################################################################################

output.lme<-data.frame()
## run the models in one go 
for(cut in unique(data.stability.facets.1$cutoff)){ 
  for (i in unique(data.stability.facets.1$stability.facets.log)){
    # i<-"log.recovery.w"; cut<-0.67
    d.e<-data.stability.facets.1%>%filter(cutoff==cut & stability.facets.log==i)%>%filter(!is.na(log.value))%>%filter(!(log.value%in%c("Inf", "-Inf")))
    
    n.sites<-length(unique(d.e$site_code))
    mod<-lme(log.value ~ trt, random= ~1|site_code/block, data=d.e)
    plot(mod)
    t<-xtable(summary(mod)$tTable)
    rs<-rsquared(mod)
    t$r2m<-rs$Marginal
    t$r2c<-rs$Conditional
    t$terms<-rownames(t)
    t$variable<-i
    t$n.sites<-n.sites
    t$cutoff<-cut
    output.lme<-rbind(output.lme, t)
  }
}
output.lme[,1:7]<-lapply(output.lme[,1:7], function(x) round(x, 2))
output.lme.sig<-output.lme%>%dplyr::rename(p1='p-value')%>%filter(terms!="(Intercept)" & p1<=0.1)

colnames(output.lme)
## calculate 95% confidence intervals 
treat2<-output.lme%>%filter(terms!="(Intercept)")%>%filter(variable!="log.invariability")%>%mutate(lower=Value-1.96*Std.Error, upper=Value+1.96*Std.Error)
# adjust variable names 
treat2<-treat2%>%mutate(stability.facets.pub=case_when(variable=="log.invariability.d"~"Temporal invariability",
                                                       variable=="log.resistance.d"~"Resistance during dry growing seasons",
                                                       variable=="log.resistance.w"~"Resistance during wet growing seasons",
                                                       variable=="log.recovery.d"~"Recovery after dry growing seasons",
                                                       TRUE~"Recovery after wet growing seasons"))
treat2$stability.facets.pub<-factor(treat2$stability.facets.pub, levels = c("Temporal invariability",
                                                                            "Resistance during dry growing seasons",
                                                                            "Resistance during wet growing seasons",
                                                                            "Recovery after dry growing seasons",
                                                                            "Recovery after wet growing seasons"))

## show results based on cutoff of 1.28 and add data points from individual sites
stability.individual.sites<-data.stability.facets.1%>%filter(stability.facets1!="invariability" & cutoff==1.28)%>%
  select(site_code, block, trt, log.value, stability.facets.log)%>%
  pivot_wider(names_from=trt, values_from=log.value)%>%mutate(dif=NPK-Control)%>%
  # unique(stability.individual.sites$stability.facets.log)
  mutate(stability.facets.pub=case_when(stability.facets.log=="log.invariability.d"~"Temporal invariability",
                                        stability.facets.log=="log.resistance.d"~"Resistance during dry growing seasons",
                                        stability.facets.log=="log.resistance.w"~"Resistance during wet growing seasons",
                                        stability.facets.log=="log.recovery.d"~"Recovery after dry growing seasons",
                                        TRUE~"Recovery after wet growing seasons"))%>% filter(!is.na(dif))

stability.individual.sites$stability.facets.pub<-factor(stability.individual.sites$stability.facets.pub, levels = c("Temporal invariability",
                                                                                                                    "Resistance during dry growing seasons",
                                                                                                                    "Resistance during wet growing seasons",
                                                                                                                    "Recovery after dry growing seasons",
                                                                                                                    "Recovery after wet growing seasons"))

(pp.stability.facets<-treat2%>%filter(cutoff==1.28)%>%ggplot(aes(stability.facets.pub, Value))+theme_bw(base_size =22)+mt+
    #facet_wrap(~cutoff, nrow=5)+
    geom_violin(stability.individual.sites, mapping = aes(stability.facets.pub, dif), color="light grey")+
    geom_point(stability.individual.sites, mapping = aes(stability.facets.pub, dif), alpha=0.15)+
    geom_point(position=pd, size=5, pch=19)+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, position=pd) +
    geom_hline(yintercept = 0, linetype ="dotted") +
    theme(axis.text.x = element_text(angle=25, hjust=0.95))+
    labs(x=NULL, y="Effect size of nutrient addition"))
ggsave(pp.stability.facets, dpi=600, width=13.3, height=6.64, file="treatment effects on stability facets for all sites cutoff of 1.28.png")

###########################################################################################
##########get species richness and add it together with functional traits
###########################################################################################
colnames(d7[,c(1:3, ncol(d7))])
rich<-d7%>%select(site_code, trt, block, year_trt, richness)%>%distinct()%>%filter(site_code %in% d.select$site_code)
range(rich$richness)
## explore the data distribution 
ggplot(rich)+geom_histogram(aes(x=richness))+theme_bw(base_size=18)+facet_wrap(~trt)## a little skewed

# load trait data 
fun2<-read.csv("fd and cwm using traits measured at NutNet sites.csv")
colnames(fun2)
# load trait data for shared species 
fun.share2<-read.csv("fd and cwm using traits measured at NutNet sites for shared species.csv")
colnames(fun.share2)
# merge all plant diversity facets and delete years that were not used for resistance and recovery
nutnet.traits1<-fun2%>%mutate(species="All")%>%select("species", "trait.name" , "site_code", "block" ,  "trt", "year_trt",  "f.dispersion" ,   "cwm1" )%>%
  bind_rows(fun.share2%>%mutate(species="Shared")%>%select("species", "trait.name" , "site_code", "block" ,  "trt", "year_trt",  "f.dispersion" ,   "cwm1" ))
colnames(nutnet.traits1)

global.traits<-read.csv("fd and cwm based on compiled traits.csv")
colnames(global.traits)
unique(global.traits$data.type1)
global.traits1<-global.traits%>%select(data.type1, site_code, block, trt, year_trt, f.dispersion, cwm1, trait.name, total_cover, total_cover_traits)%>%
  merge(nutnet.traits1%>%select(site_code, block, trt, year_trt)%>%distinct(), by=c("site_code", "block", "trt", "year_trt"))%>%
  mutate(species="All")%>%mutate(source=ifelse(data.type1=="species.and.genus", "Globaldatabases_genus", "Globaldatabases_species"))%>%
  select(species, trait.name, site_code, year_trt,  block, trt, f.dispersion,  cwm1,  source)

# add all traits together; 
colnames(global.traits1)
colnames(nutnet.traits1)
driver1<-nutnet.traits1%>% mutate(source="NutNet")%>%
  select(species, trait.name, site_code, year_trt,  block, trt, f.dispersion,  cwm1,  source)%>%bind_rows(global.traits1)%>%
  merge(rich[,c("site_code","block","trt", "year_trt", "richness")], by=c("site_code", "year_trt", "block", "trt")) %>%
  mutate(trait.name1=ifelse((trait.name %in% c("slow_fast")), "LES", trait.name))%>%
  mutate(trait.name2=ifelse((trait.name1 %in% c("slow_fast_without_K")), "LES without K", trait.name1))%>%
  mutate(trait.name3=ifelse((trait.name2 %in% c("slow_fast_with_LDMC")), "LES with LDMC", trait.name2))%>%
  mutate(predict.id=paste(source, species, trait.name3, sep = "_"))%>%filter(species %in% c("All"))

#############################################################################################
# SEM for the NPK effects on stability and diversity facets using sites have both dry and wet growing seasons
#############################################################################################
## fix the convergence error for lme 
cc<-lmeControl(opt="optim")
# cc<-lmeControl(opt = "optim", optMethod = "L-BFGS-B")
# when cutoff is 1.5, only 3 sites are available for SEM, which is too few. 
data.stability.facets.2<-data.stability.facets.1%>%filter(stability.facets1!="invariability" & cutoff==1.28)%>%
  select(site_code, block, trt, stability.facets.log, log.value)%>%
  pivot_wider(names_from="stability.facets.log", values_from="log.value")

r.sem.t<-c();r.sem.v<-c(); partition.eff.sem<-c(); r.sem.fit<-c()

for(id in unique(driver1$predict.id)){
  # id<-"NutNet_All_Leaf C"
  ## combine traits and stability facets 
  fun.temp<-driver1%>%filter(predict.id== id)%>%
    group_by(predict.id, site_code, trt, block)%>%dplyr::summarise_at(c("cwm1", "f.dispersion", "richness"), list(mean))%>%distinct()
  # ggpairs(fun.temp[, c("trt", "cwm1", "f.dispersion", "richness")], aes(color=trt))
  
  sss.dw1<-data.stability.facets.2%>%
    merge(y=fun.temp, by=c("site_code","trt","block"))%>%mutate(trt1=ifelse(trt=="Control", 0, 1))
  N.site.rr<-length(unique(sss.dw1$site_code))## 
  
  ## include CWM, FD, and richness as predictors
  mod.tt.dw<-psem(
    lme(log.invariability.d ~log.resistance.d +log.resistance.w + log.recovery.d +log.recovery.w +cwm1  +f.dispersion  +richness +trt1, random = ~1|site_code/block, control=cc, data=sss.dw1),
    lme(log.resistance.d ~cwm1  +f.dispersion  +richness +trt1, random = ~1|site_code/block, control=cc, data=sss.dw1),
    lme(log.resistance.w ~cwm1  +f.dispersion  +richness +trt1, random = ~1|site_code/block, control=cc, data=sss.dw1),
    lme(log.recovery.d ~ cwm1  +f.dispersion  +richness +trt1, random = ~1|site_code/block, control=cc, data=sss.dw1),
    lme(log.recovery.w ~ cwm1  +f.dispersion  +richness +trt1, random = ~1|site_code/block, control=cc, data=sss.dw1),
    lme(richness ~ trt1, random = ~1|site_code/block, control=cc, data=sss.dw1),
    lme(cwm1 ~ trt1, random = ~1|site_code/block, control=cc, data=sss.dw1),
    lme(f.dispersion ~ trt1, random = ~1|site_code/block, control=cc, data=sss.dw1),
    f.dispersion %~~% cwm1,
    log.recovery.w %~~% log.resistance.w,
    log.recovery.d %~~%  log.resistance.d,
    log.recovery.w %~~% log.resistance.d) # this correlated error is needed when using cutoff of 1.28
  lapply(mod.tt.dw[1:5], vif)
  
  (t.dw<-summary(mod.tt.dw))
  t.coef.dw<-xtable(t.dw$coefficients)
  t.coef.dw.sig<-t.coef.dw%>%data.frame()%>%filter(P.Value < 0.1)
  colnames(t.coef.dw)[9]<-"sig"
  t.r2.dw<-xtable(t.dw$R2)
  # t.model.fit1<-xtable(t.dw$ChiSq)
  t.model.fit2<-xtable(t.dw$Cstat)
  
  ## sort out model output 
  r.sem.t<-r.sem.t%>%bind_rows(t.coef.dw%>%data.frame()%>%mutate(predict.id=id, N.site=N.site.rr))
  r.sem.v<-r.sem.v%>%bind_rows(t.r2.dw%>%data.frame()%>%mutate(predict.id=id, N.site=N.site.rr))
  r.sem.fit<-r.sem.fit%>%bind_rows(t.model.fit2%>%data.frame()%>%mutate(predict.id=id))
  
  ############ summarize the direct and indirect effects of NPK on stability facets #############
  # [] numbers inside are positions of each variable in the SEM, from left to right, from up to down
  #############################resistance.d#################################################
                      # indirect through cwm
  ind.eff.resis.d<- (t.coef.dw$Std.Estimate[26]*t.coef.dw$Std.Estimate[9]+
                       # indirect through fd
                       t.coef.dw$Std.Estimate[27]*t.coef.dw$Std.Estimate[10]+
                       # indirect through richness
                       t.coef.dw$Std.Estimate[25]*t.coef.dw$Std.Estimate[11]) 
  # total and direct effects of trt on resistance.d:
  dir.eff.resis.d<-t.coef.dw$Std.Estimate[12]
  total.eff.resis.d<-dir.eff.resis.d+ind.eff.resis.d
  
  #############################resistance.w#################################################
                      # indirect through cwm
  ind.eff.resis.w<- (t.coef.dw$Std.Estimate[26]*t.coef.dw$Std.Estimate[13]+
                       # indirect through fd
                       t.coef.dw$Std.Estimate[27]*t.coef.dw$Std.Estimate[14]+
                       # indirect through richness
                       t.coef.dw$Std.Estimate[25]*t.coef.dw$Std.Estimate[15]) 
  # total and direct effects of trt on resistance.w:
  dir.eff.resis.w<-t.coef.dw$Std.Estimate[16]
  total.eff.resis.w<-dir.eff.resis.w+ind.eff.resis.w
  
  #############################recovery.d#################################################
                      # indirect through cwm
  ind.eff.recov.d<- (t.coef.dw$Std.Estimate[26]*t.coef.dw$Std.Estimate[17]+
                       # indirect through fd
                       t.coef.dw$Std.Estimate[27]*t.coef.dw$Std.Estimate[18]+
                       # indirect through richness
                       t.coef.dw$Std.Estimate[25]*t.coef.dw$Std.Estimate[19]) 
  # total and direct effects of trt on recovery.d:
  dir.eff.recov.d<-t.coef.dw$Std.Estimate[20]
  total.eff.recov.d<-dir.eff.recov.d+ind.eff.recov.d
  
  #############################recovery.w#################################################
                      # indirect through cwm
  ind.eff.recov.w<- (t.coef.dw$Std.Estimate[26]*t.coef.dw$Std.Estimate[21]+
                       # indirect through fd
                       t.coef.dw$Std.Estimate[27]*t.coef.dw$Std.Estimate[22]+
                       # indirect through richness
                       t.coef.dw$Std.Estimate[25]*t.coef.dw$Std.Estimate[23]) 
  # total and direct effects of trt on recovery.w:
  dir.eff.recov.w<-t.coef.dw$Std.Estimate[24]
  total.eff.recov.w<-dir.eff.recov.w+ind.eff.recov.w
  
  #############################invariability#################################################
               # indirect through cwm
  ind.eff.inv<-t.coef.dw$Std.Estimate[26]*t.coef.dw$Std.Estimate[5]+
    # indirect through fd
    t.coef.dw$Std.Estimate[27]*t.coef.dw$Std.Estimate[6]+
    # indirect through richness
    t.coef.dw$Std.Estimate[25]*t.coef.dw$Std.Estimate[7]
  # indrect through cwm, fd, richness, and then through stability facets
  ind.eff.inv.stab<- ind.eff.resis.d*t.coef.dw$Std.Estimate[1]+ind.eff.resis.w*t.coef.dw$Std.Estimate[2]+
    ind.eff.recov.d*t.coef.dw$Std.Estimate[3]+ind.eff.recov.w*t.coef.dw$Std.Estimate[4]
  # direct effects of trt on invariability:
  dir.eff.inv<-t.coef.dw$Std.Estimate[8]
  # direct through other stability facets
  dir.eff.inv.stab <- dir.eff.resis.d*t.coef.dw$Std.Estimate[1]+dir.eff.resis.w*t.coef.dw$Std.Estimate[2]+
    dir.eff.recov.d*t.coef.dw$Std.Estimate[3]+dir.eff.recov.w*t.coef.dw$Std.Estimate[4]
  # total effects
  total.eff.inv<-dir.eff.inv+dir.eff.inv.stab + ind.eff.inv+ind.eff.inv.stab
  #############################diversity facets################################################
  total.eff.rich<-t.coef.dw$Std.Estimate[25]
  total.eff.cwm<-t.coef.dw$Std.Estimate[26]
  total.eff.fd<-t.coef.dw$Std.Estimate[27]
  eff.sem<-data.frame(c(total.eff.cwm, 0, 0, 0,total.eff.cwm), c(total.eff.fd, 0, 0,0, total.eff.fd), c(total.eff.rich, 0, 0, 0,total.eff.rich),
                      c(dir.eff.resis.d, 0, ind.eff.resis.d, 0,total.eff.resis.d), c(dir.eff.resis.w, 0, ind.eff.resis.w, 0,total.eff.resis.w), 
                      c(dir.eff.recov.d, 0, ind.eff.recov.d, 0,total.eff.recov.d),  c(dir.eff.recov.w, 0, ind.eff.recov.w, 0, total.eff.recov.w), 
                      c(dir.eff.inv, dir.eff.inv.stab, ind.eff.inv, ind.eff.inv.stab, total.eff.inv))
  colnames(eff.sem)<-c("CWM", "FD", "Richness", "log.resistance.d", "log.resistance.w",  "log.recovery.d", "log.recovery.w", "log.invariability.d")
  eff.sem$effects.paths<-c("Direct effects", "Direct effects through resistance and recovery", "Indirect effects through diversity facets",  "Indirect effects through stability facets", "Total effects")
  
  partition.eff.sem<-partition.eff.sem%>%bind_rows(eff.sem%>%data.frame()%>%mutate(predict.id=id))
}

# sort the results
colnames(r.sem.t)
r.sem.t1<-r.sem.t
r.sem.t1[,3:8]<-lapply(r.sem.t1[,3:8], as.numeric)
r.sem.t1[,3:8]<-lapply(r.sem.t1[,3:8], function(x) round(x, 2))
r.sem.t2<-driver1%>%select(source, species, trait.name3, predict.id)%>%distinct()%>%merge(r.sem.t1, by=c("predict.id"))%>%
  dplyr::select(Response, Predictor,  Std.Estimate,  DF,  P.Value,  N.site, source, species, trait.name3, predict.id)
# significant effects 
r.sem.t2.sig<-r.sem.t2%>%filter(P.Value<=0.1 & species=="All" )
# save the significant test statistics
write.csv(r.sem.t2.sig, file="significant effects in SEM.csv")

colnames(r.sem.t2); colnames(r.sem.v); colnames(r.sem.fit)
r.sem.t3<-r.sem.t2%>%
  left_join(r.sem.v[,c("predict.id", "Response", "Marginal",  "Conditional")], by=c( "predict.id", "Response"))%>%
  merge(r.sem.fit%>%dplyr::rename(p=P.Value), by=c( "predict.id"))%>%
  mutate(Response1=gsub("cwm1", "CWM", Response))%>%mutate(Response2=gsub("f.dispersion", "FD", Response1))%>%
  mutate(Predictor1=gsub("cwm1", "CWM", Predictor))%>%mutate(Predictor2=gsub("f.dispersion", "FD", Predictor1))%>%mutate(Predictor3=gsub("trt1", "trt", Predictor2))%>%
  select(source, species, trait.name3, Response2, Predictor3,  Std.Estimate,  DF,  P.Value,  Marginal, Conditional, Fisher.C, df, p)

colnames(r.sem.t3)<-c("source", "Species", "Trait name", "Response", "Predictor",  "Std.Estimate",  "DF",  "P",  "R2m",    "R2c", "Fisher.C", "df", "p")
write_xlsx(r.sem.t3, path="variance explained in sem.xlsx")

colnames(output.lme)
output.lme.2<-output.lme%>%
  select("cutoff", "variable", "terms", "Value",     "Std.Error", "DF",        "t-value",   "p-value",   "r2m",       "r2c", "n.sites")
colnames(output.lme.2)<-c("Cutoff", "Variable", "Terms", "Value",  "Std.Error", "DF",  "t-value",   "p", "R2m",    "R2c", "N.sites")
write_xlsx(output.lme.2, path="effects of nutrients on stability facets using bivariate analyses.xlsx")

partition.eff.sem[partition.eff.sem==0]<-NA
colnames(partition.eff.sem)
partition.eff.sem.1<-driver1%>%select(source, species, trait.name3, predict.id)%>%distinct()%>%merge(partition.eff.sem, by=c("predict.id"))%>%
  select("source", "trait.name3", "effects.paths", "CWM",  "FD", "Richness",
         "log.resistance.d",  "log.resistance.w", 
         "log.recovery.d",  "log.recovery.w",  "log.invariability.d")

partition.eff.sem.1[,4:11]<-lapply(partition.eff.sem.1[,4:11], function(x) round(x, 2))
write_xlsx(partition.eff.sem.1, path="partition direct and indirect effects in SEM.xlsx")

# the end