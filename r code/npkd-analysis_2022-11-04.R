###### NPKD-net data cleaning and initial analysis

##### data management
library(tidyverse)

## analysis
library(nlme)
library(emmeans)
library(vegan)

#### graphics
library(cowplot)

theme_figs <- theme_classic() +
  theme(axis.title = element_text(size = 24), axis.text = element_text(size = 24),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = 'right', legend.justification = c(1, 0.4),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_text(size=18,face='bold'),legend.text = element_text(size = 18),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        strip.text =element_text(size=18,face='bold'),
        axis.line = element_line(colour = 'grey60', size = 0.35), axis.ticks = element_line(colour = 'grey60', size = 0.35)) 
theme_set(theme_figs)


mass <- read_csv('./data/npkd-full-biomass_2022-11-04.csv', na=c("NULL", "NA"))
cover <- read_csv('./data/npkd-full-cover_2022-11-04.csv', na=c("NULL", "NA"))
comb <- read_csv('./data/npkd-comb_clim_soil_2022-11-04.csv', na=c("NULL", "NA"))


comb_npk_max <- comb %>% filter(!site_code %in% c('paike.ar','yarradrt.au','freiburg.de','sonora.us')) %>% 
  filter(!substr(site_code,nchar(site_code)-1,nchar(site_code)) == 'cn') %>%
group_by(site_code) %>% mutate(max_year = max(year)) %>%
 ungroup() %>%
 mutate(max_year = if_else(site_code == 'chilcasdrt.ar', 2019, max_year)) %>%
 filter(year == max_year)


unique(comb_npk_max$site_code)
#8

comb_npk_max %>% select(live_mass) %>% print(n = 500)
comb_npk_max %>% filter(is.na(live_mass)) %>% select(site_code, year, year_trt, plot, total_mass, live_mass, dead_mass) %>% print(n=500)

mass %>% filter(site_code == 'chilcasdrt.ar' & year == 2020) %>% print(n=200)
# chilcas had covid restrictions in 2020, defining 2019 as max year above

comb_npk_max[is.na(comb_npk_max$live_mass),]$live_mass <- comb_npk_max[is.na(comb_npk_max$live_mass),]$total_mass

comb_npk_max[is.na(comb_npk_max$live_mass),] #10

comb_npk_max <- filter(comb_npk_max, !is.na(live_mass))

#### total mass analysis ####

mod.mass <- lme(log(live_mass) ~ Drought*N, random = ~1|site_code/plot, data = comb_npk_max)
summary(mod.mass)

em.mass <- data.frame(emmeans(mod.mass, ~Drought*N))

pd <- position_dodge(0.3)

gg_mass <- ggplot(em.mass, aes(x = as.factor(N), y = emmean, col = as.factor(Drought))) + geom_point(position = pd, size = 3) + 
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL),width=0.3,position=pd, size = 2) +
  #geom_jitter(comb_npk_max)
  scale_color_manual(values= c('gray25','firebrick3'),guide = guide_legend(title='Drought Treatment'),
                     labels= c('0' = 'Control','1'='Drought')) +
  ylab('Biomass') + xlab('Soil resources') + scale_x_discrete(labels = c('Ambient', '+NPK')) +
    theme_bw()

gg_mass


#### biomass across precip gradient ####

mass_lrr <- comb_npk_max %>% 
  #filter(trt == 'Control') %>% 
  group_by(site_code,block,year) %>% 
  mutate(ctrl_mass = mean(live_mass[trt == 'Control'])) %>% mutate(LRR_mass = log(live_mass/ctrl_mass)) 
# %>% 
#   filter((trt != 'Control'))

mod_mass_ppt <- lme(LRR_mass ~ MAP_v2*Drought*N, random = ~1|site_code, data = mass_lrr)
summary(mod_mass_ppt)

em_mass_ppt <- data.frame(summary(emmeans(mod_mass_ppt, pairwise ~ N * Drought | MAP_v2,
                                          at = list(MAP_v2 = seq(from = min(mass_lrr$MAP_v2), to =  max(mass_lrr$MAP_v2), by =max(mass_lrr$MAP_v2)-min(mass_lrr$MAP_v2))),
                                          cov.reduce = TRUE,
                                          type = "response")$emmeans))

em_mass_ppt

em_mass_ppt <- em_mass_ppt %>%
  mutate(trt = if_else(N == 0 & Drought == 0, 'Control',
                       if_else(N == 0 & Drought == 1, 'Drought',
                               if_else(N == 1 & Drought == 0, 'NPK',
                                       if_else(N == 1 & Drought == 1, 'NPK+Drought',NULL))))) %>% 
  filter(trt != 'Control')

col_pal <- c("#D53E4F","#6BAED6","#5E4FA2")

gg_mass_ppt <- ggplot(em_mass_ppt, aes(x = MAP_v2, y = emmean, col = trt)) +
  #geom_point(data = mass_lrr %>% filter(trt != 'Control'), aes(x=MAP_v2, y = LRR_mass, col = trt), alpha = 0.7) +
  geom_ribbon(aes(ymin = lower.CL,ymax= upper.CL,fill=trt),alpha=0.25,color=NA) +
  geom_line(aes(group = trt), size = 2, color= 'black') +
  geom_line(aes(group = trt), size = 1.8) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('LRR biomass') + xlab('Mean Annual Precipitation') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) + xlim(300,1400)

gg_mass_ppt

#### richness analysis ####

mod.rich <- lme(rich ~ Drought*N, random = ~1|site_code, data = comb_npk_max)
summary(mod.rich)

em.rich <- data.frame(emmeans(mod.rich, ~Drought*N))

pd <- position_dodge(0.3)

gg_rich <- ggplot(em.rich, aes(x = as.factor(N), y = emmean, col = as.factor(Drought))) + geom_point(position = pd, size = 3) + 
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL),width=0.3,position=pd, size = 2) +
  #geom_jitter(comb_npk_max)
scale_color_manual(values= c('gray25','firebrick3'),guide = guide_legend(title='Drought Treatment'),
                   labels= c('0' = 'Control','1'='Drought')) +
  ylab('Species Richness') + xlab('Soil resources') + scale_x_discrete(labels = c('Ambient', '+NPK')) +
  theme_bw()

gg_rich

### 

#### richness across precip gradient ####

rich_lrr <- comb_npk_max %>% 
  #filter(trt == 'Control') %>% 
  group_by(site_code,year) %>% 
  mutate(ctrl_rich = mean(rich[trt == 'Control'])) %>% mutate(LRR_rich = log(rich/ctrl_rich)) 
# %>% 
#   filter((trt != 'Control'))

mod_rich_ppt <- lme(LRR_rich ~ MAP_v2*Drought*N, random = ~1|site_code, data = rich_lrr)
summary(mod_rich_ppt)

em_rich_ppt <- data.frame(summary(emmeans(mod_rich_ppt, pairwise ~ N * Drought | MAP_v2,
                                          at = list(MAP_v2 = seq(from = min(rich_lrr$MAP_v2), to =  max(rich_lrr$MAP_v2), by =max(rich_lrr$MAP_v2)-min(rich_lrr$MAP_v2))),
                                          cov.reduce = TRUE,
                                          type = "response")$emmeans))

em_rich_ppt

em_rich_ppt <- em_rich_ppt %>%
  mutate(trt = if_else(N == 0 & Drought == 0, 'Control',
                       if_else(N == 0 & Drought == 1, 'Drought',
                               if_else(N == 1 & Drought == 0, 'NPK',
                                       if_else(N == 1 & Drought == 1, 'NPK+Drought',NULL))))) %>% 
  filter(trt != 'Control')

col_pal <- c("#D53E4F","#6BAED6","#5E4FA2")

gg_rich_ppt <- ggplot(em_rich_ppt, aes(x = MAP_v2, y = emmean, col = trt)) +
  #geom_jitter(data = dominant_pop, aes(x=year_trt, y = perc_rank),  width=0.11, alpha = 0.7, col = 'grey90') +
  geom_ribbon(aes(ymin = lower.CL,ymax= upper.CL,fill=trt),alpha=0.25,color=NA) +
  geom_line(aes(group = trt), size = 2, color= 'black') +
  geom_line(aes(group = trt), size = 1.8) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('LRR richness') + xlab('Mean Annual Precipitation') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal)

gg_rich_ppt


### let's look at functional group biomass and relative cover breakdown
#### func_group ####


unique_cat <- function(x) {length(unique(x))}

tapply(mass$category,mass$site_code,unique_cat)
#chinese sites aand netherlands only have one category (presumed 'total')

unique(mass$category)
mass[is.na(mass$category) & mass$site_code == 'bayrdrt.de',]

mass[mass$site_code == 'bayrdrt.de' & mass$plot == 1,]

gram_mass <- mass %>% filter(!(category %in% c('TOTAL','LITTER',"Pine needles","STANDING DEAD") | is.na(category))) %>% 
  group_by(site_code, plot,year) %>% 
  mutate(total_mass = sum(mass)) %>% 
  mutate(rel_mass = mass/total_mass) %>% 
  filter(category == 'GRAMINOID') %>% 
  group_by(site_code) %>% mutate(max_year = max(year)) %>%
  ungroup() %>% 
  mutate(max_year = if_else(site_code == 'chilcasdrt.ar', 2019, max_year)) %>% 
  filter(year == max_year & !site_code %in% c('paike.ar','sonora.us','haibei.cn','freiburg.de')) %>% 
  left_join(.,comb_npk_max %>% select(site_code,plot,Drought,N),by = c('site_code','plot'))

unique(gram_mass$site_code)

mod_gram_mass <- lme(rel_mass ~ Drought*N, random = ~1|site_code, data = gram_mass)
summary(mod_gram_mass)

em_gram_mass <- data.frame(emmeans(mod_gram_mass, ~Drought*N))

pd <- position_dodge(0.3)

gg_gram_mass <- ggplot(em_gram_mass, aes(x = as.factor(N), y = emmean, col = as.factor(Drought))) + geom_point(position = pd, size = 3) + 
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL),width=0.3,position=pd, size = 2) +
  #geom_jitter(comb_npk_max)
scale_color_manual(values= c('gray25','firebrick3'),guide = guide_legend(title='Drought Treatment'),
                   labels= c('0' = 'Control','1'='Drought')) +
  ylab('Relative Graminoid Biomass') + xlab('Soil resources') + scale_x_discrete(labels = c('Ambient', '+NPK')) +
  theme_bw()

gg_gram_mass


#### Relative cover graminoid ####

gram_cover<- cover %>% filter(!(functional_group %in% c('NON-LIVE') | is.na(functional_group))) %>% 
  group_by(site_code) %>% mutate(max_year = max(year)) %>%
  ungroup() %>% 
  mutate(max_year = if_else(site_code == 'chilcasdrt.ar', 2019, max_year)) %>% 
  filter(year == max_year & !site_code %in% c('paike.ar','sonora.us','haibei.cn','freiburg.de')) %>% 
  group_by(site_code, plot,year) %>% 
  mutate(total_cover = sum(max_cover)) %>% 
  mutate(rel_cover = max_cover/total_cover) %>% 
  filter(functional_group %in% c('GRAMINOID','GRASS')) %>% 
  summarize(rel_gram_cover = sum(rel_cover))%>% 
  left_join(.,comb_npk_max %>% select(site_code,plot,Drought,N),by = c('site_code','plot')) %>% 
  filter(!is.na(N))

mod_gram_cover <- lme(rel_gram_cover ~ Drought*N, random = ~1|site_code, data = gram_cover)
summary(mod_gram_cover)

em_gram_cover <- data.frame(emmeans(mod_gram_cover, ~Drought*N))

pd <- position_dodge(0.3)

gg_gram_cover <- ggplot(em_gram_cover, aes(x = as.factor(N), y = emmean, col = as.factor(Drought))) + geom_point(position = pd, size = 3) + 
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL),width=0.3,position=pd, size = 2) +
  #geom_jitter(comb_npk_max)
  scale_color_manual(values= c('gray25','firebrick3'),guide = guide_legend(title='Drought Treatment'),
                     labels= c('0' = 'Control','1'='Drought')) +
  ylab('Relative Graminoid Cover') + xlab('Soil resources') + scale_x_discrete(labels = c('Ambient', '+NPK')) +
  theme_bw()

gg_gram_cover


#### NMDS of composition ####

gg_list <- NULL

for(i in unique(gram_cover$site_code)){
  df_wide <- cover %>% filter(site_code == i & live == 1) %>% 
    #filter(year_trt == max(year_trt)) %>% 
    pivot_wider(id_cols = c('site_code','plot','year'), names_from = 'Taxon',values_from = 'max_cover', values_fill = 0)
  
  nmds <- metaMDS(df_wide[,4:ncol(df_wide)])
  nmds_scores = as.data.frame(scores(nmds)$sites)
  nms_site <- cbind(df_wide[,1:3],nmds_scores) %>% 
    left_join(.,comb %>% select(site_code,plot,year,year_trt,Drought,N),by=c('site_code','plot','year')) %>% 
    #filter(year != min(year))
    mutate(max_year = max(year)) %>% filter(year == max_year)
  # 
  gg_nms <- ggplot(nms_site,aes(x = NMDS1, y = NMDS2, col = as.factor(Drought), shape = as.factor(N))) + geom_point(size = 3) +
    scale_color_manual(values= c('gray25','firebrick3'),guide = guide_legend(title='Drought Treatment'),
                       labels= c('0' = 'Control','1'='Drought'))+
    scale_shape_manual(values = c(16,17),guide = guide_legend(title='NPK Treatment'),
                       labels= c( 'Control', 'NPK')) + theme_bw()
  
  gg_nms_empty <- ggplot(nms_site,aes(x = NMDS1, y = NMDS2, col = as.factor(Drought), shape = as.factor(N))) + geom_point(size = 3) +
    scale_color_manual(values= c('gray25','firebrick3'),guide = 'none',
                       labels= c('0' = 'Control','1'='Drought'))+
    scale_shape_manual(values = c(16,17),guide = 'none',
                       labels= c( 'Control', 'NPK')) + theme_bw() + ggtitle(label = i)
  
  
    
  gg_nms  
  gg_nms_empty
  
  assign(paste0('gg_',i),gg_nms_empty)
  }

plot_grid(gg_baddrt.de,gg_bayrdrt.de,gg_cdpt_drt.us,gg_cedarsav.us,
                       gg_chilcasdrt.ar,gg_marcdrt.ar,gg_rhijn.nl,gg_sand.us,
                      nrow = 2)
          