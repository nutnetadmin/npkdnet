###### NPKD-net data cleaning and initial analysis

##### data management
library(tidyverse)

## analysis
library(nlme)
library(emmeans)
library(vegan)

#### graphics
library(cowplot)
library(grid)
library(gridExtra)

theme_figs <- theme_classic() +
  theme(axis.title = element_text(size = 24), axis.text = element_text(size = 24),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = 'right', legend.justification = c(1, 0.4),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_text(size=24,face='bold'),legend.text = element_text(size = 24),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        strip.text =element_text(size=18,face='bold'),
        axis.line = element_line(colour = 'grey60', size = 0.35), axis.ticks = element_line(colour = 'grey60', size = 0.35)) 
theme_set(theme_figs)

summary.tablefunc <- function(mod) {  
  dat <- data.frame(summary(mod)$tTable) %>%
    tibble::rownames_to_column(var = 'Effect') %>%
    rename_with(stringr::str_replace, 
                pattern = "-", replacement = ".", 
                matches("Length")) %>% 
    dplyr::mutate(Estimate = signif(Value, digits = 3),
                  Std.Error = signif(Std.Error, digits = 2),
                  t.value = signif(t.value, digits = 2),
                  p.value = signif(p.value, digits = 2)) %>%
    dplyr::mutate(p.value = ifelse(p.value <= 0.001, '< 0.001', p.value)) %>% 
    dplyr::select(-Value) %>% 
    relocate(Estimate,.before = Std.Error)
  return(dat)
}


col_pal <- c("#D53E4F","#6BAED6","#5E4FA2")

setwd('/Users/wilf0020/Library/CloudStorage/Dropbox/NPKDNet')

mass <- read_csv('./full-biomass_2023-05-07.csv', na=c("NULL", "NA"))
cover <- read_csv('./full-cover_2023-04-27.csv', na=c("NULL", "NA"))
comb <- read_csv('./comb-by-plot_2023-05-07.csv', na=c("NULL", "NA"))
ppt <- read_csv('/Users/wilf0020/Library/Mobile Documents/com~apple~CloudDocs/Documents/NPKD Net/npkdnet/data/MSWEP-precipitation-by-treatment-date_2023-05-11.csv')

unique(mass$max_date[mass$site_code=='chilcasdrt.ar'])
mass$max_date <- gsub('2016-11-03','2016-03-11',mass$max_date)

unique(comb$trt)

#reduce treatments to common categories
comb <- comb  %>% 
  mutate(trt_use = 
                 if_else(trt %in% c('NPK','P','N',"NPK+Control"), 'Fertilization',
                    if_else(trt %in% c('NPK+Drought','P+Drought','N+Drought'), 'Fertilization+Drought', trt)))

unique(comb$trt_use)
                     
ppt %>% distinct(site_code,n_nutrient_days, n_drought_days) %>% print(n=150)
ppt <- ppt %>% filter(!(site_code == 'chilcasdrt.ar' & n_drought_days < 300)) %>% mutate(map_cv = map_sd/map)

unique(ppt$n_drought_days[ppt$site_code=='chilcasdrt.ar'])

comb_yr1 <- left_join(comb,ppt, by = c('site_code','year')) %>% filter(n_drought_days > 20) %>%
  group_by(site_code,plot,subplot) %>% 
  filter(n_nutrient_days == min(n_nutrient_days[n_nutrient_days > 20])) %>% mutate(id = paste0(site_code,plot,subplot,year))

write_csv(data.frame(unique(comb_yr1$site_code)),
          '/Users/wilf0020/Library/Mobile Documents/com~apple~CloudDocs/Documents/NPKD Net/npkdnet/data/site_list.csv')

### extract live biomass from full_biomass

### calculate long-term plot average mass for each treatment
massM <- mass %>% mutate(live = if_else(is.na(live),1,live)) %>% 
  #filter(n_nutrient_days > -365) %>% 
  group_by(site_code,live,plot,subplot,trt,year) %>% 
  summarize(sum_mass = sum(mass)) %>% 
  ungroup() %>% group_by(site_code,trt,live) %>% 
  summarize(massM = mean(sum_mass))
  
massG <- mass %>% complete(nesting(site_code,plot,block, subplot, year,trt, n_drought_yrs,n_nutrient_yrs,max_date), category, fill = list(mass =0)) %>%
  filter(!(mass == 0 & category != 'GRAMINOID')) %>%  
  mutate(live = if_else(category == 'GRAMINOID',1,live)) %>% group_by(site_code,plot,year,live) %>% 
  mutate(total_mass = sum(mass), prop_mass = mass/total_mass) %>% 
  filter(category == 'GRAMINOID', trt == 'Control') %>% 
  group_by(site_code,year) %>% 
  summarise(avg_gram_prop = mean(prop_mass))

mass_yr1 <- mass %>% mutate(id = paste0(site_code,plot,subplot,year)) %>% filter(id %in% comb_yr1$id)  %>% 
  filter(!trt %in% c("Control Infrastructure","NPK+Water_addition","Water_addition")) 

unique(mass_yr1$site_code)

live_mass_yr1 <- mass_yr1 %>% group_by(site_code, plot,subplot, year,live, trt) %>% summarize(mass = sum(mass),.groups = 'keep')  %>% filter(live == 1)

mass_use <- live_mass_yr1 %>% left_join(.,comb_yr1 %>% 
                                          select(site_code,plot,year,n_drought_yrs,n_nutrient_yrs,trt_use,
                                                 MAT_v2,MAP_v2,AI,MAP_VAR_v2,site_year_rich, precip_sum_365, precip_sum_365_drt,
                                                 map,map_sd,,map_cv,ctrl_percentile_365,drt_percentile_365,
                                                 mean90, sd90,
                                                 precip_sum_90,precip_sum_90_drt,ctrl_percentile_90, drt_percentile_90),
                                        by = c('site_code','plot','subplot','year')) %>% 
  mutate(fert = if_else(trt_use %in% c('Fertilization+Drought','Fertilization'),1,0),
         drought= if_else(trt_use %in% c('Fertilization+Drought','Drought'),1,0)) %>% filter(live == 1) %>% 
  mutate(drt_sev = (precip_sum_365_drt-map)/map,
         drt_sev_90 = (precip_sum_90_drt-mean90)/mean90) %>% 
  left_join(.,massM,
            by = c('site_code','live','trt')) %>% 
  left_join(.,massG,
            by=c('site_code','year'))


mass_use <- mass_use %>% 
  #filter(trt == 'Control') %>% 
  group_by(site_code) %>% 
  mutate(ctrl_mass = mean(mass[trt_use == 'Control']),
         fert_mass = mean(mass[trt_use == 'Fertilization']),
         drt_mass = mean(mass[trt_use == 'Drought']),
         df_mass = mean(mass[trt_use == 'Fertilization+Drought']),
         massM_ctrl = mean(massM[trt_use == 'Control'])) %>% 
  mutate(LRR_mass = log(mass/ctrl_mass)) %>% 
  mutate(ppt = if_else(precip_sum_365 < map,'Below average','Above average'),
         ppt90 = if_else(precip_sum_90 < mean90, 'Below average','Above average')) %>% 
  mutate(AR = if_else(trt_use == 'Fertilization+Drought',abs((fert_mass - mass)/massM_ctrl), abs((ctrl_mass - mass)/massM_ctrl)),
         PUEm = abs((precip_sum_365-precip_sum_365_drt)/map)) %>% 
  mutate(RM = AR/PUEm)

unique(mass_use$site_code)

# Hypothesis 1 - treatment effects on (a) biomass #####
# ~~~  (a) biomass #####


#mod_mass_yr1 <- lme(LRR_mass ~ drought*fert, random = ~1|site_code, data = mass_use)
mod_mass_yr1 <- lme(LRR_mass ~ drought*fert*ppt, random = ~1|site_code, data = mass_use)
#mod_mass_yr1_nout <- lme(LRR_mass ~ drought*fert*ppt, random = ~1|site_code, data = mass_use %>% filter(site_code != 'urat.cn'))
summary(mod_mass_yr1)
anova(mod_mass_yr1)
#confint(mod_mass_yr1)

#anova(mod_mass_yr1_nout) #without urat severe drought flips to more significant


em.mass <- data.frame(emmeans(mod_mass_yr1, ~drought*fert))

em.mass_alt <- em.mass %>% 
  mutate(trt_use = if_else(fert == 0 & drought == 0, 'Control',
                           if_else(fert == 0 & drought == 1, 'Drought',
                                   if_else(fert == 1 & drought == 0, 'Fertilization',
                                           'Fertilization+Drought')))) %>% 
  filter(trt_use != 'Control')

pd <- position_dodge(0.3)

range(mass_use$LRR_mass)

gg_mass <- ggplot(em.mass, aes(x = as.factor(fert), y = emmean, col = as.factor(drought))) + geom_hline(yintercept=log(1), linetype = 2) + 
  geom_point(position = pd, size = 3) + 
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL),width=0.3,position=pd, size = 2) +
  #geom_jitter(comb_npk_max)
  scale_color_manual(values= c('gray25','firebrick3'),guide = guide_legend(title='Drought Treatment'),
                     labels= c('0' = 'Control','1'='Drought')) +
  ylab('Biomass LRR') + xlab('') + scale_x_discrete(limits = rev(levels(factor(em.mass$fert))),labels = c('Fertilized', 'Ambient')) +coord_flip()
   #+  scale_y_continuous(breaks = c(log(0.75),log(1),log(1.25)), labels = c(0.75,1,1.25)) # to relabel axis

gg_mass

exp(-0.22028671) # 0.8022887
exp(0.26899591) # 1.30865

gg_mass_alt <- ggplot(em.mass_alt, aes(x = trt_use, y = emmean, col = as.factor(trt_use))) + geom_hline(yintercept=log(1), linetype = 2) + 
  geom_point(position = pd, size = 3) + 
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL),width=0.3,position=pd, size = 2) +
  #geom_jitter(comb_npk_max)
  scale_color_manual(values= col_pal,guide = 'none') +
  ylab('Biomass LRR') + xlab('') + scale_x_discrete(limits = rev(levels(factor(em.mass_alt$trt_use)))) +coord_flip()
#+  scale_y_continuous(breaks = c(log(0.75),log(1),log(1.25)), labels = c(0.75,1,1.25)) # to relabel axis

gg_mass_alt

figwd <- ('/Users/wilf0020/Library/Mobile Documents/com~apple~CloudDocs/Documents/NPKD Net/npkdnet/Year 1 figures')

# ggsave(paste0(figwd,'/Fig-1_mass-LRR_',Sys.Date(),'.png'),
#        gg_mass_alt,'png',width = 12, height = 6)

# ~~~  (b) sensitivity #####

#mod_mass_yr1 <- lme(LRR_mass ~ drought*fert, random = ~1|site_code, data = mass_use)
mod_rm_yr1 <- lme(RM ~ trt_use, random = ~1|site_code, data = mass_use %>% filter(trt_use %in% c('Drought','Fertilization+Drought')))
#mod_mass_yr1_nout <- lme(LRR_mass ~ drought*fert*ppt, random = ~1|site_code, data = mass_use %>% filter(site_code != 'urat.cn'))
summary(mod_rm_yr1)
anova(mod_rm_yr1)
#confint(mod_mass_yr1)

#anova(mod_mass_yr1_nout) #without urat severe drought flips to more significant


em.rm <- data.frame(emmeans(mod_rm_yr1, ~trt_use))

# em.rm_alt <- em.rm %>% 
#   mutate(trt_use = if_else(fert == 0 & drought == 0, 'Control',
#                            if_else(fert == 0 & drought == 1, 'Drought',
#                                    if_else(fert == 1 & drought == 0, 'Fertilization',
#                                            'Fertilization+Drought')))) %>% 
#   filter(trt_use != 'Control')


gg_rm <- ggplot(em.rm, aes(x = trt_use, y = emmean, col = as.factor(trt_use))) + geom_hline(yintercept=log(1), linetype = 2) + 
  geom_point(position = pd, size = 3) + 
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL),width=0.2,position=pd, size = 2) +
  #geom_jitter(comb_npk_max)
  scale_color_manual(values= c("#D53E4F","#5E4FA2"),guide = 'none') +
  ylab('Drought Sensitivity') + xlab('') #+ scale_x_discrete(limits = rev(levels(factor(em.mass$fert))),labels = c('Fertilized', 'Ambient')) +coord_flip()
#+  scale_y_continuous(breaks = c(log(0.75),log(1),log(1.25)), labels = c(0.75,1,1.25)) # to relabel axis

gg_rm


# ggsave(paste0(figwd,'/Fig-1-alt_sens-index_',Sys.Date(),'.png'),
#        gg_rm,'png',width = 9, height = 6)

# ~~~  (c) precipitation effect #####

ppt_less <- unique(mass_use[mass_use$ppt == 'Below average',]$site_code)
ppt_more <- unique(mass_use[mass_use$ppt == 'Above average',]$site_code)



em.mass_ppt <- data.frame(emmeans(mod_mass_yr1, ~drought*fert*ppt)) %>% filter((drought == 1 & fert == 0) | (fert == 1 & drought == 0)) %>% 
   mutate(trt = if_else(drought == 1,'Drought','Fertilization'))

pd <- position_dodge(0.3)


em.mass_ppt_alt <- data.frame(emmeans(mod_mass_yr1, ~drought*fert*ppt)) %>% 
  mutate(trt_use = if_else(fert == 0 & drought == 0, 'Control',
                                  if_else(fert == 0 & drought == 1, 'Drought',
                                          if_else(fert == 1 & drought == 0, 'Fertilization',
                                                  'Fertilization+Drought')))) %>% 
           filter(trt_use != 'Control')

gg_mass_ppt <- ggplot(em.mass_ppt, aes(x = trt, y = emmean, col = as.factor(ppt))) + geom_hline(yintercept=log(1), linetype = 2) + 
  geom_point(position = pd, size = 3) + 
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL),width=0.3,position=pd, size = 2) +
  #geom_jitter(comb_npk_max)
  scale_color_manual(values= c('gray25','#66CC99'),guide = guide_legend(title='Precipitation')) +
  ylab('Biomass LRR') + xlab('') + 
  #scale_x_discrete(limits = rev(levels(factor(em.mass$fert))),labels = c('Fertilized', 'Ambient')) +
  coord_flip()
#+  scale_y_continuous(breaks = c(log(0.75),log(1),log(1.25)), labels = c(0.75,1,1.25)) # to relabel axis

gg_mass_ppt


gg_mass_ppt_alt <- ggplot(em.mass_ppt_alt, aes(x = trt_use, y = emmean, col = as.factor(ppt))) + geom_hline(yintercept=log(1), linetype = 2) + 
  geom_point(position = pd, size = 3) + 
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL),width=0.3,position=pd, size = 2) +
  #geom_jitter(comb_npk_max)
  scale_color_manual(values= c('gray25','#66CC99'),guide = guide_legend(title='Precipitation')
                    ) +
  ylab('Biomass LRR') + xlab('') + 
  scale_x_discrete(limits = rev(levels(factor(em.mass_ppt_alt$trt_use)))) +
  coord_flip()
#+  scale_y_continuous(breaks = c(log(0.75),log(1),log(1.25)), labels = c(0.75,1,1.25)) # to relabel axis

gg_mass_ppt_alt

# ggsave(paste0(figwd,'/Fig-1_mass-LRR-precip_',Sys.Date(),'.png'),
#        gg_mass_ppt_alt,'png',width = 12, height = 6)


# site specific results ####

RRdata <- mass_use %>%
  filter(trt_use %in% c ('Control','Fertilization','Drought','Fertilization+Drought')) %>% 
  group_by(site_code) %>% 
  mutate(LRR_mass = if_else(trt_use == "Fertilization+Drought", 
                            LRR_mass - mean(LRR_mass[trt_use == 'Fertilization']) - mean(LRR_mass[trt_use == 'Drought']),
                            LRR_mass)) %>% 
  dplyr::group_by(site_code,trt_use) %>%
  dplyr::summarize(mean_lrr = mean(LRR_mass,na.rm=T),n=n(), sd = sd(LRR_mass,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_lrr - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_lrr + qt(1 - (0.05 / 2), n - 1) * se) %>%
  dplyr::as_tibble()

### ~stitch plots ####
unique(RRdata$site_code)

RRdata.drt <- RRdata %>% filter(trt_use == 'Drought') %>% arrange(mean_lrr) %>% 
mutate(site_code = factor(site_code, levels=unique(site_code))) %>% 
  mutate(ppt = if_else(site_code %in% ppt_less,'Extreme','Nominal'))

gg_site_drt <- ggplot(RRdata.drt,aes(x = mean_lrr,y=site_code, col = ppt))+
  geom_vline(xintercept=0,size=1.5,lty=2)+
  geom_segment(aes(x=LowerCI,y=site_code,xend=UpperCI,yend=site_code),size=1)+
  geom_point(size=6)+
  scale_color_manual(values= c('#66CC99','gray25'),guide = 'none') +
 #guide_legend(title='Precipitation'),  labels= c('Less' = 'Below Average','More'='Above Average')) +
  labs(x = "Drought Response - LRR", y="")+
  theme(axis.text = element_text(size = 15)) +
  geom_hline(yintercept=0)

gg_site_drt
  
RRdata.fert <- RRdata %>% filter(trt_use == 'Fertilization') %>% arrange(mean_lrr) %>% 
mutate(site_code = factor(site_code, levels=unique(site_code))) %>% 
  mutate(ppt = if_else(site_code %in% ppt_less,'Extreme','Nominal'))

gg_site_fert <- ggplot(RRdata.fert ,aes(x = mean_lrr,y=site_code, col = ppt))+
  geom_vline(xintercept=0,size=1.5,lty=2)+
  geom_segment(aes(x=LowerCI,y=site_code,xend=UpperCI,yend=site_code),size=1)+
  geom_point(size=6)+
  scale_color_manual(values= c('#66CC99','gray25'),guide = 'none') +
                     #   guide_legend(title='Precipitation'),
                     # labels= c('Extreme' = 'Below Average','Nominal'='Above Average')) +
  theme(axis.text = element_text(size = 15)) +
  labs(x = "Fertilization Response - LRR", y="")+
  geom_hline(yintercept=0)

gg_site_fert

RRdata.comb <- RRdata %>% filter(trt_use == 'Fertilization+Drought') %>% arrange(mean_lrr) %>% 
  mutate(site_code = factor(site_code, levels=unique(site_code))) %>% 
  mutate(ppt = if_else(site_code %in% ppt_less,'Extreme','Nominal'))

gg_site_comb <- ggplot(RRdata.comb ,aes(x = mean_lrr,y=site_code, col = ppt))+
  geom_vline(xintercept=0,size=1.5,lty=2)+
  geom_segment(aes(x=LowerCI,y=site_code,xend=UpperCI,yend=site_code),size=1)+
  geom_point(size=6)+
  scale_color_manual(values= c('#66CC99','gray25'),guide = guide_legend(title='Precipitation'),
                     labels= c('Extreme' = 'Below Average','Nominal'='Above Average')) +
  labs(x = "Drought + Fertilization Response (LRR)", y="")+
  theme(axis.text = element_text(size = 15)) +
  geom_hline(yintercept=0)

gg_site_comb 
# 
# ggsave(paste0(figwd,'/Fig-1_stitch-plots_',Sys.Date(),'.png'),
# plot_grid(gg_site_drt,gg_site_fert,gg_site_comb,ncol = 3,rel_widths = c(0.28,0.28,0.44)),
# 'png',width = 20, height = 6)

#### site summary graph ####

site_summary <- bind_rows(
  RRdata.drt %>% 
    mutate(sig = if_else(LowerCI * UpperCI >= 0, 'Yes','No'),
           effect = if_else(mean_lrr < 0, 'Negative','Positve'),
           sig_effect = paste0(sig,effect)) %>% 
    select(site_code,trt_use,mean_lrr,sig,effect,sig_effect),
  RRdata.fert %>% 
    mutate(sig = if_else(LowerCI * UpperCI >= 0, 'Yes','No'),
           effect = if_else(mean_lrr < 0, 'Negative','Positve'),
           sig_effect = paste0(sig,effect)) %>% 
    select(site_code,trt_use,mean_lrr,sig,effect,sig_effect),
  RRdata.comb %>% 
    mutate(sig = if_else(LowerCI * UpperCI >= 0, 'Yes','No'),
           effect = if_else(mean_lrr < 0, 'Negative','Positve'),
           sig_effect = paste0(sig,effect)) %>% 
    select(site_code,trt_use,mean_lrr,sig,effect,sig_effect)
) %>% 
  mutate(trt_use = factor(trt_use, levels = c('Drought','Fertilization','Fertilization+Drought')))


site_summary

gg_summary <- ggplot(site_summary,aes(x = trt_use, y = site_code, fill = mean_lrr)) + geom_tile(linewidth=.5,color = 'black') +
  scale_fill_gradient2(low = 'red',high = 'blue',guide = guide_legend(title='Direction of effect')) +
  ylab('') + xlab('') +
  coord_equal(ratio = .50) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(labels=c("Drought" = "D", "Fertilization" = "F", "Fertilization+Drought" = "D+F"))

site_summary_wide <- site_summary %>% pivot_wider(id_cols = site_code, values_from = mean_lrr,names_from = trt_use)

summary(lm(Drought ~ Fertilization, data=site_summary_wide)) #%>% filter(site_code != 'urat.cn')))
# not sig with or without urat

ggsave(paste0(figwd,'/Fig-1b-alt_site-summary_',Sys.Date(),'.png'),
  gg_summary, 'png',width = 9, height = 16)

### ~~graph a negative interaction sites ####

mass_ycn <- mass_use %>% filter(site_code == 'yanchi.cn')

lm_ycn <- lm(mass ~ drought*fert, mass_ycn)
summary(lm_ycn)

em_ycn <- data.frame(emmeans(lm_ycn, ~drought*fert))

em_ycn <- em_ycn %>% 
  mutate(trt_use = if_else(fert == 0 & drought == 0, 'Control',
                           if_else(fert == 0 & drought == 1, 'Drought',
                                   if_else(fert == 1 & drought == 0, 'Fertilization',
                                           'Fertilization+Drought'))))

em_ycn


gg_ycn <- ggplot(em_ycn, aes(x = as.factor(fert), y = emmean, col = as.factor(drought))) + 
  geom_point(position = pd, size = 3) + 
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL),width=0.3,position=pd, size = 2) +
  #geom_jitter(comb_npk_max)
  scale_color_manual(values= c('gray25','firebrick3'),guide = guide_legend(title='Drought Treatment'),
                     labels= c('0' = 'Control','1'='Drought')) +
  ylab('Biomass (g/m2)') + xlab('') +
  scale_x_discrete(labels = c('Ambient','Fertilized'))
# scale_y_continuous(breaks = c(log(0.75),log(1),log(1.25)), labels = c(0.75,1,1.25)) # to relabel axis

gg_ycn

#biomass across gradients ####

mod_mass_grad <- lme(LRR_mass ~ drought*fert*(map + site_year_rich + map_cv + AI + avg_gram_prop + drt_percentile_365) - 1, random = ~1|site_code, data = mass_use)
summary(mod_mass_grad)
anova(mod_mass_grad)

tab_mass_grad <- summary.tablefunc(mod_mass_grad)

# write_csv(tab_mass_grad,
#           paste0(figwd,'/mass-gradients-table_',
#                  Sys.Date(),'.csv'))

mass_graph <- mass_use %>% 
  group_by(site_code,trt_use) %>% 
  summarize(mean_lrr = mean(LRR_mass,na.rm=T),n=n(), sd = sd(LRR_mass,na.rm=T), se = sd/sqrt(n),
            LowerCI = mean_lrr - qt(1 - (0.05 / 2), n - 1) * se,
            UpperCI = mean_lrr + qt(1 - (0.05 / 2), n - 1) * se, .groups = 'keep') %>% 
  left_join(.,mass_use %>% distinct(site_code,ppt,map,site_year_rich,drt_percentile_365,map_cv,avg_gram_prop,AI))

# mod_mass_ppt_noout <- lme(LRR_mass ~ map*drought*fert*ppt, random = ~1|site_code, data = mass_use %>% 
#                             filter(!site_code == 'urat.cn'))
# summary(mod_mass_ppt_noout)
# anova(mod_mass_ppt_noout)

# ~~~~ precipitation ####
em_mass_map <- data.frame(summary(emmeans(mod_mass_grad, pairwise ~ fert * drought | map,
                                          at = list(map = seq(from = min(mass_use$map), to =  max(mass_use$map), by = (max(mass_use$map)-min(mass_use$map))/20)),
                                          cov.reduce = TRUE,
                                          type = "response")$emmeans))

tab_mass_grad %>% filter(grepl("map", Effect) & !grepl('map_cv', Effect))
# Fertilization is significantly negative

em_mass_map

em_mass_map <- em_mass_map %>%
  mutate(trt_use = if_else(fert == 0 & drought == 0, 'Control',
                       if_else(fert == 0 & drought == 1, 'Drought',
                               if_else(fert == 1 & drought == 0, 'Fertilization',
                                       'Fertilization+Drought')))) %>% 
  filter(trt_use != 'Control') %>% 
  mutate(sig = if_else(trt_use == 'Fertilization','Yes','No'))


gg_map <- ggplot(mass_graph %>% filter(trt_use != 'Control'), aes(x = map, y = mean_lrr, col = trt_use)) +
  geom_hline(yintercept=0, linetype = 'dotted') +
  geom_point(alpha = 0.15) + geom_segment(aes(x = map, xend = map, y=LowerCI,yend=UpperCI),size=1, alpha = 0.15) +
  geom_ribbon(data = em_mass_map,aes(x = map, y= emmean,ymin = lower.CL,ymax= upper.CL, fill = trt_use),col = NA,alpha=0.35)+
  #geom_line(data = em_mass_map,aes(x = map, y= emmean), size = 2, color= 'black') +
  geom_line(data = em_mass_map,aes(x = map, y= emmean, linetype = sig, alpha = sig), size = 1.8) +
  scale_color_manual(values= col_pal, guide = "none") +
  scale_fill_manual(values= col_pal, guide = 'none') +
  scale_linetype_manual(values = c('dashed','solid'), guide = 'none') +
  scale_alpha_manual(values = c(0.35,1), guide = 'none') +
  labs(y = "", x="Mean annual precipitation")
  
gg_map

# ~~~~ precipitation variability ####
em_mass_mapcv <- data.frame(summary(emmeans(mod_mass_grad, pairwise ~ fert * drought | map_cv,
                                          at = list(map_cv = seq(from = min(mass_use$map_cv), to =  max(mass_use$map_cv), by = (max(mass_use$map_cv)-min(mass_use$map_cv))/20)),
                                          cov.reduce = TRUE,
                                          type = "response")$emmeans))

tab_mass_grad %>% filter(grepl('map_cv', Effect))
# Drought and Fertilization are significantly negative

em_mass_mapcv

em_mass_mapcv <- em_mass_mapcv %>%
  mutate(trt_use = if_else(fert == 0 & drought == 0, 'Control',
                           if_else(fert == 0 & drought == 1, 'Drought',
                                   if_else(fert == 1 & drought == 0, 'Fertilization',
                                           'Fertilization+Drought')))) %>% 
  filter(trt_use != 'Control') %>% 
  mutate(sig = if_else(trt_use == 'Fertilization+Drought','No','Yes'))

gg_map_cv <- ggplot(mass_graph %>% filter(trt_use != 'Control'), aes(x = map_cv, y = mean_lrr, col = trt_use)) +
  geom_hline(yintercept=0, linetype = 'dotted') +
  geom_point(alpha = 0.15) + geom_segment(aes(x = map_cv, xend = map_cv, y=LowerCI,yend=UpperCI),size=1, alpha = 0.15) +
  geom_ribbon(data = em_mass_mapcv,aes(x = map_cv, y= emmean,ymin = lower.CL,ymax= upper.CL, fill = trt_use),col = NA,alpha=0.35)+
  #geom_line(data = em_mass_mapcv,aes(x = map_cv, y= emmean), size = 2, color= 'black') +
  geom_line(data = em_mass_mapcv,aes(x = map_cv, y= emmean, linetype = sig, alpha = sig), size = 1.8) +
  scale_color_manual(values= col_pal, guide = "none") +
  scale_fill_manual(values= col_pal, guide = 'none') +
  scale_linetype_manual(values = c('dashed','solid'), guide = 'none') +
  scale_alpha_manual(values = c(0.35,1), guide = 'none') +
  labs(y = "", x="CV of mean annual precipitation")

gg_map_cv

# ~~~~ site richness ####
em_mass_rich <- data.frame(summary(emmeans(mod_mass_grad, pairwise ~ fert * drought | site_year_rich,
                                            at = list(site_year_rich = seq(from = min(mass_use$site_year_rich), to =  max(mass_use$site_year_rich), by = (max(mass_use$site_year_rich)-min(mass_use$site_year_rich))/20)),
                                            cov.reduce = TRUE,
                                            type = "response")$emmeans))

tab_mass_grad %>% filter(grepl('site_year_rich', Effect))
# Fertilization is significantly positive

em_mass_rich

em_mass_rich <- em_mass_rich %>%
  mutate(trt_use = if_else(fert == 0 & drought == 0, 'Control',
                           if_else(fert == 0 & drought == 1, 'Drought',
                                   if_else(fert == 1 & drought == 0, 'Fertilization',
                                           'Fertilization+Drought')))) %>% 
  filter(trt_use != 'Control') %>% 
  mutate(sig = if_else(trt_use == 'Fertilization','Yes','No'))

gg_site_year_rich <- ggplot(mass_graph %>% filter(trt_use != 'Control'), aes(x = site_year_rich, y = mean_lrr, col = trt_use)) +
  geom_hline(yintercept=0, linetype = 'dotted') +
  geom_point(alpha = 0.15) + geom_segment(aes(x = site_year_rich, xend = site_year_rich, y=LowerCI,yend=UpperCI),size=1, alpha = 0.15) +
  geom_ribbon(data = em_mass_rich,aes(x = site_year_rich, y= emmean,ymin = lower.CL,ymax= upper.CL, fill = trt_use),col = NA,alpha=0.35)+
  #geom_line(data = em_mass_rich,aes(x = site_year_rich, y= emmean), size = 2, color= 'black') +
  geom_line(data = em_mass_rich,aes(x = site_year_rich, y= emmean, linetype = sig, alpha = sig), size = 1.8) +
  scale_color_manual(values= col_pal, guide = "none") +
  scale_fill_manual(values= col_pal, guide = 'none') +
  scale_linetype_manual(values = c('dashed','solid'), guide = 'none') +
  scale_alpha_manual(values = c(0.35,1), guide = 'none') +
  labs(y = "", x="Site richness")

gg_site_year_rich

# ~~~~ proportion of graminoids ####
em_mass_gram <- data.frame(summary(emmeans(mod_mass_grad, pairwise ~ fert * drought | avg_gram_prop,
                                           at = list(avg_gram_prop = seq(from = min(mass_use$avg_gram_prop), to =  max(mass_use$avg_gram_prop), by = (max(mass_use$avg_gram_prop)-min(mass_use$avg_gram_prop))/20)),
                                           cov.reduce = TRUE,
                                           type = "response")$emmeans))

tab_mass_grad %>% filter(grepl('avg_gram_prop', Effect))
# Drought is significantly positive

em_mass_gram

em_mass_gram <- em_mass_gram %>%
  mutate(trt_use = if_else(fert == 0 & drought == 0, 'Control',
                           if_else(fert == 0 & drought == 1, 'Drought',
                                   if_else(fert == 1 & drought == 0, 'Fertilization',
                                           'Fertilization+Drought')))) %>% 
  filter(trt_use != 'Control') %>% 
  mutate(sig = if_else(trt_use == 'Drought','Yes','No'))

gg_avg_gram_prop <- ggplot(mass_graph %>% filter(trt_use != 'Control'), aes(x = avg_gram_prop, y = mean_lrr, col = trt_use)) +
  geom_hline(yintercept=0, linetype = 'dotted') +
  geom_point(alpha = 0.15) + geom_segment(aes(x = avg_gram_prop, xend = avg_gram_prop, y=LowerCI,yend=UpperCI),size=1, alpha = 0.15) +
  geom_ribbon(data = em_mass_gram,aes(x = avg_gram_prop, y= emmean,ymin = lower.CL,ymax= upper.CL, fill = trt_use), col = NA, alpha = 0.35)+
  #geom_line(data = em_mass_gram,aes(x = avg_gram_prop, y= emmean), size = 2, color= 'black') +
  geom_line(data = em_mass_gram,aes(x = avg_gram_prop, y= emmean, linetype = sig, alpha = sig), size = 1.8) +
  scale_color_manual(values= col_pal, guide = "none") +
  scale_fill_manual(values= col_pal, guide = 'none') +
  scale_linetype_manual(values = c('dashed','solid'), guide = 'none') +
  scale_alpha_manual(values = c(0.35,1), guide = 'none') +
  labs(y = "", x="Graminoid proportion")

gg_avg_gram_prop

#### Gradient figure ####
#### stitch them together

rank_legend <- get_legend(ggplot(em_mass_map %>% filter(trt_use != 'Control'), aes(x = map, y = emmean, col = trt_use)) +
                            scale_color_manual(values= col_pal, guide = guide_legend(title='Treatment')) + geom_line(size = 2))



blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

grad_fig <- grid.arrange(rank_legend,blankPlot,
                         gg_map, gg_map_cv,
                         gg_site_year_rich,gg_avg_gram_prop,
                          ncol = 2,nrow=3,
                          widths = c(3.5,3.5), heights = c(1.3,3.5,3.5),
                          left = textGrob("Treatment LRR", rot = 90, gp=gpar(fontsize=24))
                          )


# save_plot(paste0(figwd,'/Fig-2_treatment-gradients_',Sys.Date(),'.png'),
#           grad_fig,
#           base_height = 12, base_width = 12)


# biomass across richness gradient ####

mod_mass_rich <- lme(LRR_mass ~ site_year_rich*drought*fert, random = ~1|site_code, data = mass_use)
summary(mod_mass_rich)
anova(mod_mass_rich)

em_mass_rich <- data.frame(summary(emmeans(mod_mass_rich, pairwise ~ fert * drought | site_year_rich ,
                                          at = list(site_year_rich = seq(from = min(mass_use$site_year_rich), to =  max(mass_use$site_year_rich), by =max(mass_use$site_year_rich)-min(mass_use$site_year_rich))),
                                          cov.reduce = TRUE,
                                          type = "response")$emmeans))

em_mass_rich

em_mass_rich <- em_mass_rich %>%
  mutate(trt_use = if_else(fert == 0 & drought == 0, 'Control',
                           if_else(fert == 0 & drought == 1, 'Drought',
                                   if_else(fert == 1 & drought == 0, 'Fertilization',
                                           'Fertilization+Drought')))) %>% 
  filter(trt_use != 'Control')


# gg_rich_drt <- ggplot(mass_graph %>% filter(trt_use == 'Drought'), aes(x = site_year_rich, y = mean_lrr, col = ppt)) +
#   geom_point() + geom_segment(aes(x = site_year_rich, xend = site_year_rich, y=LowerCI,yend=UpperCI),size=1, alpha = 0.5) +
#   geom_ribbon(data = em_mass_rich %>% filter(trt_use == 'Drought'),aes(x = site_year_rich, y= emmean,ymin = lower.CL,ymax= upper.CL,fill=ppt),alpha=0.25)+
#   #geom_line(data = em_mass_rich,aes(x = map, y= emmean), size = 2, color= 'black') +
#   geom_line(data = em_mass_rich %>% filter(trt_use == 'Drought'),aes(x = site_year_rich, y= emmean,col = ppt), size = 1.8) +
#   scale_color_manual(values= c('#66CC99','gray25')) + guides(colour = "none") +
#   labs(y = "Drought Response (LRR)", x="Site richness")+
#   geom_hline(yintercept=0)
# 
# gg_rich_drt

gg_rich_fert <- ggplot(mass_graph %>% filter(trt_use == 'Fertilization'), aes(x = site_year_rich, y = mean_lrr)) +
  geom_point() + geom_segment(aes(x = site_year_rich, xend = site_year_rich, y=LowerCI,yend=UpperCI),size=1, alpha = 0.5) +
  geom_ribbon(data = em_mass_rich %>% filter(trt_use == 'Fertilization'),aes(x = site_year_rich, y= emmean,ymin = lower.CL,ymax= upper.CL),alpha=0.25)+
  #geom_line(data = em_mass_rich,aes(x = map, y= emmean), size = 2, color= 'black') +
  geom_line(data = em_mass_rich %>% filter(trt_use == 'Fertilization'),aes(x = site_year_rich, y= emmean), size = 1.8) +
  scale_color_manual(values= c('#66CC99','gray25')) + guides(colour = "none") +
  labs(y = "Fertilization Response (LRR)", x="Site richness")+
  geom_hline(yintercept=0)

gg_rich_fert



#### Drought severity analysis ####


mod_mass_drt <- lme(LRR_mass ~ drt_sev*drought*fert, random = ~1|site_code, data = mass_use)
summary(mod_mass_drt)
anova(mod_mass_drt)

em_mass_drt <- data.frame(summary(emmeans(mod_mass_drt, pairwise ~ fert * drought | drt_sev,
                                           at = list(drt_sev = seq(from = min(mass_use$drt_sev), to =  max(mass_use$drt_sev), by =max(mass_use$drt_sev)-min(mass_use$drt_sev))),
                                           cov.reduce = TRUE,
                                           type = "response")$emmeans))

em_mass_drt

em_mass_drt <- em_mass_drt %>%
  mutate(trt_use = if_else(fert == 0 & drought == 0, 'Control',
                           if_else(fert == 0 & drought == 1, 'Drought',
                                   if_else(fert == 1 & drought == 0, 'Fertilization',
                                           'Fertilization+Drought')))) %>% 
  filter(trt_use != 'Control')

gg_mass_drt <- ggplot(mass_graph %>% filter(trt_use == 'Drought'), aes(x = drt_sev, y = mean_lrr, col = ppt)) +
  geom_point() + geom_segment(aes(x = drt_sev, xend = drt_sev, y=LowerCI,yend=UpperCI),size=1, alpha = 0.5) +
  geom_ribbon(data = em_mass_drt %>% filter(trt_use == 'Drought'),aes(x = drt_sev, y= emmean,ymin = lower.CL,ymax= upper.CL), col='black',alpha=0.25)+
  #geom_line(data = em_mass_rich,aes(x = map, y= emmean), size = 2, color= 'black') +
  geom_line(data = em_mass_drt %>% filter(trt_use == 'Drought'),aes(x = drt_sev, y= emmean), col='black', size = 1.8) +
  scale_color_manual(values= c('gray25','#66CC99')) + guides(colour = "none") +
  # scale_fill_manual(values= c('#66CC99','gray25'),guide = guide_legend(title='Drought effect'),
  #                   labels= c('Less' = 'Extreme','More'='Nominal')) +
  labs(y = "Drought Response (LRR)", x="Drought Severity")+
  geom_hline(yintercept=0)

gg_mass_drt

gg_mass_fert <- ggplot(mass_graph %>% filter(trt_use == 'Fertilization'), aes(x = drt_sev, y = mean_lrr, col = ppt)) +
  geom_point() + geom_segment(aes(x = drt_sev, xend = drt_sev, y=LowerCI,yend=UpperCI),size=1, alpha = 0.5) +
  geom_ribbon(data = em_mass_drt %>% filter(trt_use == 'Fertilization'),aes(x = drt_sev, y= emmean,ymin = lower.CL,ymax= upper.CL), col='black',alpha=0.25)+
  #geom_line(data = em_mass_rich,aes(x = map, y= emmean), size = 2, color= 'black') +
  geom_line(data = em_mass_drt %>% filter(trt_use == 'Fertilization'),aes(x = drt_sev, y= emmean), col='black', size = 1.8) +
  scale_color_manual(values= c('gray25','#66CC99')) + guides(colour = "none") +
  labs(y = "Fertilization Response (LRR)", x="Drought Severity")+
  geom_hline(yintercept=0)

gg_mass_fert



  
### drought sensitivity analysis ####
  #per Bondurak et al 2022 (https://doi.org/10.1111/avsc.12674)

colnames(mass_use)

mod_drt_sens <- lme(RM ~ trt_use*(map + map_cv + site_year_rich + avg_gram_prop + AI), random = ~1|site_code, data = mass_use %>% filter(trt_use %in% c('Drought','Fertilization+Drought')))
summary(mod_drt_sens)
anova(mod_drt_sens)


em_drt_sens_map <- data.frame(summary(emmeans(mod_drt_sens, pairwise ~ trt_use | map, 
                                          at = list(map = seq(from = min(mass_use$map), to =  max(mass_use$map), by = (max(mass_use$map)-min(mass_use$map))/50)),
                                  type = "response")$emmeans))

em_drt_sens_map

mass_graph_sense <- mass_use %>% 
  group_by(site_code,trt_use) %>% 
  summarize(mean_RM = mean(RM,na.rm=T),n=n(), sd = sd(RM,na.rm=T), se = sd/sqrt(n),
            LowerCI = mean_RM - qt(1 - (0.05 / 2), n - 1) * se,
            UpperCI = mean_RM + qt(1 - (0.05 / 2), n - 1) * se, .groups = 'keep') %>% 
  left_join(.,mass_use %>% distinct(site_code,trt_use,map,site_year_rich)) %>% 
  filter(trt_use %in% c('Drought','Fertilization+Drought'))

gg_drt_sens_map <- ggplot(em_drt_sens_map, aes(x = map, y = emmean, col = trt_use, fill = trt_use)) +
  geom_point(mass_graph_sense, mapping = aes(x = map, y = mean_RM, col = trt_use),position=pd,alpha = 0.25) + 
  #geom_segment(mass_graph_sense, mapping = aes(x = map, xend = map, y=LowerCI,yend=UpperCI),size=1, position=pd,alpha = 0.5) +
  geom_ribbon(aes(x = map, y= emmean,ymin = lower.CL,ymax= upper.CL),alpha=0.25)+
  #geom_line(data = em_mass_rich,aes(x = map, y= emmean), size = 2, color= 'black') +
  geom_line(size = 1.8) +
  scale_color_manual(values= c("#D53E4F","#5E4FA2")) + guides(colour = "none") +
  scale_fill_manual(values= c("#D53E4F","#5E4FA2"),guide = guide_legend(title='Treatment'),
                   ) +
  labs(y = "Drought Sensitivity", x="Mean Annual Precipitation (mm)")+
  geom_hline(yintercept=0, linetype = 'dotted')

gg_drt_sens_map

mass_graph_sense %>% filter(site_code == 'sand.us')


em_drt_sens_rich <- data.frame(summary(emmeans(mod_drt_sens, pairwise ~ trt_use | site_year_rich, 
                                              at = list(site_year_rich = seq(from = min(mass_use$site_year_rich), to =  max(mass_use$site_year_rich), 
                                                                             by = (max(mass_use$site_year_rich)-min(mass_use$site_year_rich))/50)),
                                              type = "response")$emmeans))

em_drt_sens_rich


gg_drt_sens_rich <- ggplot(em_drt_sens_rich, aes(x = site_year_rich, y = emmean, col = trt_use, fill = trt_use)) +
  geom_point(mass_graph_sense, mapping = aes(x = site_year_rich, y = mean_RM, col = trt_use),position=pd,alpha = 0.25) + 
  #geom_segment(mass_graph_sense, mapping = aes(x = map, xend = map, y=LowerCI,yend=UpperCI),size=1, position=pd,alpha = 0.5) +
  geom_ribbon(aes(x = site_year_rich, y= emmean,ymin = lower.CL,ymax= upper.CL),alpha=0.25)+
  #geom_line(data = em_mass_rich,aes(x = map, y= emmean), size = 2, color= 'black') +
  geom_line(size = 1.8) +
  scale_color_manual(values= c("#D53E4F","#5E4FA2")) + guides(colour = "none") +
  scale_fill_manual(values= c("#D53E4F","#5E4FA2"),guide = guide_legend(title='Treatment'),
  ) +
  labs(y = "Drought Sensitivity", x="Site Richness")+
  geom_hline(yintercept=0, linetype = 'dotted')

gg_drt_sens_rich
