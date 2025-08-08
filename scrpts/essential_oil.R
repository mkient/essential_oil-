
### Load packages ####
library(tidyverse)
library(ecotox)
library(emmeans)
library(reshape)
library(drc)
#library(multcomp)


### Read data #### 
setwd("C:/Farm/Papiers/ecotox/huilles_veg/")
#list.files()
df_tube <- read.table('data_df.csv', sep=",", header=TRUE, dec=".",
                      na.string=NA, quote = "\"'")

View(df_tube)

#### recode oils names ####
df_tube$code_huiles[df_tube$huiles=="CoD003"] <- "CC"
df_tube$code_huiles[df_tube$huiles=="CoD008"] <- "CC/HS"
df_tube$code_huiles[df_tube$huiles=="CoD009"] <- "HS"

## calculating mortality rate 
df_tube$total <- df_tube$dead24h + df_tube$alive24h
df_tube$mort_rate <- df_tube$dead24h/(df_tube$dead24h + df_tube$alive24h)
df_tube$dose <- as.factor(df_tube$dose)
df_tube$Kd5_rate <- df_tube$Kd5 / df_tube$total
df_tube$Kd10_rate <- df_tube$Kd10 / df_tube$total
df_tube$Kd15_rate <- df_tube$Kd15 / df_tube$total
df_tube$Kd20_rate <- df_tube$Kd20 / df_tube$total
df_tube$Kd25_rate <- df_tube$Kd25 / df_tube$total
df_tube$Kd30_rate <- df_tube$Kd30 / df_tube$total
df_tube$Kd40_rate <- df_tube$Kd40 / df_tube$total
df_tube$Kd50_rate <- df_tube$Kd50 / df_tube$total
df_tube$Kd60_rate <- df_tube$Kd60 / df_tube$total

#### Data handling ### 
df_tube <- df_tube %>% tidyr::unite(code, code_huiles, dose, sep="_",remove = FALSE) 

df_tube$code_huiles <- as.factor(df_tube$code_huiles)
df_tube$code <- as.factor(df_tube$code)
df_tube$souches <- as.factor(df_tube$souches)

### Mortality rate model ###
# Normality test
shapiro.test(df_tube$mort_rate)
qqnorm(df_tube$mort_rate,ylim=c(0,1), main = "Inactive")
qqline(df_tube$mort_rate)
# data don't follow noromal distribution

## Plot density
density(df_tube$mort_rate)

mod_mort1 <- glm(cbind(dead24h,alive24h)~code_huiles+souches+concent, data=df_tube, 
                 family =binomial())
#deviance(mod_mort1)/df.residual(mod_mort1) # sans importance
summary(mod_mort1)
car::Anova(mod_mort1,type=3)

### 
data_HC <- summary(emmeans(mod_mort1, ~code_huiles+concent), type="response")
data_HS <- summary(emmeans(mod_mort1, ~souches+code_huiles), type="response")
data_HSC <- summary(emmeans(mod_mort1, ~souches+code_huiles+concent), type="response")

## Multiples comparaison ##
#Huiles - Concent
emm_hc <- emmeans(mod_mort1, specs = pairwise~code_huiles|concent, type="response", adjust = "none")
hc_contrast <- summary(emm_hc$contrasts, infer=TRUE)%>%as.data.frame() 
hc_emm <- summary(emm_hc$emmeans, infer=TRUE)%>%as.data.frame() 

#Huiles - souches
emm_hs <- emmeans(mod_mort1, specs = pairwise~code_huiles|souches, type="response", adjust = "none")
hs_contrast <- summary(emm_hs$contrasts, infer=TRUE)%>%as.data.frame() 
hs_emm <- summary(emm_hs$emmeans, infer=TRUE)%>%as.data.frame() 

# Huiles - souches - Concent

# to be continued 

#### data table for mortality rate graphic ####
mort_tab <- df_tube%>%
  dplyr::group_by(souches, code_huiles, dose)%>%
  dplyr::summarise(mort = mean(mort_rate)*100,sd=100*sd(mort_rate)*sqrt(3/4),
                   se=(100*sd(mort_rate)*sqrt(3/4))/sqrt(4))
mort_tab$dose <- as.factor(mort_tab$dose)

## vk-labo
p1 <- ggplot(data = subset(mort_tab, souches=='VK-labo'), aes(x=dose,y=mort))+
  geom_bar(aes(fill = dose),stat="identity")+
  geom_errorbar(aes(ymin = mort - se, ymax = mort + se),
                width=0.2,)+
  facet_wrap(~code_huiles)

p1+
  geom_hline(aes(yintercept = 90,linetype = "90%"),col="darkred")+
  theme_bw()+ 
  theme(panel.grid = element_blank())+
  scale_y_continuous(breaks = seq(0,100, by=20),limits = c(0,100))+
  guides(color = 'legend')+
  scale_fill_brewer(palette = "Reds")+
  scale_linetype_manual(name="WHO Susceptibility \nthreshold", values = "dashed")+
  labs(title="Mortality rate of An. gambiae VK", x="essential oil doses", 
       y = "Mortality rate (%)")

p1+
  geom_hline(aes(yintercept = 90,linetype = "90%"),col="darkred")+
  theme_bw()+ 
  theme(panel.grid = element_blank())+
  scale_y_continuous(breaks = seq(0,100, by=20),limits = c(0,100))+
  guides(color = 'legend')+
  scale_fill_brewer(palette = "Reds")+
  scale_linetype_manual(name="WHO Susceptibility \nthreshold", values = "dashed")+
  labs(x="essential oil doses", y = "Mortality rate (%)")


## An. gambiae kisumu
p2 <- ggplot(data = subset(mort_tab, souches=='Kisumu'), aes(x=dose,y=mort))+
  geom_bar(aes(fill = dose),stat="identity")+
  geom_errorbar(aes(ymin = mort - se, ymax = mort + se),
                width=0.2,)+
  facet_wrap(~code_huiles)

p2+
  geom_hline(aes(yintercept = 90,linetype = "90%"),col="darkred")+
  theme_bw()+ 
  theme(panel.grid = element_blank())+
  scale_y_continuous(breaks = seq(0,100, by=20),limits = c(0,100))+
  guides(color = 'legend')+
  scale_fill_brewer(palette = "Reds")+
  scale_linetype_manual(name="WHO Susceptibility \nthreshold", values = "dashed")+
  labs(title="Mortality rate of An. gambiae Kisumu", x="essential oil doses", 
       y = "Mortality rate (%)")

p2+
  geom_hline(aes(yintercept = 90,linetype = "90%"),col="darkred")+
  theme_bw()+ 
  theme(panel.grid = element_blank())+
  scale_y_continuous(breaks = seq(0,100, by=20),limits = c(0,100))+
  guides(color = 'legend')+
  scale_fill_brewer(palette = "Reds")+
  scale_linetype_manual(name="WHO Susceptibility \nthreshold", values = "dashed")+
  labs(x="essential oil doses", y = "Mortality rate (%)")

#### Model with drc _ An. gambiae kisumu #####
## dose response Model of essentials oils
drc_m_kis <- drm(dead24h/total~concent, code_huiles, data = subset(df_tube, souches=='Kisumu'), 
                 fct = LL.3u(), weights = total,
              type = 'binomial', pmodels = data.frame(code_huiles,1,code_huiles))
summary(drc_m_kis)
modelFit(drc_m_kis)
mselect(drc_m_kis)


## effective dose estimation 
kis_ed1 <- ED(drc_m_kis, c(50,95),interval = 'delta')%>%
  as.data.frame()
kis_ed2 <- ED(drc_m_kis, c(10,50,95), interval = 'fls')%>%
  as.data.frame()

## comparing effective dose 
EDcomp(drc_m_kis,c(50,50), interval = 'delta', od = T)
EDcomp(drc_m_kis,c(95,95), interval = 'delta')
#comparing parameters 
compParm(drc_m_kis,'b')
compParm(drc_m_kis,'e')

### plot the lethal dose graph
plot(drc_m_kis, ylim = c(0,1), bp=-1, lty = 'solid', lwd=2, pch = c(1,6,13), cex = 1, 
     broken = T, col = c('red', 'blue', 'green'),xlab = 'Doses', 
     ylab = 'Mortality rates of An. gambiae Kisumu', legendPos = c(0.05,0.8), 
     legend = T)
abline(h=0.95, col='darkred', lty='dashed')
abline(h=0.5, col='darkred', lty='dashed')

ggplot(subset(df_tube, souches=='Kisumu'), aes(x=concent, y=dead24h))+
  geom_point() +
  stat_smooth(method = "drm", method.args = list(fct = LL.3u()),
              aes(weight = total),se = TRUE)+
  facet_wrap(~code_huiles)


#### Model with drc _ An. gambiae vk #####
## dose response Model of essentials oils
drc_m_vk <- drm(dead24h/total~concent, code_huiles, data = subset(df_tube, souches=='VK-labo'), fct = LL.3u(), weights = total,
                 type = 'binomial', pmodels = data.frame(huiles,1,huiles))
summary(drc_m_vk)
modelFit(drc_m_vk)

## effective dose estimation 
vk_ed1 <- ED(drc_m_vk, c(50,95),interval = 'delta')%>%
  as.data.frame()
vk_ed2 <- ED(drc_m_vk, c(10,50,95), interval = 'fls')%>%
  as.data.frame()

## comparing effective dose 
EDcomp(drc_m_vk,c(50,50), interval = 'delta', od = T)
EDcomp(drc_m_vk,c(95,95), interval = 'delta')
#comparing parameters 
compParm(drc_m_vk,'b')
compParm(drc_m_vk,'e')

### plot the lethal dose graph
plot(drc_m_vk, ylim = c(0,1), bp=-1, lty = 'solid', lwd=2, pch = c(1,6,13), cex = 1, 
     broken = T, col = c('red', 'blue', 'green'),xlab = 'Doses', 
     ylab = 'Mortality rates of An. gambiae vk', legendPos = c(0.05,0.8), 
     legend = T)
abline(h=0.95, col='darkred', lty='dashed')
abline(h=0.5, col='darkred', lty='dashed')


#### Lethal dose ####
## An. gambiae kisumu ##
## souches == Ag kisumu & huiles == 'CoD003' 
data_kis1 <- subset(df_tube, souches=='Kisumu'& huiles== 'CoD003' )

kis_lc1 <- LC_probit((dead24h/total)~log10(concent),weights = total,p=c(50,95),
                         data = subset(df_tube, souches=='Kisumu'& huiles== 'CoD003'),
                         conf_level = 0.95)

## souches == Ag kisumu & huiles == 'CoD008' 
kis_lc2 <- LC_probit((dead24h/total)~log10(concent),weights = total,p=c(50,95),
                     data = subset(df_tube, souches=='Kisumu'& huiles== 'CoD008'),
                     conf_level = 0.95)

## souches == Ag kisumu & huiles == 'CoD009' 
kis_lc3 <- LC_probit((dead24h/total)~log10(concent),weights = total,p=c(50,95),
                     data = subset(df_tube, souches=='Kisumu'& huiles== 'CoD009'),
                     conf_level = 0.95)

# regression graph 1
p_reg1 <- ggplot(data = subset(df_tube, souches=='Kisumu'),
       aes(x=concent, y=dead24h/total))+
  geom_point(size=1, color='darkred')+
  geom_smooth(method = "glm",
              method.args = list(family = binomial(link = "probit")),
              aes(weight = total), colour = "#0000FF", se = TRUE)+
  facet_wrap(~code_huiles)

p_reg1+
  theme_bw()+ 
  theme(panel.grid = element_blank())+
  labs(title="regression graph1 - Ag Kisumu", x="essential oil doses (%)", 
       y = "Mortality rate")

## An. gambiae VK ##
## souches == Ag VK & huiles == 'CoD003' 
vk_lc1 <- LC_probit((dead24h/total)~log10(concent),weights = total,p=c(50,95),
                     data = subset(df_tube, souches=='VK-labo'& huiles== 'CoD003'),
                     conf_level = 0.95)


## souches == Ag kisumu & huiles == 'CoD008' 
vk_lc2 <- LC_probit((dead24h/total)~log10(concent),weights = total,p=c(50,95),
                     data = subset(df_tube, souches=='VK-labo'& huiles== 'CoD008'),
                     conf_level = 0.95)

## souches == Ag kisumu & huiles == 'CoD009' 
vk_lc3 <- LC_probit((dead24h/total)~log10(concent),weights = total,p=c(50,95),
                     data = subset(df_tube, souches=='VK-labo'& huiles== 'CoD009'),
                     conf_level = 0.95)

# regression graph 2
p_reg2 <- ggplot(data = subset(df_tube, souches=='VK-labo'),
                 aes(x=concent, y=dead24h/total))+
  geom_point(size=1, color='darkred')+
  geom_smooth(method = "glm",
              method.args = list(family = binomial(link = "probit")),
              aes(weight = total), colour = "#0000FF", se = TRUE)+
  facet_wrap(~code_huiles)

p_reg2+
  theme_bw()+ 
  theme(panel.grid = element_blank())+
  labs(title="regression graph1 - Ag vk-labo", x="essential oil doses (%)", 
       y = "Mortality rate")

#### Taux de KD #### 
#data handling for KD time graphic
kd_tab <- df_tube%>%group_by(souches,code,code_huiles)%>%
  summarise(mean_KD5=mean(Kd5_rate), sd_KD5=sd(Kd5_rate)*sqrt(3/4), 
            mean_KD10=mean(Kd10_rate), sd_KD10=sd(Kd10_rate)*sqrt(3/4),
            mean_KD15=mean(Kd15_rate), sd_KD15=sd(Kd15_rate)*sqrt(3/4),
            mean_KD20=mean(Kd20_rate), sd_KD20=sd(Kd20_rate)*sqrt(3/4),
            mean_KD25=mean(Kd25_rate), sd_KD25=sd(Kd25_rate)*sqrt(3/4),
            mean_KD30=mean(Kd30_rate), sd_KD30=sd(Kd30_rate)*sqrt(3/4),
            mean_KD40=mean(Kd40_rate), sd_KD40=sd(Kd40_rate)*sqrt(3/4),
            mean_KD50=mean(Kd50_rate), sd_KD50=sd(Kd50_rate)*sqrt(3/4),
            mean_KD60=mean(Kd60_rate), sd_KD60=sd(Kd60_rate)*sqrt(3/4))%>%
  as.data.frame()

## An. gambiae kisumu 
# select data to trsnform
names(kd_tab)
#View(kd_tab[c(2,seq(4,21,2))])
kdt_kis <- subset(kd_tab, souches=='Kisumu')[c(2,seq(4,21,2))]

# transform data
data_kd_kis <- setNames(data.frame(t(kdt_kis[,-1])), kdt_kis[,1])
data_kd_kis$kd_time <- c(5,10,15,20,25,30,40,50,60)
colnames(data_kd_kis)

## melt data and add some columns
data_kd1 <- melt(data_kd_kis, id = c("kd_time"))
dose <- c(rep(c(0.01, 0.1, 1.0), each = 9, times=3))
huile <- c(rep(c('CC', 'CC/HS', 'HS'), each = 27))

data_kd1$dose <- dose
data_kd1$huile <- huile

## KD time graphic - Ag kisumu 
p_kd1 <- ggplot(data = data_kd1, aes(x=kd_time,y=100*value))+
  geom_line(aes(colour=huile))+
  facet_wrap(~dose)

p_kd1+
  geom_hline(aes(yintercept = 50,linetype = '50%',col="darkred"))+
  geom_hline(aes(yintercept = 90,linetype = '90%',col="darkred"))+
  theme_bw()+
  labs(title = 'Ag - kisumu', x = "KD time",y="KD rate (%)",color="dose")+
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 8))+
  scale_linetype_manual(name= "KD threshold", values = c("dashed", "dashed"))

p_kd1+
  geom_hline(aes(yintercept = 50,linetype = '50%'),col="darkred")+
  geom_hline(aes(yintercept = 90,linetype = '90%'),col="darkred")+
  theme_bw()+
  scale_y_continuous(breaks = seq(0,100, by=10),limits = c(0,100))+
  scale_x_continuous(breaks = seq(0,60, by=10),limits = c(0,60))+
  labs(x = "KD time",y="KD rate (%)",color="dose")+
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 8))+
  scale_linetype_manual(name= "KD threshold", values = c("dashed", "dashed"))

## An. gambiae vk labo
# select data to trsnform
#names(kd_tab)
#View(kd_tab[c(2,seq(4,21,2))])
kdt_vk <- subset(kd_tab, souches=='VK-labo')[c(2,seq(4,21,2))]

# transform data
data_kd_vk <- setNames(data.frame(t(kdt_vk[,-1])), kdt_vk[,1])
data_kd_vk$kd_time <- c(5,10,15,20,25,30,40,50,60)
colnames(data_kd_vk)

## melt data and add some columns
data_kd2 <- melt(data_kd_vk, id = c("kd_time"))
dose <- c(rep(c(0.01, 0.1, 1.0), each = 9, times=3))
huile <- c(rep(c('CC', 'CC/HS', 'HS'), each = 27))

data_kd2$dose <- dose
data_kd2$huile <- huile

## KD time graphic - Ag kisumu 
p_kd2 <- ggplot(data = data_kd2, aes(x=kd_time,y=100*value))+
  geom_line(aes(colour=huile))+
  facet_wrap(~dose)

p_kd2+
  geom_hline(aes(yintercept = 50,linetype = '50%'),col="darkred")+
  geom_hline(aes(yintercept = 90,linetype = '90%'),col="darkred")+
  theme_bw()+
  scale_y_continuous(breaks = seq(0,100, by=20),limits = c(0,100))+
  scale_x_continuous(breaks = seq(0,60, by=10),limits = c(0,60))+
  labs(title = 'Ag - VK labo', x = "KD time",y="KD rate (%)",color="dose")+
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 8))+
  scale_linetype_manual(name= "KD threshold", values = c("dashed", "dashed"))

p_kd2+
  geom_hline(aes(yintercept = 50,linetype = '50%'),col="darkred")+
  geom_hline(aes(yintercept = 90,linetype = '90%'),col="darkred")+
  theme_bw()+
  scale_y_continuous(breaks = seq(0,100, by=10),limits = c(0,100))+
  scale_x_continuous(breaks = seq(0,60, by=10),limits = c(0,60))+
  labs(x = "KD time",y="KD rate (%)",color="dose")+
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 8))+
  scale_linetype_manual(name= "KD threshold", values = c("dashed", "dashed"))


#### end of analyses ####

