# import of individual & population table (reduced from STUCTURE results):
# -> import format struct: "individual_index;Label;cluster1;cluster2;cluster3;cluster4
# -> pop format: pop_index;cluster1;cluster2;cluster3;cluster4;n_cluster

setwd('~/Desktop/')

struct <- read.table("k4.csv", sep = ' ')
pop <- read.table("pop.txt", sep = ';')

library(ggplot2)
library(reshape2)

# reordering of tabels ()

strm <- melt(struct[,c('V1','V3','V4','V5','V6')],id.vars = 1)
popm <- melt(pop[,c('V1','V2','V3','V4','V5')],id.vars = 1)

# k = 2
strm <- melt(struct[,c('V1','V3','V4')],id.vars = 1)
vals =  c('V2','V3','V4')
cols = c(rgb(0,0,.4),rgb(0.9,0,0))

# k = 3
strm <- melt(struct[,c('V1','V3','V4','V5')],id.vars = 1)
vals =  c('V2','V3','V4')
cols = c(rgb(0,.6,0),rgb(0.9,0,0),rgb(0,0,.4))

# k = 4
strm <- melt(struct[,c('V1','V3','V4','V5','V6')],id.vars = 1)
vals =  c('V2','V3','V4','V5')
cols = c(rgb(0.9,0,0),rgb(0,.6,0),rgb(.9,.6,0),rgb(0,0,.4))

# get groupsize from pop-tabel (pop$V6)
grpSize <- c(32,32,32,50,30,30,30,30,48,30,40,27,35)
grpSize_half <- grpSize/2
# calculate starting index for each group
grpInd_end <- cumsum(grpSize)
grpInd <- round(grpInd_end - grpSize_half)

# define group labels
sites <-  c('WFE1','WFE2','BR','MG','LS','RB','BS','SP','IR','QS','SM','CL','LP')
cols_label = c(rep(rgb(.9,.6,0),3),rep(rgb(0.9,0,0),3),rgb(.4,0,.4),rep(rgb(0,0,.4),4),rep(rgb(0,.6,0),2))
# ------ structure plot -----------------

pdf(file="~/Desktop/Structurek2.pdf",width=15,height=3.5,pointsize=12, bg = 'transparent')
# (yvalues normalized to 1 (minimal variance in pure sum of clusters))
ggplot(strm,aes(x= V1,y = value/(rep(tapply(strm$value,strm$V1,sum),2)), fill = variable))+
  # the actual STRUCTURE plot
  geom_bar(stat = 'identity', width=1)+
  geom_vline(xintercept= rep(grpInd_end,each=2), col = rgb(0,0,0,.9),lwd = .8)+
  # add frame
  #geom_polygon(aes(x= c(0.5,446.5,446.5,0.5), y = c(0,0,1,1)),fill=NA,colour ='black')+
  # define manual colors
  scale_fill_manual(name = 'K = 4',
                    values = cols,
                    breaks = vals)+
  # custom axes
  scale_y_continuous(expand = c(0,.001),breaks = (0:4)/4)+
  scale_x_continuous(expand = c(0.01,0.01),
                   breaks = grpInd,
                   labels = sites)+
  # general layout issues
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA,colour = NA),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text.y=element_text(size = 20),
        axis.text.x=element_text(angle = 90,vjust = 0.5,size = 22,color =cols_label))+
  guides(fill=FALSE)
dev.off()


# ----------------------------------------
# (older versions)

# population overview
ggplot(popm,aes(x= V1,y = value/(rep(tapply(popm$value,popm$V1,sum),4)), fill = variable))+
  geom_bar(stat = 'identity', width=1)+
  scale_fill_manual(values=c(rgb(.9,.6,0),rgb(0.9,0,0),rgb(0,0,.4),rgb(0,.6,0)),
                    breaks=c('V2','V3','V4','V5'))


# ------ basic structure -----------------
ggplot(strm,aes(x= V1,y = value/(rep(tapply(strm$value,strm$V1,sum),4)), fill = variable))+
  geom_bar(stat = 'identity', width=1)+
  geom_vline(xintercept= rep(grpInd,each=2), col = rgb(1,1,1,.8),lwd = .8)+
  geom_polygon(aes(x= c(0,446,446,0), y = c(0,0,1,1)),fill=NA,colour ='black')+
  scale_fill_manual(name = 'K = 4',
                    values=c(rgb(.9,.6,0),rgb(0.9,0,0),rgb(0,0,.4),rgb(0,.6,0)),
                    breaks=c('V2','V3','V4','V5'),
                    labels = c(expression('k'[1]),
                               expression('k'[2]),
                               expression('k'[3]),
                               expression('k'[4])))+
  scale_y_continuous(expand = c(0,.01),breaks = (0:4)/4)+
  scale_x_discrete(expand = c(0.01,0.01),breaks = 100*(0:4))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(), 
        axis.line = element_blank(),axis.title = element_blank(),
        axis.text=element_text(size=18,face = "bold"),
        legend.text=element_text(size=18,face = "bold"),legend.title=element_text(size=18))
