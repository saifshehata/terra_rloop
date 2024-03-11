# plot profile outmatrixdata file

profileData<-fread("~/MondalLab/IGVsubtelomereProject/deepToolsMatrices/scaleRegions.csh.mettl3kd.tab", sep="\t", header = F, skip=1, col.names = c("sample","genes",as.character(1:900)))

profileData %>% gather(.,key="bin",value="value", -sample,-genes) -> df

df %>% mutate(unscaled=if_else(sample=="Control-sh m6A",value*1.05469301340861,if_else(sample=="METTL3-KD input",value*1.3987881508079,if_else(sample=="METTL3-KD m6A",value*1.42780862879112,if_else(sample=="Control-sh input",value,value))))) -> df.unscaled

m.input<-df.unscaled %>% filter(sample=="METTL3-KD input")
m.m6a<-df.unscaled %>% filter(sample=="METTL3-KD m6A") 

c.input<-df.unscaled %>% filter(sample=="Control-sh input") 
c.m6a<-df.unscaled %>% filter(sample=="Control-sh m6A")

m.logratio<-left_join(m.input,m.m6a,by=c("bin"="bin")) %>% mutate(log2ratio=log2(value.y/value.x)) %>% mutate(sample="METTL3-KD") %>%  dplyr::select(sample,bin,log2ratio)

c.logratio<-left_join(c.input,c.m6a,by=c("bin"="bin")) %>% mutate(log2ratio=log2(value.y/value.x)) %>% mutate(sample="Control-sh") %>%  dplyr::select(sample,bin,log2ratio)

proflog2<-rbind(c.logratio,m.logratio)




proflog2 %>% mutate(bin=as.numeric(bin)) %>% ggplot(.,aes(bin,log2ratio,group=sample,color=sample)) + geom_line() + theme_cowplot() + scale_x_continuous(breaks = c(0,300,600,900), labels = c("-3.0kb","TSS","TES","+3.0kb")) + ylim(-0.5,1.5) + scale_color_manual(values=c("steelblue","salmon")) + ylab(expression(log[2](m6A/input))) + xlab("") + labs(color="")