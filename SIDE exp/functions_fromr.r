all_subnum="-"
p.meanRT=function(d,tit){
  tit=paste(tit, all_subnum)
  d0=d#%>%filter(Block %in% filterBlock) 
  dase = d0 %>% 
    # filter(Block<=4)%>%
    group_by(FileCondi, Correctness, Old,Setsize, Stimkind) %>%
    dplyr::summarize(RTm= mean(RT),
                     se =sd(RT)/sqrt(n()))
  # print(dase)
  da1se=subset(dase,Correctness==1)#;da1se
  p=ggplot(data=da1se,aes(Setsize,RTm))+
    geom_errorbar(aes(ymin=RTm-se, ymax=RTm+se), width=.2,
                  position=position_dodge(.9),alpha=0.3) +
    geom_point(aes(color=as.factor(Old),shape=as.factor(Stimkind),group=Old),size=5)+
    geom_line(aes(color=as.factor(Old),#linetype=as.factor(FileCondi),
                  group=interaction(Old,Stimkind)))+
    scale_color_manual(name="Old-New",#breaks=c(1,0),
                       labels=c(`1`="Old", `2`="New"),
                       values=c("#F23005","#CC984D"))+
    scale_shape_discrete(
      name="ProbKind",labels = c(`1`="CM",`0`="AN",`3`="VM"))+
    # scale_linetype_discrete(
    #   name="Condition")+
    ggtitle(paste("mean CorrectRT - Set Size",tit))+
    theme_bw()+
    theme(text=element_text(size=16))+
    # ylim(700,1550)+
    facet_wrap(FileCondi~.,ncol=6)
  
  return(p)
}
p.medianRT=function(d,tit){
  cols <- c(`1`=1,`0`=2,`3`=3)
  
  tit=paste(' ',all_subnum," ",tit)
  dase = d %>% #filter(Block %in% BlockFilter) %>%
    group_by(FileCondi,Correctness,  Old,Setsize, Stimkind) %>%
    dplyr::summarize(RTm= median(RT),
                     se =sd(RT)/sqrt(n()))
  # print(dase)
  da1se=subset(dase,Correctness==1)#;da1se
  p=ggplot(data=da1se,aes(Setsize,RTm))+
    geom_errorbar(aes(ymin=RTm-se, ymax=RTm+se), width=.2,
                  position=position_dodge(.9),alpha=0.3) +
    geom_point(aes(color=as.factor(Old),shape=as.factor(Stimkind),group=Old),size=5)+
    geom_line(aes(color=as.factor(Old),#linetype=as.factor(FileCondi),
                  group=interaction(Old,Stimkind)))+
    scale_color_manual(name="Old-New",#breaks=c(1,0),
                       labels=c(`1`="Old", `2`="New"),
                       values=c("#F23005","#CC984D"))+
    scale_shape_discrete(
      name="ProbKind",labels = c(`1`="CM",`0`="AN",`3`="VM"))+
    # scale_linetype_discrete(
    #   name="Condition")+
    ggtitle(paste("Median CorrectRT - Set Size",tit))+
    theme_bw()+
    theme(text=element_text(size=16))+
    ylim(500,1100)+
    facet_wrap(FileCondi~.,ncol=6)
  # scale_colour_manual(name="Error Bars",values=cols, guide = guide_legend(shape = NULL,colour = NULL)) + 
  # scale_shape_manual(name="Bar",values=cols, guide="none") 
  
  return(p)
}

p.ER=function(d,tit,BlockFilter){
  tit=paste(' ',all_subnum," ",tit)
  d0= d %>% filter(Block %in% BlockFilter)
  dase=summarySEwithin(data=d0,measurevar = "Correctness",withinvars =
                         c("FileCondi","Old","Setsize","Stimkind"))
  # print(dase)
  da1se=dase
  p=ggplot(data=da1se,aes(Setsize,1-Correctness))+
    geom_errorbar(aes(ymin=(1-Correctness)-se, ymax=(1-Correctness)+se), width=.2,
                  position=position_dodge(.9),alpha=0.5) +
    geom_point(aes(color=as.factor(Old),shape=as.factor(Stimkind),group=Old),size=5)+
    geom_line(aes(color=as.factor(Old),#linetype=as.factor(FileCondi),
                  group=interaction(Old,Stimkind)))+
    scale_color_manual(name="Old-New",#breaks=c(1,0),
                       labels=c(`1`="Old", `2`="New"),
                       values=c("#F23005","#CC984D"))+
    scale_shape_discrete(
      name="ProbKind",labels = c(`1`="CM",`0`="AN",`3`="VM"))+
    # scale_linetype_discrete(
    #   name="Condition")+
    ggtitle(paste("Error Rate - Set Size",tit))+
    theme_bw()+
    theme(text=element_text(size=16))+
    ylim(0,0.25)+
    facet_wrap(FileCondi~.,ncol=6)+
    scale_y_continuous("Probability of Error")
  return(p)
}