library(ggplot2)
output_dir<- "/sdata/images/projects/GENSCOT/1/andrew/painMDD/reports/"

graph_pain_painPGRS<-read.table("/sdata/images/projects/GENSCOT/1/andrew/painMDD/results/pgrs-ukb/summary_PAINpgrs_into_PAIN.txt", header=T, sep='\t')
graph_pain_MDDPGRS<-read.table("/sdata/images/projects/GENSCOT/1/andrew/painMDD/results/pgrs-ukb/summary_MDDpgrs_into_PAIN.txt", header=T, sep='\t')



ggfigplot1 <- ggplot(data=graph_pain_painPGRS, 
              aes(x=as.factor(pT), y=Beta, fill=Beta )) + 
              geom_bar(stat = "identity", width=0.75) + 
              facet_wrap(~Cohort) +
              theme(axis.text=element_text(size=12), 
                    axis.title=element_text(size=16),
                    strip.text.x=element_text(size=16)) +
              xlab("Polygenic threshold") +
              ylab("Standardised Beta")

jpeg(filename = paste0(output_dir,"PGRS_plot_PAINpgrsintoPAIN.jpeg"),
     width = 2400, height = 1200, units = "px", 
     quality = 150, res=120,
     bg = "white")
ggfigplot1
dev.off()

ggfigplot2 <- ggplot(data=graph_pain_MDDPGRS, 
              aes(x=as.factor(pT), y=Beta, fill=Beta )) + 
              geom_bar(stat = "identity", width=0.75) + 
              facet_wrap(~Cohort) +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16),
        axis.text.x=element_text(size=14, colour="grey20"),
        strip.text.x=element_text(size=16)) +
              xlab("Polygenic threshold") +
              ylab("Standardised Beta")

jpeg(filename = paste0(output_dir,"PGRS_plot_MDDpgrsintoPAIN.jpeg"),
     width = 2400, height = 1200, units = "px", 
     quality = 150, res=120,
     bg = "white")
ggfigplot2
dev.off()