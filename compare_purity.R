library ( ggplot2 )
library ( patchwork )
library ( dplyr )
library ( openxlsx )

rna.key  <- read.xlsx( config$spcg.key , sheet="masterkey.stable.txt")
rna.key[ is.na(rna.key) ] <- ""
rna.key = rna.key[rna.key$sample.source == "patient" , ]
rna.key=rna.key[rna.key$RNAseq.id != "None", ]
rna.key=rna.key[rna.key$patient.id != "None", ]


# open est 

est = readRDS(  paste0( "/projects/spcg/rna/SUMMARY/ESTIMATE/EST.rds"))



# get wgs purity. 

wgs_purity <- read.xlsx( "https://www.dropbox.com/s/ugw0iu0g5eih02b/SPCG-purity.xlsx?dl=1", sheet="Sheet1" )

cmp_key   =  rna.key[ rna.key$WGS.id %in% intersect ( wgs_purity$Sample, rna.key$WGS.id), ] 

est = est[ row.names ( est )  %in% cmp_key$RNAseq.id,  ]


cmp_key = cmp_key [ , c("RNAseq.id", "WGS.id",  "cancer.subtype", "treatment", "sample.classification"  )]
cmp_key = merge ( cmp_key , est, by.x="RNAseq.id", by.y="row.names" )


cmp_key = merge ( cmp_key, wgs_purity[ , c("Sample", "Purity")], by.x="WGS.id", by.y="Sample")

cmp_key$delta = abs ( cmp_key$purity - cmp_key$Purity) 
cmp_key = cmp_key[ complete.cases(cmp_key),   ]

colnames(cmp_key)[ which ( colnames(cmp_key)  %in% c("purity", "Purity") ) ] = c("rna_purity", "dna_purity")


cor ( cmp_key$rna_purity, cmp_key$dna_purity )


ggplot(cmp_key, aes(x=dna_purity, y=rna_purity)) +
  geom_point(size=10, color="steelblue", alpha=.25) + 
  geom_smooth(method=lm
              , se=TRUE
              , fullrange=TRUE
              , color = "black"
  )+
  theme_classic() +  stat_cor(method = "spearman", size=12) +
  theme(legend.position="none", legend.title=element_blank(), legend.key = element_blank(),
        
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(angle = 0, size=35, hjust=.92,vjust=0.1),
        axis.title.x = element_text(size=30),
        plot.title = element_text(size=20), 
        axis.title.y     = element_text(size=30), 
        legend.text      =element_text(size=12)
  ) 

# finally wgs_purity plot 
wkey  <- read.xlsx( config$spcg.key , sheet="masterkey.stable.txt")
wkey[ is.na(wkey) ] <- ""

wgs_purity <- read.xlsx( "https://www.dropbox.com/s/ugw0iu0g5eih02b/SPCG-purity.xlsx?dl=1", sheet="Sheet1" )

wkey <- merge ( wkey, wgs_purity[ , c("Sample", "Purity")] , by.x="WGS.id", by.y="Sample") 
wkey = wkey[ complete.cases(wkey),   ]


ggplot(wkey, aes(x=Purity)) + 
  geom_histogram(aes(y=(..count..)/sum(..count..)), colour='black' , fill="#6e9be5", alpha=.2 , binwidth=.1 )+
 # geom_density(alpha=0, fill="grey") +
  theme(legend.position="none",  legend.key = element_blank(),
        # element_blank()
        axis.text.y = element_text(size= 25 ),
        axis.text.x =  element_text(size= 25 ),
        axis.title.x = element_text(size=25),
        axis.title.y     = element_text(size=25), 
        legend.text      =element_text(size=25),
        legend.title = element_text(size=25),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5)
        # hjust centers the title
  ) + theme(panel.grid.major = element_blank()
            , panel.grid.minor = element_blank()
            ,panel.background = element_blank()
            , axis.line = element_line(colour = "black") # plot border
  ) + ylab ( "percent")



ggplot(wkey, aes(x=Purity)) + 
  geom_histogram(aes(y=(..count..)/sum(..count..), fill = treatment ), colour='black' , alpha=.2 , binwidth=.1 )+
  # geom_density(alpha=0, fill="grey") +
  theme(legend.position="bottom",  legend.key = element_blank(),
        # element_blank()
        axis.text.y = element_text(size= 25 ),
        axis.text.x =  element_text(size= 25 ),
        axis.title.x = element_text(size=25),
        axis.title.y     = element_text(size=25), 
        legend.text      =element_text(size=25),
        legend.title = element_text(size=25),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5)
        # hjust centers the title
  ) + theme(panel.grid.major = element_blank()
            , panel.grid.minor = element_blank()
            ,panel.background = element_blank()
            , axis.line = element_line(colour = "black") # plot border
  ) + ylab ( "percent")



ggplot(wkey, aes(x=Purity)) + 
  geom_histogram(aes(y=(..count..)/sum(..count..), fill = sample.classification ), colour='black' , alpha=.2 , binwidth=.1 )+
  # geom_density(alpha=0, fill="grey") +
  theme(legend.position="bottom",  legend.key = element_blank(),
        # element_blank()
        axis.text.y = element_text(size= 25 ),
        axis.text.x =  element_text(size= 25 ),
        axis.title.x = element_text(size=25),
        axis.title.y     = element_text(size=25), 
        legend.text      =element_text(size=25),
        legend.title = element_text(size=25),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5)
        # hjust centers the title
  ) + theme(panel.grid.major = element_blank()
            , panel.grid.minor = element_blank()
            ,panel.background = element_blank()
            , axis.line = element_line(colour = "black") # plot border
  ) + ylab ( "percent")



ggplot(wkey, aes(x=Purity)) + 
  geom_histogram(aes(y=(..count..)/sum(..count..), fill = cancer.subtype ), colour='black' , alpha=.2 , binwidth=.1 )+
  # geom_density(alpha=0, fill="grey") +
  theme(legend.position="bottom",  legend.key = element_blank(),
        # element_blank()
        axis.text.y = element_text(size= 25 ),
        axis.text.x =  element_text(size= 25 ),
        axis.title.x = element_text(size=25),
        axis.title.y     = element_text(size=25), 
        legend.text      =element_text(size=25),
        legend.title = element_text(size=25),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5)
        # hjust centers the title
  ) + theme(panel.grid.major = element_blank()
            , panel.grid.minor = element_blank()
            ,panel.background = element_blank()
            , axis.line = element_line(colour = "black") # plot border
  ) + ylab ( "percent")



ggplot(wkey, aes(x=Purity)) + 
  geom_histogram(aes(y=(..count..)/sum(..count..), fill = cancer.type ), colour='black' , alpha=.2 , binwidth=.1 )+
  # geom_density(alpha=0, fill="grey") +
  theme(legend.position="bottom",  legend.key = element_blank(),
        # element_blank()
        axis.text.y = element_text(size= 25 ),
        axis.text.x =  element_text(size= 25 ),
        axis.title.x = element_text(size=25),
        axis.title.y     = element_text(size=25), 
        legend.text      =element_text(size=25),
        legend.title = element_text(size=25),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5)
        # hjust centers the title
  ) + theme(panel.grid.major = element_blank()
            , panel.grid.minor = element_blank()
            ,panel.background = element_blank()
            , axis.line = element_line(colour = "black") # plot border
  ) + ylab ( "percent")


