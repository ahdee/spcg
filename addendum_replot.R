


counts = paste0 ( resource , "public/treehouse/STATIC.TPM.V11.rds" )
counts = readRDS(counts)
 # as.numeric ( tpm [ c( "TP53", "GAPDH") , "AA3" ] ) 
 # as.numeric ( counts$treehouse [ c( "TP53", "GAPDH") , "AA3" ] ) 
 # dim ( counts$treehouse )
 # dim ( tpm )
# remove all ASC 
gtex_key = counts$gtex.key
gtex = counts$gtex
v11_key  = counts$tree.key
counts = counts$treehouse


## figuring out GTEx 
gtype = unique ( sort ( gtex_key$type2))


gHeme = gtype [ grepl ( "lymp|blood", gtype, ignore.case = T)]
gCNS = gtype [ grepl ( "brain|pitu|nerv", gtype, ignore.case = T)]
gSolid = gtype [ ! gtype %in% c( gHeme, gCNS)]

# this will match cancer.type of the key 
# > unique ( rna.key$cancer.type)
# "Leukemia" "Solid"    "CNS"     
gtex_id = list()

gtex_id [["CNS"]] = gtex_key[ gtex_key$type2 %in% gCNS,    ]$sample
gtex_id [["Solid"]] = gtex_key[ gtex_key$type2 %in% gSolid,    ]$sample
gtex_id [["Leukemia"]] = gtex_key[ gtex_key$type2 %in% gHeme,    ]$sample



by.pid.newtuk <- function (ids,test.this, label, tpm.df, tpm_id, size=12, size2 = 12, sizept = 1, sizeshaperange=22, low.border = 0, 
                           sizeshape = 5, sizexy = 22, ap=.5, cpp = "#afc1de" ,xl = " log2 ( TPM + 1) " , gtitle='outlier' ){
  
  if (is.na(cpp)){
    cpp = getPalette( length ( unique (test.this)))
  }
  
  if ( !any ( names ( tpm.df) %in% label ) ) {
     tpm.df = merge ( tpm_id , tpm.df , by="row.names" )
     row.names ( tpm.df ) = tpm.df$Row.names
     tpm.df$Row.names = NULL 
  }
  
  g1 = melt ( t (  log2 (  tpm.df[ row.names (tpm.df) %in% c(test.this) ,  ids , drop=F] + 1)        )  ) 
  
  colnames ( g1) = c("tube","gene","value")
  
  # get tukey boundary 
  t1 = setDT(g1)[, get.skew.out (value, tube=unique ( g1$tube),  pidthis =label), by = gene]
  
  
  g1 = merge ( g1, t1, by="gene")
  g1$gene = as.character(g1$gene)
  # after merging value.x is the samples value while value.y is the sample of interest's value
  # make it easier to understand the names 
  g1$gene = gsub ( "CD276", "CD276 (B7-H3)", g1$gene)
  g1$gene = gsub ( "CD274", "PDL1 (B7-H1; CD274)", g1$gene)
  g1$gene = gsub ( "VTCN1", "VTCN1 (B7-H4)", g1$gene)
  g1$gene = gsub ("PDCD1LG2", "PDL2 (PDCD1LG2)", g1$gene ) 
  #  
  # 
  gene = gsub ( "CD276", "CD276 (B7-H3)", g1$gene )
  gene = gsub ( "CD274", "PDL1 (B7-H1; CD274)", g1$gene )
  gene = gsub ( "VTCN1", "VTCN1 (B7-H4)", g1$gene )
  gene = gsub ("PDCD1LG2", "PDL2 (PDCD1LG2)", g1$gene  ) 
  
  g1$gene = factor ( g1$gene, levels = rev( unique ( g1$gene) ) )
  
  
  
  
  
  g3box <- ggplot(g1, aes(y=value.x, x=gene))+ #geom_boxplot(outlier.shape = NA)+ 
    geom_point( aes(x=gene, y=max), shape="|", size=size2, col="#e6004c") +
    geom_point( aes(x=gene, y=min), shape="|", size=size2, col="#b9d159") +
    geom_point(data= g1[g1$tube %in% label,] , aes(x=gene, y=value.x), colour="black", shape=8, size=size ) +
    geom_violin(alpha = 0.5) +
    geom_point(aes(fill = gene), size = sizeshape, shape = 21, position = position_jitterdodge(), alpha = 0.15 )+
    #geom_quasirandom(aes(fill = gene, color=gene ), size = sizeshape, shape = 21, alpha = ap) +
    theme_bw() +
    ylab(xl) +
    xlab("") +
    theme(legend.position="none", legend.title=element_blank(), legend.key = element_blank(),
          
          axis.text.y = element_text(size= sizexy ),
          axis.text.x = element_text(angle = 90, size= sizexy ),
          axis.title.x = element_text(size=22),
          
          axis.title.y     = element_text(size=22), 
          legend.text      =element_text(size=12)
    ) + coord_flip()  +
    
    scale_fill_manual ( values = cpp ) +
    scale_color_manual ( values = rep("grey", length(unique(g1$gene))) ) 
  
  g3box$layers = rev( g3box$layers)
  
  # joy plot 
  jplot = ggplot(g1, aes(x=value.x, y=gene, fill = gene)) +
    geom_density_ridges(
      scale = 0.9,
      aes( point_shape = 21),
      alpha = .52, point_alpha = .2, jittered_points = TRUE
    ) +scale_fill_manual(values = c(cpp)) +
    geom_point(data= g1[g1$tube %in% label,] , aes(x=value.x, y=gene), colour="black", shape=8, size=size ) +
    geom_point( aes(x=max, y=gene), shape="|", size=size2, col="#e6004c") + 
    geom_point( aes(y=gene, x=min), shape="|", size=size2, col="#b9d159") +
    xlab("log2 ( TPM + 1) ") + theme_bw() +
    theme(legend.position="none", legend.title=element_blank(), legend.key = element_blank(),
          
          axis.text.y = element_text(size= sizexy ),
          axis.text.x = element_text(angle = 90, size= sizexy ),
          axis.title.x = element_text(size=22),
          plot.title = element_text(size = 22, face = "bold"), 
          axis.title.y     = element_text(size=22), 
          legend.text      =element_text(size=12)
    ) + ggtitle ( paste ( label, gtitle) )
  
  
  if ( low.border == 1){
    #jplot = jplot + geom_point( aes(x=min, y=gene), shape="|", size=size2, col="#b9d159") 
    # already added on top 
  }
  # + facet_grid (gene~., scales="free")
  
  
  
  ## no down outlier 
  
  g3box.no <- ggplot(g1, aes(y=value.x, x=gene))+ #geom_boxplot(outlier.shape = NA)+ 
    geom_point( aes(x=gene, y=max), shape="|", size=size2, col="#e6004c") +
    #geom_point( aes(x=gene, y=min), shape="|", size=size2, col="#b9d159") +
    geom_point(data= g1[g1$tube %in% label,] , aes(x=gene, y=value.x), colour="black", shape=8, size=size ) +
    geom_violin(alpha = 0.5) +
    geom_point(aes(fill = gene), size = sizeshape, shape = 21, position = position_jitterdodge(), alpha = 0.15 )+
    #geom_quasirandom(aes(fill = gene, color=gene ), size = sizeshapplote, shape = 21, alpha = ap) +
    theme_bw() +
    ylab(xl) +
    xlab("") +
    theme(legend.position="none", legend.title=element_blank(), legend.key = element_blank(),
          
          axis.text.y = element_text(size= sizexy ),
          axis.text.x = element_text(angle = 90, size= sizexy ),
          axis.title.x = element_text(size=22),
          
          axis.title.y     = element_text(size=22), 
          legend.text      =element_text(size=12)
    ) + coord_flip()  +
    
    scale_fill_manual ( values = cpp ) +
    scale_color_manual ( values = rep("grey", length(unique(g1$gene))) ) 
  
  g3box.no$layers = rev( g3box.no$layers)
  
  
  
  
  min1 = min ( g1$value.x, g1$value.y, g1$max) - .02
  max1 = max ( g1$value.x, g1$value.y, g1$min) + .02
  
  g4box <- ggplot(g1, aes(y=value.x, x=gene))+
    #geom_point(aes(fill = gene), size = sizeshape, shape = 21, position = position_jitterdodge(), alpha = ap )+
    geom_quasirandom(aes(fill = gene, color=gene ), size = sizeshape, shape = 21, alpha = ap) +
    theme_bw() +
    ylab(xl) +
    xlab("") +
    theme(legend.position="none", legend.title=element_blank(), legend.key = element_blank(),
          
          axis.text.y = element_text(size= sizexy ),
          axis.text.x = element_text(angle = 0, size= sizexy ),
          axis.title.x = element_text(size=22),
          
          axis.title.y     = element_text(size=22), 
          legend.text      =element_text(size=12)
    ) + coord_flip() +
    geom_point( aes(x=gene, y=max), shape="|", size=size2, col="#e6004c") +
    geom_point( aes(x=gene, y=min), shape="|", size=size2, col="#b9d159") +
    scale_fill_manual ( values = cpp ) +
    scale_color_manual ( values = rep("grey", length(unique(g1$gene))) ) +
    geom_point(data= g1[g1$tube %in% label,] , aes(x=gene, y=value.x), colour="black", shape=8, size=size )
  
  
  
  ypos = max ( density(g1$value.x)$y ) + 2   
  gden =  ggplot(g1) + 
    geom_point(aes(value.x, ypos , fill = gene, size = value.x), 
               position = position_jitterdodge(jitter.height = 2), 
               shape = 21, alpha = 0.14) + scale_size_continuous(range = c(1, sizeshaperange)) +
    theme_bw()+ #geom_histogram(aes(x = value.x ), binwidth = 1) 
    geom_density(aes(x = value.x,  fill = gene ) , alpha = 0.2) +
    theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank(),
          
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 0, size= sizexy ),
          axis.title.x = element_text(size=22),
          
          axis.title.y     = element_text(size=22), 
          legend.text      =element_text(size=12)
    ) + ylab("") + xlab( "log2 ( TPM + 1) ") +
    scale_fill_manual ( values = cpp ) + geom_point(data= g1[g1$tube %in% label,] , aes(x=value.x, ypos), colour="black", shape=8, size=size ) +
    scale_color_manual ( values = rep("grey", length(unique(g1$gene))) ) +
    annotate("rect", xmin = g1[1:1,]$max, xmax = max(g1$value.x)+.05, ymin = 0, ymax = ypos+2.1, # jitter height set to 2 so just add a bit
             alpha = .2)  
  
  
  return ( list (violin=g3box, bee=g4box, dense=gden, basic = g3box.no, jplot=jplot)  )
  
  
  
  
}




## calculates oultier boarders based on how skew the data is 
## needs min.reads, that is it will only calculate with cohort samples with minimum amount of tpm. 
## default is .05 which is very low but from looking at the data I think .01'inish is just noice 
## Note default now set to 0 which means that it will no longer consider min.reads since I think the the skewness should adjust for that. 
## when to use tukey is < 1 by default tuk1 = 1
## when to use the less stringeng < 3.7 tuk2 = 3.7 
## after that use the algo by hubert & vandervieren, 2008
# https://cran.r-project.org/web/packages/univOutl/univOutl.pdf
# https://www.notion.so/Building-outlier-analysis-76d0e5c99c76445b8f5bd5ae6c71a53f#a8ccda5ba6354d908c55ac852b281705
#### IMPORTANT, THE FIRST PID IS THE ONE ITS TESTING AGAINST! 
get.skew.out <- function ( gg, tube, min.sd=.05, tuk1=1.6, tuk2=3.7,pidthis )  {
  
  wherepid = which ( tube == pidthis)
  test = gg[ wherepid ]    # get value for the sample u are testing. 
  medg = median( gg )
  
  # get zscore 
  zs = scale ( as.numeric ( gg) )
  zs = zs[  wherepid ]
  if ( is.na ( zs ) ) {zs=0}
  # not going to rely on min.reads anymore! 
  #if ( test <= min.reads ){
  # return ( list (  skew = 1000, sk.class = "NA" , min= 0 , max = 0, value= 0, class= "NONE",diff=0, median=medg, logfc=0) )
  #}
  
  # instead will use CV however CV is sensitive when most values are zero 
  # thus if that is the case then it defaults to sd 
  # min.sd is set to .05 
  
  # instead of min.reads we will do variance instead 
  # if the mean is close to zero that is < .05 then we will use variance else default to cv 
  if ( mean ( gg) < .5){
    cv = sd ( gg )
  }else{
    cv = sd(gg)/ mean ( abs ( gg)  )
  }
  
  if ( is.na ( cv ) | cv < min.sd ){
    return ( list (  skew = 1000, sk.class = "NA" , min= 0 , max = 0, value= 0, class= "NONE",diff=0, median=medg, logfc=0, cv=cv,zscore=zs) )
  }
  
  skv  = skewness(gg)
  logfc = test - mean ( gg[2:length(gg)] ) 
  
  sk.class = "tukey"
  if (skv < tuk1 ){
    box.tuk = suppressMessages( boxB( gg, k=1.5, method='resistant', weights=NULL, id=NULL,
                                      exclude=NA, logt=FALSE) )
  }else if( skv < tuk2 ){
    
    box.tuk = suppressMessages(  boxB( gg , k=1.5, method='asymmetric', weights=NULL, id=NULL,
                                       exclude=NA, logt=FALSE) )
    sk.class = "tukey.mod"
  }else{
    
    box.tuk = suppressMessages(  boxB( gg , k=1.5, method="adjbox", weights=NULL, id=NULL,
                                       exclude=NA, logt=FALSE) )
    sk.class = "Hubert"
  }
  
  
  box.tuk = suppressMessages( boxB( gg, k=1.5, method='resistant', weights=NULL, id=NULL,
                                    exclude=NA, logt=FALSE) )
  
  
  
  exp = "NONE"
  diff = abs( medg - gg[wherepid])
  if ( gg[wherepid] >  box.tuk$fences[2] ){
    exp = "OUTLIER_HIGH"
    
  }else if ( gg[wherepid] < box.tuk$fences[1] ){
    exp = "OUTLIER_LOW"
    
  }
  
  return ( list (  skew = skv, sk.class = sk.class , min= box.tuk$fences[1]  , max = box.tuk$fences[2], value= gg[wherepid], class= exp, diff=diff, median=medg
                   , logfc=logfc
                   , cv=cv, zscore=zs ) )
  
  
}

library ( univOutl)

get.info <- function ( tube, rk ){
  rk = rk [ rk$RNAseq.id == tube, ]
  mdirs = "/projects/spcg/rna/DATA/"
  
  mdirs = paste0( "/projects/spcg/rna/DATA/", rk$patient.institution, "_", rk$patient.id, "/", rk$RNAseq.id, "/summaries/", 
                  rk$patient.id,"_",rk$RNA.parent,"_",rk$RNAseq.id,"_rna___", rk$RNAseq.id, "_outlier.rds" )
  info = readRDS ( mdirs )
  return ( info )
}




# info[["txn2.best"]], 
# or info[[tree.best]]$neighbor$Row.names


test.id = "BY1"
gene.id = "FGFR1"
type_analysis = "tree.xcell" # "tree.xcell" # txn2.best

if ( type_analysis == "tree.xcell"){
  n_near = colnames ( info[[type_analysis]]  )
}else {
  n_near =  info[[type_analysis]]  
}

id_cancertype = rna.key[ rna.key$RNAseq.id == test.id, ]$cancer.type
info = get.info ( test.id , rna.key )
tg = by.pid.newtuk( ids=as.character ( unique ( c ( test.id , n_near  )  ) )
                    ,  test.this=unique ( c ( gene.id ) ) , label=test.id , tpm.df= counts, tpm_id = tpm [ , test.id, drop=F ],
                    size=8, size2 = 15, sizeshape = 3.5, sizexy = 16, ap=.15 , low.border=1, cpp = '#91badb', gtitle="Disease") 

# not goint to use matched: info$gtex.best$neighbor$Row.names
# use solid CNA or heme instead. 
gg = by.pid.newtuk( ids= as.character ( unique ( c ( test.id, gtex_id[[ id_cancertype ]] )  ) ) 
                    ,  test.this=gene.id , label=test.id , tpm.df= gtex, tpm_id = tpm [ , test.id, drop=F ],
                    size=8, size2 = 15, sizeshape = 3.5, sizexy = 16, ap=.15 , low.border=1, cpp = '#88d134',  gtitle= "Normal")

pan = by.pid.newtuk( ids=as.character ( unique ( c ( test.id , colnames(counts) )  ) )
                    ,  test.this=unique ( c ( gene.id ) ) , label=test.id , tpm.df= counts, tpm_id = tpm [ , test.id, drop=F ],
                    size=8, size2 = 15, sizeshape = 3.5, sizexy = 16, ap=.15 , low.border=1, cpp = '#8e8c8f', gtitle= "Pan") 

tg$jplot / gg$jplot / pan$jplot

