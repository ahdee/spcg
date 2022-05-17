

## setup
spcg_fig = "./spcg_fig/"
dir.create( spcg_fig)

# ongoing excel sheet 
paper_table

addWorksheet(paper_table, 'QC')
#writeData(wb, 'QC' , df  , rowNames=T  )
# writeData(paper_table, 'QC', )
writeData(paper_table
          ,'QC'
          ,rna.key
          , rowNames=F  )



################ heatmap main 1 figure

pdf(paste0(spcg_fig, "/heatmap_clusters.pdf"), width=16, height=14.5)
draw(figs_out[['heatmap1']],annotation_legend_list = figs_out[['heatmap1_leg']], merge_legend = TRUE   )
dev.off()


#################  heatmap main1 pathays between each clusters
pdf(paste0(spcg_fig, "/heatmap_clusters_pathway.pdf"), width=18, height=13)

bar_path[["BP"]] + scale_fill_manual( values= rep ("#67A9CF", 10) ) +ylab("NES") +
  bar_path[["MF"]] + scale_fill_manual( values= rep ("#CCC7C7", 10) )+ylab("NES") + 
  bar_path[["REACTOME"]] + scale_fill_manual( values= rep ("#EF8A62", 10) )+ylab("NES")

dev.off()



################ multi-sample stack plots 
laym = "
AAAAAAAABBBBBB
AAAAAAAACCCCCC
AAAAAAAADDD###
"
pdf(paste0(spcg_fig, "/multi_sample.pdf"), width=25, height=18)
rw_only_plot + w_only_plot+r_only_plot+ multi_t1  + plot_layout(design = laym)
dev.off()

# the legend is here 
multi_legend

## treemap version to tabulate the table 
pdf(paste0(spcg_fig, "/multi_tree_tabulate.pdf"), width=25, height=18)
tree.long.tab
dev.off()





## diagnostic fusions 


## fusions with outlier 
pdf(paste0(spcg_fig, "/fusion_outlier.pdf"), width=20, height=12)
fuse_outlier
dev.off()

## fusion gene body 
pdf(paste0(spcg_fig, "/fusion_genebody.pdf"), width=12.8, height=12.3)
balloonplot_type
dev.off()

## diagnostic fusions 
pdf(paste0(spcg_fig, "/diagn_fusions.pdf"), width=9, height=6)
well.known.diag
dev.off()


### fusion 3' domains. 
pdf(paste0(spcg_fig, "/fusion_domain_3prime.pdf"), width=17, height=9)
 ( in_frame_domain$right_plot + ggtitle ("3' Domain in.frame") )   + ( out_frame_domain$right_plot + ggtitle ("3' Domain out.of.frame") ) 
dev.off()

### :: network 
in_domain_net = draw_network1(in_domain, layout="kk")
render_graph(in_domain_net, layout = "fr", height=1500, width=2000  )

out_domain_net = draw_network1(out_domain, layout="kk")
render_graph(out_domain_net, layout = "fr")

###########

### itd 

pdf(paste0(spcg_fig, "/itd.pdf"), width=6, height=6)
itd2_plot
dev.off()

## top 6 matches ( manually save here or export to pdf )
pdf(paste0(spcg_fig, "/top6.pdf"), width=6, height=6)
top6_t1
dev.off()


kable(  df6 , format = "html" , row.names = F, caption="", table.attr = "style='width:30%;'" ) %>% kable_classic(full_width = T, position = "center" ,bootstrap_options = c("striped") )


##############
# fusion tally 

# 3399

pdf(paste0(spcg_fig, "/fuse_tally.pdf"), width=12, height=6)
fuse_tally_tbl
dev.off()


fusion_percent # 3451 

# recurrent fusions 
# line 3952
fusion_left.plot | fusion_right.plot




kable( fusion_tally , format = "html" , row.names = F, caption="", table.attr = "style='width:40%;'" ) %>% kable_classic(full_width = T, position = "center" ,bootstrap_options = c("striped") )%>%
  footnote(
    symbol  = c("CS: Cancer Subtype", "UnqP: Unique by patient", "UnqCS: Unique by cancer subtype", "pct.pos: Percent of Patient samples with at least 1 in.frame fusion")) %>% row_spec(0, background = "#e6e6e6" )  %>%
  row_spec( seq ( 2, nrow (fusion_tally ), by=2) )  %>%
  column_spec( 1, border_right = T ) 


# fusion compare niave to treated
pdf(paste0(spcg_fig, "/fuse_cmp.pdf"), width=7.4, height=4)
fuse.by.class_plot + fuse.by.treatment
dev.off()


### fusion + wgs 

wgs_is.obsv = plot_wgs_fusion( df = temp.all , type = "is.obsv", legend_pos = "bottom", fs=18)
wgs_is.diag = plot_wgs_fusion( df = temp.all , type = "is.diag", legend_pos = "none" , fs=18)
wgs_is.novel = plot_wgs_fusion( df = temp.all , type = "is.novel", legend_pos = "none" , fs=18 )



pdf(paste0(spcg_fig, "/fuse_wgs.pdf"), width=12, height=8)
wgs_is.diag  +   wgs_is.obsv + wgs_is.novel 
dev.off()




### fusion drug table 

layoutd <- 
"AAAAAAA##B"
pdf(paste0(spcg_fig, "/outlier_drugtable.pdf"), width=7.4, height=4)
drugpie_t1 
dev.off()


pdf(paste0(spcg_fig, "/outlier_drug_pie.pdf"), width=12.4, height=7)
drug_per_patient
dev.off()


# druggable table

pdf(paste0(spcg_fig, "/fusion_drug_pie.pdf"), width=8, height=3)
is_drug_tbl 
dev.off()



### fusion scatter expression 




layout <- "
AAC
BBD

"


pdf(paste0(spcg_fig, "/fusion_scatter.pdf"), width=12, height=10)

(right_in$expfusion.scat.plot + ggtitle ( "3' Gene in-frame") + ylim(miny,maxy) + xlim(minx, maxx) + theme( legend.position = "none")  )+
  (right_out$expfusion.scat.plot + ggtitle ( "3' Gene out-of-frame") + ylim(miny,maxy)  + xlim(minx, maxx) + theme( legend.position = "none") ) +
 
  (right_in$pie.ratio + ggtitle ( "") + theme( legend.position = "right")  )+
  (right_out$pie.ratio + ggtitle ( "") + theme( legend.position = "right") ) +
  plot_layout(design = layout) 


dev.off()








pdf(paste0(spcg_fig, "/fusion_scatter_violin.pdf"), width=10, height=8)
(right_in$exp.fusion.plot + ggtitle ( "") + theme( legend.position = "right")  ) 
dev.off()

######## outlier pie charts


pdf(paste0(spcg_fig, "/outlier_pie.pdf"), width=18, height=15)
up_plots$pie.out2 + down_plots$pie.out2
dev.off()

########## outlier family grid 


pdf(paste0(spcg_fig, "/outlier_family.pdf"), width=8, height=8)
out_family_grid
dev.off()


########## outlier with cna 

layout <- "
AAAC
BBBD

"


pdf(paste0(spcg_fig, "/outlier_scatter.pdf"), width=12, height=10)

( expCNA.scat.plot.gain + theme ( legend.position = "none") + xlab("") ) +
  expCNA.scat.plot.loss + pie.gain + pie.loss  + plot_layout(design = layout)

dev.off()


pdf(paste0(spcg_fig, "/outlier_scatter_exprs.pdf"), width=9, height=9)

( cna.tab.gain.plot2 + ggtitle ( "GAIN: top most significant genes" ) ) /
  cna.tab.loss.plot2 + ggtitle ( "LOSS: top most significant genes" )
dev.off()



######### outlier with drugs 

l1_d = plot_drug_pie(up.drugName = up.drugName , l = "L1", ft=15, c=5)
l2_d = plot_drug_pie(up.drugName = up.drugName , l = "L2", ft=40, c=1)
l3_d = plot_drug_pie(up.drugName = up.drugName , l = "L3", ft=15 , c=4)
l1_2 =  plot_drug_pie(up.drugName = up.drugName , l = c ( "L1", "L2"), ft=15, c=5)

pdf(paste0(spcg_fig, "/outlier_drugs.pdf"), width=15, height=12)
l1_2 / l3_d  + plot_layout(heights  =c ( .4,.4 ) )
dev.off()

lay2 = "
AAAAAAAAAAAA##
AAAAAAAAAAAA##
###BB#########
###BB#########
CCCCCCCCCCCC##
CCCCCCCCCCCC##
CCCCCCCCCCCC##
CCCCCCCCCCCC##
"
pdf(paste0(spcg_fig, "/outlier_drugs.pdf"), width=17, height=12)
l1_d + l2_d + l3_d + plot_layout(design = lay2)
dev.off()

### outlier kinase 
pdf(paste0(spcg_fig, "/outlier_kinase.pdf"), width=9, height=7.7)
out_kinase_stack
dev.off()

#### raf daf main plot 
pdf(paste0(spcg_fig, "/raf_daf.pdf"), width=15, height=12)
daf_vs_raf
dev.off()

#### 

#raf1  + l2_raf_HIGH2    +  l2_raf_HIGH3  

p2ae = l2_raf_HIGH2 + scale_fill_manual(values =  c( "HIGH/MODERATE" = '#ed557e', 'LOW/MODIFIER' ='#9E6977') )
p3ae = l2_raf_HIGH3 + scale_fill_manual(values =  c( "is.cancer" = '#f54272', 'not.cancer' ='#6e5d60') )
 


pdf(paste0(spcg_fig, "/raf_daf_pies.pdf"), width=25, height=8)
 
p2ae    +  p3ae 
dev.off()

# by disease 

af_port_plot_nosig

rafdaf_cosmic_plot + rafdaf_oncTSG_plot
var_outcome


# retrace colors? 


## expressed variants

#raf1  + l2_raf_HIGH2    +  l2_raf_HIGH3    
#by_disease + by_disease_ins + plot_layout(widths  =c ( .6,.4 ) )




#dafraf_plot_all_log
#dafraf_plot_all









pdf(paste0(spcg_fig, "/raf_daf_topGenes.pdf"), width=7, height=12)
( snv_gene_cancer + xlab("")  )/ snv_gene_non_cancer + ylab ( "Frequency")
dev.off()


pdf(paste0(spcg_fig, "/raf_daf_disease.pdf"), width=7, height=9)
by_disease_ins
dev.off()



## longitudinal plots 



fuse_longout$singleg
draw ( fuse_longout$ht ) 

draw ( outUP_longout$ht )

## shared 
long_fusion_ofInt = psb_fusions [ (psb_fusions$type.x == "both" & psb_fusions$is.cancer == "yes" | psb_fusions$type.y == "is.diag") & psb_fusions$left != psb_fusions$right , ]
long_fusion_ofInt=long_fusion_ofInt[ order ( long_fusion_ofInt$type.y), ]
kable ( long_fusion_ofInt[!duplicated ( paste(long_fusion_ofInt$pid, long_fusion_ofInt$info )) , c("pid", "info","cancer.subtype", "type.y")], format="html", row.names=F, table.attr = "style='width:40%;'" ) %>%
  kable_classic(full_width = F, position = "center", )  %>%
  column_spec(1, bold = T, border_right = T) 

## secondary 

long_fusion_ofInt = psb_fusions [ (psb_fusions$type.x == "secondary_only" & psb_fusions$is.cancer == "yes" & psb_fusions$frameness == "in.frame" | psb_fusions$type.y == "is.obsv") & psb_fusions$left != psb_fusions$right , ]
long_fusion_ofInt=long_fusion_ofInt[ order ( long_fusion_ofInt$type.y), ]
kable ( long_fusion_ofInt[!duplicated ( paste(long_fusion_ofInt$pid, long_fusion_ofInt$info )) , c("pid", "info","cancer.subtype", "type.y")], format="html", row.names=F, table.attr = "style='width:40%;'" ) %>%
  kable_classic(full_width = F, position = "center", )  %>%
  column_spec(1, bold = T, border_right = T) 


############### excel 


kable( diag_table[  , c(
  "RNAseq.id", "patient.id",  "fusion" , "diagnosis" , "cancer.subtype"
)
] , format = "html" , row.names = F, caption="", table.attr = "style='width:40%;'" ) %>% kable_classic(full_width = T, position = "center", )  %>%
  column_spec(1, bold = T, border_right = T) %>%
  row_spec(rthis, background = "#edebeb")


addWorksheet(paper_table, 'Diagnostic_fusion')
writeData(paper_table
          ,'Diagnostic_fusion'
          ,diag_table[  , c(
            "RNAseq.id", "patient.id",  "fusion" , "diagnosis" , "cancer.subtype"
          )
          ] 
          , rowNames=F  )



### 
kable( is.drug , format = "html" , row.names = F, caption="", table.attr = "style='width:40%;'" ) %>% kable_classic(full_width = T, position = "center", )  %>% column_spec(1, bold = T, border_right = T)


addWorksheet(paper_table, 'Druggable_fusion')
writeData(paper_table
          ,'Druggable_fusion'
          ,is.drug
          , rowNames=F  )








### 
kable( top6_results , format = "html" , row.names = F, caption="", table.attr = "style='width:40%;'" ) %>% kable_classic(full_width = T, position = "center", )  %>% column_spec(1, bold = T, border_right = T)

plot_top6







kable( drugpie , format = "html" , row.names = F, caption="", table.attr = "style='width:40%;'" ) %>% kable_classic(full_width = F, position = "center", )  %>%
  column_spec(1, bold = T, border_right = T) %>%
  row_spec(2, background = drug.color["L1"])  %>%
  row_spec(3, background = drug.color["L2"])  %>%
  row_spec(4, background = drug.color["L3"])  

drug_per_patient = make.pie( drugpie2[, 1:2]
                             , name.first="Level"
                             , cc = drug.color
                             , leg.pos = "right"
                             , title = "Drug Levels"
                             , pietext = 30
)






out_drug_grid + theme(
  axis.text.y = element_text(size= 35 ),
  axis.text.x =  element_text(size= 20 ),
)
  
grid.arrange (  out_drug_grid2 
)










## 
tab.fusion




# stack diagnosis 

well.known.diag


## just outliers



# outlier 

out_family_grid

up_plots$pie.out2 + down_plots$pie.out2

(up_oncogene + down_tsg ) / 

(out_up_grid + out_down_grid)

top_kinase


cancer_tab_f2_plot
nocancer_tab_f2_plot


out.total.bar.up

top_drug
drugpie_t1
waffle_drugs

lay2 = "
AAAAAAAAAAAA##
AAAAAAAAAAAA##
###BB#########
###BB#########
CCCCCCCCCCCC##
CCCCCCCCCCCC##
CCCCCCCCCCCC##
CCCCCCCCCCCC##
"
l1_d + l2_d + l3_d + plot_layout(design = lay2)



( lolli.up + ggtitle ("Oncogene Up" ) )  | fusion.plot$lolli_med + ggtitle ( "Fusions (in.frame )")



# making table for henry 

hn = diag_table[ , c("RNAseq.id", "fusion")]
hn$actionable = "diagnostic"
hn$druggable = 0 
colnames ( hn )[2] = "gene"
hn$alteration = "fusion"

hn2 = is.drug[ , c("RNAseq.id", "fusion")]
hn2$actionable = '' 
hn2$druggable = 1
colnames ( hn2 )[2] = "gene"
hn2$alteration = "fusion"

hn = rbind ( hn, hn2)


d = this.up_drug_tab[ , c("RNAseq.id", "gene", "level")]
d$druggable = 1 
d$actionable = ''
d$alteration = 'overexpression'


hn2 = rbind ( hn2, d [ , names ( hn2)])
tail ( hn2 )

write.table(hn2[, c("RNAseq.id", "gene", "alteration", "actionable", "druggable")],paste0( "combined.tsv"), col.names=TRUE, append = FALSE, sep = "\t",quote=FALSE)





