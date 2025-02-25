# Unifying all of the cell type annotations
# conda activate seurat5
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/
library(Seurat)
library(tidyverse)
library(ggrastr)
set.seed(1234)
out_dir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/"
obj <- readRDS("PDAC_merge_primary_KRAS_20241205_single_SCT.rds")
prefix = "PDAC_merge_primary_KRAS_20241205_single_SCT"
# Used the 50PC UMAP reduction as input since it seems to have the fewest batch effect
pdf(paste0(out_dir,"/","Dimplots_",prefix,"_50PC_res0.3_split.by_seurat_cluster_20241212.pdf"),useDingbats = F, width=100, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", split.by = "seurat_clusters", label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
dev.off()
cell_type_ind_by_sc <- as.matrix(table(obj$cell_type_individual, obj$seurat_clusters))
write.table(cell_type_ind_by_sc, paste0(out_dir,"/Table_cell_type_ind_by_sc_50PC",prefix,".tsv"),sep='\t',quote=F)
sample_vector <- unique(obj$sample_ID)
meta.data <- obj@meta.data
dir.create(paste0(out_dir,"/","cell_barcode_csv_files"))
for (sample_ID in sample_vector) {
    tmp_df <- meta.data[(meta.data[,"sample_ID"] == sample_ID),c('original_Xenium_barcode','seurat_clusters')]
    colnames(tmp_df) <- c('cell_id', 'group')
    write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",sample_ID,"_seurat_clusters.csv"),sep=",",quote=F,row.names=F)
    colnames(tmp_df) <- c('original_Xenium_barcode', 'seurat_cluster')
    write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",sample_ID,"_seurat_clusters.tsv"),sep="\t",quote=F)
}
tmp_df <- obj@meta.data[,c('barcode',"seurat_clusters")]                                       
colnames(tmp_df) <- c('barcode', 'seurat_clusters')                                                                                                                       
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_seurat_cluster",".csv"),sep=",",quote=F,row.names=F)
print(as.matrix(table(obj$sample_ID, obj$seurat_clusters)))
#                          0      1      2      3      4      5      6      7
# HT060P1-S1R1Fp1U1   131192  43520  40736  28645   1101  21198   9261   6433
# HT061P1-S1P1A1L1U1    8580   5031   6421    596     65   3428    402   1449
# HT061P1-S1P1A1L4U1   16441  19216  10904  10288    194   5254   4897  15769
# HT125P1-S1H4A1L1U1   32072   8426  14484   3565    118   5309    240   6564
# HT125P1-S1H8A1U1      8334   7965   8770     98     76   3342   1022   7622
# HT185P1-S1H2L1U1     20621   7769   4686  10813    624   4524   7359   4108
# HT227P1-S1H1L1U1     13126  10674   6006  12909    302   2672   3344    496
# HT242P1-S1H4L4U1     10301   7557  11674     60  60461   3005  17127    517
# HT270P1-S1H1A1US2_1  14893   5778   8071   3401    111   4175   3887   1765
# 
#                          8      9     10     11     12     13     14     15
# HT060P1-S1R1Fp1U1    14282  14765  25215   1909  11652   7193   5937   1456
# HT061P1-S1P1A1L1U1    1386   7031   1142     75   1462    118   1728    399
# HT061P1-S1P1A1L4U1    3056   1104    108   9445    819    391   3496    198
# HT125P1-S1H4A1L1U1    4239   2406    524    178   3713    586   1441   2255
# HT125P1-S1H8A1U1      3475    448   3839    117    695    458    523   3289
# HT185P1-S1H2L1U1      3780   1134     70   2915    742   1232    445    225
# HT227P1-S1H1L1U1      1553    679     39    281    422   6437    784    350
# HT242P1-S1H4L4U1      2066   3181     72    240   1100   3298   3087   6278
# HT270P1-S1H1A1US2_1   3150    541    162  12769   1084    220    331    259
# 
#                         16     17     18     19     20     21     22     23
# HT060P1-S1R1Fp1U1      885   4959   2124   3352   1420   4792   2699   1697
# HT061P1-S1P1A1L1U1    8289    152    151    450    335    176     97     20
# HT061P1-S1P1A1L4U1     240    130    581    773    459    155     92     38
# HT125P1-S1H4A1L1U1     233    535   2073    616    527    753    287    246
# HT125P1-S1H8A1U1        19    176    633    516    218    365    199     76
# HT185P1-S1H2L1U1       218   2418    856    456    179    329   1326   1435
# HT227P1-S1H1L1U1       126   1302    383   1053    373    332    844    644
# HT242P1-S1H4L4U1        44    944    862   1060   4376    380   1399    855
# HT270P1-S1H1A1US2_1   2934    232    976    180    168    246    128     68
# 
#                         24     25     26     27     28     29     30     31
# HT060P1-S1R1Fp1U1     1361   2143    549    694    500    351    671    440
# HT061P1-S1P1A1L1U1     239     72      4     36     33    107      7     36
# HT061P1-S1P1A1L4U1     384    174     22     15     53     34    153     72
# HT125P1-S1H4A1L1U1     271    147     22    204     79     80     20    168
# HT125P1-S1H8A1U1       503     85     15     56     52     61    391     77
# HT185P1-S1H2L1U1       811    564     40    598    327    741     71    441
# HT227P1-S1H1L1U1       375    436      5    236    270    209     19    147
# HT242P1-S1H4L4U1       799    706   1365     88    622     78     39     96
# HT270P1-S1H1A1US2_1    286    103     30     88     46     61    201     52
# 
#                         32     33     34
# HT060P1-S1R1Fp1U1      216     16      4
# HT061P1-S1P1A1L1U1       5      1      0
# HT061P1-S1P1A1L4U1      10      0      0
# HT125P1-S1H4A1L1U1      17      3      1
# HT125P1-S1H8A1U1        11      2      1
# HT185P1-S1H2L1U1       147     18      5
# HT227P1-S1H1L1U1        31     14      4
# HT242P1-S1H4L4U1        52      9      4
# HT270P1-S1H1A1US2_1      6      2      0
test <- as.matrix(table(obj$sample_ID, obj$seurat_clusters))
test <- cbind(test, total = rowSums(test))
test/test[,"total"]*100
test[,"total"] <- NULL
print(test)
#                             0         1         2           3          4                                                                                                                  
# HT060P1-S1R1Fp1U1   33.350959 11.063432 10.355697  7.28198532  0.2798906                                                                                                                  
# HT061P1-S1P1A1L1U1  17.325283 10.158916 12.965693  1.20348121  0.1312521                                                                                                                  
# HT061P1-S1P1A1L4U1  15.663316 18.307055 10.388225  9.80136236  0.1848235                                                                                                                  
# HT125P1-S1H4A1L1U1  34.709205  9.118850 15.674985  3.85814160  0.1277029                                                                                                                  
# HT125P1-S1H8A1U1    15.569131 14.879785 16.383643  0.18307833  0.1419791                                                                                                                  
# HT185P1-S1H2L1U1    25.139283  9.471272  5.712753 13.18224487  0.7607251                                                                                                                  
# HT227P1-S1H1L1U1    19.627077 15.960644  8.980666 19.30260030  0.4515753                                                                                                                  
# HT242P1-S1H4L4U1     7.163322  5.255142  8.118107  0.04172404 42.0446169                                                                                                                  
# HT270P1-S1H1A1US2_1 22.427866  8.701283 12.154388  5.12167942  0.1671586                                                                                                                  
#                            5          6          7        8          9                                                                                                                    
# HT060P1-S1R1Fp1U1   5.388847  2.3542840  1.6353643 3.630697  3.7534827                                                                                                                    
# HT061P1-S1P1A1L1U1  6.922036  0.8117440  2.9259132 2.798700 14.1974436                                                                                                                    
# HT061P1-S1P1A1L4U1  5.005478  4.6653646 15.0231029 2.911447  1.0517792                                                                                                                    
# HT125P1-S1H4A1L1U1  5.745547  0.2597346  7.1037423 4.587563  2.6038397                                                                                                                    
# HT125P1-S1H8A1U1    6.243345  1.9092455 14.2390106 6.491808  0.8369295                                                                                                                    
# HT185P1-S1H2L1U1    5.515257  8.9714362  5.0081071 4.608239  1.3824716                                                                                                                    
# HT227P1-S1H1L1U1    3.995395  5.0002243  0.7416601 2.322174  1.0152967                                                                                                                    
# HT242P1-S1H4L4U1    2.089679 11.9101264  0.3595221 1.436698  2.2120694                                                                                                                    
# HT270P1-S1H1A1US2_1 6.287272  5.8535630  2.6579724 4.743690  0.8147100 
#                             10         11        12        13        14                                                                                                                   
# HT060P1-S1R1Fp1U1   6.41002827  0.4852962 2.9621118 1.8285677 1.5092738                                                                                                                   
# HT061P1-S1P1A1L1U1  2.30599923  0.1514448 2.9521636 0.2382731 3.4892878                                                                                                                   
# HT061P1-S1P1A1L4U1  0.10289144  8.9982375 0.7802601 0.3725051 3.3306340                                                                                                                   
# HT125P1-S1H4A1L1U1  0.56708729  0.1926365 4.0183113 0.6341854 1.5594901                                                                                                                   
# HT125P1-S1H8A1U1    7.17181341  0.2185731 1.2983616 0.8556110 0.9770405                                                                                                                   
# HT185P1-S1H2L1U1    0.08533775  3.5537079 0.9045802 1.5019445 0.5425043
# HT227P1-S1H1L1U1    0.05831601  0.4201743 0.6310092 9.6251327 1.1723014
# HT242P1-S1H4L4U1    0.05006884  0.1668961 0.7649407 2.2934312 2.1467017
# HT270P1-S1H1A1US2_1 0.24396121 19.2292633 1.6324318 0.3313053 0.4984639
#                            15          16        17        18        19
# HT060P1-S1R1Fp1U1   0.3701369  0.22498017 1.2606516 0.5399524 0.8521283
# HT061P1-S1P1A1L1U1  0.8056862 16.73767744 0.3069281 0.3049088 0.9086687
# HT061P1-S1P1A1L4U1  0.1886343  0.22864764 0.1238508 0.5535178 0.7364360
# HT125P1-S1H4A1L1U1  2.4404234  0.25215904 0.5789918 2.2434579 0.6666522
# HT125P1-S1H8A1U1    6.1443330  0.03549478 0.3287937 1.1825366 0.9639635
# HT185P1-S1H2L1U1    0.2742999  0.26576615 2.9478099 1.0435588 0.5559145
# HT227P1-S1H1L1U1    0.5233488  0.18840558 1.9468577 0.5726932 1.5745324
# HT242P1-S1H4L4U1    4.3657251  0.03059763 0.6564582 0.5994353 0.7371247
# HT270P1-S1H1A1US2_1 0.3900367  4.41840853 0.3493765 1.4697910 0.2710680
#                            20        21         22         23        24
# HT060P1-S1R1Fp1U1   0.3609851 1.2181977 0.68612597 0.43140266 0.3459865
# HT061P1-S1P1A1L1U1  0.6764534 0.3553904 0.19586859 0.04038528 0.4826040
# HT061P1-S1P1A1L4U1  0.4372886 0.1476683 0.08764826 0.03620254 0.3658362
# HT125P1-S1H4A1L1U1  0.5703340 0.8149174 0.31059934 0.26622800 0.2932837
# HT125P1-S1H8A1U1    0.4072559 0.6818734 0.37176110 0.14197911 0.9396776
# HT185P1-S1H2L1U1    0.2182208 0.4010874 1.61654090 1.74942397 0.9886988
# HT227P1-S1H1L1U1    0.5577403 0.4964338 1.26201833 0.96296186 0.5607309
# HT242P1-S1H4L4U1    3.0430731 0.2642522 0.97286547 0.59456753 0.5556251
# HT270P1-S1H1A1US2_1 0.2529968 0.3704596 0.19275947 0.10240347 0.4306969
#                            25          26         27         28         29
# HT060P1-S1R1Fp1U1   0.5447825 0.139563971 0.17642513 0.12710744 0.08922942
# HT061P1-S1P1A1L1U1  0.1453870 0.008077055 0.07269350 0.06663570 0.21606122
# HT061P1-S1P1A1L4U1  0.1657695 0.020959367 0.01429048 0.05049302 0.03239175
# HT125P1-S1H4A1L1U1  0.1590875 0.023809008 0.22077444 0.08549598 0.08657821
# HT125P1-S1H8A1U1    0.1587924 0.028022194 0.10461619 0.09714360 0.11395692
# HT185P1-S1H2L1U1    0.6875785 0.048764431 0.72902825 0.39864923 0.90336109
# HT227P1-S1H1L1U1    0.6519431 0.007476412 0.35288664 0.40372624 0.31251402
# HT242P1-S1H4L4U1    0.4909528 0.949221847 0.06119525 0.43253919 0.05424125
# HT270P1-S1H1A1US2_1 0.1551111 0.045178001 0.13252214 0.06927294 0.09186194
#                             30         31          32          33          34
# HT060P1-S1R1Fp1U1   0.17057819 0.11185455 0.054910415 0.004067438 0.001016860
# HT061P1-S1P1A1L1U1  0.01413485 0.07269350 0.010096319 0.002019264 0.000000000
# HT061P1-S1P1A1L4U1  0.14576287 0.06859429 0.009526985 0.000000000 0.000000000
# HT125P1-S1H4A1L1U1  0.02164455 0.18181425 0.018397870 0.003246683 0.001082228
# HT125P1-S1H8A1U1    0.73044518 0.14384726 0.020549609 0.003736292 0.001868146
# HT185P1-S1H2L1U1    0.08655687 0.53762785 0.179209285 0.021943994 0.006095554
# HT227P1-S1H1L1U1    0.02841037 0.21980651 0.046353754 0.020933953 0.005981130
# HT242P1-S1H4L4U1    0.02712062 0.06675846 0.036160832 0.006258606 0.002781602
# HT270P1-S1H1A1US2_1 0.30269261 0.07830854 0.009035600 0.003011867 0.000000000
test <- as.matrix(table(obj$sample_ID, obj$seurat_clusters))
test <- rbind(test, total = colSums(test))
for (col in colnames(test)) {
    test[,col] <- test[,col]/test["total",col]*100
}
print(test)
#                              0          1          2            3           4                                                                                                             
# HT060P1-S1R1Fp1U1    51.335107  37.537952  36.452144  40.70337478   1.7461778                                                                                                             
# HT061P1-S1P1A1L1U1    3.357333   4.339463   5.745758   0.84689165   0.1030895                                                                                                             
# HT061P1-S1P1A1L4U1    6.433323  16.574662   9.757320  14.61882771   0.3076825                                                                                                             
# HT125P1-S1H4A1L1U1   12.549695   7.267803  12.960842   5.06571936   0.1871471                                                                                                             
# HT125P1-S1H8A1U1      3.261074   6.870170   7.847734   0.13925400   0.1205354                                                                                                             
# HT185P1-S1H2L1U1      8.068947   6.701111   4.193214  15.36483126   0.9896593                                                                                                             
# HT227P1-S1H1L1U1      5.136172   9.206804   5.374400  18.34316163   0.4789697                                                                                                             
# HT242P1-S1H4L4U1      4.030756   6.518251  10.446345   0.08525755  95.8906934                                                                                                             
# HT270P1-S1H1A1US2_1   5.827594   4.983784   7.222242   4.83268206   0.1760452                                                                                                             
# total               100.000000 100.000000 100.000000 100.00000000 100.0000000 
#                              5           6          7          8          9                                                                                                               
# HT060P1-S1R1Fp1U1    40.066532  19.4808473  14.384098  38.613567  47.189108                                                                                                               
# HT061P1-S1P1A1L1U1    6.479294   0.8456215   3.239944   3.747263  22.471156                                                                                                               
# HT061P1-S1P1A1L4U1    9.930633  10.3010160  35.259263   8.262362   3.528397                                                                                                               
# HT125P1-S1H4A1L1U1   10.034589   0.5048487  14.677012  11.460784   7.689603                                                                                                               
# HT125P1-S1H8A1U1      6.316744   2.1498138  17.042685   9.395193   1.431813                                                                                                               
# HT185P1-S1H2L1U1      8.550853  15.4799217   9.185430  10.219807   3.624277                                                                                                               
# HT227P1-S1H1L1U1      5.050371   7.0342245   1.109049   4.198773   2.170092                                                                                                               
# HT242P1-S1H4L4U1      5.679778  36.0272618   1.156005   5.585746  10.166512                                                                                                               
# HT270P1-S1H1A1US2_1   7.891205   8.1764446   3.946515   8.516506   1.729042                                                                                                               
# total               100.000000 100.0000000 100.000000 100.000000 100.000000                                                                                                               
#                              10          11         12          13         14                                                                                                             
# HT060P1-S1R1Fp1U1    80.8924962   6.8351892  53.723085  36.0858877  33.406482                                                                                                             
# HT061P1-S1P1A1L1U1    3.6636617   0.2685381   6.740744   0.5919831   9.723160                                                                                                             
# HT061P1-S1P1A1L4U1    0.3464759  33.8178954   3.776108   1.9615713  19.671393                                                                                                             
# HT125P1-S1H4A1L1U1    1.6810497   0.6373304  17.119277   2.9398485   8.108260                                                                                                             
# HT125P1-S1H8A1U1     12.3159347   0.4189194   3.204389   2.2976973   2.942831                                                                                                             
# HT185P1-S1H2L1U1      0.2245677  10.4371800   3.421089   6.1807054   2.503939                                                                                                             
# HT227P1-S1H1L1U1      0.1251163   1.0061227   1.945687  32.2931822   4.411434                                                                                                             
# HT242P1-S1H4L4U1      0.2309839   0.8593219   5.071695  16.5454272  17.370020                                                                                                             
# HT270P1-S1H1A1US2_1   0.5197138  45.7195030   4.997925   1.1036974   1.862480                                                                                                             
# total               100.0000000 100.0000000 100.000000 100.0000000 100.000000
#                             15          16         17         18         19                                                                                                               
# HT060P1-S1R1Fp1U1     9.898701   6.8139821  45.713496  24.586179  39.640492                                                                                                               
# HT061P1-S1P1A1L1U1    2.712625  63.8204496   1.401180   1.747887   5.321665                                                                                                               
# HT061P1-S1P1A1L4U1    1.346115   1.8478596   1.198378   6.725315   9.141438                                                                                                               
# HT125P1-S1H4A1L1U1   15.330750   1.7939637   4.931785  23.995833   7.284768                                                                                                               
# HT125P1-S1H8A1U1     22.360460   0.1462889   1.622419   7.327237   6.102176                                                                                                               
# HT185P1-S1H2L1U1      1.529676   1.6784724  22.289823   9.908554   5.392621                                                                                                               
# HT227P1-S1H1L1U1      2.379496   0.9701263  12.002212   4.433383  12.452696                                                                                                               
# HT242P1-S1H4L4U1     42.681352   0.3387743   8.702065   9.978007  12.535478                                                                                                               
# HT270P1-S1H1A1US2_1   1.760827  22.5900832   2.138643  11.297604   2.128666                                                                                                               
# total               100.000000 100.0000000 100.000000 100.000000 100.000000                                                                                                               
#                             20         21         22          23         24                                                                                                               
# HT060P1-S1R1Fp1U1    17.628802  63.655685  38.169990  33.4120890  27.063034                                                                                                               
# HT061P1-S1P1A1L1U1    4.158908   2.337938   1.371800   0.3937783   4.752436                                                                                                               
# HT061P1-S1P1A1L4U1    5.698324   2.058980   1.301089   0.7481788   7.635713                                                                                                               
# HT125P1-S1H4A1L1U1    6.542520  10.002657   4.058832   4.8434731   5.388745                                                                                                               
# HT125P1-S1H8A1U1      2.706394   4.848565   2.814312   1.4963576  10.001988                                                                                                               
# HT185P1-S1H2L1U1      2.222222   4.370351  18.752652  28.2535932  16.126466                                                                                                               
# HT227P1-S1H1L1U1      4.630664   4.410202  11.936077  12.6796614   7.456751                                                                                                               
# HT242P1-S1H4L4U1     54.326505   5.047821  19.785037  16.8340224  15.887850                                                                                                               
# HT270P1-S1H1A1US2_1   2.085661   3.267800   1.810211   1.3388462   5.687015                                                                                                               
# total               100.000000 100.000000 100.000000 100.0000000 100.000000
#                             25          26          27         28         29
# HT060P1-S1R1Fp1U1    48.374718  26.7543860  34.4416873  25.227043  20.383275
# HT061P1-S1P1A1L1U1    1.625282   0.1949318   1.7866005   1.664985   6.213705
# HT061P1-S1P1A1L4U1    3.927765   1.0721248   0.7444169   2.674067   1.974448
# HT125P1-S1H4A1L1U1    3.318284   1.0721248  10.1240695   3.985873   4.645761
# HT125P1-S1H8A1U1      1.918736   0.7309942   2.7791563   2.623613   3.542393
# HT185P1-S1H2L1U1     12.731377   1.9493177  29.6774194  16.498486  43.031359
# HT227P1-S1H1L1U1      9.841986   0.2436647  11.7121588  13.622603  12.137050
# HT242P1-S1H4L4U1     15.936795  66.5204678   4.3672457  31.382442   4.529617
# HT270P1-S1H1A1US2_1   2.325056   1.4619883   4.3672457   2.320888   3.542393
# total               100.000000 100.0000000 100.0000000 100.000000 100.000000
#                              30         31         32         33         34
# HT060P1-S1R1Fp1U1    42.6844784  28.776978  43.636364  24.615385  21.052632
# HT061P1-S1P1A1L1U1    0.4452926   2.354480   1.010101   1.538462   0.000000
# HT061P1-S1P1A1L4U1    9.7328244   4.708960   2.020202   0.000000   0.000000
# HT125P1-S1H4A1L1U1    1.2722646  10.987574   3.434343   4.615385   5.263158
# HT125P1-S1H8A1U1     24.8727735   5.035971   2.222222   3.076923   5.263158
# HT185P1-S1H2L1U1      4.5165394  28.842381  29.696970  27.692308  26.315789
# HT227P1-S1H1L1U1      1.2086514   9.614127   6.262626  21.538462  21.052632
# HT242P1-S1H4L4U1      2.4809160   6.278613  10.505051  13.846154  21.052632
# HT270P1-S1H1A1US2_1  12.7862595   3.400916   1.212121   3.076923   0.000000
# total               100.0000000 100.000000 100.000000 100.000000 100.000000
summary(obj$nFeature_Xenium[obj$seurat_clusters %in% c(34)])
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.000   1.000   1.000   1.211   1.000   3.000
summary(obj$nFeature_Xenium[obj$seurat_clusters %in% c(0)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.00   25.00   37.00   37.12   48.00  123.00
summary(obj$nFeature_Xenium[obj$seurat_clusters %in% c(4)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.00   19.00   29.00   29.07   38.00  102.00
summary(obj$nFeature_Xenium[obj$seurat_clusters %in% c(7)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.00   40.00   55.00   53.86   69.00  131.00
summary(obj$nFeature_Xenium[obj$seurat_clusters %in% c(33)])
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   1.277   2.000   2.000
summary(obj$nFeature_Xenium[obj$seurat_clusters %in% c(33)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   1.277   2.000   2.000 
summary(obj$nFeature_Xenium[obj$seurat_clusters %in% c(32)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   3.000   6.000   7.947   9.000  66.000 
summary(obj$nFeature_Xenium[obj$seurat_clusters %in% c(31)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   3.000   6.000   7.924  10.000  87.000 
summary(obj$nFeature_Xenium[obj$seurat_clusters %in% c(30)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   13.00   34.00   33.95   50.00  104.00 
summary(obj$nFeature_Xenium[obj$seurat_clusters %in% c(29)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   4.000   6.000   8.738  10.000 104.000 
summary(obj$nFeature_Xenium[obj$seurat_clusters %in% c(28)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   3.000   5.000   6.247   8.000  75.000 
summary(obj$nFeature_Xenium[obj$seurat_clusters %in% c(27)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   3.000   5.000   7.201   8.000  94.000 
summary(obj$nFeature_Xenium[obj$seurat_clusters %in% c(26)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   23.00   35.00   32.55   43.00   95.00
summary(obj$nFeature_Xenium[obj$seurat_clusters %in% c(25)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    5.00    9.00   11.59   14.00  112.00
summary(obj$nFeature_Xenium[obj$seurat_clusters %in% c(24)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    5.00    9.00   12.14   16.00  100.00
summary(obj$nFeature_Xenium[obj$seurat_clusters %in% c(23)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   1.983   2.000  42.000 
summary(obj$nFeature_Xenium[obj$seurat_clusters %in% c(22)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   2.000   3.000   4.191   5.000 104.000
summary(obj$nFeature_Xenium[obj$seurat_clusters %in% c(17)])
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   3.000   4.000   6.478   8.000 107.000
summary(obj$nFeature_Xenium[obj$seurat_clusters %in% c(13)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   10.00   19.00   22.98   33.00  112.00 
table(obj$seurat_clusters)
#      0      1      2      3      4      5      6      7      8      9     10 
# 255560 115936 111752  70375  63052  52907  47539  44723  36987  31289  31171 
#     11     12     13     14     15     16     17     18     19     20     21 
#  27929  21689  19933  17772  14709  12988  10848   8639   8456   8055   7528 
#     22     23     24     25     26     27     28     29     30     31     32 
#   7071   5079   5029   4430   2052   2015   1982   1722   1572   1529    495 
#     33     34 
#     65     19
table(obj$seurat_clusters)/sum(table(obj$seurat_clusters))*100
#            0            1            2            3            4            5 
# 24.272079795 11.011143540 10.613763739  6.683939645  5.988430017  5.024897972 
#            6            7            8            9           10           11 
#  4.515066526  4.247613964  3.512879228  2.971705684  2.960498510  2.652586150 
#           12           13           14           15           16           17 
#  2.059935587  1.893157640  1.687914392  1.397002746  1.233548961  1.030300210 
#           18           19           20           21           22           23 
#  0.820498111  0.803117494  0.765032097  0.714979718  0.671575662  0.482383367 
#           24           25           26           27           28           29 
#  0.477634564  0.420743909  0.194890858  0.191376744  0.188242535  0.163548761 
#           30           31           32           33           34 
#  0.149302353  0.145218383  0.047013146  0.006173443  0.001804545
# Making sure that cell types are still consistent:
# 0 CAF (mix of myCAF and iCAF)
# 1 Myeloid
# 2 T_cell
# 3 Ductal (mix of Duct_like_1, Duct_like_2, Reactive_duct)
# 4 Tumor
# 5 Endothelial
# 6 Islet (GCG+ INS+ (equal))
# 7 Tumor
# 8 SMC
# 9 Plasma
# 10 Islet PPY+ (almost all), INS+ weak
# 11 Acinar
# 12 Mast
# 13 Mix (TNC+ (fibro), COL17A1+ (tumor), CXCL6+ (duct), GPC1+ (glial/SMC), SERPINB3+ )
# 14 Myeloid - most likely to be DC (IRF8+, FGL2+, CSF2RA+, )
# 15 B_cell
# 16 Epi - PanIN
# 17 LowCount (no positive DEGs and median nCount_Xenium is 4)
# 18 Glial
# 19 Proliferative (some Myeloid and T_cell markers here. my guess is this may split out by population if I subcluster it)
# 20 NK (CD3E+ CD3D- CD4-,GNLY+(specific),GZMA+,GZMB+(specific),KLRC1+) #https://pmc.ncbi.nlm.nih.gov/articles/PMC41376/ CD3E is the only part of the t cell receptor expressed in NK consistently
# 21 CAF (C7+ )
# 22 LowCount (no positive DEGs and median nCount_Xenium is 3)
# 23 LowCount (no positive DEGs and median nCount_Xenium is 1)
# 24 CAF_LowCount (ACTA2+ only positive DEG. sparse other myCAF markers from dotplots include THBS2+ (15%) and MFAP5+ (maybe 1%)
# 25 Myeloid_LowCount (CXCL2+ (high) and RGS16+ and nothing else) both of these are non-specific. roughly 50% of the cells are from HT060P1 and these are in immune dense regions. The dotplots do have some week signal in some myeloid genes. May try subclustering with other myeloid cells.
# 26 KRT20+SLC26A3+_Tumor (EPCAM+ (moderate), MET+ (high), COL17A1+ (high,specific), SLC26A3+ (specific) ANPEP+ (moderate) KRT20+ (specific)) localized to specific regions in HT242P1 which is 66% of this cluster.
# 27 CAF_LowCount (FBLN1+ high, SRPX+ low) this is likely fibroblasts but I can't tell what type. Will try to include with other fibroblast subclustering. May end up marking as Unknown)
# 28 Unknown (VCAN+ only DEG here. no other markers up on dotplots. VCAN is non-specific so going to mark this as Unknown)
# 29 CAF_LowCount (SFRP4+ highest not much else. this is likely fibroblasts but I can't tell what type. Will try to include with other fibroblast subclustering. May end up marking as Unknown)
# 30 Adipocyte (PLIN4+, ADIPOQ+, LPL+)
# 31 CAF_LowCount
# 32 LowCount Markers could indicate Erythrocyte AHSP+ ALAS2+ GYPA+ however upon inspections of the histology there were no visible blood cells here overlapping these cells. (CHGA+, PCSK2+, SCGN+) - no Islet hormones though. This label is accurate based on the markers but makes little sense because these cells should not be segmented since they lack nuclei and every cell in each sample in this object was only segmented based on the nuclei
# 33 LowCount MS4A4A+ and nothing else
# 34 LowCount (GNG+ and only that)
# Plan: 
# Subcluster proliferative and split it out into t_cell, myeloid, and tumor (if present).
# Cluster 13 seems to be a mixture of 4 different populations (CAF, Myeloid, Tumor, and LowCount) so will try to split it before subclustering other larger populations.
# Subcluster T_cell and split out
# Subcluster Myeloid and split out
# Subcluster CAF and split out
# Subcluster Acinar with duct_like_1, duct_like_2, and reactive_duct.
# subclustering cluster 19 independently of other clusters
library(future)
options(future.globals.maxSize= +Inf)
Idents(obj) <- "seurat_clusters"
plan(sequential)
subcluster_ident <- "merge.sub.19"
subcluster_idents <- c(19)
obj <- FindSubCluster(obj, cluster = subcluster_idents, graph.name = "SCT_snn", subcluster.name = subcluster_ident, resolution = 0.3, algorithm = 1) #based on the recommendation here https://github.com/satijalab/seurat/issues/6682 by the authors I used the SCT_snn instead of SCT_nn graph
Idents(obj) <- subcluster_ident
tmp_df <- obj@meta.data[,c('barcode',subcluster_ident)]
colnames(tmp_df) <- c('barcode',subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
# Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
# /diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
# ./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
# /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_merge.sub.19.tsv \
# merge.sub.19
# results in /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/marker_dotplots/merge.sub.19.PDAC_merge_primary_KRAS_20241205_single_SCT.Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.pdf
plan(multicore, workers = 10)
degs <- FindAllMarkers(obj, assay='SCT', slot='data', test.use='wilcox', logfc.threshold = 0, min.pct=0.1, min.diff.pct=0.1)
plan(sequential)
write.table(degs,paste0(out_dir,"/","FindAllMarkers_",subcluster_ident,"_",prefix,"_20241212.tsv"),sep="\t",quote=FALSE)
system(paste0("Rscript /diskmnt/Projects/Users/austins2/tools/markers/Shared_marker_comparison_from_seurat_degs.R -i ",paste0(out_dir,"/","FindAllMarkers_",unified_subcluster_ident,"_",prefix,"_20241212.tsv")))
pdf(paste0(out_dir,"/","Dimplot_",subcluster_ident,"_",prefix,"_res0.3_split.by_20241212.pdf"),useDingbats = F, width=100, height=6)
print(DimPlot(obj, reduction = "umap.50PC", split.by = subcluster_ident, label=TRUE, label.size=4, raster=FALSE))
dev.off()
pdf(paste0(out_dir,"/","Dimplot_",subcluster_ident,"_",prefix,"_res0.3_group.by_20241112.pdf"),useDingbats = F, width=8, height=6)
print(DimPlot(obj, reduction = "umap.50PC", group.by = subcluster_ident, label=TRUE, label.size=4, raster=FALSE))
dev.off()
library(scCustomize)
cells_to_highlight <- list()
merge_idents_highlight <- unique(obj@meta.data[,subcluster_ident])[grepl("_",unique(obj@meta.data[,subcluster_ident]))]
for (sub_ident in merge_idents_highlight) {
    print(sub_ident)
    #string_sub_ident <- paste0(sub_ident,"_")
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,subcluster_ident] == sub_ident),"barcode"]
    cells_to_highlight[[sub_ident]] <- barcodes_in_sub
    print(summary(obj@meta.data[,"nCount_Xenium"][obj@meta.data[,subcluster_ident] == sub_ident]))
}
# [1] "19_0"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.0    31.0    66.0    90.6   128.0   533.0 
# [1] "19_6"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 20.00   87.25  136.50  152.88  198.00  523.00 
# [1] "19_7"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 16.00   63.75   94.50  108.69  145.25  324.00 
# [1] "19_4"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 11.0    64.0   109.0   125.8   165.5   477.0 
# [1] "19_3"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    4.00    7.00   11.33   13.00  191.00 
# [1] "19_5"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 17.0    71.5   118.5   132.6   176.5   492.0 
# [1] "19_1"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8.0    46.0    85.0   102.2   142.0   463.0 
# [1] "19_2"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00   23.00   54.00   84.91  124.00  609.00
varibow_pal <- DiscretePalette_scCustomize(num_colors = length(cells_to_highlight), palette = "varibow")
pdf(paste0(out_dir,"/","Dimplot_",subcluster_ident,"_",prefix,"_res0.3_highlight_subclusters_20241112.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = subcluster_ident, cells.highlight = cells_to_highlight, cols.highlight = varibow_pal, sizes.highlight=0.5, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
dev.off()
print(table(obj$cell_type_individual[obj$seurat_clusters %in% subcluster_idents], obj$sample_ID[obj$seurat_clusters %in% subcluster_idents]))
print(table(obj$cell_type_individual[obj@meta.data[,c(subcluster_ident)] %in% merge_idents_highlight], obj@meta.data[,subcluster_ident][obj@meta.data[,c(subcluster_ident)] %in% merge_idents_highlight]))
table(obj@meta.data[,subcluster_ident])
#      0      1     10     11     12     13     14     15     16     17     18 
# 255560 115936  31171  27929  21689  19933  17772  14709  12988  10848   8639 
#   19_0   19_1   19_2   19_3   19_4   19_5   19_6   19_7      2     20     21 
#   2758   1429   1262   1258    591    580    362    216 111752   8055   7528 
#     22     23     24     25     26     27     28     29      3     30     31 
#   7071   5079   5029   4430   2052   2015   1982   1722  70375   1572   1529 
#     32     33     34      4      5      6      7      8      9 
#    495     65     19  63052  52907  47539  44723  36987  31289
for (sample_ID in sample_vector) {
    tmp_df <- obj@meta.data[(obj@meta.data[,"sample_ID"] == sample_ID),c('original_Xenium_barcode',subcluster_ident)]
    colnames(tmp_df) <- c('cell_id', 'group')
    write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",sample_ID,"_",prefix,"_",subcluster_ident,".csv"),sep=",",quote=F,row.names=F)
}
table(obj$sample_ID[obj@meta.data[,subcluster_ident] %in% merge_idents_highlight], obj@meta.data[,subcluster_ident][obj@meta.data[,subcluster_ident] %in% merge_idents_highlight])
#                     19_0 19_1 19_2 19_3 19_4 19_5 19_6 19_7
# HT060P1-S1R1Fp1U1    944  471  314  562  259  312  312  178
# HT061P1-S1P1A1L1U1   101  226   33   37   27   22    2    2
# HT061P1-S1P1A1L4U1   405  108  165   43   13   33    3    3
# HT125P1-S1H4A1L1U1   156  111   86   64  123   52   12   12 
# HT125P1-S1H8A1U1     145  105  118   53   67    7   19    2
# HT185P1-S1H2L1U1     124   28   96  185   11    6    1    5
# HT227P1-S1H1L1U1     549   93  187  137   65   12    8    2
# HT242P1-S1H4L4U1     267  266  228  153   12  123    3    8
# HT270P1-S1H1A1US2_1   67   21   35   24   14   13    2    4
# based on the above this is mostly proliferative immune cells (mix of T and myeloid) with a small population of proliferative fibroblasts?
# merge.sub.19
# 19_0 Myeloid
# 19_1 CD4+T_cell_proliferative
# 19_2 myCAF_proliferative # some of this might be tumor cells but it is a small proportion
# 19_3 LowCount # only a small proportion of this cluster is actually low count.
# 19_4 Myeloid
# 19_5 CD8+T_cell_proliferative
# 19_6 Myeloid
# 19_7 iCAF_proliferative
sub_tumor_individual <- obj@meta.data[,"barcode"][(obj@meta.data[,subcluster_ident] %in% c('19_0','19_2')) & (grepl("Tumor",obj$cell_type_individual))]
length(sub_tumor_individual)
# [1] 756
table(obj$sample_ID[(obj@meta.data[,subcluster_ident] %in% merge_idents_highlight) & (obj@meta.data[,'cell_type_individual'] %in% unique(obj$cell_type_individual)[grepl("Tumor",unique(obj$cell_type_individual))])], 
      obj@meta.data[,subcluster_ident][(obj@meta.data[,subcluster_ident] %in% merge_idents_highlight) & (obj@meta.data[,'cell_type_individual'] %in% unique(obj$cell_type_individual)[grepl("Tumor",unique(obj$cell_type_individual))])])
#                    19_0 19_1 19_2 19_3 19_4 19_5 19_6
# HT060P1-S1R1Fp1U1     1    3    1    1    2    0    0
# HT125P1-S1H4A1L1U1    0    2    1    0    1    0    0
# HT125P1-S1H8A1U1      1    1    1    0    0    0    0
# HT185P1-S1H2L1U1      0    0    3    3    0    0    0
# HT227P1-S1H1L1U1    387    9  160   85   51    2    6 
# HT242P1-S1H4L4U1     39    3  162  124    0    0    1
# based on a qualitative visual inspection of these cells overlayed on the H&E in 
# Xenium Explorer these are cells that are in close proximity to tumor cells 
# but are unlikely to be actual tumor cells themselves.
pdf(paste0(out_dir,"/","Dimplot_",subcluster_ident,"_",prefix,"_res0.3_highlight_19_0__19_2_tumor_20241112.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = "merge.sub.19", cells.highlight = sub_tumor_individual, cols.highlight = "green", sizes.highlight=0.3, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
dev.off()
for (sample_ID in sample_vector) {
    tmp_df <- obj@meta.data[(obj@meta.data[,"sample_ID"] == sample_ID),c('original_Xenium_barcode',subcluster_ident)]
    colnames(tmp_df) <- c('cell_id', 'group')
    write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",sample_ID,"_",subcluster_ident,".csv"),sep=",",quote=F,row.names=F)
    colnames(tmp_df) <- c('original_Xenium_barcode', subcluster_ident)
    write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",sample_ID,"_",subcluster_ident,".tsv"),sep="\t",quote=F)
}
# subclustering cluster 13 independently of other clusters
Idents(obj) <- "seurat_clusters"
plan(sequential)
subcluster_ident <- "merge.sub.13"
subcluster_idents <- c(13)
obj <- FindSubCluster(obj, cluster = subcluster_idents, graph.name = "SCT_snn", subcluster.name = subcluster_ident, resolution = 0.3, algorithm = 1) #based on the recommendation here https://github.com/satijalab/seurat/issues/6682 by the authors I used the SCT_snn instead of SCT_nn graph
Idents(obj) <- subcluster_ident
tmp_df <- obj@meta.data[,c('barcode',subcluster_ident)]
colnames(tmp_df) <- c('barcode',subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
# Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
# /diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
# ./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
# /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_merge.sub.13.tsv \
# merge.sub.13
# results in /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only//marker_dotplots/merge.sub.13.PDAC_merge_primary_KRAS_20241205_single_SCT.Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.pdf
plan(multicore, workers = 10)
degs <- FindAllMarkers(obj, assay='SCT', slot='data', test.use='wilcox', logfc.threshold = 0, min.pct=0.1, min.diff.pct=0.1)
plan(sequential)
write.table(degs,paste0(out_dir,"/","FindAllMarkers_",subcluster_ident,"_",prefix,"_20241212.tsv"),sep="\t",quote=FALSE)
system(paste0("Rscript /diskmnt/Projects/Users/austins2/tools/markers/Shared_marker_comparison_from_seurat_degs.R -i ",paste0(out_dir,"/","FindAllMarkers_",unified_subcluster_ident,"_",prefix,"_20241212.tsv")))
pdf(paste0(out_dir,"/","Dimplot_",subcluster_ident,"_",prefix,"_res0.3_split.by_20241212.pdf"),useDingbats = F, width=100, height=6)
print(DimPlot(obj, reduction = "umap.50PC", split.by = subcluster_ident, label=TRUE, label.size=4, raster=FALSE))
dev.off()
pdf(paste0(out_dir,"/","Dimplot_",subcluster_ident,"_",prefix,"_res0.3_group.by_20241112.pdf"),useDingbats = F, width=8, height=6)
print(DimPlot(obj, reduction = "umap.50PC", group.by = subcluster_ident, label=TRUE, label.size=4, raster=FALSE))
dev.off()
library(scCustomize)
cells_to_highlight <- list()
merge_idents_highlight <- unique(obj@meta.data[,subcluster_ident])[grepl("_",unique(obj@meta.data[,subcluster_ident]))]
for (sub_ident in merge_idents_highlight) {
    print(sub_ident)
    #string_sub_ident <- paste0(sub_ident,"_")
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,subcluster_ident] == sub_ident),"barcode"]
    cells_to_highlight[[sub_ident]] <- barcodes_in_sub
    print(summary(obj@meta.data[,"nCount_Xenium"][obj@meta.data[,subcluster_ident] == sub_ident]))
}
# [1] "13_2"                                                                                                                                           
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                      
# 1.00    7.00   13.00   23.58   22.00  404.00                                                                                                      
# [1] "13_0"                                                                                                                                           
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                      
# 2.00   26.00   52.00   68.95   92.00  424.00 
# [1] "13_12"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    5.00    9.00   12.74   14.00  133.00 
# [1] "13_3"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    8.00   14.00   24.12   25.00  342.00 
# [1] "13_1"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.00   27.00   44.00   54.79   69.00  588.00 
# [1] "13_6"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    5.00   10.00   15.28   20.00  165.00 
# [1] "13_4"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00    5.00    9.00   12.47   16.00  167.00 
# [1] "13_5"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.00   47.00   72.00   84.64  106.00  795.00 
# [1] "13_9"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00    7.00   10.00   13.06   16.00   87.00 
# [1] "13_7"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00    7.00   12.00   15.31   17.00  230.00 
# [1] "13_8"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.0     7.0    12.0    12.5    16.0    51.0 
# [1] "13_10"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00    8.00   13.00   14.48   18.00   68.00 
# [1] "13_11"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   1.382   2.000   3.000
varibow_pal <- DiscretePalette_scCustomize(num_colors = length(cells_to_highlight), palette = "varibow")
pdf(paste0(out_dir,"/","Dimplot_",subcluster_ident,"_",prefix,"_res0.3_highlight_subclusters_20241112.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = subcluster_ident, cells.highlight = cells_to_highlight, cols.highlight = varibow_pal, sizes.highlight=0.1, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
dev.off()
print(table(obj$cell_type_individual[obj$seurat_clusters %in% subcluster_idents], obj$sample_ID[obj$seurat_clusters %in% subcluster_idents]))
print(table(obj$cell_type_individual[obj@meta.data[,c(subcluster_ident)] %in% merge_idents_highlight],
            obj@meta.data[,subcluster_ident][obj@meta.data[,c(subcluster_ident)] %in% merge_idents_highlight]))
print(table(obj@meta.data[,subcluster_ident]))
#      0      1     10     11     12   13_0   13_1  13_10  13_11  13_12   13_2 
# 255560 115936  31171  27929  21689   4474   3389    140    123    118   2693 
#   13_3   13_4   13_5   13_6   13_7   13_8   13_9     14     15     16     17 
#   2656   2556   2245    547    393    313    286  17772  14709  12988  10848 
#     18     19      2     20     21     22     23     24     25     26     27 
#   8639   8456 111752   8055   7528   7071   5079   5029   4430   2052   2015 
#     28     29      3     30     31     32     33     34      4      5      6 
#   1982   1722  70375   1572   1529    495     65     19  63052  52907  47539 
#      7      8      9 
#  44723  36987  31289
for (sample_ID in sample_vector) {
    tmp_df <- obj@meta.data[(obj@meta.data[,"sample_ID"] == sample_ID),c('original_Xenium_barcode',subcluster_ident)]
    colnames(tmp_df) <- c('cell_id', 'group')
    write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",sample_ID,"_",prefix,"_",subcluster_ident,".csv"),sep=",",quote=F,row.names=F)
}
table(obj$sample_ID[obj@meta.data[,subcluster_ident] %in% merge_idents_highlight], obj@meta.data[,subcluster_ident][obj@meta.data[,subcluster_ident] %in% merge_idents_highlight])
#                     13_0 13_1 13_10 13_11 13_12 13_2 13_3 13_4 13_5 13_6 13_7
# HT060P1-S1R1Fp1U1   2030  172    90    37    74 1310 1146 1420   97  260  192
# HT061P1-S1P1A1L1U1    25    5     3     0     0   19   26    8   18   13    1
# HT061P1-S1P1A1L4U1   100   14     1     0     2  180   31   19   10   19   10
# HT125P1-S1H4A1L1U1   179   24     9     3     6  152  110   50    3   35    7
# HT125P1-S1H8A1U1     109   36     7     0     7  105   73   27   68   20    5
# HT185P1-S1H2L1U1     135   49    14    29    18  292  204  301    3   76   43
# HT227P1-S1H1L1U1     135 2858    10    26     1  387  394  374 2036   82   83
# HT242P1-S1H4L4U1    1718  224     4    27     4  179  608  348    8   27   50
# HT270P1-S1H1A1US2_1   43    7     2     1     6   69   64    9    2   15    2
# 
#                     13_8 13_9
# HT060P1-S1R1Fp1U1    252  113
# HT061P1-S1P1A1L1U1     0    0
# HT061P1-S1P1A1L4U1     1    4
# HT125P1-S1H4A1L1U1     4    4
# HT125P1-S1H8A1U1       0    1
# HT185P1-S1H2L1U1      26   42
# HT227P1-S1H1L1U1       9   42
# HT242P1-S1H4L4U1      21   80
# HT270P1-S1H1A1US2_1    0    0
print(table(obj@meta.data[,subcluster_ident][obj@meta.data[,c(subcluster_ident)] %in% merge_idents_highlight]))
# 13_0  13_1 13_10 13_11 13_12  13_2  13_3  13_4  13_5  13_6  13_7  13_8  13_9 
# 4474  3389   140   123   118  2693  2656  2556  2245   547   393   313   286
# 13_0 myCAF TNC+ ACTA2+ COL5A2+ THBS2+ FBN1+ VCAN+ RGS16+ C7- PTGDS- CXCR4-
# 13_1 Tumor_proliferative (mod_high_counts weekly proliferative (20% TOP2A+ CDK1+ MKI67+ CENPF+), MET+ KRT7+ COL17A1+
# 13_2 CXCL6+ and no other DEGs EPCAM- EHF- FXYD2- (not epithelial) mostly HT060P1. and in dense immume, fibroblast region with minimal visible tumor. cells are not over tumor cells
# 13_3 myCAF_low_count  TNC+ (highest)+ PDGFRB+ (low) PDGFRA- C7- PTGDS- THBS2+ (low) COL5A2+ spatially this is in immune dense regions.
# 13_4 LowCount COL17A1+ but that is it. mostly HT060P1. and in dense immume, fibroblast region with minimal visible tumor. cells are not over tumor cells
# 13_5 Tumor (MET+, KRT7+, COL17A1+, (low), EHF+, TFF2- LGR5- GPX2+ (sparse, low))
# 13_6 LowCount 
# 13_7 LowCount COL17A1+ but that is it. mostly HT060P1. and in dense immume, fibroblast region with minimal visible tumor. cells are not over tumor cells
# 13_8 LowCount COL17A1+ but that is it. mostly HT060P1. and in dense immume, fibroblast region with minimal visible tumor. cells are not over tumor cells
# 13_9 LowCount COL17A1+ but that is it. mostly HT060P1. and in dense immume, fibroblast region with minimal visible tumor. cells are not over tumor cells
# 13_10 LowCount COL17A1+ but that is it. mostly HT060P1. and in dense immume, fibroblast region with minimal visible tumor. cells are not over tumor cells
# 13_11 LowCount COL17A1+ but that is it. mostly HT060P1. and in dense immume, fibroblast region with minimal visible tumor. cells are not over tumor cells
# 13_12 LowCount
# given 13_3 and 13_0 are likely myCAFs based on their expression profiles it is not unlikely that the other 13_* LowCount clusters are also myCAF however, I cannot tell for certain.
# subclustering likely T_cell populations independently of other clusters
Idents(obj) <- "merge.sub.19"
plan(sequential)
subcluster_ident <- "merge.sub.2.20.19_1.19_5"
subcluster_idents <- c('2','20','19_1','19_5')
obj <- FindSubCluster(obj, cluster = subcluster_idents, graph.name = "SCT_snn", subcluster.name = subcluster_ident, resolution = 0.3, algorithm = 1) #based on the recommendation here https://github.com/satijalab/seurat/issues/6682 by the authors I used the SCT_snn instead of SCT_nn graph
cells_to_highlight <- list()
merge_idents_highlight <- c()
for (new_ident in subcluster_idents) {
    merge_idents_highlight_single <- unique(obj@meta.data[,subcluster_ident])[grepl(paste0(new_ident,"_"),unique(obj@meta.data[,subcluster_ident]))]
    merge_idents_highlight <- unique(c(merge_idents_highlight, merge_idents_highlight_single))
}
unified_subcluster_ident <- paste0("unified.",subcluster_ident)
obj@meta.data[,unified_subcluster_ident] <- obj@meta.data[,subcluster_ident]
unified_cells_to_highlight <- list()
merge_unified_subcluster_label <- c()
for (sub_ident in merge_idents_highlight) {
    print(sub_ident)
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,subcluster_ident] == sub_ident),"barcode"]
    cells_to_highlight[[sub_ident]] <- barcodes_in_sub
    new_subcluster <- unlist(strsplit(sub_ident, "_"))[length(unlist(strsplit(sub_ident, "_")))]
    new_subcluster_label <- paste0("sub_",new_subcluster)
    obj@meta.data[,unified_subcluster_ident][obj@meta.data[,subcluster_ident] == sub_ident] <- new_subcluster_label
    merge_unified_subcluster_label <- unique(c(merge_unified_subcluster_label, new_subcluster_label))
}
print(unique(obj@meta.data[,unified_subcluster_ident]))
unified_cells_to_highlight <- list()
for (unified_subcluster_label in merge_unified_subcluster_label) {
    print(unified_subcluster_label)
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,unified_subcluster_ident] == unified_subcluster_label),"barcode"]
    unified_cells_to_highlight[[unified_subcluster_label]] <- barcodes_in_sub
    print(summary(obj@meta.data[,"nCount_Xenium"][obj@meta.data[,unified_subcluster_ident] == unified_subcluster_label]))
}
# [1] "sub_2"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.0    40.0    71.0    79.5   109.0   647.0                                                                                                                     
# [1] "sub_0"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.00   42.00   73.00   82.45  112.00  526.00                                                                                                                     
# [1] "sub_3"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.00   13.00   31.00   53.09   76.00  423.00                                                                                                                     
# [1] "sub_1"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.00   37.00   67.00   72.41  100.00  382.00                                                                                                                     
# [1] "sub_5"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.00   37.00   68.00   76.83  105.00  447.00                                                                                                                     
# [1] "sub_10"                                                                                                                                                        
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 7.00   52.00   79.00   87.39  115.00  315.00                                                                                                                     
# [1] "sub_4"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 1.00   28.00   50.00   62.12   84.00  427.00
# [1] "sub_9"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 11.00   60.00   91.00   96.92  128.00  344.00 
# [1] "sub_11"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 9.0    57.0    92.0   101.5   134.0   365.0 
# [1] "sub_12"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 15      81     121     123     160     399 
# [1] "sub_8"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8.0    51.5    94.0   109.9   151.0   492.0 
# [1] "sub_7"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    7.00   13.00   20.12   24.00  227.00 
# [1] "sub_6"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.00   51.00   90.00   95.76  130.00  365.00 
# [1] "sub_13"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    5.00   10.00   18.93   21.00  257.00 
# [1] "sub_14"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   1.344   2.000   3.000
Idents(obj) <- unified_subcluster_ident
tmp_df <- obj@meta.data[,c('barcode',subcluster_ident)]
colnames(tmp_df) <- c('barcode',subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
tmp_df <- obj@meta.data[,c('barcode',unified_subcluster_ident)]
colnames(tmp_df) <- c('barcode',unified_subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",unified_subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
# Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
# /diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
# ./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
# /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_unified.merge.sub.2.20.19_1.19_5.tsv \
# unified.merge.sub.2.20.19_1.19_5
# results in /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only//marker_dotplots/unified.merge.sub.2.20.19_1.19_5.PDAC_merge_primary_KRAS_20241205_single_SCT.Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.pdf
plan(multicore, workers = 10)
degs <- FindAllMarkers(obj, assay='SCT', slot='data', test.use='wilcox', logfc.threshold = 0, min.pct=0.1, min.diff.pct=0.1)
plan(sequential)
write.table(degs,paste0(out_dir,"/","FindAllMarkers_",unified_subcluster_ident,"_",prefix,"_20241212.tsv"),sep="\t",quote=FALSE)
system(paste0("Rscript /diskmnt/Projects/Users/austins2/tools/markers/Shared_marker_comparison_from_seurat_degs.R -i ",paste0(out_dir,"/","FindAllMarkers_",unified_subcluster_ident,"_",prefix,"_20241212.tsv")))
pdf(paste0(out_dir,"/","Dimplot_",unified_subcluster_ident,"_",prefix,"_res0.3_split.by_20241212.pdf"),useDingbats = F, width=100, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", split.by = unified_subcluster_ident, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
dev.off()
pdf(paste0(out_dir,"/","Dimplot_",unified_subcluster_ident,"_",prefix,"_res0.3_group.by_20241112.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
dev.off()
pdf(paste0(out_dir,"/","Dimplot_",subcluster_ident,"_",prefix,"_res0.3_split.by_20241212.pdf"),useDingbats = F, width=100, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", split.by = subcluster_ident, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
dev.off()
pdf(paste0(out_dir,"/","Dimplot_",subcluster_ident,"_",prefix,"_res0.3_group.by_20241112.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = subcluster_ident, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
dev.off()
library(scCustomize)
varibow_pal <- DiscretePalette_scCustomize(num_colors = length(unified_cells_to_highlight), palette = "varibow")
ordered_labels_for_plotting <- c()
for (i in 0:(length(merge_unified_subcluster_label)-1)) {
    unified_subcluster_label <- paste0("sub_",i)
    ordered_labels_for_plotting <- c(ordered_labels_for_plotting, unified_subcluster_label)
}
unified_cells_to_highlight <- unified_cells_to_highlight[ordered_labels_for_plotting]
pdf(paste0(out_dir,"/","Dimplot_",unified_subcluster_ident,"_",prefix,"_res0.3_highlight_subclusters_20241212.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, cells.highlight = unified_cells_to_highlight, cols.highlight = varibow_pal, sizes.highlight=0.1, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
for (i in 1:length(ordered_labels_for_plotting)) {
    unified_subcluster_label <- ordered_labels_for_plotting[i]
    print(unified_subcluster_label)
    unified_cells_to_highlight_sep <- unified_cells_to_highlight[[unified_subcluster_label]]
    col_highlight <- varibow_pal[i]
    print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, cells.highlight = unified_cells_to_highlight_sep, cols.highlight = col_highlight, sizes.highlight=0.1, label=TRUE, label.size=4, raster=FALSE) + labs(title = paste0("highlight ",unified_subcluster_label)), layers='Point', dpi=300))
}
dev.off()
print(table(obj$cell_type_individual[obj$merge.sub.19 %in% subcluster_idents], obj$sample_ID[obj$merge.sub.19 %in% subcluster_idents]))
print(table(obj$cell_type_individual[obj@meta.data[,c(unified_subcluster_ident)] %in% merge_unified_subcluster_label],
            obj@meta.data[,unified_subcluster_ident][obj@meta.data[,c(unified_subcluster_ident)] %in% merge_unified_subcluster_label]))
print(table(obj@meta.data[,unified_subcluster_ident]))
#      0      1     10     11     12     13     14     15     16     17     18 
# 255560 115936  31171  27929  21689  19933  17772  14709  12988  10848   8639 
#   19_0   19_2   19_3   19_4   19_6   19_7     21     22     23     24     25 
#   2758   1262   1258    591    362    216   7528   7071   5079   5029   4430 
#     26     27     28     29      3     30     31     32     33     34      4 
#   2052   2015   1982   1722  70375   1572   1529    495     65     19  63052 
#      5      6      7      8      9  sub_0  sub_1 sub_10 sub_11 sub_12 sub_13 
#  52907  47539  44723  36987  31289  35033  21409   1140    842    341    121 
# sub_14  sub_2  sub_3  sub_4  sub_5  sub_6  sub_7  sub_8  sub_9 
#     96  19719  17668   9159   7024   3488   2301   2055   1420 
for (sample_ID in sample_vector) {
    tmp_df <- obj@meta.data[(obj@meta.data[,"sample_ID"] == sample_ID),c('original_Xenium_barcode',unified_subcluster_ident)]
    colnames(tmp_df) <- c('cell_id', 'group')
    write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",sample_ID,"_",prefix,"_",unified_subcluster_ident,".csv"),sep=",",quote=F,row.names=F)
}
table(obj$sample_ID[obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting], obj@meta.data[,unified_subcluster_ident][obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting])
#                     sub_0 sub_1 sub_10 sub_11 sub_12 sub_13 sub_14 sub_2
# HT060P1-S1R1Fp1U1    8647  7641    520    358    135     75     35  9624
# HT061P1-S1P1A1L1U1   2205   734    131     88     19      6      0  1034
# HT061P1-S1P1A1L4U1   5353  1157     77     79    103      3      0  1431
# HT125P1-S1H4A1L1U1   5117  2424    222    103     62      8      2  3182
# HT125P1-S1H8A1U1     3385  1833     78     76      0      3      0  1466
# HT185P1-S1H2L1U1      927  1093      9     14      8      8     20   353
# HT227P1-S1H1L1U1      951  1830     12     13     12      8     29   866
# HT242P1-S1H4L4U1     4963  2213     49     60      0      6      9  1410
# HT270P1-S1H1A1US2_1  3485  2484     42     51      2      4      1   353
#                     sub_3 sub_4 sub_5 sub_6 sub_7 sub_8 sub_9
# HT060P1-S1R1Fp1U1    8143  1607  3398   919   670   796   371
# HT061P1-S1P1A1L1U1    901   369   162   480   107   254   514
# HT061P1-S1P1A1L4U1   1287   561   209   949    51   142   102
# HT125P1-S1H4A1L1U1   1344   660  1055   531   106   164   194
# HT125P1-S1H8A1U1      905   272   708   149    77   113    35
# HT185P1-S1H2L1U1     1494   203   496    36   200    35     3
# HT227P1-S1H1L1U1     1612   403   398    36   199   108     7
# HT242P1-S1H4L4U1      712  4868   425   357   832   408   127
# HT270P1-S1H1A1US2_1  1270   216   173    31    59    35    67
print(table(obj@meta.data[,unified_subcluster_ident][obj@meta.data[,c(unified_subcluster_ident)] %in% ordered_labels_for_plotting]))
# sub_0  sub_1 sub_10 sub_11 sub_12 sub_13 sub_14  sub_2  sub_3  sub_4  sub_5 
# 35033  21409   1140    842    341    121     96  19719  17668   9159   7024 
# sub_6  sub_7  sub_8  sub_9 
# 3 488   2301   2055   1420
# sub_0 CD8+T_cell CCL5+ TRAC+ CXCR4+ CD2+ PTPRC+ CD3E+ CD8A+ CD3D+ GZMA+ GZMK+
# sub_1 CD4+T_cell IL7R+ TRAC+ CXCR4+ CD2+ PTPRC+ CD3D+ CD4+ SELL+
# sub_2 Tregs FOXP3+ CTLA4+ strongest here
# sub_3 CD4+T_cell
# sub_4 NK_cell CD3E+ CD3D- GNLY+ (specific) GZMA+ GZMB+ (specific -only other expression is sub_12) CD247+ FCGR3A+
# sub_5 CD4+T_cell
# sub_6 CD4+T_cell maybe some myeloid (this is much weaker than T cell signals)
# sub_7 CD8+T_cell cD8A+ CD4- 
# sub_8 T_cell_proliferative CD4+ CD8A+ CD3D+ CD3E+ 
# sub_9 Naive_T_cell CD27+ (highest), CD28+ (lower) CD3E+ CD3D+ CD79a+ BANK1- MS4A1- Naive T cells express CD27 and CD28 which decreases as T cells mature. mixed_Plasma
# sub_10 CD4+T_cell CPA3+ KIT+ TRAC+ CD2+ CD3E+ CD274-. This paper indicates that CD117+ (gene: KIT+) T cells are early an T cell progenitor stages: https://stemcellres.biomedcentral.com/articles/10.1186/s13287-017-0495-4. Naive T cells express CD27 and CD28 which decreases as T cells mature.
# sub_11 CD4+T_cell some small endothelial markers (VWF, PECAM1, CD34, STC1) extravasation being picked up by imperfect segmentation?
# sub_12 CD8+T_cell Epi markers here (EHF+, EPCAM+ - T_cells close to epithelial structures?). Also plenty of CD8+
# sub_13 LowCount the markers here indicate it could be a small populations of Erythrocytes A(HSP+ ALAS2+ GYPA+ CD3D- CD3E- (not T_cells)), however, there are no visible red blood cells in the histology of the sample overlapping these cells. Also mature erythrocytes don't have nuclei?
# sub_14 LowCount

# subclustering cluster Subcluster Myeloid and split out independently of other clusters
Idents(obj) <- "merge.sub.19"
plan(sequential)
subcluster_ident <- "merge.sub.1.14.25.19_0.19_4.19_6"
subcluster_idents <- c("1","14","25",'19_0','19_4','19_6')
obj <- FindSubCluster(obj, cluster = subcluster_idents, graph.name = "SCT_snn", subcluster.name = subcluster_ident, resolution = 0.3, algorithm = 1) #based on the recommendation here https://github.com/satijalab/seurat/issues/6682 by the authors I used the SCT_snn instead of SCT_nn graph
cells_to_highlight <- list()
merge_idents_highlight <- c()
for (new_ident in subcluster_idents) {
    merge_idents_highlight_single <- unique(obj@meta.data[,subcluster_ident])[grepl(paste0(new_ident,"_"),unique(obj@meta.data[,subcluster_ident]))]
    merge_idents_highlight <- unique(c(merge_idents_highlight, merge_idents_highlight_single))
}
unified_subcluster_ident <- paste0("unified.",subcluster_ident)
obj@meta.data[,unified_subcluster_ident] <- obj@meta.data[,subcluster_ident]
unified_cells_to_highlight <- list()
merge_unified_subcluster_label <- c()
for (sub_ident in merge_idents_highlight) {
    print(sub_ident)
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,subcluster_ident] == sub_ident),"barcode"]
    cells_to_highlight[[sub_ident]] <- barcodes_in_sub
    new_subcluster <- unlist(strsplit(sub_ident, "_"))[length(unlist(strsplit(sub_ident, "_")))]
    new_subcluster_label <- paste0("sub_",new_subcluster)
    obj@meta.data[,unified_subcluster_ident][obj@meta.data[,subcluster_ident] == sub_ident] <- new_subcluster_label
    merge_unified_subcluster_label <- unique(c(merge_unified_subcluster_label, new_subcluster_label))
}
print(unique(obj@meta.data[,unified_subcluster_ident]))
unified_cells_to_highlight <- list()
for (unified_subcluster_label in merge_unified_subcluster_label) {
    print(unified_subcluster_label)
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,unified_subcluster_ident] == unified_subcluster_label),"barcode"]
    unified_cells_to_highlight[[unified_subcluster_label]] <- barcodes_in_sub
    print(summary(obj@meta.data[,"nCount_Xenium"][obj@meta.data[,unified_subcluster_ident] == unified_subcluster_label]))
}
# [1] "sub_0"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 3.00   40.00   73.00   87.65  120.00  579.00                                                                                                                     
# [1] "sub_12"                                                                                                                                                        
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.00   34.00   66.00   76.95  106.00  364.00                                                                                                                     
# [1] "sub_2"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.00   21.00   48.00   68.72   96.00  556.00                                                                                                                     
# [1] "sub_1"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 1.0    45.0    86.0   100.1   139.0   663.0                                                                                                                     
# [1] "sub_8"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.0    24.0    51.0    67.1    91.0   886.0                                                                                                                     
# [1] "sub_7"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 1.00   29.00   67.00   85.57  122.00  582.00                                                                                                                     
# [1] "sub_3"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 7.00   48.00   82.00   96.17  129.00  525.00
# [1] "sub_5"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.00    7.00   14.00   24.63   28.00  420.00                                                                                                                     
# [1] "sub_4"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 1.00   33.00   65.00   83.82  117.00  765.00                                                                                                                     
# [1] "sub_9"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 5.0    38.0    80.0   101.7   143.0   533.0                                                                                                                     
# [1] "sub_6"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 1.00    4.00    7.00   12.34   14.00  266.00                                                                                                                     
# [1] "sub_11"                                                                                                                                                        
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 1.00    6.00   11.00   14.81   17.00  250.00                                                                                                                     
# [1] "sub_10"                                                                                                                                                        
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.00   10.75   21.00   30.91   40.00  335.00                                                                                                                     
# [1] "sub_14"                                                                                                                                                        
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 1.00    5.00   10.00   13.59   17.00  101.00                                                                                                                     
# [1] "sub_13"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    5.00    9.00   13.91   13.00  278.00 
# [1] "sub_20"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   1.258   1.000   3.000 
# [1] "sub_16"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   2.000   1.852   2.000   8.000 
# [1] "sub_15"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   1.274   2.000   3.000 
# [1] "sub_19"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   2.000   1.939   2.000  13.000 
# [1] "sub_17"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   2.000   2.585   3.000  10.000 
# [1] "sub_18"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0     1.0     1.0     1.5     2.0     4.0 
# [1] "sub_21"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 96.00   98.25  100.50  100.50  102.75  105.00
Idents(obj) <- unified_subcluster_ident
tmp_df <- obj@meta.data[,c('barcode',subcluster_ident)]
colnames(tmp_df) <- c('barcode',subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
tmp_df <- obj@meta.data[,c('barcode',unified_subcluster_ident)]
colnames(tmp_df) <- c('barcode',unified_subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",unified_subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
# Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
# /diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
# ./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
# /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_unified.merge.sub.1.14.25.19_0.19_4.19_6.tsv \
# unified.merge.sub.1.14.25.19_0.19_4.19_6
# results in /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only//marker_dotplots/unified.merge.sub.1.14.25.19_0.19_4.19_6.PDAC_merge_primary_KRAS_20241205_single_SCT.Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.pdf
plan(multicore, workers = 10)
degs <- FindAllMarkers(obj, assay='SCT', slot='data', test.use='wilcox', logfc.threshold = 0, min.pct=0.1, min.diff.pct=0.1)
plan(sequential)
write.table(degs,paste0(out_dir,"/","FindAllMarkers_",unified_subcluster_ident,"_",prefix,"_20241212.tsv"),sep="\t",quote=FALSE)
system(paste0("Rscript /diskmnt/Projects/Users/austins2/tools/markers/Shared_marker_comparison_from_seurat_degs.R -i ",paste0(out_dir,"/","FindAllMarkers_",unified_subcluster_ident,"_",prefix,"_20241212.tsv")))
pdf(paste0(out_dir,"/","Dimplot_",unified_subcluster_ident,"_",prefix,"_res0.3_split.by_20241212.pdf"),useDingbats = F, width=100, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", split.by = unified_subcluster_ident, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
dev.off()
pdf(paste0(out_dir,"/","Dimplot_",unified_subcluster_ident,"_",prefix,"_res0.3_group.by_20241112.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
dev.off()
pdf(paste0(out_dir,"/","Dimplot_",subcluster_ident,"_",prefix,"_res0.3_split.by_20241212.pdf"),useDingbats = F, width=100, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", split.by = subcluster_ident, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
dev.off()
pdf(paste0(out_dir,"/","Dimplot_",subcluster_ident,"_",prefix,"_res0.3_group.by_20241112.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = subcluster_ident, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
dev.off()
library(scCustomize)
varibow_pal <- DiscretePalette_scCustomize(num_colors = length(unified_cells_to_highlight), palette = "varibow")
ordered_labels_for_plotting <- c()
for (i in 0:(length(merge_unified_subcluster_label)-1)) {
    unified_subcluster_label <- paste0("sub_",i)
    ordered_labels_for_plotting <- c(ordered_labels_for_plotting, unified_subcluster_label)
}
unified_cells_to_highlight <- unified_cells_to_highlight[ordered_labels_for_plotting]
pdf(paste0(out_dir,"/","Dimplot_",unified_subcluster_ident,"_",prefix,"_res0.3_highlight_subclusters_20241212.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, cells.highlight = unified_cells_to_highlight, cols.highlight = varibow_pal, sizes.highlight=0.1, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
for (i in 1:length(ordered_labels_for_plotting)) {
    unified_subcluster_label <- ordered_labels_for_plotting[i]
    print(unified_subcluster_label)
    unified_cells_to_highlight_sep <- unified_cells_to_highlight[[unified_subcluster_label]]
    col_highlight <- varibow_pal[i]
    print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, cells.highlight = unified_cells_to_highlight_sep, cols.highlight = col_highlight, sizes.highlight=0.1, label=TRUE, label.size=4, raster=FALSE) + labs(title = paste0("highlight ",unified_subcluster_label)), layers='Point', dpi=300))
}
dev.off()
print(table(obj$cell_type_individual[obj$merge.sub.19 %in% subcluster_idents], obj$sample_ID[obj$merge.sub.19 %in% subcluster_idents]))
print(table(obj$cell_type_individual[obj@meta.data[,c(unified_subcluster_ident)] %in% merge_unified_subcluster_label],
            obj@meta.data[,unified_subcluster_ident][obj@meta.data[,c(unified_subcluster_ident)] %in% merge_unified_subcluster_label]))
print(table(obj@meta.data[,unified_subcluster_ident]))
#      0     10     11     12     13     15     16     17     18   19_1   19_2 
# 255560  31171  27929  21689  19933  14709  12988  10848   8639   1429   1262 
#   19_3   19_5   19_7      2     20     21     22     23     24     26     27 
#   1258    580    216 111752   8055   7528   7071   5079   5029   2052   2015 
#     28     29      3     30     31     32     33     34      4      5      6 
#   1982   1722  70375   1572   1529    495     65     19  63052  52907  47539 
#      7      8      9  sub_0  sub_1 sub_10 sub_11 sub_12 sub_13 sub_14 sub_15 
#  44723  36987  31289  39908  16293   3432   2711   1184   1029    741    168 
# sub_16 sub_17 sub_18 sub_19  sub_2 sub_20 sub_21  sub_3  sub_4  sub_5  sub_6 
#    128     94     68     66  15452     62      2  14762   9519   9221   8628 
#  sub_7  sub_8  sub_9 
#   8275   6352   3754
for (sample_ID in sample_vector) {
    tmp_df <- obj@meta.data[(obj@meta.data[,"sample_ID"] == sample_ID),c('original_Xenium_barcode',unified_subcluster_ident)]
    colnames(tmp_df) <- c('cell_id', 'group')
    write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",sample_ID,"_",prefix,"_",unified_subcluster_ident,".csv"),sep=",",quote=F,row.names=F)
}
table(obj$sample_ID[obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting], obj@meta.data[,unified_subcluster_ident][obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting])
#                     sub_0 sub_1 sub_10 sub_11 sub_12 sub_13 sub_14 sub_15                                                                                         
# HT060P1-S1R1Fp1U1   12155 11268    595   1351    624    378    436     30                                                                                         
# HT061P1-S1P1A1L1U1   2391   301     95     39    125     20     15      1                                                                                         
# HT061P1-S1P1A1L4U1   9553  1192    199     84     84     90      7      4                                                                                         
# HT125P1-S1H4A1L1U1   2909   911     81     95    127     41     12      6                                                                                         
# HT125P1-S1H8A1U1     2307  1253    142     57     56     21      8      5                                                                                         
# HT185P1-S1H2L1U1     1499   349    336    300     34    178     89     48                                                                                         
# HT227P1-S1H1L1U1     3129   386    671    253     40    152     34     34                                                                                         
# HT242P1-S1H4L4U1     3137   460   1173    492     43    103    121     36                                                                                         
# HT270P1-S1H1A1US2_1  2828   173    140     40     51     46     19      4
#                     sub_16 sub_17 sub_18 sub_19 sub_2 sub_20 sub_21 sub_3
# HT060P1-S1R1Fp1U1       31     47     16     21  4753     19      0  5242
# HT061P1-S1P1A1L1U1       2      0      1      3   808      0      0   609
# HT061P1-S1P1A1L4U1       1      1      4      0  4422      0      0   769
# HT125P1-S1H4A1L1U1       5      4      4      3   783      1      0  2657
# HT125P1-S1H8A1U1         2      0      1      5  1039      1      0  1889
# HT185P1-S1H2L1U1        62     22     18     16  1478     23      0   946
# HT227P1-S1H1L1U1        16      3      7      7   956     12      0  1421
# HT242P1-S1H4L4U1         6     16     17     10   374      6      1   605
# HT270P1-S1H1A1US2_1      3      1      0      1   839      0      1   624
#                     sub_4 sub_5 sub_6 sub_7 sub_8 sub_9
# HT060P1-S1R1Fp1U1    3149  2016  3236  2796  3428  1524
# HT061P1-S1P1A1L1U1    875   232   266   853   196   129
# HT061P1-S1P1A1L4U1   1494  2398   323  2009   246   427
# HT125P1-S1H4A1L1U1    710   303   268   732   360   293
# HT125P1-S1H8A1U1      275   594   385   250   278   236
# HT185P1-S1H2L1U1      177   937  1492   268   502   140
# HT227P1-S1H1L1U1      516  1697  1659   269   613   641
# HT242P1-S1H4L4U1     2100   737   735   990   189   281
# HT270P1-S1H1A1US2_1   223   307   264   108   540    83
print(table(obj@meta.data[,unified_subcluster_ident][obj@meta.data[,c(unified_subcluster_ident)] %in% ordered_labels_for_plotting]))
# sub_0  sub_1 sub_10 sub_11 sub_12 sub_13 sub_14 sub_15 sub_16 sub_17 sub_18 
# 39908  16293   3432   2711   1184   1029    741    168    128     94     68 
# sub_19  sub_2 sub_20 sub_21  sub_3  sub_4  sub_5  sub_6  sub_7  sub_8  sub_9 
# 66  15452     62      2  14762   9519   9221   8628   8275   6352   3754
# Pretty much all of the macrophages here have an M2 phenotype.
# this is consistent with what is reported in this article:
# https://doi.org/10.1038/s41392-021-00769-z
# and this review:
# https://doi.org/10.1038/s41698-024-00522-z
# If there are different macrophage M1 vs M2 populations I'm not able to distinguish them in the low gene-space
# sub_0 Macrophage_M2
# sub_1 Macrophage_M2
# sub_2 Macrophage_M2
# sub_3 dendritic_cell CSF2RA+ CD1C+ IL1R2+ SERPINB9+ IRF8+ (higher in dendritic cells)
# sub_4 Macrophage_M2
# sub_5 Macrophage_M2
# sub_6 Macrophage_M2
# sub_7 Macrophage_M2
# sub_8 Macrophage_M2 subcluster this further there may fibroblasts mixed in here.
# sub_9 Myeloid_proliferative
# sub_10 Macrophage_M2
# sub_11 Macrophage_M2
# sub_12 Macrophage_M2 MS4A6A+ CD163+ CD68+ MS4A4A+ AIF1+ HPGDS+ KIT+ 
# sub_13 LowCount  GATM+ is literally the only marker here in the DEGs and on the DotPlot (clusters alongside other acinar cells but no other acinar markers expressed besides GATM)
# sub_14 LowCount RGS16+ is literally the only positive DEG and it is non-specific
# sub_15 LowCount
# sub_16 LowCount
# sub_17 LowCount
# sub_18 LowCount
# sub_19 LowCount
# sub_20 LowCount
# sub_21 iCAF FBN1+ PTGDS+ C7+ SFRP4+ SFRP2+ FBLN1+ VCAN+ bro its literally 2 cells how.

# subclustering cluste Fibroblasts and split out independently of other clusters
# setting up the input column that will be used for subclustering. 
# I am using the output of the Myeloid subclustering column but compressing all of the Myeloid subclusters into a single label except for sub_8
obj$merge.sub.1.14.25.19_0.19_4.19_6_merge.sub.13 <- obj$unified.merge.sub.1.14.25.19_0.19_4.19_6
sub.13.labels <- unique(obj$merge.sub.13)[grepl("13_",unique(obj$merge.sub.13))]
for (sub.13.label in sub.13.labels) {
    obj$merge.sub.1.14.25.19_0.19_4.19_6_merge.sub.13[obj$merge.sub.13 == sub.13.label] = sub.13.label
}
for (ordered_label in ordered_labels_for_plotting){
    obj$merge.sub.1.14.25.19_0.19_4.19_6_merge.sub.13[obj$unified.merge.sub.1.14.25.19_0.19_4.19_6 == ordered_label] <- ordered_label
}
obj$merge.sub.1.14.25.19_0.19_4.19_6_merge.sub.13[(obj$merge.sub.1.14.25.19_0.19_4.19_6_merge.sub.13 %in% ordered_labels_for_plotting) & !(obj$merge.sub.1.14.25.19_0.19_4.19_6_merge.sub.13 == "sub_8")] <- "Myeloid"
obj$merge.sub.1.14.25.19_0.19_4.19_6_merge.sub.13[(obj$merge.sub.1.14.25.19_0.19_4.19_6_merge.sub.13 %in% ordered_labels_for_plotting) & (obj$merge.sub.1.14.25.19_0.19_4.19_6_merge.sub.13 == "sub_8")] <- "Msub_8"
# doing the subclustering of Fibroblasts and split out for interrogation
Idents(obj) <- "merge.sub.1.14.25.19_0.19_4.19_6_merge.sub.13"
plan(sequential)
subcluster_ident <- "merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8"
subcluster_idents <- c('0','21','19_2','19_7','13_3','13_0','Msub_8')
obj <- FindSubCluster(obj, cluster = subcluster_idents, graph.name = "SCT_snn", subcluster.name = subcluster_ident, resolution = 0.3, algorithm = 1) #based on the recommendation here https://github.com/satijalab/seurat/issues/6682 by the authors I used the SCT_snn instead of SCT_nn graph
cells_to_highlight <- list()
merge_idents_highlight <- c()
for (new_ident in subcluster_idents) {
    merge_idents_highlight_single <- unique(obj@meta.data[,subcluster_ident])[grepl(paste0(new_ident,"_"),unique(obj@meta.data[,subcluster_ident]))]
    merge_idents_highlight <- unique(c(merge_idents_highlight, merge_idents_highlight_single))
}
unified_subcluster_ident <- paste0("unified.",subcluster_ident)
obj@meta.data[,unified_subcluster_ident] <- obj@meta.data[,subcluster_ident]
unified_cells_to_highlight <- list()
merge_unified_subcluster_label <- c()
for (sub_ident in merge_idents_highlight) {
    print(sub_ident)
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,subcluster_ident] == sub_ident),"barcode"]
    cells_to_highlight[[sub_ident]] <- barcodes_in_sub
    new_subcluster <- unlist(strsplit(sub_ident, "_"))[length(unlist(strsplit(sub_ident, "_")))]
    new_subcluster_label <- paste0("sub_",new_subcluster)
    obj@meta.data[,unified_subcluster_ident][obj@meta.data[,subcluster_ident] == sub_ident] <- new_subcluster_label
    merge_unified_subcluster_label <- unique(c(merge_unified_subcluster_label, new_subcluster_label))
}
print(unique(obj@meta.data[,unified_subcluster_ident]))
unified_cells_to_highlight <- list()
for (unified_subcluster_label in merge_unified_subcluster_label) {
    print(unified_subcluster_label)
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,unified_subcluster_ident] == unified_subcluster_label),"barcode"]
    unified_cells_to_highlight[[unified_subcluster_label]] <- barcodes_in_sub
    print(summary(obj@meta.data[,"nCount_Xenium"][obj@meta.data[,unified_subcluster_ident] == unified_subcluster_label]))
}
# [1] "sub_0"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.0    60.0    94.0   107.3   140.0   730.0                                                                                                                     
# [1] "sub_2"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.0    40.0    79.0   100.9   136.0   883.0                                                                                                                     
# [1] "sub_4"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.00   21.00   44.00   65.38   88.00  672.00                                                                                                                     
# [1] "sub_5"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 1.00   17.00   43.00   60.88   86.00  886.00
# [1] "sub_6"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 1      46      86      94     130     514                                                                                                                     
# [1] "sub_10"                                                                                                                                                        
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 19.0    75.0   110.0   118.6   148.2   373.0                                                                                                                     
# [1] "sub_8"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 1.00   15.00   32.00   52.37   70.00  424.00                                                                                                                     
# [1] "sub_9"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.00   27.00   62.50   88.24  127.75  609.00                                                                                                                     
# [1] "sub_7"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.00   32.00   61.00   77.57  106.00  778.00 
# [1] "sub_1"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.0    54.0    94.0   107.4   144.0   721.0 
# [1] "sub_3"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00   21.00   48.00   66.85   95.00  604.00
Idents(obj) <- unified_subcluster_ident
tmp_df <- obj@meta.data[,c('barcode',subcluster_ident)]
colnames(tmp_df) <- c('barcode',subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
tmp_df <- obj@meta.data[,c('barcode',unified_subcluster_ident)]
colnames(tmp_df) <- c('barcode',unified_subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",unified_subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
# Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
# /diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
# ./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
# /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8.tsv \
# unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8
# results in /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only//marker_dotplots/unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8.PDAC_merge_primary_KRAS_20241205_single_SCT.Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.pdf
plan(multicore, workers = 10)
degs <- FindAllMarkers(obj, assay='SCT', slot='data', test.use='wilcox', logfc.threshold = 0, min.pct=0.1, min.diff.pct=0.1)
plan(sequential)
write.table(degs,paste0(out_dir,"/","FindAllMarkers_",unified_subcluster_ident,"_",prefix,"_20241212.tsv"),sep="\t",quote=FALSE)
pdf(paste0(out_dir,"/","Dimplot_",unified_subcluster_ident,"_",prefix,"_res0.3_split.by_20241212.pdf"),useDingbats = F, width=100, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", split.by = unified_subcluster_ident, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
dev.off()
pdf(paste0(out_dir,"/","Dimplot_",unified_subcluster_ident,"_",prefix,"_res0.3_group.by_20241112.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
dev.off()
pdf(paste0(out_dir,"/","Dimplot_",subcluster_ident,"_",prefix,"_res0.3_split.by_20241212.pdf"),useDingbats = F, width=100, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", split.by = subcluster_ident, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
dev.off()
pdf(paste0(out_dir,"/","Dimplot_",subcluster_ident,"_",prefix,"_res0.3_group.by_20241112.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = subcluster_ident, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
dev.off()
library(scCustomize)
varibow_pal <- DiscretePalette_scCustomize(num_colors = length(unified_cells_to_highlight), palette = "varibow")
ordered_labels_for_plotting <- c()
for (i in 0:(length(merge_unified_subcluster_label)-1)) {
    unified_subcluster_label <- paste0("sub_",i)
    ordered_labels_for_plotting <- c(ordered_labels_for_plotting, unified_subcluster_label)
}
unified_cells_to_highlight <- unified_cells_to_highlight[ordered_labels_for_plotting]
pdf(paste0(out_dir,"/","Dimplot_",unified_subcluster_ident,"_",prefix,"_res0.3_highlight_subclusters_20241212.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, cells.highlight = unified_cells_to_highlight, cols.highlight = varibow_pal, sizes.highlight=0.1, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
for (i in 1:length(ordered_labels_for_plotting)) {
    unified_subcluster_label <- ordered_labels_for_plotting[i]
    print(unified_subcluster_label)
    unified_cells_to_highlight_sep <- unified_cells_to_highlight[[unified_subcluster_label]]
    col_highlight <- varibow_pal[i]
    print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, cells.highlight = unified_cells_to_highlight_sep, cols.highlight = col_highlight, sizes.highlight=0.1, label=TRUE, label.size=4, raster=FALSE) + labs(title = paste0("highlight ",unified_subcluster_label)), layers='Point', dpi=300))
}
dev.off()
print(table(obj$cell_type_individual[obj$merge.sub.19 %in% subcluster_idents], obj$sample_ID[obj$merge.sub.19 %in% subcluster_idents]))
print(table(obj$cell_type_individual[obj@meta.data[,c(unified_subcluster_ident)] %in% merge_unified_subcluster_label],
            obj@meta.data[,unified_subcluster_ident][obj@meta.data[,c(unified_subcluster_ident)] %in% merge_unified_subcluster_label]))
print(table(obj@meta.data[,unified_subcluster_ident]))
#    10      11      12    13_1   13_10   13_11   13_12    13_2    13_4    13_5
# 31171   27929   21689    3389     140     123     118    2693    2556    2245
#  13_6    13_7    13_8    13_9      15      16      17      18    19_1    19_3
#   547     393     313     286   14709   12988   10848    8639    1429    1258
#  19_5       2      20      22      23      24      26      27      28      29
#   580  111752    8055    7071    5079    5029    2052    2015    1982    1722
#     3      30      31      32      33      34       4       5       6       7
# 70375    1572    1529     495      65      19   63052   52907   47539   44723
#     8       9 Myeloid   sub_0   sub_1  sub_10   sub_2   sub_3   sub_4   sub_5
# 36987   31289  135497   91050   53492     872   49943   24627   17776   14515
# sub_6   sub_7   sub_8   sub_9
#  8570    8547    7138    1518
for (sample_ID in sample_vector) {
    tmp_df <- obj@meta.data[(obj@meta.data[,"sample_ID"] == sample_ID),c('original_Xenium_barcode',unified_subcluster_ident)]
    colnames(tmp_df) <- c('cell_id', 'group')
    write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",sample_ID,"_",prefix,"_",unified_subcluster_ident,".csv"),sep=",",quote=F,row.names=F)
}
table(obj$sample_ID[obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting], obj@meta.data[,unified_subcluster_ident][obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting])
#                     sub_0 sub_1 sub_10 sub_2 sub_3 sub_4 sub_5 sub_6 sub_7
# HT060P1-S1R1Fp1U1   58440 35577    689 12220  7090  5797  7790  5409  6386
# HT061P1-S1P1A1L1U1   2255   835      7  3428   717  1016   384   219    85
# HT061P1-S1P1A1L4U1   1908  2879     19  8308  1671  1097   540   190   223
# HT125P1-S1H4A1L1U1  10176  7990     53  6521  5049  1475   971   852    90
# HT125P1-S1H8A1U1     1107   755      9  4755  1058   351   450   451    42
# HT185P1-S1H2L1U1     6005  2408     24  3666  2224  5094  1242   375   412
# HT227P1-S1H1L1U1     2590  1007     10  2885  4058  1211  1695   372   239
# HT242P1-S1H4L4U1     1737   778      7  4624  1540   576   366   420   815
# HT270P1-S1H1A1US2_1  6832  1263     54  3536  1220  1159  1077   282   255
#                     sub_8 sub_9
# HT060P1-S1R1Fp1U1    3178   504
# HT061P1-S1P1A1L1U1     55    37
# HT061P1-S1P1A1L4U1    135   171
# HT125P1-S1H4A1L1U1    293   102
# HT125P1-S1H8A1U1      180   121
# HT185P1-S1H2L1U1      337   105
# HT227P1-S1H1L1U1      529   193
# HT242P1-S1H4L4U1     2323   246
# HT270P1-S1H1A1US2_1   108    39
# dotplot calls of subclusters
# sub_0 iCAF
# sub_1 myCAF
# sub_2 myCAF
# sub_3 myCAF
# sub_4 mix CAF
# sub_5 iCAF_myeloid_mix??? this also has macrophages from the input Msub_8 mixed in here. so this cluster needs to be split in.
# sub_6 iCAF
# sub_7 myCAF
# sub_8 myCAF 
# sub_9 myCAF
# sub_10 iCAF
# Trying out how using gene modules to assign cell types might alter cell classification.
iCAF_markers <- list(c('PTGDS','C7','SFRP4','SFRP2','VCAN','PDGFRA','PDGFRB'))
myCAF_markers <- list(c('THBS2','PDGFRB','COL5A2','SFRP4','VCAN','ASPN','ACTA2','TNC','MFAP5'))
obj <- AddModuleScore(obj, features = iCAF_markers, assay="SCT", name="iCAF_score", ctrl = 15, seed=1234)
obj <- AddModuleScore(obj, features = myCAF_markers, assay="SCT", name="myCAF_score", ctrl = 15, seed=1234)
obj$sub_4_iCAF_score <- 0
obj$sub_4_iCAF_score[obj@meta.data[,unified_subcluster_ident] == "sub_4"] <- obj$iCAF_score1[obj@meta.data[,unified_subcluster_ident] == "sub_4"]
obj$sub_4_myCAF_score <- 0
obj$sub_4_myCAF_score[obj@meta.data[,unified_subcluster_ident] == "sub_4"] <- obj$myCAF_score1[obj@meta.data[,unified_subcluster_ident] == "sub_4"]
scores_to_plot <- c("sub_4_iCAF_score","sub_4_myCAF_score")
pdf(paste0(out_dir,"/","Dimplot_",unified_subcluster_ident,"_",prefix,"_res0.3_sub_4_CAF_scores_20241212.pdf"),useDingbats = F, width=8, height=6)
for (score in scores_to_plot) {
    print(rasterize(FeaturePlot(obj, reduction = "umap.50PC", features = score, label=T, label.size=1, raster=FALSE,order = T) + scale_colour_gradient2(low = "blue", mid = "lightgrey", high = "red"), layers='Point', dpi=300))
}
dev.off()
obj$sub_4_CAF_call <- NA
obj$sub_4_CAF_call[(obj@meta.data[,unified_subcluster_ident] == "sub_4") & (obj@meta.data[,"sub_4_iCAF_score"] >= obj@meta.data[,"sub_4_myCAF_score"])] <- "iCAF"
obj$sub_4_CAF_call[(obj@meta.data[,unified_subcluster_ident] == "sub_4") & (obj@meta.data[,"sub_4_iCAF_score"] < obj@meta.data[,"sub_4_myCAF_score"])] <- "myCAF"
alphas <- rep(1, times = length(obj$sub_4_CAF_call))
alphas[is.na(obj$sub_4_CAF_call)] <- 0
pdf(paste0(out_dir,"/","Dimplot_",unified_subcluster_ident,"_",prefix,"_res0.3_sub_4_CAF_call_20241212.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = "sub_4_CAF_call", alpha = alphas, label=TRUE, label.size=2, raster=FALSE), layers='Point', dpi=300))
scores_to_plot <- c("sub_4_iCAF_score","sub_4_myCAF_score")
for (score in scores_to_plot) {
    print(rasterize(FeaturePlot(obj, reduction = "umap.50PC", features = score, label=T, label.size=1, raster=FALSE, order = T) + scale_colour_gradient2(low = "blue", mid = "lightgrey", high = "red"), layers='Point', dpi=300))
}
dev.off()
obj$sub_iCAF_score <- 0
obj$sub_iCAF_score[obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting] <- obj$iCAF_score1[obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting]
obj$sub_myCAF_score <- 0
obj$sub_myCAF_score[obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting] <- obj$myCAF_score1[obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting]
scores_to_plot <- c("sub_iCAF_score","sub_myCAF_score")
pdf(paste0(out_dir,"/","Dimplot_",unified_subcluster_ident,"_",prefix,"_res0.3_sub_CAF_scores_20241212.pdf"),useDingbats = F, width=8, height=6)
for (score in scores_to_plot) {
    print(rasterize(FeaturePlot(obj, reduction = "umap.50PC", features = score, label=T, label.size=1, raster=FALSE,order = T) + scale_colour_gradient2(low = "blue", mid = "lightgrey", high = "red"), layers='Point', dpi=300))
}
dev.off()
obj$sub_CAF_call <- NA
obj$sub_CAF_call[(obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting) & (obj@meta.data[,"sub_iCAF_score"] >= obj@meta.data[,"sub_myCAF_score"])] <- "iCAF"
obj$sub_CAF_call[(obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting) & (obj@meta.data[,"sub_iCAF_score"] < obj@meta.data[,"sub_myCAF_score"])] <- "myCAF"
alphas <- rep(1, times = length(obj$sub_CAF_call))
alphas[is.na(obj$sub_CAF_call)] <- 0
pdf(paste0(out_dir,"/","Dimplot_",unified_subcluster_ident,"_",prefix,"_res0.3_sub_CAF_call_20241212.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = "sub_CAF_call", alpha = alphas, label=TRUE, label.size=2, raster=FALSE), layers='Point', dpi=300))
scores_to_plot <- c("sub_iCAF_score","sub_myCAF_score")
for (score in scores_to_plot) {
    print(rasterize(FeaturePlot(obj, reduction = "umap.50PC", features = score, label=T, label.size=1, raster=FALSE, order = T) + scale_colour_gradient2(low = "blue", mid = "lightgrey", high = "red"), layers='Point', dpi=300))
}
dev.off()
table(obj@meta.data[,unified_subcluster_ident][obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting], obj@meta.data[,"sub_CAF_call"][obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting])
#         iCAF myCAF
# sub_0  72332 18718
# sub_1  11915 41577
# sub_2   3971 45972
# sub_3   4389 20238
# sub_4   8678  9098
# sub_5  11210  3305
# sub_6   7523  1047
# sub_7   1719  6828
# sub_8    150  6988
# sub_9    416  1102 
# sub_10   806    66 
# sub_0 iCAF
# sub_1 myCAF
# sub_2 myCAF
# sub_3 myCAF
# sub_4 mix CAF
# sub_5 iCAF_myeloid_mix??? this also has macrophages from the input Msub_8 mixed in here. so this cluster needs to be split in.
# sub_6 iCAF
# sub_7 myCAF
# sub_8 myCAF 
# sub_9 myCAF
# sub_10 iCAF
# 18718 + 11915 + 3971 + 4389 + 3305 + 1047 + 1719 + 150 + 416 + 66 = 45696
# 72332 +18718 + 11915+ 41577 +  3971+ 45972+   4389+ 20238+  11210+  3305+   7523+  1047 +  1719 + 6828+    150 + 6988  +  416 + 1102  +  806  +  66 = 260272
# This means if I rely on the GeneModuleScores instead of the clusters and dotplots then (45696/260272)*100 = 17% of the cells 
# will be assigned to different categories. Since the clusters should represent cell that are more similar I'm going to go with them.
# Then I'm going to sub cluster sub_4 at low resolution since it is appears mixed and assign as either myCAF or iCAF.
# for sub_5 which looks like a mix of myeloid cells and CAFs I'm going to subcluster at low resolution to try to pull out the Myeloid population
obj$CAF_sub_labels <- obj@meta.data[,unified_subcluster_ident]
obj@meta.data[,"CAF_sub_labels"][(obj@meta.data[,unified_subcluster_ident] %in% c("sub_0","sub_6","sub_10")) ] <- "iCAF"
obj@meta.data[,"CAF_sub_labels"][(obj@meta.data[,unified_subcluster_ident] %in% c("sub_1","sub_2","sub_3","sub_7","sub_8","sub_9")) ] <- "myCAF"
obj@meta.data[,"CAF_sub_labels"][(obj@meta.data[,unified_subcluster_ident] %in% c("sub_4")) ] <- "Fsub_4"
obj@meta.data[,"CAF_sub_labels"][(obj@meta.data[,unified_subcluster_ident] %in% c("sub_5")) ] <- "Fsub_5"
Idents(obj) <- "CAF_sub_labels"
plan(sequential)
subcluster_ident <- "merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8.Fsub_4"
subcluster_idents <- c('Fsub_4')
obj <- FindSubCluster(obj, cluster = subcluster_idents, graph.name = "SCT_snn", subcluster.name = subcluster_ident, resolution = 0.1, algorithm = 1) #based on the recommendation here https://github.com/satijalab/seurat/issues/6682 by the authors I used the SCT_snn instead of SCT_nn graph
cells_to_highlight <- list()
merge_idents_highlight <- c()
for (new_ident in subcluster_idents) {
    merge_idents_highlight_single <- unique(obj@meta.data[,subcluster_ident])[grepl(paste0(new_ident,"_"),unique(obj@meta.data[,subcluster_ident]))]
    merge_idents_highlight <- unique(c(merge_idents_highlight, merge_idents_highlight_single))
}
unified_subcluster_ident <- paste0("unified.",subcluster_ident)
obj@meta.data[,unified_subcluster_ident] <- obj@meta.data[,subcluster_ident]
unified_cells_to_highlight <- list()
merge_unified_subcluster_label <- c()
for (sub_ident in merge_idents_highlight) {
    print(sub_ident)
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,subcluster_ident] == sub_ident),"barcode"]
    cells_to_highlight[[sub_ident]] <- barcodes_in_sub
    new_subcluster <- unlist(strsplit(sub_ident, "_"))[length(unlist(strsplit(sub_ident, "_")))]
    new_subcluster_label <- paste0("sub_",new_subcluster)
    obj@meta.data[,unified_subcluster_ident][obj@meta.data[,subcluster_ident] == sub_ident] <- new_subcluster_label
    merge_unified_subcluster_label <- unique(c(merge_unified_subcluster_label, new_subcluster_label))
}
print(unique(obj@meta.data[,unified_subcluster_ident]))
unified_cells_to_highlight <- list()
for (unified_subcluster_label in merge_unified_subcluster_label) {
    print(unified_subcluster_label)
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,unified_subcluster_ident] == unified_subcluster_label),"barcode"]
    unified_cells_to_highlight[[unified_subcluster_label]] <- barcodes_in_sub
    print(summary(obj@meta.data[,"nCount_Xenium"][obj@meta.data[,unified_subcluster_ident] == unified_subcluster_label]))
}
# [1] "sub_0"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00   31.00   60.50   80.82  109.00  672.00 
# [1] "sub_1"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00   13.00   26.00   41.17   52.00  426.00 
# [1] "sub_2"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00   12.00   20.00   29.88   37.00  279.00 
# [1] "sub_3"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.00   32.00   66.50   83.24  118.00  329.00
Idents(obj) <- unified_subcluster_ident
tmp_df <- obj@meta.data[,c('barcode',subcluster_ident)]
colnames(tmp_df) <- c('barcode',subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
tmp_df <- obj@meta.data[,c('barcode',unified_subcluster_ident)]
colnames(tmp_df) <- c('barcode',unified_subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",unified_subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
# Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
# /diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
# ./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
# /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8.Fsub_4.tsv \
# unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8.Fsub_4
library(scCustomize)
varibow_pal <- DiscretePalette_scCustomize(num_colors = length(unified_cells_to_highlight), palette = "varibow")
ordered_labels_for_plotting <- c()
for (i in 0:(length(merge_unified_subcluster_label)-1)) {
    unified_subcluster_label <- paste0("sub_",i)
    ordered_labels_for_plotting <- c(ordered_labels_for_plotting, unified_subcluster_label)
}
unified_cells_to_highlight <- unified_cells_to_highlight[ordered_labels_for_plotting]
pdf(paste0(out_dir,"/","Dimplot_",unified_subcluster_ident,"_",prefix,"_res0.3_highlight_subclusters_20241212.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, cells.highlight = unified_cells_to_highlight, cols.highlight = varibow_pal, sizes.highlight=0.1, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
for (i in 1:length(ordered_labels_for_plotting)) {
    unified_subcluster_label <- ordered_labels_for_plotting[i]
    print(unified_subcluster_label)
    unified_cells_to_highlight_sep <- unified_cells_to_highlight[[unified_subcluster_label]]
    col_highlight <- varibow_pal[i]
    print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, cells.highlight = unified_cells_to_highlight_sep, cols.highlight = col_highlight, sizes.highlight=0.1, label=TRUE, label.size=4, raster=FALSE) + labs(title = paste0("highlight ",unified_subcluster_label)), layers='Point', dpi=300))
}
dev.off()
print(table(obj$cell_type_individual[obj@meta.data[,c(unified_subcluster_ident)] %in% merge_unified_subcluster_label],
            obj@meta.data[,unified_subcluster_ident][obj@meta.data[,c(unified_subcluster_ident)] %in% merge_unified_subcluster_label]))
print(table(obj@meta.data[,unified_subcluster_ident]))
#    10      11      12    13_1   13_10   13_11   13_12    13_2    13_4    13_5 
# 31171   27929   21689    3389     140     123     118    2693    2556    2245 
#  13_6    13_7    13_8    13_9      15      16      17      18    19_1    19_3 
#   547     393     313     286   14709   12988   10848    8639    1429    1258 
#  19_5       2      20      22      23      24      26      27      28      29 
#   580  111752    8055    7071    5079    5029    2052    2015    1982    1722 
#     3      30      31      32      33      34       4       5       6       7 
# 70375    1572    1529     495      65      19   63052   52907   47539   44723 
#     8       9  Fsub_5    iCAF   myCAF Myeloid   sub_0   sub_1   sub_2   sub_3 
# 36987   31289   14515  100492  145265  135497   10992    4090    2226     468
# sub_0 myCAF
# sub_1 iCAF
# sub_2 iCAF
# sub_3 iCAF
Idents(obj) <- "CAF_sub_labels"
plan(sequential)
subcluster_ident <- "merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8.Fsub_5"
subcluster_idents <- c('Fsub_5')
obj <- FindSubCluster(obj, cluster = subcluster_idents, graph.name = "SCT_snn", subcluster.name = subcluster_ident, resolution = 0.07, algorithm = 1) #based on the recommendation here https://github.com/satijalab/seurat/issues/6682 by the authors I used the SCT_snn instead of SCT_nn graph
cells_to_highlight <- list()
merge_idents_highlight <- c()
for (new_ident in subcluster_idents) {
    merge_idents_highlight_single <- unique(obj@meta.data[,subcluster_ident])[grepl(paste0(new_ident,"_"),unique(obj@meta.data[,subcluster_ident]))]
    merge_idents_highlight <- unique(c(merge_idents_highlight, merge_idents_highlight_single))
}
unified_subcluster_ident <- paste0("unified.",subcluster_ident)
obj@meta.data[,unified_subcluster_ident] <- obj@meta.data[,subcluster_ident]
unified_cells_to_highlight <- list()
merge_unified_subcluster_label <- c()
for (sub_ident in merge_idents_highlight) {
    print(sub_ident)
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,subcluster_ident] == sub_ident),"barcode"]
    cells_to_highlight[[sub_ident]] <- barcodes_in_sub
    new_subcluster <- unlist(strsplit(sub_ident, "_"))[length(unlist(strsplit(sub_ident, "_")))]
    new_subcluster_label <- paste0("sub_",new_subcluster)
    obj@meta.data[,unified_subcluster_ident][obj@meta.data[,subcluster_ident] == sub_ident] <- new_subcluster_label
    merge_unified_subcluster_label <- unique(c(merge_unified_subcluster_label, new_subcluster_label))
}
print(unique(obj@meta.data[,unified_subcluster_ident]))
unified_cells_to_highlight <- list()
for (unified_subcluster_label in merge_unified_subcluster_label) {
    print(unified_subcluster_label)
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,unified_subcluster_ident] == unified_subcluster_label),"barcode"]
    unified_cells_to_highlight[[unified_subcluster_label]] <- barcodes_in_sub
    print(summary(obj@meta.data[,"nCount_Xenium"][obj@meta.data[,unified_subcluster_ident] == unified_subcluster_label]))
}
# [1] "sub_0"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00   37.00   68.00   82.42  112.00  886.00 
# [1] "sub_1"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0     7.0    14.0    23.4    29.0   277.0 
# [1] "sub_2"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.0    11.0    21.0    30.9    38.0   388.0 
# [1] "sub_3"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00   10.00   16.00   25.56   31.00  203.00
Idents(obj) <- unified_subcluster_ident
tmp_df <- obj@meta.data[,c('barcode',subcluster_ident)]
colnames(tmp_df) <- c('barcode',subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
tmp_df <- obj@meta.data[,c('barcode',unified_subcluster_ident)]
colnames(tmp_df) <- c('barcode',unified_subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",unified_subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
print("plots")
# Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
# /diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
# ./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
# /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8.Fsub_4.tsv \
# unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8.Fsub_5
library(scCustomize)
varibow_pal <- DiscretePalette_scCustomize(num_colors = length(unified_cells_to_highlight), palette = "varibow")
ordered_labels_for_plotting <- c()
for (i in 0:(length(merge_unified_subcluster_label)-1)) {
    unified_subcluster_label <- paste0("sub_",i)
    ordered_labels_for_plotting <- c(ordered_labels_for_plotting, unified_subcluster_label)
}
unified_cells_to_highlight <- unified_cells_to_highlight[ordered_labels_for_plotting]
pdf(paste0(out_dir,"/","Dimplot_",unified_subcluster_ident,"_",prefix,"_res0.3_highlight_subclusters_20241212.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, cells.highlight = unified_cells_to_highlight, cols.highlight = varibow_pal, sizes.highlight=0.1, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
for (i in 1:length(ordered_labels_for_plotting)) {
    unified_subcluster_label <- ordered_labels_for_plotting[i]
    print(unified_subcluster_label)
    unified_cells_to_highlight_sep <- unified_cells_to_highlight[[unified_subcluster_label]]
    col_highlight <- varibow_pal[i]
    print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, cells.highlight = unified_cells_to_highlight_sep, cols.highlight = col_highlight, sizes.highlight=0.1, label=TRUE, label.size=4, raster=FALSE) + labs(title = paste0("highlight ",unified_subcluster_label)), layers='Point', dpi=300))
}
dev.off()
# I have tried 0.1, 0.07, and 0.05 in the resolution parameter for FindSubCluster
# after reviewing all of the above I have come to the conclusion that by including 
# Msub_8 in the subclustering of Fibroblasts in an attempt to subset out the small 
# amount of fibroblasts it contained I instead ended up being unable separate out the myeloid cells from the fibroblasts during the fibroblast subclustering
# so I'm just going to redo the subclustering for fibroblasts and leave Msub_8 out 
Idents(obj) <- "merge.sub.1.14.25.19_0.19_4.19_6_merge.sub.13"
plan(sequential)
subcluster_ident <- "merge.sub.0.21.19_2.19_7.13_3.13_0"
subcluster_idents <- c('0','21','19_2','19_7','13_3','13_0')
obj <- FindSubCluster(obj, cluster = subcluster_idents, graph.name = "SCT_snn", subcluster.name = subcluster_ident, resolution = 0.3, algorithm = 1) #based on the recommendation here https://github.com/satijalab/seurat/issues/6682 by the authors I used the SCT_snn instead of SCT_nn graph
cells_to_highlight <- list()
merge_idents_highlight <- c()
for (new_ident in subcluster_idents) {
    merge_idents_highlight_single <- unique(obj@meta.data[,subcluster_ident])[grepl(paste0(new_ident,"_"),unique(obj@meta.data[,subcluster_ident]))]
    merge_idents_highlight <- unique(c(merge_idents_highlight, merge_idents_highlight_single))
}
unified_subcluster_ident <- paste0("unified.",subcluster_ident)
obj@meta.data[,unified_subcluster_ident] <- obj@meta.data[,subcluster_ident]
unified_cells_to_highlight <- list()
merge_unified_subcluster_label <- c()
for (sub_ident in merge_idents_highlight) {
    print(sub_ident)
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,subcluster_ident] == sub_ident),"barcode"]
    cells_to_highlight[[sub_ident]] <- barcodes_in_sub
    new_subcluster <- unlist(strsplit(sub_ident, "_"))[length(unlist(strsplit(sub_ident, "_")))]
    new_subcluster_label <- paste0("sub_",new_subcluster)
    obj@meta.data[,unified_subcluster_ident][obj@meta.data[,subcluster_ident] == sub_ident] <- new_subcluster_label
    merge_unified_subcluster_label <- unique(c(merge_unified_subcluster_label, new_subcluster_label))
}
print(unique(obj@meta.data[,unified_subcluster_ident]))
unified_cells_to_highlight <- list()
for (unified_subcluster_label in merge_unified_subcluster_label) {
    print(unified_subcluster_label)
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,unified_subcluster_ident] == unified_subcluster_label),"barcode"]
    unified_cells_to_highlight[[unified_subcluster_label]] <- barcodes_in_sub
    print(summary(obj@meta.data[,"nCount_Xenium"][obj@meta.data[,unified_subcluster_ident] == unified_subcluster_label]))
}
# [1] "sub_0"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.0    60.0    94.0   106.6   139.0   721.0                                                                                                                     
# [1] "sub_2"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.00   35.00   73.00   94.09  128.00  842.00                                                                                                                     
# [1] "sub_6"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   47.00   87.00   95.02  131.00  514.00 
# [1] "sub_4"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   20.00   45.00   65.43   89.00  730.00 
# [1] "sub_8"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00    9.00   16.00   22.79   29.00  317.00 
# [1] "sub_7"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   15.00   32.00   52.49   71.00  424.00 
# [1] "sub_1"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.0    66.0   103.0   118.9   154.0   791.0 
# [1] "sub_5"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00   37.00   70.00   86.34  118.00  778.00 
# [1] "sub_3"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00   28.00   60.00   78.74  110.00  883.00 
# [1] "sub_9"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00   26.00   61.00   87.53  125.75  609.00
Idents(obj) <- unified_subcluster_ident
tmp_df <- obj@meta.data[,c('barcode',subcluster_ident)]
colnames(tmp_df) <- c('barcode',subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
tmp_df <- obj@meta.data[,c('barcode',unified_subcluster_ident)]
colnames(tmp_df) <- c('barcode',unified_subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",unified_subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
# Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
# /diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
# ./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
# /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_unified.merge.sub.0.21.19_2.19_7.13_3.13_0.tsv \
# unified.merge.sub.0.21.19_2.19_7.13_3.13_0
# results are in /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/marker_dotplots/unified.merge.sub.0.21.19_2.19_7.13_3.13_0.PDAC_merge_primary_KRAS_20241205_single_SCT.Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.pdf
print("plots")
library(scCustomize)
varibow_pal <- DiscretePalette_scCustomize(num_colors = length(unified_cells_to_highlight), palette = "varibow")
ordered_labels_for_plotting <- c()
for (i in 0:(length(merge_unified_subcluster_label)-1)) {
    unified_subcluster_label <- paste0("sub_",i)
    ordered_labels_for_plotting <- c(ordered_labels_for_plotting, unified_subcluster_label)
}
unified_cells_to_highlight <- unified_cells_to_highlight[ordered_labels_for_plotting]
pdf(paste0(out_dir,"/","Dimplot_",unified_subcluster_ident,"_",prefix,"_res0.3_highlight_subclusters_20241212.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, cells.highlight = unified_cells_to_highlight, cols.highlight = varibow_pal, sizes.highlight=0.1, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
for (i in 1:length(ordered_labels_for_plotting)) {
    unified_subcluster_label <- ordered_labels_for_plotting[i]
    print(unified_subcluster_label)
    unified_cells_to_highlight_sep <- unified_cells_to_highlight[[unified_subcluster_label]]
    col_highlight <- varibow_pal[i]
    print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, cells.highlight = unified_cells_to_highlight_sep, cols.highlight = col_highlight, sizes.highlight=0.1, label=TRUE, label.size=4, raster=FALSE) + labs(title = paste0("highlight ",unified_subcluster_label)), layers='Point', dpi=300))
}
dev.off()
print(table(obj$cell_type_individual[obj$merge.sub.19 %in% subcluster_idents], obj$sample_ID[obj$merge.sub.19 %in% subcluster_idents]))
print(table(obj$cell_type_individual[obj@meta.data[,c(unified_subcluster_ident)] %in% merge_unified_subcluster_label],
            obj@meta.data[,unified_subcluster_ident][obj@meta.data[,c(unified_subcluster_ident)] %in% merge_unified_subcluster_label]))
print(table(obj@meta.data[,unified_subcluster_ident]))
for (sample_ID in sample_vector) {
    tmp_df <- obj@meta.data[(obj@meta.data[,"sample_ID"] == sample_ID),c('original_Xenium_barcode',unified_subcluster_ident)]
    colnames(tmp_df) <- c('cell_id', 'group')
    write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",sample_ID,"_",prefix,"_",unified_subcluster_ident,".csv"),sep=",",quote=F,row.names=F)
}
table(obj$sample_ID[obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting], obj@meta.data[,unified_subcluster_ident][obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting])
#                     sub_0 sub_1 sub_2 sub_3 sub_4 sub_5 sub_6 sub_7 sub_8
# HT060P1-S1R1Fp1U1   51544 34461 15082  9630  9844  8462  5109  3181  1835
# HT061P1-S1P1A1L1U1   1900   938  4011  1115   413   116   183    53    78
# HT061P1-S1P1A1L4U1   1351  5191  3592  5054   749   346   167   133   144
# HT125P1-S1H4A1L1U1   9622  7365  6659  6297  1449    97   811   306   505
# HT125P1-S1H8A1U1     1099   934  1675  3810   604    53   382   184   140
# HT185P1-S1H2L1U1     4687  1406  7918  3087  2044   440   342   343  1020
# HT227P1-S1H1L1U1     1183   493  2486  3493  4682   241   341   528   534
# HT242P1-S1H4L4U1     1089   657  1422  4966   884  1014   402  2343   229
# HT270P1-S1H1A1US2_1  6946   802  3220  2213  1215   283   274   109   184
# 
#                     sub_9
# HT060P1-S1R1Fp1U1     504
# HT061P1-S1P1A1L1U1     35
# HT061P1-S1P1A1L4U1    168
# HT125P1-S1H4A1L1U1    101
# HT125P1-S1H8A1U1      120
# HT185P1-S1H2L1U1      103
# HT227P1-S1H1L1U1      195
# HT242P1-S1H4L4U1      237
# HT270P1-S1H1A1US2_1    39
# dotplot calls of subclusters
# sub_0 iCAF
# sub_1 myCAF
# sub_2 myCAF
# sub_3 myCAF
# sub_4 iCAF with some myCAF
# sub_5 myCAF
# sub_6 iCAF
# sub_7 mycAF
# sub_8 iCAF
# sub_9 CAF_proliferative
# sub_4 is a mix of different CAFs so I'm going to sublcuster just that population
obj$CAF_v2_sub_labels <- obj@meta.data[,unified_subcluster_ident]
obj@meta.data[,"CAF_v2_sub_labels"][(obj@meta.data[,unified_subcluster_ident] %in% c("sub_0","sub_6","sub_8")) ] <- "iCAF"
obj@meta.data[,"CAF_v2_sub_labels"][(obj@meta.data[,unified_subcluster_ident] %in% c("sub_1","sub_2","sub_3","sub_5","sub_7")) ] <- "myCAF"
obj@meta.data[,"CAF_v2_sub_labels"][(obj@meta.data[,unified_subcluster_ident] %in% c("sub_9")) ] <- "CAF_proliferative"
obj@meta.data[,"CAF_v2_sub_labels"][(obj@meta.data[,unified_subcluster_ident] %in% c("sub_4")) ] <- "Fsub_4"
Idents(obj) <- "CAF_v2_sub_labels"
plan(sequential)
subcluster_ident <- "merge.sub.0.21.19_2.19_7.13_3.13_0.Fsub_4"
subcluster_idents <- c('Fsub_4')
obj <- FindSubCluster(obj, cluster = subcluster_idents, graph.name = "SCT_snn", subcluster.name = subcluster_ident, resolution = 0.1, algorithm = 1) #based on the recommendation here https://github.com/satijalab/seurat/issues/6682 by the authors I used the SCT_snn instead of SCT_nn graph
cells_to_highlight <- list()
merge_idents_highlight <- c()
for (new_ident in subcluster_idents) {
    merge_idents_highlight_single <- unique(obj@meta.data[,subcluster_ident])[grepl(paste0(new_ident,"_"),unique(obj@meta.data[,subcluster_ident]))]
    merge_idents_highlight <- unique(c(merge_idents_highlight, merge_idents_highlight_single))
}
unified_subcluster_ident <- paste0("unified.",subcluster_ident)
obj@meta.data[,unified_subcluster_ident] <- obj@meta.data[,subcluster_ident]
unified_cells_to_highlight <- list()
merge_unified_subcluster_label <- c()
for (sub_ident in merge_idents_highlight) {
    print(sub_ident)
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,subcluster_ident] == sub_ident),"barcode"]
    cells_to_highlight[[sub_ident]] <- barcodes_in_sub
    new_subcluster <- unlist(strsplit(sub_ident, "_"))[length(unlist(strsplit(sub_ident, "_")))]
    new_subcluster_label <- paste0("sub_",new_subcluster)
    obj@meta.data[,unified_subcluster_ident][obj@meta.data[,subcluster_ident] == sub_ident] <- new_subcluster_label
    merge_unified_subcluster_label <- unique(c(merge_unified_subcluster_label, new_subcluster_label))
}
print(unique(obj@meta.data[,unified_subcluster_ident]))
unified_cells_to_highlight <- list()
for (unified_subcluster_label in merge_unified_subcluster_label) {
    print(unified_subcluster_label)
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,unified_subcluster_ident] == unified_subcluster_label),"barcode"]
    unified_cells_to_highlight[[unified_subcluster_label]] <- barcodes_in_sub
    print(summary(obj@meta.data[,"nCount_Xenium"][obj@meta.data[,unified_subcluster_ident] == unified_subcluster_label]))
}
# [1] "sub_1"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    9.00   17.00   25.86   33.00  388.00 
# [1] "sub_0"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.0    30.0    58.0    77.8   105.0   730.0 
Idents(obj) <- unified_subcluster_ident
tmp_df <- obj@meta.data[,c('barcode',subcluster_ident)]
colnames(tmp_df) <- c('barcode',subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
tmp_df <- obj@meta.data[,c('barcode',unified_subcluster_ident)]
colnames(tmp_df) <- c('barcode',unified_subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",unified_subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
# Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
# /diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
# ./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
# /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Fsub_4.tsv \
# unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Fsub_4
print("plots")
library(scCustomize)
varibow_pal <- DiscretePalette_scCustomize(num_colors = length(unified_cells_to_highlight), palette = "varibow")
ordered_labels_for_plotting <- c()
for (i in 0:(length(merge_unified_subcluster_label)-1)) {
    unified_subcluster_label <- paste0("sub_",i)
    ordered_labels_for_plotting <- c(ordered_labels_for_plotting, unified_subcluster_label)
}
unified_cells_to_highlight <- unified_cells_to_highlight[ordered_labels_for_plotting]
pdf(paste0(out_dir,"/","Dimplot_",unified_subcluster_ident,"_",prefix,"_res0.3_highlight_subclusters_20241212.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, cells.highlight = unified_cells_to_highlight, cols.highlight = varibow_pal, sizes.highlight=0.1, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
for (i in 1:length(ordered_labels_for_plotting)) {
    unified_subcluster_label <- ordered_labels_for_plotting[i]
    print(unified_subcluster_label)
    unified_cells_to_highlight_sep <- unified_cells_to_highlight[[unified_subcluster_label]]
    col_highlight <- varibow_pal[i]
    print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, cells.highlight = unified_cells_to_highlight_sep, cols.highlight = col_highlight, sizes.highlight=0.1, label=TRUE, label.size=4, raster=FALSE) + labs(title = paste0("highlight ",unified_subcluster_label)), layers='Point', dpi=300))
}
dev.off()
print(table(obj$cell_type_individual[obj$merge.sub.19 %in% subcluster_idents], obj$sample_ID[obj$merge.sub.19 %in% subcluster_idents]))
print(table(obj$cell_type_individual[obj@meta.data[,c(unified_subcluster_ident)] %in% merge_unified_subcluster_label],
            obj@meta.data[,unified_subcluster_ident][obj@meta.data[,c(unified_subcluster_ident)] %in% merge_unified_subcluster_label]))
print(table(obj@meta.data[,unified_subcluster_ident]))
#     10                11                12              13_1                                                                                             
#  31171             27929             21689              3389                                                                                             
#  13_10             13_11             13_12              13_2                                                                                             
#    140               123               118              2693                                                                                             
#   13_4              13_5              13_6              13_7                                                                                             
#   2556              2245               547               393                                                                                             
#   13_8              13_9                15                16                                                                                             
#    313               286             14709             12988                                                                                             
#     17                18              19_1              19_3                                                                                             
#  10848              8639              1429              1258                                                                                             
#   19_5                 2                20                22                                                                                             
#    580            111752              8055              7071                                                                                             
#     23                24                26                27                                                                                             
#   5079              5029              2052              2015                                                                                             
#     28                29                 3                30 
#   1982              1722             70375              1572 
#     31                32                33                34 
#   1529               495                65                19 
#      4                 5                 6                 7 
#  63052             52907             47539             44723 
#      8                 9 CAF_proliferative              iCAF 
#  36987             31289              1502             92101 
# Msub_8             myCAF           Myeloid             sub_0 
#   6352            156209            135497             16672 
#  sub_1 
#   5212
for (sample_ID in sample_vector) {
    tmp_df <- obj@meta.data[(obj@meta.data[,"sample_ID"] == sample_ID),c('original_Xenium_barcode',unified_subcluster_ident)]
    colnames(tmp_df) <- c('cell_id', 'group')
    write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",sample_ID,"_",prefix,"_",unified_subcluster_ident,".csv"),sep=",",quote=F,row.names=F)
}
table(obj$sample_ID[obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting], obj@meta.data[,unified_subcluster_ident][obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting])
#                     sub_0 sub_1
# HT060P1-S1R1Fp1U1    7623  2221
# HT061P1-S1P1A1L1U1    322    91
# HT061P1-S1P1A1L4U1    635   114
# HT125P1-S1H4A1L1U1   1092   357
# HT125P1-S1H8A1U1      435   169
# HT185P1-S1H2L1U1     1329   715
# HT227P1-S1H1L1U1     3664  1018
# HT242P1-S1H4L4U1      682   202
# HT270P1-S1H1A1US2_1   890   325
# dotplot calls of subclusters
# sub_0 iCAF
# sub_1 iCAF
# this doesn't separate out any myCAF-like cells so I'm just going to call them iCAFs.

# I am going to subcluster the ductal and acinar cells in order to make sure the acinar, 
# reactive ductal, duct_like_2, and duct_like_1 cells are accurately and consistently assigned. 
# I can see in the file Table_cell_type_ind_by_sc_50PCPDAC_merge_primary_KRAS_20241205_single_SCT.tsv
# that the majority of the cells that were called as PanIN when doing the individual cell typing (cell_type_individual)
# have clustered separately in seurat cluster 16. However, there are some cells that were called as PanIN that are now part of other clusters
# and some other cells that are now listed as PanIN. In order to see if these are actually PanIN (which should be confirmed histologically)
# I am going to combine the seurat cluster by cell_type_individual labels to one column and then import that 
# into Xenium Explorer to manually inspect.
obj$unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Fsub_4.by.cell_type_individual <- paste0(obj$unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Fsub_4,"__",obj$cell_type_individual)
for (sample_ID in sample_vector) {
    tmp_df <- obj@meta.data[(obj@meta.data[,"sample_ID"] == sample_ID),c('original_Xenium_barcode',"unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Fsub_4.by.cell_type_individual")]
    colnames(tmp_df) <- c('cell_id', 'group')
    write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",sample_ID,"_",prefix,"_","unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Fsub_4.by.cell_type_individual",".csv"),sep=",",quote=F,row.names=F)
}
# upon inspection of the clusters in the histology:
# HT060P1-S1R1Fp1U1: in cluster 16 are unlikely to be PanIN and seem to be more a part of normal reactive ductal structures than anything.
# HT061P1-S1P1A1L1U1: the original tumor vs PanIN labels are probably more accurate.
# HT061P1-S1P1A1L4U1: the original tumor vs panin labels are probably more accurate
# HT125P1-S1H4A1L1U1: the original tumor vs panin labels are probably more accurate
# HT125P1-S1H8A1U1: the original tumor vs panin labels are probably more accurate
# HT270P1-S1H1A1US2_1: the original tumor vs panin labels are probably just as accurate as the merged seurat cluster labels.
# based on this I am going to just re-cluster the normal duct and acinar clusters as a result.
# I am going to leave the cluster 16 labels as is besides,
# just updating them so the names are consistent with the rest of the labels.
# For Seurat Cluster 4 and 7 I am going to also leave the labels as is since when I transfer them to the Xenium Explorer application
# I can see that some cells in this cluster are cells adjacent to tumor cells that are unlikely to actually be tumor cells.
print(unified_subcluster_ident)
# [1] "unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Fsub_4"
# Don't need to include any of clusters 19 or 13 in this so I'm just going to use the base Seurat Clusters
Idents(obj) <- "seurat_clusters"
plan(sequential)
subcluster_ident <- "merge.sub.3.11.res_0.3"
subcluster_idents <- c('3','11')
obj <- FindSubCluster(obj, cluster = subcluster_idents, graph.name = "SCT_snn", subcluster.name = subcluster_ident, resolution = 0.3, algorithm = 1) #based on the recommendation here https://github.com/satijalab/seurat/issues/6682 by the authors I used the SCT_snn instead of SCT_nn graph
cells_to_highlight <- list()
merge_idents_highlight <- c()
for (new_ident in subcluster_idents) {
    merge_idents_highlight_single <- unique(obj@meta.data[,subcluster_ident])[grepl(paste0(new_ident,"_"),unique(obj@meta.data[,subcluster_ident]))]
    merge_idents_highlight <- unique(c(merge_idents_highlight, merge_idents_highlight_single))
}
unified_subcluster_ident <- paste0("unified.",subcluster_ident)
obj@meta.data[,unified_subcluster_ident] <- obj@meta.data[,subcluster_ident]
unified_cells_to_highlight <- list()
merge_unified_subcluster_label <- c()
for (sub_ident in merge_idents_highlight) {
    print(sub_ident)
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,subcluster_ident] == sub_ident),"barcode"]
    cells_to_highlight[[sub_ident]] <- barcodes_in_sub
    new_subcluster <- unlist(strsplit(sub_ident, "_"))[length(unlist(strsplit(sub_ident, "_")))]
    new_subcluster_label <- paste0("sub_",new_subcluster)
    obj@meta.data[,unified_subcluster_ident][obj@meta.data[,subcluster_ident] == sub_ident] <- new_subcluster_label
    merge_unified_subcluster_label <- unique(c(merge_unified_subcluster_label, new_subcluster_label))
}
print(unique(obj@meta.data[,unified_subcluster_ident]))
#  [1] "10"     "0"      "2"      "9"      "5"      "12"     "1"      "8"                                                                                             
#  [9] "sub_3"  "18"     "sub_2"  "sub_10" "sub_0"  "27"     "20"     "sub_1"                                                                                         
# [17] "6"      "13"     "21"     "24"     "19"     "31"     "sub_6"  "sub_12"                                                                                        
# [25] "25"     "sub_5"  "sub_4"  "4"      "17"     "14"     "15"     "30"                                                                                            
# [33] "28"     "16"     "7"      "22"     "29"     "26"     "23"     "32"                                                                                            
# [41] "sub_7"  "sub_11" "sub_9"  "33"     "34"     "sub_8"
unified_cells_to_highlight <- list()
for (unified_subcluster_label in merge_unified_subcluster_label) {
    print(unified_subcluster_label)
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,unified_subcluster_ident] == unified_subcluster_label),"barcode"]
    unified_cells_to_highlight[[unified_subcluster_label]] <- barcodes_in_sub
    print(summary(obj@meta.data[,"nCount_Xenium"][obj@meta.data[,unified_subcluster_ident] == unified_subcluster_label]))
}
# [1] "sub_3"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 4.0    59.0   108.0   116.1   160.0   614.0                                                                                                                     
# [1] "sub_0"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 5.00   45.00   79.00   87.94  122.00  481.00                                                                                                                     
# [1] "sub_2"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 1.00   16.00   35.00   54.27   79.00  488.00                                                                                                                     
# [1] "sub_1"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 1.0    82.0   128.0   141.5   186.0   727.0                                                                                                                     
# [1] "sub_12"                                                                                                                                                        
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 29.0   113.8   155.0   172.8   216.5   487.0
# [1] "sub_10"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.0    78.0   113.0   126.4   162.0   467.0 
# [1] "sub_6"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   10.00   19.00   44.12   53.00  528.00 
# [1] "sub_5"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.00   46.00   74.00   87.83  119.00  387.00 
# [1] "sub_7"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.0    90.0   126.0   128.5   163.8   455.0 
# [1] "sub_4"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.00   59.00   81.00   84.42  106.00  335.00 
# [1] "sub_9"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 10.0    56.0    88.0   105.9   136.0   476.0 
# [1] "sub_11"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 12.0    78.0   114.0   120.1   156.0   302.0 
# [1] "sub_8"
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.00   32.00   59.00   76.05  102.00  442.00
Idents(obj) <- unified_subcluster_ident
tmp_df <- obj@meta.data[,c('barcode',subcluster_ident)]
colnames(tmp_df) <- c('barcode',subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
tmp_df <- obj@meta.data[,c('barcode',unified_subcluster_ident)]
colnames(tmp_df) <- c('barcode',unified_subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",unified_subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
# Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
# /diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
# ./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
# /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_unified.merge.sub.3.11.res_0.3.tsv \
# unified.merge.sub.3.11.res_0.3
print("plots")
library(scCustomize)
varibow_pal <- DiscretePalette_scCustomize(num_colors = length(unified_cells_to_highlight), palette = "varibow")
ordered_labels_for_plotting <- c()
for (i in 0:(length(merge_unified_subcluster_label)-1)) {
    unified_subcluster_label <- paste0("sub_",i)
    ordered_labels_for_plotting <- c(ordered_labels_for_plotting, unified_subcluster_label)
}
unified_cells_to_highlight <- unified_cells_to_highlight[ordered_labels_for_plotting]
pdf(paste0(out_dir,"/","Dimplot_",unified_subcluster_ident,"_",prefix,"_highlight_subclusters_20241212.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, cells.highlight = unified_cells_to_highlight, cols.highlight = varibow_pal, sizes.highlight=0.1, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
for (i in 1:length(ordered_labels_for_plotting)) {
    unified_subcluster_label <- ordered_labels_for_plotting[i]
    print(unified_subcluster_label)
    unified_cells_to_highlight_sep <- unified_cells_to_highlight[[unified_subcluster_label]]
    col_highlight <- varibow_pal[i]
    print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, cells.highlight = unified_cells_to_highlight_sep, cols.highlight = col_highlight, sizes.highlight=0.1, label=TRUE, label.size=4, raster=FALSE) + labs(title = paste0("highlight ",unified_subcluster_label)), layers='Point', dpi=300))
}
dev.off()
print(table(obj$cell_type_individual[obj$merge.sub.19 %in% subcluster_idents], obj$sample_ID[obj$merge.sub.19 %in% subcluster_idents]))
print(table(obj$cell_type_individual[obj@meta.data[,c(unified_subcluster_ident)] %in% merge_unified_subcluster_label],
            obj@meta.data[,unified_subcluster_ident][obj@meta.data[,c(unified_subcluster_ident)] %in% merge_unified_subcluster_label]))
print(table(obj@meta.data[,unified_subcluster_ident]))
#      0      1     10     12     13     14     15     16     17     18     19 
# 255560 115936  31171  21689  19933  17772  14709  12988  10848   8639   8456 
#      2     20     21     22     23     24     25     26     27     28     29 
# 111752   8055   7528   7071   5079   5029   4430   2052   2015   1982   1722 
#     30     31     32     33     34      4      5      6      7      8      9 
#   1572   1529    495     65     19  63052  52907  47539  44723  36987  31289 
#  sub_0  sub_1 sub_10 sub_11 sub_12  sub_2  sub_3  sub_4  sub_5  sub_6  sub_7 
#  23359  16043    885    459    236  14988   9600   9063   7659   7081   6642 
#  sub_8  sub_9 
#   1338    951 
for (sample_ID in sample_vector) {
    tmp_df <- obj@meta.data[(obj@meta.data[,"sample_ID"] == sample_ID),c('original_Xenium_barcode',unified_subcluster_ident)]
    colnames(tmp_df) <- c('cell_id', 'group')
    write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",sample_ID,"_",prefix,"_",unified_subcluster_ident,".csv"),sep=",",quote=F,row.names=F)
}
table(obj$sample_ID[obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting], obj@meta.data[,unified_subcluster_ident][obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting])
#                     sub_0 sub_1 sub_10 sub_11 sub_12 sub_2 sub_3 sub_4 sub_5
# HT060P1-S1R1Fp1U1   10865  8486    864      3    177  5297  2656   128   212
# HT061P1-S1P1A1L1U1     42   454      1      0      0    89    10     0     0
# HT061P1-S1P1A1L4U1   1917  1206     11    431      2  1162  3950   185  2673
# HT125P1-S1H4A1L1U1     37  3126      0      0      5   302    94     0     0
# HT125P1-S1H8A1U1       11    17      0      0      2    35    31     0     1
# HT185P1-S1H2L1U1     5198  1644      2      6     35  2191   695   958   860
# HT227P1-S1H1L1U1     4290   581      4      0     12  5325  2111    24   200
# HT242P1-S1H4L4U1        0     6      0      0      0    48     4     0     0
# HT270P1-S1H1A1US2_1   999   523      3     19      3   539    49  7768  3713
#                     sub_6 sub_7 sub_8 sub_9
# HT060P1-S1R1Fp1U1    1640    63   146    17
# HT061P1-S1P1A1L1U1     74     0     0     1
# HT061P1-S1P1A1L4U1   1681  5951   235   329
# HT125P1-S1H4A1L1U1    179     0     0     0
# HT125P1-S1H8A1U1      118     0     0     0
# HT185P1-S1H2L1U1     1215   299   584    41
# HT227P1-S1H1L1U1      233    41   366     3
# HT242P1-S1H4L4U1      240     1     1     0
# HT270P1-S1H1A1US2_1  1701   287     6   560
# sub_0 Duct_like_1 CFTR+ FXYD2+ PROX1+ CXCL6- (very low and ~20%)
# sub_1 Ductal_reactive
# sub_2 Duct_like_1 CFTR+ FXYD2+ PROX1+ CXCL6- (very low)
# sub_3 Duct_like_2 CFTR+ FXYD2+ PROX1+ CXCL6+ (CXCL6+ is higher)
# sub_4 Acinar AMY2A+ GATM+ AQP8+
# sub_5 Duct_like_1 CFTR+ FXYD2+ PROX1+ AMY2A+ GATM+ These are ductal cells directly adjacent to acinar cells (maybe centro-acinar?)
# sub_6 iCAF #This cluster has low expression of Acinar markers however looking at the H&E histology many appear to be adjacent to the actual acinar cells. Most are PDGFRA+ PDGFRB+ COL5A2+ VCAN+ FBN1+. Those that are not are PTPRC+. Since the majority are expressing CAF markers and are C7+ and PTGDS+ I'm going to call them iCAFs
# sub_7 Acinar
# sub_8 Duct_like_1 CFTR+ FXYD2+ PROX1+ also GCG+ but upon looking at the spatial location of these cells it seems they are just normal ductal cells adjacent to islets.
# sub_9 Endothelial VWF+ PECAM1+ CD34+ STC1+ AMY2A+ GATM+ (low). these are endothelial cells in the middle of the exocrine lobules. The nuclear expansion model for cell segmentation is resulting in background expression.
# sub_10 Duct_like_2 CFTR+ FXYD2+ CXCL6+ PROX1+ (lower than sub_0 and sub_5), PPY+ but this population is strictly adjacent to PPY+ endocrine cells and almost exclusively comes from 1 sample.
# sub_11 CD8+T_cell AMY2A+ GATM+ AQP8+ (low) PTPRC+ (90-100% on dot plot) CD3D+ CD3E+ CD8A+ CD4+ CD8A is higher thand CD4+. Spatially this looks like T cells mixed in the exocrine (acinar/ductal) lobules 
# sub_12 Ductal_reactive
# Regarding the label for sub_10, this cluster predominantly comes from sample HT060P1-S1R1Fp1U1 which,
# has a number of abnormal shaped structures that are all strongly positive for PPY. 
# sub_10 often but not always is positive for both PPY and duct_like_1/duct_like_2 markers. It is not positive
# for reactive or ampullary ductal markers ACE2, GPX2, TFF2, SPDEF and has moderate expression 
# of EHF which is more highly expressed in the reactive ductal clusters. This paper https://doi.org/10.1002/path.6295
# talks about a mouse model in which aberant expression of SV40 large T antigen (TAg) in PPY islet 
# cells results in the formation of PDAC-like tumors. This group 
# https://doi.org/10.1053/j.gastro.2024.07.016 reported an increase in the PPY+ endocrine cells in the
# context of PDAC and they suspected transdifferentiation of GCG+ endocrine cells into hormone-negative 
# endocrine cells and then PPY+ cells. Upon close inspection of the histology of the sub_10 cells and the cells
# around them I am inclined to think that these are just ducts that are in close proximity to the PPY islet cells
# that make up most of cluster 10. As a result I am going to label them as duct_like_2 based on the increased
# expression of CXCL6 alongside CFTR and FXYD2. PROX1 is also expressed here but it is not as strong as sub_0 or sub_5.
# I originally wasn't going to try to subtype the islet cells but given the interesting PPY+ endocrine expansion we see in HT060P1 I'm going to do so.
Idents(obj) <- "seurat_clusters"
plan(sequential)
subcluster_ident <- "merge.sub.6.10.res_0.3"
subcluster_idents <- c('6','10')
obj <- FindSubCluster(obj, cluster = subcluster_idents, graph.name = "SCT_snn", subcluster.name = subcluster_ident, resolution = 0.3, algorithm = 1) #based on the recommendation here https://github.com/satijalab/seurat/issues/6682 by the authors I used the SCT_snn instead of SCT_nn graph
cells_to_highlight <- list()
merge_idents_highlight <- c()
for (new_ident in subcluster_idents) {
    merge_idents_highlight_single <- unique(obj@meta.data[,subcluster_ident])[grepl(paste0(new_ident,"_"),unique(obj@meta.data[,subcluster_ident]))]
    merge_idents_highlight <- unique(c(merge_idents_highlight, merge_idents_highlight_single))
}
unified_subcluster_ident <- paste0("unified.",subcluster_ident)
obj@meta.data[,unified_subcluster_ident] <- obj@meta.data[,subcluster_ident]
unified_cells_to_highlight <- list()
merge_unified_subcluster_label <- c()
for (sub_ident in merge_idents_highlight) {
    print(sub_ident)
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,subcluster_ident] == sub_ident),"barcode"]
    cells_to_highlight[[sub_ident]] <- barcodes_in_sub
    new_subcluster <- unlist(strsplit(sub_ident, "_"))[length(unlist(strsplit(sub_ident, "_")))]
    new_subcluster_label <- paste0("sub_",new_subcluster)
    obj@meta.data[,unified_subcluster_ident][obj@meta.data[,subcluster_ident] == sub_ident] <- new_subcluster_label
    merge_unified_subcluster_label <- unique(c(merge_unified_subcluster_label, new_subcluster_label))
}
print(unique(obj@meta.data[,unified_subcluster_ident]))
#  [1] "sub_12" "0"      "sub_2"  "2"      "9"      "sub_10" "sub_3"  "5"                                                                                             
#  [9] "12"     "sub_5"  "sub_11" "1"      "8"      "3"      "18"     "27"                                                                                            
# [17] "20"     "sub_8"  "sub_1"  "13"     "21"     "24"     "19"     "31"                                                                                            
# [25] "11"     "25"     "4"      "sub_4"  "17"     "14"     "sub_9"  "sub_6"                                                                                         
# [33] "sub_0"  "sub_13" "sub_7"  "15"     "30"     "28"     "16"     "7"                                                                                             
# [41] "22"     "29"     "26"     "23"     "sub_14" "32"     "33"     "34"
unified_cells_to_highlight <- list()
for (unified_subcluster_label in merge_unified_subcluster_label) {
    print(unified_subcluster_label)
    barcodes_in_sub <- obj@meta.data[(obj@meta.data[,unified_subcluster_ident] == unified_subcluster_label),"barcode"]
    unified_cells_to_highlight[[unified_subcluster_label]] <- barcodes_in_sub
    print(summary(obj@meta.data[,"nCount_Xenium"][obj@meta.data[,unified_subcluster_ident] == unified_subcluster_label]))
}
# [1] "sub_12"                                                                                                                                                        
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 7.00   57.75   89.00   95.21  123.25  339.00
# [1] "sub_2"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.0    55.0    93.0   102.8   138.0   565.0                                                                                                                     
# [1] "sub_10"                                                                                                                                                        
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 5.00   48.00   75.00   86.45  112.00  401.00                                                                                                                     
# [1] "sub_11"                                                                                                                                                        
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 1.00    9.00   15.00   23.19   25.00  262.00                                                                                                                     
# [1] "sub_3"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 8.0    67.0    94.0   101.1   127.0   499.0                                                                                                                     
# [1] "sub_5"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 7.0    74.0   106.0   113.9   143.0   514.0                                                                                                                     
# [1] "sub_1"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 11.0    75.0   106.0   111.9   142.0   429.0                                                                                                                     
# [1] "sub_4"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 1.00   10.00   20.00   40.63   54.00  575.00
# [1] "sub_9"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 1.00    7.00   11.00   22.17   19.00  366.00                                                                                                                     
# [1] "sub_6"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 1.00   12.00   24.00   42.35   55.00  507.00                                                                                                                     
# [1] "sub_0"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 7.0    64.0    96.0   104.9   134.0   619.0                                                                                                                     
# [1] "sub_8"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.0    21.0    59.0    75.1   110.2   423.0                                                                                                                     
# [1] "sub_13"                                                                                                                                                        
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 6.0    59.0    95.0   108.1   143.0   398.0                                                                                                                     
# [1] "sub_14"                                                                                                                                                        
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 9.00   49.00   76.00   80.98  108.00  376.00                                                                                                                     
# [1] "sub_7"                                                                                                                                                         
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.                                                                                                                     
# 2.00   20.00   46.00   61.96   86.00  475.00
Idents(obj) <- unified_subcluster_ident
tmp_df <- obj@meta.data[,c('barcode',subcluster_ident)]
colnames(tmp_df) <- c('barcode',subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
tmp_df <- obj@meta.data[,c('barcode',unified_subcluster_ident)]
colnames(tmp_df) <- c('barcode',unified_subcluster_ident)
write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_",unified_subcluster_ident,".tsv"),sep="\t",quote=F,row.names=F)
# Rscript $austin/tools/markers/dotplot_new_metadata_PDAC_xen_377_snv_v20241113.R \
# /diskmnt/Projects/Users/austins2/tools/markers/Cell_state_markers_PDAC_xen_377_hmulti_with_SNV_v20241113.txt \
# ./marker_dotplots/ ./PDAC_merge_primary_KRAS_20241205_single_SCT.rds \
# /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/cell_barcode_csv_files/PDAC_merge_primary_KRAS_20241205_single_SCT_unified.merge.sub.6.10.res_0.3.tsv \
# unified.merge.sub.6.10.res_0.3
print("plots")
library(scCustomize)
varibow_pal <- DiscretePalette_scCustomize(num_colors = length(unified_cells_to_highlight), palette = "varibow")
ordered_labels_for_plotting <- c()
for (i in 0:(length(merge_unified_subcluster_label)-1)) {
    unified_subcluster_label <- paste0("sub_",i)
    ordered_labels_for_plotting <- c(ordered_labels_for_plotting, unified_subcluster_label)
}
unified_cells_to_highlight <- unified_cells_to_highlight[ordered_labels_for_plotting]
pdf(paste0(out_dir,"/","Dimplot_",unified_subcluster_ident,"_",prefix,"_highlight_subclusters_20241212.pdf"),useDingbats = F, width=8, height=6)
print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, cells.highlight = unified_cells_to_highlight, cols.highlight = varibow_pal, sizes.highlight=0.1, label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
for (i in 1:length(ordered_labels_for_plotting)) {
    unified_subcluster_label <- ordered_labels_for_plotting[i]
    print(unified_subcluster_label)
    unified_cells_to_highlight_sep <- unified_cells_to_highlight[[unified_subcluster_label]]
    col_highlight <- varibow_pal[i]
    print(rasterize(DimPlot(obj, reduction = "umap.50PC", group.by = unified_subcluster_ident, cells.highlight = unified_cells_to_highlight_sep, cols.highlight = col_highlight, sizes.highlight=0.1, label=TRUE, label.size=4, raster=FALSE) + labs(title = paste0("highlight ",unified_subcluster_label)), layers='Point', dpi=300))
}
dev.off()
print(table(obj$cell_type_individual[obj$merge.sub.19 %in% subcluster_idents], obj$sample_ID[obj$merge.sub.19 %in% subcluster_idents]))
print(table(obj$cell_type_individual[obj@meta.data[,c(unified_subcluster_ident)] %in% merge_unified_subcluster_label],
            obj@meta.data[,unified_subcluster_ident][obj@meta.data[,c(unified_subcluster_ident)] %in% merge_unified_subcluster_label]))
print(table(obj@meta.data[,unified_subcluster_ident]))
#      0      1     11     12     13     14     15     16     17     18     19 
# 255560 115936  27929  21689  19933  17772  14709  12988  10848   8639   8456 
#      2     20     21     22     23     24     25     26     27     28     29 
# 111752   8055   7528   7071   5079   5029   4430   2052   2015   1982   1722 
#      3     30     31     32     33     34      4      5      7      8      9 
#  70375   1572   1529    495     65     19  63052  52907  44723  36987  31289 
#  sub_0  sub_1 sub_10 sub_11 sub_12 sub_13 sub_14  sub_2  sub_3  sub_4  sub_5 
#  13356  11604   3167   2888   1960   1829   1043  10246   7906   5190   5078 
#  sub_6  sub_7  sub_8  sub_9 
#   3889   3687   3584   3283
for (sample_ID in sample_vector) {                                                                                                                                
    tmp_df <- obj@meta.data[(obj@meta.data[,"sample_ID"] == sample_ID),c('original_Xenium_barcode',unified_subcluster_ident)]                                       
    colnames(tmp_df) <- c('cell_id', 'group')                                                                                                                       
    write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",sample_ID,"_",prefix,"_",unified_subcluster_ident,".csv"),sep=",",quote=F,row.names=F)      
}
table(obj$sample_ID[obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting], obj@meta.data[,unified_subcluster_ident][obj@meta.data[,unified_subcluster_ident] %in% ordered_labels_for_plotting])
plan(multicore, workers = 40)
degs <- FindAllMarkers(obj, assay='SCT', slot='data', test.use='wilcox', logfc.threshold = 0, min.pct=0.1, min.diff.pct=0.1)
plan(sequential)
write.table(degs,paste0(out_dir,"/","FindAllMarkers_",unified_subcluster_ident,"_",prefix,"_20241212.tsv"),sep="\t",quote=FALSE)
system(paste0("Rscript /diskmnt/Projects/Users/austins2/tools/markers/Shared_marker_comparison_from_seurat_degs.R -i ",paste0(out_dir,"/","FindAllMarkers_",unified_subcluster_ident,"_",prefix,"_20241212.tsv")))
#                     sub_0 sub_1 sub_10 sub_11 sub_12 sub_13 sub_14 sub_2
# HT060P1-S1R1Fp1U1    2760  1288   2890   2046   1134    250    180  8888
# HT061P1-S1P1A1L1U1     42     2    116    179     56      8      9   273
# HT061P1-S1P1A1L4U1   1651  1475      0      9      3    218     11    11
# HT125P1-S1H4A1L1U1     19    19     16    157     40      2      1   150
# HT125P1-S1H8A1U1      132   205    126    381    696     47     31   855
# HT185P1-S1H2L1U1     1799   534      3     49      3    191    192     7
# HT227P1-S1H1L1U1      915   492      4     14      2    104    113     7
# HT242P1-S1H4L4U1     5254  6320      7     19      3    750    458    11
# HT270P1-S1H1A1US2_1   784  1269      5     34     23    259     48    44
#                     sub_3 sub_4 sub_5 sub_6 sub_7 sub_8 sub_9
# HT060P1-S1R1Fp1U1    6645  1019  3593  1471   537   960   815 # I'm not sure why this sample has such a large expansion of PPY+ cells.
# HT061P1-S1P1A1L1U1    397   110   142    98    43    21    48 
# HT061P1-S1P1A1L4U1     46   591    57    74   375   415    69 
# HT125P1-S1H4A1L1U1     80    32    86    43     5    58    56 
# HT125P1-S1H8A1U1      699    53  1121   204    14   179   118 
# HT185P1-S1H2L1U1        1  1184     8   837   717   691  1213 
# HT227P1-S1H1L1U1        7   692     1   223   324   208   277 
# HT242P1-S1H4L4U1        9  1159    30   671  1312   728   468 
# HT270P1-S1H1A1US2_1    22   350    40   268   360   324   219 
print(unified_subcluster_ident)
# unified.merge.sub.6.10.res_0.3
# sub_0 Islet_INS+GCG+
# sub_1 Islet_INS+GCG+
# sub_2 iCAF 
# sub_3 Islet_gamma # mostly from HT060P1-S1R1Fp1U1
# sub_4 Islet_alpha
# sub_5 Islet_beta
# sub_6 Islet_delta
# sub_7 Islet_alpha
# sub_8 Islet_beta
# sub_9 Islet_beta
# sub_10 Islet_gamma # mostly from HT060P1-S1R1Fp1U1
# sub_11 Islet_gamma # mostly from HT060P1-S1R1Fp1U1
# sub_12 CD8+T_cell CD3E+ CD3D+ PTPRC+ CD8A+ CD4+ (lower than CD8A) also PPY+. This is predominantly from sample HT060P1-S1R1Fp1U1. these are T cells adjacent to the PPY+ endocrine cells.
# sub_13 Endothelial VWF+ PECAM1+ these are endothelial cells adjacent to islets
# sub_14 Islet_delta
write.table(obj@meta.data, paste0(out_dir,"/PDAC_merge_primary_KRAS_single_SCT_meta.data_20241216.tsv"),sep="\t",quote=FALSE)
saveRDS(obj, paste0(out_dir,"/",prefix,"_20241216.rds"))
# Unifying all of the cell type annotations
# conda activate seurat5
# cd /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/
library(Seurat)
library(tidyverse)
library(ggrastr)
set.seed(1234)
out_dir = "/diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/snvs_project/merged_analysis/v3_freeze/merge/primary_only/"
obj <- readRDS("PDAC_merge_primary_KRAS_20241205_single_SCT_20241216.rds")
prefix = "PDAC_merge_primary_KRAS_20241205_single_SCT"
meta.data <- read.table("PDAC_merge_primary_KRAS_single_SCT_meta.data_20241216.tsv",sep='\t',header=T)
# merge.sub.19 - c(19) #done
# merge.sub.13 - c(13) #done
# unified.merge.sub.2.20.19_1.19_5 #T_cell
# unified.merge.sub.1.14.25.19_0.19_4.19_6 #Myeloid
# unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8 # do not use
# unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8.Fsub_4 # do not use
# unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Msub_8.Fsub_5 # do not use
# unified.merge.sub.0.21.19_2.19_7.13_3.13_0 #Fibroblast
# unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Fsub_4
# unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Fsub_4.by.cell_type_individual
# unified.merge.sub.3.11.res_0.3 # ductal/acinar
# unified.merge.sub.6.10.res_0.3 #Islet
# 24 is myCAF
meta.data$cell_type_merge_v1 <- NA
meta.data$cell_type_merge_v1[ meta.data$seurat_clusters %in% c(17,22,23,27,28,29,31,32,33,34) | 
                                  meta.data$merge.sub.19 %in% c("19_3") |
                                  meta.data$merge.sub.13 %in% c("13_2","13_4","13_6","13_7","13_8","13_9","13_10","13_11","13_12") |
                                  meta.data$unified.merge.sub.2.20.19_1.19_5 %in% c("sub_13","sub_14") |
                                  meta.data$unified.merge.sub.1.14.25.19_0.19_4.19_6 %in% c("sub_13","sub_14","sub_15","sub_16","sub_17","sub_18","sub_19","sub_20") ] <- "LowCount"
meta.data$cell_type_merge_v1[ meta.data$seurat_clusters %in% c(18) ] <- "Glial"
meta.data$cell_type_merge_v1[ meta.data$seurat_clusters %in% c(8) ] <- "SMC"
meta.data$cell_type_merge_v1[ meta.data$seurat_clusters %in% c(9) ] <- "Plasma"
meta.data$cell_type_merge_v1[ meta.data$seurat_clusters %in% c(15) ] <- "B_cell"
meta.data$cell_type_merge_v1[ meta.data$seurat_clusters %in% c(12) ] <- "Mast"
meta.data$cell_type_merge_v1[ meta.data$seurat_clusters %in% c(5) ] <- "Endothelial"
meta.data$cell_type_merge_v1[ meta.data$seurat_clusters %in% c(30) ] <- "Adipocyte"
# clusters 4, 26 and 7 are almost entirely tumor. The cells in cluster 16 are PanIN. For these cells we will rely on the cell_type_individual labels.
meta.data$cell_type_merge_v1[ meta.data$seurat_clusters %in% c(4,7,26,16) ] <- meta.data$cell_type_individual[ meta.data$seurat_clusters %in% c(4,7,26,16) ]
CAF
meta.data$cell_type_merge_v1[ meta.data$merge.sub.19 %in% c("19_2") ] <- "myCAF_proliferative"
meta.data$cell_type_merge_v1[ meta.data$merge.sub.19 %in% c("19_7") ] <- "iCAF_proliferative"

meta.data$cell_type_merge_v1[ meta.data$merge.sub.13 %in% c("13_1") ] <- "Tumor_proliferative"
meta.data$cell_type_merge_v1[ meta.data$merge.sub.13 %in% c("13_5") ] <- "Tumor"

meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.2.20.19_1.19_5 %in% c("sub_0","sub_7","sub_12") ] <- "CD8+T_cell" 
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.2.20.19_1.19_5 %in% c("sub_1","sub_3","sub_5","sub_6","sub_10","sub_11") ] <- "CD4+T_cell"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.2.20.19_1.19_5 %in% c("sub_2") ] <- "Treg"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.2.20.19_1.19_5 %in% c("sub_4") ] <- "NK_cell"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.2.20.19_1.19_5 %in% c("sub_9") ] <- "Naive_T_cell"
meta.data$cell_type_merge_v1[ meta.data$merge.sub.19 %in% c("19_1") ] <- "CD4+T_cell_proliferative"
meta.data$cell_type_merge_v1[ meta.data$merge.sub.19 %in% c("19_5") ] <- "CD8+T_cell_proliferative"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.2.20.19_1.19_5 %in% c("sub_8") & !(meta.data$merge.sub.19 %in% c("19_1","19_5")) ] <- "CD8+T_cell_proliferative"

meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.1.14.25.19_0.19_4.19_6 %in% c("sub_3") ] <- "Dendritic_cell"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.1.14.25.19_0.19_4.19_6 %in% c("sub_9") ] <- "Macrophage_proliferative"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.1.14.25.19_0.19_4.19_6 %in% c("sub_0","sub_1","sub_2","sub_4","sub_5","sub_6","sub_7","sub_8","sub_10","sub_11","sub_12") ] <- "Macrophage"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.1.14.25.19_0.19_4.19_6 %in% c("sub_21") ] <- "iCAF"

meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.0.21.19_2.19_7.13_3.13_0 %in% c("sub_0","sub_4","sub_6","sub_8") ] <- "iCAF" # subclustering in unified.merge.sub.0.21.19_2.19_7.13_3.13_0.Fsub_4 found that sub_4 was mostly iCAFs
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.0.21.19_2.19_7.13_3.13_0 %in% c("sub_1","sub_2","sub_3","sub_5","sub_7") | meta.data$seurat_clusters %in% c(24)] <- "myCAF"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.0.21.19_2.19_7.13_3.13_0 %in% c("sub_9") ] <- "CAF_proliferative"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.3.11.res_0.3 %in% c("sub_0","sub_2","sub_5","sub_8") ] <- "Duct_like_1"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.3.11.res_0.3 %in% c("sub_1","sub_12") ] <- "Ductal_reactive"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.3.11.res_0.3 %in% c("sub_3","sub_10") ] <- "Duct_like_2"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.3.11.res_0.3 %in% c("sub_4","sub_7") ] <- "Acinar"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.3.11.res_0.3 %in% c("sub_6") ] <- "iCAF"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.3.11.res_0.3 %in% c("sub_9") ] <- "Endothelial"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.3.11.res_0.3 %in% c("sub_11") ] <- "CD8+T_cell"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.6.10.res_0.3 %in% c("sub_0","sub_1") ] <- "Islet_INS+GCG+"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.6.10.res_0.3 %in% c("sub_2") ] <- "iCAF"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.6.10.res_0.3 %in% c("sub_3","sub_10","sub_11") ] <- "Islet_gamma"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.6.10.res_0.3 %in% c("sub_4","sub_7") ] <- "Islet_alpha"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.6.10.res_0.3 %in% c("sub_5","sub_8","sub_9") ] <- "Islet_beta"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.6.10.res_0.3 %in% c("sub_12") ] <- "CD8+T_cell"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.6.10.res_0.3 %in% c("sub_13") ] <- "Endothelial"
meta.data$cell_type_merge_v1[ meta.data$unified.merge.sub.6.10.res_0.3 %in% c("sub_6","sub_14") ] <- "Islet_delta"

print(unique(meta.data$seurat_clusters[is.na(meta.data$cell_type_merge_v1)]))
# integer(0)
print(unique(meta.data$cell_type_merge_v1))
#                   Acinar                Adipocyte                   B_cell                                                                                               
#                    15715                     1572                    14720                                                                                               
#                      CAF        CAF_proliferative              CD4+_T_cell                                                                                               
#                     2672                     1502                       15                                                                                               
#               CD4+T_cell CD4+T_cell_proliferative          CD8+_T_cells_NK                                                                                               
#                    51571                     1429                        7                                                                                               
#               CD8+T_cell CD8+T_cell_proliferative    CD8+T_near_epithelial                                                                                               
#                    40094                      630                        7                                                                                               
#       CXCL10+_Macrophage     CXCL10+CXCL9+Myeloid                       DC                                                                                               
#                      317                       16                       98                                                                                               
#           Dendritic_cell              Duct_like_1              Duct_like_2                                                                                               
#                    14762                    47557                    10585                                                                                               
#          Ductal_reactive              Endothelial               Fibroblast                                                                                               
#                    16279                    55901                      114                                                                                               
#                    Glial                     iCAF     Immune_proliferative                                                                                               
#                     8675                   131469                       64
#                    Islet              Islet_alpha               Islet_beta 
#                      151                     8877                    11945 
#              Islet_delta              Islet_gamma           Islet_INS+GCG+ 
#                     4973                    13971                    24960 
#                 LowCount         LowCount_Myeloid               Macrophage 
#                    43608                      204                   120975 
# Macrophage_proliferative                     Mast                    myCAF 
#                     3754                    21828                   163531 
#                  Myeloid             Naive_T_cell                       NK 
#                     2520                     1418                      155 
#                  NK_cell                    PanIN                   Plasma 
#                     9157                     9850                    32004 
#            reactive_duct            Reactive_duct                      SMC 
#                      208                      376                    37410 
#                   T_cell     T_cell_proliferative                T_cell/NK 
#                      551                        8                       36 
#                T_NK_cell                     Treg                    Tumor 
#                        1                    19719                    91072 
#      Tumor_proliferative      Tumor_Proliferative                  Unknown 
#                     5748                     5817                     2222 
#   Vascular_smooth_muscle                     vSMC 
#                       10                       67
# fixing the labeling to be consistent since I included some labels from the individual objects for clusters 4, 26, 7, and 
write.table(meta.data[,c('barcode','cell_type_merge_v1')], paste0(out_dir,"/","cell_barcode_csv_files","/",prefix,"_","cell_type_merge",".tsv"),sep="\t",quote=F,row.names=F)
meta.data$cell_type_merge_v2 <- meta.data$cell_type_merge_v1
meta.data$cell_type_merge_v2[ meta.data$cell_type_merge_v1 %in% c("CD4+_T_cell")] <- "CD4+T_cell"
meta.data$cell_type_merge_v2[ meta.data$cell_type_merge_v1 %in% c("T_NK_cell","T_cell/NK","CD8+_T_cells_NK","CD8+T_near_epithelial","T_cell")] <- "CD8+T_cell"
meta.data$cell_type_merge_v2[ meta.data$cell_type_merge_v1 %in% c("CXCL10+_Macrophage","CXCL10+CXCL9+Myeloid")] <- "Macrophage"
meta.data$cell_type_merge_v2[ meta.data$cell_type_merge_v1 %in% c("DC")] <- "Dendritic_cell"
meta.data$cell_type_merge_v2[ meta.data$cell_type_merge_v1 %in% c("Immune_proliferative","T_cell_proliferative")] <- "CD4+T_cell_proliferative"
meta.data$cell_type_merge_v2[ meta.data$cell_type_merge_v1 %in% c("Islet")] <- "Islet_gamma"
meta.data$cell_type_merge_v2[ meta.data$cell_type_merge_v1 %in% c("LowCount_Myeloid")] <- "LowCount"
meta.data$cell_type_merge_v2[ meta.data$cell_type_merge_v1 %in% c("Myeloid")] <- "Macrophage"
meta.data$cell_type_merge_v2[ meta.data$cell_type_merge_v1 %in% c("NK")] <- "NK_cell"
meta.data$cell_type_merge_v2[ meta.data$cell_type_merge_v1 %in% c("Reactive_duct","reactive_duct")] <- "Ductal_reactive"
meta.data$cell_type_merge_v2[ meta.data$cell_type_merge_v1 %in% c("Tumor_Proliferative")] <- "Tumor_proliferative"
meta.data$cell_type_merge_v2[ meta.data$cell_type_merge_v1 %in% c("Unknown")] <- "LowCount"
meta.data$cell_type_merge_v2[ meta.data$cell_type_merge_v1 %in% c("NK")] <- "NK_cell"
meta.data$cell_type_merge_v2[ meta.data$cell_type_merge_v1 %in% c("Vascular_smooth_muscle","vSMC")] <- "SMC"
meta.data$cell_type_merge_v2[ meta.data$cell_type_merge_v1 %in% c("CAF","Fibroblast")] <- "myCAF"
meta.data$cell_type_merge_v2[ meta.data$cell_type_merge_v1 %in% c("CAF_proliferative")] <- "myCAF_proliferative"
# 'Acinar'
# 'Adipocyte'
# 'B_cell'
# 'CAF' <- 'myCAF'
# 'CAF_proliferative' -> "myCAF_proliferative"
# 'CD4+_T_cell' -> 'CD4+T_cell' #
# 'CD4+T_cell'
# 'CD4+T_cell_proliferative'
# 'CD8+_T_cells_NK' -> 'CD8+T_cells' #
# 'CD8+T_cell'
# 'CD8+T_cell_proliferative'
# 'CD8+T_near_epithelial' -> 'CD8+T_cell'# 
# 'CXCL10+_Macrophage' -> 'Macrophage'
summary(meta.data$nFeature_Xenium[meta.data$cell_type_merge_v1 %in% c('CXCL10+_Macrophage')])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5.00   21.00   29.00   30.47   38.00   80.00 
# 'CXCL10+CXCL9+Myeloid' -> 'Macrophage'
summary(meta.data$nFeature_Xenium[meta.data$cell_type_merge_v1 %in% c('CXCL10+CXCL9+Myeloid')])
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 15.00   21.00   26.00   29.88   34.50   72.00
# 'DC' -> 'Dendritic_cell'
# 'Dendritic_cell'
# 'Duct_like_1'
# 'Duct_like_2'
# 'Ductal_reactive'
# 'Endothelial'
# 'Fibroblast' -> 
# 'Glial'
# 'iCAF'
# 'Immune_proliferative' -> 'CD4+T_cell_proliferative'
# 'Islet' -> 'Islet_gamma'
# 'Islet_alpha'
# 'Islet_beta'
# 'Islet_delta'
# 'Islet_gamma'
# 'Islet_INS+GCG+'
# 'LowCount'
# 'LowCount_Myeloid' -> 'LowCount'
summary(meta.data$nFeature_Xenium[meta.data$cell_type_merge_v1 %in% c('LowCount_Myeloid')])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00    7.00   10.50   14.39   19.00   63.00
# 'Macrophage'
# 'Macrophage_proliferative'
# 'Mast'
# 'myCAF'
# 'Myeloid' -> "Macrophage"
# 'Naive_T_cell'
# 'NK' -> 'NK_cell'
# 'NK_cell'
# 'PanIN'
# 'Plasma'
# 'Reactive_duct' -> 'Ductal_reactive'
# 'reactive_duct' -> 'Ductal_reactive'
# 'SMC'
# 'T_cell' -> 'CD8+T_cell'
# 'T_cell_proliferative' -> 'CD4+T_cell_proliferative'
# 'T_cell/NK' -> 'CD8+T_cell'
# 'T_NK_cell' -> 'CD8+T_cell'
# 'Treg'
# 'Tumor'
# 'Tumor_proliferative'
# 'Tumor_Proliferative' -> 'Tumor_proliferative'
# 'Unknown' -> "LowCount"
# 'Vascular_smooth_muscle' -> 'SMC'
# 'vSMC' -> 'SMC'
unique(meta.data$cell_type_merge_v2)
# [1] "CD8+T_cell"               "iCAF"                    
# [3] "myCAF"                    "Treg"                    
# [5] "Plasma"                   "CD4+T_cell"              
# [7] "Islet_gamma"              "Endothelial"             
# [9] "Mast"                     "Islet_beta"              
# [11] "Macrophage"               "SMC"                     
# [13] "Duct_like_2"              "Dendritic_cell"          
# [15] "Glial"                    "Duct_like_1"             
# [17] "LowCount"                 "NK_cell"                 
# [19] "Ductal_reactive"          "Islet_INS+GCG+"          
# [21] "Macrophage_proliferative" "myCAF_proliferative"     
# [23] "Acinar"                   "Islet_alpha"             
# [25] "Islet_delta"              "B_cell"                  
# [27] "Naive_T_cell"             "Adipocyte"               
# [29] "CD8+T_cell_proliferative" "CD4+T_cell_proliferative"
# [31] "Tumor_proliferative"      "Tumor"                   
# [33] "PanIN"
meta.data$cell_type_merge_v3 <- meta.data$cell_type_merge_v2
meta.data$cell_type_merge_v3[meta.data$cell_type_merge_v2 %in% c("CD8+T_cell_proliferative")] <- "CD8+T_cell"
meta.data$cell_type_merge_v3[meta.data$cell_type_merge_v2 %in% c("CD4+T_cell_proliferative")] <- "CD4+T_cell"
meta.data$cell_type_merge_v3[meta.data$cell_type_merge_v2 %in% c("Macrophage_proliferative")] <- "Macrophage"
meta.data$cell_type_merge_v3[meta.data$cell_type_merge_v2 %in% c('myCAF_proliferative')] <- "myCAF"
meta.data$cell_type_merge_v3[meta.data$cell_type_merge_v2 %in% c('Tumor_proliferative')] <- "Tumor"
unique(meta.data$cell_type_merge_v3)
#  [1] "CD8+T_cell"      "iCAF"            "myCAF"           "Treg"           
#  [5] "Plasma"          "CD4+T_cell"      "Islet_gamma"     "Endothelial"    
#  [9] "Mast"            "Islet_beta"      "Macrophage"      "SMC"            
# [13] "Duct_like_2"     "Dendritic_cell"  "Glial"           "Duct_like_1"    
# [17] "LowCount"        "NK_cell"         "Ductal_reactive" "Islet_INS+GCG+" 
# [21] "Acinar"          "Islet_alpha"     "Islet_delta"     "B_cell"         
# [25] "Naive_T_cell"    "Adipocyte"       "Tumor"           "PanIN"
tmp_df <- meta.data[,c("cell_type_merge_v1","cell_type_merge_v2","cell_type_merge_v3")]
obj <- AddMetaData(obj, metadata = tmp_df)
print(all(table(obj$cell_type_merge_v2) == table(meta.data$cell_type_merge_v2)))
print(all(table(obj$cell_type_merge_v3) == table(meta.data$cell_type_merge_v3)))
sample_vector <- unique(obj$sample_ID)
for (sample_ID in sample_vector) {                                                                                                                                
    tmp_df <- obj@meta.data[(obj@meta.data[,"sample_ID"] == sample_ID),c('original_Xenium_barcode',"cell_type_merge_v3")]                                       
    colnames(tmp_df) <- c('cell_id', 'group')                                                                                                                       
    write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",sample_ID,"_",prefix,"_","cell_type_merge_v3",".csv"),sep=",",quote=F,row.names=F)
    write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",sample_ID,"_",prefix,"_","cell_type_merge_v3",".tsv"),sep="\t",quote=F,row.names=F)
    tmp_df <- obj@meta.data[(obj@meta.data[,"sample_ID"] == sample_ID),c('original_Xenium_barcode',"cell_type_merge_v2")]                                       
    colnames(tmp_df) <- c('cell_id', 'group')                                                                                                                       
    write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",sample_ID,"_",prefix,"_","cell_type_merge_v2",".csv"),sep=",",quote=F,row.names=F)
    write.table(tmp_df, paste0(out_dir,"/","cell_barcode_csv_files","/",sample_ID,"_",prefix,"_","cell_type_merge_v2",".tsv"),sep="\t",quote=F,row.names=F)
}
write.table(obj@meta.data, paste0(out_dir,"/PDAC_merge_primary_KRAS_single_SCT_meta.data_20240102.tsv"),sep="\t",quote=FALSE)
saveRDS(obj, paste0(out_dir,"/",prefix,"_20240102.rds"))
# cell_type_merge_v3 were the cell types used in the analysis of the paper.