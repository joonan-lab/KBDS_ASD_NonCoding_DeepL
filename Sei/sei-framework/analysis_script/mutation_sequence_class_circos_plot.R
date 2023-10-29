require(data.table)
library(shades)
require(ggrepel)
require(plyr)
library(scales)


## load plotdata
plotdata = fread("/data2/sei-framework/result/Korean_ASD_WGS/hg38/data/Korean_ASD_WGS_proband_DNMs_strong_prediction_distribution_plotdata.tsv")
plotdata$group_label = factor(plotdata$group_label,levels=c('E','CTCF','P','TF','HET','PC','TN','L'))


## load plotdata_group and assign types for columns
plotdata_group = fread("/data2/sei-framework/result/Korean_ASD_WGS/hg38/data/Korean_ASD_WGS_proband_DNMs_strong_prediction_distribution_plotdata_group.tsv")
plotdata_group$index = factor(plotdata_group$index, levels=plotdata_group$index)
plotdata_group$group_label = factor(plotdata_group$group_label,levels=c('E','CTCF','P','TF','HET','PC','TN','L'))


## color list
scale_fill_Publication <- function(...){
      #discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#984ea3","#662506","#a6cee3","#fb9a99","#dddddd")), ...)
      discrete_scale("fill","Publication",manual_pal(values = c("#984ef3","#fdb462","#ef3b2c","#386cb0","#662506","#a6cee3","#fb9a99","#dddddd")), ...)
}
scale_color_Publication <- function(...){
      #discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#984ea3","#662506","#a6cee3","#fb9a99","#dddddd")), ...)
      discrete_scale("color","Publication",manual_pal(values = brightness(c("#984ef3","#fdb462","#ef3b2c","#386cb0","#662506","#a6cee3","#fb9a99","#dddddd"),0.6)), ...)
}

max_score = max(abs(plotdata$seqclass_pred_score))
bar_w = ceiling(max_score)+3
x_min = (-1) * (bar_w+2)
x_max = bar_w+2

## Set arguments for the plot
options(repr.plot.width = 15, repr.plot.height = 15, repr.plot.res = 200)

## Draw and save the plot
## 수정해야되는 부분
# 1. second geom_bar의 width -> maximum score * 2
# 2. xlim( , ) -> maximum score의 abs+2의 -, +
# 3. 마지막 geom_text의 x -> xlim의 max 값
# 4. 마지막 geom_text의 size -> 3
# 5. geom_point의 size -> abs(seqclass_pred_score) ** 0.7
ggplot() +
  geom_bar(aes(x = 0, y = Freq, group = index, fill=group_label), width=2.2, stat = "identity", color = "white", data=plotdata_group) +
  geom_bar(aes(x = 0, y = Freq, group = index, fill=group_label), width=bar_w*2, alpha=I(0.1),stat="identity", color = "white", data=plotdata_group) +
  coord_polar(theta = "y", start = 0)+
  theme_void()+
  xlim(x_min, x_max)+
geom_point(aes(y=newind, x=-seqclass_pred_score,color=group_label,size=I(abs(seqclass_pred_score)**0.7)), data=plotdata[abs(seqclass_pred_score)>0.,])+
geom_text_repel(aes(y=newind, x=-(seqclass_pred_score),label=""), nudge_x = -0.3, data=plotdata[seqclass_pred_score>1.3,],size=4.,force=6)+
geom_text_repel(aes(y=newind, x=-(seqclass_pred_score),label=""), nudge_x = 0.3, data=plotdata[seqclass_pred_score< -1.3,],size=4.,force=6)+
geom_text(aes(x = 0, y = rev(cumsum(rev(Freq)))-Freq/2, label = label), size=3, alpha=0.8, color = "white", fontface = "bold", data=plotdata_group)+
geom_text(aes(x = x_max, y = rev(cumsum(rev(Freq)))-Freq/2, label = names), size=5, alpha=0.8, color = "black", fontface = "bold", data=plotdata_group)+
scale_fill_Publication()+scale_color_Publication()

ggsave('/data2/sei-framework/result/Korean_ASD_WGS/hg38/figures/Korean_ASD_WGS_proband_DNMs_strong_prediction_distribution.pdf', width=40, height = 20, dpi = 200, device=cairo_pdf)