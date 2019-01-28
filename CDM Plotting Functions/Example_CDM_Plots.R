#### Example CDM Plots
## Created by: Reuben Biel
## Created on: May 31, 2018
## Last Modified: May 31, 2018

#### Load Functions
source('CDM_plotting_functions.R')

#### Directories
parent_dir <- getwd() ## change as necessary
sub_dir <- 'Example_Batch_CDM_Simulation' ## change as necessary
dir <- paste(parent_dir, sub_dir, sep='/')


#### 3D animation ####
p <- CDM.plotly(dir = paste(dir, 'model_iter4', sep='/'))
p
htmlwidgets::saveWidget(as_widget(p), file.path(parent_dir, "CDM_animation.html"))

CDM.plotly.timepoint(dir = paste(dir, 'model_iter1', sep='/'), timepoint='100000', xrange=c(1,100))

#### figures ####
## Timeseries Plots
# Max Dune Height
batch_CDM_df<-batch.CDM.df(dir, pattern = '^h.*')
n<-names(batch_CDM_df)
p1 <- batch.CDM.summary.plot(batch_CDM_df, n[6], ylab = "Maximum Dune Elevation (m)")
p1

# Veg Cover
batch_CDM_df<-batch.CDM.df(dir, pattern = '^veget_x.*')
n<-names(batch_CDM_df)
p2 <- batch.CDM.summary.plot(batch_CDM_df, n[6], ylab = expression("Maximum Proportional Cover"))
p2

## Final Timepoint Plots
# Dune Height
batch_timepoint_df <- batch.timepoint.df(dir, timepoint = '100000', pattern = "^h.")
n<-names(batch_timepoint_df)
p3 <- batch.CDM.timepoint.plot(batch_timepoint_df, n[6])
p3

# Veg Cover
batch_timepoint_df <- batch.timepoint.df(dir, timepoint = '100000', pattern = "^veget_x.")
n<-names(batch_timepoint_df)
p4 <- batch.CDM.timepoint.plot(batch_timepoint_df, n[6], ylab = "Vegetation Proportional Cover")
p4

p5 <- grid.arrange(p1, p2, p3, p4, nrow = 2)
p5

ggsave(file.path(dir, "h.timeseries.png"), 
       plot = p1, # or give ggplot object name as in myPlot,
       width = 5, height = 4, 
       units = "in", # other options c("in", "cm", "mm"), 
       dpi = 200, scale = 1.5)
ggsave(file.path(dir, "veget_x.timeseries.png"), 
       plot = p2, # or give ggplot object name as in myPlot,
       width = 5, height = 4, 
       units = "in", # other options c("in", "cm", "mm"), 
       dpi = 200, scale = 1.5)
ggsave(file.path(dir, "h.final.png"), 
       plot = p3, # or give ggplot object name as in myPlot,
       width = 5, height = 4, 
       units = "in", # other options c("in", "cm", "mm"), 
       dpi = 200, scale = 1.5)
ggsave(file.path(dir, "veget_x.final.png"), 
       plot = p4, # or give ggplot object name as in myPlot,
       width = 5, height = 4, 
       units = "in", # other options c("in", "cm", "mm"), 
       dpi = 200, scale = 1.5)
ggsave(file.path(dir, "summary_plots.png"), 
       plot = p5, # or give ggplot object name as in myPlot,
       width = 10, height = 8, 
       units = "in", # other options c("in", "cm", "mm"), 
       dpi = 200, scale = 1.5)

save(p1, p2, p3, p4, file=file.path(dir, "plots.rda"))

