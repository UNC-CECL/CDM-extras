#### CDM Plotting Functions
## Created by: Reuben Biel
## Created on: May 23, 2018
## Last Modified: May 31, 2018

## Load libraries
library(plotly)
library(ggplot2)
library(viridis)
library(plyr)
library(gridExtra)

#### CDM Plotly 3D Surface Animation ####
## creates a list of data frames for each timepoint in the time series
## for time series with positive rates of shoreline change, this function pads out the dataframe to create a spatially stationary time series
timeseries.stationary.list <- function(file.list, metadata) {
  for(ii in 1:length(file.list)) {
    # read elevation data
    dat<-read.table(file.list[ii])
    if(ii==1)
      xshore_transect_ct=ncol(dat)
    names(dat)<-paste0('xshore_trans',1:xshore_transect_ct)
    time=as.numeric(strsplit(file.list[ii], split="[.]")[[1]][2])
    
    # reposition by shoreline change distance
    x_range=range(metadata[,'scr'])
    dat_reposition=as.data.frame(matrix(0, ncol=ncol(dat), nrow=nrow(dat)+abs(diff(x_range))))
    names(dat_reposition)=names(dat)
    if (time>max(metadata[,'iter'])) {
      x_offset = metadata[which.max(metadata[,'iter']),'scr']
      time_yr = metadata[which.max(metadata[,'iter']),'time_yr']
    } else {
      x_offset = metadata[time == metadata[,'iter'],'scr']
      time_yr = metadata[time == metadata[,'iter'],'time_yr']
    }
    
    dat_reposition[(x_offset-x_range[1]+1):(nrow(dat_reposition)+x_offset-x_range[2]),]=dat
    
    # list of dataframes for each timepoint
    if (ii==1) {
      dat_complete<-list()
      timepts<-c()
      timepts_yr<-c()
    }
    dat_complete[[ii]] <- dat_reposition
    timepts[ii]<-time
    timepts_yr[ii]<-time_yr
    
    
    rm(dat, dat_reposition)
  }
  
  names(dat_complete)<-formatC( round( timepts_yr, 2 ), format='f', digits=2 )
  dat_complete[order(timepts)] # sorted by timepoints
}


## function to create colormap scalar within a single dataframe
wsv_scaling<-function(h_dat, v_dat, water1 = 0.0, water2 = 0.5) {
  v=v_dat
  h=h_dat
  xmax=nrow(h_dat)
  ymax=ncol(h_dat)
  
  for (j in 1:ymax) {
    idx<-which(h_dat[,j] < water1)
    idx2<-which(h_dat[,j] >= water1 & h[,j] < water2)
    v[idx,j] <- -1
    v[idx2,j] <- -1*(h[idx2,j] - water2)/(water1 - water2)
  }
  v
}

## function to create list of colormap scalars
wsv_scaling_list<-function(elev_complete, veget_complete, water1 = 0.0, water2=0.5) {
  wsv_scaling_complete <- elev_complete
  for(ii in seq_along(wsv_scaling_complete)) {
    wsv_scaling_complete[[ii]] <- wsv_scaling(elev_complete[[ii]], veget_complete[[ii]], water1, water2)
  }
  wsv_scaling_complete
}


## plotly animation
CDM.plotly <- function(dir = getwd()) {
  data_dir = paste0(dir, '/DATA')
  ## load metadata
  metadata=read.table(file.path(data_dir, 'time.dat'), skip=9)
  names(metadata)<-c('iter', 'time_yr', 'max_ht', 'max_cover', 'sand_volume', 'scr') 
  iter.year=metadata[nrow(metadata),'iter']/metadata[nrow(metadata),'time_yr']
  
  # elevation list
  file.list<-list.files(path=data_dir, pattern = '^h.*', full.names = T)
  elev_complete<-timeseries.stationary.list(file.list = file.list, metadata = metadata)
  
  # vegetation list
  file.list2<-list.files(path=data_dir, pattern = '^veget.x.*', full.names=T)
  veget_complete<-timeseries.stationary.list(file.list = file.list2, metadata = metadata)
  
  # colormap scalar list
  wsv_scaling_complete<-wsv_scaling_list(elev_complete, veget_complete)
  
  # colormap
  cm1<-colorRamp(c('blue','tan'), bias=1)( seq(0,1,length.out=50) )
  cm2<-colorRamp(c('tan','forestgreen'), bias=2)( seq(0,1,length.out=50) )
  colormap<-rbind(cm1,cm2)
  
  ## plotly animation
  z = elev_complete
  
  for (ii in seq_along(z)){
    if(ii==1) {
      p <- plot_ly(showscale = FALSE)
      r= range(z)
    }
    
    p <- p %>%
      layout(scene = list(aspectmode = 'manual', aspectratio = list(x=0.4, y=2, z=0.2),
                          camera = list(eye = list(x = 1.25, y = -0.5, z = 0.5)),
                          xaxis=list(title='alongshore location (m)'),
                          yaxis=list(title='cross-shore location (m)'),
                          zaxis=list(title='elevation (m)', range = r))) %>% 
      add_surface(as.matrix(z[[ii]]), frame=as.numeric(names(z)[ii]), surfacecolor=as.matrix(wsv_scaling_complete[[ii]]), 
                  colors=colormap, showscale=F, cauto = F, cmin = -1, cmax = 1) %>%
      animation_opts(1000, redraw = FALSE)
  }
  
  p <- p %>% 
    animation_button( x = 1, xanchor = "right", y = 0, yanchor = "bottom") %>%
    animation_slider( currentvalue = list(prefix = "Year ", font = list(color="red")) )
  
  # display animation
  p
}


### Example ####
## set directories
# parent_dir <- paste0('/home/rgbiel/Documents/Enhanced_CDM/3_Stochastic_Seeding/')
# sub_dir='Comparison2'
# dir = paste0(parent_dir,sub_dir,'/model_iter2')
# setwd(dir)
# 
# p <- CDM.plotly(dir = dir)
# p
# htmlwidgets::saveWidget(as_widget(p), file.path(parent_dir, "graph.html"))


#### ggplot2 ####
plot.CDM <- function(dat_complete, ylab, p) {
  l1<-do.call('rbind',lapply(dat_complete, function(x) {apply(x, 2, max)}))
  l2<-data.frame(year = as.numeric(row.names(l1)),median = apply(l1,1, median),
                 lower = apply(l1,1, function(x) quantile(x, probs=0.25)),
                 upper = apply(l1,1, function(x) quantile(x, probs=0.75)))
  
  p <- ggplot(l2, aes(x=year, y=median)) + 
    scale_x_continuous(breaks = pretty(l2$year, n = 5)) +
    scale_y_continuous(breaks = pretty(l2$upper, n = 5)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill = "grey70") +
    geom_line(size = 1) + 
    xlab("Year") + ylab(ylab) + 
    theme_light(base_size = 18)  +
    scale_color_viridis(option="plasma",discrete=TRUE) + 
    scale_fill_viridis(option="plasma",discrete=TRUE)

  p
}

## plot max vegetation cover over time
plot.vegmax <- function(dat_complete, ylab = "Maximum Proportional Vegetation Cover") {
  plot.CDM(dat_complete, ylab)
}

## plot max dune height over time
plot.htmax <- function(dat_complete, ylab = "Maximum Dune Elevation") {
  plot.CDM(dat_complete, ylab)
}

# plot.vegmax(veget_complete)

##### Batch Summary Plots ####
# dir1


## sim_dir summary df function
CDM.summary.df <- function(sim_dir, pattern = '^h.*', iter_metadata) {
  
  data_dir = paste0(sim_dir, '/DATA')
  file.list<-list.files(path=data_dir, pattern = pattern, full.names = T)
  metadata=read.table(file.path(data_dir,'time.dat'), skip=9)
  names(metadata)<-c('iter', 'time_yr', 'max_ht', 'max_cover', 'sand_volume', 'scr') 
  
  dat_complete<-timeseries.stationary.list(file.list = file.list, metadata = metadata)
  
  l1<-do.call('rbind',lapply(dat_complete, function(x) {apply(x, 2, max)}))
  l2<-data.frame(model_iter = as.numeric(gsub('model_iter', "", basename(sim_dir))), 
                 year = as.numeric(row.names(l1)),median = apply(l1,1, median),
                 lower = apply(l1,1, function(x) quantile(x, probs=0.25)),
                 upper = apply(l1,1, function(x) quantile(x, probs=0.75)))
  l3 <- merge(l2, iter_metadata, by="model_iter")
  l3
}


#for (sim_dir in sim_dirs)
batch.CDM.df <- function(dir, pattern = '^h.*') {
  # find all simulation folders
  folder_list=list.files(path=dir, pattern='model_iter')
  sim_dirs = paste(dir, folder_list, sep='/')
  
  # load iter metadata
  iter_metadata=read.table(file.path(dir,'iter_metadata.txt'), sep="\t", header=T)
  
  batch_CDM_df <- ldply(sim_dirs, CDM.summary.df, pattern, iter_metadata)
  
  # return
  batch_CDM_df
}

`+.uneval` <- function(a,b) { ## from: https://stackoverflow.com/questions/28777626/how-do-i-combine-aes-and-aes-string-options
  `class<-`(modifyList(a,b), "uneval")
}

batch.CDM.summary.plot <- function(batch_CDM_df, variable, ylab = "Maximum Dune Elevation (m)") {
  batch_CDM_df[[variable]] <- as.factor(batch_CDM_df[[variable]])
  
  ggplot(batch_CDM_df, aes_string(x="year", y="median", color=variable)) + 
    scale_x_continuous(breaks = pretty(batch_CDM_df$year, n = 5)) +
    scale_y_continuous(breaks = pretty(batch_CDM_df$upper, n = 5)) +
    geom_ribbon(aes(ymin=lower, ymax=upper,linetype=NA) + aes_string(fill=variable), alpha = 0.3) +
    geom_line(size = 1) + 
    xlab("Year") + ylab(ylab) + 
    theme_light(base_size = 18) +
    theme(legend.position="bottom", legend.box = "horizontal")  +
    scale_color_viridis(option="plasma",discrete=TRUE) + 
    scale_fill_viridis(option="plasma",discrete=TRUE)
}


### example code
# parent_dir <- paste0('/home/rgbiel/Documents/Enhanced_CDM/3_Stochastic_Seeding/')
# sub_dir='Comparison1'
# dir <- paste0(parent_dir, sub_dir)
# batch_CDM_df<-batch.CDM.df(dir, pattern = '^h.*') 
# names(batch_CDM_df)
# batch.CDM.summary.plot(subset(batch_CDM_df, subset=replicate==1, ylab=), probabilistic.propagule.pressure, ylab = "Maximum Dune Elevation (m)")
# 
# batch_CDM_df<-batch.CDM.df(dir, pattern = '^veget_x.*') 
# batch.CDM.summary.plot(subset(batch_CDM_df, subset=replicate==1), probabilistic.propagule.pressure, ylab = expression("Maximum Proportional Cover"))


#### Batch Final Timepoint Plots ####
## creates a df for CDM simulation at timepoint
CDM.timepoint.df <- function(sim_dir, timepoint = '170000', pattern = '^h.', iter_metadata) {
  
  data_dir = paste0(sim_dir, '/DATA')
  file.list<-list.files(path=data_dir, pattern = paste0(pattern, timepoint), full.names = T)
  metadata=read.table(file.path(data_dir,'time.dat'), skip=9)
  names(metadata)<-c('iter', 'time_yr', 'max_ht', 'max_cover', 'sand_volume', 'scr') 
  
  dat_complete<-read.table(file.list)
  
  l2<-data.frame(model_iter = as.numeric(gsub('model_iter', "", basename(sim_dir))), xshore_pos = 1:nrow(dat_complete),median = apply(dat_complete,1, median),
                 lower = apply(dat_complete,1, function(x) quantile(x, probs=0.25)),
                 upper = apply(dat_complete,1, function(x) quantile(x, probs=0.75)))
  l3 <- merge(l2, iter_metadata, by="model_iter")
  l3
}

## creates data frame for all CDM simulations in dir at timepoint
batch.timepoint.df <- function(dir, timepoint = '170000', pattern = '^h.*') {
  # find all simulation folders
  folder_list=list.files(path=dir, pattern='model_iter')
  sim_dirs = paste(dir, folder_list, sep='/')
  
  # load iter metadata
  iter_metadata=read.table(file.path(dir,'iter_metadata.txt'), sep="\t", header=T)
  
  batch_timepoint_df <- ldply(sim_dirs, CDM.timepoint.df, timepoint, pattern, iter_metadata)
  
  # return
  batch_timepoint_df
}

## plots data from batch.timepoint.df
batch.CDM.timepoint.plot <- function(batch_timepoint_df, variable, ylab = "Elevation (m)") {
  batch_timepoint_df[[variable]] <- as.factor(batch_timepoint_df[[variable]])
  
  ggplot(batch_timepoint_df, aes_string(x="xshore_pos", y="median", color=variable)) + 
    scale_x_continuous(breaks = pretty(batch_timepoint_df$xshore_pos, n = 5)) +
    scale_y_continuous(breaks = pretty(batch_timepoint_df$upper, n = 5)) +
    geom_ribbon(aes(ymin=lower, ymax=upper,linetype=NA) + aes_string(fill=variable), alpha = 0.3) +
    geom_line(size = 1) + 
    xlab("Cross-Shore Position (m)") + ylab(ylab) + 
    theme_light(base_size = 18) +
    theme(legend.position="bottom", legend.box = "horizontal") +
    scale_color_viridis(option="plasma",discrete=TRUE) + 
    scale_fill_viridis(option="plasma",discrete=TRUE)
}

#### Example ####
# parent_dir <- paste0('/home/rgbiel/Documents/Enhanced_CDM/3_Stochastic_Seeding/')
# sub_dir='Comparison3'
# dir <- paste0(parent_dir, sub_dir)
# 
# batch_timepoint_df <- batch.timepoint.df(dir, timepoint = '170000')
# names(batch_timepoint_df)
# batch.CDM.timepoint.plot(subset(batch_timepoint_df, subset=replicate==1), SCR)
# 
# 
# batch_timepoint_df <- batch.timepoint.df(dir, timepoint = '170000', pattern = "^veget_x.")
# names(batch_timepoint_df)
# batch.CDM.timepoint.plot(subset(batch_timepoint_df, subset=replicate==1), SCR, ylab = "Vegetation Proportional Cover")
