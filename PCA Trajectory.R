# PCA visualization
# Yongzhe Wang

Dt <- read.csv('Simulated Biomarker Dataset.csv')
head(Dt)

# Extract columns for RNA-seq
X <- Dt[, -c(1:3)]
# Perform PCA 
pca <- prcomp(x = X, scale = T)
# Extract first 3 PCs for each point in X
X_pca <- pca$x[, c('PC1', 'PC2', 'PC3')]
# Combine the value of rotated data with time and treatment info
Dt.Plot <- cbind(Dt[, 1:3], X_pca)
head(Dt.Plot)

require(ggplot2)
p_2D_PCA <- 
  ggplot(Dt.Plot, aes(x = PC1, y = PC2, col = Treatment)) + 
  geom_point(alpha = 0.8) + 
  theme_bw() +                                                      # dark-on-light theme
  theme(panel.border = element_blank(),
        panel.background = element_blank(),                    
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(),                               # these two are for the axis line
        axis.line.y = element_line(),
        axis.text.x = element_text(colour = "black"),               # there two are for texts in axes
        axis.text.y = element_text(colour = "black"),
        axis.ticks.x = element_line(),                              # these two are for ticks in axes
        axis.ticks.y = element_line(),
        axis.title.x = element_text(colour = "black", face = 'bold', vjust = -1),                              
        axis.title.y = element_text(colour = "black", face = 'bold'),
        legend.title = element_text(colour = "black", face = 'bold'),
        legend.text = element_text(colour = "black"))  
p_2D_PCA

p_2D_PCA_ID <- 
  ggplot(Dt.Plot, aes(x = PC1, y = PC2, col = Treatment)) + 
  geom_point(alpha = 0.8) + 
  facet_grid(cols = vars(Patient_ID)) + 
  theme_bw() +                                                      # dark-on-light theme
  theme(panel.border = element_blank(),
        panel.background = element_blank(),                    
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(),                               # these two are for the axis line
        axis.line.y = element_line(),
        axis.text.x = element_text(colour = "black"),               # there two are for texts in axes
        axis.text.y = element_text(colour = "black"),
        axis.ticks.x = element_line(),                              # these two are for ticks in axes
        axis.ticks.y = element_line(),
        axis.title.x = element_text(colour = "black", face = 'bold', vjust = -1),                              
        axis.title.y = element_text(colour = "black", face = 'bold'),
        legend.title = element_text(colour = "black", face = 'bold'),
        legend.text = element_text(colour = "black"))  
p_2D_PCA_ID

require(plotly)
p_3D_PCA <- plot_ly(data = Dt.Plot,                                         # dataset
                    x = ~ PC1, y = ~ PC2, z = ~ PC3,                        # specify x/y/z axes
                    color = ~ Treatment, colors = c('#D6604D', '#4393C3'),  # colors
                    marker = list(size = 5))                                # point size
p_3D_PCA <- p_3D_PCA %>% add_markers()                                      # specify scatter plot
p_3D_PCA <- p_3D_PCA %>% layout(scene = list(xaxis = list(title = 'PC1'),   # axes title
                                             yaxis = list(title = 'PC2'), 
                                             zaxis = list(title = 'PC3')))
p_3D_PCA


require(tidyr)
require(dplyr)
# Aggregated 
Dt.Plot.Aggregated <- Dt.Plot %>% group_by(Time) %>% summarise(PC1_Aggregated = mean(PC1), 
                                                               PC2_Aggregated = mean(PC2), 
                                                               PC3_Aggregated = mean(PC3))
Dt.Plot.Aggregated$'Treatment' <- factor(c(rep('No', 12), rep('Yes', 12)), 
                                         levels = c('No', 'Yes'))
head(Dt.Plot.Aggregated)

# Interprolated points
Time <- Dt.Plot.Aggregated$Time
Interprolate <- seq(from = min(Time), to = max(Time), len = 200)
# Spline regression on each PC dimension
PC1 <- splinefun(Time, Dt.Plot.Aggregated$PC1_Aggregated)(Interprolate)
PC2 <- splinefun(Time, Dt.Plot.Aggregated$PC2_Aggregated)(Interprolate)
PC3 <- splinefun(Time, Dt.Plot.Aggregated$PC3_Aggregated)(Interprolate)
# Combine with treatment variable
Treatment <- factor(ifelse(Interprolate <= 12, 'No', 'Yes'),  # 1-12 month: no treatment, 13-24 month: treatment
                    levels = c('No', 'Yes'))
Dt.Plot.Traj <- data.frame(Treatment, 'PC1' = PC1, 'PC2' = PC2, 'PC3' = PC3)
head(Dt.Plot.Traj)

p_2D_PCA_Traj <- 
  ggplot(Dt.Plot.Traj, aes(x = PC1, y = PC2, col = Treatment)) + 
  geom_path(arrow = arrow(angle = 10, 
                          ends = "last",                            # arrow at the end of a curve
                          type = "closed")) +                       # closed triangle for arrow head
  theme_bw() +                                                      # dark-on-light theme
  theme(axis.line.x = element_line(),                               # these two are for the axis line
        axis.line.y = element_line(),
        axis.text.x = element_text(colour = "black"),               # there two are for texts in axes
        axis.text.y = element_text(colour = "black"),
        axis.ticks.x = element_line(),                              # these two are for ticks in axes
        axis.ticks.y = element_line(),
        axis.title.x = element_text(colour = "black", face = 'bold', vjust = -1),                              
        axis.title.y = element_text(colour = "black", face = 'bold'),
        legend.title = element_text(colour = "black", face = 'bold'),
        legend.text = element_text(colour = "black"))  
p_2D_PCA_Traj

p_3D_PCA_Traj <- plot_ly(Dt.Plot.Traj, 
                         x = ~ PC1, y = ~ PC2, z = ~ PC3, 
                         color = ~ Treatment, colors = c('#D6604D', '#4393C3'), 
                         type = 'scatter3d',                                         # trace is 'scatter3d'
                         mode = 'lines',                                             # show lines
                         line = list(width = 4))                                     # Control line width
p_3D_PCA_Traj <- p_3D_PCA_Traj %>% layout(scene = list(xaxis = list(title = 'PC1'), 
                                                       yaxis = list(title = 'PC2'), 
                                                       zaxis = list(title = 'PC3')))
p_3D_PCA_Traj

id <- unique(Dt.Plot$Patient_ID)
n_id <- length(id)
Dt.Plot.Subject <- data.frame()
for(i in 1:n_id){
  Loca <- which(Dt.Plot$Patient_ID == id[i])
  Trans <- Dt.Plot[Loca, ]
  # Interprolated points
  Time <- Trans$Time
  Interprolate <- seq(from = min(Time), to = max(Time), len = 200)
  # Spline regression on each PC dimension
  PC1 <- splinefun(Time, Trans$PC1)(Interprolate)
  PC2 <- splinefun(Time, Trans$PC2)(Interprolate)
  PC3 <- splinefun(Time, Trans$PC3)(Interprolate)
  # Combine treatment variable
  Treatment <- factor(ifelse(Interprolate <= 12, 'No', 'Yes'),  # 1-12 month: no treatment, 13-24 month: treatment
                      levels = c('No', 'Yes'))
  # Result
  Result <- data.frame('Patient_ID' = rep(id[i], 200), 
                       'Treatment'= Treatment,
                       'PC1' = PC1, 
                       'PC2' = PC2, 
                       'PC3' = PC3)
  Dt.Plot.Subject <- rbind(Dt.Plot.Subject, Result)
}

p_3D_PCA_Traj <- plot_ly(Dt.Plot.Subject, 
                         x = ~ PC1, y = ~ PC2, z = ~ PC3, 
                         color = ~ Patient_ID, 
                         type = 'scatter3d',                                         # trace is 'scatter3d'
                         mode = 'lines',                                             # show lines
                         line = list(width = 4))                                     # Control line width
p_3D_PCA_Traj <- p_3D_PCA_Traj %>% layout(scene = list(xaxis = list(title = 'PC1'), 
                                                       yaxis = list(title = 'PC2'), 
                                                       zaxis = list(title = 'PC3')))
p_3D_PCA_Traj

