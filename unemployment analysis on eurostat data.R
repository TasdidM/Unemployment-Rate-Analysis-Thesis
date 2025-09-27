library(tidyverse)
library(ggplot2)
library(zoo)
library(fda)
library(funFEM)
library(giscoR)
library(countrycode)
library(sf)

# read data (collected from eurostat)
unemploy_rate <- read.csv("unemployment_in_europe.csv")
head(unemploy_rate)
nrow(unemploy_rate)
# 30 countries under observation
ncol(unemploy_rate)
# 133 - 1(first column) = 132 data points for individual observations

# missing data check
anyNA(unemploy_rate) # no data are missing

# changing time format
# removes x in front of dates
time_clean <- sub("^X", "", colnames(unemploy_rate[1,-1]))
# transform dates into month-year format
ym_dates <- as.yearmon(time_clean, "%Y.%m")
# transform dates into year-month-date format
formatted_dates <- format(ym_dates, "%Y-%m-%d")
# changes column names to formatted date
cols <- c(colnames(unemploy_rate[1]), formatted_dates)
colnames(unemploy_rate) <- cols
head(unemploy_rate)

### data visualisation

# convert wide data to long data
long_data <- unemploy_rate %>% 
  pivot_longer(cols = -Reference_area,
               names_to = "Time_Period",
               values_to = "Unemployment")
# character dates to date type
long_data$Time_Period <- as.Date(long_data$Time_Period)

# Figure 1
# highest and lowest average unemployment value over the year 2013 to 2023

# average rates over 11 years observation period
avg_rate <- aggregate(Unemployment ~ Reference_area, long_data, mean)
# standard deviation of unemployment rates over 11 years observation period
sd_error <- aggregate(Unemployment ~ Reference_area, long_data, sd)
# create a data frame with reference name, mean value and standard error columns
data_avg_sd <- data.frame(Reference_area = avg_rate$Reference_area,
                          Average = avg_rate$Unemployment,
                          Sderror = sd_error$Unemployment)
# plot point graph with addition and subtraction of 1 standard error wicks
windows()
ggplot(data_avg_sd, aes(x = reorder(Reference_area, Average), 
                     y = Average)) +
  geom_point(size = 5, color = "dodgerblue") +
  geom_errorbar(aes(ymin = Average-Sderror, ymax = Average+Sderror),
                width = 0.2, size = 1.2, color = "red")+
  labs(title = "Mean and Standard Error of Unemployment Rates by Country",
       x = "Country", y = "Unemployment rate") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
        legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        axis.text.y = element_text(size = 25),
        plot.title = element_text(size = 26))

# higher average unemployment, a potential indicator of economic challenges in 
# those regions.
# lower average unemployment, possibly reflecting stronger or more stable labor markets.

#figure 2
#heatmap on rates over each year
# extracting years 
year_data <- long_data %>% 
  mutate(Year = format(Time_Period, "%Y"))

# mean rates over each years
mean_year_data <- aggregate(Unemployment ~ Reference_area + Year, year_data, mean) %>% 
  arrange(Reference_area)

# plotting heat map with each square representing each year
windows()
ggplot(mean_year_data, aes(x = reorder(Reference_area, Unemployment), 
                           y = Year, fill = Unemployment)) +
  geom_tile() +
  scale_fill_viridis_c(option = "cividis")+
  labs(title = "Mean Unemployment Rate Heatmap for Each Year",
       x = "Country", y = "Year", fill = "Unemployment Rate (%)") +
  theme_minimal()+
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    axis.text.x = element_text(size = 23, angle = 45, hjust =1),
    plot.title = element_text(size = 26),
    text = element_text(size = 25))


# Figure 3
# identify values that are significantly greater than others
overall_median <- median(long_data$Unemployment)
Q1 <- quantile(long_data$Unemployment, 0.25)
Q3 <- quantile(long_data$Unemployment, 0.75)
IQR <- Q3 - Q1
upper_fence <- Q3 + 1.5 * IQR
long_data$sig_high <- long_data$Unemployment > upper_fence

# unemployment rate over time line plot for each country
windows()
ggplot(long_data, aes(x = Time_Period, y = Unemployment, 
                      color = sig_high, group = Reference_area)) +
  labs(title = paste0("Unemployment Rates Over Time by Country (overall median = ",
                      overall_median, "%)"),
       x = "Time Period", y = "Unemployment Rates") +
  geom_smooth(method = lm, se = F, color = "gray20")+
  theme_minimal() + 
  geom_point(aes(color = sig_high), alpha = 0.5) +
  scale_colour_manual(values = c("FALSE" = "dodgerblue4", "TRUE" = "red")) +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(size = 14, angle = 45),
    plot.title = element_text(size = 26),
    text = element_text(size = 25))+
  facet_wrap(~Reference_area, scales = "free_y") +
  theme(legend.position = "none")


#figure 4
#rough unemployment rate over at starting and ending of the time frame
unemploy_trends <- long_data %>%
  group_by(Reference_area) %>%
  summarize(
    first_rate = first(Unemployment),
    last_rate = last(Unemployment),
    rate_change = last(Unemployment) - first(Unemployment)
  )
# Sort by rate change
unemploy_trends <- unemploy_trends %>%
  arrange(rate_change)

#improving and worsening countries
windows()
ggplot(unemploy_trends, aes(x = reorder(Reference_area, rate_change, decreasing = T), 
                        y = rate_change, fill = rate_change > 0)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "green"),
                    labels = c("Improved", "Worsened")) +
  coord_flip() +
  labs(
    title = "Improving and Worsening Countries in Unemployment Rates",
    x = "Country",
    y = "Rate Change"
  ) +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank())+
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 26),
    text = element_text(size = 25))




# Functional Data Analysis
ncol(unemploy_rate)
basis <- create.bspline.basis(1:132, nbasis = 134)
data_matrix <- unemploy_rate[,-1]

# function to select lambda value using GCV method
gcv_lambda_select <- function(data_matrix){
    #GCV to choose lambda for absolute data
    loglam <- seq(-3, 3, 0.25)
    nlam <- length(loglam)
    dfsave <- rep(NA, nlam)
    gcvsave <- rep(NA, nlam)
    
    # create linear differential operator of Lfd class
    lfd_unemploy <- int2Lfd(2) # second derivative as penalty term
    
    #function to check gcv for all provided log lambda
    for(ilam in 1:nlam){
      cat(paste("log10 Lambda =", loglam[ilam], "\n"))
      lambda = 10^loglam[ilam]
      fdParobj = fdPar(basis, lfd_unemploy, lambda)
      smoothlist = smooth.basis(1:132, data_matrix,
                                fdParobj)
      dfsave[ilam] = smoothlist$df
      gcvsave[ilam] = sum(smoothlist$gcv)
    }
    
    # store lambda which has minimum gcv value
    log_lambda <- loglam[which.min(gcvsave)]
    lambda <- 10^log_lambda
    df_lambda <- dfsave[which.min(gcvsave)]
    
    return(list(
      lambda = lambda,
      df_lambda = df_lambda,
      gcvsave = gcvsave,
      loglam = loglam,
      dfsave = dfsave,
      log_lambda = log_lambda,
      df_lambda = df_lambda
    ))
}

# plot of gcv and df curve over selected lambdas
gcv_lambda <- gcv_lambda_select(t(data_matrix))
windows()
par(mar = c(5,5,4,2),mfrow=c(1,2))
plot(gcv_lambda$loglam, gcv_lambda$gcvsave, type = "l", lwd = 2, 
     ylab = "GCV(lambda)", xlab = "log(lambda)",
     main = "GCV vs log(lambda)",
     cex.axis = 2, cex.lab = 2,cex.main = 3)
points(gcv_lambda$loglam, gcv_lambda$gcvsave, cex = 1, pch = 16)
abline(v = gcv_lambda$log_lambda, col = "red", lwd = 2)
plot(gcv_lambda$loglam, gcv_lambda$dfsave, type = "l", lwd = 2, 
     ylab = "Degree of Freedom", xlab = "log(lambda)",
     main = "df vs log(lambda)",
     cex.axis = 2, cex.lab = 2,cex.main = 3)
abline(v =gcv_lambda$log_lambda, h = gcv_lambda$df_lambda, col = "red", lwd = 2)
par(mfrow = c(1,1))

# start and end date for the given data
start_date <- as.Date("2013-01-01")
end_date <- as.Date("2023-12-01")
date_seq <- seq(start_date, end_date, by = "month")
form_date <- format(date_seq, "%b-%Y")

# smooth splines with the selected lambda (all curves)
windows()
par(mar = c(5,5,4,2))
lfd_unemploy <- int2Lfd(2)
fdpar_unemploy <- fdPar(basis, lfd_unemploy, gcv_lambda$lambda)
fd_unemploy <- smooth.basis(1:132, t(data_matrix), fdpar_unemploy)$fd
plot(fd_unemploy, ylab = "Unemployment Rate (in %)",
     xlab = "time(2013-Jan to 2023-Dec)",
     xaxt = "n", lwd = 2,
     main = "Smooth Curves for Unemployment Rates", 
     cex.axis = 2, cex.lab = 2, cex.main = 3)
abline(h = 11, lwd = 5, col = "skyblue")
text(80, 10.5, "Threshold = 11", col = "skyblue", cex = 1.5, pos = 4)
m_unemploy <- mean.fd(fd_unemploy)
lines(m_unemploy, col = 2, lwd = 8)
text(132, eval.fd(132, m_unemploy), "Mean", col = 2, cex = 1.5, pos = 4)
axis(1, cex.axis= 2, at = seq(1,132, by = 12), labels = form_date[c(1, seq(13,132, by = 12))])

# magnitude outliers detection
windows()
par(mar = c(5,5,4,2))
smooth_data_matrix <- eval.fd(1:132, fd_unemploy)
boxplot_fdata <- fbplot(smooth_data_matrix, 1:132, method = "MBD", xlab = "Time", 
                        ylab = "Unemployment Rate",
                        main = "Functional Boxplot of Absolute Values",
                        cex.axis = 2, cex.lab = 2, cex.main = 2)
#this plot shows curves that are further away from central region
outliers <- boxplot_fdata$outpoint
for(i in outliers){
  text(x = 70, y = data_matrix[i, 70],
       labels = unemploy_rate[i, 1], pos = 3, col = "red", cex = 1.5)
}
med_curve <- boxplot_fdata$medcurve
text(x = 70, y = data_matrix[med_curve, 70],
     labels = unemploy_rate[med_curve, 1], pos = 3, col = "black", cex = 1.5)


# Functional PCA
fpca <- pca.fd(fd_unemploy, nharm = 4, centerfns = TRUE)


# plot cumulative variance explained
windows()
par(mar = c(5,5,4,2), mfrow = c(1,2))
# plot the explained cumulative percentage of total variations
plot(cumsum(fpca$values[1:10])/sum(fpca$values),
     xlab = "Number of Components",
     ylab = "Cumulative Variance Explained",
     main = "Variance Explained by FPCs",
     cex.axis = 2, cex.lab = 2,
     cex.main = 3, pch = 16, cex = 2)
abline(h = 0.95)
# two pcs are enough to explain more than 95% variances

# plot FPC curves
harm <- fpca$harmonics
harmvals <- eval.fd(1:132, harm)
matplot(1:132, harmvals[,1:2], type = "l", col = 1:3, lty = 1, lwd = 2.5,
        ylab = "Value of PC components", xlab = "Time", xaxt = "n",
        main = "Functional Principal Components",
        cex.axis = 2, cex.lab = 2, cex.main = 3)
axis(1, cex.axis= 2, at = seq(1,132, by = 12), labels = form_date[c(1, seq(13,132, by = 12))])
abline(h=0, lty = 2, col = 1, lwd = 2)
legend("bottomright", c(paste("FPC1 (",round(fpca$varprop[1],2)*100, "%)"), 
                        paste("FPC2 (",round(fpca$varprop[2],2)*100, "%)")),
       col = 1:2, lty = 1, cex = 2, bty = "n", lwd = 2)
par(mfrow=c(1,1))
# fpc 1 explains countries that had higher unemployment rate but decreased over time
# fpc 2 positive value explains rates were high between 2016 and 2020 respect to other years in observation 
# and negative value explains higher unemployment rates at  2013 to 2016



#plot fpc scores
windows()
par(mar = c(5,5,4,2))
plot(fpca$scores[,1:2], cex = 0.01, xlab = "FPC 1",
     ylab = "FPC 2", cex.axis = 2, cex.lab = 2,
     cex.main = 3, main = "FPC Scores of Individual Countries")
abline(h = 0, v = 0, lwd = 2)
text(fpca$scores[,1], fpca$scores[,2], 
     labels = unemploy_rate$Reference_area, cex = 2)


# clustering using DFM model
set.seed(420)
femmodels <- c("DkBk", "DkB", "DBk",
               "DB", "AkjBk", "AkjB", "AkB", "AkBk", "AjBk", "AjB", "ABk",
               "AB")
nmodels <- length(femmodels)
femresults <- list()
bestk <- bestbic <- numeric(0)
K=2:10
fembic <- matrix(NA,nrow=nmodels,ncol=max(K))
for (i in 1:nmodels){
  try({print(femmodels[i])
      femresults[[i]] <- funFEM(fd_unemploy,model=femmodels[i],K=K)
      fembic[i,K] <- femresults[[i]]$allCriterions$bic
      bestk[i] <- which(fembic[i,]==max(fembic[i,K],na.rm=TRUE))
      bestbic[i] <- max(fembic[i,K],na.rm=TRUE)})
}
besti <- which(bestbic==max(bestbic,na.rm=TRUE))
# DFM model selected using BIC index
femmodels[besti] #AjBk
# number of clusters that returns higher BIC index for the chosen DFM model
bestk[besti] #8
femresult <- femresults[[besti]]

# print countries in clusters
for(i in 1:femresult$K){
  print(i)
  print(unemploy_rate[femresult$cls==i, 1])
}

# BIC plot for all models and K
windows()
par(mar = c(5,5,4,2), mfrow = c(1,2))
i <- 1
plot(1:max(K), fembic[i,], col = i, pch = i, 
     ylim = c(min(fembic, na.rm = T), max(fembic, na.rm = T)), type = "n",
     xlab = "Number of clusters", ylab = "BIC index",
     cex.axis = 2, cex.lab = 2, cex.main = 3)
for(i in 1:nmodels){
  text(1:max(K), fembic[i,], femmodels[i], col = i, cex = 2)
}

#Curve plot for each individual clusters
clmeans <- fd_unemploy
clmeans$coefs <- t(femresult$prms$my)
plot(clmeans, lwd = 2, xaxt="n", ylab = "",
     cex.axis = 2, cex.lab = 2, cex.main = 3,
     main = "Mean Curve of Individual Group")
axis(1, cex.axis= 2, at = seq(1,132, by = 12), 
     labels = form_date[c(1, seq(13,132, by = 12))])
legend(100,25,legend=1:8,col=c(1:6,1:2),lty=c(1:5,1:3),
       lwd = 2, cex = 1.5)


# Colour scheme for the clusters
cluster_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                    "#FF7F00", "#A65628", "#F781BF", "#00CED1")

# cluster on functional principal components
windows()
par(mar = c(5,5,4,2))
plot(fpca$scores[,1:2], cex = 0.01, xlab = "FPC 1",
     ylab = "FPC 2", cex.axis = 2, cex.lab = 2,
     cex.main = 3, main = "FPC Scores of Individual Countries")
abline(h = 0, v = 0, lwd = 2)
text(fpca$scores[,1], fpca$scores[,2], 
     labels = unemploy_rate$Reference_area, cex = 2,
     col = cluster_colors[femresult$cls])


# Mapping of clusters

# getting map polygon values
europe <- gisco_get_nuts(year = 2024, nuts_level = 0)
# extract country code that are relevant to analysis
eu_countries <- countrycode(unemploy_rate[,1], origin = "country.name", 
                            destination = 'iso2c',
                            custom_match = c("Greece" = "EL"))

# Data frame with Country code and the cluster that individual country belongs to
cluster_data <- data.frame(CNTR_CODE = eu_countries,
                           cluster = femresult$cls)

# Join cluster info to spatial data
europe_clusters <- europe %>%
  filter(CNTR_CODE %in% cluster_data$CNTR_CODE) %>%
  left_join(cluster_data, by = "CNTR_CODE")

# Plot map colored by cluster
windows()
ggplot(europe_clusters) +
  geom_sf(aes(fill = factor(cluster)), color = "black") +
  scale_fill_manual(values = cluster_colors, name = "Clusters") +
  coord_sf(xlim = c(-25,35), ylim = c(25, 73), expand = T) +
  theme_minimal() +
  labs(title = "European Countries Coloured by Clusters") +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        text = element_text(size = 20))










