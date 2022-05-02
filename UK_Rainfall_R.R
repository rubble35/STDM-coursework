# Install packages and libraries
install.packages("maptools")
install.packages("forecast")
install.packages("rnn")


library(sp)
library(maptools)
library(lattice)
library(spdep)
library(rgdal)
library(tmap)
library(ggplot2)
library(gridExtra)
library(scatterplot3d)
library(plot3D)
library(reshape)
library(rgl)
library(gstat)
library(OpenStreetMap)
library(raster)
library(spacetime)
library(forecast)
library(nnet)
library(caret)
library(rnn)




## 1 IMPORT THE DATA

# 1.1 Monthly UK Rainfall Data by Region

# Dataset being used is monthly rainfall records for the UK sourced from the UK Met Office.
# This data was recently updated as a result of the crowdsourcing project 'Rainfall Rescue'.
# The period covered is January 1836 to December 2021.

# Read geodatabase and set CRS using a proj4string
uk_rain_raw <- readOGR(dsn="Data/UK Rainfall.gdb", layer="uk_rain_all_districts", 
                        p4s = CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs")@projargs)

# Dataframe has 10 elements/features: 1 for each district region of the UK as defined by the MetOffice.
# Dataframe has 2232 columns/fields of data: one for each month for 186 years.

summary(uk_rain_raw)

# Remove the Shape_Area and Shape_Length fields from dataframe to avoid skewing rainfall data with these large values
uk_rain <- subset(uk_rain_raw, select = -c(Shape_Length, Shape_Area))

# Create a matrix and add row names
rain_matrix <- data.matrix(uk_rain@data[,-c(1)])
rownames(rain_matrix) <- uk_rain@data[,"REGION"]

# Visualise the data for one region by way of example
# Set parameters so that axis labels show
par(mar = c(5, 5, 4, 2) + 0.1)
plot(rain_matrix[1,], type="l", xaxt="n", xlab="Year (monthly data)", ylab="Monthly rainfall (mm)", main="Monthly Rainfall for Scotland N region, 1836 to 2021")
axis(1, at=seq(49, 2209, 240), labels=seq(1840, 2020, 20))
# Very long time series, so hard to see detail of any patterns or trends.
# Look at a shorter time period
plot(rain_matrix[1, 1969:2232], type="l", xaxt="n", xlab="Year (monthly data)", ylab="Monthly rainfall (mm)", main="Monthly Rainfall for Scotland N region, Jan 2000 to Dec 2021")
axis(1, at=seq(1, 264, 24), labels=seq(2000, 2020, 2))
# Would expect a seasonal trend, and it looks like there is one, but it is noisy
# Have a look at another region
plot(rain_matrix[6, 1969:2232], type="l", xaxt="n", xlab="Year (monthly data)", ylab="Monthly rainfall (mm)", main="Monthly Rainfall for East Anglia region, Jan 2000 to Dec 2021")
axis(1, at=seq(1, 264, 24), labels=seq(2000, 2020, 2))
# Also looks quite noisy, with significant variation between max and min and from year to year

# 1.2 Monthly UK Rainfall Records for Historic UK Weather Stations

# This dataset comprises historic monthly rainfall records for 37 individual weather stations across the UK sourced from the UK Met Office.
# The stations have latitude, longitude and altitude data, i.e. they are point locations rather than areas
# Data goes back to 1855 for some stations, but most stations have full records going back only to 1960's.
# Data imported for analysis is monthly for the period Sep 1964 to Aug 2016 (624 months), but just for 29 stations that have a full record in this timeframe.
uk_station <- read.csv("Data/rainfall_by_station.csv")

# Rename first column to remove strange characters at front
colnames(uk_station)[1] <- "Station"

# Create a matrix excluding the station name, lat, long and altitude columns
station_matrix <- data.matrix(uk_station[,-c(1:4)])

# Create a vector containing the station names and add to matrix
station <- as.vector(uk_station[,c(1)], mode="any")
rownames(station_matrix) <- station

# Visualise the data for one weather station by way of example
plot(station_matrix[1,], type="l", xaxt="n", xlab="Year (monthly data)", ylab="Rainfall (mm)", main="Average Monthly Rainfall for Aberporth, 1965 to 2015")
axis(1, at=seq(5, 605, 60), labels=seq(1965, 2015, 5))

# Also bring in file containing annual average rainfall for the same 29 weather stations for period 1965 to 2015
# And rename columns and create a data matrix
uk_station_annual <- read.csv("Data/station_annual_ave.csv")
colnames(uk_station_annual)[1] <- "Station"
colnames(uk_station_annual)[5:ncol(uk_station_annual)] <- as.character(c(1965:2015))
station_annual_matrix <- data.matrix(uk_station_annual[,-c(1:4)])
rownames(station_annual_matrix) <- station

# Visualise the data for one weather station by way of example
plot(station_annual_matrix[1,], type="l", xaxt="n", xlab="Year", ylab="Rainfall (mm)", main="Average Annual Rainfall for Aberporth, 1965 to 2015")
axis(1, at=seq(1, 51, 5), labels=seq(1965, 2015, 5))
# Very variable from one month to the next. Possibly seasonal as there are a series of peaks and troughs

# Create column averages and plot to see any trends in average annual rainfall for whole of UK
plot(colMeans(station_annual_matrix), type="l", xaxt="n", xlab="Year", ylab="Rainfall (mm)", main="Average Annual Rainfall for all UK Weather Stations, 1965 to 2015")
axis(1, at=seq(1, 51, 5), labels=seq(1965, 2015, 5))
# No obvious long-term trend in average annual rainfall


# 1.3 Examine non-spatio-temporal data characteristics

# Mean and standard deviation for regional data
mean_rain <- mean(rain_matrix)
mean_scotland_n <- mean(rain_matrix[1,])
mean_scotland_e <- mean(rain_matrix[2,])
mean_scotland_w <- mean(rain_matrix[3,])
mean_england_e_ne <- mean(rain_matrix[4,])
mean_england_nw_n_wales <- mean(rain_matrix[5,])
mean_midlands <- mean(rain_matrix[6,])
mean_e_anglia <- mean(rain_matrix[7,])
mean_s_wales_england_sw <- mean(rain_matrix[8,])
mean_england_se_central_s <- mean(rain_matrix[9,])
mean_n_ireland <- mean(rain_matrix[10,])
sd_rain <- sd(rain_matrix)
mean_rain   # Mean monthly rainfall for all UK regions across all years
sd_rain

# Mean and standard deviation for weather station data
mean_station <- mean(station_matrix)
sd_station <- sd(station_matrix)
mean_station    # Mean monthly rainfall recorded at UK selected weather stations across all years
sd_station

# Histograms
hist(rain_matrix, col="lightblue", xlab="Monthly Rainfall in mm", main="Histogram of UK Regional Monthly Rainfall")
abline(v=mean_rain, col="red", lwd=3)
text(120, 3849, "Mean = 88.55")
# Shows positive skew, bounded by zero (cannot have -ve rainfall)

hist(station_matrix, col="lightgreen", xlab="Monthly Rainfall in mm", main="Histogram of UK Weather Station Monthly Rainfall")
abline(v=mean_station, col="red", lwd=3)
text(120, 7650, "Mean = 70.20")
# Shows positive skew, bounded by zero (cannot have -ve rainfall)

# Q-Q Plots
qqnorm(rain_matrix)
qqline(rain_matrix, col="red", lwd=3)
qqnorm(station_matrix)
qqline(station_matrix, col="red", lwd=3)


## 2 EXPLORATORY SPATIO-TEMPORAL DATA ANALYSIS AND VISUALISATION

# 2.1. Examining Temporal Characteristics

# First examine the average monthly rainfall, both on the regional dataset and the individual station dataset
# Regional dataset, monthly granularity, 186 years, average rainfall across all regions
plot(colMeans(rain_matrix), xlab="Year (monthly data)", ylab="Rainfall in mm", type="l", xaxt="n", main="Average UK Monthly Rainfall 1836-2021")
axis(1, at=seq(49, 2209, 240), labels=seq(1840, 2020, 20))
# No obvious trends - data too dense to really see - significant variation across the months within a year
# Suggest using annual averages?

# Station dataset, monthly granularity, 52 years, average rainfall across all weather stations
plot(colMeans(station_matrix), xlab="Year (monthly data)", ylab="Rainfall (mm)", type="l", xaxt="n", main="Average Monthly Rainfall for Selected UK Weather Stations 1965-2015")
axis(1, at=seq(5, 605, 120), labels=seq(1965, 2015, 10))
# No obvious trends - data too dense to really see - significant variation across the months within a year
# Hint of a seasonal pattern

# Try looking at a shorter time range: 2005-2015
plot(colMeans(station_matrix[,504:624]), xlab="Year (monthly data)", ylab="Rainfall (mm)", type="l", xaxt="n", main="Average Monthly Rainfall for Selected UK Weather Stations 2005-2015")
axis(1, at=seq(1, 121, 12), labels=seq(2005, 2015, 1))
# No obvious long-term trend - significant variation across the months within a year
# Looks like a potentially seasonal pattern

# Suggest using annual averages:
plot(colMeans(station_annual_matrix), xlab="Year", ylab="Rainfall (mm)", type="l", xaxt="n", main="Average Annual Rainfall for Selected UK Weather Stations 1965-2015")
axis(1, at=seq(1,51,5), labels=seq(1965, 2015, 5))
# Possible trend of increased rainfall with time, but extent of variation from year to year dominates

# Create Lattice Plot for 10 selected stations that cover all parts of the UK
# Rainfall to be dependent variable
# Month of the year (MMM.YY) the independent variable
station_melt <- melt(uk_station, id.vars=1:4, measure.vars = 5:ncol(uk_station))
colnames(station_melt)[5:6] <- c("MMMYY", "Rainfall")
station.chosen=c("Aberporth","Armagh","Bradford","Braemar","Cambridge","Oxford","Yeovilton","Lowestoft","Lerwick","Durham")
s <- station_melt[station %in% station.chosen,]
xyplot(Rainfall ~ MMMYY | Station, xlab = "MMM.YY", type ="l",
       layout = c(5,2),
       data = s,
       main = "Monthly Rainfall at Selected UK Weather Stations Sep 1965 to Aug 2016")
# No obvious trends here either, although possible increase over time seen in Lerwick?
# Data possibly too granular

# Use annual averages instead
station_annual_melt <- melt(uk_station_annual, id.vars=1:4, measure.vars = 5:ncol(uk_station_annual))
colnames(station_annual_melt)[5:6] <- c("Year", "Rainfall")
sam <- station_annual_melt[station %in% station.chosen,]
xyplot(Rainfall ~ Year | Station, xlab = "Year", type ="l",
       layout = c(5,2),
       data = sam,
       main = "Annual Average Rainfall at Selected UK Weather Stations 1965-2015")
# Data easier to interpret than monthly data, but still no consistent temporal trends
# Some hint at slightly increased rainfall over time, e.g.Lerwick and Braemar


# 2.3 Examining Spatial Characteristics

# Scatterplot Matrix for UK Weather Stations
pairs(~LONG+LAT+ALT+rowMeans(station_matrix),data=uk_station,main="Simple Scatterplot Matrix for Rainfall at UK Weather Stations")
# bottom left plot hints at a relationship between longitude and rainfall - further west, higher rainfall
# furthest right plot on second row also hints at a relationship between latitude and rainfall - further south, lower rainfall
# but trends are not clear

# 3D Scatterplots
scatterplot3d(x=uk_station$LONG, y=uk_station$LAT, z=rowMeans(station_matrix), main="3D Scatterplot of Rainfall at UK Weather Stations by Longitude and Latitude", xlab="Longitude", ylab="Latitude", zlab="Average Monthly Rainfall in mm")
scatter3D(uk_station$LONG, uk_station$LAT, rowMeans(station_matrix), main="3D Scatter of Rainfall at UK Weather Stations", xlab="Longitude", ylab="Latitude", zlab="Average Monthly Rainfall in mm")
plot3d(uk_station$LONG, uk_station$LAT, rowMeans(station_matrix), main="Dynamic 3D Plot of Rainfall at UK Weather Stations", xlab="Longitude", ylab="Latitude", zlab="Average Monthly Rainfall in mm")
# 3D plot confirms that there appears to be a relationship between latitude, longitude and rainfall, in that the further north and west, the higher the average monthly rainfall. But pattern is not clear

# Heatmap
heatmap(station_matrix,Rowv=NA,Colv=NA, col=cm.colors(256), scale="column", margins=c(5,3), xlab="MMM.YY", ylab="Station", main="Heatmap of Rainfall at UK Weather Stations, Unordered", cexCol=1.1, y.scale.components.subticks(n=10))
# Shows that there is not much temporal variation, as each row shows limited colour variation
# Try ordering stations by latitude to see if that shows anything
station_latorder <- uk_station[order(uk_station$LAT, decreasing=FALSE),]
station_latorder_matrix <- data.matrix(station_latorder[,5:ncol(uk_station)])
heatmap(station_latorder_matrix,Rowv=NA,Colv=NA, col=cm.colors(256), scale="column", margins=c(5,3), xlab="MMM.YY", ylab="Station (N at the top, S at the bottom)", main="Heatmap of Rainfall at UK Weather Stations, Ordered by Latitude", cexCol=1.1, y.scale.components.subticks(n=10))
# Definitely shows pinker colours (most rainfall) at the top (highest latitudes), i.e. furthest north
# Try ordering by longitude as well
station_longorder <- uk_station[order(uk_station$LONG, decreasing=TRUE),]
station_longorder_matrix <- data.matrix(station_longorder[,5:ncol(uk_station)])
heatmap(station_longorder_matrix,Rowv=NA,Colv=NA, col=cm.colors(256), scale="column", margins=c(5,3), xlab="MMM.YY", ylab="Station (W at the top, E at the bottom)", main="Heatmap of Rainfall at UK Weather Stations, Ordered by Longitude",cexCol=1.1, y.scale.components.subticks(n=10))
# Less obvious than latitude heat map, but suggests pinker colours (most rainfall) at the bottom (lowest longitudes), i.e. furthest west
# Station 12 is a bit of an anomaly (this is Lerwick, the furthest north by far in the Shetlands)
# Suggests latitude is more of a factor that longitude, but this might be expected as UK has a greater N/S extent than E/W

# Plot annual average rainfall for weather stations on a map
ann_ave <- cbind(uk_station_annual[1:4], rowMeans(station_annual_matrix))
colnames(ann_ave)[5] <- "Ave.Rainfall"
ann_ave[,2:3] <- projectMercator(ann_ave$LAT, ann_ave$LONG)
map <- openmap(c(49,-11), c(61,3), type='esri-topo')
autoplot.OpenStreetMap(map) +
  geom_point(data = ann_ave, aes(x = LONG, y = LAT, color = Ave.Rainfall, size = Ave.Rainfall)) +
  ggtitle("Annual Average Rainfall in UK, all years") +
  scale_color_gradient(low="lightblue", high="darkblue")
# Shows wettest areas are in the north and west; drier to south and east and away from coasts
# Eskdalemuir records most rain annually out of the 29 weather stations in the dataset


# Going back to the regional data
# Create breaks that can be used consistently in each map based on entire data set
brks=quantile(as.numeric(unlist(uk_rain@data[,-c(1)])), seq(0,1,1/5))
# Produce a choropleth map of UK rainfall by region for March 1990, which based on analysis in Excel is the month that shows the greatest contrast between wettest and driest region
tm_shape(uk_rain) +
  tm_fill("Mar_1990", style="fixed", palette="Blues", breaks=brks) +
  tm_borders("white") +
  tm_compass(position=c("left", "bottom")) +
  tm_legend(position=c("left", "top")) +
  tm_scale_bar()
# Very wet in Scotland, particularly in the N and W; dry elsewhere

# Produce a choropleth map of UK rainfall by region for Nov 1907, which based on analysis in Excel is the month that is closest to the overall average rainfall for each region
tm_shape(uk_rain) +
  tm_fill("Nov_1907", style="fixed", palette="Blues", breaks=brks) +
  tm_borders("white") +
  tm_compass(position=c("left", "bottom")) +
  tm_legend(position=c("left", "top")) +
  tm_scale_bar()
# Wettest in the West, driest in the East

# Whilst temperature follows fairly distinct seasonal patterns (warmer in summer, colder in winter), rainfall in the UK may not
# Produce choropleth maps for each month in 1984, which is identified as the year that is closest to the overall annual average rainfall for each region to see if there are any hints of seasonal patterns
year84 <- subset(uk_rain, select = c(Jan_1984:Dec_1984))
jan84 <- tm_shape(year84) +
  tm_fill("Jan_1984", style="fixed", palette="Blues", breaks=brks) +
  tm_borders("white") +
  tm_compass(position=c("right", "top"), size=1) +
  tm_legend(position=c("left", "top"), text.size=0.4) +
  tm_scale_bar()

feb84 <- tm_shape(year84) +
  tm_fill("Feb_1984", style="fixed", palette="Blues", breaks=brks) +
  tm_borders("white") +
  tm_legend(position=c("left", "top"), text.size=0.4)

mar84 <- tm_shape(year84) +
  tm_fill("Mar_1984", style="fixed", palette="Blues", breaks=brks) +
  tm_borders("white") +
  tm_legend(position=c("left", "top"), text.size=0.4)

apr84 <- tm_shape(year84) +
  tm_fill("Apr_1984", style="fixed", palette="Blues", breaks=brks) +
  tm_borders("white") +
  tm_legend(position=c("left", "top"), text.size=0.4)

may84 <- tm_shape(year84) +
  tm_fill("May_1984", style="fixed", palette="Blues", breaks=brks) +
  tm_borders("white") +
  tm_legend(position=c("left", "top"), text.size=0.4)

jun84 <- tm_shape(year84) +
  tm_fill("Jun_1984", style="fixed", palette="Blues", breaks=brks) +
  tm_borders("white") +
  tm_legend(position=c("left", "top"), text.size=0.4)

jul84 <- tm_shape(year84) +
  tm_fill("Jul_1984", style="fixed", palette="Blues", breaks=brks) +
  tm_borders("white") +
  tm_legend(position=c("left", "top"), text.size=0.4)

aug84 <- tm_shape(year84) +
  tm_fill("Aug_1984", style="fixed", palette="Blues", breaks=brks) +
  tm_borders("white") +
  tm_legend(position=c("left", "top"), text.size=0.4)

sep84 <- tm_shape(year84) +
  tm_fill("Sep_1984", style="fixed", palette="Blues", breaks=brks) +
  tm_borders("white") +
  tm_legend(position=c("left", "top"), text.size=0.4)

oct84 <- tm_shape(year84) +
  tm_fill("Oct_1984", style="fixed", palette="Blues", breaks=brks) +
  tm_borders("white") +
  tm_legend(position=c("left", "top"), text.size=0.4)

nov84 <- tm_shape(year84) +
  tm_fill("Nov_1984", style="fixed", palette="Blues", breaks=brks) +
  tm_borders("white") +
  tm_legend(position=c("left", "top"), text.size=0.4)

dec84 <- tm_shape(year84) +
  tm_fill("Dec_1984", style="fixed", palette="Blues", breaks=brks) +
  tm_borders("white") +
  tm_legend(position=c("left", "top"), text.size=0.4)

tmap_arrange(jan84, feb84, mar84, apr84, may84, jun84, jul84, aug84, sep84, oct84, nov84, dec84) # This takes a long time to plot

# These maps clearly show that (for this single year) the western and northern regions of the UK are the wettest, although N.Ireland is relatively drier than nearby W Scotland
# The driest regions are towards the east and the south
# This is in line with expectations based on knowledge of the south-westerly prevailing weather pattern in the UK, which brings rain in from the Atlantic and the central hills and mountains of England and Scotland create a rain shadow effect to the east
# Looking just at 1984, there is an indication of seasonality with more rain in the autumn/winter than spring/summer.
# Needs investigating further to see if this is a pattern that bears out across more/all years in the study.


# 2.3 Spatial and Temporal Autocorrelation

# 2.3.1. Spatial Autocorrelation

# Build a row standardised spatial weight matrix for the UK Regions (using Queen's case)
W <- nb2listw(poly2nb(uk_rain, queen=TRUE))
# Error due to empty neighbour sets found
# Test whether this is due to N.Ireland, which is not connected to rest of UK
W <- nb2listw(poly2nb(uk_rain[1:9,], queen=TRUE))
# Success - it was due to N.Ireland
W

# 2.3.1.1 Global Spatial Autocorrelation

# Calculate Moran's I

# This is a purely spatial index so can just take average rainfall for whole period under review
# Need to create subset of uk_rain that excludes N.Ireland due to empty neighbour sets issue above
gb_rain <- uk_rain[1:9,] # excluding N.Ireland
gb_rain_matrix <- data.matrix(gb_rain@data[,-c(1)])
gb_rain_ave <- rowMeans(gb_rain_matrix) # excluding N.Ireland
moran(gb_rain_ave, W, n=length(W$neighbours), S0=Szero(W))
# Gives a Moran's I value of 0.516

# Test min and max values to see whether reasonable on an absolute scale
moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(W)
# Moran range is approx -0.6 to 1.0, with mid point of approx 0.2.
# Therefore, Moran's I of 0.516 appears significant

# Test whether it is statistically significant using moran.test and moran.mc
moran.test(x=gb_rain_ave, listw=W)
moran.mc(x=gb_rain_ave, listw=W, nsim=9999)
# Both tests give p-values <0.05, so can conclude that there is significant spatial autocorrelation between rainfall levels across the regions of the UK at this spatial order

# 2.3.1.2 Local Spatial Autocorrelation

# Calculate Local Moran's I for Regional data
local_m <- localmoran(x=gb_rain_ave, listw=W)
gb_rain$Ii <- local_m[,"Ii"]
gb_rain$Iip_unadj <- local_m[,"Pr(z != E(Ii))"]
gb_rain$Unadjusted_significance <- "nonsignificant"
gb_rain$Unadjusted_significance[which(gb_rain$Iip_unadj < 0.05)] <- "significant"
tm_shape(gb_rain) + 
  tm_polygons(col="Unadjusted_significance", palette="-RdBu", style="quantile", lwd=0.1, border.alpha=0.1) +
  tm_layout(main.title="Local Moran's I by Region", main.title.size=1.25, title="P-values unadjusted", title.size=1,legend.position=c("left", "top")) +
  tm_scale_bar(width=0.25, position=c("right", "bottom")) + 
  tm_compass(position=c("right", "top"))
# "S Wales & England SW" and "England NW & N Wales" are the regions with the lowest local Moran values (close to zero)
# These regions are wet regions, being in the West, but they border the driest regions, i.e. those in the south and east
# "East Anglia" is a dry region bordering other dry regions, hence high local Moran
# "Scotland N" is a wet region bordering other wet regions, hence high local Moran
# "East Anglia" and "Scotland W" have significant local Moran's I based on unadjusted p-values

local_m_adj <- localmoran(x=gb_rain_ave, listw=W, p.adjust.method="bonferroni")
gb_rain$Iip_adj <- local_m_adj[,"Pr(z != E(Ii))"]
gb_rain$Adjusted_significance <- "nonsignificant"
gb_rain$Adjusted_significance[which(gb_rain$Iip_adj <0.05)] <- "significant"
tm_shape(gb_rain) + 
  tm_polygons(col="Adjusted_significance", palette="-RdBu", style="quantile", lwd=0.1, border.alpha=0.1) +
  tm_layout(main.title="Local Moran's I by Region", main.title.size=1.25, title="P-values adjusted", title.size=1, legend.position=c("left", "top")) +
  tm_scale_bar(width=0.25, position=c("right", "bottom")) + 
  tm_compass(position=c("right", "top"))
# No regions with significant local Moran's I based on adjusted p-values

# Measure autocorrelation in weather station point data using semivariogram
# Use average rainfall across all time periods under review as temporal variation not taken into account in semivariogram
coords = list(projectMercator(uk_station[,3], uk_station[,2]))
plot(variogram(list(rowMeans(station_matrix)), location=coords), cex=2, main="Semivariogram for UK Weather Station Rainfall", col="blue")
# Result is a bit scattered, but can conclude that rainfall levels at weather stations that are closer are more similar than those further away

# Consider different directions
plot(variogram(list(rowMeans(station_matrix)), location=coords, alpha=c(0,45,90,135)), cex=1.5, main="Directional semivariograms for UK Weather Station Rainfall", col="blue")
# Not very clear results, but hints of anisotropy in that not all semivariograms are the same


# 2.3.2 Temporal Autocorrelation

# Using the monthly rainfall data by region
# Create dataframes containing lagged variables at lag of 1 month
# Just describing the time intervals as 'time_in_months' and giving them a numerical index as cannot find a way to convert the format in the original data to a date that the formula will accept
Scotland_N_lagged <- data.frame(time_in_months = 1:2231, t=rain_matrix[1,][2:(ncol(rain_matrix))], t_minus_1=rain_matrix[1,][1:(ncol(rain_matrix)-1)])
Scotland_E_lagged <- data.frame(time_in_months = 1:2231, t=rain_matrix[2,][2:(ncol(rain_matrix))], t_minus_1=rain_matrix[2,][1:(ncol(rain_matrix)-1)])
Scotland_W_lagged <- data.frame(time_in_months = 1:2231, t=rain_matrix[3,][2:(ncol(rain_matrix))], t_minus_1=rain_matrix[3,][1:(ncol(rain_matrix)-1)])
England_E_NE_lagged <- data.frame(time_in_months = 1:2231, t=rain_matrix[4,][2:(ncol(rain_matrix))], t_minus_1=rain_matrix[4,][1:(ncol(rain_matrix)-1)])
England_NW_N_Wales_lagged <- data.frame(time_in_months = 1:2231, t=rain_matrix[5,][2:(ncol(rain_matrix))], t_minus_1=rain_matrix[5,][1:(ncol(rain_matrix)-1)])
Midlands_lagged <- data.frame(time_in_months = 1:2231, t=rain_matrix[6,][2:(ncol(rain_matrix))], t_minus_1=rain_matrix[6,][1:(ncol(rain_matrix)-1)])
East_Anglia_lagged <- data.frame(time_in_months = 1:2231, t=rain_matrix[7,][2:(ncol(rain_matrix))], t_minus_1=rain_matrix[7,][1:(ncol(rain_matrix)-1)])
S_Wales_England_SW_lagged <- data.frame(time_in_months = 1:2231, t=rain_matrix[8,][2:(ncol(rain_matrix))], t_minus_1=rain_matrix[8,][1:(ncol(rain_matrix)-1)])
England_SE_lagged <- data.frame(time_in_months = 1:2231, t=rain_matrix[9,][2:(ncol(rain_matrix))], t_minus_1=rain_matrix[9,][1:(ncol(rain_matrix)-1)])
N_Ireland_lagged <- data.frame(time_in_months = 1:2231, t=rain_matrix[10,][2:(ncol(rain_matrix))], t_minus_1=rain_matrix[10,][1:(ncol(rain_matrix)-1)])

# Calculate PMCC based on 1 month lags
SN_r <- round(cor(Scotland_N_lagged$t, Scotland_N_lagged$t_minus_1), 3) # gives 0.306
SE_r <- round(cor(Scotland_E_lagged$t, Scotland_E_lagged$t_minus_1), 3) # gives 0.152
SW_r <- round(cor(Scotland_W_lagged$t, Scotland_W_lagged$t_minus_1), 3) # gives 0.280
EENE_r <- round(cor(England_E_NE_lagged$t, England_E_NE_lagged$t_minus_1), 3) # gives 0.091
EWNW_r <- round(cor(England_NW_N_Wales_lagged$t, England_NW_N_Wales_lagged$t_minus_1), 3) # gives 0.180
M_r <- round(cor(Midlands_lagged$t, Midlands_lagged$t_minus_1), 3)  # gives 0.081
EA_r <- round(cor(East_Anglia_lagged$t, East_Anglia_lagged$t_minus_1), 3) # gives 0.103
EWSW_r <- round(cor(S_Wales_England_SW_lagged$t, S_Wales_England_SW_lagged$t_minus_1), 3) # gives 0.1178
ESE_r <- round(cor(England_SE_lagged$t, England_SE_lagged$t_minus_1), 3)  # gives 0.127
NI_r <- round(cor(N_Ireland_lagged$t, N_Ireland_lagged$t_minus_1), 3) # gives 0.135

# Most of these PMCCs show only weak or very weak correlation at 1 month lag interval
# Scotland N shows the highest correlation at 0.306
# Plot to visualise
p1 <- ggplot(Scotland_N_lagged, aes(x=time_in_months, y=t)) + geom_line()
p2 <- ggplot(Scotland_N_lagged, aes(x=t, y=t_minus_1)) +
  geom_point() +
  labs(y="t-1") +
  geom_smooth(method="lm") +
  annotate("text", 30, 320, label=paste("r =", SN_r))

grid.arrange(p1,p2, nrow=1)

# Now just plot Scotland N PMCC
ggplot(Scotland_N_lagged, aes(x=t, y=t_minus_1)) +
  geom_point() +
  labs(y="t-1") +
  geom_smooth(method="lm") +
  annotate("text", 30, 320, label=paste("r =", SN_r))

# Also try looking at 12month lag
Scotland_N_lagged12 <- data.frame(time_in_months = 12:2231, t=rain_matrix[1,][13:(ncol(rain_matrix))], t_minus_12=rain_matrix[1,][1:(ncol(rain_matrix)-12)])
SN_r_minus12 <- round(cor(Scotland_N_lagged12$t, Scotland_N_lagged12$t_minus_12), 3)
ggplot(Scotland_N_lagged12, aes(x=t, y=t_minus_12)) +
  geom_point() +
  labs(y="t-12") +
  geom_smooth(method="lm") +
  annotate("text", 30, 320, label=paste("r =", SN_r_minus12))
# PMCC of 0.331 implies stronger temporal correlation at 12 month lag than 1 month lag
# Suggests rainfall in a given month is more similar to the same month in the previous year than in the previous month.
# Suggests rainfall varies more month-to-month than year-to-year in Scotland N region

# Next, examining annual averages by weather station
# First create a dataset of annual averages across all weather stations
uk_station_annual_aves <- colMeans(station_annual_matrix)
uk_annual_lagged <- data.frame(year = 1965:2014, t=uk_station_annual_aves[2:(length(uk_station_annual_aves))], t_minus_1=uk_station_annual_aves[1:(length(uk_station_annual_aves)-1)])
uk_r <- round(cor(uk_annual_lagged$t, uk_annual_lagged$t_minus_1), 3)

# Plot to visualise
p3 <- ggplot(uk_annual_lagged, aes(x=year, y=t)) + geom_line()
p4 <- ggplot(uk_annual_lagged, aes(x=t, y=t_minus_1)) +
  geom_point() +
  labs(y="t-1") +
  geom_smooth(method="lm") +
  annotate("text", 60, 90, label=paste("r =", uk_r))

grid.arrange(p3,p4, nrow=1)
# Shows that annual average rainfall is weakly correlated at a lag of one year (PMCC=0.121)
# This is not surprising as no reason to suggest that rainfall in one year is related to rainfall in previous year


# 2.3.3 Spatio-Temporal Autocorrelation

# Import starima_package
source("starima_package.R")

# Attempt space-time semivariogram using UK weather station point data
# Project points to get spatial lag in metres
points <- SpatialPoints(uk_station_annual[,2:3], proj4string = CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
# Convert years to correct date format
years <- seq(as.Date("1965-01-01"), length=51, by="year")
# Create spatio-temporal data frame
stfdf <- STFDF(points, years, data.frame(as.vector(t(station_annual_matrix))))
names(stfdf@data) <- "Rainfall"

# Calculate space-time semivariogram and plot
station_ST_var <- variogram(Rainfall~1, stfdf, width=50, cutoff=500, tlags=0:10)
plot(station_ST_var)
plot(station_ST_var, wireframe=T)
# What does this show?
# Semivariogram increases much more rapidly with increasing temporal lag than with increasing spatial lag
# But time and space are not the same thing so cannot directly compare these 'distances'


## 3 STATISTICAL MODELLING OF TIME SERIES AND SPATIO-TEMPORAL SERIES

# 3.1 Time Series Decomposition

# Use STL to decompose the Scotland N data
# Choosing this dataset as it showed the highest levels of temporal autocorrelation
# First, separate the Scotland N data from the rest of the UK rain data and convert to a time-series
Scot_N_ts <- ts(rain_matrix[1,], frequency = 12, start=c(1836,1))

# Then run stl function (t.window parameter should be odd or NULL)
decom <- stl(Scot_N_ts, t.window=NULL, s.window="periodic")
autoplot(decom)
# Hard to see due to length of timeseries

# Try on shorter timeseries (1990-2021)
Scot_N_sts <- ts(rain_matrix[1,1849:2232], frequency = 12, start=c(1990,1))
decom <- stl(Scot_N_sts, t.window=25, s.window="periodic")
autoplot(decom)
# Various different values of t.window used to extract a trend. Settled for t.window=25 as remainders smaller
# Trend component not clear - may not need nonseasonal differencing; may only need to use seasonal differencing
# Clear seasonality extracted through the stl decomposition process, although seasonal component only ranges +/-60mm
# Remainders still significant in size, typically ranging +/-100mm - larger than the seasonal component (which is actually smaller than the trend component based on the bars to the side of the plots)


# 3.2 ARIMA

# Use the Box-Jenkins approach to ARIMA modelling

# 3.2.1 Exploratory Data Analysis

# Have already performed this above
# UK Regional rainfall data shows seasonality as it is presented as monthly, with no clear trend up or down over time
# Annual average rainfall for UK Weather Stations does not show seasonality as it is annual, but no clear trend up or down over time: it appears stationary

# Lag plot for Scotland N region
lag.plot(rain_matrix[1,], lags=12, do.lines=FALSE)
# No clear patterns identified, which is a bit surprising considering seasonality already identified

# 3.2.2 Differencing, AutoCorrelation and Partial Autocorrelation Factors

# Examine the Scotland N time series and associated ACF
autoplot(Scot_N_ts) + ggtitle("Scotland N timeseries 1836 to 2021")
acf(rain_matrix[1,], lag.max=36, xlab="Lag", ylab="ACF", main="ACF Plot of Scotland N timeseries")
# Using rain_matrix[1,] in the acf instead of the Scot_N_ts timeseries because timeseries frequency is 12 (months) and lags indices on x-axis therefore show for years, not months
# ACF shows strong seasonal pattern with positive peaks at lags 0 (as per default), 12, 24, 36, 48 and negative peaks at lags 6, 18, 30, 42
# Suggests seasonal differencing is required with an order of 12, which makes sense as this is monthly data

# Seasonal differencing
# Perform first seasonal differencing at lag 12 and re-run the ACF
SN_s_diff <- diff(rain_matrix[1,], lag=12, differences=1)
acf(SN_s_diff, lag.max=36, xlab="Lag", ylab="ACF", main="ACF Plot of Seasonally Differenced Scotland N timeseries")
# Seasonal differencing has pretty much removed the seasonal pattern after lag 12
# Fits description of "one or more spikes, rest essentially zero", which implies MA model
# Significant autocorrelation remains at lags 1 and 12 (seasonal lag 1)

# Now do PACF for seasonally differenced data
pacf(SN_s_diff, lag.max=36, xlab="Lag", ylab="PACF", main="PACF Plot of Seasonally Differenced Scotland N timeseries")
# Significant PACF values seen at lags 1, 12 (seasonal lag 1), 24 (seasonal lag 2), 25, 36 (seasonal lag 3), 48 (seasonal lag 4)

# Try non-seasonal differencing and re-run ACF and PACF
SN_ns_diff <- diff(rain_matrix[1,])
acf(SN_ns_diff, lag.max=36, xlab="Lag", ylab="ACF", main="ACF Plot of Non-Seasonally Differenced Scotland N timeseries")
pacf(SN_ns_diff, lag.max=36, xlab="Lag", ylab="PACF", main="PACF Plot of Non-Seasonally Differenced Scotland N timeseries")
# ACF shows significant -ve ACF at lag 1 and +ve ACF at lags 12, 24, 36 which suggests seasonal pattern remains
# PACF shows significant -ve PACF at lag 1, which slowly decays to insignificant values

# Try seasonal and nonseasonal differencing together
SN_sns_diff <- diff(diff(rain_matrix[1,], lag=12), lag=1)
acf(SN_sns_diff, lag.max=36, xlab="Lag", ylab="ACF", main="ACF Plot of Non-Seasonally and Seasonally Differenced Scotland N timeseries")
pacf(SN_sns_diff, lag.max=36, xlab="Lag", ylab="PACF", main="PACF Plot of NonSeasonally and Seasonally Differenced Monthly Rainfall for Scotland N Region")




# Examine the ACF and PACF plots together
p5 <- autoplot(acf(SN_s_diff, lag.max=36, plot=FALSE)) + ggtitle("ACF, Scotland N Seasonally Differenced")
p6 <- autoplot(pacf(SN_s_diff, lag.max=36, plot=FALSE)) + ggtitle("PACF, Scotland N Seasonally Differenced")
p7 <- autoplot(acf(SN_ns_diff, lag.max=36, plot=FALSE)) + ggtitle("ACF, Scotland N NonSeasonally Differenced")
p8 <- autoplot(pacf(SN_ns_diff, lag.max=36, plot=FALSE)) + ggtitle("PACF, Scotland N NonSeasonally Differenced")
p9 <- autoplot(acf(SN_sns_diff, lag.max=36, plot=FALSE)) + ggtitle("ACF, Scotland N Seasonally and NonSeasonally Differenced")
p10 <- autoplot(pacf(SN_sns_diff, lag.max=36, plot=FALSE)) + ggtitle("PACF, Scotland N Seasonally and NonSeasonally Differenced")
grid.arrange(p5,p6,p7,p8,p9,p10)


pacf(uk_station_annual_aves, lag.max=36, main="PACF, Annual Average Rainfall for Selected UK Weather Stations")
# No statistically significant PACF at any lags. Implies the data is random from one year to the next


# 3.2.2.1 Spatio-Temporal ACF and PACF

# Calculate spatio-temporal autocorrelation factors (stacf) using UK regional (areal) data (excluding N.Ireland)
weight_matrix <- listw2mat(W)
stacf(t(gb_rain_matrix), weight_matrix, 48)
# Shows seasonal pattern with peaks at lag 1, 13, 25, 37
# and troughs at lag 7, 19, 30, 43

# Calculate spatio-temporal partial autocorrelation factors (stpacf) using UK regional (areal) data (excluding N.Ireland)
stpacf(t(gb_rain_matrix), weight_matrix, 12)
# stpacf is insignificant for all temporal lags - implies the data are essentially random


# 3.2.3. Parameter Estimation and Fitting

# ACF for Scotland N showed seasonality with a seasonal component of order 12, being the 12 months of the year
acf(rain_matrix[1,], lag.max=50, main="AutoCorrelation Plot of Monthly Rainfall for Scotland N Region")

# Perform seasonal differencing and re-run the ACF
SN_diff <- diff(rain_matrix[1,], lag=12, differences=1)
acf(SN_diff, lag.max=50, xlab="Lag", ylab="ACF", main="Seasonally Differenced ACF Plot of Monthly Rainfall for Scotland N Region")
# Seasonal differencing has pretty much removed the seasonal pattern after lag 12
# Fits description of "one or more spikes, rest essentially zero"
# Significant autocorrelation remains at lags 1 and 12
# MA order identified by where plot becomes zero. Only have 1 significant lag at lag=1, therefore q = 1

# This implies a model of ARIMA(0,0,1)(0,1,1)12
# q=1 for non-seasonal MA due to only 1 significant lag in seasonally differenced ACF
# Q=1 for seasonal MA as only have 1 significant lag at lag=12 and nothing at e.g. lag=24 or 36
# D=1 as only used 1 difference in seasonal differencing
# m=12 as monthly data

# Now do PACF for seasonally differenced data
pacf(SN_diff, lag.max=50, xlab="Lag", ylab="PACF", main="Seasonally Differenced PACF Plot of Monthly Rainfall for Scotland N Region")
# Significant PACF values seen at lags 1, 12, 24, 25, 36, 48
# Suggests 2 or 3 or more seasonal AR terms? Try 2 and see what happens
# This would give ARIMA(1,0,1)(2,1,1)12# Train over first 176 years, test over last 10 years
fit.ar <- arima(rain_matrix[1,1:2112], order=c(0,1,2), seasonal=list(order=c(0,1,1), period=12))
fit.ar
# Gives sigma^2 of 2405, log likelihood of -11169.6, aic of 22351.21
# Running other models does not give better results
# In fact this seems to be the only model with a seasonal component that runs without error
# Running ARIMA(1,0,1)(3,1,1)12 gives an error
# Running ARIMA(1,0,1)(1,1,1)12 gives NaNs
# Running ARIMA(2,0,1)(2,1,1)12 gives NaNs

NRMSE_fit <- NRMSE(res=fit.ar$residuals, obs=rain_matrix[1,1:2112])
NRMSE_fit
# gives figure of 0.8485827

tsdiag(fit.ar)

pre.ar <- predict(fit.ar, n.ahead=12)
matplot(1:12, cbind(rain_matrix[1,2113:2124], pre.ar$pred), type="l", main="", xlab="Month", ylab="Rainfall in mm")
# Black line is actual, red-dashed line is prediction

fit.Ar <- Arima(rain_matrix[1,1:2112], order=c(1,0,1), seasonal=list(order=c(2,1,1), period=12))
fit.Ar

pre.Ar <- Arima(rain_matrix[1,2113:ncol(rain_matrix)], model=fit.Ar)
matplot(cbind(pre.Ar$fitted, pre.Ar$x), type="l")


# Have a look at auto.ar to see what it gives as a suggested ARIMA model and compare with my model
# Adapted from https://medium.com/analytics-vidhya/sarima-forecasting-seasonal-data-with-python-and-r-2e7472dfad83
fit.auto.ar <- auto.arima(rain_matrix[1,1:2112], trace=TRUE, test="kpss", ic="bic")
fit.auto.ar
# Auto.ar gives "ARIMA(1,1,0) with drift"
# Interestingly this didn't try any seasonal ARIMA terms

# Check what forecast versus actuals looks like
pre.auto.ar <- Arima(rain_matrix[1,2113:ncol(rain_matrix)], model=fit.auto.ar)
matplot(cbind(pre.auto.ar$fitted, pre.auto.ar$x), type="l")
# Black line is prediction, red-dashed line is actuals
# Gives a close fit in terms of shape, but prediction seems to be 1 or 2 steps behind actual??



## 5 ARTIFICIAL NEURAL NETWORKS

# Use rain_matrix created earlier
X <- t(as.matrix(rain_matrix))
Y <- as.matrix(X[-1,])

# Training on first 80% of data (0.8*2232 = 1785.6 = 1786)
rain.nnet <- nnet(X[1:1786, 1:10], Y[1:1786, 1:10], decay=1e-6, linout=TRUE, size=2)
rain.nnet["fitted.values"]
# All values seem to be constant or one of two fixed values within a narrow range - no seasonality

# Run the prediction on the remaining 20% of data (445 months)
# Y data set finishes at 2231 as there is one less month of data in Y compared to X
rain.pred <-predict(rain.nnet, Y[1787:2231, 1:10])


# Look at results for the first region: Scotland N
rain.pred

matplot(cbind(Y[1787:2231,1], rain.pred[,1]), ylab="Monthly average rainfall", xlab="Time (in months)", main="Scotland N", type="l")
# Looks bad
# Possible explanations - wrong parameters - have tried many different combinations of decay and size, to no avail - tried mygrid approach from tutorial but suspect that only works with regression rather than time-series
# Try expand.grid from (regression) tutorial to determine optimal parameter settings


mygrid <- expand.grid(.decay=c(0.5, 0.1, 0.01, 1e-3, 1e-4, 1e-5, 1e-6), .size=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

rain.nnet.parameters <- train(Y, data=X, method="nnet", metric="Rsquared", maxit=100, tuneGrid=mygrid, trace=F)

?rnn

train <- trainr(Y[1:1786, 1], X[1:1786, 1], learningrate = 0.1)
