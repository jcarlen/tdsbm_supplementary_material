# Code to reproduce the examples in our TDSBM paper 

# Jane Carlen
# Created: 4-5-18

#Set to the path to your tdsbm_supplementary_material folder:
setwd("~/Documents/285J/tdsbm_supplementary_material")

################################################################################################
# 0. Setup ####

install_list = c("devtools", "ggplot2", "dplyr", "tidyr", "ggmap","reshape2",
                 "lubridate","rgdal","raster","broom", "RColorBrewer","scales",
                 "cowplot","stringr", "foreign","xtable", "latex2exp", "mapproj",
                 "sbmt")

# install packages as necessary
for (pkg in install_list) {
  if(!pkg %in% names(installed.packages()[,1])) {
    if (pkg == "sbmt") {#install from github
      devtools::install_github("jcarlen/sbm", subdir = "sbmt") 
    } else {install.packages(pkg)} #install from CRAN
  }
}

package_list = c("ggplot2", "dplyr", "tidyr", "ggmap","reshape2",
                 "lubridate","rgdal", "broom", "RColorBrewer","scales",
                 "cowplot","stringr", "foreign", "latex2exp",
                 "sbmt") #don't load devtools, raster or xtable

# load all packages
lapply(package_list, require, character.only = TRUE)

# ggplot2 theme
theme_set(theme_minimal())

################################################################################################
# 1. Clean Data ####
################################################################################################
#   1a. LA - October - December 2016  ----

#Unzip if not already done:
if (!file.exists("data/LA/metro_station_table.csv")) {
  unzip("data/LA/metro_station_table.zip", exdir = "data/LA")
}
if (!file.exists("data/LA/Metro_trips_Q4_2016.csv")) {
  unzip("data/LA/Metro_trips_Q4_2016.zip", exdir = "data/LA")
}

la_station_info = read.csv("data/LA/metro_station_table.csv")
la = read.csv("data/LA/Metro_trips_Q4_2016.csv"); nrow(la) #43202 entries
la$count = 1
# Clean 

#   Clean ####
##    change duration to minutes ----
la$duration = round(la$duration/60,2) 

##    remove outlier trips ----
# cost is $3.50 per each 30min
la = subset(la, duration > 2); nrow(la) # now 41771 entries
la = subset(la, duration < 90); nrow(la) # now 40327 entries

##    remove test station ----
la = subset(la, !end_station_id %in% c(3000, 4108,"\\N"), drop = T) #remove virtual station and warehouse trips
la = subset(la, !start_station_id %in% c(4108, "\\N"), drop = T) 
la = droplevels(la)
la_station_info = la_station_info %>% filter(Status == " Active  " & Station.ID !=3000)
la_station_info$Station.ID = as.character(la_station_info$Station.ID)
nrow(la) #40130

##    only keep stations that have at least one start and end of a trip (doesn't affect la) ----
# In LA, all stations retained to this point have at least 40 starts and ends
identical(sort(unique(la$start_station_id)),  sort(unique(la$end_station_id)))
head(sort(table(la$start_station_id)))
head(sort(table(la$end_station_id)))

##    DON'T remove self-trips ----
#la = filter(la, start_station_id != end_station_id)
nrow(filter(la, start_station_id == end_station_id)) #3084 self-edges

##    add time fields ----
la$hour = sapply(la$start_time, function(x) {unclass(strptime(x, format = "%m/%d/%Y %H:%M"))$hour})
la$dayofweek = sapply(la$start_time, function(x) {unclass(strptime(x, format = "%m/%d/%Y %H:%M"))$wday})
la$weekend = (la$dayofweek == 0 | la$dayofweek == 6)
mean(la$weekend) # 26.63% weekend trips
la$count = 1
nrow(la) #now 40130 entries

#   Save (for reproducibility) ---- ----

#final size
dim(la)
n_distinct(la$start_station_id); n_distinct(la$end_station_id)

# cleaning removed 3072 entries, from 43202 to to 40130 (7% of entries)
write.csv(la[,c("start_time", "end_time", "start_station_id", "end_station_id",
                "start_lat", "start_lon", "end_lat", "end_lon")],
          "data/cleaned/LA16_cleaned_final.csv", row.names = F)

#--------------------------------------------------------------------------------------------
#  1b. SF - September 2015 - Aug 2016 ----

#Unzip if not already done:
if (!file.exists("data/SF/201608_station_data.csv")) {
  unzip("data/SF/201608_station_data.zip", exdir = "data/SF")
}
if (!file.exists("data/SF/201608_trip_data.csv")) {
  unzip("data/SF/201608_trip_data.zip", exdir = "data/SF")
}

bay = read.csv("data/SF/201608_trip_data.csv")
bay_stations = read.csv("data/SF/201608_station_data.csv")
bay_stations = bay_stations[1:67,]

#  Clean ####

##    reformat times ----
bay$Start.Date = as.POSIXct(bay$Start.Date, format = "%m/%d/%Y %H:%M")
bay$End.Date = as.POSIXct(bay$End.Date, format = "%m/%d/%Y %H:%M")
range(bay$Start.Date); range(bay$End.Date)
#dim(subset(bay, bay$Start.Date > "2016-01-01"))

##    remove non-sf stations ----
sf_station = subset(bay_stations, bay_stations$landmark == "San Francisco")
sf_station_name = subset(bay_stations, bay_stations$landmark == "San Francisco")$name
sf = filter(bay, bay$End.Station%in%sf_station_name); nrow(sf)
sf = filter(sf, Start.Station != "Mountain View City Hall") # <- removes just one trip
nrow(sf); # 280078 entries

##    change duration to minutes -----
sf$Duration = sf$Duration/60

##    remove outlier trips -----
# cost is $3 per first 30 min, per 15 min after that.
sf = filter(sf, sf$Duration < 90); nrow(sf) #now 276116 entreis
sf = filter(sf, sf$Duration > 2); nrow(sf) #now 274004 entries
sf = droplevels(sf)

##    only keep stations that have at least one start and end of a trip (removes 2 stations) ----

identical(sort(unique(sf$Start.Station)),  sort(unique(sf$End.Station))) #FALSE
unique(sf$End.Station)[!unique(sf$End.Station) %in% unique(sf$Start.Station)] #0
unique(sf$Start.Station)[!unique(sf$Start.Station) %in% unique(sf$End.Station)] #2
sf = filter(sf, sf$Start.Station %in% unique(sf$End.Station))
dim(sf) #267433     11
n_distinct(sf$Start.Station) # now 35 stations (2 were removed)
identical(sort(unique(sf$Start.Station)),  sort(unique(sf$End.Station))) #now TRUE

sf = left_join(sf, sf_station, by = c("Start.Station" = "name"))
sf = left_join(sf, sf_station, by = c("End.Station" = "name"))
names(sf)[c(12,13,14,18,19,20)] = c("start_station_id", "start_lat",
                                    "start_lon","end_station_id",
                                    "end_lat","end_lon")
sf$count = 1

##    DON'T remove self-trips -----
# sf = filter(sf, Start.Station != End.Station)
nrow(filter(sf, Start.Station == End.Station)) # 3787 self-edges

##    add time fields ----
sf$hour = sapply(sf$Start.Date, function(x) {unclass(strptime(x, format = "%Y-%m-%d %H:%M:%S"))$hour})
sf$dayofweek = sapply(sf$Start.Date, function(x) {unclass(strptime(x, format = "%Y-%m-%d %H:%M:%S"))$wday})
sf$weekend = (sf$dayofweek == 0 | sf$dayofweek == 6)
sf = filter(sf, !is.na(hour)) #removed 21 trips with no starting time. they're short and end just after midnight -- maybe started exactly at midnight? Or they're system testing/reboot?
mean(sf$weekend) # 7.95% weekend tripss
nrow(sf) # now 267412 entries

#  Save (for reproducibility) ----

#final size - 267142 trips, 35 stations
dim(sf)
n_distinct(sf$start_station_id); n_distinct(sf$end_station_id)

# cleaning remvoed 12666 entries, from 280078 to to 267412 (5% of entries)
write.csv(sf[,c("Start.Date", "End.Date", "start_station_id", "end_station_id",
                "start_lat", "start_lon", "end_lat", "end_lon")],
          "data/cleaned/SF16_cleaned_final.csv", row.names = F)

#--------------------------------------------------------------------------------------------
#  1c. NY - October 2016 ----

if (!file.exists("data/NY/201610-citibike-tripdata.csv")) {
  unzip("data/NY/201610-citibike-tripdata.zip", exdir = "data/NY")
}

ny1610 = read.csv("data/NY/201610-citibike-tripdata.csv"); nrow(ny1610) #1573872 entries
ny1610$count = rep(1, nrow(ny1610))

#  Clean----

##    change duration to minutes ----
names(ny1610)[1] = "tripduration"
ny1610$tripduration = round(ny1610$tripduration/60, 2)

##    remove outlier trips ----
# for subscribers it was 45 min free, $2.50 per additional 15 min 
ny1610 = filter(ny1610, tripduration > 2); nrow(ny1610) # now 1556306 entries
ny1610 = filter(ny1610, tripduration < 120); nrow(ny1610) # now 1552215 entries

##    remove test and maintenance stations ----
nystationnames = unique(c(as.character(ny1610$Start.Station.Name), as.character(ny1610$End.Station.Name)))
nydepots = nystationnames[grepl(nystationnames, pattern = "Depot|Kiosk|SSP")]
ny1610 = filter(ny1610, !((Start.Station.Name %in% nydepots) | (End.Station.Name %in% nydepots)))
ny1610 =  droplevels(ny1610)
nrow(ny1610) # now 1551709 entries

##    only keep stations that have at least one start and end of a trip (removes 6 stations) ----
length(table(ny1610$Start.Station.Name)); length(table(ny1610$End.Station.Name))
unique(ny1610$Start.Station.Name)[!unique(ny1610$Start.Station.Name) %in% unique(ny1610$End.Station.Name)] #0
unique(ny1610$End.Station.Name)[!unique(ny1610$End.Station.Name) %in% unique(ny1610$Start.Station.Name)] #6

removestations = as.character(unique(ny1610$End.Station.Name)[!unique(ny1610$End.Station.Name) %in% unique(ny1610$Start.Station.Name)])
removetrips = which(ny1610$End.Station.Name%in%removestations | ny1610$Start.Station.Name%in%removestations)
ny1610 = ny1610[-removetrips,]
ny1610 =  droplevels(ny1610)
nrow(ny1610) # now 1551701 entries

identical(sort(unique(ny1610$Start.Station.ID)), sort(unique(ny1610$End.Station.ID)))
n_distinct(ny1610$Start.Station.ID) #now 602 stations
head(sort(table(ny1610$End.Station.Name)))
head(sort(table(ny1610$Start.Station.Name)))

# if using with KLOptimization.cpp change station ids to 1-N scale,

##    also remove Soissons Landing on governor' island (only 9 trips to or from it, mostly self-trips) ----
filter(ny1610, End.Station.Name =="Soissons Landing" | Start.Station.Name =="Soissons Landing")
ny1610 = filter(ny1610, !(End.Station.Name =="Soissons Landing" | Start.Station.Name =="Soissons Landing"))
nrow(ny1610) # now 1551692 entries

##    DON'T remove self-trips ----
# ny1610 = filter(ny1610, Start.Station.ID != End.Station.ID)
nrow(filter(ny1610, Start.Station.ID == End.Station.ID)) # 22894 self-edges

##    add time fields ----
ny1610$hour = sapply(ny1610$Start.Time, function(x) {unclass(strptime(x, format = "%Y-%m-%d %H:%M:%S"))$hour})
ny1610$dayofweek = sapply(ny1610$Start.Time, function(x) {unclass(strptime(x, format = "%Y-%m-%d %H:%M:%S"))$wday})
ny1610$weekend = (ny1610$dayofweek == 0 | ny1610$dayofweek == 6)
mean(ny1610$weekend) # 24.4% weekend tripss
#ny1610 = filter(ny1610, weekend == FALSE)
nrow(ny1610) # now 1551692 entries

#  Save (for reproducibility)----

#final size - 1551692 trips, 601 stations
dim(ny1610)
n_distinct(ny1610$Start.Station.ID); n_distinct(ny1610$End.Station.ID)

# cleaning removed 22180 entries, from 1573872 to 1551692 (1.4% of entries)

write.csv(ny1610[,c("Start.Time", "Stop.Time", "Start.Station.ID", "End.Station.ID",
                    "Start.Station.Latitude", "Start.Station.Longitude", 
                    "End.Station.Latitude", "End.Station.Longitude")],
          "data/cleaned/Oct16_nyc_cleaned_final.csv", row.names = F)


#--------------------------------------------------------------------------------------------
#  1d. Aggregate station info (all) ####
##    LA  ----
la.station.all = data.frame(id = aggregate(count ~ start_station_id, data = la, sum)$start_station_id,
                        in.degree = aggregate(count ~ end_station_id, data = la, sum)$count,
                        out.degree = aggregate(count ~ start_station_id, data = la, sum)$count) %>%
  mutate(degree = in.degree + out.degree) %>%  
  left_join(la[, c("start_lat", "start_lon", "start_station_id")] %>% 
              filter(!duplicated(start_station_id)), by = c("id" = "start_station_id"))
names(la.station.all)[match(c("start_lat", "start_lon"), names(la.station.all))] = c("lat","lon")

##    SF  ----
sf.station.all = data.frame(id = aggregate(count ~ start_station_id, data = sf, sum)$start_station_id,
                        in.degree = aggregate(count ~ end_station_id, data = sf, sum)$count,
                        out.degree = aggregate(count ~ start_station_id, data = sf, sum)$count) %>%
  mutate(degree = in.degree + out.degree) %>%  
  left_join(sf[, c("start_lat", "start_lon", "start_station_id")] %>% 
              filter(!duplicated(start_station_id)), by = c("id" = "start_station_id"))

names(sf.station.all)[match(c("start_lat", "start_lon"), names(sf.station.all))] = c("lat","lon")


##    NY  ----
ny1610.station.all = data.frame(id = aggregate(count ~ Start.Station.ID, data = ny1610, sum)$Start.Station.ID,
                            in.degree = aggregate(count ~ End.Station.ID, data = ny1610, sum)$count,
                            out.degree = aggregate(count ~ Start.Station.ID, data = ny1610, sum)$count) %>%
  mutate(degree = in.degree + out.degree) %>%  
  left_join(ny1610[, c("Start.Station.Latitude", "Start.Station.Longitude", "Start.Station.ID")] %>% 
              filter(!duplicated(Start.Station.ID)), by = c("id" = "Start.Station.ID"))
names(ny1610.station.all)[match(c("Start.Station.Latitude", "Start.Station.Longitude"), names(ny1610.station.all))] = c("lat","lon")

#Useful for identifying stations in result output:
ny1610.station.names = (ny1610 %>% group_by(Start.Station.ID, ny1610$Start.Station.Name) %>% summarize(ID = first(Start.Station.ID), name = first(Start.Station.Name)))[,c("ID","name")]


#--------------------------------------------------------------------------------------------
#  1e. Trip count by hour (before removing weekends) ####
#   LA ####

la.overall = la %>% filter(start_station_id != end_station_id) %>%
  group_by(hour) %>% summarize(day_type = "overall", count = sum(count))

la.weekday =  la %>% filter(start_station_id != end_station_id & weekend == FALSE) %>%
  group_by(hour) %>% summarize(day_type = "weekday", count = sum(count))

la.weekend =  la %>% filter(start_station_id != end_station_id & weekend == TRUE) %>%
  group_by(hour) %>% summarize(day_type = "weekend", count = sum(count))

la.byhour = bind_rows(la.overall, la.weekday, la.weekend)

la_trip_count_by_hour <- ggplot(la.byhour, aes(x=hour, y = count, group = day_type, color = day_type)) + 
  geom_line(aes(linetype = day_type)) +
  ggtitle("Los Angeles") +
  theme_classic() +
  ylab("") +
  guides(color = "none", linetype = "none") + 
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 12)) +
  scale_x_continuous(breaks= seq(0, 23, by = 4), expand = expand_scale(0,0)) +
  scale_y_continuous(expand = expand_scale(0, 0))

ggsave("IMG/la_trip_count_by_hour.png", plot = la_trip_count_by_hour)

#   SF ####

sf.overall = sf %>% filter(start_station_id != end_station_id) %>%
  group_by(hour) %>% summarize(day_type = "overall", count = sum(count))

sf.weekday =  sf %>% filter(start_station_id != end_station_id & weekend == FALSE) %>%
  group_by(hour) %>% summarize(day_type = "weekday", count = sum(count))

sf.weekend =  sf %>% filter(start_station_id != end_station_id & weekend == TRUE) %>%
  group_by(hour) %>% summarize(day_type = "weekend", count = sum(count))

sf.byhour = bind_rows(sf.overall, sf.weekday, sf.weekend)

sf_trip_count_by_hour <- ggplot(sf.byhour, aes(x=hour, y = count, group = day_type, color = day_type)) + 
  geom_line(aes(linetype = day_type)) +
  ggtitle("San Francisco") +
  theme_classic() +
  guides(color = "none", linetype = "none") + 
  ylab("") +
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 12)) +
  scale_x_continuous(breaks= seq(0, 23, by = 4), expand = c(0, 0)) +
  scale_y_continuous(expand = expand_scale(0, 0))

ggsave("IMG/sf_trip_count_by_hour.png", plot = sf_trip_count_by_hour)


#   NY ####

ny.overall = ny1610 %>% filter(Start.Station.ID != End.Station.ID) %>%
  group_by(hour) %>% summarize(day_type = "  overall", trips = sum(count))

ny.weekday =  ny1610 %>% filter(Start.Station.ID != End.Station.ID & weekend == FALSE) %>%
  group_by(hour) %>% summarize(day_type = "  weekday", trips = sum(count))

ny.weekend =  ny1610 %>% filter(Start.Station.ID != End.Station.ID & weekend == TRUE) %>%
  group_by(hour) %>% summarize(day_type = "  weekend", trips = sum(count))

ny.byhour = bind_rows(ny.overall, ny.weekday, ny.weekend)

ny_trip_count_by_hour <- ggplot(ny.byhour, aes(x=hour, y = trips, group = day_type, color = day_type)) + 
  geom_line(aes(linetype = day_type)) +
  ggtitle("New York City") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "left",
        legend.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text = element_text(size = 14)) +
  scale_x_discrete(limits= seq(0, 23, by = 4), expand = c(0, 0)) +
  scale_y_continuous(expand = expand_scale(0, 0))

ggsave("IMG/ny_trip_count_by_hour.png", plot = ny_trip_count_by_hour)

#   Combined plot ####
combined_by_hour_plot = plot_grid(ny_trip_count_by_hour, sf_trip_count_by_hour, la_trip_count_by_hour,
                                  nrow = 1, rel_widths = c(.5,.3,.3))
ggsave("IMG/trip_count_by_hour.png", plot = combined_by_hour_plot, width = 10, height = 4)


#--------------------------------------------------------------------------------------------
#  1e. REMOVE WEEKENDS  ####
la.weekend = filter(la, weekend == TRUE)
sf.weekend = filter(sf, weekend == TRUE)
ny1610.weekend = filter(ny1610, weekend == TRUE) #10687, 21255, 378760
nrow(la.weekend); nrow(sf.weekend); nrow(ny1610.weekend) #
la = filter(la, weekend == FALSE)
sf = filter(sf, weekend == FALSE)
ny1610 = filter(ny1610, weekend == FALSE)
nrow(la); nrow(sf); nrow(ny1610) # 29443, 246157, 1172932
#----------------------------------------------------------------------------------------------
#  1f. Aggregate station info (no weekends) ####
##    LA  ----
la.station = data.frame(id = aggregate(count ~ start_station_id, data = la, sum)$start_station_id,
                        in.degree = aggregate(count ~ end_station_id, data = la, sum)$count,
                        out.degree = aggregate(count ~ start_station_id, data = la, sum)$count) %>%
  mutate(degree = in.degree + out.degree) %>%  
  left_join(la[, c("start_lat", "start_lon", "start_station_id")] %>% 
              filter(!duplicated(start_station_id)), by = c("id" = "start_station_id"))
names(la.station)[match(c("start_lat", "start_lon"), names(la.station))] = c("lat","lon")

##    SF  ----
sf.station = data.frame(id = aggregate(count ~ start_station_id, data = sf, sum)$start_station_id,
                        in.degree = aggregate(count ~ end_station_id, data = sf, sum)$count,
                        out.degree = aggregate(count ~ start_station_id, data = sf, sum)$count) %>%
  mutate(degree = in.degree + out.degree) %>%  
  left_join(sf[, c("start_lat", "start_lon", "start_station_id")] %>% 
              filter(!duplicated(start_station_id)), by = c("id" = "start_station_id"))

names(sf.station)[match(c("start_lat", "start_lon"), names(sf.station))] = c("lat","lon")


##    NY  ----
ny1610.station = data.frame(id = aggregate(count ~ Start.Station.ID, data = ny1610, sum)$Start.Station.ID,
                            in.degree = aggregate(count ~ End.Station.ID, data = ny1610, sum)$count,
                            out.degree = aggregate(count ~ Start.Station.ID, data = ny1610, sum)$count) %>%
  mutate(degree = in.degree + out.degree) %>%  
  left_join(ny1610[, c("Start.Station.Latitude", "Start.Station.Longitude", "Start.Station.ID")] %>% 
              filter(!duplicated(Start.Station.ID)), by = c("id" = "Start.Station.ID"))
names(ny1610.station)[match(c("Start.Station.Latitude", "Start.Station.Longitude"), names(ny1610.station))] = c("lat","lon")

#----------------------------------------------------------------------------------------------
#  1g.  Make time-sliced edgelist (no weekends) ----
##    LA ----
la_byhour = list()
for (i in 0:23) {
  la_hour = subset(la, hour == i)
  la_hour = aggregate(la_hour[,c("count")], #"duration"
                      by = data.frame(la_hour[,c("start_station_id", "end_station_id")]), #"start_lat", "end_lat", "start_lon", "end_lon"
                      FUN = sum)
  la_byhour[[i+1]] = la_hour
}
##    SF ----

sf_byhour = list()
for (i in 0:23) {
  sf_hour = subset(sf, hour == i)
  sf_hour = aggregate(sf_hour[,c("count")],
                      by = data.frame(sf_hour[,c("start_station_id", "end_station_id")]),
                      FUN = sum)
  sf_byhour[[i+1]] = sf_hour
}
##    NY ----
ny_byhour = list()
for (i in 0:23) {
  ny_hour = subset(ny1610, hour == i)
  ny_hour = aggregate(ny_hour[,c("count")], #"duration"
                      by = data.frame(ny_hour[,c("Start.Station.ID", "End.Station.ID")]), #"start_lat", "end_lat", "start_lon", "end_lon"
                      FUN = sum)
  ny_byhour[[i+1]] = ny_hour
}

#----------------------------------------------------------------------------------------------
# 1h. Save without weekends (for reproducibility and use by mixed model) ---- ----

write.csv(la[,c("start_time", "end_time", "start_station_id", "end_station_id",
                "start_lat", "start_lon", "end_lat", "end_lon")], "data/cleaned/LA16_cleaned_final_no_weekend.csv", row.names = F)

write.csv(sf[,c("Start.Date", "End.Date", "start_station_id", "end_station_id",
                "start_lat", "start_lon", "end_lat", "end_lon")], "data/cleaned/SF16_cleaned_final_no_weekend.csv", row.names = F)

write.csv(ny1610[,c("Start.Time", "Stop.Time", "Start.Station.ID", "End.Station.ID",
                    "Start.Station.Latitude", "Start.Station.Longitude", 
                    "End.Station.Latitude", "End.Station.Longitude")], "data/cleaned/Oct16_nyc_cleaned_final_no_weekend.csv", row.names = F)


################################################################################################
# 2. Data exploratory Analysis
#----------------------------------------------------------------------------------------------
# 2b. Correlations of in- and out-degree (all cities) ####
cor(la.station$in.degree, la.station$out.degree) #0.911692 
cor(sf.station$in.degree, sf.station$out.degree) #0.9799804
cor(ny1610.station$in.degree, ny1610.station$out.degree) #0.9929619

#----------------------------------------------------------------------------------------------
# 2b. SVD ####

# Helper function (svdplot2) ####
svdplot2 <- function(city.in, city.out, n.svd = 2, time = 24, title = "Singular Vectors") {
  city.in = city.in[,-1]
  city.out = city.out[,-1]
  city.svd.in  = svd(city.in)
  cat("\npercent variation explained (in):", round((city.svd.in$d)^2/sum(city.svd.in$d^2), 3)[1:n.svd])
  city.svd.out = svd(city.out)
  cat("\npercent variation explained (out):", round((city.svd.out$d)^2/sum(city.svd.out$d^2), 3)[1:n.svd])
  city.svd = data.frame("line" = rep(c(paste("in degree", 1:n.svd), paste("out degree", 1:n.svd)), each = time),
                        "vector" = rep(c(paste(" ", 1:n.svd), paste(" ", 1:n.svd)), each = time),
                        "hour" = rep(1:time, n.svd*2),
                        "value" = -c(city.svd.in$v[,1:n.svd], city.svd.out$v[,1:n.svd]),
                        "direction" = rep(c("in-degree", "out-degree"), each = time*n.svd))
  ggplot(city.svd, aes(x = hour, y = value, group = line, color = vector, lty = direction)) + 
    xlim(0,23) + 
    geom_line() + ggtitle(title) + facet_wrap(~vector) + theme_minimal()#theme(legend.title = element_blank())
}

#   LA ####

la.in = la %>% #filter(start_station_id != end_station_id) %>%
  group_by(start_station_id, hour) %>% summarize(count = n()) %>%
  spread(key = hour, value = count, fill = 0) 

la.out = la %>% #filter(start_station_id != end_station_id) %>%
  group_by(end_station_id, hour) %>% summarize(count = n()) %>%
  spread(key = hour, value = count, fill = 0) 

svdplot2(la.in, la.out, 2)

#   SF ####

sf.in = sf %>% #filter(start_station_id != end_station_id) %>%
  group_by(start_station_id, hour) %>% summarize(count = n()) %>%
  spread(key = hour, value = count, fill = 0) 

sf.out = sf %>% #filter(start_station_id != end_station_id) %>%
  group_by(end_station_id, hour) %>% summarize(count = n()) %>%
  spread(key = hour, value = count, fill = 0) 

svdplot2(sf.in, sf.out, 2)

#   NY ####
ny.out = ny1610 %>% #filter(Start.Station.ID != End.Station.ID) %>%
  group_by(Start.Station.ID, hour) %>% summarize(count = sum(count)) %>%
  spread(key = hour, value = count, fill = 0) 

ny.in =  ny1610 %>% #filter(Start.Station.ID != End.Station.ID) %>%
  group_by(End.Station.ID, hour) %>% summarize(count = sum(count)) %>%
  spread(key = hour, value = count, fill = 0) 

svdplot2(ny.in, ny.out, 2)

#   Plots in one place ####

ggsave("IMG/ny_svd.png", svdplot2(ny.in, ny.out, 2, title = "Singular Vectors — New York City"),
       device = "png", width = 10, height = 5)
ggsave("IMG/la_svd.png", svdplot2(la.in, la.out, 2, title = "Singular Vectors — Los Angeles"),
       device = "png", width = 10, height = 5)
ggsave("IMG/sf_svd.png", svdplot2(sf.in, sf.out, 2, title = "Singular Vectors — San Francisco"),
       device = "png", width = 10, height = 5)

################################################################################################
# 3. Run models
#----------------------------------------------------------------------------------------------
# 3a. Run mixed models ----
# Because some models models (particularly for new york) are time consuming to fit,
# some code used to fit the models is commented out and results are loaded below.
# Note we use python3, so on some systems this call would be modified e.g. 
# to ("python3 LosAngeles.py | tee LAoutput.txt")

# setwd("mixed_model_implementation")
# system("python LosAngeles.py | tee LAoutput.txt")
# system("python SanFrancisco.py | tee SFoutput.txt")
# system("python NewYork.py | tee NYoutput.txt") #very slow
# setwd("..")

#----------------------------------------------------------------------------------------------
# Run discrete models
#   LA ----

  # 2 blocks
  la2.3 = sbmt(la_byhour,  degreeCorrect = 3, directed = T, klPerNetwork = 50, maxComms = 2, seed = 1)
  saveRDS(la2.3, "discrete_model_results/la2_3.RDS")
  # 3 blocks
  la3.3 = sbmt(la_byhour,  degreeCorrect = 3, directed = T, klPerNetwork = 50, maxComms = 3, seed = 1)
  saveRDS(la3.3, "discrete_model_results/la3_3.RDS")
  # 4 blocks
  la4.3 = sbmt(la_byhour,  degreeCorrect = 3, directed = T, klPerNetwork = 50, maxComms = 4, seed = 1)
  saveRDS(la4.3, "discrete_model_results/la4_3.RDS")
  # 5 blocks
  la5.3 = sbmt(la_byhour,  degreeCorrect = 3, directed = T, klPerNetwork = 50, maxComms = 5, seed = 1)
  saveRDS(la5.3, "discrete_model_results/la5_3.RDS")


#   SF ----
  
 # 2 blocks
 sf2.3 = sbmt(sf_byhour,  degreeCorrect = 3, directed = T, klPerNetwork = 50, maxComms = 2, seed = 1)
 saveRDS(sf2.3, "discrete_model_results/sf2_3.RDS")
 # 3 blocks
 sf3.3 = sbmt(sf_byhour,  degreeCorrect = 3, directed = T, klPerNetwork = 50, maxComms = 3, seed = 1)
 saveRDS(sf3.3, "discrete_model_results/sf3_3.RDS")
 # 4 blocks
 sf4.3 = sbmt(sf_byhour,  degreeCorrect = 3, directed = T, klPerNetwork = 50, maxComms = 4, seed = 1)
 saveRDS(sf4.3, "discrete_model_results/sf4_3.RDS")

#   NY ----
 
 # slow so loaded instead of run
 
# # 2 blocks
# ny2.start.time = Sys.time()
# ny2.3 = sbmt(ny_byhour,  degreeCorrect = 3, directed = T, klPerNetwork = 50, maxComms = 2, seed = 1)
# saveRDS(ny2.3, "discrete_model_results/ny2_3.RDS")
# ny2.end.time = Sys.time()
# ny2.time = ny2.end.time - ny2.start.time # 7.370873 hours, 6.956191 hours, 5.146391 hours (3 trials), 5.58 another trial
# 
# # 3 blocks
# ny3.start.time = Sys.time()
# ny3.3 = sbmt(ny_byhour,  degreeCorrect = 3, directed = T, klPerNetwork = 50, maxComms = 3, seed = 1)
# saveRDS(ny3.3, "discrete_model_results/ny3_3.RDS")
# ny3.end.time = Sys.time()
# ny3.time = ny3.end.time - ny3.start.time #7.068 hours, 6.303427 hours (2 trials), found best score by iter 9, 6.63 another trial
# 
# # 4 blocks
# ny4.start.time = Sys.time()
# ny4.3 = sbmt(ny_byhour,  degreeCorrect = 3, directed = T, klPerNetwork = 50, maxComms = 4, seed = 1)
# saveRDS(ny4.3, "discrete_model_results/ny4_3.RDS")
# ny4.end.time = Sys.time()
# ny4.time = ny4.end.time - ny4.start.time #8.1524 hrs, 6.2345 hrs (2 trials), finds best by iter 12, 5.67 another trial
#    
# # 6 blocks
# ny6.start.time = Sys.time()
# ny6.3 = sbmt(ny_byhour,  degreeCorrect = 3, directed = T, klPerNetwork = 50, maxComms = 6, seed = 1)
# saveRDS(ny6.3, "discrete_model_results/ny6_3.RDS")
# ny6.end.time = Sys.time()
# ny6.time = ny6.end.time - ny6.start.time # 8.279234 hours, 6.052419 hours 2 runs, 6.49 another trial
 
################################################################################################
# 4. Load model results ####
#----------------------------------------------------------------------------------------------
#   LA ----

#     + continuous (AKA mixed) ----
# load gradient descent estimate
la_continuous.roles = read.csv("mixed_model_results/LA_2_roles.csv", sep = ",", header = T)
la_continuous.roles$X = as.character(la_continuous.roles$X)
la_continuous = left_join(la.station, la_continuous.roles, by = c("id"="X"))[,c("id", "lon", "lat", "X0", "X1")]
la_continuous.omega = lapply(read.csv("mixed_model_results/LA_2_omega.csv", sep = ",", header = T), matrix, 2, 2) # reconstruct result in same form as discrete

#     + discrete ----
# 2 blocks
la_discrete.results = readRDS("discrete_model_results/la2_3.RDS")
la_discrete = left_join(la.station, data.frame(id =  names(la_discrete.results$FoundComms), 
                                               role  =  (la_discrete.results$FoundComms)))
#write.csv(la_discrete, file = "discrete_model_results/la2_3_roles.csv", row.names = F)

# 3 blocks
la_discrete.results.3 = readRDS("discrete_model_results/la3_3.RDS")
la_discrete.3 = left_join(la.station, data.frame(id =  names(la_discrete.results.3$FoundComms), 
                                               role  =  (la_discrete.results.3$FoundComms)))
# 4 blocks
la_discrete.results.4 = readRDS("discrete_model_results/la4_3.RDS")
la_discrete.4 = left_join(la.station, data.frame(id =  names(la_discrete.results.4$FoundComms), 
                                                 role  =  (la_discrete.results.4$FoundComms)))
# 5 blocks
la_discrete.results.5 = readRDS("discrete_model_results/la5_3.RDS")
la_discrete.5 = left_join(la.station, data.frame(id =  names(la_discrete.results.5$FoundComms), 
                                                 role  =  (la_discrete.results.5$FoundComms)))

#   SF ----
#     + continuous (AKA mixed) ----

sf_continuous.roles = read.csv("mixed_model_results/SF_2_roles.csv", sep = ",", header = T)
sf_continuous = left_join(sf.station, sf_continuous.roles, by = c("id"="X"))[,c("id", "lon", "lat", "X0", "X1")]
sf_continuous.omega = lapply(read.csv("mixed_model_results/SF_2_omega.csv", sep = ",", header = T), matrix, 2, 2)

#     + discrete  ----
sf_discrete.results = readRDS("discrete_model_results/sf2_3.RDS")
sf_discrete = left_join(sf.station, data.frame(id =  as.integer(names(sf_discrete.results$FoundComms)), 
                                    role  =  (sf_discrete.results$FoundComms)), stringsAsFactors = F) 
write.csv(sf_discrete, file = "discrete_model_results/sf2_3_roles.csv", row.names = F)

#   NY ----
#     + continuous (AKA mixed) ----

# 2 blocks
ny_continuous.roles.2 = read.csv("mixed_model_results/NY_2_roles.csv", sep = ",", header = T)
ny_continuous.2 = left_join(ny1610.station, ny_continuous.roles.2, by = c("id"="X"))[,c("id", "lon", "lat", "X0", "X1")]
ny_continuous.omega.2 = lapply(read.csv("mixed_model_results/NY_2_omega.csv", sep = ",", header = T), matrix, 2, 2)

# 3 blocks
ny_continuous.roles.3 = read.csv("mixed_model_results/NY_3_roles.csv", sep = ",", header = T)
ny_continuous.3 = left_join(ny1610.station, ny_continuous.roles.3, by = c("id"="X"))[,c("id", "lon", "lat", "X0", "X1", "X2")]
ny_continuous.omega.3 = lapply(read.csv("mixed_model_results/NY_3_omega.csv", sep = ",", header = T), matrix, 3, 3)
  
# 4 blocks
ny_continuous.roles.4 = read.csv("mixed_model_results/NY_4_roles.csv", sep = ",", header = T)
ny_continuous.4 = left_join(ny1610.station, ny_continuous.roles.4, by = c("id"="X"))[,c("id", "lon", "lat", "X0", "X1", "X2", "X3")]
ny_continuous.omega.4 = lapply(read.csv("mixed_model_results/NY_4_omega.csv", sep = ",", header = T), matrix, 4, 4)



#     + discrete ----

# 2 blocks
ny_discrete.results.2 <- readRDS("discrete_model_results/ny2_3.RDS")
ny_discrete.2 = left_join(ny1610.station, data.frame(id =  as.integer(names(ny_discrete.results.2$FoundComms)), role  =  (ny_discrete.results.2$FoundComms), stringsAsFactors = F))
#write.csv(ny_discrete.2, file = "discrete_model_results/ny2_3_roles.csv", row.names = F)

# 3 blocks
ny_discrete.results.3 <- readRDS("discrete_model_results/ny3_3.RDS")
ny_discrete.3 = left_join(ny1610.station, data.frame(id =  as.integer(names(ny_discrete.results.3$FoundComms)), role  =  (ny_discrete.results.3$FoundComms), stringsAsFactors = F))
#write.csv(ny_discrete, file = "discrete_model_results/ny3_3_roles.csv", row.names = F)

# 4 blocks
ny_discrete.results.4 <- readRDS("discrete_model_results/ny4_3.RDS")
ny_discrete.4 = left_join(ny1610.station, data.frame(id =  as.integer(names(ny_discrete.results.4$FoundComms)), role  =  (ny_discrete.results.4$FoundComms), stringsAsFactors = F))
#write.csv(ny_discrete.4, file = "discrete_model_results/ny4_3_roles.csv", row.names = F)

# 6 blocks
ny_discrete.results.6 <- readRDS("discrete_model_results/ny6_3.RDS")
ny_discrete.6 = left_join(ny1610.station, data.frame(id =  as.integer(names(ny_discrete.results.6$FoundComms)), role  =  (ny_discrete.results.6$FoundComms), stringsAsFactors = F))
#write.csv(ny_discrete.6, file = "discrete_model_results/ny6_3_roles.csv", row.names = F)


################################################################################################
# 5. Plots ####
#----------------------------------------------------------------------------------------------
# 5a. Maps for backgrounds ####

# Code to generate the maps loaded below requires an individual google api key
# register_google(key = "") # <- your key here
# 
# ny_background = get_map(location = c(-74.25617, 40.93617, -73.70015, 40.51111), #min(nyzd.spdf$long), max(nyzd.spdf$lat)+.02, max(nyzd.spdf$long), min(nyzd.spdf$lat)+.02
#                         zoom = 12, maptype = "roadmap")
# 
# la_background = get_map(location = c(-118.2529, 34.0486), zoom = 14,
#                         maptype = "roadmap", scale = 1)
# 
# sf_background = get_map(location = c(-122.4037, 37.7860), zoom = 14,
#                         maptype = "roadmap", scale = 1)
#saveRDS(ny_background, "ny_background.RDS")
#saveRDS(la_background, "la_background.RDS")
#saveRDS(sf_background, "sf_background.RDS")

ny_background = readRDS("data/maps/ny_background.RDS")
la_background = readRDS("data/maps/la_background.RDS")
sf_background = readRDS("data/maps/sf_background.RDS")

#----------------------------------------------------------------------------------------------
#   LA ----
#       pre-process ####

# continuous 

names(la_continuous) = c("id", "lon", "lat", "X1", "X2")
la_continuous$role = la_continuous$X1/(la_continuous$X1 + la_continuous$X2)
la_continuous$C_total = la_continuous$X1 + la_continuous$X2
# assign roles
left_join(la_continuous, la_station_info, by = c("id"="Station.ID")) %>% arrange(X1) %>% select(X1, Station.Name)
la_continuous_labels = c("work", "home")

# discrete

#should include colnames c("id", "in.degree", "out.degree", "degree", "lon", "lat", "role")
la_discrete$role = as.factor(la_discrete$role)
left_join(la_discrete, la_station_info, by = c("id"="Station.ID")) %>% arrange(role) %>% select(role, Station.Name)
# assign roles
la_discrete_labels = c("work", "home")
levels(la_discrete$role) = la_discrete_labels

#     + TDD and TDMM over map ----
#         discrete  ####

plot_la_discrete = ggmap(la_background) +
  geom_point(data=la_discrete, shape = 21, aes(x=lon, y=lat, fill = role,
                  size = degree/max(degree)), stroke = .2) +
  ggtitle("Los Angeles, Discrete Membership") +
  scale_fill_manual(values = c("black", "white"),
                    guide = guide_legend(reverse=TRUE, order = 1)) +
  scale_size(range = c(1,5), TeX(paste0("degree","\\textbf{/}","max(degree)"))) + # TeX("$\\frac{degree}{max(degree)}$") 
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())

plot_la_discrete 

#         continuous  ####
plot_la_continuous = ggmap(la_background) +
  geom_point(data=la_continuous, shape = 21, aes(x=lon, y=lat, fill=role, size = C_total), 
             show.legend = T, stroke = .2) +
  ggtitle("Los Angeles, Mixed Membership") +
  scale_fill_gradient(low = "white", high = "black", limits = c(0,1), breaks = c(.1,.9),
                      labels=c("home","work"), 
                      guide = guide_colorbar(ticks = FALSE, reverse = T)) + 
  scale_size(range = c(1,5), name = "C total") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())
plot_la_continuous

#         save joint plot ####

plot_la_mixed_discrete = plot_grid(plot_la_continuous, plot_la_discrete, rel_widths = c(4,4.85))

ggsave("IMG/LA_mixed_discrete.png", plot_la_mixed_discrete, width = 9, height = 3.5)


#     +  Zoning ----

if (!file.exists("data/zoning/LA_Zoning")) {
  unzip("data/zoning/LA_Zoning.zip", exdir = "data/zoning")
}

lazd.shp <- readOGR("data/zoning/LA_Zoning/Zoning.shp")
#lazd.dbf <- readOGR("data/zoning/LA_Zoning/Zoning.dbf")

raster::crs(lazd.shp)

#         convert to a standard coordinate reference system with longlat ####

lazd.latlong <- spTransform(lazd.shp, CRS("+proj=longlat +datum=WGS84"))
lazd.spdf <- broom:::tidy.SpatialPolygonsDataFrame(lazd.latlong)

lazd.spdf$id = as.numeric(lazd.spdf$id) + 1 

#         what are the zoning types? ####

lazd.spdf$ZONEDIST = factor(lazd.spdf$id)
levels(lazd.spdf$ZONEDIST) = lazd.latlong$ZONE_CMPLT
lazd.spdf$zone = as.character(lazd.spdf$ZONEDIST)

#         subset to our bounding box of interest ####

lazd.spdf = subset(lazd.spdf, long > -118.28 & long < -118.2)
lazd.spdf = subset(lazd.spdf, lat > 34 & lat < 34.7)
lazd.spdf$ZONEDIST = droplevels(lazd.spdf$ZONEDIST)

lazd.spdf$zone = as.character(lazd.spdf$ZONEDIST)
tmp = str_extract(lazd.spdf$zone, "[A-Z]{1}[A-Z0-9]{1}")
tmp[which(is.na(tmp))] = "P" 
lazd.spdf$zone = tmp
lazd.spdf$zone_type = sapply(lazd.spdf$zone, function(x) {
  if (startsWith(x ,"C")) {return("residential/\ncommercial")}
  else if (startsWith(x ,"M1") | startsWith(x ,"M2")) {return("indust. light")}
  else if (startsWith(x ,"M")) {return("indust. heavy")}
  else if (startsWith(x ,"PF")) {return("public facilities")}
  else if (startsWith(x ,"R")) {return("residential")}
  #else if (startsWith(x ,"A1")) {return("park")}
  #else if (startsWith(x ,"A2")) {return("park")}
  else if (startsWith(x ,"AD")) {return("Union Station")}
  else {return("other")}
})

#         zone colors: https://www.r-bloggers.com/palettes-in-r/ ----

# M1/2 , M3/R , other ,  public facilities , res, res/com,  union station 
lapal = brewer.pal(7,"RdYlBu"); lapal = rev(lapal) 
lapal[3] = "#D9D1CE"; lapal[4] = "green4"; lapal[5] = "red"; lapal[6] = "violet"; lapal[7] = "yellow" 
#closer to ny palette:
#lapal = c(brewer.pal(4,"BuGn")[c(1,3)], "#668cb2", "#D9D1CE", brewer.pal(4,"BuGn")[c(2,4)], "#003366")

#         plot (discrete) ####

la_zone_discrete = ggplot() + geom_polygon(data = lazd.spdf,
                                           aes(x = long, y = lat, group = group, fill = zone_type)) +
  coord_map(xlim= c(-118.2748, -118.225), ylim = c(34.02, 34.07)) +
  scale_fill_manual(values = alpha(lapal,.5), name = "zone type") +
  ggtitle("Discrete Roles versus LA Zoning") +
  geom_point(data = la_discrete, aes(x = lon, y = lat, color = role, size = degree/max(degree))) +
  scale_color_manual(values = c("black", "white")) +
  geom_point(data = la_discrete, shape = 1, aes(x = lon, y = lat, size = degree/max(degree))) + 
  guides(colour = guide_legend(override.aes = list(shape = 21, fill = c("white", "black"), 
                                                   color = "black"), order = 1, reverse = TRUE)) +
  scale_size(guide = 'none') + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())
la_zone_discrete

#         plot (continuous) ####

la_zone_continuous = ggplot() + geom_polygon(data = lazd.spdf,
                                             aes(x = long, y = lat, group = group, fill = zone_type)) +
  coord_map(xlim= c(-118.2748, -118.225),  ylim = c(34.02, 34.07)) +
  scale_fill_manual(values = alpha(lapal,.5), name = "zone type") +
  geom_point(data = la_continuous, aes(x = lon, y = lat, color = role, size = C_total)) +
  ggtitle("Mixed-Membership Roles versus LA Zoning") +
  scale_color_gradient(low = "white", high = "black", limits = c(0,1), breaks = c(.1,.9),
                       labels=c("home","work"),
                       guide = guide_colorbar(ticks = FALSE, reverse = TRUE, order = 1)) +
  geom_point(data = la_continuous, shape = 1, aes(x = lon, y = lat, size = C_total)) + #add for outline
  scale_size(range = c(2,7), guide = 'none') +
  #xlab("longitude") + ylab("latitude") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())
la_zone_continuous

ggsave("IMG/la_continuous_zones.png", plot_grid(la_zone_continuous), height = 6, width = 6)

#     + Omegas ----
# time-dependent parameters

la.plot.order = c("home to home", "home to work", "work to home", "work to work")

# discrete
results = data.frame(t(sapply(la_discrete.results$EdgeMatrix, unlist)))
names(results) = c(outer(la_discrete_labels, la_discrete_labels, paste, sep = " to "))
results$hour = 0:(nrow(results)-1)
results = melt(results, id.vars = "hour", value.name = "count")
# reorder levels for plot
results$variable = factor(results$variable, levels = la.plot.order)

la_omega_discrete = ggplot(results, aes(x = hour, y = count, color = 1)) + geom_line() +
  facet_wrap(.~variable) +
  xlim(0,23) +
  guides(color = "none") +
  ggtitle("Time-Dependent Block Parameters for LA Discrete Model") +
  theme_bw() + 
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, size = 14), strip.text = element_text(size = 12))
la_omega_discrete

ggsave("IMG/la_omega_discrete.png", plot_grid(la_omega_discrete), width = 6, height = 6)

# continuous
results = data.frame(t(sapply(la_continuous.omega, unlist)))
names(results) = c(outer(la_continuous_labels, la_continuous_labels, paste, sep = " to "))
results$hour = 0:(nrow(results)-1)
results = melt(results, id.vars = "hour", value.name = "count")
results$variable = factor(results$variable, levels = la.plot.order)

la_omega_continuous = ggplot(results, aes(x = hour, y = count, color = 1)) + geom_line() + facet_wrap(.~variable) +
  guides(color = "none") +
  ggtitle("Time-Dependent Block Parameters for LA Mixed Model") +
  theme_bw() + 
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, size = 14), strip.text = element_text(size = 12))
la_omega_continuous

ggsave("IMG/la_omega.png", plot_grid(la_omega_continuous, la_omega_discrete), width = 12, height = 6)

#   SF ----
#       pre-process ####

# continuous
names(sf_continuous) = c("id", "lon", "lat", "X1", "X2")
sf_continuous$role = sf_continuous$X2/(sf_continuous$X1 + sf_continuous$X2)
sf_continuous$C_total = sf_continuous$X1 + sf_continuous$X2
# assign roles
sum(sf_continuous$X1 > sf_continuous$X2) # Use the fact that there are more work stations
sf_continuous_labels = c("home", "work")

# discrete
#should include colnames c("id", "in.degree", "out.degree", "degree", "lon", "lat", "role")
sf_discrete$role = as.factor(sf_discrete$role)
# assign roles
table(sf_discrete$role) # Use the fact that there are more work stations
sf_discrete_labels = c("work", "home") 
levels(sf_discrete$role) = sf_discrete_labels

#     + TDD and TDMM over map ----
#         discrete ####

plot_sf_discrete = ggmap(sf_background) +
  geom_point(data=sf_discrete, shape=21, aes(x=lon, y=lat, fill = role, pch = role,
                                             size = degree/max(degree)), stroke = .2) +
  ggtitle("San Francisco, Discrete Membership") +
  scale_fill_manual(values = c("black","white"),
                    guide = guide_legend(reverse=TRUE, order = 1)) +
  scale_size(range = c(1,5), TeX(paste0("degree","\\textbf{/}","max(degree)"))) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())

#         continuous ####

plot_sf_continuous = ggmap(sf_background) +
  geom_point(data=sf_continuous, shape = 21,  aes(x=lon, y=lat, fill = role, size = C_total), stroke = .2) +
  ggtitle("San Francisco, Mixed Membership") +
  scale_fill_gradient(low = "white", high = "black", limits = c(0,1), breaks = c(.1,.9),
                      labels=c("home","work"), 
                      guide = guide_colorbar(reverse = TRUE, ticks = FALSE, order = 1)) + 
  scale_size(range = c(1,5), name = "C total") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())
plot_sf_continuous

#         save joint plot ####

plot_sf_mixed_discrete = plot_grid(plot_sf_continuous, plot_sf_discrete, rel_widths = c(4,4.75))

ggsave("IMG/SF_mixed_discrete.png", plot_sf_mixed_discrete, width = 9, height  = 3.5)

#     + Omegas ----

sf.plot.order = c("home to home", "home to work", "work to home", "work to work")

# discrete
results = data.frame(t(sapply(sf_discrete.results$EdgeMatrix, unlist)))
names(results) = c(outer(sf_discrete_labels, sf_discrete_labels, paste, sep = " to "))
results$hour = 0:(nrow(results)-1)
results = melt(results, id.vars = "hour", value.name = "count")
results$variable = factor(results$variable, levels = sf.plot.order)

sf_omega_discrete = ggplot(results, aes(x = hour, y = count, color = 1)) + 
  geom_line() + facet_wrap(.~variable) +
  guides(color = "none") +
  ggtitle("Time-Dependent Block Parameters for SF Discrete Model") +
  theme_bw() +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5), 
        strip.text = element_text(size = 12))
sf_omega_discrete

# continuous
results = data.frame(t(sapply(sf_continuous.omega, unlist)))
names(results) = c(outer(sf_continuous_labels, sf_continuous_labels, paste, sep = " to "))
results$hour = 0:(nrow(results)-1)
results = melt(results, id.vars = "hour", value.name = "count")
results$variable = factor(results$variable, levels = la.plot.order)

sf_omega_continuous = ggplot(results, aes(x = hour, y = count, color = 1)) + geom_line() + facet_wrap(.~variable) +
  guides(color = "none") +
  ggtitle("Time-Dependent Block Parameters for SF Mixed Model") +
  theme_bw() +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5), strip.text = element_text(size = 12))
sf_omega_continuous

ggsave("IMG/sf_omega.png", plot_grid(sf_omega_continuous, sf_omega_discrete), width = 12, height = 6)

#   NY ----
#      pre-process ####
#         discrete ----

# 2 blocks
ny_discrete.2$role = as.factor(ny_discrete.2$role)
levels(ny_discrete.2$role) = c("Manhattan", "Brooklyn")

# 3 blocks
#should include colnames c("id", "in.degree", "out.degree", "degree", "lon", "lat", "role")
ny_discrete.3$role = as.factor(ny_discrete.3$role)
# assign roles
left_join(ny_discrete.3, ny1610.station.names, by = c("id"="ID"))
ny_discrete_labels.3 = c("Manhattan (work)", "Brooklyn", "Manhattan (home)")
levels(ny_discrete.3$role) = ny_discrete_labels.3

# 4 blocks
ny_discrete.4$role = as.factor(ny_discrete.4$role)
levels(ny_discrete.4$role) = c("Brooklyn", "Manhattan (home)", "Manhattan (midtown)", "Manhattan (work)")

# 6 blocks
ny_discrete.6$role = as.factor(ny_discrete.6$role)

# assign roles
left_join(ny_discrete.6, ny1610.station.names, by = c("id"="ID"))
ny_discrete_labels.6 = c(as.character(0:5)) #TBD
levels(ny_discrete.6$role) = ny_discrete_labels.6

#         continuous ----
names(ny_continuous.3) = c("id", "lon", "lat", "X0", "X1","X2")
ny_continuous.3$C_total = ny_continuous.3$X0 + ny_continuous.3$X1 + ny_continuous.3$X2
ny_continuous.3$role = rgb(ny_continuous.3[,c("X0", "X1", "X2")]/ny_continuous.3$C_total)
# Could also use this to assign shape as dominant class:
role.colors.3 = c("#fc8d59", "#ffffbf", "#91bfdb")
ny_continuous.3$discretize = role.colors.3[apply(ny_continuous.3[,c("X0", "X1", "X2")]/ny_continuous.3$C_total, 1, which.max)]

# continuous (2)
names(ny_continuous.2) = c("id", "lon", "lat", "X0", "X1")
ny_continuous.2$C_total = ny_continuous.2$X0 + ny_continuous.2$X1
ny_continuous.2$role =  ny_continuous.2$X1/(ny_continuous.2$X1 + ny_continuous.2$X0)

# continuous (4)
names(ny_continuous.4) = c("id", "lon", "lat", "X0", "X1", "X2", "X3")
ny_continuous.4$C_total = rowSums(ny_continuous.4[,7:4])
ny_continuous.4$role = rgb(ny_continuous.4[,7:4]/ny_continuous.4$C_total)
ny_continuous.4$discretize = apply(ny_continuous.4[,c("X0", "X1", "X2","X3")]/ny_continuous.4$C_total, 1, which.max)


#    + TDD and TDMM over map ----

# 2 blocks - splits into manhattan, brooklyn
# 4 blocks - splits into "Brooklyn", "Manhattan (midtown)", "Manhattan (home)", "Manhattan (work)"

#         ny discrete ####

#               2 blocks ----

plot_ny_discrete.2 = ggmap(ny_background) +
  geom_point(data = ny_discrete.2, aes(x = lon, y = lat, fill = role,
                                       shape = role, size = degree/max(degree)), color = "black", alpha = .7) +
  ggtitle("New York, Discrete Membership") +
  xlim(-74.02,-73.93) + ylim(40.65, 40.8) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(size = guide_legend(order = 2)) +
  scale_fill_manual(name = "role", 
                    values = c("white", "black"),
                    guide = guide_legend(order = 1)) + 
  scale_shape_manual(name = "role", 
                     values = c(22, 24), 
                     guide = guide_legend(order = 1)) +
  scale_size(range = c(1,4), name = TeX(paste0("degree","\\textbf{/}","max(degree)"))) +
  #xlab("longitude") + ylab("latitude")
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())

plot_ny_discrete.2

ggsave("IMG/ny_discrete_2.png", plot_ny_discrete.2, width = 10, height = 6, units = "in")


#               3 blocks ----
#change level order for plot
ny_discrete.3$role = factor(ny_discrete.3$role, levels = ny_discrete_labels.3[c(3,1,2)])

plot_ny_discrete.3 = ggmap(ny_background) +
  geom_point(data = ny_discrete.3, aes(x = lon, y = lat, fill = role,
                                     shape = role, size = degree/max(degree)), color = "black", alpha = .7) +
  ggtitle("New York, Discrete Membership") +
  xlim(-74.02,-73.93) + ylim(40.65, 40.8) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(size = guide_legend(order = 2)) +
  scale_fill_manual(name = "role", 
                    values = c("white", "black", "gray"),
                    guide = guide_legend(order = 1)) + 
  scale_shape_manual(name = "role", 
                     values = c(22, 24, 21), guide = guide_legend(order = 1)) +
  scale_size(range = c(1,4), name = TeX(paste0("degree","\\textbf{/}","max(degree)"))) +
  #xlab("longitude") + ylab("latitude")
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())

plot_ny_discrete.3

#               4 blocks ----

plot_ny_discrete.4 = ggmap(ny_background) +
  geom_point(data = ny_discrete.4, aes(x = lon, y = lat, fill = role,
                                       shape = role, size = degree/max(degree)), color = "black", alpha = .7) +
  ggtitle("New York, Discrete Membership") +
  xlim(-74.02,-73.93) + ylim(40.65, 40.8) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(size = guide_legend(order = 2)) +
  scale_fill_manual(name = "role", 
                    values = c("red", "blue", "green", "black"),
                    labels = c("Manhattan (midtown)", "Brooklyn", "Manhattan (work)", "Manhattan (home)"),
                    guide = guide_legend(order = 1)) + 
  scale_shape_manual(name = "role", 
                     values = c(22, 21, 23, 24), 
                     labels = c("Manhattan (midtown)", "Brooklyn", "Manhattan (work)", "Manhattan (home)"),
                     guide = guide_legend(order = 1)) +
  scale_size(range = c(1,4), name = TeX(paste0("degree","\\textbf{/}","max(degree)"))) +
  #xlab("longitude") + ylab("latitude")
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())

plot_ny_discrete.4


#               6 blocks ----

#change level order for plot
ny_discrete.6$role = factor(ny_discrete.6$role, levels = ny_discrete_labels.6[c(1,2,3,4,5,6)])

plot_ny_discrete.6 = ggmap(ny_background) +
  geom_point(data = ny_discrete.6, aes(x = lon, y = lat, fill = role,
                                       shape = role, size = degree/max(degree)), color = "black", alpha = .7) +
  ggtitle("New York, Discrete Membership") +
  xlim(-74.02,-73.93) + ylim(40.65, 40.8) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(size = guide_legend(order = 2)) +
  scale_size(range = c(1,4), name = TeX(paste0("degree","\\textbf{/}","max(degree)"))) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())

plot_ny_discrete.6

#         ny continuous ####

#               2 blocks ----

plot_ny_continuous = ggmap(ny_background) +
  geom_point(data=ny_continuous.2, aes(x=lon, y=lat, fill = role, size = C_total), shape = 21) +
  ggtitle("New York, Mixed Membership") +
  scale_fill_gradient(low = "black", high = "white", limits = c(0,1), breaks = c(.1,.9),
                      labels=c("Brooklyn","Manhattan"),
                      guide = guide_colorbar(ticks = FALSE)) + 
  scale_size(range = c(1,4), name = "C total") +
  xlim(-74.02,-73.93) + ylim(40.65, 40.8) +
  theme_classic() +
  #xlab("longitude") + ylab("latitude")
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())
plot_ny_continuous

ggsave("IMG/ny_mixed_2.png", plot_ny_continuous, width = 10, height = 6, units = "in")

#               3 blocks ----

role.colors.3 = c("#fc8d59", "#ffffbf", "#91bfdb")# <- selected using divergent, colorsafe scheme: http://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=3
plot_ny_continuous.3 = ggmap(ny_background) +
  geom_point(data=ny_continuous.3, shape = 21,
             aes(x=lon, y=lat, alpha = X0/C_total, size = C_total, fill = I(role.colors.3[1]) )) +
  geom_point(data=ny_continuous.3, shape = 21,
             aes(x=lon, y=lat, alpha = X1/C_total, size = C_total, fill = I(role.colors.3[2]) )) +
  geom_point(data=ny_continuous.3, shape = 21,
             aes(x=lon, y=lat, alpha = X2/C_total, size = C_total, fill = I(role.colors.3[3]) )) +
  guides(size = guide_legend(override.aes=list(fill = "black"), order = 2), 
         fill = guide_legend(order = 1, reverse = TRUE)) +
  #guides(alpha = guide_legend(override.aes=list(fill = "white")))) +
  scale_fill_manual(values = role.colors.3,
                    labels = c("Brooklyn", "Manhattan (work)", "Manhattan (home)"),
                    name = "role classes") + 
  scale_alpha_continuous(guide=FALSE) + 

  ggtitle("New York, Mixed Membership") +
  scale_size(range = c(1,5), name = "C total") +
  xlim(-74.02,-73.93) + ylim(40.65, 40.8) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())
plot_ny_continuous.3

#               4 blocks ----
role.colors.4 = hue_pal()(4)
plot_ny_continuous.4 = ggmap(ny_background) +
#  geom_point(data=ny_continuous.4, shape = 21, aes(x=lon, y=lat, fill = role.colors.4[1], alpha = X0/C_total, size = C_total)) +
#  geom_point(data=ny_continuous.4, shape = 21, aes(x=lon, y=lat, fill = role.colors.4[2], alpha = X1/C_total, size = C_total)) +
#  geom_point(data=ny_continuous.4, shape = 21, aes(x=lon, y=lat, fill = role.colors.4[3], alpha = X2/C_total, size = C_total)) +
#  geom_point(data=ny_continuous.4, shape = 21, aes(x=lon, y=lat, fill = role.colors.4[4], alpha = X3/C_total, size = C_total)) +
  geom_point(data=ny_continuous.4, aes(x=lon, y=lat, size = C_total,
                                     color  = I(role), fill = as.factor(discretize),
             text = paste(round(X0,3), round(X1,3), round(X2,3), round(X3,3), sep = ",")), 
             shape = 16, alpha = .8) +
  scale_fill_manual(values = c("green","blue","red","black"),
                    labels = c("Manhattan (home)", "Manhattann (work)", "Manhattan (home, lower)", "Brooklyn"))  +
  guides(fill = guide_legend(override.aes=list(shape = 21, stroke = 0, size = 3), 
                             title = "role", order = 1, reverse = TRUE)) +
  ggtitle("New York, Mixed Membership") +
  scale_size(range = c(1,5), name = "C total") +
  xlim(-74.02,-73.93) + ylim(40.65, 40.8) +
  theme_classic() +
  #xlab("longitude") + ylab("latitude")
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())
plot_ny_continuous.4

ggsave("IMG/ny_mixed_4.png", plot_ny_continuous.4, width = 6, height = 10)

#         save joint plot ####

plot_ny_mixed_discrete = plot_grid(plot_ny_continuous.3, plot_ny_discrete.3, rel_widths = c(4.5,5))

ggsave("IMG/NY_mixed_discrete.png", plot_ny_mixed_discrete, width = 10, height = 6, units = "in")

#    + Zoning ----
#         load zoning maps & convert to lat long coords ####

# new york zoning districts released 2018: https://www1.nyc.gov/site/planning/data-maps/open-data/dwn-gis-zoning.page
# co files are for commercial overlays. "commercial overlays are located in residential zoning districts. It means that you can have a commercial use in a residential area." https://jorgefontan.com/commercial-overlay-nyc-zoning/
# CRS = coordinate reference system 
# LCC = lambert conformal conic projection

#         load ####
if (!file.exists("data/zoning/nycgiszoningfeatures_201712shp")) {
  unzip("data/zoning/nycgiszoningfeatures_201712shp.zip", exdir = "data/zoning")
}

nyzd.shp <- readOGR("data/zoning/nycgiszoningfeatures_201712shp/nyzd.shp")
nyzd.dbf = read.dbf("data/zoning/nycgiszoningfeatures_201712shp/nyzd.dbf") #has ZONEDIST types
# look at given projection system (this site helpful: http://rspatial.org/spatial/rst/6-crs.html)
raster::crs(nyzd.shp)
#         convert to a standard coordinate reference system with longlat ####
nyzd.latlong <- spTransform(nyzd.shp, CRS("+proj=longlat +datum=WGS84"))

nyzd.spdf <- broom:::tidy.SpatialPolygonsDataFrame(nyzd.latlong)
head(nyzd.spdf)
#IMPORANT TO DO THIS NEXT ID STEP. not exactly sure why. 
nyzd.spdf$id = as.numeric(nyzd.spdf$id) + 1 

#         what are the zoning types? ####
#Fortunately, the ordering of "id" is the same as that stored in shape@data$region: https://stackoverflow.com/questions/40576457/keep-region-names-when-tidying-a-map-using-broom-package
nyzd.spdf$ZONEDIST = factor(nyzd.spdf$id)
levels(nyzd.spdf$ZONEDIST) = nyzd.latlong$ZONEDIST
nyzd.spdf$zone = as.character(nyzd.spdf$ZONEDIST)
nyzd.spdf$zone = sapply(nyzd.spdf$zone, function(x) {
  if (startsWith(x ,"BPC")) return("Battery Park City")
  if (startsWith(x ,"C")) return("commercial")
  if (startsWith(x ,"M")) return("manufaturing")
  if (startsWith(x ,"P")) return("park")
  if (startsWith(x ,"R")) return("residential")
})
# BPC, C, M, R, PARK - Battery park city, commercial, manufacturing, residential, park

# nyzd.co should give information about commercial overlays

#         pre-process
#         plot (discrete) ####

# create fill2 to have two non-contradictory fill scales
fill2 = (ny_discrete.3 %>% mutate(fill2 = c("white", "black", "gray")[role]))$fill2

plot_ny_discrete_zone = ggplot() + 
  geom_polygon(data = nyzd.spdf, aes(x = long, y = lat, group = group, fill = zone)) +
  scale_fill_brewer(type = "seq", palette = 2, direction = -1, name = "zone type") + 
  theme_classic() + 
  geom_point(data = ny_discrete.3, 
             aes(x = lon, y = lat, shape = role, size = degree/max(degree)),
             fill = fill2, alpha = .7) +
  ggtitle("Discrete Roles versus NYC Zoning") +
  scale_shape_manual(values = c(22, 24, 21)) +
  guides(shape = guide_legend(override.aes = list(fill = c("white", "black", "gray")), order = 1),
         size = "none") +
  scale_size(range = c(.5,2.5), name = TeX(paste0("degree","\\textbf{/}","max(degree)"))) +
  xlim(-74.02, -73.925) + ylim(40.65, 40.81) +
  coord_map() + #<- takes away angle
  #xlab("longitude") + ylab("latitude") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())

plot_ny_discrete_zone

ggsave("IMG/ny_discrete_zones.png", plot_ny_discrete_zone, width = 5, height = 8) 

# show same on regular plot for comparison

#         plot (continuous) ####

plot_ny_continuous_zone = ggplot() + geom_polygon(data = nyzd.spdf, aes(x = long, y = lat, group = group, fill = zone)) +
  scale_fill_brewer(type = "seq", palette = 2, direction = -1, name = "zone type") + 
  theme_classic() + 
  geom_point(data = ny_continuous.3, aes(x = lon, y = lat), color = ny_continuous.3$role) +
  ggtitle("Continuous Roles versus NYC Zoning Types") +
  #guides(color = "none") +
  scale_size(range = c(1,5)) +
  xlim(-74.02, -73.925) + ylim(40.65, 40.85) +
  coord_map() + #<- takes away angle
  #xlab("longitude") + ylab("latitude") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())
plot_ny_continuous_zone

ggsave("IMG/ny_continuous_zones.png", plot_ny_continuous_zone) 
#     + Omegas ----

# 3 - block

ny.plot.order = c("M (home) to M (home)", "M (home) to M (work)",  "M (home) to BK", "M (work) to M (home)",
                    "M (work) to M (work)", "M (work) to BK" , "BK to M (home)", "BK to M (work)", "BK to BK")

# discrete
  results = data.frame(t(sapply(ny_discrete.results.3$EdgeMatrix, unlist)))
  results$hour = 0:(nrow(results)-1)
  # Had to abbreviate names from ny_discrete_labels.3
  colnames(results) = c(outer(c("M (work)", "BK", "M (home)"), c("M (work)", "BK", "M (home)"), 
                              paste, sep = " to "), "hour")
  # ny_discrete %>% group_by(role) %>% summarize(sum(degree)) # a check
  results = melt(results, id.vars = "hour", value.name = "count")
  results$variable = factor(results$variable, levels = ny.plot.order)

  ny_omega_discrete = ggplot(results, aes(x = hour, y = count, color = 1)) + 
    geom_line() + facet_wrap(.~variable) + 
    guides(color = "none") +
    ggtitle("Time-Dependent Block Parameters for NY Discrete Model") +
    theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5), 
                     strip.text = element_text(size = 12))
  ny_omega_discrete
  ggsave("IMG/ny_omega_discrete.png", ny_omega_discrete, width =12, height = 8)

# continuous 
  
  results = data.frame(t(sapply(ny_continuous.omega.3, unlist)))
  results$hour = 0:(nrow(results)-1)
  # Had to abbreviate names from ny_con 
  names(results) = c(outer(c("M (work)", "M (home)","BK"),  c("M (work)", "M (home)","BK"), 
                           paste, sep = " to "), "hour")
  results = melt(results, id.vars = "hour", value.name = "count")
  # reorder for plot
  results$variable = factor(results$variable, levels = ny.plot.order)

  ny_omega_continuous = ggplot(results, aes(x = hour, y = count, color = 1)) +
    geom_line() + facet_wrap(.~variable) +
    guides(color = "none") +
    ggtitle("Time-Dependent Block Parameters for NY Mixed Model") +
    theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5), 
                       strip.text = element_text(size = 12))

  ny_omega_continuous
  ggsave("IMG/ny_omega.png", plot_grid(ny_omega_continuous, ny_omega_discrete), width = 13, height = 6)

################################################################################################
# 6d. STATIC models (time-independent) ----

#     run the static models (load saved output) ####

#      pre-process all cities 
la.tmp = data.frame(la %>% group_by(start_station_id, end_station_id) %>% summarize(count = sum(count)), stringsAsFactors = F)
sf.tmp = data.frame(sf %>% group_by(start_station_id, end_station_id) %>% summarize(count = sum(count)), stringsAsFactors = F)
ny.tmp = data.frame(ny1610 %>% group_by(Start.Station.ID, End.Station.ID) %>% summarize(count = sum(count)), stringsAsFactors = F)

la2.T  = sbmt(list(la.tmp), maxComms = 2, degreeCorrect = 3, directed = T, klPerNetwork = 50, seed = 1)
sf2.T  = sbmt(list(sf.tmp), maxComms = 2, degreeCorrect = 3, directed = T, klPerNetwork = 50, seed = 1)
#ny2.T  = sbmt(list(ny.tmp), maxComms = 2, degreeCorrect = 3, directed = T, klPerNetwork = 50, seed = 1) <- slow, load pre-run
#ny3.T  = sbmt(list(ny.tmp), maxComms = 3, degreeCorrect = 3, directed = T, klPerNetwork = 50, seed = 1)
  saveRDS(la2.T, "discrete_model_results/la2_T.RDS")
  saveRDS(sf2.T, "discrete_model_results/sf2_T.RDS")
  #saveRDS(ny2.T, "discrete_model_results/ny2_T.RDS")
  #saveRDS(ny3.T, "discrete_model_results/ny3_T.RDS")

la2.T = readRDS("discrete_model_results/la2_T.RDS")
sf2.T = readRDS("discrete_model_results/sf2_T.RDS")
ny3.T = readRDS("discrete_model_results/ny3_T.RDS")

la.static = left_join(la.station, data.frame(id =  names(la2.T$FoundComms), 
                                             block  =  as.factor((la2.T$FoundComms))))
sf.static = left_join(sf.station, data.frame(id =  as.integer(names(sf2.T$FoundComms)), 
                                             block  =  as.factor((sf2.T$FoundComms))), stringsAsFactors = F)
ny.static = left_join(ny1610.station, data.frame(id =  as.integer(names(ny3.T$FoundComms)), 
                                                 block  =  as.factor((ny3.T$FoundComms)), stringsAsFactors = F))

#      plots ####
#       LA ----

plot_la_static = ggmap(la_background) +
  geom_point(data=la.static, shape = 21, aes(x=lon, y=lat, fill = block), size = 4) +
  ggtitle("Los Angeles, Static Discrete Membership") +
  scale_fill_manual(values = c("black", "white")) +
  #scale_size(guide = "none") +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())
plot_la_static

ggsave(filename = "IMG/la_static_discrete.png", plot_la_static)

#       SF ----
plot_sf_static = ggmap(sf_background) +
  geom_point(data=sf.static, shape = 21, aes(x=lon, y=lat, fill = block), size = 4) +
  ggtitle("San Francisco, Static Discrete Membership") +
  scale_fill_manual(values = c("black", "white")) +
  #scale_size(guide = "none") +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())
plot_sf_static

#       NY ----

plot_ny_static = ggmap(ny_background) +
  geom_point(data = ny.static, aes(x = lon, y = lat, color = block, fill = block,
                                   shape = block), color = "black", alpha = .7) +
  ggtitle("New York, Static Discrete Membership") +
  xlim(-74.02,-73.93) + ylim(40.65, 40.8) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(size = guide_legend(order = 2)) +
  scale_fill_manual(name = "role", 
                    #labels = c("Brooklyn", "Manhattan (home)", "Manhattan (work)"), 
                    values = c("white", "gray", "black"), guide = guide_legend(order = 1)) + 
  scale_shape_manual(name = "role", 
                     #labels = c("Brooklyn", "Manhattan (home)", "Manhattan (work)"), 
                     values = c(24, 22, 21), guide = guide_legend(order = 1)) +
  scale_size(range = c(1,4)) +
  #xlab("longitude") + ylab("latitude")
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())
plot_ny_static

#     Final plots ----

ggsave(filename = "IMG/la_sf_static_discrete.png", plot_grid(plot_la_static, plot_sf_static, ncol = 2, align = "h"), width = 9, height = 4)

ggsave(filename = "IMG/ny_static_discrete.png", plot_ny_static)
################################################################################################
# 7. NEW YORK (MANHATTAN) SUBSET  ----
#   Subset -- choose the smallest one for convenience ----

table(ny_discrete.3$role)

home_manhattan_id = filter(ny_discrete.3, role == "Manhattan (home)")$id #make sure to do the pre-processing first
ny1610_hm = filter(ny1610, Start.Station.ID %in% home_manhattan_id & End.Station.ID %in% home_manhattan_id) #324292
ny1610_hm %>% group_by(End.Station.ID) %>% summarize(n()) %>% select("n()")
n_distinct(ny1610_hm$Start.Station.ID); n_distinct(ny1610_hm$End.Station.ID) # check equal
to.remove = unique(ny1610_hm$End.Station.ID)[which(!(unique(ny1610_hm$End.Station.ID) %in% unique(ny1610_hm$Start.Station.ID)))]
ny1610_hm = filter(ny1610_hm, !End.Station.ID %in% to.remove)
# don't include degrees because they reflect whole-system degrees
ny1610_hm.station = filter(ny1610.station %>% select(id, lat, lon), id %in% unique(ny1610_hm$Start.Station.ID))

write.csv(ny1610_hm,"data/cleaned/ny1610_hm_no_weekend.csv", row.names = F)

#   Make time-sliced ---
ny_byhour_hm = list()
for (i in 0:23) {
  ny_hour_hm = subset(ny1610_hm, hour == i)
  ny_hour_hm = aggregate(ny_hour_hm[,c("count")], #"duration"
                         by = data.frame(ny_hour_hm[,c("Start.Station.ID", "End.Station.ID")]), 
                         FUN = sum)
  ny_byhour_hm[[i+1]] = ny_hour_hm
}

# saveRDS(object = ny_byhour_hm, "data/ny_byhour_hm.RDS")
#----------------------------------------------------------------------------------------------
#   Run subset models ----

# Time-dependent - discrete

  for (i in 2:6) {
    ny_hm = sbmt(ny_byhour_hm,  degreeCorrect = 3, directed = T, klPerNetwork = 50, maxComms = i, seed = 1)
    saveRDS(ny_hm, paste0("discrete_model_results/ny_hm", i, "_3.RDS"))
  }


ny_hm.tmp = data.frame(ny1610_hm %>% group_by(Start.Station.ID, End.Station.ID) %>% summarize(count = sum(count)), stringsAsFactors = F)

# Static - discrete

  for (i in 2:6) {
    ny_hm_T = sbmt(list(ny_hm.tmp), maxComms = i, degreeCorrect = 3, directed = T, klPerNetwork = 50, seed = 1)
    saveRDS(ny_hm_T, paste0("discrete_model_results/ny_hm", i, "_T.RDS"))
  }

# Time-dependent - Continuous

# setwd("mixed_model_implementation")
# system("python LosAngeles.py | tee LAoutput.txt")
# setwd("..")

# 
#----------------------------------------------------------------------------------------------
#   Load model results ----

#     Discrete ----
ny_hm2.3 = readRDS("discrete_model_results/ny_hm2_3.RDS")
ny_hm3.3 = readRDS("discrete_model_results/ny_hm3_3.RDS")
ny_hm4.3 = readRDS("discrete_model_results/ny_hm4_3.RDS")
ny_hm5.3 = readRDS("discrete_model_results/ny_hm5_3.RDS")
ny_hm6.3 = readRDS("discrete_model_results/ny_hm6_3.RDS")


ny_hm_discrete.2 = left_join(ny1610_hm.station, data.frame(id =  as.integer(names(ny_hm2.3$FoundComms)),
                                                           role  =  ny_hm2.3$FoundComms, stringsAsFactors = F))
ny_hm_discrete.3 = left_join(ny1610_hm.station, data.frame(id =  as.integer(names(ny_hm3.3$FoundComms)),
                                                           role  =  ny_hm3.3$FoundComms, stringsAsFactors = F))
ny_hm_discrete.4 = left_join(ny1610_hm.station, data.frame(id =  as.integer(names(ny_hm4.3$FoundComms)),
                                                           role  =  ny_hm4.3$FoundComms, stringsAsFactors = F))
ny_hm_discrete.5 = left_join(ny1610_hm.station, data.frame(id =  as.integer(names(ny_hm5.3$FoundComms)),
                                                           role  =  ny_hm5.3$FoundComms , stringsAsFactors = F))
ny_hm_discrete.6 = left_join(ny1610_hm.station, data.frame(id =  as.integer(names(ny_hm6.3$FoundComms)),
                                                           role  =  ny_hm6.3$FoundComms, stringsAsFactors = F))

## assign labels

left_join(ny1610.station.names, ny_hm_discrete.2, by = c("ID"="id")) %>% filter(!is.na(lat)) %>% select(ID, name, role)
ny_hm_discrete.2$role = as.factor(ny_hm_discrete.2$role)
levels(ny_hm_discrete.2$role) = c("work", "home")

# 3 blocks

left_join(ny1610.station.names, ny_hm_discrete.3, by = c("ID"="id")) %>% filter(!is.na(lat)) %>% select(ID, name, role)
ny_hm_discrete.3$role = as.factor(ny_hm_discrete.3$role)
levels(ny_hm_discrete.3$role) = c("work", "home", "mixed")

# 4 blocks

left_join(ny1610.station.names, ny_hm_discrete.4, by = c("ID"="id")) %>% filter(!is.na(lat)) %>% select(ID, name, role)
ny_hm_discrete.4$role = as.factor(ny_hm_discrete.4$role)
levels(ny_hm_discrete.4$role) = c("work", "park", "home", "mixed")

# 5 blocks

left_join(ny1610.station.names, ny_hm_discrete.5, by = c("ID"="id")) %>% filter(!is.na(lat)) %>% select(ID, name, role)
ny_hm_discrete.5$role = as.factor(ny_hm_discrete.5$role)
levels(ny_hm_discrete.5$role) = c("work", "park", "mixed", "home (east)", "home (west)")

#     Static discrete ----
ny_hm2.T = readRDS("discrete_model_results/ny_hm2_T.RDS")
ny_hm3.T = readRDS("discrete_model_results/ny_hm3_T.RDS")
ny_hm4.T = readRDS("discrete_model_results/ny_hm4_T.RDS")
ny_hm5.T = readRDS("discrete_model_results/ny_hm5_T.RDS")
ny_hm6.T = readRDS("discrete_model_results/ny_hm6_T.RDS")

ny_hm_static.2 = left_join(ny1610_hm.station, data.frame(id =  as.integer(names(ny_hm2.T$FoundComms)),
                                                           role  =  (ny_hm2.T$FoundComms), stringsAsFactors = F))
ny_hm_static.3 = left_join(ny1610_hm.station, data.frame(id =  as.integer(names(ny_hm3.T$F)),
                                                         role  =  (ny_hm3.T$FoundComms), stringsAsFactors = F))
ny_hm_static.4 = left_join(ny1610_hm.station, data.frame(id =  as.integer(names(ny_hm4.T$FoundComms)),
                                                         role  =  (ny_hm4.T$FoundComms), stringsAsFactors = F))
ny_hm_static.5 = left_join(ny1610_hm.station, data.frame(id =  as.integer(names(ny_hm5.T$FoundComms)),
                                                         role  =  (ny_hm5.T$FoundComms), stringsAsFactors = F))
ny_hm_static.6 = left_join(ny1610_hm.station, data.frame(id =  as.integer(names(ny_hm6.T$FoundComms)),
                                                         role  =  (ny_hm6.T$FoundComms), stringsAsFactors = F))


#     Continuous (aka Mixed) ----


# 2 blocks 
ny_hm_continuous.roles.2 = read.csv("mixed_model_results/Manhattan_2_roles.csv", sep = ",", header = T)
ny_hm_continuous.2 = left_join(ny1610_hm.station, ny_hm_continuous.roles.2, 
                               by = c("id"="X"))
ny_hm_continuous.omega.2 = lapply(read.csv("mixed_model_results/Manhattan_2_omega.csv", sep = ",", header = T), matrix, 2, 2)

# pre-process 
ny_hm_continuous.2$C_total = rowSums(ny_hm_continuous.2[,c("X0","X1")])
ny_hm_continuous.2$role = ny_hm_continuous.2[,c("X0")]/ny_hm_continuous.2$C_total
ny_hm_continuous.2$discretize = apply(ny_hm_continuous.2[,c("X0", "X1")]/ny_hm_continuous.2$C_total, 1, which.max)
# assign labels
left_join(ny1610.station.names, ny_hm_continuous.2, by = c("ID"="id")) %>% filter(!is.na(lat)) %>% select(ID, name, discretize)
ny_hm_continuous_labels.2 = c("work", "home")



# 3 blocks 
ny_hm_continuous.roles.3 = read.csv("mixed_model_results/Manhattan_3_roles.csv", sep = ",", header = T)
ny_hm_continuous.3 = left_join(ny1610_hm.station, ny_hm_continuous.roles.3, 
                               by = c("id"="X"))[,c("id", "lon", "lat", "X0", "X1", "X2")]
ny_hm_continuous.omega.3 = lapply(read.csv("mixed_model_results/Manhattan_3_omega.csv", sep = ",", header = T), matrix, 3, 3)
# pre-process 
ny_hm_continuous.3$C_total = rowSums(ny_hm_continuous.3[,c("X0","X1","X2")])
ny_hm_continuous.3$role = rgb(ny_hm_continuous.3[,c("X0","X1","X2")]/ny_hm_continuous.3$C_total)
ny_hm_continuous.3$discretize = apply(ny_hm_continuous.3[,c("X2", "X1", "X0")]/ny_hm_continuous.3$C_total, 1, which.max)
# assign labels
left_join(ny1610.station.names, ny_hm_continuous.3, by = c("ID"="id")) %>% filter(!is.na(lat)) %>% select(ID, name, discretize)
ny_hm_continuous_labels.3 = c("home", "work", "mixed")

# 4 blocks 
ny_hm_continuous.roles.4 = read.csv("mixed_model_results/Manhattan_4_roles.csv", sep = ",", header = T)
ny_hm_continuous.4 = left_join(ny1610_hm.station, ny_hm_continuous.roles.4, 
                               by = c("id"="X"))[,c("id", "lon", "lat", "X0", "X1", "X2", "X3")]
ny_hm_continuous.omega.4 = lapply(read.csv("mixed_model_results/Manhattan_4_omega.csv", sep = ",", header = T), matrix, 4, 4)

# pre-process 
ny_hm_continuous.4$C_total = rowSums(ny_hm_continuous.4[,c("X0","X1","X2","X3")])
ny_hm_continuous.4$role = rgb(ny_hm_continuous.4[,c("X0","X3","X2","X1")]/ny_hm_continuous.4$C_total)
ny_hm_continuous.4$discretize = apply(ny_hm_continuous.4[,c("X0", "X1", "X2", "X3")]/ny_hm_continuous.4$C_total, 1, which.max)
# assign labels
left_join(ny1610.station.names, ny_hm_continuous.4, by = c("ID"="id")) %>% filter(!is.na(lat)) %>% select(ID, name, discretize)
ny_hm_continuous_labels.4 = c("home (west)", "work", "home (east)", "mixed")


# 5 blocks 
ny_hm_continuous.roles.5 = read.csv("mixed_model_results/Manhattan_5_roles.csv", sep = ",", header = T)
ny_hm_continuous.5 = left_join(ny1610_hm.station, ny_hm_continuous.roles.5, 
                               by = c("id"="X"))[,c("id", "lon", "lat", "X0", "X1", "X2", "X3", "X4")]
ny_hm_continuous.omega.5 = lapply(read.csv("mixed_model_results/Manhattan_5_omega.csv", sep = ",", header = T), matrix, 5, 5)

# pre-process 
ny_hm_continuous.5$C_total = rowSums(ny_hm_continuous.5[,c("X0","X1","X2","X3","X4")])

# enough roles that we can't use rgb anymore
role.colors.5 = matrix(
  c( 215,25,28,
     253,174,97,
     255,255,191,
     171,221,164,
     43,131,186), nrow = 5, ncol = 3, byrow = T) # <- selected using divergent, colorsafe scheme: http://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=3. Hex: role.colors.5 = c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba") 
ny_hm_continuous.5$role  = rgb(t(apply(ny_hm_continuous.5[,c("X0", "X1", "X2", "X3","X4")],
                                 1, function(x) {colSums(diag(x/sum(x)) %*% role.colors.5)})), maxColorValue = 255)
ny_hm_continuous.5$discretize = apply(ny_hm_continuous.5[,c("X0", "X1", "X2", "X3","X4")]/ny_hm_continuous.5$C_total, 1, which.max)
# assign labels
left_join(ny1610.station.names, ny_hm_continuous.5, by = c("ID"="id")) %>% filter(!is.na(lat)) %>% select(ID, name, discretize)
ny_hm_continuous_labels.5 = c("mixed", "home (west)", "park", "home (east)", "work")


# 6 blocks 
ny_hm_continuous.roles.6 = read.csv("mixed_model_results/Manhattan_6_roles.csv", sep = ",", header = T)
ny_hm_continuous.6 = left_join(ny1610_hm.station, ny_hm_continuous.roles.6, 
                               by = c("id"="X"))[,c("id", "lon", "lat", "X0", "X1", "X2", "X3", "X4", "X5")]
ny_hm_continuous.omega.6 = lapply(read.csv("mixed_model_results/Manhattan_6_omega.csv", sep = ",", header = T), matrix, 6, 6)

# pre-process 
ny_hm_continuous.6$C_total = rowSums(ny_hm_continuous.6[,c("X0","X1","X2","X3","X4","X5")])
ny_hm_continuous.6$role = rgb(ny_hm_continuous.6[,c("X0","X3","X2","X1","X4","X5")]/ny_hm_continuous.6$C_total)
ny_hm_continuous.6$discretize = apply(ny_hm_continuous.6[,c("X0", "X1", "X2", "X3","X4","X5")]/ny_hm_continuous.6$C_total, 1, which.max)
# assign labels
left_join(ny1610.station.names, ny_hm_continuous.6, by = c("ID"="id")) %>% filter(!is.na(lat)) %>% select(ID, name, discretize)
ny_hm_continuous_labels.6 = c("park", "home (east)", "mixed (east)", "home (west)", "work", "mixed (west)")


#     likelihood ####
    
  A = edgelist_to_adj(ny_byhour_hm)
  N = nrow(ny1610_hm.station)
  Time = 24
  
  ny_hm_lk_table = data.frame(
    
    "Blocks" = c(2,3,4,5,6),
    
    "Mixed # par." = c(tdmm_n_param(N, 2, Time),  tdmm_n_param(N, 3, Time), tdmm_n_param(N, 4, Time),
                     tdmm_n_param(N, 5, Time), tdmm_n_param(N, 6, Time)),
    
    "Mixed log-lik" = 
      round(  c(tdmm_sbm_llik(A, data.frame(ny_hm_continuous.roles.2[,-1], 
                row.names = ny_hm_continuous.roles.2[,1]), ny_hm_continuous.omega.2),
      tdmm_sbm_llik(A, data.frame(ny_hm_continuous.roles.3[,-1], row.names = ny_hm_continuous.roles.3[,1]), 
                    ny_hm_continuous.omega.3),
      tdmm_sbm_llik(A, data.frame(ny_hm_continuous.roles.4[,-1], row.names = ny_hm_continuous.roles.4[,1]), 
                    ny_hm_continuous.omega.4),
      tdmm_sbm_llik(A, data.frame(ny_hm_continuous.roles.5[,-1], row.names = ny_hm_continuous.roles.5[,1]), 
                    ny_hm_continuous.omega.5),
      tdmm_sbm_llik(A, data.frame(ny_hm_continuous.roles.6[,-1], row.names = ny_hm_continuous.roles.6[,1]), 
                    ny_hm_continuous.omega.6))),
  
    "Discrete par." = c(tdd_n_param(N, 2, Time),  tdd_n_param(N, 3, Time), tdd_n_param(N, 4, Time),
                        tdd_n_param(N, 5, Time), tdd_n_param(N, 6, Time)),
    
    "Discrete log-lik" = round(c(tdd_sbm_llik(A, ny_hm2.3$FoundComms, ny_hm2.3$EdgeMatrix),
                                 tdd_sbm_llik(A, ny_hm3.3$FoundComms, ny_hm3.3$EdgeMatrix),
                                 tdd_sbm_llik(A, ny_hm4.3$FoundComms, ny_hm4.3$EdgeMatrix),
                                 tdd_sbm_llik(A, ny_hm5.3$FoundComms, ny_hm5.3$EdgeMatrix),
                                 tdd_sbm_llik(A, ny_hm6.3$FoundComms, ny_hm6.3$EdgeMatrix)
                               ))
    )
    
  ny_hm_lk_table  

  xtable::xtable(ny_hm_lk_table, caption = "Model likelihood comparisons",
                 label = "ny_hm_llik_table", digits = 0)

#Checked they agree with tdmm-sbm llik
#Best log-likelihoods[-260625.18052058 -235162.42763762 -212295.31117002 -198489.91494118, -189670.85830275]
  
#----------------------------------------------------------------------------------------------
#   Plots ----

#     position plot ----
#       2 blocks ----

ny_hm_compare.2 = cbind(rbind(ny_hm_discrete.2[,c("id", "lon", "lat")],
                              ny_hm_static.2[,c("id", "lon", "lat")], 
                              ny_hm_continuous.2[,c("id", "lon", "lat")]),
                        
                        model = c(rep("TDD-SBM (2)", nrow(ny_hm_discrete.2)), 
                                  rep("SBM (2)", nrow(ny_hm_static.2)), 
                                  rep("TDMM-SBM (2)", nrow(ny_hm_continuous.2))),
                        
                        role_tdmm = c(gray(c(0:1)[as.numeric(ny_hm_discrete.2$role)]),
                                      gray(c(0:1)[ny_hm_static.2$role+1]), 
                                      gray(1-ny_hm_continuous.2$role))) %>% 
    
  mutate(role_tdd = as.factor(c(as.numeric(ny_hm_discrete.2$role), 
                                c(1,1)[ny_hm_static.2$role+1],
                                rep(1, nrow(ny_hm_continuous.2)))))

ny_hm_compare.2$model = factor(ny_hm_compare.2$model, levels = levels(ny_hm_compare.2$model)[c(3,2,1)])  #reorder to put sbm on right

ny_plot_hm.2 = ggplot(ny_hm_compare.2) +
  geom_point(aes(x = lon, y = lat, shape = role_tdd, fill = I(as.character(ny_hm_compare.2$role_tdmm))), size = 4) + 
  facet_grid(~model) +
  #scale_fill_continuous(low = "black", high = "white", breaks = seq(0,1,by=0.5), 
   #                     labels = c("work", "", "home")) + 
  scale_shape_manual(values = c(21, 22), labels = c("work", "home"), name = "TD role") +
  guides(shape = guide_legend(override.aes = list(fill = c("black","white")))) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())
ny_plot_hm.2

ggsave("IMG/ny_hm_2.png", ny_plot_hm.2, width = 12, height = 6)

#       3 blocks ----

ny_hm_static.3$role  =  (0:2)[ny_hm_static.3$role +1] #re-order static to match tdd coloring

ny_hm_compare.3 = cbind(rbind(ny_hm_discrete.3[,c("id", "lon", "lat")],
                              ny_hm_static.3[,c("id", "lon", "lat")], 
                              ny_hm_continuous.3[,c("id", "lon", "lat")]),
                        model = c(rep("TDD-SBM (3)", nrow(ny_hm_discrete.3)), 
                                  rep("SBM (3)", nrow(ny_hm_static.3)), 
                                  rep("TDMM-SBM (3)", nrow(ny_hm_continuous.3))),
                        role_tdmm = c(gray((as.numeric(ny_hm_discrete.3$role)-1)/(max(as.numeric(ny_hm_discrete.3$role))-1)),
                                      gray(ny_hm_static.3$role/max(ny_hm_static.3$role)), 
                                      ny_hm_continuous.3$role)) %>% 
  mutate(role_tdd = as.factor(c(as.numeric(ny_hm_discrete.3$role), 
                                c(1,1,1)[ny_hm_static.3$role+1],
                                rep(1, nrow(ny_hm_continuous.3)))))

#order facets
ny_hm_compare.3$model = factor(ny_hm_compare.3$model, levels = c("TDMM-SBM (3)", "TDD-SBM (3)", "SBM (3)")) 

ny_plot_hm.3 = ggplot(ny_hm_compare.3) +
  geom_point(aes(x = lon, y = lat, shape = role_tdd, fill = I(as.character(ny_hm_compare.3$role_tdmm))), size = 4) + 
  facet_grid(~model) +
  scale_shape_manual(values = c(21, 22, 23), labels = c("work", "home", "mixed"), name = "role (TDD-SBM)") +
  guides(shape = guide_legend(override.aes = list(fill = c("black","gray","white")))) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())
ny_plot_hm.3

#       4 blocks ----

ny_hm_compare.4 = cbind(rbind(ny_hm_discrete.4[,c("id", "lon", "lat")],
                            ny_hm_static.4[,c("id", "lon", "lat")], 
                            ny_hm_continuous.4[,c("id", "lon", "lat")]),
                      model = c(rep("TDD-SBM (4)", nrow(ny_hm_discrete.4)), 
                                rep("SBM (4)", nrow(ny_hm_static.4)), 
                                rep("TDMM-SBM (4)", nrow(ny_hm_continuous.4))),
                      role_tdmm = c(gray(as.numeric(ny_hm_discrete.4$role)/max(as.numeric(ny_hm_discrete.4$role))),
                                    gray(ny_hm_static.4$role/max(ny_hm_static.4$role)), 
                                    ny_hm_continuous.4$role)) %>% 
  mutate(role_tdd = as.factor(c(c(4,3,2,1)[as.numeric(ny_hm_discrete.4$role)], 
                                c(1,1,1,1)[ny_hm_static.4$role+1],
                      rep(1, nrow(ny_hm_continuous.4)))))
                      
ny_hm_compare.4$model = factor(ny_hm_compare.4$model, levels = c("TDMM-SBM (4)", "TDD-SBM (4)" , "SBM (4)"))  #reorder to put sbm on right


ny_plot_hm.4 = ggplot(ny_hm_compare.4) +
  geom_point(aes(x = lon, y = lat, shape = role_tdd), fill = as.character(ny_hm_compare.4$role_tdmm), size = 4) + 
  facet_grid(~model) +
  scale_shape_manual(values = c(21, 22, 23, 24), labels = c("home", "mixed", "home|park", "work")) +
  guides(shape = guide_legend(title = "TD role", override.aes = list(fill = c("blue", "green", "red", "black")))) +
  #geom_text(aes(label = id), col = "red") +
  #geom_point(x = -73.9935, y = 40.7506, color = "green", size = 3, shape = 13) + #penn station
  #geom_point(x = -73.9903, y = 40.7569, color = "green", size = 3, shape = 13) + #port authority
  #scale_fill_manual(values = c("white","black","gray")) + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())
ny_plot_hm.4

ggsave("IMG/ny_hm_4.png", ny_plot_hm.4, width = 12, height = 6)

#       5 blocks ----

ny_hm_compare.5 = cbind(rbind(ny_hm_discrete.5[,c("id", "lon", "lat")],
                              ny_hm_static.5[,c("id", "lon", "lat")],
                              ny_hm_continuous.5[,c("id", "lon", "lat")]),
                        model = c(rep("TDD-SBM (5)", nrow(ny_hm_discrete.5)),
                                  rep("SBM (5)", nrow(ny_hm_static.5)),
                                  rep("TDMM-SBM (5)", nrow(ny_hm_continuous.5))),
                        role_tdmm = c(#gray(as.numeric(ny_hm_discrete.5$role)/max(as.numeric(ny_hm_discrete.5$role))),
                          (rgb(role.colors.5, maxColorValue = 255)[c(5,3,1,4,2)])[as.numeric(ny_hm_discrete.5$role)],
                          gray(ny_hm_static.5$role/max(ny_hm_static.5$role)),
                          ny_hm_continuous.5$role)) %>%
  mutate(role_tdd = as.factor(c(c(1,2,3,4,5)[as.numeric(ny_hm_discrete.5$role)],
                                c(1,1,1,1,1)[ny_hm_static.5$role+1],
                                #c(1,3,2,4,5)[ny_hm_static.5$role+1],
                                rep(1, nrow(ny_hm_continuous.5)))))

ny_hm_compare.5$model = factor(ny_hm_compare.5$model, levels = c("TDMM-SBM (5)", "TDD-SBM (5)" , "SBM (5)"))  #reorder to put sbm on right

ny_plot_hm.5 = ggplot(ny_hm_compare.5) +
  geom_point(aes(x = lon, y = lat, shape = role_tdd), 
             fill = as.character(ny_hm_compare.5$role_tdmm), size = 4) + 
  facet_grid(~model) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25), labels = levels(ny_hm_discrete.5$role)) +
  guides(shape = guide_legend(title = "TD role", 
                              override.aes = list(fill = rgb(role.colors.5, maxColorValue = 255)[c(5,3,1,4,2)]))) +
  #geom_text(aes(label = id), col = "red") +
  #geom_point(x = -73.9935, y = 40.7506, color = "green", size = 3, shape = 13) + #penn station
  #geom_point(x = -73.9903, y = 40.7569, color = "green", size = 3, shape = 13) + #port authority
  #scale_fill_manual(values = c("white","black","gray")) + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())

ggsave("IMG/ny_hm_5.png", ny_plot_hm.5, width = 12, height = 6)

#      final plot (5 blocks) ----

ny_plot_hm = ggplot(ny_hm_compare.5 %>% filter(model != "SBM (5)")) +
  geom_point(aes(x = lon, y = lat, shape = role_tdd), 
             fill = as.character(ny_hm_compare.5$role_tdmm[!ny_hm_compare.5$model=="SBM (5)"]), size = 4) + 
  facet_grid(~model) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25), labels = levels(ny_hm_discrete.5$role)) +
  guides(shape = guide_legend(title = "TD role", 
                              override.aes = list(fill = rgb(role.colors.5, maxColorValue = 255)[c(5,3,1,4,2)]))) +
  #geom_text(aes(label = id), col = "red") +
  #geom_point(x = -73.9935, y = 40.7506, color = "green", size = 3, shape = 13) + #penn station
  #geom_point(x = -73.9903, y = 40.7569, color = "green", size = 3, shape = 13) + #port authority
  #scale_fill_manual(values = c("white","black","gray")) + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.line.y = element_blank())
ny_plot_hm

ggsave("IMG/NY_HM_mixed_discrete.png", ny_plot_hm, width = 12, height = 6)

#       check on penn station and port authority - 3230 (Penn), 3236 (Port) found by annotating plot_ny_hm ----
grep(ny1610.station.names$name, pattern = "Penn"); ny1610.station.names[444,] 
grep(ny1610.station.names$name, pattern = "42"); ny1610.station.names[449,]; 

# departure in morning
filter(ny1610, Start.Station.ID == 3236) %>% group_by(hour) %>% 
  summarize( hour1 = first(hour), sum = sum(count)) %>% 
  select(hour1, sum); plot(.Last.value, type = "l")

# arrivals in evening
filter(ny1610, End.Station.ID == 3236) %>% group_by(hour) %>% 
  summarize( hour1 = first(hour), sum = sum(count)) %>% 
  select(hour1, sum); plot(.Last.value, type = "l")

#invesigate lenox hill hospital
filter(ny1610, Start.Station.ID == "3154") %>% group_by(hour) %>% summarize(n()) %>% plot(type = "l")
filter(ny1610, End.Station.ID == "3154") %>% group_by(hour) %>% summarize(n()) %>% plot(type = "l")

#     omega ----

#       discrete ----


#Figure out labels

ny_dis = ny_hm5.3

left_join(data.frame(ny_dis$FoundComms, id = as.numeric(names(ny_dis$FoundComms))), ny1610.station.names, by = c('id' = "ID"))

tmp = ny1610_hm 
tmp$Start.Station.role= ny_dis$FoundComms[as.character(tmp$Start.Station.ID)]
tmp$End.Station.role= ny_dis$FoundComms[as.character(tmp$End.Station.ID)]
tmp %>% group_by(hour, Start.Station.role, End.Station.role) %>% summarize(s = sum(count)) %>% 
  mutate(group = interaction(Start.Station.role, End.Station.role)) %>%
  ggplot() + geom_line(aes(x = hour, y = s, color = group)) +
  facet_wrap(.~group)

# 2 blocks

ny.hm.2.plot.order = c(t(outer(c("home", "work"),
                             c("home", "work"), paste, sep =" to ")))

K = nrow(ny_hm2.3$EdgeMatrix[[1]])
results =  data.frame(t(sapply(ny_hm2.3$EdgeMatrix, unlist)))
names(results) = c(outer(levels(ny_hm_discrete.2$role),
                           levels(ny_hm_discrete.2$role), paste, sep =" to "))
results$hour = 0:(nrow(results)-1)
results = melt(results, id.vars = "hour", value.name = "count")
results$variable  = factor(results$variable, levels = ny.hm.2.plot.order)

# 5 blocks

ny.hm.5.plot.order = c(t(outer(c("home W","home E", "mixed","park", "work"),
                             c("home W","home E", "mixed","park", "work"), paste, sep =" to ")))

K = nrow(ny_hm5.3$EdgeMatrix[[1]])
results =  data.frame(t(sapply(ny_hm5.3$EdgeMatrix, unlist)))
# levels(ny_hm_discrete.5$role)
names(results) = c(outer(c("work","park","mixed", "home E", "home W"), 
                         c("work","park","mixed", "home E", "home W"), paste, sep =" to "))
results$hour = 0:(nrow(results)-1)
results = melt(results, id.vars = "hour", value.name = "count")
results$variable  = factor(results$variable, levels = ny.hm.5.plot.order)

# plot either

ny_hm_omega_discrete = ggplot(results, aes(x = hour, y = count, color = 1)) + geom_line() + 
  facet_wrap(.~variable, scales = "fixed") +
  guides(color = "none") +
  ggtitle("Time-Dependent Block Parameters for Manhattan Subnetwork Discrete Model") +
  theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
ny_hm_omega_discrete

ggsave("IMG/ny_hm_omega_discrete.png", ny_hm_omega_discrete, width = 12, height = 8)

#       continuous ----

# 2 block
results = data.frame(t(sapply(ny_hm_continuous.omega.2, unlist)))
results$hour = 0:(nrow(results)-1)
names(results) = c(as.vector(outer(ny_hm_continuous_labels.2, ny_hm_continuous_labels.2, 
                                   paste, sep = " to ")), "hour")
results = melt(results, id.vars = "hour", value.name = "count")
results$variable  = factor(results$variable, levels = la.plot.order) #can use la's order

ny_hm_omega_continuous.2 = ggplot(results, aes(x = hour, y = count)) +
  geom_line() + facet_wrap(.~variable) +
  guides(color = "none") +
  ggtitle("Time-Dependent Block Parameters for Manhattan Subnetwork, NY Mixed Model") +
  theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5), 
                     strip.text = element_text(size = 12))

ny_hm_omega_continuous.2

# 3 block
results = data.frame(t(sapply(ny_hm_continuous.omega.3, unlist)))
names3 = ny_hm_continuous_labels.3
results$hour = 0:(nrow(results)-1)
names(results) = c(as.vector(outer(names3, names3, paste, sep = " to ")), "hour")
results = melt(results, id.vars = "hour", value.name = "count")

ny_hm_omega_continuous.3 = ggplot(results, aes(x = hour, y = count)) +
  geom_line() + facet_wrap(.~variable, scales = "fixed", dir = "v") +
  guides(color = "none") +
  ggtitle("Time-Dependent Block Parameters for Manhattan Subnetwork, NY Mixed Model") +
  theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5), 
                     strip.text = element_text(size = 12))

ny_hm_omega_continuous.3

# 4 block
results = data.frame(t(sapply(ny_hm_continuous.omega.4, unlist)))
names4 =  ny_hm_continuous_labels.4
results$hour = 0:(nrow(results)-1)
names(results) = c(as.vector(outer(names4, names4, paste, sep = " to ")), "hour")
results = melt(results, id.vars = "hour", value.name = "count")

ny_hm_omega_continuous.4 = ggplot(results, aes(x = hour, y = count)) +
  geom_line(color = "black") + facet_wrap(.~variable, scales = "fixed", dir = "v") +
  guides(color = "none") +
  ggtitle("Time-Dependent Block Parameters for Manhattan Subnetwork Mixed Model") +
  theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5), 
                     strip.text = element_text(size = 12))

ny_hm_omega_continuous.4

# 5 block
results = data.frame(t(sapply(ny_hm_continuous.omega.5, unlist)))
names5 =  c("mixed","home W", "park", "home E", "work") # based on ny_hm_continuous_labels.5
results$hour = 0:(nrow(results)-1)
names(results) = c(as.vector(outer(names5, names5, paste, sep = " to ")), "hour")
results = melt(results, id.vars = "hour", value.name = "count")
results$variable  = factor(results$variable, levels = ny.hm.5.plot.order)

ny_hm_omega_continuous.5 = ggplot(results, aes(x = hour, y = count)) +
  geom_line(color = "black") + facet_wrap(.~variable, scales = "fixed", dir = "h") +
  guides(color = "none") +
  ggtitle("Time-Dependent Block Parameters for Manhattan Subnetwork Mixed Model") +
  theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

ny_hm_omega_continuous.5

#       Final plot ####
ggsave("IMG/ny_hm_omega.png", plot_grid(ny_hm_omega_continuous.5, ny_hm_omega_discrete), width = 14, height = 7)

################################################################################################
# 8. Comparison of role distribution for two-block models (not in paper) ####

ggplot() + geom_density(data = la_continuous, aes(x = role), col = "blue", n  = 1024) + 
           geom_density(data = sf_continuous, aes(x = role), col = "green", n  = 1024) + 
           geom_density(data = ny_continuous.2, aes(x = role), col = "red", n  = 1024) +
           geom_density(data = ny_hm_continuous.2, aes(x = role), col = "black", n  = 1024) 
# notice new york city nodes are  less mixed in general

round(prop.table(table(la_discrete$role)),2)
round(prop.table(table(sf_discrete$role)),2)
round(prop.table(table(ny_discrete.2$role)),2)
round(prop.table(table(ny_hm_discrete.2$role)),2)



