README

FILE LIST
1) 201608_status_data.csv – approx. 36 million records of bike and dock availability by minute
2) 201608_station_data.csv – 67 records – station ID, name, latitude, longitude, dockcount, city, installation date
3) 201608_trip_data.csv – approx. 314,000 records of individual trips
4) 201608_weather_data.csv – 1,830 records of daily weather by city

Files contain data from 9/1/15 to 8/31/16. This is the third year of Bay Area Bike Share's operation. The first and second year of data can be downloaded separately from the website.

1) STATUS DATA
FILE = "201608_status_data.csv"
-station_id: station ID number (use "201608_station_data.csv" to find corresponding station information)
-bikes_available: number of available bikes
-docks_available: number of available docks
-time: date and time, PST

Note: On 6/30/16, Service in Redwood City was discontinued due to low usage. This included 7 stations: 21, 22, 23, 24, 25, and 26.

Four of these stations have since been moved to either San Francisco or San Jose. (Stations 23, 24, 25 and 26 have become stations 88, 89, 90 and 91 respectively). Although these stations were promptly re-named, there was a delay in assigning them new station IDs. Full details:

On 7/5/16, Station 23, "San Mateo County Center," was renamed to be "5th S. at E. San Salvador St.” On 8/24/16, the station was reassigned to Station 88.
On 7/5/16, Station 24, "Redwood City Public Library," was renamed to be "S. Market St at Park Ave.” On 8/24/16, the station was reassigned to Station 89.
On 8/4/16, Station 25, "Stanford in Redwood City," was renamed to be "Cyril Magnin St at Ellis St.” On 8/24/16, the station was reassigned to Station 91.
On 8/4/16, Station 26, "Kaiser Hospital," was renamed to be "5th St at Folsom St.” On 8/24/16, the station was reassigned to Station 90.

2) STATION INFORMATION
FILE = "201608_station_data.csv"
-station_id: station ID number (corresponds to "station_id" in "201608_status_data.csv")
-name: name of station
-lat: latitude
-long: longitude
-dockcount: number of total docks at station
-landmark: city (San Francisco, Redwood City, Palo Alto, Mountain View, San Jose)
-installation: original date that station was installed. If station was moved, it is noted below.

Note 1: Station names and locations listed on "201608_station_data.csv" represent data that was collected on 8/31/16. However, please note that during 9/1/15 and 8/31/16, the following moves and name changes took place:

Station 21: On 9/16/15, this station was renamed from "Franklin at Maple" to "Sequoia Hospital" and moved to (37.479303,-122.253755)
Station 26: On 9/16/15, this station was renamed from "Redwood City Medical Center" to "Kaiser Hospital" and moved to (37.489704,-122.224728)
Station 30: On 9/28/15, this station was renamed from "Evelyn Park and Ride" to "Middlefield Light Rail Station" and moved to (37.395337,-122.052476)
Station 33: On 9/16/15, this station was renamed from "Rengstorff Avenue / California Street" to "Charleston Park/ North Bayshore Area" and moved to (37.420909,-122.080623)
Station 73: Moved twice. From 3/14/16 – 5/19/16, this station was located at (37.797746, -122.407073). From 5/19/16 to 8/31/16, the station was located at (37.7979, -122.405942). The station name stayed the same for all moves. 
Station 83: On 9/16/15, this station was renamed from "Mezes Park" to "Mezes" and moved to (37.491405,-122.233051)

Note 2: On 6/30/16, Service in Redwood City was discontinued due to low usage. This included 7 stations: 21, 22, 23, 24, 25, and 26.

Four of these stations have since been moved to either San Francisco or San Jose. (Stations 23, 24, 25 and 26 have become stations 88, 89, 90 and 91 respectively). Although these stations were promptly re-named, there was a delay in assigning them new station IDs. Full details:

On 7/5/16, Station 23, "San Mateo County Center," was renamed to be "5th S. at E. San Salvador St.” On 8/24/16, the station was reassigned to Station 88.
On 7/5/16, Station 24, "Redwood City Public Library," was renamed to be "S. Market St at Park Ave.” On 8/24/16, the station was reassigned to Station 89.
On 8/4/16, Station 25, "Stanford in Redwood City," was renamed to be "Cyril Magnin St at Ellis St.” On 8/24/16, the station was reassigned to Station 91.
On 8/4/16, Station 26, "Kaiser Hospital," was renamed to be "5th St at Folsom St.” On 8/24/16, the station was reassigned to Station 90.

3) TRIP DATA
FILE = "201608_trip_data.csv"
-Trip ID: numeric ID of bike trip
-Duration: time of trip in seconds  (trips <1 min and >24 hours are excluded)
-Start Date: start date of trip with date and time, in PST
-Start Station: station name of start station
-Start Terminal: numeric reference for start station
-End Date: end date of trip with date and time, in PST
-End Station: station name for end station
-End Terminal: numeric reference for end station
-Bike #: ID of bike used
-Subscription Type: Subscriber = annual or 30-day member; Customer = 24-hour or 3-day member
-Zip Code: Home zip code of subscriber (customers can choose to manually enter zip at kiosk however data is unreliable) 


4) WEATHER DATA
FILE = "201608_weather_data.csv"
Daily weather information per service area, provided from Weather Underground in PST. Weather is listed from north to south (San Francisco, Redwood City, Palo Alto, Mountain View, San Jose).
        
-Precipitation_In         "numeric, in form x.xx but alpha ""T""= trace when amount less than .01 inch"        
-Cloud_Cover         "scale of 0-8, 0=clear"        
-Zip: 94107=San Francisco, 94063=Redwood City, 94301=Palo Alto, 94041=Mountain View, 95113= San Jose"