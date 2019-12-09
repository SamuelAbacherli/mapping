# ------------------------------------------------------------------------------
# GDP PER CAPITAL PPS OF EUROPE
# ------------------------------------------------------------------------------

shp.europe <- readOGR(dsn = "/Users/SamuelAbacherli/Downloads/CNTR_RG_60M_2016_3857.shp/",
                      layer = "CNTR_RG_60M_2016_3857")
class(shp.europe)

sf.europe <- st_as_sf(shp.europe)
class(sf.europe)
head(sf.europe)

df.gdp <- read.table(file='/Users/SamuelAbacherli/Downloads/tec00114.tsv', sep = '\t', header = TRUE, fill = TRUE)
View(head(df.gdp))

df.gdp$CNTR_ID <- df.gdp$na_item.ppp_cat.geo.time %>% 
  substr(nchar(as.character(df.gdp$na_item.ppp_cat.geo.time)) - 1, 
         nchar(as.character(df.gdp$na_item.ppp_cat.geo.time)))

df.ec <- read.table(file='/Users/SamuelAbacherli/Downloads/tec00013.tsv', sep = '\t', header = TRUE, fill = TRUE)
head(df.ec)

df.ec$CNTR_ID <- df.ec$na_item.unit.geo.time %>% 
  substr(nchar(as.character(df.ec$na_item.unit.geo.time)) - 1, 
         nchar(as.character(df.ec$na_item.unit.geo.time)))

sf.europe <- sf.europe %>% 
  right_join(df.gdp %>% select(CNTR_ID, GDP = X2018), by = "CNTR_ID", suffix = c("", ".gdp")) %>% 
  right_join(df.ec %>% select(CNTR_ID, EC = X2018), by = "CNTR_ID", suffix = c("", ".ec"))

tmap_mode("view")
tm_basemap(leaflet::providers$Stamen.TonerBackground) +
  tm_shape(sf.europe) +
  tm_borders(col = "white", lwd = 1) +
  tm_polygons(col = "GDP", style = "cont", palette = viridis::magma(20, 
                                                      begin = 0.1, 
                                                      end = 0.9, 
                                                      direction = -1)) +
  tm_text(text = "GDP", col = "white", size = "AREA")


# ------------------------------------------------------------------------------
# AREA OF NORTH CAROLINA OVERLAYED
# ------------------------------------------------------------------------------

demo(nc, ask = FALSE, echo = FALSE)
class(nc)

tmap_mode("view")
tm_basemap(leaflet::providers$Stamen.TonerBackground) +
  tm_shape(nc) +
  tm_fill(col = "AREA", alpha = 0.3, palette = viridis::magma(nrow(nc), 
                                                 begin = 0.2, 
                                                 end = 0.9, 
                                                 direction = -1)) +
  tm_borders(col = "white", lwd = 1)
