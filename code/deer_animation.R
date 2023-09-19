library(moveVis)
library(move)
library(raster)

# The dat_move object is created in the deer_analysis script.
# create time filter to include only one month.
tfilter <- timestamps(dat_move)>as.POSIXct("2023-4-30") & timestamps(dat_move)<as.POSIXct("2023-6-01") 
# align move_data to a uniform time scale, here every 2 hours.
m <- align_move(dat_move[tfilter], res = 2, unit = "hours", spaceMethod = "euclidean")
# m <- spTransform(m, crs("+init=epsg:3857"))

# create spatial frames
frames <- frames_spatial(m, path_colours = hcl.colors(5, "Dark 3")) %>% 
  add_labels(x = "Longitude", y = "Latitude") %>% # add some customizations, such as axis labels
  add_northarrow() %>% 
  add_scalebar() %>% 
  add_timestamps(type = "label") %>% 
  add_progress()

# animate frames
animate_frames(frames, out_file = "moveVis.gif")