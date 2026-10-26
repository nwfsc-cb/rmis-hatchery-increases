library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridisLite)
library(ggspatial)

sf_use_s2(FALSE)

# helper function to make sure polygons are actually closed up
close_poly <- function(m) {
  if (!all(m[1, ] == m[nrow(m), ])) {
    m <- rbind(m, m[1, ])
  }
  st_polygon(list(m))
}

# Get base land geometry
land <- ne_states(country = c("United States of America", "Canada"), returnclass = "sf")

target_crs <- 3005 # BC Albers, since BC is central
land_proj <- st_transform(land, target_crs) |>
  st_make_valid()
land_union <- st_union(land_proj) |>
  st_make_valid()

coastlines <- st_cast(land_union, "MULTILINESTRING")

# Make SRKW range from Ford et al. 2006
srkw_north_coords <- matrix(c(
  -133.5, 54.40,
  -132.5, 54.55,
  -131.6, 54.40,
  -131.1, 53.60,
  -130.8, 52.80,
  -130.3, 52.10,
  -129.5, 51.65,
  -128.6, 51.25,
  -127.6, 50.95,
  -126.8, 50.60,
  -125.2, 50.10,
  -124.2, 49.50,
  -123.3, 48.90,
  -122.5, 48.20,
  -122.2, 47.60,
  -122.6, 47.10,
  -123.2, 47.15,
  -123.8, 48.15,
  -124.7, 48.35,
  -125.4, 48.35,
  -126.5, 49.00,
  -127.8, 49.80,
  -128.8, 50.60,
  -129.8, 51.10,
  -130.6, 51.45,
  -131.4, 51.80,
  -132.3, 52.50,
  -133.1, 53.30,
  -133.8, 54.10,
  -133.5, 54.40
), ncol = 2, byrow = TRUE)

srkw_north <- close_poly(srkw_north_coords) |>
  st_sfc(crs = 4326) |>
  st_transform(target_crs) |>
  st_make_valid()

us_coast_box <- st_bbox(c(xmin = -126.0, ymin = 36.3, xmax = -121.5, ymax = 48.45), crs = 4326) |>
  st_as_sfc() |>
  st_transform(target_crs)

us_outer_coast <- st_intersection(coastlines, us_coast_box)
us_coast_strip <- st_buffer(us_outer_coast, dist = 38000, endCapStyle = "FLAT")

srkw_combined <- st_union(srkw_north, us_coast_strip) |>
  st_make_valid()

srkw_sf <- st_sf(
  layer = factor("SRKW Range"),
  geometry = srkw_combined
)

# Make fishery management zones (50 nm Offshore)
# Exact shared maritime lines:
lat_ak_bc <- 54.67 # Dixon Entrance / A-B Line
lat_wa_or <- 46.25 # Columbia River mouth
lat_or_ca <- 42.00 # 42nd parallel

# Shared international boundary through Juan de Fuca (WA / BC)
jdf_line <- matrix(c(
  -120.00, 49.00,
  -123.00, 49.00,
  -123.25, 48.65,
  -123.50, 48.30,
  -124.75, 48.50,
  -126.15, 48.15 # 50 nm offshore terminus
), ncol = 2, byrow = TRUE)

# Shared AK / BC boundary line across Dixon Entrance
ak_bc_line <- matrix(c(
  -129.50, lat_ak_bc,
  -131.00, lat_ak_bc,
  -133.50, lat_ak_bc,
  -135.00, lat_ak_bc # 50 nm offshore terminus
), ncol = 2, byrow = TRUE)

# A. ALASKA (50 nm offshore SEAK down to Dixon Entrance)
ak_coords <- rbind(
  matrix(c(
    -129.50, 56.50,
    -137.00, 56.50,
    -136.20, 55.50
  ), ncol = 2, byrow = TRUE),
  ak_bc_line[nrow(ak_bc_line):1, ],
  matrix(c(-129.50, 56.50), ncol = 2, byrow = TRUE)
)

# B. BRITISH COLUMBIA (Traces true 50-nm shelf: Haida Gwaii, QC Sound waist, VI)
bc_coords <- rbind(
  ak_bc_line,
  matrix(c(
    -135.00, 54.10, # 50 nm off west Haida Gwaii
    -134.20, 53.20,
    -132.80, 52.00, # 50 nm off Cape St. James
    -130.60, 51.10, # Follows QC Sound natural shelf waist
    -129.60, 50.40, # 50 nm off Brooks Peninsula
    -127.40, 49.00, # 50 nm off SW Vancouver Island
    -126.15, 48.15 # Welds directly to Juan de Fuca offshore terminus
  ), ncol = 2, byrow = TRUE),
  jdf_line[nrow(jdf_line):1, ],
  matrix(c(
    -120.00, 54.67,
    -129.50, 54.67
  ), ncol = 2, byrow = TRUE)
)

# C. WASHINGTON (50 nm offshore Cape Flattery down to Columbia River)
wa_coords <- rbind(
  jdf_line,
  matrix(c(
    -125.80, 47.30, # 50 nm outer shelf
    -125.60, lat_wa_or, # 50 nm offshore Columbia River
    -124.06, lat_wa_or, # Coast
    -120.00, lat_wa_or,
    -120.00, 49.00
  ), ncol = 2, byrow = TRUE)
)

# D. OREGON (50 nm offshore Columbia River to 42°N)
or_coords <- matrix(c(
  -125.60, lat_wa_or,
  -124.06, lat_wa_or,
  -120.00, lat_wa_or,
  -120.00, lat_or_ca,
  -124.21, lat_or_ca,
  -125.70, lat_or_ca,
  -125.80, 44.50,
  -125.60, lat_wa_or
), ncol = 2, byrow = TRUE)

# E. CALIFORNIA (50 nm offshore 42°N to 34°N)
ca_coords <- matrix(c(
  -125.70, lat_or_ca,
  -124.21, lat_or_ca,
  -118.00, lat_or_ca,
  -118.00, 34.00,
  -120.50, 34.00,
  -122.00, 34.00,
  -123.15, 36.35, # 50 nm off Point Sur
  -124.10, 38.00, # 50 nm off Point Reyes
  -125.60, 40.30, # 50 nm off Cape Mendocino
  -125.70, lat_or_ca
), ncol = 2, byrow = TRUE)

fishery_zones <- st_sf(
  region = factor(c("Alaska", "British Columbia", "Washington", "Oregon", "California"),
    levels = c("Alaska", "British Columbia", "Washington", "Oregon", "California")
  ),
  geometry = st_sfc(list(
    close_poly(ak_coords),
    close_poly(bc_coords),
    close_poly(wa_coords),
    close_poly(or_coords),
    close_poly(ca_coords)
  ), crs = 4326)
) |>
  st_transform(target_crs) |>
  st_make_valid()

# Clean boundary divider lines that end at the 50 nm mark
dividers <- st_sf(geometry = st_sfc(list(
  st_linestring(matrix(c(-135.00, lat_ak_bc, -130.00, lat_ak_bc), ncol = 2, byrow = TRUE)),
  st_linestring(matrix(c(-126.15, 48.15, -124.75, 48.50), ncol = 2, byrow = TRUE)),
  st_linestring(matrix(c(-125.60, lat_wa_or, -124.06, lat_wa_or), ncol = 2, byrow = TRUE)),
  st_linestring(matrix(c(-125.70, lat_or_ca, -124.21, lat_or_ca), ncol = 2, byrow = TRUE))
), crs = 4326)) |> st_transform(target_crs)

# Transform for ggplot
srkw_plot <- st_transform(srkw_sf, 4326)
fishery_plot <- st_transform(fishery_zones, 4326)
land_plot <- st_transform(land_proj, 4326)
dividers_plot <- st_transform(dividers, 4326)

mako_palette <- mako(n = 7, begin = 0.22, end = 0.90)[1:5]
names(mako_palette) <- levels(fishery_zones$region)

# Make map
ggplot() +
  # Continuous 50-nm Fishery Zones
  geom_sf(
    data = fishery_plot,
    aes(fill = region),
    color = NA,
    alpha = 0.50
  ) +
  scale_fill_manual(
    name = "Troll Fishery Region",
    values = mako_palette
  ) +

  # Shared jurisdictional divider lines
  geom_sf(
    data = dividers_plot,
    color = "grey20",
    linewidth = 0.45
  ) +

  # Solid SRKW Range Polygon
  geom_sf(
    data = srkw_plot,
    aes(color = layer),
    fill = "#e63946",
    linewidth = 0.55,
    alpha = 0.50
  ) +
  scale_color_manual(
    name = NULL,
    values = c("SRKW Range" = "#6a040f"),
    guide = guide_legend(
      override.aes = list(
        fill = alpha("#e63946", 0.50),
        linewidth = 0.55
      )
    )
  ) +

  # Solid grey land on top (masks inland areas cleanly)
  geom_sf(
    data = land_plot,
    fill = "#bdbdbd",
    color = "#2e2e2e",
    linewidth = 0.3
  ) +

  # Compass Rose
  annotation_north_arrow(
    location = "bl",
    which_north = "true",
    pad_x = unit(0.4, "in"),
    pad_y = unit(0.4, "in"),
    style = north_arrow_fancy_orienteering(
      fill = c("grey30", "white"),
      line_col = "grey20"
    ),
    height = unit(1.1, "cm"),
    width = unit(1.1, "cm")
  ) +

  # Coordinate framing with 4° axis breaks
  coord_sf(
    xlim = c(-136.0, -119.5),
    ylim = c(35.0, 56.0),
    expand = FALSE,
    label_axes = list(bottom = "E", left = "N")
  ) +
  scale_x_continuous(breaks = seq(-136, -120, by = 4)) +
  scale_y_continuous(breaks = seq(36, 56, by = 4)) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
    axis.text = element_text(size = 9, color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.title = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.15, "cm"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key = element_rect(color = "grey40"),
    plot.margin = margin(8, 8, 8, 8)
  )
ggsave("figures/range_map.png", width = 7, height = 8)
