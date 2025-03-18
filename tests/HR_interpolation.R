library(ggforce)
library(tidyverse)
library(ggtext)
library(patchwork)
data("san_andreas")
data("plates") # load plate boundary data set
data("nuvel1")

por <- subset(nuvel1, nuvel1$plate.rot == "na")
san_andreas.res <-
  bind_cols(
    san_andreas,
    PoR_shmax(san_andreas, por, type = "right")
  ) |>
  filter(
    between(lat, 31.5, 90),
    between(lon, -124, -114.5)
  )


 n_lines <- function(d){
   180/d
 }
trajectories <- eulerpole_loxodromes(por, n_lines(2), cw = FALSE)

ggplot() +
  geom_sf(
    data = plates,
    color = "red",
    lwd = 2,
    alpha = .5
  ) +

  geom_sf(
    data = trajectories,
    lty = 2
  ) +
  geom_spoke(
    data = san_andreas.res,
    aes(
      x = lon,
      y = lat,
      angle = deg2rad(90 - azi),
      color = deviation_norm(dev),
      alpha = quality
    ),
    radius = .3,
    position = "center_spoke",
    na.rm = TRUE
  ) +
  scale_color_continuous(
    type = "viridis",
    limits = c(NA, NA),
    name = "|Deviation| in (\u00B0)",
    breaks = seq(0, 90, 10)
  ) +
  scale_alpha_discrete(name = "Quality rank", range = c(1, 0.2)) +
  scale_x_continuous(breaks = seq(-360, 360, 2)) +
  scale_y_continuous(breaks = seq(-360, 360, 2)) +

  coord_sf(
    xlim = range(san_andreas.res$lon),
    ylim = range(san_andreas.res$lat)
  ) +
  # theme_classic() +
  theme(
  )


gridsize <- .25

san_andreas_int <- PoR_stress2grid_stats(san_andreas.res, por, gridsize = gridsize,  R_range = seq(25, 250, 25)) |>
  compact_grid() |>
  dplyr::mutate(
    dev = mean.PoR - 135,
    cdist = circular_distance(mean.PoR, 135)
  )

san_andreas_kdsip <-
  san_andreas.res |>
  select(-azi) |>
  rename(azi = azi.PoR) |>
  kernel_dispersion(gridsize = gridsize, R_range = seq(25, 500, 25)) |>
  compact_grid('dispersion')




ggplot(san_andreas_int) +
  geom_spoke(
    data = san_andreas_int,
    aes(
      x = lon,
      y = lat,
      angle = deg2rad(90 - mean),
      # color = mdr,
      alpha = sd
    ),
    # color = 'white',
    lwd = .2,
    radius = .15,
    position = "center_spoke",
    na.rm = TRUE
  ) +
  scale_alpha(name = "Sd. (\u00B0)", range = c(1, 0.5)) +


ggplot(san_andreas_int) +
  ggforce::geom_voronoi_tile(
    aes(lon, lat, fill = R),
    max.radius = gridsize, normalize = FALSE
  ) +
  scale_fill_viridis_b('Wavelength (km)', option='A', begin = .2) +

ggplot(san_andreas_int) +
  ggforce::geom_voronoi_tile(
    aes(lon, lat, fill = dev),
    max.radius = gridsize, normalize = FALSE
  ) +
  # scale_fill_viridis_c("Angular distance", limits = c(0, 1)) +
  scico::scale_fill_scico("Angle (\u00B0)", palette = 'vik', limits = c(-90, 90), breaks = seq(-90, 90, 22.5)) +




ggplot(san_andreas_int) +
  ggforce::geom_voronoi_tile(
    aes(lon, lat, fill = skewness),
    max.radius = gridsize, normalize = FALSE
  ) +
  scico::scale_fill_scico("Skewness", palette = 'bam', limits = c(-max(abs(san_andreas_int$skewness)), max(abs(san_andreas_int$skewness)))/2, oob = scales::squish) +


  ggplot(san_andreas_int) +
  ggforce::geom_voronoi_tile(
    aes(lon, lat, fill = abs(kurtosis)),
    max.radius = gridsize, normalize = FALSE
  ) +
  scale_fill_viridis_c("Kurtosis", option = 'E', transform = 'log10', breaks = scales::breaks_log(), labels = scales::label_log()) +


  ggplot(san_andreas_kdsip) +
  ggforce::geom_voronoi_tile(
    aes(lon, lat, fill = stat),
    max.radius = gridsize, normalize = FALSE
  ) +
  scale_fill_viridis_c("Dispersion", breaks = seq(0, 1, .2)) +

  plot_layout() &
  geom_sf(
    data = plates,
    color = "red",
    lwd = 2,
    alpha = .5
  ) &
  geom_sf(data = san_andreas, size = .25) &
  geom_sf(
    data = trajectories,
    lty = 2
  ) &

  # scale_alpha_discrete(name = "Quality rank", range = c(1, 0.2)) +
  scale_x_continuous(breaks = seq(-360, 360, 2)) &
  scale_y_continuous(breaks = seq(-360, 360, 2)) &
  labs(x=NULL, y = NULL) &
  coord_sf(
    xlim = range(san_andreas.res$lon),
    ylim = range(san_andreas.res$lat),
    expand = FALSE
  ) &
  theme_classic() +
  theme(
  )
