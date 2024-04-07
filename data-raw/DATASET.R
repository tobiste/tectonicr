## code to prepare `san_andreas` dataset goes here  -------------

wsm2016 <-
  dplyr::mutate(wsm2016,
    quality = forcats::fct_relevel(quality, "A", "B", "C", "D", "E"),
    regime = ifelse(regime == "NF", "N", regime),
    regime = ifelse(regime == "TF", "T", regime),
    regime = ifelse(regime == "SS", "S", regime),
    regime = ifelse(regime == "U", NA, regime),
    quality.quant = tectonicr::quantise_wsm_quality(quality),
    unc = ifelse(is.na(sd), quality.quant, sd),
    unc = ifelse(unc > quality.quant, quality.quant, unc),
    unc = ifelse(unc == 0, 15, unc),
  ) |>
  dplyr::arrange(quality, unc) |>
  dplyr::filter(quality != "E") |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = sf::st_crs("WGS84"), remove = FALSE) |>
  dplyr::select(id, lat, lon, azi, unc, type, depth, quality, regime)

san_andreas <- filter(
  wsm2016,
  between(lat, 23, 40) &
    between(lon, -126, -108)
)
usethis::use_data(san_andreas, overwrite = TRUE, ascii = TRUE)

frame_iceland <- readRDS("~/GIT_repos/europe-tectonics/data-raw/iceland_frame.rds")
iceland <- sf::st_intersection(wsm2016, frame_iceland)
usethis::use_data(iceland, overwrite = TRUE, ascii = TRUE)

zoom_asia <- readRDS("~/GIT_repos/europe-tectonics/data-raw/asia_zoom.rds")
tibet <- sf::st_intersection(wsm2016, zoom_asia)
usethis::use_data(tibet, overwrite = TRUE, ascii = TRUE)


# ggplot(san_andreas) +
#   geom_spoke(
#     aes(
#       x = lon,
#       y = lat,
#       radius = 1,
#       angle = deg2rad(90-azi),
#       alpha = san_andreas$quality
#     ),
#     position = "center_spoke",
#     na.rm = TRUE
#   ) +
#   coord_sf(
#     xlim = range(san_andreas.res$lon),
#     ylim = range(san_andreas.res$lat),
#     expand = FALSE,
#     default_crs = "WGS84"
#   ) +
#   scale_alpha_discrete(name = "Quality rank", range = c(1, 0.1))


# plate boundaries --------------------------------------------------------

## PB2002 plates --------------------------------------------------------
# data(PB2002)
# pb2002_plates <- PB2002 |> sf::st_make_valid() |>
#  rename(name = Name, source = Source, plateA = PlateA, plateB = PlateB, type = Type) |>
#  select(-LAYER, name, plateA, plateB, type, source)
# usethis::use_data(pb2002_plates, overwrite = TRUE)

plate.names <- readxl::read_excel("../europe-tectonics/data/euler/recent_plate_motion.xlsx", sheet = "Bird") |>
  dplyr::select(plate.rot, plate.name) |>
  dplyr::mutate(plate.rot = tolower(plate.rot))

pbty.def <- data.frame(
  pb.type = c("divergent", "convergent", "transform_L", "transform_R"),
  pb.type.stress = c("out", "in", "left", "right")
)
plates <- sf::read_sf("../europe-tectonics/data/gis/PB2002_mod.shp") |>
  sf::st_set_crs("EPSG:4326") |>
  sf::st_make_valid() |>
  sf::st_wrap_dateline() |>
  dplyr::rename(plate.pair = plat_pr, pb.type = type) |>
  dplyr::mutate(pb = paste0(plate.pair, "_", pb.type)) |>
  dplyr::left_join(pbty.def) |>
  tidyr::drop_na(pb.type) |>
  sf::st_as_sf() |>
  tidyr::separate(plate.pair, into = c("PlateA", "PlateB"), sep = "-", remove = FALSE) |>
  dplyr::left_join(plate.names, by = c("PlateA" = "plate.rot")) |>
  dplyr::left_join(plate.names, by = c("PlateB" = "plate.rot"), suffix = c("A", "B")) |>
  dplyr::rename(type = pb.type, pair = plate.pair, plateA = PlateA, plateB = PlateB, name = pb, displacement = pb.type.stress, nameA = plate.nameA, nameB = plate.nameB) |>
  dplyr::select(pair, plateA, plateB, type, displacement, name, nameA, nameB) |>
  dplyr::arrange(pair, displacement) |>
  dplyr::group_by(pair)
# plot(plates)

usethis::use_data(plates, overwrite = TRUE, ascii = TRUE)


## NUVEL1 plates --------------------------------------------------------
# nuvel1_plates <- sf::st_read("E:/Global_data/Plate Boundaries/NUVEL1/nuvel1_plates.shp", quiet = TRUE) |>
#   sf::st_make_valid() |>
#   rename(plateA = PlateA, plateB = PlateB)
# data("nuvel1_plates")
# plot(nuvel1_plates)
# usethis::use_data(nuvel1_plates, overwrite = TRUE, ascii = TRUE)

## GSRM2 plates --------------------------------------------------------
# dat <- read.table("../cordillera-stress/data/gsrm2/GSRM_plate_outlines.gmt") |>
#   dplyr::mutate(plate = NA)
# for(i in seq_along(dat$V1)){
#   if(dat$V1[i] == ">") {
#     dat$plate[i] <-  dat$V2[i]
#   } else {
#     dat$plate[i] <-  dat$plate[i-1]
#   }
# }
# gsrm2_plates <- dat |> dplyr::filter(V1 != ">") |>
#   dplyr::mutate(lon = as.numeric(V1), lat = as.numeric(V2)) |>
#   dplyr::select(-c(V1, V2)) |>
#   dplyr::as_tibble() |>
#   sf::st_as_sf(coords = c("lon", "lat"), crs = "WGS84", remove = FALSE) |>
#   dplyr::group_by(plate) |>
#   dplyr::summarize(plate = dplyr::first(plate),do_union=FALSE) |>
#   sf::st_cast("LINESTRING") |>
#   sf::st_wrap_dateline()
#
# usethis::use_data(gsrm2_plates, overwrite = TRUE, ascii = TRUE)

## Rotation parameters --------------
### NUVEL 1
data("nuvel1")
nuvel1$plate.name <- stringi::stri_enc_toascii(nuvel1$plate.name)
nuvel1$plate.rot <- stringi::stri_enc_toascii(nuvel1$plate.rot)
nuvel1$plate.fix <- stringi::stri_enc_toascii(nuvel1$plate.fix)
nuvel1$source <- stringi::stri_enc_toascii(nuvel1$source)
usethis::use_data(nuvel1, overwrite = TRUE, ascii = TRUE)


pb2002 <-
  readxl::read_excel("../europe-tectonics/data/euler/recent_plate_motion.xlsx", sheet = "Bird") |>
  dplyr::mutate(plate.rot = tolower(plate.rot), source = ifelse(source == "this paper", "Bird [2003]", source), model = "PB2002") |>
  dplyr::rename(area = Area) |>
  dplyr::select(plate.name, plate.rot, lon, lat, angle, plate.fix, model) |>
  dplyr::mutate(
    plate.name = stringi::stri_enc_toascii(plate.name),
    plate.rot = stringi::stri_enc_toascii(plate.rot),
    plate.rot = stringi::stri_enc_toascii(plate.rot),
    plate.fix = stringi::stri_enc_toascii(plate.fix),
    model = stringi::stri_enc_toascii(model)
  ) |>
  as.data.frame()
usethis::use_data(pb2002, overwrite = TRUE, ascii = TRUE)

morvel56 <-
  readxl::read_excel("../europe-tectonics/data/euler/recent_plate_motion.xlsx", sheet = "NNR-MORVEL56") |>
  dplyr::mutate(plate.rot = tolower(plate.rot), model = "NNR-MORVEL56") |>
  # rename(area = Area) |>
  dplyr::select(plate.name, plate.rot, lon, lat, angle, plate.fix, model)


gsrm2 <-
  readxl::read_excel("../europe-tectonics/data/euler/recent_plate_motion.xlsx", sheet = "GSRM") |>
  dplyr::mutate(
    plate.rot = tolower(plate.rot),
    lat = lat_NNR,
    lon = lon_NNR,
    angle = angle_NNR,
    plate.fix = "NNR",
    model = "GSRM2.1"
  ) |>
  dplyr::select(plate.name, plate.rot, lon, lat, angle, plate.fix, model)

hsnuvel1a <- nuvel.hs <-
  readxl::read_excel("../europe-tectonics/data/euler/recent_plate_motion.xlsx", sheet = "HS3-NUVEL1A") |>
  dplyr::mutate(
    plate.rot = tolower(plate.rot),
    angle = rate,
    plate.fix = "hs",
    model = "HS3-NUVEL1A"
  ) |>
  dplyr::select(plate.name, plate.rot, lon, lat, angle, plate.fix, model)

revel <-
  readxl::read_excel("../europe-tectonics/data/euler/recent_plate_motion.xlsx", sheet = "REVEL") |>
  dplyr::mutate(
    plate.rot = tolower(plate.rot),
    plate.fix = "itrf97",
    model = "REVEL"
  ) |>
  dplyr::select(plate.name, plate.rot, lon, lat, angle, plate.fix, model)

cpm_models <- rbind(
  nuvel1 |> dplyr::mutate(model = "NNR-NUVEL1A") |> sdplyr::elect(-source),
  morvel56,
  gsrm2,
  hsnuvel1a,
  revel,
  pb2002
) # |>
# mutate(plate.name = stringi::stri_enc_toascii(plate.name),
#        plate.rot = stringi::stri_enc_toascii(plate.rot),
#        plate.rot = stringi::stri_enc_toascii(plate.rot),
#        plate.fix = stringi::stri_enc_toascii(plate.fix),
#        model = stringi::stri_enc_toascii(model)
# ) #|> group_by(model)
usethis::use_data(cpm_models, overwrite = TRUE, ascii = TRUE)

borders <- rnaturalearth::ne_download(returnclass = "sf") |> sf::st_geometry()
usethis::use_data(borders, overwrite = TRUE, ascii = TRUE, compress = "xz")
