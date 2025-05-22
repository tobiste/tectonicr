## code to prepare `san_andreas` dataset goes here  -------------
library(tidyverse)

wsm2025 <- download_WSM(version = '2016') |>
  dplyr::select(id, lat, lon, azi, unc, type, depth, quality, regime)

all(validUTF8(wsm2025$id))
all(validUTF8(wsm2025$type))
all(validUTF8(wsm2025$regime))

san_andreas <- filter(
  wsm2025,
  between(lat, 23, 40),
  between(lon, -126, -108),
  !is.na(azi)
)
usethis::use_data(san_andreas, overwrite = TRUE, ascii = TRUE)

frame_iceland <- readRDS("../europe-tectonics/data-raw/iceland_frame.rds")
iceland <- sf::st_intersection(wsm2025, frame_iceland) |>
  filter(!is.na(azi)) |>
  select(-c(Name, Description))
usethis::use_data(iceland, overwrite = TRUE, ascii = TRUE)

zoom_asia <- readRDS("../europe-tectonics/data-raw/asia_zoom.rds")
tibet <- sf::st_intersection(wsm2025, zoom_asia) |>
  filter(!is.na(azi)) |>
  select(-c(Name, Description))
usethis::use_data(tibet, overwrite = TRUE, ascii = TRUE)


# ggplot(san_andreas) +
#   geom_spoke(
#     aes(
#       x = lon,
#       y = lat,
#       radius = 1,
#       angle = deg2rad(90-azi),
#       alpha = quality
#     ),
#     position = "center_spoke",
#     na.rm = TRUE
#   ) +
#   coord_sf(
#     xlim = range(san_andreas$lon),
#     ylim = range(san_andreas$lat),
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

plate.names <- readxl::read_excel("inst/recent_plate_motion.xlsx", sheet = "Bird") |>
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

all(validUTF8(plates$displacement))

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
# nuvel1$plate.name <- stringi::stri_enc_toascii(nuvel1$plate.name)
# nuvel1$plate.rot <- stringi::stri_enc_toascii(nuvel1$plate.rot)
# nuvel1$plate.fix <- stringi::stri_enc_toascii(nuvel1$plate.fix)
# nuvel1$source <- stringi::stri_enc_toascii(nuvel1$source)
# all(validUTF8(nuvel1$plate.name))
# usethis::use_data(nuvel1, overwrite = TRUE, ascii = TRUE)
nnr.nuvel1a <- nuvel1

pb2002 <-
  readxl::read_excel("inst/recent_plate_motion.xlsx", sheet = "Bird") |>
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
all(validUTF8(pb2002$plate.name))
usethis::use_data(pb2002, overwrite = TRUE, ascii = TRUE)

morvel56 <-
  readxl::read_excel("inst/recent_plate_motion.xlsx", sheet = "NNR-MORVEL56") |>
  dplyr::mutate(plate.rot = tolower(plate.rot), model = "NNR-MORVEL56", plate.name = str_squish(plate.name)) |>
  # rename(area = Area) |>
  dplyr::select(plate.name, plate.rot, lon, lat, angle, plate.fix, model)
all(validUTF8(morvel56$plate.name))


gsrm2 <-
  readxl::read_excel("inst/recent_plate_motion.xlsx", sheet = "GSRM") |>
  dplyr::mutate(
    plate.rot = tolower(plate.rot),
    lat = lat_NNR,
    lon = lon_NNR,
    angle = angle_NNR,
    plate.fix = "NNR",
    model = "GSRM2.1"
  ) |>
  dplyr::select(plate.name, plate.rot, lon, lat, angle, plate.fix, model)

nuvel1 <- readxl::read_excel("inst/recent_plate_motion.xlsx", sheet = "NUVEL1") |>
  dplyr::mutate(
    plate.rot = tolower(plate.rot),
    angle = rate,
    # plate.fix = "pa",
    model = "NUVEL1"
  ) |>
  dplyr::select(plate.name, plate.rot, lon, lat, angle, plate.fix, model)

hs3nuvel1a <- nuvel.hs <-
  readxl::read_excel("inst/recent_plate_motion.xlsx", sheet = "HS3-NUVEL1A") |>
  dplyr::mutate(
    plate.rot = tolower(plate.rot),
    angle = rate,
    # plate.fix = "hs",
    model = "HS3-NUVEL1A"
  ) |>
  dplyr::select(plate.name, plate.rot, lon, lat, angle, plate.fix, model)

hs2nuvel1 <-
  readxl::read_excel("inst/recent_plate_motion.xlsx", sheet = "HS2-NUVEL1") |>
  dplyr::mutate(
    plate.rot = tolower(plate.rot),
    angle = rate,
    # plate.fix = "hs",
    model = "HS2-NUVEL1"
  ) |>
  dplyr::select(plate.name, plate.rot, lon, lat, angle, plate.fix, model)

p073 <- readxl::read_excel("inst/recent_plate_motion.xlsx", sheet = "P073") |>
  dplyr::mutate(
    plate.rot = tolower(plate.rot),
    angle = rate, #* 10^(-7),
    # plate.fix = "hs",
    model = "P073",
  ) |>
  dplyr::select(plate.name, plate.rot, lon, lat, angle, plate.fix, model)

am <- readxl::read_excel("inst/recent_plate_motion.xlsx", sheet = "AM1-2") |>
  dplyr::mutate(
    plate.rot = tolower(plate.rot),
    angle = rate, #* 10^(-7),
    # plate.fix = "hs",
    model = "AM1-2",
  ) |>
  dplyr::select(plate.name, plate.rot, lon, lat, angle, plate.fix, model)

revel <-
  readxl::read_excel("inst/recent_plate_motion.xlsx", sheet = "REVEL") |>
  dplyr::mutate(
    plate.rot = tolower(plate.rot),
    # plate.fix = "itrf97",
    model = "REVEL"
  ) |>
  dplyr::select(plate.name, plate.rot, lon, lat, angle, plate.fix, model)

ITRF2020 <- readxl::read_excel("inst/recent_plate_motion.xlsx", sheet = "ITRF2020-PMM")
ITRF2020_cart <- euler::euler_cart2geo(ITRF2020$x_deg, ITRF2020$y_deg, ITRF2020$z_deg)


ITRF2020_PMM <- ITRF2020 |>
  select(plate.name, plate.rot, plate.fix) |>
  bind_cols(ITRF2020_cart) |>
  rename(angle = mag) |>
  dplyr::mutate(
    lon = ifelse(is.nan(lon), 0, lon),
    lat = ifelse(is.nan(lat), 0, lat),
    plate.rot = tolower(plate.rot),
    # plate.fix = "itrf97",
    model = "ITRF2020-PMM"
  ) |>
  dplyr::select(plate.name, plate.rot, lon, lat, angle, plate.fix, model)



cpm_models_df <- rbind(
  nnr.nuvel1a |> dplyr::mutate(model = "NNR-NUVEL1A") |> dplyr::select(-source),
  nuvel1,
  am,
  morvel56,
  gsrm2,
  hs3nuvel1a,
  hs2nuvel1,
  revel,
  pb2002,
  p073,
  ITRF2020_PMM
) |>
  mutate(
    plate.name = stringi::stri_enc_toascii(plate.name),
    plate.rot = stringi::stri_enc_toascii(plate.rot),
    plate.rot = stringi::stri_enc_toascii(plate.rot),
    plate.fix = stringi::stri_enc_toascii(plate.fix),
    model = stringi::stri_enc_toascii(model)
  ) |>
  arrange(model)

cpm_models <- cpm_models_df |>
  group_by(model) |>
  group_split()
names(cpm_models) <- unique(cpm_models_df$model)

usethis::use_data(cpm_models, overwrite = TRUE, ascii = TRUE)
# borders <- rnaturalearth::ne_download(returnclass = "sf") |> dplyr::select() |> sf::st_geometry()
# usethis::use_data(borders, overwrite = TRUE, ascii = F, compress = "xz",version = 3)

# check all the files:
tools::checkRdaFiles("data")
