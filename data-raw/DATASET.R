## code to prepare `san_andreas` dataset goes here  -------------
library(dplyr)

data("wsm2016", package = "ptrotR")
san_andreas <- filter(
  wsm2016,
  between(lat, 23, 40) &
    between(lon, -126, -108)
) %>%
  mutate(
    x = lon, y = lat,
    quality = forcats::fct_relevel(quality, "E", "D", "C", "B", "A"),
    regime = ifelse(regime == "NF", "N", regime),
    regime = ifelse(regime == "TF", "T", regime),
    regime = ifelse(regime == "SS", "S", regime),
    regime = ifelse(regime == "U", NA, regime),
    quality.quant = tectonicr::quantise_wsm_quality(quality),
    unc = ifelse(is.na(sd), quality.quant, sd),
    unc = ifelse(unc > quality.quant, quality.quant, unc),
    unc = ifelse(unc == 0, 15, unc),

    # id = stringi::stri_enc_toascii(id),
    # type = stringi::stri_enc_toascii(type),
    # quality = stringi::stri_enc_toascii(quality),
    # regime = stringi::stri_enc_toascii(regime)
  ) %>%
  arrange(quality, unc) %>%
  filter(quality != "E") %>%
  sf::st_as_sf(coords = c("x", "y"), crs = sf::st_crs("WGS84")) %>%
  select(id, lat, lon, azi, unc, type, depth, quality, regime)
usethis::use_data(san_andreas, overwrite = TRUE, ascii = TRUE)


## PB2002 plates --------------
# data(PB2002)
# pb2002_plates <- PB2002 %>% sf::st_make_valid() %>%
#  rename(name = Name, source = Source, plateA = PlateA, plateB = PlateB, type = Type) %>%
#  select(-LAYER, name, plateA, plateB, type, source)
# usethis::use_data(pb2002_plates, overwrite = TRUE)

plate.names <- readxl::read_excel("../europe-tectonics/data/euler/recent_plate_motion.xlsx", sheet = "Bird") %>%
  select(plate.rot, plate.name) %>%
  mutate(plate.rot = tolower(plate.rot))

pbty.def <- data.frame(
  pb.type = c("divergent", "convergent", "transform_L", "transform_R"),
  pb.type.stress = c("out", "in", "left", "right")
)
plates <- sf::read_sf("../europe-tectonics/data/gis/PB2002_mod.shp") %>%
  sf::st_set_crs("EPSG:4326") %>%
  sf::st_make_valid() %>%
  sf::st_wrap_dateline() %>%
  rename(plate.pair = plat_pr, pb.type = type) %>%
  mutate(pb = paste0(plate.pair, "_", pb.type)) %>%
  left_join(pbty.def) %>%
  tidyr::drop_na(pb.type) %>%
  sf::st_as_sf() %>%
  tidyr::separate(plate.pair, into = c("PlateA", "PlateB"), sep = "-", remove = FALSE) %>%
  left_join(plate.names, by = c("PlateA" = "plate.rot")) %>%
  left_join(plate.names, by = c("PlateB" = "plate.rot"), suffix = c("A", "B")) %>%
  rename(type = pb.type, pair = plate.pair, plateA = PlateA, plateB = PlateB, name = pb, displacement = pb.type.stress, nameA = plate.nameA, nameB = plate.nameB) %>%
  select(pair, plateA, plateB, type, displacement, name, nameA, nameB) %>%
  arrange(pair, displacement) %>%
  group_by(pair)
# plot(plates)

usethis::use_data(plates, overwrite = TRUE, ascii = TRUE)


## NUVEL1 plates -------------
# nuvel1_plates <- sf::st_read("E:/Global_data/Plate Boundaries/NUVEL1/nuvel1_plates.shp", quiet = TRUE) %>%
#   sf::st_make_valid() %>%
#   rename(plateA = PlateA, plateB = PlateB)
# #data("nuvel1_plates")
# #plot(nuvel1_plates)
# usethis::use_data(nuvel1_plates, overwrite = TRUE, ascii = TRUE)

## Rotation parameters --------------
### NUVEL 1
data("nuvel1")
nuvel1$plate.name <- stringi::stri_enc_toascii(nuvel1$plate.name)
nuvel1$plate.rot <- stringi::stri_enc_toascii(nuvel1$plate.rot)
nuvel1$plate.fix <- stringi::stri_enc_toascii(nuvel1$plate.fix)
nuvel1$source <- stringi::stri_enc_toascii(nuvel1$source)
usethis::use_data(nuvel1, overwrite = TRUE, ascii = TRUE)


pb2002 <-
  readxl::read_excel("../europe-tectonics/data/euler/recent_plate_motion.xlsx", sheet = "Bird") %>%
  mutate(plate.rot = tolower(plate.rot), source = ifelse(source == "this paper", "Bird [2003]", source), model = "PB2002") %>%
  rename(area = Area) %>%
  select(plate.name, plate.rot, lon, lat, angle, plate.fix, model) %>%
  mutate(
    plate.name = stringi::stri_enc_toascii(plate.name),
    plate.rot = stringi::stri_enc_toascii(plate.rot),
    plate.rot = stringi::stri_enc_toascii(plate.rot),
    plate.fix = stringi::stri_enc_toascii(plate.fix),
    model = stringi::stri_enc_toascii(model)
  ) %>%
  as.data.frame()
usethis::use_data(pb2002, overwrite = TRUE, ascii = TRUE)

morvel56 <-
  readxl::read_excel("../europe-tectonics/data/euler/recent_plate_motion.xlsx", sheet = "NNR-MORVEL56") %>%
  mutate(plate.rot = tolower(plate.rot), model = "NNR-MORVEL56") %>%
  # rename(area = Area) %>%
  select(plate.name, plate.rot, lon, lat, angle, plate.fix, model)


gsrm2 <-
  readxl::read_excel("../europe-tectonics/data/euler/recent_plate_motion.xlsx", sheet = "GSRM") %>%
  mutate(
    plate.rot = tolower(plate.rot),
    lat = lat_NNR,
    lon = lon_NNR,
    angle = angle_NNR,
    plate.fix = "NNR",
    model = "GSRM2.1"
  ) %>%
  select(plate.name, plate.rot, lon, lat, angle, plate.fix, model)

hsnuvel1a <- nuvel.hs <-
  readxl::read_excel("../europe-tectonics/data/euler/recent_plate_motion.xlsx", sheet = "HS3-NUVEL1A") %>%
  mutate(
    plate.rot = tolower(plate.rot),
    angle = rate,
    plate.fix = "hs",
    model = "HS3-NUVEL1A"
  ) %>%
  select(plate.name, plate.rot, lon, lat, angle, plate.fix, model)

revel <-
  readxl::read_excel("../europe-tectonics/data/euler/recent_plate_motion.xlsx", sheet = "REVEL") %>%
  mutate(
    plate.rot = tolower(plate.rot),
    plate.fix = "itrf97",
    model = "REVEL"
  ) %>%
  select(plate.name, plate.rot, lon, lat, angle, plate.fix, model)

cpm_models <- rbind(
  nuvel1 %>% mutate(model = "NNR-NUVEL1A") %>% select(-source),
  morvel56,
  gsrm2,
  hsnuvel1a,
  revel,
  pb2002
) # %>%
# mutate(plate.name = stringi::stri_enc_toascii(plate.name),
#        plate.rot = stringi::stri_enc_toascii(plate.rot),
#        plate.rot = stringi::stri_enc_toascii(plate.rot),
#        plate.fix = stringi::stri_enc_toascii(plate.fix),
#        model = stringi::stri_enc_toascii(model)
# ) #%>% group_by(model)
usethis::use_data(cpm_models, overwrite = TRUE, ascii = TRUE)
