## code to prepare `san_andreas` dataset goes here
data("wsm2016")
library(dplyr)
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
    quality.quant = tectonicr::quantise_wsm_quality(quality),
    unc = ifelse(is.na(sd), quality.quant, sd),
    unc = ifelse(unc > quality.quant, quality.quant, unc),
    unc = ifelse(unc == 0, 15, unc)
  ) %>%
  arrange(quality, unc) %>%
  filter(quality != "E") %>%
  sf::st_as_sf(coords = c('x', 'y'), crs = "WSG84") %>%
  select(id, lat, lon, azi, unc, type, depth, quality, regime)

usethis::use_data(san_andreas, overwrite = TRUE)
