lapply(c("dplyr", "ggplot"), require, character.only = T)

#Mortality assessment

trap.raw <- read.csv("Trapping - Data.csv") %>%
  select(Bird.ID = AlumBand, Date, Trans.Type, Sex) %>%
  mutate(Fate = "L") %>%
  filter(Trans.Type == "GPS_Back" & Sex == "F") %>%
  select(Bird.ID, Date, Fate)

telem.raw <- read.csv("Telemetry_Data - Telemetry.csv") %>%
  select(Bird.ID = AlumBand, Date, Fate) %>%
  filter(Bird.ID %in% trap.raw$Bird.ID)

all.checks <- rbind(trap.raw, telem.raw) %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>%
  arrange(Bird.ID, Date)

first.L <- all.checks %>%
  filter(Fate == "L") %>%
  group_by(Bird.ID) %>%
  slice(1L) %>%
  ungroup() %>%
  rename(Cap = Date) %>%
  select(-Fate)

all.D <- all.checks %>%
  filter(Fate == "D") %>%
  group_by(Bird.ID) %>%
  slice(1L) %>%
  ungroup() %>%
  rename(Death = Date) %>%
  select(-Fate)

first.LD <- merge(first.L, all.D, all.x = T, by = "Bird.ID") %>%
  arrange(Bird.ID) %>%
  mutate(Days.Alive = Death-Cap) %>%
  mutate(CapMort = ifelse(Days.Alive < 14, 1, 0)) %>%
  mutate(CapMort = ifelse(is.na(CapMort), 0, CapMort)) %>%
  filter(CapMort != 1) %>%
  mutate(Winter = as.Date(paste(lubridate::year(Cap), "-03-15", sep = ""))) %>%
  mutate(WintDeath = ifelse(Winter > Death, 1, 0)) %>%
  mutate(WintDeath = ifelse(is.na(WintDeath), 0, WintDeath))
  
