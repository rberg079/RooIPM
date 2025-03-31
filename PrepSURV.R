# 25 March 2025
# Prep survival data

library(readxl)
library(tidyverse)

# setwd("C:/Code/RooIPM")

## Load & clean up -------------------------------------------------------------

surv <- read_excel("data/PromSurvivalOct24.xlsx", sheet = "YEARLY SURV")
env  <- read_csv("data/Env_Mar25.csv")

# surv <- SURV
surv <- surv %>% 
  rename(Dead = "Found dead") %>%
  rename_with(~gsub("(\\d{2})to(\\d{2})", "in20\\2", .), everything()) %>% 
  # split "08to09" columns into "in2008" & "in2009" columns
  # 1s in "in20XX" should be read as "seen in Aug-Nov 20XX"
  mutate(in2008 = ifelse(!is.na(in2009), 1, NA)) %>%
  mutate(across(in2009:in2023, ~ifelse(is.na(.x) & !is.na(
    get(sub("\\d{4}", as.integer(gsub("\\D", "", cur_column())) + 1, cur_column()))), 1, .x))) %>%
  select(ID, Sex, Dead, in2008, matches("^in20\\d{2}$"), -in2025, matches("^Age\\d{2}")) %>%
  filter(Sex == 2) # females


## obs data --------------------------------------------------------------------

obs <- surv %>% 
  select(in2008:in2024) %>% 
  mutate_all(~ifelse(.> 1 | is.na(.), 0, .))

# write csv
# write_csv(obs, "obsF.csv")

# roos observed outside of the study area
# as well as roadkilled or poached roos
# are treated as unobserved


## state & id data -------------------------------------------------------------

state <- surv %>% 
  mutate(HRDead = as.numeric(apply(select(., starts_with("in")) == 2, 1, any)),
         HRDead = replace_na(HRDead, 0)) %>% 
  select(ID, Dead, HRDead, in2008:in2024) %>% 
  mutate_at(vars(in2008:in2024), ~case_when(. == 2 ~ 0,   # roadkills are dead
                                            . == 3 ~ 1,   # observed emigrants are alive
                                            . == 4 ~ NA,  # unobserved emigrants are unknown
                                            TRUE ~ .))

# roos that were unobserved one year
# but observed the next were alive
state <- state %>% 
  mutate(across(in2008:in2023, ~ifelse(. == 0 & !is.na(
    get(sub("\\d{4}", as.integer(gsub("\\D", "", cur_column())) + 1, cur_column()))),
    1, .)))

# split into roos found dead vs vanished
# NAs become 0s for roos found dead,
# but remain NAs for vanished roos
state.found <- state %>% 
  filter(!is.na(Dead))

for(i in which(is.na(state.found[,20]))){
  state.found[i,(max(which(state.found[i,] == 1)) +1):ncol(state.found)] <- 0
}

state.vanished <- state %>% 
  filter(is.na(Dead)) %>% 
  mutate_at(vars(in2008:in2024), ~na_if(., 0))

# join all roos again
# create first & last variables
state <- bind_rows(state.found, state.vanished) %>%
  arrange(., ID) %>% 
  rowwise() %>% 
  mutate(first = min(which(c_across(4:20) == 1)),
         last = min(which(c_across(4:20) == 0)),
         last = ifelse(HRDead == 1, last-1, last),
         last = ifelse(last == Inf, ncol(obs), last))

# set aside id info
id <- state %>% select(ID, first, last)

# write state csv
state <- state %>% select(in2008:in2024)
# write_csv(state, "stateF.csv")
remove(state.found, state.vanished)


## age data --------------------------------------------------------------------

age <- surv %>% 
  select(Age08:Age24) %>% 
  mutate_at(vars(Age08:Age24), ~ifelse(. == "A", NA, .)) %>% 
  mutate_all(~as.numeric(.))

fill_ages <- function(row) {
  if(all(is.na(row))) return(row)  # if all NAs, return as is
  first <- which(!is.na(row))[1]   # first non-NA value
  
  # fill forward from first known age
  row[first:length(row)] <- seq(from = row[first], by = 1, length.out = length(row) - first + 1)
  
  # fill backward if needed
  if(first > 1) {
    row[1:(first - 1)] <- seq(from = row[first] - first + 1, by = 1, length.out = first - 1)
  }
  return(row)
}

age <- as.data.frame(t(apply(age, 1, fill_ages)))
# write_csv(age, "ageF.csv")

# add uka & write id csv
id$uka <- as.logical(is.na(age[,ncol(age)]))
# write_csv(id, "idF.csv")


## env data --------------------------------------------------------------------

env <- env %>% 
  filter(Date > "2008-03-31") %>% 
  mutate(Year = ifelse(Month < 10, Year-1, Year)) %>%
  distinct(Date, Year, Month, Day, Dens, DensSE, Veg, VegSE)

# create yearly dens, veg & veg/roo
# propagate uncertainty throughout
# additions: if Q = A + B, error becomes eQ = sqrt((eA)^2 + (eB)^2)
# multiplications: if Q = A*B, error eQ = sqrt((eA/A)^2 + (eB/B)^2)

# density estimates
dens <- env %>% 
  filter(!is.na(Dens)) %>% 
  distinct(Year, Dens, DensSE) %>% 
  mutate(DensSE = DensSE^2) %>% 
  group_by(Year) %>% 
  mutate(Dens = mean(Dens),
         DensSE = sqrt(sum(DensSE)) / 5) %>% 
  ungroup() %>% 
  distinct(Year, Dens, DensSE)

# vegetation estimates
veg <- env %>% 
  filter(!is.na(Veg)) %>% 
  select(Date, Year, Veg, VegSE) %>% 
  group_by(Year) %>% 
  mutate(Veg = sum(Veg),
         VegSE = sqrt(sum(VegSE^2)),
         across(c(Veg, VegSE), ~replace(., Year == 2008, NA)),
         across(c(Veg, VegSE), ~replace(., Year == 2024, NA))) %>% 
  ungroup() %>% 
  distinct(Year, Veg, VegSE)

# veg/roo estimates
env <- veg %>% 
  left_join(dens) %>% 
  mutate(VegRoo = Veg / Dens,
         VegRooSE = abs(VegRoo)*sqrt((VegSE / Veg)^2+(DensSE / Dens)^2))

# # plot all of this for fun
# g1 <- ggplot(env, aes(x = Year)) +
#   geom_point(aes(y = Dens)) +
#   geom_ribbon(aes(ymin = Dens-DensSE, ymax = Dens+DensSE), alpha = 0.2) +
#   theme_bw(); g1
# 
# g2 <- ggplot(env, aes(x = Year)) +
#   geom_point(aes(y = Veg), colour = '#738844') +
#   geom_ribbon(aes(ymin = Veg-VegSE, ymax = Veg+VegSE), fill = "#738844", alpha = 0.2) +
#   theme_bw(); g2
# 
# g3 <- ggplot(env, aes(x = Year)) +
#   geom_point(aes(y = VegRoo), colour = '#A44200') +
#   geom_ribbon(aes(ymin = VegRoo-VegRooSE, ymax = VegRoo+VegRooSE), fill = "#A44200", alpha = 0.2) +
#   theme_bw(); g3
# 
# library(cowplot)
# plot_grid(g1 + labs(y = "Population density (roos/ha)"),
#           g2 + labs(y = "Available forage (g/m^2/year)"),
#           g3 + labs(y = "Available forage per capita"), nrow = 1)

# # see what a total area decrease following roughly
# # area = -0.026(year)+0.75 means for veg & veg/roo
# env <- env %>% 
#   mutate(Area = as.numeric(factor(Year)),
#          Area = (-0.026*Area+0.75)*100,
#          Roos = Dens*Area,
#          VegArea = Veg*10000*Area,
#          VegRoos = VegArea/Roos)
# 
# g4 <- ggplot(env, aes(x = Year)) +
#   geom_point(aes(y = Roos)) +
#   geom_path(aes(y = Roos), alpha = 0.4) +
#   theme_bw(); g4
# 
# g5 <- ggplot(env, aes(x = Year)) +
#   geom_point(aes(y = VegArea), colour = '#738844') +
#   geom_path(aes(y = VegArea), colour = "#738844", alpha = 0.4) +
#   theme_bw(); g5
# 
# g6 <- ggplot(env, aes(x = Year)) +
#   geom_point(aes(y = VegRoos), colour = '#A44200') +
#   geom_path(aes(y = VegRoos), colour = "#A44200", alpha = 0.4) +
#   theme_bw(); g6
# 
# plot_grid(g4 + labs(y = "Population size (roos)"),
#           g5 + labs(y = "Available forage over whole area (g)"),
#           g6 + labs(y = "Available forage per capita (g/roo)"), nrow = 1)


## nimble ----------------------------------------------------------------------

# assemble nimble lists
nind   <- nrow(state)
ntimes <- ncol(state)

ageC   <- c(1,2,2,3,3,3,3,4,4,4, rep(5,60))
nAge   <- max(ageC, na.rm = T)
noAge  <- which(is.na(age[,ncol(age)]))
nNoAge <- length(noAge)

first <- as.numeric(id$first)
last  <- as.numeric(id$last)

obs   <- as.matrix(obs) %>% unname()
state <- as.matrix(state) %>% unname()
age   <- as.matrix(age)+1 %>% unname()

veg  <- as.numeric(scale(env$Veg))
dens <- as.numeric(scale(env$Dens))

vegE  <- as.numeric(ifelse(is.na(env$VegSE), 2, env$VegSE/sd(env$Veg, na.rm = T)))
densE <- as.numeric(ifelse(is.na(env$DensSE), 2, env$DensSE/sd(env$Dens, na.rm = T)))

nNoVeg  <- sum(is.na(veg))
nNoDens <- sum(is.na(dens))

mydata  <- list(obs = obs, state = state, age = age, ageC = ageC,
                dens = dens, densE = densE, veg = veg, vegE = vegE)

myconst <- list(nind = nind, ntimes = ntimes, 
                nAge = nAge, noAge = noAge, nNoAge = nNoAge,
                first = first, last = last, W = diag(nAge), DF = nAge,
                nNoVeg = nNoVeg, nNoDens = nNoDens)

