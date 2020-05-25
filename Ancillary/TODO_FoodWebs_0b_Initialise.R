# Loads the required model data

### Get the environmental data
# enviro_data <- read_rds(paste0(data_dir,"enviro200_20200209.RDS")) # Load in environmental data

TestGroups <- read_csv("~/Nextcloud/MME2Work/ZooMSS/_LatestModel/20200212_Control_Full_UNSW/TestGroups.csv") # Load in functional group information

# Lets add these colours to TestGroup
TestGroups <- TestGroups %>%
  mutate(Colour = case_when(species == "Flagellates" ~ "seagreen3",
                            species == "Ciliates" ~ "seagreen4",
                            species == "Larvaceans" ~ "lightgoldenrod4",
                            species == "OmniCopepods" ~ "cornflowerblue",
                            species == "CarnCopepods" ~ "darkred",
                            species == "Euphausiids" ~ "darkblue",
                            species == "Chaetognaths" ~ "firebrick1",
                            species == "Salps" ~ "lightgoldenrod3",
                            species == "Jellyfish" ~ "orange",
                            species == "Fish_Small" ~ "lightcoral",
                            species == "Fish_Med" ~ "coral3",
                            species == "Fish_Large" ~ "coral4"),
         CommonName = species,
         CommonName = case_when(species == "Euphausids" ~ "Krill",
                                species != "Euphausids" ~ species),
         Abbrev = case_when(species == "Flagellates" ~ "Het. Flag.",
                            species == "Ciliates" ~ "Het. Cil.",
                            species == "Larvaceans" ~ "Larv.",
                            species == "OmniCopepods" ~ "Omni.Cop.",
                            species == "CarnCopepods" ~ "Carn.Cop.",
                            species == "Euphausiids" ~ "Euph.",
                            species == "Chaetognaths" ~ "Chaet.",
                            species == "Salps" ~ "Salps",
                            species == "Jellyfish" ~ "Jelly.",
                            species == "Fish_Small" ~ "Fish_Small",
                            species == "Fish_Med" ~ "Fish_Med",
                            species == "Fish_Large" ~ "Fish_Large"),
         Feeding = case_when(species == "Flagellates" ~ "Heterotroph",
                             species == "Ciliates" ~ "Heterotroph",
                             species == "Larvaceans" ~ "FilterFeeder",
                             species == "OmniCopepods" ~ "Omnivore",
                             species == "CarnCopepods" ~ "Carnivore",
                             species == "Euphausiids" ~ "Omnivore",
                             species == "Chaetognaths" ~ "Carnivore",
                             species == "Salps" ~ "FilterFeeder",
                             species == "Jellyfish" ~ "Carnivore",
                             species == "Fish_Small" ~ "Fish",
                             species == "Fish_Med" ~ "Fish",
                             species == "Fish_Large" ~ "Fish"
         ),
         FeedingNo = case_when(species == "Flagellates" ~ 0,
                             species == "Ciliates" ~ 0,
                             species == "Larvaceans" ~ 1,
                             species == "OmniCopepods" ~ 2,
                             species == "CarnCopepods" ~ 3,
                             species == "Euphausiids" ~ 2,
                             species == "Chaetognaths" ~ 3,
                             species == "Salps" ~ 1,
                             species == "Jellyfish" ~ 3,
                             species == "Fish_Small" ~ 4,
                             species == "Fish_Med" ~ 4,
                             species == "Fish_Large" ~ 4
         )
  )
