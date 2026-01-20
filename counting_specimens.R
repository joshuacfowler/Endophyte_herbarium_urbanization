#how many specimens? from which herbaria? let's find out
library(dplyr)


herb_info <- read.csv(file = "C:/Users/malpa/Downloads/Endo_Herbarium - Herbarium specimen info.csv")%>%
  mutate(Sample_id = Specimen_id)
endo_herb <- read.csv("C:/Users/malpa/OneDrive/Documents/Endophyte_herbarium_urbanization/endo_herb_georef.csv")



endo_herb_info <- left_join(endo_herb_georef, herb_info, by = "Sample_id")

aghybrit <- length(which(endo_herb_info$Herb_code=="BRIT" & endo_herb_info$Spp_code.x == "AGHY"))
elvibrit <- length(which(endo_herb_info$Herb_code=="BRIT" & endo_herb_info$Spp_code.x == "ELVI"))
agpebrit <- length(which(endo_herb_info$Herb_code=="BRIT" & endo_herb_info$Spp_code.x == "AGPE"))


aghytex <- length(which(endo_herb_info$Herb_code=="LL" & endo_herb_info$Spp_code.x == "AGHY"))
elvitex <- length(which(endo_herb_info$Herb_code=="LL" & endo_herb_info$Spp_code.x == "ELVI"))
agpetex <- length(which(endo_herb_info$Herb_code=="LL" & endo_herb_info$Spp_code.x == "AGPE"))

aghyokl <- length(which(endo_herb_info$Herb_code=="OKL" & endo_herb_info$Spp_code.x == "AGHY"))
elviokl <- length(which(endo_herb_info$Herb_code=="OKL" & endo_herb_info$Spp_code.x == "ELVI"))
agpeokl <- length(which(endo_herb_info$Herb_code=="OKL" & endo_herb_info$Spp_code.x == "AGPE"))

aghyokla <- length(which(endo_herb_info$Herb_code=="OKLA" & endo_herb_info$Spp_code.x == "AGHY"))
elviokla <- length(which(endo_herb_info$Herb_code=="OKLA" & endo_herb_info$Spp_code.x == "ELVI"))
agpeokla <- length(which(endo_herb_info$Herb_code=="OKLA" & endo_herb_info$Spp_code.x == "AGPE"))

aghymo <- length(which(endo_herb_info$Herb_code=="MO" & endo_herb_info$Spp_code.x == "AGHY"))
elvimo <- length(which(endo_herb_info$Herb_code=="MO" & endo_herb_info$Spp_code.x == "ELVI"))
agpemo <- length(which(endo_herb_info$Herb_code=="MO" & endo_herb_info$Spp_code.x == "AGPE"))

aghykanu <- length(which(endo_herb_info$Herb_code=="KANU" & endo_herb_info$Spp_code.x == "AGHY"))
elvikanu <- length(which(endo_herb_info$Herb_code=="KANU" & endo_herb_info$Spp_code.x == "ELVI"))
agpekanu <- length(which(endo_herb_info$Herb_code=="KANU" & endo_herb_info$Spp_code.x == "AGPE"))

aghylsu <- length(which(endo_herb_info$Herb_code=="LSU" & endo_herb_info$Spp_code.x == "AGHY"))
elvilsu <- length(which(endo_herb_info$Herb_code=="LSU" & endo_herb_info$Spp_code.x == "ELVI"))
agpelsu <- length(which(endo_herb_info$Herb_code=="LSU" & endo_herb_info$Spp_code.x == "AGPE"))

aghyam <- length(which(endo_herb_info$Herb_code=="AM" & endo_herb_info$Spp_code.x == "AGHY"))
elviam <- length(which(endo_herb_info$Herb_code=="AM" & endo_herb_info$Spp_code.x == "ELVI"))
agpeam <- length(which(endo_herb_info$Herb_code=="AM" & endo_herb_info$Spp_code.x == "AGPE"))

aghymerc <- length(which(endo_herb_info$Herb_code=="MERCA" & endo_herb_info$Spp_code.x == "AGHY"))
elvimerc <- length(which(endo_herb_info$Herb_code=="MERCA" & endo_herb_info$Spp_code.x == "ELVI"))
agpemerc <- length(which(endo_herb_info$Herb_code=="MERCA" & endo_herb_info$Spp_code.x == "AGPE"))


