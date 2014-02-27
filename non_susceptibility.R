# Script for classifying OPIEs as resistant or sensitive for thorax2
# Assumes worst-case; for an OPIE, if any sample is non-susceptible for an antibiotic, 
#   then the OPIE as a whole is considered resistant. 

# general approach is to mark records as tested against ab, then as separate column, as non-susceptible to ab. 

# Meticillin ####
abx <- subset(df, select = c(opie.id, generic.antimicrobial.name))
#head(abx)
abx$met.tested <- 0
abx$met.tested[abx$generic.antimicrobial.name == "METHICILLIN"] <- 1
abx <- subset(abx, met.tested==1)
abx <- subset(abx, !duplicated(abx$opie.id))
abx <- subset(abx, select = c(opie.id, met.tested))
#head(abx)
df <- join(df, abx, type = "left", match = "all")

abx <- subset(df, 
              generic.antimicrobial.name == "METHICILLIN" & susceptibility.result.description != "" & susceptibility.result.description != "SUSCEPTIBLE",
              select = c(opie.id, generic.antimicrobial.name, susceptibility.result.description))
#head(abx)
table(abx$generic.antimicrobial.name, useNA = "ifany")
table(abx$susceptibility.result.description, useNA = "ifany")
abx$met.resistant <- 1
abx <- subset(abx, !duplicated(opie.id), select= c(opie.id, met.resistant))
df <- join(df, abx, type = "left", match = "all")
#head(df)

#head(df[!is.na(df$met.resistant), ])
#df[df$opie.id=="10ZVE034104/01", ]
df$met.tested[df$uniq!=1] <- 0
df$met.resistant[df$uniq!=1] <- 0

# Tetracycline ####
df$cycline <- grepl("CYCLINE", df$antimicrobial.name)
abx <- subset(df, cycline=="TRUE", select = c(opie.id, generic.antimicrobial.name))
#head(abx)
abx <- subset(abx, !duplicated(opie.id))
abx$tet.tested <- 1
abx <- abx[,c(1,3)]
df <- join(df, abx, type = "left", match = "all")

abx <- subset(df, 
              generic.antimicrobial.name == "TETRACYCLINE" & susceptibility.result.description != "" & susceptibility.result.description != "SUSCEPTIBLE",
              select = c(opie.id, generic.antimicrobial.name, susceptibility.result.description))
#head(abx)
abx$tet.resistant <- 1
abx <- subset(abx, !duplicated(opie.id), select= c(opie.id, tet.resistant))
df <- join(df, abx, type = "left", match = "all")
#head(df)
df$tet.tested[df$uniq!=1] <- 0
df$tet.resistant[df$uniq!=1] <- 0

# co-amox ####
abx <- subset(df, generic.antimicrobial.name=="AMOXYCILLIN/CLAVULANATE", select = c(opie.id, generic.antimicrobial.name))
#head(abx)
abx <- subset(abx, !duplicated(opie.id))
abx$amoxclav.tested <- 1
abx <- abx[,c(1,3)]
df <- join(df, abx, type = "left", match = "all")

abx <- subset(df, 
              generic.antimicrobial.name == "AMOXYCILLIN/CLAVULANATE" & susceptibility.result.description != "" & susceptibility.result.description != "SUSCEPTIBLE",
              select = c(opie.id, generic.antimicrobial.name, susceptibility.result.description))
#head(abx)
abx$amoxclav.resistant <- 1
abx <- subset(abx, !duplicated(opie.id), select= c(opie.id, amoxclav.resistant))
df <- join(df, abx, type = "left", match = "all")
#head(df)
df$amoxclav.tested[df$uniq!=1] <- 0
df$amoxclav.resistant[df$uniq!=1] <- 0

# cipro ####
abx <- subset(df, generic.antimicrobial.name=="CIPROFLOXACIN", select = c(opie.id, generic.antimicrobial.name))
#head(abx)
abx <- subset(abx, !duplicated(opie.id))
abx$cip.tested <- 1
abx <- abx[,c(1,3)]
df <- join(df, abx, type = "left", match = "all")

abx <- subset(df, 
              generic.antimicrobial.name == "CIPROFLOXACIN" & susceptibility.result.description != "" & susceptibility.result.description != "SUSCEPTIBLE",
              select = c(opie.id, generic.antimicrobial.name, susceptibility.result.description))
#head(abx)
abx$cip.resistant <- 1
abx <- subset(abx, !duplicated(opie.id), select= c(opie.id, cip.resistant))
df <- join(df, abx, type = "left", match = "all")
#head(df)
df$cip.tested[df$uniq!=1] <- 0
df$cip.resistant[df$uniq!=1] <- 0

# ampamox ####
abx <- subset(df, generic.antimicrobial.name=="AMPICILLIN/AMOXYCILLIN", select = c(opie.id, generic.antimicrobial.name))
#head(abx)
abx <- subset(abx, !duplicated(opie.id))
abx$ampamox.tested <- 1
abx <- abx[,c(1,3)]
df <- join(df, abx, type = "left", match = "all")

abx <- subset(df, 
              generic.antimicrobial.name == "AMPICILLIN/AMOXYCILLIN" & susceptibility.result.description != "" & susceptibility.result.description != "SUSCEPTIBLE",
              select = c(opie.id, generic.antimicrobial.name, susceptibility.result.description))
#head(abx)
abx$ampamox.resistant <- 1
abx <- subset(abx, !duplicated(opie.id), select= c(opie.id, ampamox.resistant))
df <- join(df, abx, type = "left", match = "all")
#head(df)
df$ampamox.tested[df$uniq!=1] <- 0
df$ampamox.resistant[df$uniq!=1] <- 0

# azi ####
abx <- subset(df, generic.antimicrobial.name=="AZITHROMYCIN", select = c(opie.id, generic.antimicrobial.name))
#head(abx)
abx <- subset(abx, !duplicated(opie.id))
abx$azi.tested <- 1
abx <- abx[,c(1,3)]
df <- join(df, abx, type = "left", match = "all")

abx <- subset(df, 
              generic.antimicrobial.name == "AZITHROMYCIN" & susceptibility.result.description != "" & susceptibility.result.description != "SUSCEPTIBLE",
              select = c(opie.id, generic.antimicrobial.name, susceptibility.result.description))
head(abx)
abx$azi.resistant <- 1
abx <- subset(abx, !duplicated(opie.id), select= c(opie.id, azi.resistant))
df <- join(df, abx, type = "left", match = "all")
#head(df)
df$azi.tested[df$uniq!=1] <- 0
df$azi.resistant[df$uniq!=1] <- 0

# ery ####
abx <- subset(df, generic.antimicrobial.name=="ERYTHROMYCIN", select = c(opie.id, generic.antimicrobial.name))
#head(abx)
abx <- subset(abx, !duplicated(opie.id))
abx$ery.tested <- 1
abx <- abx[,c(1,3)]
df <- join(df, abx, type = "left", match = "all")

abx <- subset(df, 
              generic.antimicrobial.name == "ERYTHROMYCIN" & susceptibility.result.description != "" & susceptibility.result.description != "SUSCEPTIBLE",
              select = c(opie.id, generic.antimicrobial.name, susceptibility.result.description))
#head(abx)
abx$ery.resistant <- 1
abx <- subset(abx, !duplicated(opie.id), select= c(opie.id, ery.resistant))
df <- join(df, abx, type = "left", match = "all")
#head(df)
df$ery.tested[df$uniq!=1] <- 0
df$ery.resistant[df$uniq!=1] <- 0

# pen ####
abx <- subset(df, generic.antimicrobial.name=="PENICILLIN", select = c(opie.id, generic.antimicrobial.name))
#head(abx)
abx <- subset(abx, !duplicated(opie.id))
abx$pen.tested <- 1
abx <- abx[,c(1,3)]
df <- join(df, abx, type = "left", match = "all")

abx <- subset(df, 
              generic.antimicrobial.name == "PENICILLIN" & susceptibility.result.description != "" & susceptibility.result.description != "SUSCEPTIBLE",
              select = c(opie.id, generic.antimicrobial.name, susceptibility.result.description))
#head(abx)
abx$pen.resistant <- 1
abx <- subset(abx, !duplicated(opie.id), select= c(opie.id, pen.resistant))
df <- join(df, abx, type = "left", match = "all")
#head(df)
df$pen.tested[df$uniq!=1] <- 0
df$pen.resistant[df$uniq!=1] <- 0

# cla ####
abx <- subset(df, generic.antimicrobial.name=="CLARITHROMYCIN", select = c(opie.id, generic.antimicrobial.name))
#head(abx)
abx <- subset(abx, !duplicated(opie.id))
abx$cla.tested <- 1
abx <- abx[,c(1,3)]
df <- join(df, abx, type = "left", match = "all")

abx <- subset(df, 
              generic.antimicrobial.name == "CLARITHROMYCIN" & susceptibility.result.description != "SUSCEPTIBLE" & susceptibility.result.description != "",
              select = c(opie.id, generic.antimicrobial.name, susceptibility.result.description))
#head(abx)
abx$cla.resistant <- 1
abx <- subset(abx, !duplicated(opie.id), select= c(opie.id, cla.resistant))
df <- join(df, abx, type = "left", match = "all")
#head(df)
df$cla.tested[df$uniq!=1] <- 0
df$cla.resistant[df$uniq!=1] <- 0

# macrolides ####

abx <- subset(df, generic.antimicrobial.name=="CLARITHROMYCIN" | generic.antimicrobial.name=="AZITHROMYCIN" | 
                generic.antimicrobial.name=="ERYTHROMYCIN", select = c(opie.id, generic.antimicrobial.name))
#head(abx)
abx <- subset(abx, !duplicated(opie.id))
abx$mac.tested <- 1
sum(abx$mac.tested)
abx <- abx[,c(1,3)]
df <- join(df, abx, type = "left", match = "all")

abx <- subset(df, generic.antimicrobial.name=="CLARITHROMYCIN" | generic.antimicrobial.name=="AZITHROMYCIN" | 
                generic.antimicrobial.name=="ERYTHROMYCIN", select = c(opie.id, generic.antimicrobial.name, susceptibility.result.description))

abx <- subset(abx, 
              susceptibility.result.description != "SUSCEPTIBLE" & susceptibility.result.description != "",
              select = c(opie.id, generic.antimicrobial.name, susceptibility.result.description))
#head(abx)
abx$mac.resistant <- 1
# abx$uniq <- 0
# abx$uniq[!duplicated(abx$opie.id)] <- 1
abx <- subset(abx, !duplicated(opie.id), select= c(opie.id, mac.resistant))
df <- join(df, abx, type = "left", match = "all")
#head(df)
df$mac.tested[df$uniq!=1] <- 0
df$mac.resistant[df$uniq!=1] <- 0
rm(abx)

# doxy ####

abx <- subset(df, antimicrobial.name=="DOXYCYCLINE" , select = c(opie.id, antimicrobial.name))
#head(abx)
abx <- subset(abx, !duplicated(opie.id))
abx$dox.tested <- 1
sum(abx$dox.tested)
abx <- abx[,c(1,3)]
df <- join(df, abx, type = "left", match = "all")

abx <- subset(df, antimicrobial.name=="DOXYCYCLINE", select = c(opie.id, antimicrobial.name, susceptibility.result.description))

abx <- subset(abx, 
              susceptibility.result.description != "SUSCEPTIBLE" & susceptibility.result.description != "",
              select = c(opie.id, antimicrobial.name, susceptibility.result.description))
#head(abx)
abx$dox.resistant <- 1
# abx$uniq <- 0
# abx$uniq[!duplicated(abx$opie.id)] <- 1
abx <- subset(abx, !duplicated(opie.id), select= c(opie.id, dox.resistant))
df <- join(df, abx, type = "left", match = "all")
#head(df)
df$dox.tested[df$uniq!=1] <- 0
df$dox.resistant[df$uniq!=1] <- 0
rm(abx)
gc()
