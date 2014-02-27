# Figure 2:  trend in resistances by organism

table(df$ampamox.tested, useNA = "ifany")
table(df$uniq, df$ampamox.tested, useNA = "ifany", dnn = c("uniq", "ampamox tested"))
table(df$uniq, df$cla.tested, useNA = "ifany", dnn = c("uniq", "cla tested"))
table(df$uniq, df$dox.tested, useNA = "ifany", dnn = c("uniq", "dox tested"))

test.trend <- ddply(df, .(year, organism.name), summarise, n = sum(uniq, na.rm = TRUE),
                  ampamox.tested = sum(ampamox.tested, na.rm = TRUE), 
                  cla.tested = sum(cla.tested, na.rm = TRUE),
                  dox.tested = sum(dox.tested, na.rm = TRUE)
                  #                  rec_cef.tested = sum(rec.cef.tested, na.rm = TRUE), 
                  #                  rec_cef.resistant = sum(rec.cef.resistant, na.rm = TRUE)
)

head(test.trend)

m.test.trend <- melt(test.trend, id.vars = c("year", "organism.name", "n")) # edited to add n as id.var
m.test.trend$abx <- sub("\\.\\w+$","", m.test.trend$variable)
# don't need to extract tested or resistant any more. 
# m.test.trend$tested <- grepl("tested", m.test.trend$variable)
# m.test.trend$tested[m.test.trend$tested == "TRUE"] <- "tested"
# m.test.trend$resistant <- grepl("resistant", m.test.trend$variable)
# m.test.trend$resistant[m.test.trend$resistant == "TRUE"] <- "resistant"
m.test.trend$variable <- NULL
m.test.trend$pc.tested <- round(binom.confint(m.test.trend$value, m.test.trend$n, methods = "exact")$mean*100, 2)
m.test.trend$lci <- round(binom.confint(m.test.trend$value, m.test.trend$n, methods = "exact")$lower*100, 2)
m.test.trend$uci <- round(binom.confint(m.test.trend$value, m.test.trend$n, methods = "exact")$upper*100, 2)
m.test.trend$organism.name <- simpleCap(m.test.trend$organism.name)

head(m.test.trend)

m.test.trend$Antibiotic[m.test.trend$abx == "ampamox"] <- "Ampicillin/\namoxicillin"
m.test.trend$Antibiotic[m.test.trend$abx == "cla"] <- "Clarithromycin"
m.test.trend$Antibiotic[m.test.trend$abx == "dox"] <- "Doxycycline"

names(m.test.trend)[1] <- "Year"
names(m.test.trend)[2] <- "Organism"

m.test.trend <- arrange(m.test.trend, Year, Organism)

u <- ggplot(m.test.trend, aes(x = Year, y = pc.tested, group = Antibiotic)) + facet_wrap(~Organism) + 
  geom_line(aes(colour = Antibiotic)) + 
  geom_errorbar(aes(ymax = uci, ymin = lci, colour = Antibiotic), width = 0.25)
u

png("F:\\antimicrobial resistance\\Simon\\thorax2\\FES2014\\presentation\\trend_testing.png", 
    width = 3060, height = 2295, res = 300)
u + theme_grey(base_size = 16) + 
  theme(strip.text.x = element_text(face = 'italic'), legend.position = c(0.84,0.86), 
        legend.background = element_rect(fill = "#ffffffaa", colour = NA), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_y_continuous("Per cent tested", limits = c(0,100))
dev.off()

rm(m.test.trend, u, test.trend)
