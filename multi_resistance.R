# trends in resistance to any, and trends in resistance to all. 

df$any.t <- 0
df$any.t[df$ampamox.tested == 1 | df$cla.tested == 1 | df$dox.tested == 1] <- 1
df$any.t[df$uniq == 0] <- 0

df$any.r <- 0
df$any.r[df$ampamox.resistant == 1 | df$cla.resistant == 1 | df$dox.resistant == 1] <- 1
df$any.r[df$uniq == 0] <- 0

table(df$any.t, df$any.r, dnn = c("tested", "resistant"))

df$all.t <- 0
df$all.t[df$ampamox.tested == 1 & df$cla.tested == 1 & df$dox.tested == 1] <- 1
df$all.t[df$uniq == 0] <- 0

df$all.r <- 0
df$all.r[df$ampamox.resistant == 1 & df$cla.resistant == 1 & df$dox.resistant == 1] <- 1
df$all.r[df$uniq == 0] <- 0

table(df$all.t, df$all.r, dnn = c("tested", "resistant"))

df$ampmac.tested <- 0
df$ampmac.tested[df$ampamox.tested == 1 & df$mac.tested == 1 & df$uniq == 1] <- 1
table(df$ampmac.tested[df$uniq == 1])

df$ampmac.resistant <- 0
df$ampmac.resistant[df$ampamox.resistant == 1 & df$mac.resistant == 1 & df$uniq == 1] <- 1
table(df$ampmac.resistant[df$uniq == 1 & df$ampmac.tested == 1])

# trend in resistance to front line and second line. 
df$all.r <- NULL
df$ant.r <- NULL
df$ant.t <- NULL
df$all.t <- NULL
gc()

trend_multi <- ddply(df, .(year, organism.name), summarise, 
                     ampmac.tested = sum(ampmac.tested, na.rm = TRUE),
                     ampmac.resistant = sum(ampmac.resistant, na.rm = TRUE))

head(trend_multi)

trend_multi$pc.resistant <- round(binom.confint(trend_multi$ampmac.resistant, trend_multi$ampmac.tested, methods = "exact")$mean*100, 2)
trend_multi$lci <- round(binom.confint(trend_multi$ampmac.resistant, trend_multi$ampmac.tested, methods = "exact")$lower*100, 2)
trend_multi$uci <- round(binom.confint(trend_multi$ampmac.resistant, trend_multi$ampmac.tested, methods = "exact")$upper*100, 2)

names(trend_multi)[1] <- "Year"
names(trend_multi)[2] <- "Organism"

trend_multi$Organism <- simpleCap(trend_multi$Organism)

o <- ggplot(trend_multi, aes(x = Year, y = pc.resistant, group = Organism)) + 
  geom_line(aes(colour = Organism)) + 
  geom_errorbar(aes(ymax = uci, ymin = lci, colour = Organism), width = 0.15) + 
  scale_y_continuous("Per cent resistant", limits = c(0, 50)) + 
  theme(legend.text = element_text(face = 'italic'),
        legend.position = c(0.95,0.7), legend.direction = "vertical", legend.justification = c(1,0))

png("F:\\antimicrobial resistance\\Simon\\thorax2\\FES2014\\presentation\\trend_multi.png", 
    width = 3060, height = 2000, res = 300)
o
dev.off()
