# age structure figure

df2 <- subset(df, uniq==1)
df2$organism.name <- simpleCap(df2$organism.name)

levels(factor(df2$age.group))
df2$age.group[df2$age.group == "75+ years"] <- "> 75 years"
Encoding(df2$age.group) <- "UTF-8"
df2$age.group <- factor(df2$age.group, levels = c("<1 month", "1-11 months", "1-4 years", "5-9 years", 
                                                  "10-14 years", "15-44 years", "45-64 years", "65-74 years", 
                                                  "> 75 years", "Unknown"))
df2$sex[df2$sex == "M"] <- "Male"
df2$sex[df2$sex == "F"] <- "Female"
df2$sex[df2$sex == "U"] <- "Unknown"

p <- ggplot(df2, aes(x = age.group)) + geom_histogram(aes(fill = sex)) + facet_wrap(~organism.name, ncol = 3, scales = "free_y") + 
  theme_grey(base_size = 14) + theme(strip.text.x = element_text(face = 'italic'), 
                                    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  scale_x_discrete("Age group", labels = c("<1 month", "1-11 months", "1-4 years", "5-9 years", 
                                           "10-14 years", "15-44 years", "45-64 years", "65-74 years", 
                                           expression(phantom(x) >=75), "Unknown")) + 
  scale_y_continuous("n. isolates") + 
  scale_fill_brewer("Sex", type = "qual", palette = 3)

png("F:\\antimicrobial resistance\\Simon\\thorax2\\FES2014\\presentation\\age_fig.png", 
    width = 3060, height = 1530, res = 300)
p
dev.off()

age.trend <- ddply(df2, .(organism.name, age.group, yrqtr), summarise, n = sum(uniq))
age.trend$yrqtr <- factor(age.trend$yrqtr, 
                          levels = c("2008-1", "2008-2", "2008-3", "2008-4", "2009-1", "2009-2", "2009-3", "2009-4", 
                                     "2010-1", "2010-2", "2010-3", "2010-4", "2011-1", "2011-2", "2011-3", "2011-4", 
                                     "2012-1", "2012-2", "2012-3", "2012-4", "2013-1", "2013-2"))
head(age.trend)
# age.trend <- ddply(df2, .(organism.name, age.group, year), summarise, n = sum(uniq))
# age.trend <- age.trend[age.trend$year != 2013, ]

age.trend$yrqtr2 <- as.numeric(age.trend$yrqtr)
trend.lm <- lm(n ~ yrqtr2 * factor(organism.name) + age.group, data = age.trend )
age.trend$predicted <- predict(trend.lm, type = "response")
age.trend$se <- predict(trend.lm, type = "response", se.fit = TRUE)$se.fit
age.trend$lci <- age.trend$predicted - 1.96*age.trend$se
age.trend$uci <- age.trend$predicted + 1.96*age.trend$se

cols <- c("Trend" = "#000000", "95 % CI" = "#000000")

q <- ggplot(age.trend, aes(x = yrqtr, y = n, group = age.group)) + geom_line(aes(colour = age.group)) + 
  facet_wrap(~organism.name, ncol = 3, scales = "free_y") + geom_line(aes(y = predicted), colour = "Black") + 
  geom_line(aes(y = lci, colour = "95 % CI"), linetype = 2) + geom_line(aes(y = uci, colour = "95 % CI"), linetype = 2)
q
