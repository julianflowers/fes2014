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

# regression
ampmac <- subset(df, ampmac.tested == 1 & uniq == 1)
# complete cases for organism, quarter, country, specimen source, age2 and male sex. 
ampmac <- ampmac[complete.cases(ampmac[,c(2, 18, 14, 10, 46, 47)]), ]

ampmac.m <- glm(ampmac.resistant ~ year * organism.name + age2 + factor(quarter) + postcode.derived.region.name +
                  specimen.source.type.description + male, data = ampmac, family = "poisson")

ampmac.m2 <- glm(ampmac.resistant ~ year * organism.name * age2 + factor(quarter) + postcode.derived.region.name +
  specimen.source.type.description + male, data = ampmac, family = "poisson")

epicalc::lrtest(ampmac.m, ampmac.m2) # No evidence for interaction. 

# functions for extracting IRRs ####
orCalc <- function(outcome, riskf){
  # x should be outcome, y stratifying variable. 
  # cribbed from mhodds in epicalc package
  # with reference to p157 + p164 of Kirkewood and Sterne, Essential Medical Statistics. 2nd Ed
  tab <- table(riskf, outcome)
  #print(tab)
  or <- c(1:dim(tab)[1]) # create vector of same length as table rows. 
  se.log.or <- c(1:dim(tab)[1])
  lci <- c(1:dim(tab)[1])
  uci <- c(1:dim(tab)[1])
  z <- c(1:dim(tab)[1])
  p <- c(1:dim(tab)[1])
  for (i in 1:dim(tab)[1]){
    or[i] <- (tab[1,1]*tab[i,2])/(tab[1,2]*tab[i,1])
    se.log.or[i] <-sqrt(1/tab[1,1] + 1/tab[1,2] + 1/tab[i,1] + 1/tab[i,2])
    lci[i] <- or[i]/exp(1.96 * se.log.or[i])
    uci[i] <- or[i]*exp(1.96 * se.log.or[i])
    z[i] <- log(or[i])/se.log.or[i]
    p[i] <- 2*(1-pnorm(abs(z[i])))
  }
  m <- as.matrix(cbind(or, se.log.or, lci, uci, z, p))
  return(m)
}

lincom <- function(svycontrast_object){
  require(survey)
  if (class(svycontrast_object)=="svystat"){
    or <- exp(svycontrast_object[1])
    lci <- or / exp(1.96*sqrt(attributes(svycontrast_object)$var))
    uci <- or * exp(1.96*sqrt(attributes(svycontrast_object)$var))
    return(as.list(c(or, lci, uci)))
  } else {
    print("Requires object of class svystat")
  }
}

sehac<-function(fit,vcov=sandwich){ #Convenience function for robust standard errors
  coeftest(fit,vcov)
}

funinteff<-function(mod,var,vcov=sandwich){ #mod is an lm() object, var is the name of the main effect that was interacted, vcov is the type of variance covariance method you want to use 
  #Extract Coefficient names create 'beta names' to feed to deltaMethod()
  cnames<-coef(mod)
  pnams<-data.frame('b'=paste('b',0:(length(cnames)-1),sep=""),'est'=cnames) #assign parameter names so that deltaMethod does not throw an error
  
  #Extract the specific parameters of interest
  vars<-grep(var,names(cnames),value=T)
  var1<-vars[1]
  intvars<-vars[2:length(vars)]
  bi<-pnams[var1,'b']
  
  #--Create Data Frame to store Main Effect
  int<-sehac(mod,vcov=vcov)[var1,c('Estimate','Std. Error')]
  int<-as.data.frame(t(int))
  names(int)<-c('Estimate','SE')
  row.names(int)<-var1
  
  #Loop through and store the results in a data.frame
  for(i in 1:length(intvars)){
    bint<-pnams[intvars[i],'b']
    eq<-paste(bi,bint,sep="+")
    interac<-deltaMethod(mod,eq,parameterNames=pnams[,1],vcov=vcov)
    row.names(interac)<-intvars[i]
    
    int<-rbind(int,interac)
  }
  return(int)
}


ampmac.robust <- funinteff(ampmac.m, "year")
ampmac.robust$z <- ampmac.robust$Estimate/ampmac.robust$SE
ampmac.robust$p <- 2*pnorm(-abs(ampmac.robust$z))
ampmac.robust

ampmac.robust$irr <- sprintf("%.2f", exp(ampmac.robust$Estimate) )
ampmac.robust$ci <- paste("(", sprintf("%.2f", exp(ampmac.robust$Estimate)/exp(1.96*ampmac.robust$SE)), " - ",
                          sprintf("%.2f", exp(ampmac.robust$Estimate)*exp(1.96*ampmac.robust$SE)), ")", 
                          sep = "")
ampmac.robust$organism <- c("\\textit{H. influenzae}", "","", "", "\\textit{S. aureus}", 
                            "\\textit{S. pneumoniae}" )
ampmac.robust <- ampmac.robust[ampmac.robust$organism != "" ,]
ampmac.robust$Estimate <- NULL
ampmac.robust$SE <- NULL
ampmac.robust$z <- NULL
ampmac.robust$p <- NULL
row.names(ampmac.robust) <- c(1:3)
ampmac.robust <- ampmac.robust[,c(3,1,2)]
names(ampmac.robust) <- c("Organism", "Incident rate ratio", "95 \\% confidence interval")
out <- xtable(ampmac.robust, caption = "")
sink("multiresistance_table.txt")
print(out, include.rownames = FALSE, booktabs = TRUE, sanitize.text.function = function(x){x}, )
sink()
