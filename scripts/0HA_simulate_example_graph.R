library(pertInv)
library(GGally)

n <- 200
dts <- list()

e <- 1
X5 <- rnorm(n)
X4 <-  (2*X5 + rnorm(n)) / sqrt(5)
X2 <- rnorm(n)
Y <- (-2*X2 + 2*X4 + rnorm(n)) /sqrt(9)
X3 <- (2*Y + rnorm(n)) / sqrt(5)
dts[[e]] <- data.table(e, Y, X2, X3, X4, X5)

e <- 2
X5 <- rnorm(n)
X4 <-  (2*X5 + rnorm(n)) / sqrt(5)
X2 <- rnorm(n)
Y <- (-2*X2 + 2*X4 + rnorm(n)) /sqrt(9)
X3 <- rnorm(n)
dts[[e]] <- data.table(e, Y, X2, X3, X4, X5)

e <- 3
X5 <- rnorm(n)
X4 <-  rnorm(n)
X2 <- rnorm(n)
Y <- (-2*X2 + 2*X4 + rnorm(n)) /sqrt(9)
X3 <- (2*Y + rnorm(n)) / sqrt(5)
dts[[e]] <- data.table(e, Y, X2, X3, X4, X5)

dt <- rbindlist(dts)
dt[,e:=as.factor(e)]

#pairs(dt)
ggpairs(dt, aes(colour = e, alpha = 0.4), axisLabels="none", columns=2:ncol(dt))

dt2 <- melt(dt,id.vars = c("e","Y"))[,Y_name:="Y"]

for (j in 1:3) {
  dt3 <-  dt2[as.integer(e)<=j]
  p <- ggplot(dt3, aes(y=value,x=Y))+
    geom_blank(data=dt2)+
    facet_grid(variable~Y_name, switch="both") +
    geom_point(aes(fill=e,color=NA),shape=21, stroke=0, alpha=0.2) +
    geom_smooth(aes(color=e),method="lm") + theme_void() +
    theme(legend.justification = c(1, 1),
          legend.position = c(1, 1),
          strip.background = element_rect(fill="gray")
          )
  print(figure("ICP-Visualization: pairs plot", p, sub_title=paste0("e<=",j),width=2, height=8))
}

dt2[,e_str:=sprintf("e=%i",e)]
ggplot(dt2, aes(y=value,x=Y))+
  geom_blank(data=dt2)+
  facet_grid(variable~e_str, switch="both") +
  geom_point() +
  geom_smooth(method="lm",color="red") + theme_void() +
  theme(legend.justification = c(1, 1),
        legend.position = c(1, 1),
        strip.background = element_rect(fill="gray",color=NA)
  )

