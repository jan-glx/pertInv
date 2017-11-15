library(pertInv)
library(GGally)

n <- 300
dts <- list()

e <- 1
X_5 <- rnorm(n)
X_4 <-  (2*X_5 + rnorm(n)) / sqrt(5)
X_2 <- rnorm(n)
Y <- (-2*X_2 + 2*X_4 + rnorm(n)) /sqrt(9)
X_3 <- (2*Y + rnorm(n)) / sqrt(5)
dts[[e]] <- data.table(e, Y, X_2, X_3, X_4, X_5)

e <- 2
X_5 <- rnorm(n)
X_4 <-  (2*X_5 + rnorm(n)) / sqrt(5)
X_2 <- rnorm(n) + 1
Y <- (-2*X_2 + 2*X_4 + rnorm(n)) /sqrt(9)
X_3 <- rnorm(n) + 1
dts[[e]] <- data.table(e, Y, X_2, X_3, X_4, X_5)

e <- 3
X_5 <- rnorm(n)
X_4 <-  rnorm(n) + 1
X_2 <- rnorm(n)
Y <- (-2*X_2 + 2*X_4 + rnorm(n)) /sqrt(9)
X_3 <- (2*Y + rnorm(n)) / sqrt(5)
dts[[e]] <- data.table(e, Y, X_2, X_3, X_4, X_5)

dt <- rbindlist(dts)
dt[,e:=as.factor(e)]

#pairs(dt)
ggpairs(dt, aes(colour = e, alpha = 0.4), axisLabels="none", columns=2:ncol(dt))
dt[,i:=.I]
dt[,environment:=sprintf("e=%i",e)]
dt2 <- melt(dt,id.vars = c("i","e","environment","Y"))[,Y_name:="Y"]

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


figure("ICP-Visualization: pairs plot -- individual environments",
  ggplot(dt2, aes(y=value,x=Y,color=environment))+
    geom_blank(data=dt2)+
    facet_grid(variable~environment, switch="both") +
    geom_point() +
    geom_smooth(method="lm", fullrange=TRUE) + theme_void() +
    theme(#legend.justification = c(1, 1),
          #legend.position = c(1, 1),
          strip.background = element_rect(fill="gray",color=NA)
    ) +   scale_color_manual(values=cbPalette[-1]),
  width=5, height=8
)

dt2[,dist:="Y"]

figure("ICP-Visualization: pairs plot--merged",
  ggplot(dt2, aes(y=value,x=Y,color=environment))+
    geom_blank(data=dt2)+
    facet_grid(variable~dist, switch="both") +
    geom_point() +
    geom_smooth(method="lm", fullrange=TRUE) + theme_void() +
    theme(#legend.justification = c(1, 1),
          #legend.position = c(1.1, 1.1),
          strip.background = element_rect(fill="gray",color=NA)
    ) +   scale_color_manual(values=cbPalette[-1]),
  width=3, height=8
)


#ggplot(dt2, aes(color=Y,x=value,y=0)) +
#  facet_grid(variable~environment)+geom_point()


#dt3 <- dt2[dt2, .(S=paste0("S={", variable, ",", i.variable, "}"),
#                  S1_name = variable, S2_name = i.variable,
#                  S1 = value, S2 = i.value,
#                  Y = Y, environment, e),
#           on = .(e, i), allow.cartesian = TRUE]


#ggplot(dt3[as.character(S1_name) < as.character(S2_name)], aes(color=Y, x=S1, y=S2)) +
#  facet_grid(S~environment) + geom_point()
#
dt2[,Y_adj:=residuals(lm(Y~value)), by=.(e,variable)]

dt3 <- dt2[dt2, .(S=i.variable,Y=i.Y,Y_adj=i.Y_adj,variable,value,environment,e),on=.(e,i),allow.cartesian=TRUE]
dt3[,dist:=sprintf("Y|%s",S)]

figure("ICP-Visualization: pairs plot--adj",
       ggplot(dt3, aes(y=value,x=Y_adj,color=environment))+
         geom_blank(data=dt3)+
         facet_grid(variable~dist, switch="both") +
         geom_point() +
         geom_smooth(method="lm", fullrange=TRUE) + theme_void() +
         theme(#legend.justification = c(1, 1),
           #legend.position = c(1.1, 1.1),
           strip.background = element_rect(fill="gray",color=NA)
         ) +   scale_color_manual(values=cbPalette[-1]),
       width=7, height=8
)
library(latex2exp)

dt2[,Y_adj:=residuals(lm(Y~value)), by=.(e,variable)]

dt2[, variable_int := as.integer(variable)]
dt3 <- dt2[dt2, .(S=paste0("S=\\{", variable, ",", i.variable, "\\}\\} "),
                  S1_name = variable, S2_name = i.variable,
                  S1 = value, S2 = i.value,
                  Y_adj_name = sprintf("Y-E[Y|%s,%s]",variable, i.variable),
                  Y = Y, environment, e),
           on = .(e, i, variable_int<variable_int),,nomatch=0, allow.cartesian = TRUE][,
             Y_adj := residuals(lm(Y~S1+S2)), by=S
             ]

dt2[,':='(S=paste0("S=\\{", variable,"\\}\\} "),
      S1_name = variable,
      S1 = value,
      Y_adj_name = sprintf("Y-E[Y|%s]",variable),
      Y_adj = residuals(lm(Y~value))), by= variable]
dt[,':='(S="S=\\{ \\}\\} ",
      Y_adj_name = "Y-E[Y]",
      Y_adj = residuals(lm(Y~1)))]


figure("ICP-Visualization: |S|=2",
       ggplot(dt3, aes(y=Y_adj, x="", color=environment))+
         facet_grid(environment~latex2exp::TeX(S), scales="free_y") +
         geom_jitter_normal(alpha=0.5, stroke=0, width=0.2) +
         geom_boxplot(fill=NA, color="black") +
         ylab(expression(Y-hat(Y)[S])) + geom_abline(slope=0, linetype="dashed",color=cbPalette[1]) +
         theme(axis.text.x=element_blank(),
               axis.title.y=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.line.y= element_blank(),
               #legend.justification = c(1, 1),
           #legend.position = c(1.1, 1.1),
           legend.position="none",
           strip.background = element_rect(fill="gray",color=NA)
         ) +   scale_color_manual(values=cbPalette[-1]) + coord_flip() ,
       width=8, height=3
)


dt4 <- rbind(dt,dt2,dt3,fill=TRUE)
dt4[, S:= factor(S, levels= unique(S), labels=sapply(unique(S), latex2exp::TeX, output="character"))]
dt4[, environment:= factor(environment, levels= unique(environment), labels=sapply(unique(environment), latex2exp::TeX, output="character"))]



figure("ICP-Visualization: |S|<=2",
       ggplot(dt4, aes(y=Y_adj, x="", color=environment))+
         facet_grid(environment~S, scales="free_y",labeller = label_parsed) +
         geom_jitter_normal(alpha=0.5, stroke=0, width=0.2) +
         geom_boxplot(fill=NA, color="black") +
         ylab(expression(Y-hat(Y)[S])) + geom_abline(slope=0, linetype="dashed",color=cbPalette[1]) +
         theme(axis.text.x=element_blank(),
               axis.title.y=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               axis.line.y= element_blank(),
               #legend.justification = c(1, 1),
               #legend.position = c(1.1, 1.1),
               legend.position="none",
               strip.background = element_rect(fill="gray",color=NA)
         ) +   scale_color_manual(values=cbPalette[-1]) + coord_flip() ,
       width=10, height=3
)
