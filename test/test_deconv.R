
# --------------
x=c(rgamma(10000,9),rgamma(10000,30,1))/10+3
#hist(x,breaks = 30)
E=rgamma(length(x),4)
E<- E/mean(E)
#hist(E,breaks = 30)
y=rpois(length(x),lambda = x*E)
library(data.table)
fwrite(data.table(x=x,E=E,y=y), "deconv.csv")

#hist(y,breaks = 30)

f <-  deconv2(tau=seq(0,10,0.1)[-1],X=E,y=y, family ="PoissonExp", n=200,ignoreZero = FALSE,c0=10,pDegree=10)
f <-  deconv2(tau=seq(0,10,0.1)[-1],X=E,y=y, family ="PoissonExp", n=200,ignoreZero = FALSE,c0=1, aStart = f$mle,pDegree=10)
ggplot(data = as.data.table(f$stats)) +
  geom_smooth(stat="identity", mapping = aes(x = theta, y = g, ymin = g-SE.g, ymax =g+ SE.g )) +
  labs(x = expression(theta), y = expression(hat(g)))
