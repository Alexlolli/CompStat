##### Microbenchmark bandwidth selection #####



myKernel_ep <- myKernel(K = expression(3 / 4 * (1 - x^2)), 
                        rng = c(-1, 1))

bandwidth_bench <- microbenchmark(bandwidth.myKernel(obj = myKernel_ep, 
                                                     data = log(F12), 
                                                     method = "silverman"),
                                  bandwidth.myKernel(obj = myKernel_ep, 
                                                     data = log(F12), 
                                                     method = "plug-in",
                                                     pilot = "gauss"),
                                  bandwidth.myKernel(obj = myKernel_ep, 
                                                     data = log(F12), 
                                                     method = "plug-in",
                                                     pilot = "ep"),
                                  bandwidth.myKernel(obj = myKernel_ep, 
                                                     data = log(F12), 
                                                     method = "cv"))

levels(bandwidth_bench$expr) <- c("silverman", "plug-in gaussian", 
                                  "plug-in Epanechnikov", "LOOCV")


autoplot(bandwidth_bench) + 
  geom_jitter(position = position_jitter(0.2, 0), 
              aes(color = expr), alpha = 0.4) + 
  aes(fill = I("gray")) +
  theme(legend.position = "none")


set.seed(1234)
Nsim <- 2^20
par_1 <- c( 0, 1)
par_2 <- c( 6, 2)

MM_var<-rbinom(Nsim, size = 1,p=0.5) #Bernouilli variables

x <- (1-MM_var)* rnorm(Nsim, mean = par_1[1], sd = par_1[2]) +
  MM_var * rnorm(Nsim, mean = par_2[1], sd = par_2[2])

conf <- expand.grid(
  n = 2^(5:11),
  method = c("silverman", "plug-in\", pilot = \"ep", 
             "plug-in\", pilot = \"gauss", "cv")
)

calls <- paste0("bandwidth.myKernel(myKernel_ep, x[1:", conf[ , 1], "], 
                method = \"", conf[ , 2], "\")")
expr_list <- lapply(calls, function(x) parse(text = x)[[1]])

bandwidth_bench_n <- microbenchmark(list = expr_list, times = 40L)

levels(conf$method)<-c(levels(conf$method),"Plug_in_epan","Plug_in_gauss")

conf$method[conf$method=='plug-in", pilot = "ep']<-"Plug_in_epan"
conf$method[conf$method=='plug-in", pilot = "gauss']<-"Plug_in_gauss"

bench_data <- as.data.frame(cbind(conf, time = summary(bandwidth_bench_n)$median))

ggplot(bench_data, aes(x = n, y = time, color = method)) + 
  geom_abline(intercept = 7.5, slope = 1, color = "gray", linetype = 2) +
  stat_summary(fun.y = "median", geom = "line") + 
  stat_summary(fun.y = "median", geom = "point") + 
  facet_wrap(~ "Bandwidth calculation time") + 
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous("Time (s)", trans = "log2", 
                     breaks = c(1e3, 1e4, 1e5, 1e6, 1e7, 1e8), 
                     labels = c("0.001", "0.01", "0.1", "1", "10", "100")) +
  scale_color_discrete("Method:") + 
  theme(legend.position="top")


###### Profiling ######

library("profvis")


getwd()
setwd('~/Dropbox/GitHub/CompStat/Assignment_1/Assignment-1-profiling-debbie')
source('~/Dropbox/GitHub/CompStat/Assignment_1/Assignment-1-profiling-debbie/Assignment-1-profiling-debbie.R', echo=TRUE)

x <- rnorm(2^10)

profFile <- "prof.html"
Rprof(profFile)
{bandwidth(x, method = "plug-in", pilot = "gauss")
  bandwidth(x, method = "plug-in", pilot = "ep")
  bandwidth(x, method = "silverman")
  bandwidth(x, method = "cv")}
Rprof(NULL)

summaryRprof("prof.Rprof")


profvis({bandwidth(x, method = "plug-in", pilot = "gauss")
  bandwidth(x, method = "plug-in", pilot = "ep")
  bandwidth(x, method = "silverman")
  bandwidth(x, method = "cv")}, interval = 0.005)
