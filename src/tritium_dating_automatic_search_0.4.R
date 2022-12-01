
#########################
########## script with which
########## to fit a lumped-parameter model
########## to tritium and stable isotope measurements (monthly values)
#########################

rm(list = ls())

##### written in August 2022 by Julien Farlin

# wk.dir=c("C:/Users/Lutz0003/OneDrive - Universiteit Utrecht/WATSON/modell_vergleich/h3_h2_combined_weierbach/analysis")

data.raw.dir <- c("./data/raw")

data.proces.dir <- c("./data/processed")

graph.dir <- c("./results/figures")

mod.dir <- c("./src/modules")

output.dir <- c("./results/output")

#########################
# loading the data
#########################

getwd()
#setwd(data.raw.dir)

# tritium input
# data=read.table("tritium_Vienna_1961_2005.csv",sep=",",header=T)
data <- read.table(paste(data.raw.dir,"/tritium_trier_vienna_1961_2019.csv",sep=""), sep = ",", header = T)
# tritium input function (station Vienna). Monthly values.
h3.input <- data$H3
h3.input.date <- as.Date(data$Time, "%Y-%m-%d")

# stable isotope input

sw <- 2 #### switch, sw=1 calculates the time series for 18O, sw=2 for deuterium !!!!!! User defined.
# data=read.csv("stable_Ondrasova_1990_2021.csv",sep=",",header=TRUE)
data <- read.csv(paste(data.raw.dir,"/Weierbach_H2_rainfall_2009-2019.csv",sep=""), sep = ";", header = TRUE, na.strings = "no data")
stable.input <- data$isotope
stable.input.date <- as.Date(data$Time, "%m/%d/%Y")
# stable.input.date=as.Date(data$Time,"%Y-%m-%d")
rm(data)

# tritium output
#setwd(data.proces.dir)
h3.output <- read.table(paste(data.proces.dir,"/Weierbach_tritium_2011-2017.csv", sep=""),
                        header = T, sep = ";", na.strings = "no data") # tritium measurements at the outlet
h3.output$date <- as.Date(h3.output$date, format = "%m/%d/%Y")

# stable isotope output

stable.output <- read.table(paste(data.proces.dir,"/Weierbach_H2_streamwater_2009-2019.csv",sep=""),
                            header = T, sep = ";", na.strings = "no data") # stable isotope measurements in spring water
# stable.output$date=as.Date(stable.output$date,format="%Y-%m-%d") #example of a different date format. See "as.Date" for details.
stable.output$date <- as.Date(stable.output$date, format = "%m/%d/%Y")

obs.date <- format(stable.output$date, format = "%Y-%m") # calculates the monthly mean
obs.val <- stable.output$isotope
dummy <- aggregate(obs.val, list(obs.date), median) # calculates the median value for each month in case more than one observation is available over that time period.
obs.date <- as.character(dummy$Group.1)
stable.input.date <- as.Date(paste(obs.date, "-15", sep = ""))
stable.input <- dummy$x
rm(obs.date, obs.val)

#########################
# interpolate missing values
# of the input functions
#########################

# Tritium
x <- seq(1, length(h3.input), by = 1)
prov <- approx(x, h3.input, method = "linear", n = length(x))
h3.input <- round(prov$y)

# Stable isotope
x <- seq(1, length(stable.input), by = 1)
prov <- approx(x, stable.input, method = "linear", n = length(x))
stable.input <- round(prov$y)

rm(x, prov)

#########################
# calculates the mean deuterium value of input
#########################

input <- stable.input

# Annual mean+winter and summer means
n.year <- ceiling(length(input) / 12) # number of years with monthly observations
dummy <- c(4, 5, 6, 7, 8, 9) # summer month of a year (from April until September)
dummy2 <- c(4, 5, 6, 7, 8, 9)
for (ee in 2:n.year) # creates a vector of the summer month entries
{
  dummy <- dummy + 12
  dummy2 <- c(dummy2, dummy)
}
mean.summer.input <- mean(input[dummy2]) # mean isotopic value for the summer months
mean.winter.input <- mean(input[-dummy2]) # mean isotopic value for the winter months
mean.year.input <- mean(input) # mean isotopic value (summer AND winter months)
# mean.h2=mean(stable.output$h2,na.rm=T)
# mean.o18=mean(stable.output$o18,na.rm=T)

# moving six months input means (summer and winter)
obs.val <- stable.input

dummy1 <- c(1, 1, 1, 1, 1, 1)
dummy <- dummy1
for (ee in 2:200) # creates a vector marking each successive period of six months with an different, increasing number.
{
  dummy2 <- dummy1 * ee
  dummy <- c(dummy, dummy2)
}
rm(dummy1, dummy2)
dummy <- dummy[4:length(dummy)] # crops the "time" vector so that the first three months of the first year are marked separately from the next six (January, February and March belong to the first winter period, while April, May, etc. to the first summer period).
dummy <- dummy[1:length(obs.val)] # crops the vector to the length of the stable isotope input.

dummy2 <- aggregate(obs.val, list(dummy), mean) # calculates the median value for each six-months period.
winter.summer.stable <- dummy2$x # moving mean value of each summer and winter period of the input.

dummy1 <- as.numeric(format(stable.input.date[1], "%Y"))
dummy2 <- as.numeric(format(stable.input.date[length(stable.input.date)], "%Y"))
dummy3 <- as.Date(paste(dummy1, "-04-01", sep = ""))
dummy4 <- as.Date(paste(dummy2, "-04-01", sep = ""))
dummy5 <- as.Date(paste(dummy1, "-10-01", sep = ""))
dummy6 <- as.Date(paste(dummy2, "-10-01", sep = ""))

summer.dates <- seq(dummy3, dummy4, by = "year") # output 2
winter.dates <- seq(dummy5, dummy6, by = "year") # output 2
dummy <- c(summer.dates, winter.dates)
winter.summer.dates <- c(as.Date(paste(format(stable.input.date[1], "%Y"), "-01-01", sep = "")), dummy[order(dummy)]) # final time steps of the bi-annual winter/summer means.

rm(dummy1, dummy2, dummy3, dummy4, dummy, dummy5, dummy6)

#########################
# lengthening the input function beyond the last measurements
#########################

# aa=length(h3.input)-59 #grabs the last 5 years of the tritium time series
# bb=length(h3.input)
# dummy=rep.int(h3.input[aa:bb],12)
# h3.input=c(h3.input,dummy)
# l.obs=length(h3.input)
# t.vec=seq.Date(h3.input.date[1],length.out=l.obs,by="month")
# h3.input.date=t.vec # adds dates to equal the lengthened input vector
# rm(aa,bb,dummy,l.obs,t.vec)

#########################
# visualising the input functions
#########################

o18lab <- expression(paste("", delta^18, "", O, "[\u2030]"))
h2lab <- expression(paste("", delta^2, "", H, "[\u2030]"))
x.inf <- as.Date("01.01.1960", format = "%d.%m.%Y") # date of the lower x-axis limit
x.sup <- as.Date("01.01.2020", format = "%d.%m.%Y") # date of the upper x-axis limit

#not working right now

#setwd(graph.dir)
bitmap(paste(graph.dir,"/input_functions.png",sep=""), type = "png16m", height = 20, width = 60, res = 150, encoding = "WinAnsi.enc") # the WinAnsi encoding is necessary to use the "per mille" sign (using Linux). Might not work on a windows computer. Please find a work around, or think about changing to linux, which is an open source, community supported, efficient, secure operating system.
source(paste(mod.dir,"/visualisation_inputs.Rhistory",sep=""))
dev.off()

#########################
# parameters and variables
#########################

DT.trit <- 12.3 # tritium half-life [y]
DT.trit <- DT.trit * 12 # tritium half-life [m]
mean.pf.vec <- seq(0, 30, by = 1) # piston-flow residence time [y]. User defined. Please do not use a priori a range too narrow, as this would lead to missing out on parameter combinations yielding perfectly acceptable fits. The search for the best fit or fits is done in the "post-processing" section below.
mean.exp.vec <- seq(1, 20, by = 1) # exponential residence time [y]. User defined. The same remark as for the piston-flow component applies here too.
n.vec <- seq(0.1, 0.9, by = 0.1) # parameter controling the infiltration of rainwater in summer (modified from Grabczak et al., Catena, 1984). 1=tracer input in the summer months is equal to precipitation amount, 0=tracer input in the summer months is zero.
alt.effect <- 0 # altitude effect for deuterium (in per mille per 100 metres). Literature values might be inadequate for the study site ! If possible, calculate the local gradient from available measurement stations.

#########################
# calculating the tritium and stable isotope outputs
#########################

library(signal) # for the "conv" function

# parameters that are reinitialised at each run of the loop
params <- do.call(expand.grid, list(n.vec, mean.pf.vec, mean.exp.vec))
names(params) <- c("n", "pf", "exp") # creates a matrix with all possible parameter combinations
rtdf.out <- matrix(ncol = 1200, nrow = dim(params)[1]) # stores the calculated transfer functions. Each vector is one hundred years long (12 months times a hundred years=1200 entries)
output1 <- matrix(ncol = 1200, nrow = dim(params)[1]) # stores the calculted output for input 1 (normally, tritium).
output2 <- matrix(ncol = 1200, nrow = dim(params)[1]) # stores the calculted output for input 2 (normally, a stable isotope).
gw.iso <- vector(mode = "numeric", length = dim(params)[1]) # stores the calculated mean isotopic content in the output as a function of the parameter controlling the summer to winter infiltration ratio.

for (i in 1:dim(params)[1])
{
  ###################################
  # setting the parameters (all calculations are done on a monthly basis, hence the unit of time is also in months)
  ###################################

  print(paste(i, "/", dim(params)[1], sep = ""))
  n <- params$n[i]
  mean.pf <- params$pf[i] * 12 # piston flow component [m]
  mean.exp <- params$exp[i] * 12 # exponential component [m]

  ###################################
  # pre-processing the input function
  ###################################

  input1 <- h3.input
  input2 <- stable.input
  n.year1 <- ceiling(length(input1) / 12) # number of years with monthly observations (tritium)
  # ann.mark1=rep(seq(1,n.year1),each=12)
  # ann.mark1=ann.mark1[1:length(input1)]
  # ann.mark2=rep(seq(1,n.year2),each=12)
  # ann.mark2=ann.mark2[1:length(input2)]
  # trit.med=tapply(trit,ann.mark,mean) #annual mean isotopic ratio

  # input 1
  dummy <- c(4, 5, 6, 7, 8, 9)
  for (ee in 1:n.year1) # reduces the amplitude of the sommer months by the factor "n" (tritium)
  {
    input1[dummy] <- input1[dummy] * n
    dummy <- dummy + 12
  }
  input1 <- input1[1:length(h3.input)]
  rm(dummy, ee)
  ix <- which(is.na(input1))
  input1[ix] <- 0

  # input 2

  gw.comp <- (mean.summer.input * n + mean.winter.input) / (1 + n) # mean isotopic value in groundwater as a function of the summer to winter infiltration ratio (see Grabczak et al., Catena 1984 for details)
  ratio <- abs(gw.comp) / abs(mean.year.input)
  input2 <- input2 * ratio # shifts the input function to match the mean isotopic value calculated at the outlet from the summer to winter infiltration ratio.
  # rm(dummy,dummy2,ee,ratio)

  #########################
  # calculating the outputs
  #########################

  #setwd(mod.dir)
  # computes and stores the transfer function
  source(paste(mod.dir,"/rtdf_EPM.Rhistory",sep=""))
  rtdf.out[i, ] <- rtdf.tot[1:1200]

  # calculates the convolution and stores the output for the first input (i.e. tritium)
  input <- input1
  source(paste(mod.dir,"/convolution_EPM.Rhistory",sep=""))
  output1[i, 1:1200] <- convolution.out[1:1200] # stores the predicted output (input 1)

  # calculates the convolution and stores the output for the second input (i.e. stable isotope)
  input <- c(input2, input2) # doubles the input length in order to use the first as warm-up period. Without warm-up, the beginning of the predicted output will progressively sink to the negative values typical of deuterium or oxygen-18, which is a computational artefact. Using a warm-up period reduces the influence of the unknown input, but it is still better (if possible) to make sure that the input time period precedes by at least one mean transit time the first output measurements.
  sub1 <- length(input2) + 1
  source(paste(mod.dir,"/convolution_EPM_stable.Rhistory",sep=""))
  sub2 <- sub1 + 1199
  output2[i, 1:1200] <- convolution.out[sub1:sub2] # stores the predicted output (input 2)
  gw.iso[i] <- gw.comp # stores the predicted mean isotopic value in the output (input 2)
}

rtdf.out[which(is.na(rtdf.out))] <- 0

#################################
# post-processing
#################################

### Important note: the post-processing workflow allows the user to concentrate the search on
### a subspace of the entire parameter space defined in the section "parameters and variables".
### Think carefully before reducing a priori the search of the best fit to a narrower
### parameter combination. Using different sub-spaces (by repeating the procedure below multiple
###  times changing the range of piston-flow and exponential mean transit times, and saving graphs
### and tables under different names) might work better than selecting a single area, as different
### parameter combinations might yield equally good fits, but for largely different
### combinations of mean transit times. In that case, verify whether the water column in storage is compatible
### with storage space (T=nH/R or T=nH/Q, T=exponential mean transit time, n=total porosity, H=water column in storage,
### R=annual recharge rate and Q=mean outflow rate).

## step 1: calculation of the goodness of fit  ##

input.dates <- format(seq(h3.input.date[1], by = "month", length.out = 1200), format = "%Y-%m") # first for output 1
obs.date <- format(h3.output$date, format = "%Y-%m")
obs.val <- h3.output$tritium
ix <- which(input.dates %in% obs.date)

for (i in 1:dim(params)[1])
{
  # goodness of fit
  pred.val <- output1[i, ix]
  err <- sum(abs(pred.val - obs.val)) / length(obs.val) # mean absolute error
  params$error1[i] <- err # adds to the parameter object the prediction error associated with outpu1
}

# This section selects part of the ouput time series to be used to calculate the goodness of fit for the stable isotopes
first.date <- format(as.Date("01/01/2015", "%m/%d/%Y"), format = "%Y-%m") # User defined
last.date <- format(as.Date("01/01/2021", "%m/%d/%Y"), format = "%Y-%m") # User defined
input.dates <- format(seq(stable.input.date[1], by = "month", length.out = 1200), format = "%Y-%m") # output 2
obs.date <- format(stable.output$date, format = "%Y-%m")
obs.val <- stable.output$isotope
ix <- which(obs.date > first.date & obs.date < last.date)
obs.date <- obs.date[ix]
obs.val <- obs.val[ix]

dummy <- aggregate(obs.val, list(obs.date), median) # calculates the median value for each month in case more than one observation is available over that time period.
obs.date <- dummy$Group.1
obs.val <- dummy$x
rm(dummy)
ix <- which(input.dates %in% obs.date)

for (i in 1:dim(params)[1])
{
  # goodness of fit
  pred.val <- output2[i, ix]
  err <- sum(abs(pred.val - obs.val)) / length(obs.val) # mean absolute error
  params$error2[i] <- err # adds to the parameter object the prediction error associated with output 2
}

# step 2: search for the best fit within a given range of parameters ##

pf.range <- c(9, 30) # range of piston-flow transit times [y]. User defined
exp.range <- c(1, 5) # range of piston-flow transit times [y]. User defined
n.range <- c(0, 1) # range of infiltration coefficient [-]. User defined

# ix=which(params$pf>=pf.range[1] & params$pf<=pf.range[2] & params$exp>=exp.range[1] & params$exp<=exp.range[2] & gw.iso<mean(stable.output$h2)+1 & gw.iso>mean(stable.output$h2)-1) # subsets the parameter space
# ix=which(params$pf>=pf.range[1] & params$pf<=pf.range[2] & params$exp>=exp.range[1] & params$exp<=exp.range[2]) # subsets the parameter space
ix <- which(params$pf >= pf.range[1] & params$pf <= pf.range[2] & params$exp >= exp.range[1] & params$exp <= exp.range[2] & params$n >= n.range[1] & params$n <= n.range[2]) # subsets the parameter space
params.sub <- params[ix, ]
output1.sub <- output1[ix, ]
output2.sub <- output2[ix, ]

# ix=which(params.sub$error2==min(params.sub$error2)) # finds the best parameter combination to the stable isotope output function
# ix=which(params.sub$error1==min(params.sub$error1)) # finds the best parametern to the tritium output function
error.mult <- params.sub$error1 * params.sub$error2 # calculates the product of both measures of fit
# error.mult=params.sub$error2 # calculates the product of both measures of fit
ix <- which(error.mult == min(error.mult)) # finds the best parametern to both output functions

# step 3: save the results as graph ##

xmin.h3 <- as.Date("01.01.2000", format = "%d.%m.%Y") # date of the lower x-axis limit for tritium. User defined
xmax.h3 <- as.Date("01.01.2020", format = "%d.%m.%Y") # date of the upper x-axis limit for tritium. User defined

#TO CHANGE
setwd(graph.dir)
bitmap("fitting_0.1.png", type = "png16m", height = 10, width = 20, res = 200, encoding = "WinAnsi.enc") # the WinAnsi encoding is necessary to use the "per mille" sign (in Linux). Might not work on a windows computer. Please find a work around, or think about changing to linux, which is an open source, community supported, efficient, and secure operating system.
source(paste(mod.dir,"/visualisation_fitting.Rhistory",sep=""))
dev.off()

## step 4: save the results as table ##

#TO CHANGE
setwd(output.dir)

res <- params.sub[ix, ]
res[4:5] <- round(res[4:5], digits = 2)
write.table(res, "best_fit_1.csv", sep = ",", row.names = F)
