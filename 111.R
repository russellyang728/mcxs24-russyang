library(readabs)
library(readrba)
library(xts)
library(ggplot2)
library(patchwork)
library(tseries)
library(knitr)
library(zoo)

#download data


cr_dwld <- readrba::read_rba(series_id = "FIRMMCRTD")

cr <- cr_dwld[, c("date", "value")]

cr <- xts::xts(cr$value,cr$date)

cr <- xts::to.quarterly(cr, OHLC = FALSE)

start_date <- start(cr)
end_date1 <- as.Date("2023-12-31")
end_date2 <- as.yearqtr("2023 Q4")

cr = window(cr, start=start_date, end=end_date2)

cr_p = autoplot(cr) +
  
  theme_classic()+
  
  labs(title = "Cash Rate Target")+
  
  theme(axis.title.x = element_blank(),
        
        plot.title = element_text(hjust = 0.5,
                                  
                                  face = "bold"))



#unemployment rate data
unemp_dwnld <- readabs::read_abs(series_id = 'A84423050A') 


unemp <- xts(unemp_dwnld$value, order.by=as.Date(unemp_dwnld$date))

# Covert monthly unemployment data to quarterly for analysis compatibility. 

unemp <- to.quarterly(unemp,OHLC = FALSE)

unemp <- window(unemp, start=start_date, end=end_date2)

unemp_p = autoplot(unemp) +
  
  theme_classic()+
  
  scale_x_yearqtr(format = "%Y")+
  
  labs(title = "Unemployment Rate")+
  
  theme(axis.title.x = element_blank(),
        
        plot.title = element_text(hjust = 0.5,
                                  
                                  face = "bold"))


#government spending data

govspend_dwnld <- readabs::read_abs(series_id = 'A2304080V')

govspend <- xts(govspend_dwnld$value, order.by=as.Date(govspend_dwnld$date))

govspend <- window(govspend, start=start_date, end=end_date1)

govspend_p = autoplot(govspend) +
  
  theme_classic()+
  
  labs(title = "Government Spending")+
  
  theme(axis.title.x = element_blank(),
        
        plot.title = element_text(hjust = 0.5,
                                  
                                  face = "bold"))


#tax revenue data

taxrev_dwnld <- readabs::read_abs(series_id = 'A2302794K')

taxrev <- xts(taxrev_dwnld$value, order.by=as.Date(taxrev_dwnld$date))

taxrev <- window(taxrev, start=start_date, end=end_date1)

taxrev_p = autoplot(taxrev) +
  
  theme_classic()+
  
  labs(title = "Tax Revenue")+
  
  theme(axis.title.x = element_blank(),
        
        plot.title = element_text(hjust = 0.5,
                                  
                                  face = "bold"))


#inflation data

infl_dwnld <- readabs::read_abs(series_id = 'A2325850V')

infl <- xts(infl_dwnld$value, order.by=as.Date(infl_dwnld$date))

infl <- window(infl, start=start_date, end=end_date1)

infl_p = autoplot(infl) +
  
  theme_classic()+
  
  labs(title = "Inflation")+
  
  theme(axis.title.x = element_blank(),
        
        plot.title = element_text(hjust = 0.5,
                                  
                                  face = "bold"))


# real gdp data
rgdp_dwld <- readrba::read_rba(series_id = "GGDPCVGDP")

rgdp <- rgdp_dwld[, c("date", "value")]

rgdp$quarter <- zoo::as.yearqtr(rgdp$date)

rgdp <- xts::xts(rgdp$value,rgdp$quarter)

rgdp <- window(rgdp, start=start_date, end=end_date2)

rgdp_p = autoplot(rgdp) +
  
  theme_classic()+
  
  labs(title = "Real GDP")+
  
  theme(axis.title.x = element_blank(),
        
        plot.title = element_text(hjust = 0.5,
                                  
                                  face = "bold"))

length(unemp)
length(govspend)
length(taxrev)
length(infl)
length(rgdp)
length(cr)

data_df <- data.frame(
  unemp = coredata(unemp),
  govspend = coredata(govspend),
  taxrev = coredata(taxrev),
  infl = coredata(infl),
  rgdp = coredata(rgdp),
  cr = coredata(cr)
)
data_df

variable = c('unemp', 'govspend', 'taxrev','infl', 'rgdp','cr')
N <- rep(136, length(variable))
Mean <- sapply(data_df, mean)
SD <- sapply(data_df, sd)
Min <- sapply(data_df, min)
Max <- sapply(data_df, max)
table1 =data.frame(N, Mean,SD, Min, Max)


knitr::kable(table1, caption = "Summary Statistics", digits = 2)


Sign = c('/','+','-','/','+','/')

table2 = data.frame(variable,Sign)
knitr::kable(table2, caption = "Sign Restriction", digits = 2)




set.seed(2024)
RW1 <- arima.sim(model= list(order = c(0, 1, 0)), n=1000, mean=0, sd=1)
plot.ts(RW1,main="Random Walk 1", col=4,xlab="")


RW2 <- arima.sim(model= list(order = c(0, 1, 0)), n=1000, mean=0, sd=1)
plot.ts(RW2,main="Random Walk 2", col=4,xlab="")

RW  <- cbind(RW1,RW2)
RW
