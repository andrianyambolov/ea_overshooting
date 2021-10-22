rm(list=ls())
cat("\014")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library('OECD')
library('countrycode')
library("readxl")


get_data <- function(ind, countries, start_year, end_year) {
  df <- data.frame(time = seq(start_year, end_year + 11/12, 1/12))
  for (cc in countries){
    envi <- new.env()
    assign("flag", FALSE, env=envi)
    print(cc)
    Sys.sleep(0.5)
    
    tryCatch({
      fetch_data <- get_dataset("MEI_FIN",
                                filter = paste(ind, cc, 'M', sep = '.'),
                                start_time = T0, end_time = Tend, pre_formatted = TRUE)
      dat <- cbind(fetch_data$obsTime, fetch_data$obsValue)
      dat[,2] <- as.numeric(dat[,2])
      dat[,1] <- as.numeric(substr(dat[,1], 1, 4)) + (as.numeric(substr(dat[,1], 6, 8)) - 1)/12
      dat <- data.frame(dat)
      if (cc != 'EA19'){
        current_code <-  countrycode(cc, origin = 'iso3c', destination = 'iso2c')
      }
      else{
        current_code <- 'EA'
      }
      colnames(dat) <- c('time',  current_code)
      df <- merge(df, dat, by = 'time', all.x = TRUE)
    }, error = function(msg) {
      assign("flag", TRUE, env=envi)
    })
    
    if (get("flag", env=envi) == TRUE){
      df[cc] <- NA
    }
    
  }
  return(df)
}

T0 <- 1999
Tend <- 2019
country_codes <- c('Australia', 'Canada', 'United Kingdom', 'Japan', 'South Korea', 'Sweden')
country_codes <-  countrycode(country_codes, origin = 'country.name', destination = 'iso3c')
var_codes <- c('CCUS')
var_names <- c('fx')

out <- get_data('CCUS', 'EA19', T0, Tend)

for (ii in 1:length(var_codes)){
  print(var_names[ii])
  aux <- get_data(var_codes[ii], country_codes, T0, Tend)
  out <- merge(aux, out,by="time")
}

fx <- out

for (ii in 2:(length(fx)-1)){
  fx[,ii] <- as.numeric(out[,ii])/as.numeric(out[,length(out)])
}
names(fx) <- c('time', paste('fx', tolower(colnames(fx[,-1])), sep="_"))
fx[,length(fx)] <-as.numeric(fx[,length(fx)])^(-1)
names(fx)[length(fx)] <- 'fx_us'
write.csv(fx, 'output/fx.csv', row.names = FALSE)