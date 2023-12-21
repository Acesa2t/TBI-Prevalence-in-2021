

#' This script contains all the functions used by the script titled
#' 'TB infection and reactivation rate creation.R'.


#' Fills down the values in a character variable
filltheblanks <- function(x, missing = ""){
  
  rle <- rle(as.character(x))
  empty <- which(rle$value == missing)
  rle$values[empty] <- rle$value[empty - 1] 
  inverse.rle(rle)
  
}


#' Prepares the wide census data and put it in long format
CleanseCensus <- function(dt) {

  #' Make sure the dataset is a datatable
  dt <- as.data.table(dt)
  
  #' Replace the first two column names with age and year of arrival
  setnames(dt, 1, "AGEP") 
  setnames(dt, 2 ,"YARP")
  
  #' Clean the year of arrival column
  #' Fill down the blank rows in the age column
  dt[ , AGEP := filltheblanks(AGEP)]
  #' Remove " years" if it is there
  dt[ , AGEP := gsub(" years", "", AGEP)]
  dt[ , AGEP := gsub(" and over", "", AGEP)]
  #' Get rid of rows that aren't integer
  dt <- dt[!is.na(as.integer(as.character(dt$AGEP))),]
  dt[, AGEP := as.numeric(AGEP)]
  
  #' Clean the year of arrival column
  dt[YARP == "Arrived 1 Jan 2011 - 9 August 2011", YARP := "2011"]
  dt[ , YARP := gsub("Arrived ", "", YARP)]
  dt[YARP == "Not stated", YARP := "1111"]
  dt[YARP == "Not applicable", YARP := "1111"]
  dt <- dt[!is.na(as.integer(as.character(dt$YARP))),]
  dt[, YARP := as.numeric(YARP)]
  dt[YARP == 1111, YARP := NA]
  
  #' Find the highest year of arrival (i.e.2006/2011/2016) 
  censusyear <- max(dt$YARP, na.rm = T)
  
  #' Reshape, of country of birth becomes one column
  dt <- melt(dt, id = c("AGEP", "YARP"))
  
  #' Renaming country of birth (cob) variable
  setnames(dt, "variable", "cob") 
  setnames(dt, "value", "NUMP") 
  dt <- subset(dt, cob != "Total")
  
  #' Create year of birth (yob) column
  dt <- as.data.table(dt)
  dt[, CNSY := as.numeric(censusyear)]
  
  #' Create year of birth (yob) column
  dt[, YOBP := CNSY - AGEP]
  
  #' Remove rows with no population
  dt <- subset(dt, NUMP != 0)
  
  #' Creating a column of ISO3 codes 
  dt[, ISO3 := countrycode(cob, "country.name", "iso3c")]
  
  #' Fix the countries that didn't convert properly 
  dt <- ISO3fixfunc(dt)
  
  #' Remove country of birth column and aggregate results
  dt[, cob := NULL]
  dt <- dt[, list(NUMP = sum(NUMP)), 
           by = c("AGEP", "YARP", "YOBP", "ISO3", "CNSY")]
  
  dt
  
}


#' Prepares the TB data and puts 
CleanseTBdata <- function(dt) {
  
  #' Make sure the dataset is a datatable
  dt <- as.data.table(dt)
  
  #' Replace column names 
  setnames(dt, "Diagnosis.Year.Month", "notiyear")
  setnames(dt, "Age", "AGEP")
  dt[, AGEP := as.numeric(AGEP)]
  setnames(dt, "Sex", "SEXP")
  setnames(dt, "Country.of.birth", "cob")
  setnames(dt, "Year.of.first.arrival", "YARP")
  setnames(dt, "Clinical.Presentation", "mani")
  setnames(dt, "State", "state")
  setnames(dt, "RA", "gcc")
  
  dt[dt$YARP == "Not applicable","YARP" := NA]
  dt[dt$YARP == "Unknown","YARP" := NA]
  dt[, YARP := as.character(YARP)]
  dt[, YARP := as.numeric(YARP)]
  
  #' Creating a column of ISO3 codes 
  dt$ISO3 <- countrycode(dt$cob, "country.name", "iso3c")
  
  #' Fixing ISO3 codes for those that didn't convert, and changing those
  #' that don't appear in Houben and Dodd's dataset
  dt <- ISO3fixfunc(dt)
  
  ####### Remove PNG cases from Queensland #####################
  dt[, p := 0]
  dt[ISO3 == "PNG" & state == "QLD", p := 1]
  dt <- dt[dt$p != 1,]
  dt[, p := NULL]
  ##########################################################
  
  #' Sorting the month and year columns (i.e. originally called notiyear)
  revsubstr <- function(x, start, stop) {
    x <- strsplit(x, "")
    sapply(x, 
           function(x) paste(rev(rev(x)[start:stop]), collapse = ""), 
           USE.NAMES = FALSE)
  }
  dt[, notiyear := as.character(notiyear)]
  dt[, year := substr(notiyear, start = 1, stop = 4)]
  dt[, poptb := 1]
  
  #' Get rid of columns we don't need for now: cob, mani, gcc
  dt <- dt[, c("year", "AGEP", "YARP", "ISO3", "poptb")] %>%
    group_by(ISO3, YARP, AGEP, year) %>% # aggregating  
    summarise_all(list(sum))
  sum(dt$poptb, na.rm = TRUE)
  
  dt
  
}


#' Creates the master look- up table.
CreateProbTables <-  function() {
  
  # Use census and tb tables to get a list of unique ISO codes.
  x <- unique(c(census$ISO3, tb$ISO3))
  x <- x[!is.na(x)]
  
  # Use census table to get a list of unique ISO codes.
  census.year <- unique(census$CNSY)
  
  # Create two sequences for year of birth and year of arrival.
  y <- 1900:census.year
  z <- 1900:census.year
  
  # Pass all three variable to expand() to create a data table.
  prob.Inf <- as.data.table(expand.grid(ISO3 = x, YARP = y, YOBP = z, 
                                        CNSY = census.year, stringsAsFactors = FALSE))
  
  # Removing rows with 'birth after arriving' and censusralia.
  prob.Inf <- prob.Inf[!(YOBP > YARP), ]
  prob.Inf <- prob.Inf[ISO3 != census.iso3, ]
  
  # Creating rows for censusralia with year of arrival as 'NA'.
  prob.InfAus <- as.data.table(expand.grid(ISO3 = census.iso3, 
                                           YARP = 1L, YOBP = z, 
                                           CNSY = census.year, stringsAsFactors = FALSE)) 
  prob.InfAus$YARP <- NA
  
  # Combining and reordering the data.table.
  prob.Inf <-  rbind(prob.Inf, prob.InfAus)
  prob.Inf <- setorder(prob.Inf, ISO3, YOBP, YARP) 
  
  # Creating additional key columns used for subsetting data.tables .
  # YOBPp1 - year of birth plus one year.
  # YARPm1 - year of arrival minus one year.
  # YARPp1 - year of arrival plus one year.
  
  prob.Inf[, c("YOBPp1", "YARPm1", "YARPp1") := .(as.integer(YOBP + 1), 
                                                  as.integer(YARP - 1), 
                                                  as.integer(YARP + 1))]
  
  return (prob.Inf)
  
}

#' Prepares Houben & Dodd's datasets for calculations
tbhazprep.function <- function(dt){
  
  #' Converts the log annual risk of infection estimates (lari)
  #' to hazards (FOIs).
  dt[, FOI := exp(lari)]
  
  dt[, partFOI := 1]
  dt[, cumhaz := 1]
  dt[, halfFOI := 1]
  dt[, quarterFOI := 1]

  #' #' Create data for 2015 & 2016 (by just duplicating
  #' #' the 2014 data)
  # dt2015 <- dt[year == 2014]
  # dt2016 <- dt2015
  # 
  # dt2015$year <- 2015
  # dt2016$year <- 2016
  # 
  # dt <- rbind(dt2016, dt2015, dt)

  #' Now create pre 1934 data (by just duplicating
  #' the 1934 data)
  dt1934 <- subset(dt, dt$year == 1934)
  
  df <- NULL # initialise a container
  for(i in 1:45) {
    dfnew <- dt1934
    dfnew$year <- 1934 - i
    if (is.null(df)) {
      df <- dfnew
    }
    else {
      df <- rbind(dfnew, df)
    }
  }
  
  #' Merge the pre 1934 data table to the other one
  dt <- rbind(df, dt) 
  
  #' Make sure the data table is ordered
  dt <- setorder(dt, iso3, -replicate, year)
  rownames(dt) <- seq(length = nrow(dt))
  
  #' If the year is a census year then replace the
  #' FOI with a partial FOI, becuase the census was 
  #' performed in August (still not really sure 
  #' whether we should have done this)
  dt[, partFOI := FOI * (220/365.2425)]
  dt[year == census.year, FOI := partFOI]
  dt[, partFOI := NULL]
  
  #' Create a column with half the FOI (for years of birth and migration)
  dt[, halfFOI := FOI * 0.5]
  
  #' Create a column with a quarter of the FOI
  dt[, quarterFOI := FOI * 0.25]
  
  #' Create the cumulative sum column
  dt[, cumhaz := cumsum(FOI), by = list(iso3, replicate)]
  setkey(dt, iso3)
  
  return(dt)
  
}


ISO3fixfunc <- function(DT){
  
  DT <- as.data.table(DT)
  DT[, cob := gsub("..", ".", cob, fixed = TRUE)]
  DT[, cob := gsub(",", ".", cob, fixed = TRUE)]
  DT[, cob := gsub(".", " ", cob, fixed = TRUE)]
  DT[, cob := gsub("-", "", cob, fixed = TRUE)]
  DT[, cob := gsub("nec", "", cob, fixed = TRUE)]
  DT[, cob := gsub("nfd", "", cob, fixed = TRUE)]
  DT[, cob := gsub("nf", "", cob, fixed = TRUE)]
  DT[, cob := gsub("  ", " ", cob, fixed = TRUE)]
  
  trim <- function (x) gsub("^\\s+|\\s+$", "", x)
  DT[, cob := trim(cob)]
  DT[cob == "At sea", ISO3 := NA]
  DT[cob == "Not stated", ISO3 := NA]
  DT[cob == "Not applicable", ISO3 := NA]
  DT[cob == "Not elsewhere classified", ISO3 := NA]
  DT[cob == "Overseas visitor", ISO3 := NA]
  DT[cob == "Inadequately described", ISO3 := NA]
  DT[cob == "No data available", ISO3 := NA]
  DT[cob == "Overseas  Country unknown", ISO3 := NA]
  DT[cob == "Unknown", ISO3 := NA]
  
  #' Australian regions
  DT[cob == "Australia includes External Territories", ISO3 := "AUS"] 
  DT[cob == "Norfolk Island", ISO3 := "AUS"]
  DT[cob == "Australian External Territories", ISO3 := "AUS"]
  DT[cob == "Oceania and Antarctica", ISO3 := "NZL"]
  
  #' Asian countries
  DT[cob == "China (excludes SARs and Taiwan)", ISO3 := "CHN"]
  DT[cob == "China excludes SARs and Taiwan Province", ISO3 := "CHN"]
  DT[cob == "Japan and the Koreas", ISO3 := "JPN"]
  DT[cob == "China excludes SARs and Taiwan", ISO3 := "CHN"]
  DT[ISO3 == "TWN", ISO3 := "CHN"] #Taiwan
  
  #' Regions
  DT[cob == "Middle East", ISO3 := "SAU"]
  DT[cob == "North Africa", ISO3 := "MAR"]
  DT[cob == "North Africa and the Middle East", ISO3 := "EGY"]
  DT[cob == "Southern and East Africa", ISO3 := "SOM"]
  DT[cob == "Africa", ISO3 := "SOM"]
  DT[cob == "Spanish North Africa", ISO3 := "ESP"]
  DT[cob == "Sub Saharan Africa", ISO3 := "SOM"]
  DT[cob == "Sub-Saharan Africa, ", ISO3 := "SOM"]
  DT[cob == "SubSaharan Africa", ISO3 := "SOM"] 
  DT[cob == "Central and West Africa", ISO3 := "CMR"]
  
  DT[cob == "Eastern Europe", ISO3 := "POL"]
  DT[cob == "North West Europe", ISO3 := "FRA"]
  DT[cob == "NorthWest Europe", ISO3 := "FRA"]
  DT[cob == "Southern Europe", ISO3 := "GRC"]
  DT[cob == "South Eastern Europe", ISO3 := "MKD"]
  DT[cob == "Southern and Eastern Europe", ISO3 := "MKD"]
  DT[cob == "Northern Europe", ISO3 := "DNK"]
  DT[cob == "Western Europe", ISO3 := "FRA"]
  
  DT[cob == "Americas", ISO3 := "USA"]
  DT[cob == "Polynesia excludes Hawaii", ISO3 := "WSM"]
  DT[cob == "Polynesia (excludes Hawaii)", ISO3 := "WSM"]
  
  DT[cob == "South America", ISO3 := "BRA"]
  DT[cob == "Northern America", ISO3 := "USA"]
  DT[cob == "Central America", ISO3 := "CRI"]
  DT[cob == "Virgin Islands United States", ISO3 := "USA"]
  DT[cob == "Caribbean", ISO3 := "BHS"]
  DT[cob == "Virgin Islands United States", ISO3 := "USA"]
  
  DT[cob == "Central Asia", ISO3 := "KAZ"]
  DT[cob == "South East Asia", ISO3 := "THA"]
  DT[cob == "SouthEast Asia", ISO3 := "THA"]
  DT[cob == "Southern Asia", ISO3 := "IND"]
  DT[cob == "North East Asia", ISO3 := "MNG"]
  DT[cob == "NorthEast Asia", ISO3 := "MNG"]
  
  DT[cob == "Mainland South East Asia", ISO3 := "THA"]
  DT[cob == "Mainland SouthEast Asia", ISO3 := "THA"]
  DT[cob == "Maritime South East Asia", ISO3 := "PHL"]
  DT[cob == "Maritime SouthEast Asia", ISO3 := "PHL"]
  DT[cob == "Melanesia", ISO3 := "VUT"]
  DT[cob == "Micronesia", ISO3 := "KIR"]
  DT[cob == "Southern and Central Asia", ISO3 := "AFG"]
  
  #' Great Britain and surrounding
  DT[cob == "Northern Ireland", ISO3 := "GBR"]
  DT[cob == "Scotland", ISO3 := "GBR"]
  DT[cob == "Wales", ISO3 := "GBR"]
  DT[cob == "Pitcairn Islands", ISO3 := "GBR"]  
  DT[cob == "England", ISO3 := "GBR"]
  DT[cob == "United Kingdom Channels Islands and Isle of Man", ISO3 := "GBR"]  
  DT[cob == "United Kingdom Channel Islands and Isle of Man", ISO3 := "GBR"]  
  DT[cob == "United Kingdom  Channel Islands and Isle of Man", ISO3 := "GBR"] 
  DT[cob == "United Kingdom, Channel Islands and Isle of Man", ISO3 := "GBR"] 
  DT[cob == "Channel Islands and Isle of Man", ISO3 := "GBR"]
  DT[cob == "Channel Islands", ISO3 := "GBR"]
  DT[cob == "United Kingdom", ISO3 := "GBR"]
  DT[cob == "British Antarctic Territory", ISO3 := "GBR"]
  
  #' Other
  DT[cob == "Kosovo", ISO3 := "SRB"]
  DT[ISO3 == "ESH", ISO3 := "MRT"]
  DT[cob == "Former Yugoslavia", ISO3 := "MKD"]
  DT[cob == "Yugoslavia Federal Republic of", ISO3 := "MKD"]
  DT[cob == "SFR Yugoslavia", ISO3 := "MKD"]
  DT[ISO3 == "YUG", ISO3 := "MRT"]
  DT[cob == "Yugoslavia, Federal Republic of", ISO3 := "MKD"]
  DT[cob == "Sao Tom? and Principe", ISO3 := "GAB"] 
  
  
  #' European territories
  DT[cob == "RÃ union", ISO3 := "FRA"]
  DT[cob == "Réunion", ISO3 := "FRA"] 
  DT[cob == "R? union", ISO3 := "FRA"]
  DT[cob == "St Barthelemy", ISO3 := "FRA"]
  DT[cob == "St Martin French part", ISO3 := "FRA"]
  DT[cob == "Bonaire Sint Eustatius and Saba", ISO3 := "NLD"]
  DT[cob == "Aland Islands", ISO3 := "SWE"]
  DT[cob == "Netherlands Antilles", ISO3 := "NLD"]
  
  #' Removing YARP for Norfolk Islanders since I've made their ISO3 Australia anyway
  DT[cob == "Norfolk Island", YARP := NA]
  DT[ISO3 == "IMN", ISO3 := "GBR" ] #Isle of Man
  DT[ISO3 == "GGY", ISO3 := "GBR" ] #Guernsey
  DT[ISO3 == "JEY", ISO3 := "GBR" ] #Jersey
  DT[ISO3 == "GIB", ISO3 := "ESP" ] #Gibraltar
  DT[ISO3 == "TWN", ISO3 := "CHN" ] #Taiwan
  DT[ISO3 == "SHN", ISO3 := "GBR" ] #Saint Helena
  DT[ISO3 == "FLK", ISO3 := "GBR" ] #Falkland Islands (Malvinas)
  DT[ISO3 == "REU", ISO3 := "FRA" ] #Reunion
  DT[ISO3 == "NFK", ISO3 := "AUS" ] #Norfolk Island
  DT[ISO3 == "ATA", ISO3 := "GBR" ] #Antarctica
  DT[ISO3 == "LIE", ISO3 := "AUT" ] #Liechtenstein
  DT[ISO3 == "FRO", ISO3 := "NOR" ] #Faroe Islands
  DT[ISO3 == "VAT", ISO3 := "ITA" ] #Holy See (Vatican City State)
  DT[ISO3 == "ESH", ISO3 := "MAR" ] #Western Sahara
  DT[ISO3 == "SPM", ISO3 := "FRA" ] #Saint Pierre and Miquelon
  DT[ISO3 == "TWN", ISO3 := "CHN" ] #Taiwan, Province of China
  DT[ISO3 == "GLP", ISO3 := "FRA" ] #Guadeloupe
  DT[ISO3 == "GUF", ISO3 := "FRA" ] #French Guiana
  DT[ISO3 == "MTQ", ISO3 := "FRA" ] #Martinique
  DT[ISO3 == "MYT", ISO3 := "FRA" ] #Mayotte
  DT[ISO3 == "SHN", ISO3 := "GBR" ] #Saint Helena, Ascension and Tristan da Cunha
  DT
  
}

#' Calculating the hazard for the local born
TBhazard.calc.function.local.born <- function(DT, haz.ref) {
  
  #' Calculating the hazard for the local born
  
  #' This places a list of all FOIs (all replicates) for the year of birth, 
  #' into the Ha1 field. The FOIs are also halved, because we assume that 
  #' the population birth dates were evenly distributed throughout 
  #' their birth year. See further below for an adjustment to account 
  #' for those in the population who are born in the census year.
  DT[haz.ref[iso3 == census.iso3, .(list(halfFOI)), by = year],
     Ha1 := V1,
     on = c(YOBP = "year")]
  
  #' This takes the full list of cumulative hazards from the earliest year 
  #' to one year after the year of birth (all replicates), and places them in 
  #' the Ha2a field.
  DT[haz.ref[iso3 == census.iso3, .(list(cumhaz)), by = year],
     Ha2a := V1,
     on = c(YOBPp1 = "year")]
  
  #' This takes the full list of cumulative hazards from the earliest year 
  #' to the census year (all replicates), and places them in 
  #' the Ha2b field...so the values are actually the same in each row.
  DT[haz.ref[iso3 == census.iso3, .(list(cumhaz)), by = year],
     Ha2b := V1,
     on = c(CNSY = "year")]
  
  #' This finds the difference between the cumulative hazards in the above
  #' two variables/fields (i.e. from the year of birth plus one to the 
  #' census year). So this then this gives the cumulative 
  #' hazard from the year after the year of birth until the census year.
  DT[, Ha2 := .(mapply('-',Ha2b, Ha2a, SIMPLIFY = F)),]
  
  #' Getting rid of variable/s we no longer need.
  DT[, c("Ha2a", "Ha2b") := .(NULL, NULL), ]
  
  #' Now we can add the hazards for the year of birth, census year and the 
  #' intervening period together.
  DT[, H := .(mapply('+', Ha1, Ha2, SIMPLIFY = F)),]
  
  #' Getting rid of variable/s we no longer need.
  DT[, c("Ha1", "Ha2") := .(NULL), ]
  
  #' If the year of birth is the same as the census year, then the FOIs
  #' calculated above are wrong (and negative), so the calculation below  
  #' adjusts the FOI for these individuals to simply be half of the
  #' FOI in the census year, which has already been apportioned. 
  DT[haz.ref[iso3 == census.iso3, .(list(halfFOI)), by = year],
     Ha.half.census := V1,
     on = c(CNSY = "year")]
  
  DT[YOBP == CNSY, H := Ha.half.census]
  
  #' If the year of birth is the year before the census year, then the FOIs
  #' calculated above are too small (the FOI in the census year is not
  #' included), so the calculation below sorts this out
  DT[haz.ref[iso3 == census.iso3, .(list(FOI)), by = year],
     Ha.census := V1,
     on = c(CNSY = "year")]
  
  DT[YOBP == CNSY - 1, H := .(mapply('+', H, Ha.census, SIMPLIFY = F)),]
  
  #' This finds the median of all the hazards
  #DT[,  H.med := .(mapply(median, H, SIMPLIFY = T)),]
  #DT[is.na(H.med), H.med := 0, ]
  
  #' quantile.func <- function (x) {
  #'   quantile(x, probs = c(low.percentile))
  #' }
  #' 
  #' #' This finds the desired low percentile of the same and puts it in H.low
  #' DT[,  H.low := .(mapply(quantile.func, H, SIMPLIFY = T)),]
  #' DT[is.na(H.low), H.low := 0, ]
  #' 
  #' quantile.func <- function (x) {
  #'   quantile(x, probs = c(high.percentile))
  #' }
  
  #' This finds the desired high percentile of the same and puts it in H.high
  # DT[,  H.high := .(mapply(quantile.func, H, SIMPLIFY = T)),]
  # DT[is.na(H.high), H.high := 0, ]
  
  #DT[, H := .(NULL), ]
  
  #' Getting rid of variable/s we no longer need.
  DT[, c("YOBPp1", "YARPm1", "Ha.half.census",
         "YARPp1", "Ha.census") := 
       .(NULL, NULL, NULL, NULL, NULL), ]
}


#' Calculating the hazard for the overseas born
TBhazard.calc.function.overseas.born <- function(DT, haz.ref){
  #' 5000 set - Overseas components
  #' This places a list of all FOIs (all replicates) for the year of birth, 
  #' into the Ho1 field. The FOIs are also halved, because we assume that 
  #' the population birth dates were evenly distributed throughout 
  #' their birth year. See further below for an adjustment to account 
  #' for those born in the year of migration or the census year.
  DT[haz.ref[, .(list(halfFOI)), by = .(year,iso3)], 
     Ho1 := V1,
     on = c(YOBP = "year", ISO3 = "iso3")]
  
  #' This takes the full list of cumulative hazards from the earliest year 
  #' to one year after the year of birth (all replicates), and places them in 
  #' the Ho2a field.
  DT[haz.ref[, .(list(cumhaz)), by = .(year, iso3)],
     Ho2a := V1,
     on = c(YOBPp1 = "year", ISO3 = "iso3")]

  #' This takes the full list of cumulative hazards from the earliest year 
  #' to one year before the year of migration (all replicates), and places 
  #' them in the Ho2b field.
  DT[haz.ref[, .(list(cumhaz)), by = .(year, iso3)],
     Ho2b := V1,
     on = c(YARPm1 = "year", ISO3 = "iso3")]

  #' This finds the difference between the cumulative hazards in the above
  #' two variables/fields (i.e. from the year of birth plus one to the 
  #' year before migration). So this then this gives the cumulative 
  #' hazard from the year after the year of birth until the year prior
  #' to migration.
  DT[, Ho2 := .(mapply('-', Ho2b, Ho2a, SIMPLIFY = F)),]
  
  #' Getting rid of variable/s we no longer need.
  DT[, c("Ho2a", "Ho2b") := .(NULL, NULL), ]
  
  #' If the year of migration is the year after the year of 
  #' birth the above calculation creates a negative number.
  #' Therefore, for these people Ho2 needs to just be
  #' adjusted to zero
  DT[YOBP == YARP - 1, Ho2 := 0, ]
  
  gc() # Frees up memory, supposedly, and reports the memory usage.
  
  #' This places a list of all FOIs (all replicates) for the year of migration
  #' in the country of birth into the Ho3 field. The FOIs are also halved, 
  #' because we assume that the population migration dates were evenly 
  #' distributed throughout their migration year, and so they will also 
  #' experience the risk in the country to which they are migrating.
  DT[haz.ref[, .(list(halfFOI)), by = .(year, iso3)],
     Ho3 := V1,
     on = c(YARP = "year", ISO3 = "iso3")]
  
  #' If the year of birth is the same as the year of arrival, then
  #' only a quarter of the FOI in that year for the country of birth
  #' should be applied.
  DT[haz.ref[, .(list(quarterFOI)), by = .(year, iso3)],
     Ho3.quart := V1,
     on = c(YARP = "year", ISO3 = "iso3")]

  DT[YOBP == YARP, Ho3 := Ho3.quart, ]
  
  #' If the year of birth and migration is also the
  #' census year, in this case we just apply a quarter of the hazard
  #' in both the country of origin and destination as well.
  #' Not exactly right, but close enough.
  #' This will be achieved with the calculation
  #' above and below anyway.
  
  #' I also need to get rid of the FOI calculated above in the 
  #' year of birth, and subsequently, otherwise it will be counted twice
  DT[YOBP == YARP, Ho1 := 0, ]
  DT[YOBP == YARP, Ho2 := 0, ]
  
  #' Getting rid of variable/s we no longer need.
  DT[, c("Ho3.quart") := .(NULL), ]
  
  #' If the year of birth and migration is also the
  #' census year, in this case we just apply a quarter of the hazard
  #' in both the country of origin and destination as well.
  #' Not exactly right, but close enough.
  #' This will be achieved with the calculation
  #' above and below anyway.
  
  #' This adds together all of the overseas hazards together
  DT[, Ho := .(mapply('+', Ho1, Ho2, SIMPLIFY = F)),]
  DT[, Ho := .(mapply('+', Ho, Ho3, SIMPLIFY = F)),]
  
  #' Getting rid of variable/s we no longer need.
  DT[, c("Ho1", "Ho2", "Ho3") := .(NULL, NULL, NULL), ]
  
  
  gc() # Frees up memory, supposedly, and reports the memory usage.  
  
  #' 200 set - Australian components
  
  #' This places a list of all FOIs (all replicates) for the year of migration
  #' in the destination country into the Ha1 field. The FOIs are also halved, 
  #' because we assume that the population migration dates were evenly 
  #' distributed throughout their migration year, and so they will also 
  #' experience the risk in the country from which they are migrating.
  DT[haz.ref[iso3 == census.iso3, .(list(halfFOI)), by = year],
     Ha1 := V1,
     on = c(YARP = "year")]
  
  #' This takes the full list of cumulative hazards from the earliest year 
  #' to one year after the year of migration (all replicates), and places them in 
  #' the Ha2a field.
  DT[haz.ref[iso3 == census.iso3, .(list(cumhaz)), by = year],
     Ha2a := V1,
     on = c(YARPp1 = "year")]
  
  #' This takes the full list of cumulative hazards from the earliest year 
  #' to the census year (all replicates), and places them in 
  #' the Ha2b field...so the values are actually the same in each row.
  DT[haz.ref[iso3 == census.iso3, .(list(cumhaz)), by = year],
     Ha2b := V1,
     on = c(CNSY = "year")]
  
  gc() # Frees up memory, supposedly, and reports the memory usage.  
  
  #' This finds the difference between the cumulative hazards in the above
  #' two variables/fields (i.e. from the year of migration plus one to the 
  #' census year). So this then this gives the cumulative 
  #' hazard from the year after the year of migration until the census year.
  DT[, Ha2 := .(mapply('-',Ha2b, Ha2a, SIMPLIFY = F)), ]
  
  #' Getting rid of variable/s we no longer need.
  DT[, c("Ha2a", "Ha2b") := .(NULL, NULL), ]
  
  #' If the year of migration is the census year 
  #' the above calculation will be zero, but that's
  #' okay, because the hazard in that year will be
  #' captured in the Ha1 field.
  
  #' If the year of arrival is the year before the census
  #' the above calculation creates a negative number.
  #' Therefore, for these people, Ha2 needs to just be
  #' adjusted to simply equal the FOI in the census year
  DT[haz.ref[iso3 == census.iso3, .(list(halfFOI)), by = year],
     Ha.half.census := V1,
     on = c(CNSY = "year")]
  
  DT[YOBP == YARP - 1, Ha2 := Ha.half.census, ]
  
  gc() # Frees up memory, supposedly, and reports the memory usage.  
  
  #' Now we can add the hazards for the year of migration, census year and the 
  #' intervening period together.
  DT[, Ha := .(mapply('+', Ha1, Ha2, SIMPLIFY = F)),]
  
  #' Getting rid of variable/s we no longer need.
  DT[, c("Ha1", "Ha2", "Ha.half.census") := .(NULL), ]
  
  # str(DT)
  
  gc() # Frees up memory and report the memory usage.
  
  #' If the year of migration is also the
  #' census year, we'll just apply a quarter of the hazard
  #' in both the country of origin and destination as well.
  #' The overseas component has already been sorted. Below
  #' sorts the destination component
  DT[haz.ref[iso3 == census.iso3, .(list(quarterFOI)), by = year],
     Ha.quarter.census := V1,
     on = c(CNSY = "year")]
  
  DT[YARP == CNSY, Ha := Ha.quarter.census, ]
  
  #' Getting rid of variable/s we no longer need.
  DT[, c("Ha.quarter.census") := .(NULL), ]
  
  #' Now we can add the overseas and Australia hazards together.
  DT[, H := .(mapply('+', Ho, Ha, SIMPLIFY = F)),]
  
  #' Getting rid of variable/s we no longer need.
  DT[, c("Ho", "Ha") := .(NULL, NULL), ]
  
  gc() # Frees up memory, supposedly, and reports the memory usage.  
  
  #' This finds the median of all the hazard replicates
  #DT[, H := .(mapply(median, H, SIMPLIFY = T)),]

  
  #' This finds the median of all the hazards
  #' DT[,  H.med := .(mapply(median, H, SIMPLIFY = T)),]
  #' DT[is.na(H.med), H.med := 0, ]
  #' 
  #' quantile.func <- function (x) {
  #'   quantile(x, probs = c(low.percentile))
  #' }
  #' 
  #' #' This finds the desired low percentile of the same and puts it in H.low
  #' DT[,  H.low := .(mapply(quantile.func, H, SIMPLIFY = T)),]
  #' DT[is.na(H.low), H.low := 0, ]
  #' 
  #' quantile.func <- function (x) {
  #'   quantile(x, probs = c(high.percentile))
  #' }
  #' 
  #' This finds the desired high percentile of the same and puts it in H.high
  # DT[,  H.high := .(mapply(quantile.func, H, SIMPLIFY = T)),]
  # DT[is.na(H.high), H.high := 0, ]
  # 
  # DT[, H := .(NULL), ]
  
  #' Getting rid of variable/s we no longer need.
  DT[, c("YOBPp1", "YARPm1", "YARPp1") := 
       .(NULL, NULL, NULL), ]
  #View(DT)
}
