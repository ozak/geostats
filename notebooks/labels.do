* File assignes labels to most of the variables in the dataset
set maxvar 32000
cd "../data/"

use "Country-Stats.dta", clear

foreach stat in gini mean std max min sum count median majority{
    di "`stat'"
    foreach var of varlist *`stat' {
        local lbl : variable label `var'
        if "`lbl'" == "" {
            local mylabel = ""
            local name = substr("`var'", 1, length("`var'") - length("`stat'"))
            di "`name'"

            if strpos("`name'", "post1500AverageCal") > 0 {
                local mylabel = "Average Crop Yield (post-1500CE)"
            } 
            else if strpos("`name'", "post1500TotalCal") > 0 {
                local mylabel = "Total Crop Yield (post-1500CE)"
            }
            else if strpos("`name'", "post1500MaximumCal") > 0 {
                local mylabel = "Maximum Crop Yield (post-1500CE)"
            }
            else if strpos("`name'", "pre1500AverageCal") > 0 {
                local mylabel = "Average Crop Yield (pre-1500CE)"
            }
            else if strpos("`name'", "pre1500TotalCalories") > 0 {
                local mylabel = "Total Crop Yield (pre-1500CE)"
            }
            else if strpos("`name'", "pre1500MaximumCal") > 0 {
                local mylabel = "Maximum Crop Yield (pre-1500CE)"
            }
            else if strpos("`name'", "difAverageCal") > 0 {
                local mylabel = "Change Average Crop Yield"
            }
            else if strpos("`name'", "difTotalCal") > 0 {
                local mylabel = "Change Total Crop Yield"
            }
            else if strpos("`name'", "dif2MaximumCal") > 0 {
                local mylabel = "Change Maximum Crop Yield"
            }
            else if strpos("`name'", "dif2AverageCal") > 0 {
                local mylabel = "Change Average Crop Yield"
            }
            else if strpos("`name'", "dif2TotalCal") > 0 {
                local mylabel = "Change Total Crop Yield"
            }
            else if strpos("`name'", "dif2MaximumCal") > 0 {
                local mylabel = "Change Maximum Crop Yield"
            }
            else if substr("`name'", 1, 11) == "plowpospost" {
                local mylabel = "Plow Positive Crop CSI (post-1500CE)"
            }
            else if substr("`name'", 1, 10) == "plowpospre" {
                local mylabel = "Plow Positive Crop CSI (pre-1500CE)"
            }
            else if substr("`name'", 1, 11) == "plownegpost" {
                local mylabel = "Plow Negative Crop CSI (post-1500CE)"
            }
            else if substr("`name'", 1, 10) == "plownegpre" {
                local mylabel = "Plow Negative Crop CSI (pre-1500CE)"
            }
            else if substr("`name'", 1, 11) == "plowpotpost" {
                local mylabel = "Plow Potential CSI (post-1500CE)"
            }
            else if substr("`name'", 1, 10) == "plowpotpre" {
                local mylabel = "Plow Potential CSI (pre-1500CE)"
            }
            else if substr("`name'", 1, 10) == "plowposdif" {
                local mylabel = "Change Plow Positive Crop CSI"
            }
            else if substr("`name'", 1, 10) == "plownegdif" {
                local mylabel = "Change Plow Negative Crop CSI"
            }
            else if substr("`name'", 1, 10) == "plowpotdif" {
                local mylabel = "Change Plow Potential CSI"
            }
            else if strpos("`name'", "cldmean") > 0 {
                local mylabel = "Average Cloud Cover (%)"
            }
            else if strpos("`name'", "premean") > 0 {
                local mylabel = "Average Precipitation"
            }
            else if strpos("`name'", "tmpmean") > 0 {
                local mylabel = "Average Mean Temperature"
            }
            else if strpos("`name'", "tmnmean") > 0 {
                local mylabel = "Average Minimum Temperature"
            }
            else if strpos("`name'", "tmxmean") > 0 {
                local mylabel = "Average Maximum Temperature"
            }
            else if strpos("`name'", "dtrmean") > 0 {
                local mylabel = "Average Diurnal Temperature Range"
            }
            else if strpos("`name'", "rehmean") > 0 {
                local mylabel = "Average Relative Humidity (%)"
            }
            else if strpos("`name'", "petmean") > 0 {
                local mylabel = "Average Potential Evapotranspiration (mm/month)"
            }
            else if strpos("`name'", "frsmean") > 0 {
                local mylabel = "Average Frequency of Frost Days (days/month)"
            }
            else if strpos("`name'", "vapmean") > 0 {
                local mylabel = "Average Vapor Pressure (hPa)"
            }
            else if strpos("`name'", "wetmean") > 0 {
                local mylabel = "Average Wet Days (days/month)"
            }
            else if strpos("`name'", "cldspatcorr") > 0 {
                local mylabel = "Spatial Correlation Cloud Cover (%)"
            }
            else if strpos("`name'", "prespatcorr") > 0 {
                local mylabel = "Spatial Correlation Precipitation"
            }
            else if strpos("`name'", "tmpspatcorr") > 0 {
                local mylabel = "Spatial Correlation spatcorr Temperature"
            }
            else if strpos("`name'", "tmnspatcorr") > 0 {
                local mylabel = "Spatial Correlation Minimum Temperature"
            }
            else if strpos("`name'", "tmxspatcorr") > 0 {
                local mylabel = "Spatial Correlation Maximum Temperature"
            }
            else if strpos("`name'", "dtrspatcorr") > 0 {
                local mylabel = "Spatial Correlation Diurnal Temperature Range"
            }
            else if strpos("`name'", "rehspatcorr") > 0 {
                local mylabel = "Spatial Correlation Relative Humidity (%)"
            }
            else if strpos("`name'", "petspatcorr") > 0 {
                local mylabel = "Spatial Correlation Potential Evapotranspiration (mm/month)"
            }
            else if strpos("`name'", "frsspatcorr") > 0 {
                local mylabel = "Spatial Correlation Frequency of Frost Days (days/month)"
            }
            else if strpos("`name'", "vapspatcorr") > 0 {
                local mylabel = "Spatial Correlation Vapor Pressure (hPa)"
            }
            else if strpos("`name'", "wetspatcorr") > 0 {
                local mylabel = "Spatial Correlation Wet Days (days/month)"
            }
            else if strpos("`name'", "cldvolatility") > 0 {
                local mylabel = "Volatility Cloud Cover (%)"
            }
            else if strpos("`name'", "prevolatility") > 0 {
                local mylabel = "Volatility Precipitation"
            }
            else if strpos("`name'", "tmpvolatility") > 0 {
                local mylabel = "Volatility volatility Temperature"
            }
            else if strpos("`name'", "tmnvolatility") > 0 {
                local mylabel = "Volatility Minimum Temperature"
            }
            else if strpos("`name'", "tmxvolatility") > 0 {
                local mylabel = "Volatility Maximum Temperature"
            }
            else if strpos("`name'", "dtrvolatility") > 0 {
                local mylabel = "Volatility Diurnal Temperature Range"
            }
            else if strpos("`name'", "rehvolatility") > 0 {
                local mylabel = "Volatility Relative Humidity (%)"
            }
            else if strpos("`name'", "petvolatility") > 0 {
                local mylabel = "Volatility Potential Evapotranspiration (mm/month)"
            }
            else if strpos("`name'", "frsvolatility") > 0 {
                local mylabel = "Volatility Frequency of Frost Days (days/month)"
            }
            else if strpos("`name'", "vapvolatility") > 0 {
                local mylabel = "Volatility Vapor Pressure (hPa)"
            }
            else if strpos("`name'", "wetvolatility") > 0 {
                local mylabel = "Volatility Wet Days (days/month)"
            }
            else if strpos("`name'", "CSI") > 0 {
                local ending = substr("`name'", length("`name'")-4, 5)
                if strpos("`ending'", "med")>0{
                    local mypos = length("`name'")-3
                }
                else{
                    local mypos = length("`name'")-2
                }
                local crop = strtitle(substr("`name'", 4, `mypos'-3))
                local crop2 = substr("`name'", 4, `mypos'-3)
                local mylabel = "Caloric Suitability `crop'"
                capture label var `crop2'lo`stat' "`crop' Yield (`stat')"
                capture label var `crop2'med`stat' "`crop' Yield (`stat')"
                capture label var `crop2'hi`stat' "`crop' Yield (`stat')"
            }
            else if strpos("`name'", "rix") > 0 {
                local mylabel = "Ruggedness"
            }
            else if strpos("`name'", "etopo") > 0 {
                local mylabel = "Elevation"
            }
            else if strpos("`name'", "globe") > 0 {
                local mylabel = "Elevation"
            }
            else if strpos("`name'", "HLD") > 0 {
                local mylabel = "Harmonized Lights"
            }
            else if substr("`name'", 1, 2) == "F1" {
                local year = substr("`name'", 4, 4)
                local mylabel = "Lights (`year')"
            }
            else if substr("`name'", 1, 2) == "ls" {
                local mylabel = "Population Count (Landscan)"
            }
            else if substr("`name'", 1, 6) == "afpopd" {
                local mylabel = "Population Density (UNEP)"
            }
            else if substr("`name'", 1, 9) == "gpwv4popd" {
                local year = substr("`name'", 10, 13)
                local mylabel = "Population Density (`year', GPWv4)"
            }
            else if substr("`name'", 1, 9) == "gpwv4popc" {
                local year = substr("`name'", 10, 13)
                local mylabel = "Population Count (`year', GPWv4)"
            }
            else if substr("`name'", 1, 4) == "popc" {
                local mylabel = "Population Count (HYDE)"
            }
            else if substr("`name'", 1, 4) == "popd" {
                local mylabel = "Population Density (HYDE)"
            }
            else if substr("`name'", 1, 4) == "rurc" {
                local mylabel = "Rural Population Count (HYDE)"
            }
            else if substr("`name'", 1, 4) == "urbc" {
                local mylabel = "Urban Population Count (HYDE)"
            }
            else if strpos("`name'", "stxv") > 0 {
                local mylabel = "Malaria Suitability"
            }
            else if strpos("`name'", "tse") > 0 {
                local mylabel = "Tse-tse Suitability"
            }
            else if strpos("`name'", "Fusca") > 0 {
                local mylabel = "(Fusca) Tse-tse Suitability"
            }
            else if strpos("`name'", "Morsitans") > 0 {
                local mylabel = "(Morsitans) Tse-tse Suitability"
            }
            else if strpos("`name'", "Palpalis") > 0 {
                local mylabel = "(Palpalis) Tse-tse Suitability"
            }
            else if strpos("`name'", "climfac") > 0 {
                local mylabel = "Agricultural Suitability (Climatic)"
            }
            else if strpos("`name'", "climsuit") > 0 {
                local mylabel = "Agricultural Suitability (Climatic)"
            }
            else if strpos("`name'", "soilfac") > 0 {
                local mylabel = "Agricultural Suitability (Soil)"
            }
            else if strpos("`name'", "suit") > 0 {
                local mylabel = "Agricultural Suitability (Ramankutty et al.)"
            }
            else if strpos("`name'", "ecodiversity") > 0 {
                local mylabel = "Ecological Diversity"
            }
            else if strpos("`name'", "ecopolarization") > 0 {
                local mylabel = "Ecological Polarization"
            }
            else if strpos("`name'", "ecodiversityLGM") > 0 {
                local mylabel = "Ecological Diversity (LGM)"
            }
            else if strpos("`name'", "ecopolarizationLGM") > 0 {
                local mylabel = "Ecological Polarization (LGM)"
            }
            else if "`name'"=="sea100" {
                local mylabel = "Land Area within 100km of Sea"
            }
            else if "`name'"=="sea100" {
                local mylabel = "Share of Land Area within 100km of Sea"
            }

            if length("`var'") >= 30 {
                local var2 = subinstr("`var'", "Calories", "Cal", .)
                local var2 = subinstr("`var2'", "Total", "Tot", .)
                rename `var' `var2'
                label var `var2' "`mylabel' (`stat')"
            }
            else {
                label var `var' "`mylabel' (`stat')"
            }
        }
    }
}
compress
save "Country-Stats.dta", replace
