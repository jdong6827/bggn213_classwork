Class 17 COVID-19 Homework
================

## Coronavirus

Here we analyze infection data for the 2019 novel Coronavirus COVID-19
(2019-nCoV) epidemic. The raw data is pulled from the Johns Hopkins
University Center for Systems Science and Engineering (JHU CCSE)
Coronavirus repository.

A CSV file is available here
<https://github.com/RamiKrispin/coronavirus-csv>

``` r
url <- "https://tinyurl.com/COVID-2019"
virus <- read.csv(url)
tail(virus)
```

    ##      Province.State      Country.Region     Lat     Long       date cases
    ## 3622       Shanghai      Mainland China 31.2020 121.4491 2020-03-10     4
    ## 3623         Shanxi      Mainland China 37.5777 112.2922 2020-03-10     4
    ## 3624        Sichuan      Mainland China 30.6171 102.7103 2020-03-10    12
    ## 3625         Taiwan Taipei and environs 23.7000 121.0000 2020-03-10     2
    ## 3626        Tianjin      Mainland China 39.3054 117.3230 2020-03-10     1
    ## 3627       Zhejiang      Mainland China 29.1832 120.0934 2020-03-10    15
    ##           type
    ## 3622 recovered
    ## 3623 recovered
    ## 3624 recovered
    ## 3625 recovered
    ## 3626 recovered
    ## 3627 recovered

> Q1. How many total infected cases are there around the world?

``` r
sum(virus$cases)
```

    ## [1] 187075

``` r
table(virus$type)
```

    ## 
    ## confirmed     death recovered 
    ##      2112       274      1241

> Q2. How many deaths linke to infected cases have there been?

``` r
inds<- virus$type=="death"
sum(virus$cases[inds]) # sum(virus[inds,"case"])
```

    ## [1] 4262

> Q3. What is the overall death rate?

``` r
sum(virus$case[inds])/sum(virus$case)*100
```

    ## [1] 2.278231

> Q4. What is the death rate in “Mainland China”?

``` r
country <- virus$Country.Region=="Mainland China"

total_case <- sum(virus$cases[country])

country_death <- virus$type=="death" & virus$Country.Region=="Mainland China"

total_death <- sum(virus$cases[country_death])

death_rate <- total_death/total_case*100

round(death_rate,2)
```

    ## [1] 2.18

> Q5. What is the death rate in Italy, Iran and the US?

``` r
# Define a function to calculate death rate for countries

covidDR<-function(x){country <- virus$Country.Region==x

total_case <- sum(virus$cases[country])

country_death <- virus$type=="death" & virus$Country.Region==x

total_death <- sum(virus$cases[country_death])

death_rate <- total_death/total_case*100

round(death_rate,2)
}
```

``` r
# Calculate the death rates for Italy, Iran and the US
covidDR("Italy")
```

    ## [1] 5.49

``` r
covidDR("Iran")
```

    ## [1] NaN

``` r
covidDR("US")
```

    ## [1] 3.45
