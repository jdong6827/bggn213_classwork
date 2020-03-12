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

    ##      Province.State Country.Region     Lat     Long       date cases      type
    ## 2675         Shanxi Mainland China 37.5777 112.2922 2020-03-03     5 recovered
    ## 2676        Sichuan Mainland China 30.6171 102.7103 2020-03-03     8 recovered
    ## 2677        Tianjin Mainland China 39.3054 117.3230 2020-03-03    13 recovered
    ## 2678       Xinjiang Mainland China 41.1129  85.2401 2020-03-03     2 recovered
    ## 2679         Yunnan Mainland China 24.9740 101.4870 2020-03-03     1 recovered
    ## 2680       Zhejiang Mainland China 29.1832 120.0934 2020-03-03    24 recovered

> Q1. How many total infected cases are there around the world?

``` r
sum(virus$cases)
```

    ## [1] 144233

``` r
table(virus$type)
```

    ## 
    ## confirmed     death recovered 
    ##      1461       194      1025

> Q2. How many deaths linke to infected cases have there been?

``` r
inds<- virus$type=="death"
sum(virus$cases[inds]) # sum(virus[inds,"case"])
```

    ## [1] 3160

> Q3. What is the overall death rate?

``` r
sum(virus$case[inds])/sum(virus$case)*100
```

    ## [1] 2.190899

> Q4. What is the death rate in “Mainland China”?

``` r
country <- virus$Country.Region=="Mainland China"

total_case <- sum(virus$cases[country])

country_death <- virus$type=="death" & virus$Country.Region=="Mainland China"

total_death <- sum(virus$cases[country_death])

death_rate <- total_death/total_case*100

round(death_rate,2)
```

    ## [1] 2.26

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

    ## [1] 2.88

``` r
covidDR("Iran")
```

    ## [1] 2.85

``` r
covidDR("US")
```

    ## [1] 5.11
