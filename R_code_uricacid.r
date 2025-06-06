## Purpose: This script runs descriptive statistics and statistical tests including Fisher's exact test, Kruskal-Wallis test, analysis of covariance, and piecewise linear regression model using jmv packages for R.
## Author: Naoki Omori                                             
## naoki.shimane.medical@gmail.com                                      
## Department of Neurology, Shimane University
## Izumo-shi, Shimane, Japan

#. Install and load packages

install.packages("jmv")
library("jmv")

#. Load data

data = read.csv ("uric_acid_data.csv", header = TRUE, sep = ",") ##choose your own directory

# Descriptive statistics
jmv::descriptives(
    formula = BMI ~ UA3,
    data = data)

jmv::descriptives(
    formula = age + sex + BMI + HT + DM + DL + eGFR + `etat clibre` ~ UA3,
    data = data)   

#. Fisher's exact test
jmv::contTables(
    formula = ~ sex:UA3,
    data = data,
    chiSq = FALSE,
    fisher = TRUE)

jmv::contTables(
    formula = ~ HT:UA3,
    data = data,
    chiSq = FALSE,
    fisher = TRUE)

jmv::contTables(
    formula = ~ DM:UA3,
    data = data,
    chiSq = FALSE,
    fisher = TRUE)

jmv::contTables(
    formula = ~ DL:UA3,
    data = data,
    chiSq = FALSE,
    fisher = TRUE)

jmv::contTables(
    formula = ~ `etat clibre`:UA3,
    data = data,
    chiSq = FALSE,
    fisher = TRUE)

#. Kruskal-Wallis test
jmv::anovaNP(
    formula = age ~ UA3,
    data = data)

jmv::anovaNP(
    formula = BMI ~ UA3,
    data = data)

jmv::anovaNP(
    formula = eGFR ~ UA3,
    data = data)

#. Analysis of covariance
jmv::ancova(
    formula = `CaudateHead-left_roi` ~ UA3 + age + sex + BMI + eGFR + HT + DM + DL + `etat clibre`,
    data = data,
    postHoc = ~ UA3)

jmv::ancova(
    formula = `CaudateHead-right_roi` ~ UA3 + age + sex + BMI + eGFR + HT + DM + DL + `etat clibre`,
    data = data,
    postHoc = ~ UA3)

jmv::ancova(
    formula = `CaudateTail-left_roi` ~ UA3 + age + sex + BMI + eGFR + HT + DM + DL + `etat clibre`,
    data = data,
    postHoc = ~ UA3)

jmv::ancova(
    formula = `CaudateTail-right_roi` ~ UA3 + age + sex + BMI + eGFR + HT + DM + DL + `etat clibre`,
    data = data,
    postHoc = ~ UA3)

jmv::ancova(
    formula = `LateralGlobusPallidus-left_roi` ~ UA3 + age + sex + BMI + eGFR + HT + DM + DL + `etat clibre`,
    data = data,
    postHoc = ~ UA3)

jmv::ancova(
    formula = `LateralGlobusPallidus-right_roi` ~ UA3 + age + sex + BMI + eGFR + HT + DM + DL + `etat clibre`,
    data = data,
    postHoc = ~ UA3)

jmv::ancova(
    formula = `MedicalGlobusPallidus-left_roi` ~ UA3 + age + sex + BMI + eGFR + HT + DM + DL + `etat clibre`,
    data = data,
    postHoc = ~ UA3)

jmv::ancova(
    formula = `MedicalGlobusPallidus-right_roi` ~ UA3 + age + sex + BMI + eGFR + HT + DM + DL + `etat clibre`,
    data = data,
    postHoc = ~ UA3)

jmv::ancova(
    formula = `Putamen-left_roi` ~ UA3 + age + sex + BMI + eGFR + HT + DM + DL + `etat clibre`,
    data = data,
    postHoc = ~ UA3)

jmv::ancova(
    formula = `Putamen-right_roi` ~ UA3 + age + sex + BMI + eGFR + HT + DM + DL + `etat clibre`,
    data = data,
    postHoc = ~ UA3)

jmv::ancova(
    formula = `SubstantiaNigra-left_roi` ~ UA3 + age + sex + BMI + eGFR + HT + DM + DL + `etat clibre`,
    data = data,
    postHoc = ~ UA3)

jmv::ancova(
    formula = `SubstantiaNigra-right_roi` ~ UA3 + age + sex + BMI + eGFR + HT + DM + DL + `etat clibre`,
    data = data,
    postHoc = ~ UA3)

jmv::ancova(
    formula = `SubthalamicNucleus-left_roi` ~ UA3 + age + sex + BMI + eGFR + HT + DM + DL + `etat clibre`,
    data = data,
    postHoc = ~ UA3)

jmv::ancova(
    formula = `SubthalamicNucleus-right_roi` ~ UA3 + age + sex + BMI + eGFR + HT + DM + DL + `etat clibre`,
    data = data,
    postHoc = ~ UA3)

#. Piecewise linear regression
jmv::linReg(
    data = data,
    dep = CaudateHead-left_roi,
    covs = vars(UA, age, BMI, eGFR),
    factors = vars(sex, HT, DM, DL),
    blocks = list(
        list(
            "UA",
            "age",
            "sex",
            "BMI",
            "eGFR",
            "HT",
            "DM",
            "DL")),
    refLevels = list(
        list(
            var="sex",
            ref="0"),
        list(
            var="HT",
            ref="0"),
        list(
            var="DM",
            ref="0"),
        list(
            var="DL",
            ref="0")))

jmv::linReg(
    data = data,
    dep = CaudateHead-right_roi,
    covs = vars(UA, age, BMI, eGFR),
    factors = vars(sex, HT, DM, DL),
    blocks = list(
        list(
            "UA",
            "age",
            "sex",
            "BMI",
            "eGFR",
            "HT",
            "DM",
            "DL")),
    refLevels = list(
        list(
            var="sex",
            ref="0"),
        list(
            var="HT",
            ref="0"),
        list(
            var="DM",
            ref="0"),
        list(
            var="DL",
            ref="0")))

jmv::linReg(
    data = data,
    dep = Putamen-left_roi,
    covs = vars(UA, age, BMI, eGFR),
    factors = vars(sex, HT, DM, DL),
    blocks = list(
        list(
            "UA",
            "age",
            "sex",
            "BMI",
            "eGFR",
            "HT",
            "DM",
            "DL")),
    refLevels = list(
        list(
            var="sex",
            ref="0"),
        list(
            var="HT",
            ref="0"),
        list(
            var="DM",
            ref="0"),
        list(
            var="DL",
            ref="0")))

jmv::linReg(
    data = data,
    dep = Putamen-right_roi,
    covs = vars(UA, age, BMI, eGFR),
    factors = vars(sex, HT, DM, DL),
    blocks = list(
        list(
            "UA",
            "age",
            "sex",
            "BMI",
            "eGFR",
            "HT",
            "DM",
            "DL")),
    refLevels = list(
        list(
            var="sex",
            ref="0"),
        list(
            var="HT",
            ref="0"),
        list(
            var="DM",
            ref="0"),
        list(
            var="DL",
            ref="0")))
