# Analysis of capability of detection of Extensive Air Showers (EAS) using simple scintillator dotectors
This repository contains all files connected to my work at CREDO project. It means programs for analysis and results estimation, plots, all recent progress etc. Below you will find a brief description of everything that has been done so far. However, for more detailed informations or suggestions do not hasitate to contact me directly. It is always a "work in progress" project but I hope that as time goes it will become more developed, more oragnised and more interesting :).

# Introduction

One of the main objectives of the CREDO project is to search for so-called Cosmic-Ray Ensembles (CRE). To confirm the existence of such phenomena a massive scale observation of even relatively low energy Extensive Air Showers (EAS) and an analysis of their correlations in time must be performed. To make such observations possible, an infrastructure of widely spread detectors connected in a global network should be developed using low-cost devices capable of collecting data for a long period of time. They should detect even relatively low energy EAS with high reliability and with a time accuracy as high as possible. A candidate for such a network is a system of several Cosmic Watches, small scintillator detectors, connected in a coincidence circuit. Determination of the probability of detection of an EAS is an important part in the analysis of the data collected by a variety of systems. The standard approach based on detailed and extensive simulations is not possible for many systems, thus a faster method is developed. Knowing the characteristics of EAS from more general simulations any required probability is calculated using appropriate parameterization taking into account EAS spectrum, energy dependence of particle density, zenith angle dependence and many others. This allows to estimate expected number of EAS events measured by a set of small detectors. Comparing these results with number of expected background events one can find if the system can reliably detect EAS.

# Metodology

First step in the analysis is to perform some shower simulations, which are the main source of information about EAS in this work. Those data is later analysed to find various properties of the cascades that can be described in a quantitative way by mathematical functions.

Energy spectrum of currently used simulations:
- From 1 TeV to 4000 TeV,
- Proton as primary particle,
- 18 different energies (evenly distributet in logarithmic scale),
- From 1000 up to 20 000 simulated cascades for each energy.

Currently used simulations for angle distibution analysis:
- One fixed energy: 100 TeV,
- Proton as primary particle,
- 7 different angles (up to 70 degrees every 10 degrees).

Before starting any fitting or calculations some assumptions about EAS has to be taken. Current list of those is as follows:
- Showers are circularly simetrical.
- All showers are like produced by protons. 
- Primary particles energy spectrum and frequancies are known.
- Time blur of particles from the cascade is close to none.

As w now the rules we must follow, the next step is to extract information from simulations and turn them into set of rules formuleted as mathematical funtions. It is nessesary to do so because we characterize every EAS using its particle density function \rho. In this work it is a function of many parameters such as:

- E - energy of primary particle, 
- r - distance from its center,
- \theta - primary particle approach direction zenith angle,
- N_part - total number of produced particles,
- maybe more in the future...

Final \rho function is a combination of many factors that depend on previously written parameters. The construction of particle density funtion was meant to be eimple and intuitive, and connetion between each factor and its physical meaning clear. All of those components were found by fitting some pre-defined functions to the data from simulations.

- \rho(r) - The "basis", most important factor. It represents relation between particle density and distance from the center of the shower. This function is fitted for certain vertical cascade (it is obviously an average from many cascades) of chosen energy to which other foctors are normalised. 
- F(E,r) - This factor scales the density and modifies the "shape" (reletion with r) as energy change.
- F(N_part, r) - This factor scales the density and modifies the "shape" (reletion with r) as total number of produced particles changes. Here, the total number of produced particles is strongly connected with the altitude at which the cascade started to form. However, N_part is a more convenient parameter to work with and you can always imagine "smaller" and "bigger" showers.
- F(\theta) - This factor scales the density as zenith angle of primary particle approach direction changes. Here a simplification was made. Namely, the shape of shower footprint is not sphericaly simetrical as this angle increase. However, EAS comes from all directions with equall probability and number of cascades that should occur in the time of measurement is huge. Thus, neglecting this effect is perfectly justified.




