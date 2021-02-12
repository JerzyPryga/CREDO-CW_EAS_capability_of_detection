# Analysis of capability of detection of Extensive Air Showers (EAS) using simple scintillator dotectors
This repository contains all files connected to my work at CREDO project. It means programs for analysis and results estimation, plots, all recent progress etc. It is always a "work in progress" project but I hope that as time goes it will become more developed, more oragnised and more interesting :).

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


