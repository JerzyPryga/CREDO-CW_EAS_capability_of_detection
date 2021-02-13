# Analysis of capability of detection of Extensive Air Showers (EAS) using simple scintillator dotectors
This repository contains all files connected to my work at CREDO project. It means programs for analysis and results estimation, plots, all of current results and recent progress etc. Below you will find a brief description of everything that has been done so far. However, for more detailed informations or suggestions do not hasitate to contact me directly. It is always a "work in progress" project but I hope that as time goes it will become more developed, more oragnised and more interesting :).

# Abstract

One of the main objectives of the CREDO project is to search for so-called Cosmic-Ray Ensembles (CRE). To confirm the existence of such phenomena a massive scale observation of even relatively low energy Extensive Air Showers (EAS) and an analysis of their correlations in time must be performed. To make such observations possible, an infrastructure of widely spread detectors connected in a global network should be developed using low-cost devices capable of collecting data for a long period of time. They should detect even relatively low energy EAS with high reliability and with a time accuracy as high as possible. A candidate for such a network is a system of several Cosmic Watches, small scintillator detectors, connected in a coincidence circuit. Determination of the probability of detection of an EAS is an important part in the analysis of the data collected by a variety of systems. The standard approach based on detailed and extensive simulations is not possible for many systems, thus a faster method is developed. Knowing the characteristics of EAS from more general simulations any required probability is calculated using appropriate parameterization taking into account EAS spectrum, energy dependence of particle density, zenith angle dependence and many others. This allows to estimate expected number of EAS events measured by a set of small detectors. Comparing these results with number of expected background events one can find if the system can reliably detect EAS.

# Idea and method

An ideal of how to detect EAS which just few small devices is very simple. As the cosmic-rays showers manifest themselves (among other effects in the atmosphere that require completely different methods of detection) as incresed particles density on the groud that last for a very short period in time. Thus, one can imagine that when there are several detectors placed close to each other and a cascade occuring near such system. At that moment, several devices may detect particles that ariginate form this exact EAS and give signals at almoast the same time. Probability of such situation to be cause by two uncorrelated cosmic rays or other effects should be much smaller so both types of events could be distinguished. It is possible to evaluate expected number of both cascades and background events using current knowledge about EAS as well as simulations of them and information about detectors in such system. This is the exact purpose of this work. Obtained results should be very useful for CREDO project as it will give some information like:
- How to interpret data gathered by such systems.
- If certain system works as it should.
- What are pros and cons of such systems.
- Some hints how improve the design such systems. 

# Analysis

First step in the analysis is to perform some shower simulations, which are the main source of information about EAS in this work. All simulations was performed using CORSIKA software. Those data is later analysed to find various properties of the cascades that can be described in a quantitative way by mathematical functions.

Energy spectrum of currently used simulations:
- From 1 TeV to 4000 TeV,
- Proton as primary particle,
- 18 different energies (evenly distributet in logarithmic scale),
- From 1000 up to 20 000 simulated cascades for each energy.

Currently used simulations for angle distibution analysis:
- One fixed energy: 100 TeV,
- Proton as primary particle,
- 7 different angles (up to 70 degrees every 10 degrees).

Staritng point in the analysis is to evaluate flux of bacground particles, and probability of fake signals i.e. no caused by EAS. In this work, term "background" means uncorelated, single cosmic-rays which comes from all direction all the time. Also natural radiation and noise from electronics itself my give signals. Taken assumptions are as follows:
- Only certain types of particles with specific energies gives signal.
- Background flux is known and constant in time.

Next step is to calculate probability of a fake signal to happen. Thus, again some assumptions about used devices has to be taken:
- All detectors in the system are identical and behave the same way.
- Detectors behaviour doeas not change in time.
- Particles can be detected from all directions.
- Only certain types of particles with specific energies gives signal.
- Interval in which signals from different detectors are treated as coincidance is \delta T.

Time \delta T is a unit of time in which we calculate the probability of signal to happen. Below formula shows how this is done:

\\wzór

- \delta T - coincidence time [s].
- \ni - detectors efficiency [%].
- A - surface area of the detector [cm^2].
- I_bg - background particles flux [1/s cm^2].
- f_bg - frequancy of non cosmic background signals [1/s].

After calculating this, chance of coincidance signals is computed as well as expected number of such events in the time of measurement. Results of these calculations are presented later in table and are compared with expected number of signals caused by the background.

Before starting any fitting or calculations some assumptions about EAS has to be taken. Current list of crucial ones is as follows:
- Showers are circularly simetrical.
- All showers are like produced by protons. 
- Primary particles energy spectrum and frequancies are known and constant in time.
- Time blur of particles from the cascade is close to none.

As w now the rules we must follow, the next step is to extract information from simulations and turn them into set of rules formuleted as mathematical funtions. It is nessesary to do so because we characterize every EAS using its particle density function \rho. In this work it is a function of many parameters such as:

- E - energy of primary particle [TeV], 
- r - distance from its center [m],
- \theta - primary particle approach direction zenith angle [radians],
- N_part - total number of produced particles,
- maybe more in the future...

Final \rho function is a combination of many factors that depend on previously written parameters. The construction of particle density funtion was meant to be eimple and intuitive, and connetion between each factor and its physical meaning clear. All of those components were found by fitting some pre-defined functions to the data from simulations.

\\wzór

- \rho(r) - The "basis", most important factor. It represents relation between particle density and distance from the center of the shower. This function is fitted for certain vertical cascade (it is obviously an average from many cascades) of chosen energy to which other foctors are normalised. 
- F(E,r) - This factor scales the density and modifies the "shape" (reletion with r) as energy change.
- F(N_part, r) - This factor scales the density and modifies the "shape" (reletion with r) as total number of produced particles changes. Here, the total number of produced particles is strongly connected with the altitude at which the cascade started to form. However, N_part is a more convenient parameter to work with and you can always imagine "smaller" and "bigger" showers.
- F(\theta) - This factor scales the density as zenith angle of primary particle approach direction changes. Here a simplification was made. Namely, the shape of shower footprint is not sphericaly simetrical as this angle increase. However, EAS comes from all directions with equall probability and number of cascades that should occur in the time of measurement is huge. Thus, neglecting this effect is perfectly justified.

Next step is to calculate probability of a signal to happen when there is a cascade near the detector. Assumptions about detectors are the same as described earlier. Formula for single signal probability is as follows:

\\wzór

To evaluate expected number of events coused by EAS on integration over all energies of primary particles, area around the devices and sperical angle of the sky
must be performed. Formula below shows how it is computed (notice that with assumption that cascades are circulaly symetrical reduces integral ofver horizontal angle to the factor of 2 \pi r):

\\wzór

- Q(n,k,P) - probability of k signals in system of n detectors.
- r - distance to the center of the shower [m].
- j(E) - cascades frequancy function [1/s m^2].
- E - energy of primary particle [TeV].
- T - time of measurement [s].
- \Omega - spherical angle [steradians].

# Comparison with simpler method

To check if the performed analysis was done right another approach to the problem was tested. As muons have been studied for many years, an approximate function for its density on the ground level can be found in the literature. It depends on total number of produced muons which is a function of the energy of primary particle and distance from the center of the shower. So the parameters taht charecterise the showers are the same as in analysis described earlier. In this work it was only modified by the scale factor of primary particle zenith angle. Thus, used formula is as follows:

//wzór

- r - distance to the center of the shower [m].
- F(\theta) - This factor scales the density as zenith angle of primary particle approach direction changes (as described in formula []).
- N_part - Total umber of produced muons as a function of primary particle energy.
- E - energy of primary particle [TeV].

Further steps in the analysis are the same as previously described. This comparison should help to judge if analysis described earlier, which is more general and gives more information about the behaviour of the system, is reasonable and does not consist any unrealistic assumptions.

# Current results

Here are presented results of above analysis calculated for exemplary system. Properties of the system are listed below:
- 4 devices (Cosmic Watches).
- 

# Plans for the future
