# Analysis of capability of detection of Extensive Air Showers (EAS) using simple scintillator detectors
This repository contains all files connected to my work at CREDO project. It means programs for analysis and results estimation, plots, all of current results and recent progress etc. Below you will find a brief description of everything that has been done so far. However, for more detailed information or suggestions do not hesitate to contact me directly. It is always a "work in progress" project but I hope that as time goes it will become more developed, more organised and more interesting :).

## Abstract

One of the main objectives of the CREDO project is to search for so-called Cosmic-Ray Ensembles (CRE). To confirm the existence of such phenomena a massive scale observation of even relatively low energy Extensive Air Showers (EAS) and an analysis of their correlations in time must be performed. To make such observations possible, an infrastructure of widely spread detectors connected in a global network should be developed using low-cost devices capable of collecting data for a long period of time. They should detect even relatively low energy EAS with high reliability and with a time accuracy as high as possible. A candidate for such a network is a system of several Cosmic Watches, small scintillator detectors, connected in a coincidence circuit. Determination of the probability of detection of an EAS is an important part in the analysis of the data collected by a variety of systems. The standard approach based on detailed and extensive simulations is not possible for many systems, thus a faster method is developed. Knowing the characteristics of EAS from more general simulations any required probability is calculated using appropriate parameterisation taking into account EAS spectrum, energy dependence of particle density, zenith angle dependence and many others. This allows to estimate expected number of EAS events measured by a set of small detectors. Comparing these results with number of expected background events one can find if the system can reliably detect EAS.

## Idea and method

An ideal of how to detect EAS which just few small devices is very simple. As the cosmic-rays showers manifest themselves (among other effects in the atmosphere that require completely different methods of detection) as increased particles density on the ground that last for a very short period in time. Thus, one can imagine that when there are several detectors placed close to each other and a cascade occurring near such system. At that moment, several devices may detect particles that originate form this exact EAS and give signals at almost the same time. Probability of such situation to be cause by two uncorrelated cosmic rays or other effects should be much smaller so both types of events could be distinguished. It is possible to evaluate expected number of both cascades and background events using current knowledge about EAS as well as simulations of them and information about detectors in such system. This is the exact purpose of this work. Obtained results should be very useful for CREDO project as it will give some information like:
- How to interpret data gathered by such systems.
- If certain system works as it should.
- What are pros and cons of such systems.
- Some hints how improve the design such systems. 

## Analysis

#### Simulations
First step in the analysis is to perform some shower simulations, which are the main source of information about EAS in this work. All simulations was performed using CORSIKA software. Those data is later analysed to find various properties of the cascades that can be described in a quantitative way by mathematical functions.

Energy spectrum of currently used simulations:
- From 1 TeV to 4000 TeV,
- Proton as primary particle,
- 18 different energies (evenly distributed in logarithmic scale),
- From 1000 up to 20 000 simulated cascades for each energy.

Currently used simulations for angle distribution analysis:
- One fixed energy: 100 TeV,
- Proton as primary particle,
- 7 different angles (up to <a href="https://www.codecogs.com/eqnedit.php?latex=70^o" target="_blank"><img src="https://latex.codecogs.com/gif.latex?70^o" title="70^o" /></a> with a step of <a href="https://www.codecogs.com/eqnedit.php?latex=10^o" target="_blank"><img src="https://latex.codecogs.com/gif.latex?10^o" title="10^o" /></a>).

#### Background
Starting point in the analysis is to evaluate flux of background particles, and probability of fake signals i.e. no caused by EAS. In this work, term "background" means uncorrelated, single cosmic-rays which comes from all direction all the time. Also natural radiation and noise from electronics itself my give signals. Taken assumptions are as follows:
- Only certain types of particles with specific energies gives signal.
- Background flux is known and constant in time.

Next step is to calculate probability of a fake signal to happen. Thus, again some assumptions about used devices has to be taken:
- All detectors in the system are identical and behave the same way.
- Detectors behaviour does not change in time.
- Particles can be detected from all directions.
- Only certain types of particles with specific energies gives signal.
- Interval in which signals from different detectors are treated as coincidence is <a href="https://www.codecogs.com/eqnedit.php?latex=\delta&space;T" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\delta&space;T" title="\delta T" /></a>.

Time <a href="https://www.codecogs.com/eqnedit.php?latex=\delta&space;T" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\delta&space;T" title="\delta T" /></a> is a unit of time in which we calculate the probability of signal to happen. Below formula shows how this is done:

<a href="https://www.codecogs.com/eqnedit.php?latex=P_{bg}&space;=&space;1&space;-&space;exp(\delta&space;T&space;(\eta&space;\cdot&space;A&space;\cdot&space;I_{bg})&plus;f_{bg})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?P_{bg}&space;=&space;1&space;-&space;exp(\delta&space;T&space;(\eta&space;\cdot&space;A&space;\cdot&space;I_{bg})&plus;f_{bg})" title="P_{bg} = 1 - exp(\delta T (\eta \cdot A \cdot I_{bg})+f_{bg})" /></a>

- <a href="https://www.codecogs.com/eqnedit.php?latex=\delta&space;T" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\delta&space;T" title="\delta T" /></a> - coincidence time [s].
- <a href="https://www.codecogs.com/eqnedit.php?latex=\eta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\eta" title="\eta" /></a> - detectors efficiency [%].
- <a href="https://www.codecogs.com/eqnedit.php?latex=A" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A" title="A" /></a> - surface area of the detector [<a href="https://www.codecogs.com/eqnedit.php?latex=cm^2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?cm^2" title="cm^2" /></a>].
- <a href="https://www.codecogs.com/eqnedit.php?latex=I_{bg}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?I_{bg}" title="I_{bg}" /></a> - background particles flux <a href="https://www.codecogs.com/eqnedit.php?latex=[\frac{1}{cm^2s}]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?[\frac{1}{cm^2s}]" title="[\frac{1}{cm^2s}]" /></a>.
- <a href="https://www.codecogs.com/eqnedit.php?latex=f_{bg}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f_{bg}" title="f_{bg}" /></a> - frequency of non cosmic background signals <a href="https://www.codecogs.com/eqnedit.php?latex=[\frac{1}{s}]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?[\frac{1}{s}]" title="[\frac{1}{s}]" /></a>.

After calculating this, chance of coincidence signals is computed as well as expected number of such events in the time of measurement using following formula:

<a href="https://www.codecogs.com/eqnedit.php?latex=\langle&space;N_{bg}(k)&space;\rangle&space;=&space;Q(n,k,P_{bg})&space;\cdot&space;\frac{T}{\delta&space;T}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\langle&space;N_{bg}(k)&space;\rangle&space;=&space;Q(n,k,P_{bg})&space;\cdot&space;\frac{T}{\delta&space;T}" title="\langle N_{bg}(k) \rangle = Q(n,k,P_{bg}) \cdot \frac{T}{\delta T}" /></a>

Where:
- <a href="https://www.codecogs.com/eqnedit.php?latex=Q(n,k,P_{bg)})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Q(n,k,P_{bg)})" title="Q(n,k,P_{bg)})" /></a> - probability of k signals in system of n detectors (given by binomial distribution).
- T - time of measurement [s].
- <a href="https://www.codecogs.com/eqnedit.php?latex=\delta&space;T" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\delta&space;T" title="\delta T" /></a> - time of coincidence [s].

Results of those calculations are presented later in table and are compared with expected number of signals caused by the background.

#### Signals from the cascades
Before starting any fitting or calculations some assumptions about EAS has to be taken. Current list of crucial ones is as follows:
- Showers are circularly symmetrical.
- All showers are like produced by protons. 
- Primary particles energy spectrum and frequencies are known and constant in time.
- Time blur of particles from the cascade is close to none.

As w now the rules we must follow, the next step is to extract information from simulations and turn them into set of rules formulated as mathematical functions. It is necessary to do so because we characterise every EAS using its particle density function <a href="https://www.codecogs.com/eqnedit.php?latex=\rho" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\rho" title="\rho" /></a>. In this work it is a function of many parameters such as:

- E - energy of primary particle [TeV], 
- r - distance from its centre [m],
- <a href="https://www.codecogs.com/eqnedit.php?latex=\theta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta" title="\theta" /></a> - primary particle approach direction zenith angle [radians],
- <a href="https://www.codecogs.com/eqnedit.php?latex=N_{part}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?N_{part}" title="N_{part}" /></a> - total number of produced particles,
- maybe more in the future...

Final <a href="https://www.codecogs.com/eqnedit.php?latex=\rho" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\rho" title="\rho" /></a> function is a combination of many factors that depend on previously written parameters. The construction of particle density function was meant to be simple and intuitive, and connection between each factor and its physical meaning clear. All of those components were found by fitting some pre-defined functions to the data from simulations.

<a href="https://www.codecogs.com/eqnedit.php?latex=\rho(r,&space;\theta,&space;E,&space;N_{part})&space;=&space;\rho_{norm}(r)&space;\cdot&space;F_{\theta}(\theta)&space;\cdot&space;F_E(E,r)&space;\cdot&space;F_N(N_{part},r)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\rho(r,&space;\theta,&space;E,&space;N_{part})&space;=&space;\rho_{norm}(r)&space;\cdot&space;F_{\theta}(\theta)&space;\cdot&space;F_E(E,r)&space;\cdot&space;F_N(N_{part},r)" title="\rho(r, \theta, E, N_{part}) = \rho_{norm}(r) \cdot F_{\theta}(\theta) \cdot F_E(E,r) \cdot F_N(N_{part},r)" /></a>

- <a href="https://www.codecogs.com/eqnedit.php?latex=\rho_{norm}(r)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\rho_{norm}(r)" title="\rho_{norm}(r)" /></a> - The "basis", most important factor. It represents relation between particle density and distance from the centre of the shower. This function is fitted for certain vertical cascade (it is obviously an average from many cascades) of chosen energy to which other factors are normalised. 
- <a href="https://www.codecogs.com/eqnedit.php?latex=F_E(E,r)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?F_E(E,r)" title="F_E(E,r)" /></a> - This factor scales the density and modifies the "shape" (relation with r) as energy change.
- <a href="https://www.codecogs.com/eqnedit.php?latex=F_N(N_{part},r)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?F_N(N_{part},r)" title="F_N(N_{part},r)" /></a> - This factor scales the density and modifies the "shape" (relation with r) as total number of produced particles changes. Here, the total number of produced particles is strongly connected with the altitude at which the cascade started to form. However, <a href="https://www.codecogs.com/eqnedit.php?latex=N_{part}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?N_{part}" title="N_{part}" /></a> is a more convenient parameter to work with and you can always imagine "smaller" and "bigger" showers.
- <a href="https://www.codecogs.com/eqnedit.php?latex=F_{\theta}(\theta)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?F_{\theta}(\theta)" title="F_{\theta}(\theta)" /></a> - This factor scales the density as zenith angle of primary particle approach direction changes. Here a simplification was made. Namely, the shape of shower footprint is not circularly symmetrical as this angle increase. However, EAS comes from all directions with equal probability and number of cascades that should occur in the time of measurement is huge. Thus, neglecting this effect is perfectly justified.

Next step is to calculate probability of a signal to happen when there is a cascade near the detector. Assumptions about detectors are the same as described earlier. Formula for single signal probability is as follows:

<a href="https://www.codecogs.com/eqnedit.php?latex=P&space;=&space;1&space;-&space;exp(\eta&space;\cdot&space;A&space;\cdot&space;\rho(r,&space;\theta,&space;E,&space;N_{part}))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?P&space;=&space;1&space;-&space;exp(\eta&space;\cdot&space;A&space;\cdot&space;\rho(r,&space;\theta,&space;E,&space;N_{part}))" title="P = 1 - exp(\eta \cdot A \cdot \rho(r, \theta, E, N_{part}))" /></a>

- <a href="https://www.codecogs.com/eqnedit.php?latex=\eta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\eta" title="\eta" /></a> - detectors efficiency [%].
- <a href="https://www.codecogs.com/eqnedit.php?latex=A" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A" title="A" /></a> - surface area of the detector [<a href="https://www.codecogs.com/eqnedit.php?latex=cm^2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?cm^2" title="cm^2" /></a>].
- <a href="https://www.codecogs.com/eqnedit.php?latex=\rho(r,&space;\theta,&space;E,&space;N_{part})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\rho(r,&space;\theta,&space;E,&space;N_{part})" title="\rho(r, \theta, E, N_{part})" /></a> - particles density function <a href="https://www.codecogs.com/eqnedit.php?latex=[\frac{1}{cm^2}]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?[\frac{1}{cm^2}]" title="[\frac{1}{cm^2}]" /></a>.

There is also another effect in the cascades that could have impact on the probability of coincidence signals. It is currently under investigation, but early results suggest that when one particle that originates from the shower is detected, then chance that another one reach the ground in close region of such event should increase. This "clustering" effect could be most likely explained by fact, that particles in the shower are produced in at least pairs. Thus, when particles are produced together not very high in the atmosphere they should hit the ground close to each other as they may not have time to travel further. In this work this results are included in calculations. The function which scale the probability of multiple signals is a function of particles density and it impacts coincidence signals of 2 and more.

To evaluate expected number of events caused by EAS on integration over all energies of primary particles, area around the devices and spherical angle of the sky must be performed. Formula below shows how it is computed (notice that with assumption that cascades are circularly symmetrical reduces integral over horizontal angle to the factor of <a href="https://www.codecogs.com/eqnedit.php?latex=2&space;\pi&space;r" target="_blank"><img src="https://latex.codecogs.com/gif.latex?2&space;\pi&space;r" title="2 \pi r" /></a>):

<a href="https://www.codecogs.com/eqnedit.php?latex=\langle&space;N(k)&space;\rangle&space;=&space;\int_{0}^{r_{max}}\int_{E_{min}}^{E_{max}}\int_{0}^{\frac{\pi}{2}}Q(n,k,P)&space;2\pi&space;r&space;j(E)&space;T&space;d\Omega&space;dE&space;dr" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\langle&space;N(k)&space;\rangle&space;=&space;\int_{0}^{r_{max}}\int_{E_{min}}^{E_{max}}\int_{0}^{\frac{\pi}{2}}Q(n,k,P)&space;2\pi&space;r&space;j(E)&space;T&space;d\Omega&space;dE&space;dr" title="\langle N(k) \rangle = \int_{0}^{r_{max}}\int_{E_{min}}^{E_{max}}\int_{0}^{\frac{\pi}{2}}Q(n,k,P) 2\pi r j(E) T d\Omega dE dr" /></a>

- Q(n,k,P) - probability of k signals in system of n detectors (given by binomial distribution).
- r - distance to the centre of the shower [m].
- j(E) - cascades frequency function <a href="https://www.codecogs.com/eqnedit.php?latex=[\frac{1}{m^2s}]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?[\frac{1}{m^2s}]" title="[\frac{1}{m^2s}]" /></a>.
- E - energy of primary particle [TeV].
- T - time of measurement [s].
- <a href="https://www.codecogs.com/eqnedit.php?latex=\Omega" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Omega" title="\Omega" /></a> - spherical angle [steradians].

Of course the integration must be performed in certain limits, thus the primary cosmic-rays energy range and maximal distance from the shower must be chosen. Integration up to pi/2 over zenith angle is obvious as it was assumed that particles from all directions can be detected. Criteria of choice of the maximal distance may be different and in this work two were considered. Namely:
- <a href="https://www.codecogs.com/eqnedit.php?latex=R_{prc}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?R_{prc}" title="R_{prc}" /></a> - Radius in which certain percent (in this case 95% was chosen) of particles produced by the cascade are included.
- <a href="https://www.codecogs.com/eqnedit.php?latex=R_{rho}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?R_{rho}" title="R_{rho}" /></a> - Radius in which shower particles density is greater than the avereage background.
Both of this quantities were analysed and characterised as a function of energy of primary particle.

## Comparison with simpler method

To check if the performed analysis was done right another approach to the problem was tested. As muons have been studied for many years, an approximate function for its density on the ground level can be found in the literature. It depends on total number of produced muons which is a function of the energy of primary particle and distance from the centre of the shower. So the parameters that characterise the showers are the same as in analysis described earlier. In this work it was only modified by the scale factor of primary particle zenith angle. Thus, used formula is as follows:

//wzór

- r - distance to the centre of the shower [m].
- <a href="https://www.codecogs.com/eqnedit.php?latex=F_{\theta}(\theta)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?F_{\theta}(\theta)" title="F_{\theta}(\theta)" /></a> - This factor scales the density as zenith angle of primary particle approach direction changes (as described in formula []).
- <a href="https://www.codecogs.com/eqnedit.php?latex=N_{part}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?N_{part}" title="N_{part}" /></a> - Total umber of produced muons as a function of primary particle energy.
- E - energy of primary particle [TeV].

Further steps in the analysis are the same as previously described. This comparison should help to judge if analysis described earlier, which is more general and gives more information about the behaviour of the system, is reasonable and does not consist any unrealistic assumptions.

## Current results

#### Parameters
Here are presented results of above analysis calculated for exemplary system. Assumed properties of the system are listed below:
- 4 devices (Cosmic Watches).
- Time of coincidence:

<a href="https://www.codecogs.com/eqnedit.php?latex=\delta&space;T&space;=&space;200\&space;ns" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\delta&space;T&space;=&space;200\&space;ns" title="\delta T = 200\ ns" /></a>.
- Area of the surface:

<a href="https://www.codecogs.com/eqnedit.php?latex=A&space;=&space;25\&space;cm^2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A&space;=&space;25\&space;cm^2" title="A = 25\ cm^2" /></a>.
- Efficiency of the detector:

<a href="https://www.codecogs.com/eqnedit.php?latex=\eta&space;=&space;95\&space;%" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\eta&space;=&space;95\&space;%" title="\eta = 95\ %" /></a>.
- Particle that gives signal:

<a href="https://www.codecogs.com/eqnedit.php?latex=muons\&space;with\&space;E&space;\geqslant&space;1&space;GeV" target="_blank"><img src="https://latex.codecogs.com/gif.latex?muons\&space;with\&space;E&space;\geqslant&space;1&space;GeV" title="muons\ with\ E \geqslant 1 GeV" /></a>.
- Time of measurement t = 7 days.

The choice of above parameters are not arbitrary. They are the same as in the system tested by prof. Tadeusz Wibig in his work (only efficiency was arbitrary chosen), however there is no certainty about which particles gives signal in the detector. However, according to information given by designers of Cosmic Watch detector it should be muons. These assumptions and choice of parameters yield one of the following properties of the background:
- Background particles flux:

<a href="https://www.codecogs.com/eqnedit.php?latex=I_{bg}&space;=&space;109.96\&space;[\frac{1}{m^2s}]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?I_{bg}&space;=&space;109.96\&space;[\frac{1}{m^2s}]" title="I_{bg} = 109.96\ [\frac{1}{m^2s}]" /></a>.
- Frequency of non cosmic background signals (arbitrary chosen):

<a href="https://www.codecogs.com/eqnedit.php?latex=f_{bg}&space;=&space;0.1\&space;[\frac{1}{s}]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f_{bg}&space;=&space;0.1\&space;[\frac{1}{s}]" title="f_{bg} = 0.1\ [\frac{1}{s}]" /></a>.

Choice of integration limits:
- Energy range:

<a href="https://www.codecogs.com/eqnedit.php?latex=1\&space;TeV&space;-&space;10^5\&space;TeV" target="_blank"><img src="https://latex.codecogs.com/gif.latex?1\&space;TeV&space;-&space;10^5\&space;TeV" title="1\ TeV - 10^5\ TeV" /></a>.
- Distance from the centre:

<a href="https://www.codecogs.com/eqnedit.php?latex=0\&space;m&space;-&space;R_{prc}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?0\&space;m&space;-&space;R_{prc}" title="0\ m - R_{prc}" /></a>.

In further calculations, only moun component of the EAS has to be taken into account. After carrying all steps of the analysis the final results are as follows:

//tableka

## Conclusions

As one can see, the results of the analysis and measurement differ significantly. It is most likely due to assumptions and simplifications that are still very far from reality, such as:
- All showers are treated as originating from protons while they make up only about 74% of primary cosmic-rays.
- Only muons in the showers are taken into account while it is most likely that electromagnetic component also gives signal.
- Fluctuations in number of produced particles and their distribution over distance was not yet fully characterised and sometimes an average was taken.
- There are many unclear effects that may have impact on the measurement like productions of particles in the upper parts of building in which it was performed.
- Background level might be under-evaluated.

However, as predicted the average number of coincidence signals caused by EAS is significantly higher than for the background. Thus, the level of confidence that certain event indicates occurrence of the cascade in a close surrounding of the system should be high. 

#### Advantages
This type of analysis has also another advantage, as it can even give some information about the energy of the primary particle that caused it. It is easy to see when the expected number of events for different number of coincidence signals is plotted over energy spectrum:

//obrazki

Thus, one can approximate the energy of the shower which is a useful information when searching for cosmic-rays ensembles.


# Plans for the future

Obviously, there are many way in which this work can be improved and developed. For now the main goal is to include particles other than muons, helium and probably iron nuclei as another primary particles, and to study the distribution of particles within the shower with more details.

Another plan is to build our own system of detectors. It would be just like the example presented in the results section. Than, acquired knowledge and possible future tests of behaviour of such system could help to take some more realistic assumptions in the analysis. Also with such tool working, new results of measurement could be constantly compared with computations which could help to improve this work. As written earlier, everything is still in progress.

Final goal of this mini project is to confirm that such method is a reliable way of detecting EAS.
