# DigitalPicoTools
DigiPicoTools

DIGITAL ANALYSIS OF PICOGRAM QUANTITIES OF DNA

Developed in the Ovarian Cancer Cell Laboratory, University of Oxford, United Kingdom.

This package offers a set of tools to process and analyse DNA Long Fragment Read sequencing data obtained from a DigiPico Experiment.

To Download and install the package.

A- Using devtools :

1- If you dont already have devtools installed, proceed to its installation :

install.packages("devtools")

library(devtools)

2- Download, install and launch OncoPhase :

install_github("OvarianCancerCell/DigiPicoTools", auth_token ="ce43XXXXXXXXXXXXXX")

As it is a private repository, an autnetification token will be requiered. Token by be obtained by contacting the Ovarian Cancer Lab.

library(DigiPicoTools)

help(DigiPicoTools)
