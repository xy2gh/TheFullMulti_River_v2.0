# TheFullMulti_River_v2.0
New version of The Full Multi River model after correction of several bugs found in the original model code [TheFullMulti_RIVER](https://github.com/Nano2PlastProject/TheFullMulti_RIVER)

The corrections include:

- Calculation of mixing rates between the surface, flowing and stagnant water compartments has been corrected to reflect the configuration described in Domercq et al. 2022 (Table S2). Mixing rates between stagnant water and sediment layer have been removed.

- Biofouling of already biofouled particles (i.e. "Biofouled-MP" and "heteroagg-Biofouled-MP") has been removed.

- Changes in the estimation of the fragmentation rates of biofouled particles have been introduced to only take the internal diameter of the biofouled particle (i.e. excluding the biofilm layer).

- Bug on the built up of the interactions dataframe have been corrected to properly assign sediment transport values as well as advective flow interactions



### Authors
===========

Repository and Jupyter Notebooks: Dr. Maria del Prado Domercq (@PradoDomercq)

Coding of python files and functions: Dr. Antonia Praetorius (@apraetorius) and Dr. Maria del Prado Domercq (@PradoDomercq)

Contact:

fullMultiModel@aces.su.se
