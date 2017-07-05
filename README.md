# Table of Contents
1. [Demo Summary](README.md#demo-summary)
2. [Step-by-Step](README.md#step-by-step)
3. [Dependencies](README.md#dependencies)

# Demo Summary
In this snippet I show how to integrate a two-compoenent sersic profile, given the coefficients of best-fit. 
Using those results, I convert the calculated total magnitude within a specified radius to a stellar mass using the stellar sysnthesis population models generated with [EZGAL](https://www.baryons.com/ezgal).
Finally, I pull this information together with two other samples of galaxy clusters with stellar mass measurements and look at the stellar mass as a fuction of total cluster mass, M_500.
I fit the relation ![eqn](https://github.com/wickedalchemist/StellarMass_Sersics_Demo/blob/master/CodeCogsEqn.png) with an orthogonal distance regression (ODR), which takes into account uncertainties in both axes.


# Step by Step

1. Decide on the annulus within with you want to know the stellar mass of the brightest cluster galaxy and extended halo of intracluster light (BCG+ICL). These are set in m500_stellarmass_ring.py with the variables term1 and term2 (term2 = outer radius, term1 = inner radius in units of kpc)
2. On command line run `/astro/data/siesta1/soft/anaconda_2.2.0/bin/python m500_stellarmass_ring.py`  
  * This program first uses the function AG2005_sersics_2dev() defined in integrate_sersics.py to find the total magnitude (brightness) of the BCG+ICL within the specified radiual range for a sample of low-redshift (nearby) galaxy clusters from the paper by Gonzalez et al. (2005). 
  * Then using the stellar population synthesis models provided in the Python module EZGAL (see above for link to package), this total magnitude is converted to a total stellar mass.
  * I have two additional samples of galaxy clusters for which I have previosly measured the stellar mass based on Hubble Space Telescope imaging data that I have reduced and analyzed (i.e. not simply used published profiles to derive these quantities). I have stored the stellar mass and total cluster mass for each galaxy cluster in these additional samples in a Python pickle file.
  * The relation of interest between the stellar mass of the BCG+ICL and total cluster mass that hosts the BCG+ICL, is thought to go as ![eqn](https://github.com/wickedalchemist/StellarMass_Sersics_Demo/blob/master/CodeCogsEqn.png). I fit this function using my mid-redshift sample (from the pickle file) using an orthogonal distance regression, which takes into consideration uncertaintines the values of both stellar mass and cluster mass.
  * Using the resulting best-fit slope, alpha, from the above fit, I then find the best-fit relation with the same functional form for the high-redshift sample (also from the pickle).
3. Plot! The stellar mass as a function of cluster mass is plotted for all systems in the composite sample. The best fit functions found above are overlaid. This plot is displayed as well as saved as a png file, with a name that changes to reflect the choice of radial range probed from step 1 above. An example of the plot produced: ![stellar-mass-cluster-mass](https://github.com/wickedalchemist/StellarMass_Sersics_Demo/blob/master/sm100-10_m500.png)

# Dependencies
To run this code you will need the Python package [EZGAL](https://www.baryons.com/ezgal), a tool that takes stellar population synthesis models and evolves them with time to calculate magnitudes, mass-to-light ratios, passband corrections, etc in your desired filter. 
Other Python packages used include: glob, pickle, scipy, mpmath, and matplotlib.pypot.
