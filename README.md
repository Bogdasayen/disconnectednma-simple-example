# disconnectednma-simple-example
Simple example of aggregate level matching and reference prediction on network with 1 covariate
All code runs in the R statistical programming language linked to OpenBUGS using the R2OpenBUGS library. 
The script disconnected.nma.main.R calls data loading and model scripts to run independent baselines, ALM and reference prediction on an example evidence network. The evidence network is one draw from the simulation study with 5 connected RCTs with 100 patients per arm on treatments 1 and 2 and 3, 5 disconnected RCTs on treatments 4 and 5, and 10 single-arm studies with 5 on treatments 4 and 5 each. Data and code are described with in-line comments.  
