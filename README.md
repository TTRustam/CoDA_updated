# CoDA updated data, functions and results.
This file contains some new data and rewritten and restructured scripts. This was done to ease the reading and testing of the obtained results. Code contains comments to explain what is going on there. 
I have deliberately excluded parts of the code that aim for data manipulation of LT calculation only. I have these scripts locally and will add them to Git as soon as I will be sure that the analysis is correct. I have thus created several data files to simply read them into R, which in my idea will simplify the review process. 

# Contents
**exposure** folder contains the mid-year population estimates for Russia (Rus) /n
mortality contains the total Mx, cause specific Mx, and cause-specific cancer mortality.
**life_table** is Rus life table (Calculated following method described in Preston et al.)
**RMSE** folder contains two data files with observed dx to compare with the fitted values.
**coda_functions** folder contains several files: 
**auxillary_functions** are the functions called by other scripts. They mostly do prefill other functions or are used for data cleaning and manipulation. 
**inner_functions** are the functions used to calculate the SVD and resulting LT after forecasting by the main functions.
**main_functions** are the 4 compositional models.
**to_plot** (a rather ugly name, I know and shall rename this file slightly later) contains some functions used for plotting and data manipulation process that were pulled from the main analysis to ease the review process. **It also contains the RMSE calculation function** which is very important for us, since this is how we select the fitting period and best model.
**rmse_results** - is a markdown file containing some short description and tables for RMSE checks on different fitting period and also shows the best one.
**results** Is the main file that calls all other files to produce the results. It prints the most important figures and tables and has some comments (Note that what is being printed starts in line 220)
**soren_ppaper** - contains the Sorens paper itself, as well as its statistical supplement that went onluie.
**old_version_of_the_paper** Contains the old version of the draf prepared during EDSD.

# Short description of results, ideas and timeline
I made a mistake in a code (not in a CoDA models, but in further calculations), so some previous conclusions were partially incorrect. It seems, that decreasing the fitting period only works if we decrease it by 10 years. If more, the RMSE actually increases, and if decrease fitting period to < 1992, models go mental and VECM absolutely fails. This is an interesting finding since LC model do not seem to have such a heavy restriction. Conclusion here is that the CoDA models should be used only if the long enough time series is available.
Even with taking the 1975:2018 perspective, males still have high error values, but at this point I officialy have only one guess on how to improve this. I tried different arima models, with negligible effect, width and number of age groups also do not rely affect the model performance, any type of LT related problem also seems to be invalid in the end. MAYBE if we change the number of causes and their grouping by assigning a meaningful (with uniform dynamics) groups to other cancers and even to those causes that are not being considered explicitly, we can improve the model performance. This is the only thing I can think about and it has no grantee to work.
One additional possibility is that models mainly use rank 3 approximation for SVD, while 2s-CoDA uses rank 1 approximation, so maybe changing it to rank 3 will yield better results, but to do so, I will have to rewrite the code, as I should confess, I did not add this option to the current code version.

Other than that I think we pulled everything we could from the models. In terms of text itself, I think that maybe the part of the paper that ascribes the dynamics of age-standardized cancer Mx should be removed as it (in my opinion) do not add any particular value to the research. Maybe it will be better to substitute it with age specific Mx dynamics. I also think that adding at least one country for comparison might be beneficial. 
 
# Difference in RMSE calculations
Soren calculates the RMSE by each age group and then simply takes the average. I combine all the age groups and then calculate RMSE. We both use rolling bases for the calculation. He also calculates RMSE for LC model in terms om Mx, while I use dx for this task. I do not think that any of these affects results significantly, but is still worth mentioning.


# Conclusions
I do not see any methodological contribution here, so I would frame the conclusion in terms of the Russian case vs. the chosen comparator + some ideas on why the models perform poorly in males. This reasons are:
Models depend on the long enough and stable time series of mortality.
The unstable k(t) - the time pattern of mortality resulting from significant fluctuations in Mx.
POSSIBLY: Choice of grops of causes of death considered in model affect its result.

# PS
I would highly appreciate any comments on code and calculations. If there is no mistakes, I can easily finalize the paper in 2-3 days, and we can proceed to the application process. However if there is something else we could pull off the models or there is any mistake, I would really like to find it beforehand.






