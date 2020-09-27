# CoDA updated data, functions and results.
This file contains some new data and rewritten and restructured R scripts. Restructurization was done to ease the review process and testing of the obtained results. Code contains comments to explain what is going on there.\
I have deliberately excluded parts of the code that aim for data manipulation or LT calculations. I have these scripts locally and will add them to Git as soon as I will be sure that the analysis is correct. I have thus created several data files to simply read them into R, which in my idea will simplify the review process.\

# Contents
**exposure** folder contains the mid-year population estimates for Russia (Rus)\
**mortality** contains the total Mx, cause specific Mx, and cause-specific cancer mortality\
**life_table** is Rus life table (Calculated following method described in Preston et al.)\
**RMSE** folder contains two data files with observed dx to compare with the fitted values\
**coda_functions** folder contains several files:\
**auxillary_functions** are the functions called by other scripts. They mostly do prefill other functions or are used for data cleaning and manipulation\
**inner_functions** are the functions used to calculate the SVD and resulting LT after forecasting by the main functions\
**main_functions** are the 4 compositional models\
**to_plot** (a rather ugly name, I know and shall rename this file slightly later) contains some functions used for plotting and data manipulation process that were pulled from the main analysis to ease the review process. **It also contains the RMSE calculation function** which is very important for us, since this is how we select the fitting period and best model\
**rmse_results** - is a markdown file containing some short description and tables for RMSE checks on different fitting period and also shows the best one\
**results** Is the main file that calls all other files to produce the results. It prints the most important figures and tables and has some comments (Note that what is being printed starts in line 220)\
**soren_ppaper** - contains the Sorens paper itself, as well as its statistical supplement that went onluie.\
**old_version_of_the_paper** Contains the old version of the draf prepared during EDSD.\

# Short description of results, ideas and timeline
I made a mistake in a code (not in a CoDA models, but in further calculations), so some previous conclusions were partially incorrect. It seems, that decreasing of the fitting period only works if it does not exceed 10-20 years. If more, the RMSE actually increases, and if decrease fitting period to 1992+, models go mental and VECM-model absolutely fails and sometimes even not working. This is an interesting finding since LC model do not seem to have such a heavy restriction. Conclusion here is that the CoDA models should be used only if the long enough time series of mortality is available.
Even with taking the best possible fitting period - 1975:2018 perspective, males still have high error values, but at this point I officialy have only one guess on how to this can be approached. I tried different arima models with negligible effect, width and number of age groups also do not rely affect the model performance, any type of LT related problem also seems to be of no importance in the end. MAYBE if we change the number of causes and their grouping by assigning a meaningful (with uniform mortality dynamics) groups to other cancers and even to those causes that are not being considered explicitly (CVD, External), rather than collapsing them in two quantities, we can improve the model performance. This is the only thing I can think about and it has no grantee to work. Does it even worth trying this, or it should have no effect due to any theory that I miss?\
One additional possibility is that most CoDA models use rank 3 approximation for SVD, while 2s-CoDA uses rank 1 approximation, so maybe changing it to rank 3 will yield better results for this particular model, but to do so, I will have to rewrite the code substantially, as I should confess, I did not add this option to the current code version.\

Other than that I think we pulled everything we could from the models. In terms of the paper itself, I think that maybe the part of the paper that ascribes the dynamics of age-standardized cancer Mx should be removed as it (in my opinion) do not add any particular value to the research. Maybe it will be better to substitute it with age specific Mx dynamics or just cut this part, as we already have some compositional dx distributions. I also think that adding at least one comparator country might be beneficial in terms of general robustness of interpretaion. SOmething like England and Wales might work, as it has close to average European epidemiologic chcracteristics of cancer.\
 
# Difference in RMSE calculations
There is some difference in ,ethodlogy of RMSE calculation that I use compared to that of Soren et al. They calculates the RMSE for each age group and then simply take the average. I combine all the age groups and then calculate RMSE. We both use rolling bases for the calculation. He also calculates RMSE for LC model in terms om Mx, while I use dx for this task. I do not think that any of these affects results significantly, but is still worth mentioning.\
NOTE: RMSE listed in old version of paper are essentially wrong, so best to reffer to _rmse_results_ file.\

# Conclusions
I do not see any methodological contribution made by me, other than the observational ones. I would frame the conclusion in terms of the Russian case vs. the chosen comparator + some ideas on why the models perform poorly in males. This reasons are:\
Models depend on the long enough time series of mortality with preferably uniform k(t) (the time pattern of mortality) dynamics and no significant fluctuations + POSSIBLY certainly groupped causes of death.\
LC at the same time, has less restrictions, so one of the conclusions should be that LC can be preffered in case of very unstable dynamics or short time series. Although in both cases it might be the best possible option, but still a bad one.\
Also, it seems that forecasting 20 years ahead given the observed dynamics and its uncertainty is a false approach for Russia, and we should focus of 10 years ahead forecast.\
But again, all these conclusions make sence if I made no mistakes in a code. I cheked it many times and found no other mistakes, but I think that review would be great still. I would highly appreciate any comments on code and calculations. If there is no mistakes, I can easily finalize the paper in 2-3 days, and we can proceed to the journal application process itsef. However if there is something else we could pull off the models or there is any mistake, Iam very much eager to fix it.\
