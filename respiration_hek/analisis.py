"""Script to analyze data from Oroboros Oximeter."""

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm 
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd                                             
from statsmodels.stats.multicomp import MultiComparison 
from statsmodels.stats.libqsturng import psturng


#load csv data
#respiration
oroboros = pd.read_csv("datos.csv")

#number of cells

cell_number = pd.read_csv("Cell_number.csv")

#select and slice the data for use
#Starting points

index_list = oroboros[oroboros["Event"].notnull()].index[:] 

#lists for new data frame
event = []
slopeA = []
rsquaredA = []
slopeB = []
rsquaredB = []
experiment = []
condition = []

#Chamber A data
#Calculates the slope for each condition from minute 1 to minute 4 after treatment

for i in index_list:
	
	if oroboros.loc[i]["Chamber"] == "Both" or oroboros.loc[i]["Chamber"] == "Right":
		event.append(oroboros.loc[i]["Event"])
		experiment.append(oroboros.loc[i]["Experiment"])
		condition.append(oroboros.loc[i]["Condition"])
		x = oroboros.loc[i+20:i+80]["Minute"]
		y = oroboros.loc[i+20: i+80]["A"]
		slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
		slopeA.append(slope)
		rsquaredA.append(r_value ** 2)
	
#Chamber B data

for i in index_list:
	
	if oroboros.loc[i]["Chamber"] == "Both" or oroboros.loc[i]["Chamber"] == "Left":
		x = oroboros.loc[i+20:i+80]["Minute"]
		y = oroboros.loc[i+20 : i+80]["B"]
		slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
		slopeB.append(slope)
		rsquaredB.append(r_value** 2)

#Add cell number to df

#create df from previous lists 

tables = [experiment, condition, event, slopeA, slopeB, rsquaredA, rsquaredB]
colnames = ["Experiment","Condition", "Event", "SlopeA", "SlopeB", "RsquaredA", "RsquaredB"]
new_df = pd.DataFrame(tables)
new_df = new_df.transpose()
new_df.columns = colnames

#merge slope data frame with cell number df

df = pd.merge(new_df, cell_number)

#normalize slopes per million cells

a = (df.SlopeA/df.CellsA)*1000000
b = (df.SlopeB/df.CellsB)*1000000

#add normalized values columns to df

df["NormA"] = a
df["NormB"] = b

#Extract data frame and clean it
df1 = df.iloc[:,0:3]
df2 = df.iloc[:, -2:]
df_1 = pd.concat([df1, df2], axis = 1)
df_2 = pd.melt(df_1, id_vars=["Experiment", "Condition", "Event"], value_vars= ["NormA", "NormB"], var_name="Chamber", value_name= "Slope")


#Change Names of events
df_2.loc[df_2.Event.str.startswith("HEK293_IP10"), "Event"] = "Pre-incubation 10uM"
df_2.loc[df_2.Event.str.startswith("HEK293_IP50"), "Event"] = "Pre-incubation 50uM"
df_2.loc[df_2.Event.str.startswith("HEK"), "Event"] = "Cells"
df_2.loc[df_2.Event.str.startswith("DM"), "Event"] = "Medium"

#Create a final array with selected events only

final_df = df_2.loc[(df_2.Event != 'ANTI')&(df_2.Event != 'END')&(df_2.Event != 'Medium')&(df_2.Event !='OLIGO')]

treatments=list(np.unique(final_df.Event))
treatments.remove('Cells')


#plot
#Create multipanel figure with 3 rows and 2 columns
fig, ax = plt.subplots(3,2)

#iterate over treatments to compare each one vs Control (cells)

for ind, panel in enumerate(ax.flatten()):
	sns.stripplot(x="Event", y="Slope", jitter=True, data=final_df[(final_df.Event=='Cells')|(final_df.Event==treatments[ind])], ax=panel)
	
plt.savefig('respiration.svg')
plt.show()

#statistical analysis

final_df['Slope'] = final_df.Slope.astype(np.float)

#fit linear model
model = ols('Slope ~ Event + Experiment', data=final_df).fit()
print(model.summary())

print('ANOVA analysis')

aov_table= sm.stats.anova_lm(model, type=2)
print(aov_table)

print('Post hoc tukey')

mc = MultiComparison(final_df['Slope'], final_df['Event'])
mc_results = mc.tukeyhsd()
print(mc_results.summary())

p_values = psturng(np.abs(mc_results.meandiffs / mc_results.std_pairs), len(mc_results.groupsunique), mc_results.df_total)

print('p_values: ', p_values)