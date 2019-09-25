import pandas as pd
import seaborn as sns
import numpy as np
from scipy.stats import wilcoxon
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt

## Load and reshape data
data = pd.read_csv('atp_ip10.csv', header=0)

### correct baseline by subtracting the first measurement of every control to the corresponding control and IP10
### Then calulate the log difference in both groups.

controls = data[['Control 1', 'Control 2', 'Control 3', 'Control 4']]
peptides = data[['IP10 1', 'IP10 2', 'IP10 3', 'IP10 4']]

base = controls.loc[0].values

controls = np.log(controls) - np.log(base)
peptides = np.log(peptides) - np.log(base)

### assemble data frame and reshape it top a tidy format

data = pd.concat([data['Time(s)'], controls, peptides], axis=1)
long_data = pd.melt(data, id_vars='Time(s)', value_vars=['Control 1', 
                                                        'IP10 1', 'Control 2', 'IP10 2',
                                                        'Control 3', 'IP10 3', 'Control 4',
                                                        'IP10 4'])
## Create new variable from strings in column and concatenate with long_data DF

long_data = pd.concat([long_data,
                       long_data.variable.str.split(' ', expand=True).rename({0:'treatment', 1:'repetition'}, axis=1)],
                       axis=1)
long_data.drop('variable', axis=1, inplace=True)

## Create boxplot for comparing distributions

f, ax = plt.subplots(figsize=(7, 6))                                                                                                          

sns.boxplot(x='treatment', y='value', data=long_data, palette='GnBu_d')                                                                       
#sns.swarmplot(x='treatment', y='value', data=long_data, color='black', alpha=0.3, jitter=True)                                                                     
ax.set(ylabel='Total relative change\n Fold decrease')
ax.yaxis.label.set_size(20)
ax.set(xlabel='Treatment')
ax.xaxis.label.set_size(20)
plt.savefig('boxplot.svg')                                                                                                        


## Create lineplots
plt.figure(figsize=(4, 3))
sns.lineplot(x='Time(s)', y='value', hue='treatment', style='treatment', data=long_data)
plt.savefig('mean_line.svg')
plt.figure(figsize=(4, 3))
sns.lineplot(x='Time(s)', y='value', hue='treatment', style='repetition', data=long_data)
plt.savefig('all_line.svg')

## statistical analysis
# Mann Whiteney U

man_s, man_p = mannwhitneyu(long_data.value[long_data['treatment'] == 'Control'], long_data.value[long_data['treatment'] == 'IP10'])
w_s, w_p = wilcoxon(long_data.value[long_data['treatment'] == 'Control'], long_data.value[long_data['treatment'] == 'IP10'])

alpha = 0.05

if w_p < alpha:
    print(f'The differences are statistically significant with a p value of {w_p}, we reject H0')
else:
    print('Both groups come from a population with the same distribution, we accept H0') 

print('Done master!! :)')
