#
# Copyright (c) 2025      High Performance Computing Center Stuttgart,
#                         University of Stuttgart.  All rights reserved.
#
# Authors: Kingshuk Haldar <kingshuk.haldar@hlrs.de>
#

import os
import subprocess
import re
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
import seaborn as sns
import seaborn.objects as so

print(sys.argv[1])
mon_fn= sys.argv[1]
#effs= [<lb>, <ser>, <trf>]
colnames= ['time', 'ideal_max', 'ideal_avg', 'useful_max', 'useful_avg', 'elapsed', 'ideal_loc', 'nevts_min']
df= pd.read_csv(mon_fn, sep=r'\s+', names= colnames, skiprows= 1)
df['time']= df.time* 1.0e-9

df['LB']= df.useful_avg/ df.useful_max
df['LB'].values[df['LB']> 1.0]= 1.0

df['Ser']= df.useful_max/ df.ideal_max
df['Ser'].values[df['Ser']> 1.0]= 1.0

df['Trf']= df.ideal_max/ df.elapsed
df['Trf'].values[df['Trf']> 1.0]= 1.0

row0= df.loc[0].values.flatten().tolist()
row0[0]= 0.0
df.loc[-1]= row0
df.index+= 1
df= df.sort_index()

ax= sns.lineplot(df, x= 'time', y= 'LB', drawstyle='steps-pre', color= '#f8766d', linewidth= 2, label= '$LB$', legend=True)
ax= sns.lineplot(df, x= 'time', y= 'Ser', drawstyle='steps-pre', color= '#00ba38', linewidth= 2, label= '$Ser$', legend=True)
ax= sns.lineplot(df, x= 'time', y= 'Trf', drawstyle='steps-pre', color= '#619cff', linewidth= 2, label= '$Trf$', legend=True)

ax.set_xlabel('Timeline [s]', fontsize=20)
ax.yaxis.label.set_visible(False)
ax.set_ylabel('Metrics values', fontsize=20)
x_max= plt.xlim()[1]
xshift= -15

#uncomment next 6 lines when the overall efficiencies are known
#ax.hlines(effs[0], xmin= 0, xmax=x_max, color='#f8766d', ls='dashed')
#ax.text(x_max+ xshift, effs[0]- 0.01,'$LB_{full}$',rotation=0)
#ax.hlines(effs[1], xmin= 0, xmax=x_max, color='#00ba38', ls='dashed')
#ax.text(x_max+ xshift, effs[1]- 0.01,'$Ser_{full}$',rotation=0)
#ax.hlines(effs[2], xmin= 0, xmax=x_max, color='#619cff', ls='dashed')
#ax.text(x_max+ xshift, effs[2]- 0.01,'$Trf_{full}$',rotation=0)

plt.ylim(-0.05,1.05)
#ax.set_aspect('0.5')
ax.grid()
plt.legend(loc= 0, fontsize=18)
plt.title('Time-resolved metrics', fontsize=20)
plt.savefig(os.path.splitext(mon_fn)[0]+ '.tiff', format= 'tiff', bbox_inches= 'tight')
plt.show()
