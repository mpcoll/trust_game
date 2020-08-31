#-*- coding: utf-8  -*-
"""
Author: michel-pierre.coll
Date: 2020-07-20
Description: Import data for the multiround trust game and exports to a
.bin file for model fitting.
# TODO Double check if rôle 1 == invest
# TODO Double check exceptions at line 70 and 74
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os
from os.path import join as opj
import shutil
import matplotlib.pyplot as plt
import seaborn as sns

# Change working dir for docker
os.chdir(os.path.dirname(os.path.realpath(__file__)))

###################################################################
# Import data in long format
###################################################################
# Read the data
wide_data = pd.read_csv('data/trust_raw_data.csv')

# Drop participants with missing values
wide_data = wide_data[wide_data['code'] != 999]

# Drop participants with type de pair == 1
wide_data = wide_data[wide_data['Type_de_paire'] != 1]


# Loop groups and find pairing so we get a number for each pair across all groups
# Probably a better way to do it but here find people who played exactly de same game
game = []
for row in wide_data.iterrows():
    strgame = ''
    for tr in range(1, 11):
        strgame = strgame + str(row[1]['inv_' + str(tr)]) + str(row[1]['remise_' + str(tr)])

    game.append(strgame + str(row[1]['Groupe']))
wide_data['game'] = game

# Drop participants with no game data
wide_data = wide_data[wide_data['game'] != str(999)*20]

# Assign a code to each pair
pairs = [999]*len(list(wide_data['game']))
isalone = [999]*len(list(wide_data['game']))
for idx, ga in enumerate(list(set(wide_data['game']))):

    if len(np.where(np.asarray(wide_data['game']) == ga)[0]) > 2:
           print('More than two identical games found for these IDs')
           print(np.asarray(wide_data['ID'])[np.where(np.asarray(wide_data['game']) == ga)[0]])

    pairs = np.where(np.asarray(wide_data['game']) == ga, idx, pairs)
    # Find participants not part of a pair
    if len(np.where(np.asarray(wide_data['game']) == ga)[0]) < 2:
        isalone = np.where(np.asarray(wide_data['game']) == ga, 1, isalone)
    else:
        isalone = np.where(np.asarray(wide_data['game']) == ga, 0, isalone)

wide_data['isalone'] = isalone
wide_data['pairs'] = pairs


# _________________________________________________________________
# Keep only one line/pair

# For this pair, one role == 999, assumed it is the missing one.
wide_data['rôle'][wide_data['ID'] == 110] = 2

# For this pair, both role were == 2, assumed it the first one == 1
wide_data['rôle'][wide_data['ID'] == 125] = 1

# Drop participants that are alone
wide_data = wide_data[wide_data['isalone'] == 0]

wide_data_pairs = pd.DataFrame()
for p in list(set(pairs)):

    inv_row = wide_data[(wide_data['pairs'] == p) & (wide_data['rôle'] == 1)]
    inv_row.columns = ['i_' + c for c in list(inv_row.columns)]

    tru_row = wide_data[(wide_data['pairs'] == p) & (wide_data['rôle'] == 2)]
    tru_row.columns = ['t_' + c for c in list(tru_row.columns)]

    pairdata = pd.concat([inv_row.reset_index(), tru_row.reset_index()],
                         axis=1)
    wide_data_pairs = pd.concat([wide_data_pairs, pairdata.reset_index()],
                                axis=0)

wide_data_pairs['pairs_valid'] = range(0, len(wide_data_pairs))

wide_data_pairs.to_csv('data/trust_wide_data.csv')


# _________________________________________________________________
# Transform data to long format

# Init list to collect data frames
pairs_df = []

# Loop participants
for idx, p in enumerate(list(set(wide_data_pairs['pairs_valid']))):

    # keep data for this pair only
    pairdat = wide_data_pairs[wide_data_pairs['pairs_valid'] == p]

    # Create a new long dataset for this part with one row per trial
    pairdat_new = pd.DataFrame(data={'pair': [int(pairdat.pairs_valid)]*10,
                                     'trials': range(1, 11),
                                     'bin_id': idx})

    # Add all the other data to each row
    for c in list(wide_data_pairs.columns):
        pairdat_new[c] = [list(pairdat[c])[0]]*10

    # Unpack the investment and return at each trial
    # We could add the other values but we just need investment and return,
    # we can calulate the other stuff on the fly
    invests, returns = [], []
    for tr in range(1, 11):
        inv = int(pairdat['i_inv_' + str(tr)])
        ret = int(pairdat['i_remise_' + str(tr)])
        invests.append(inv)
        returns.append(ret)

    pairdat_new['invest'] = invests
    pairdat_new['returns'] = returns

    # Add to list of data frames
    pairs_df.append(pairdat_new)

# Concatenate all subs in a single long data frame of shape ntrials * nparticipants
longdata = pd.concat(pairs_df)

print(str(len(set(longdata.i_pairs))) + " valid pairs")


###################################################################
# Plot average behaviour
###################################################################

longdata['perc_inv'] = longdata['invest'] / 20 * 100
longdata['perc_ret'] = np.where(longdata['invest'] == 0, 0,
                                longdata['returns']
                                / (longdata['invest']*3)*100)
longdata['prop_inv'] = longdata['invest'] /20
longdata['prop_ret'] = np.where(longdata['invest'] == 0, 0,
                                longdata['returns'] / (longdata['invest']*3))


# Plot percent invested and percent returned across trials
# (ie investor and trustee, respectively(?))
fig, ax = plt.subplots(facecolor='white')
sns.pointplot(x='trials', y='perc_inv', data=longdata,
              ci=68, label='investment', ax=ax) #NOTE CI = 68 shows sem
sns.pointplot(x='trials', y='perc_ret', data=longdata,
              ci=68, color='g', label='returns', ax=ax)
ax.set_ylabel('% of amount received', fontsize=25)
ax.set_xlabel('Trials', fontsize=25)
ax.tick_params('both', labelsize=20)
ax.legend(handles=[ax.lines[::11][0], ax.lines[::-1][0]],
          labels=["Investor", "Trustee"], loc='lower left',
          fontsize=20, frameon=False)

# Plot...actual amounts?
fig, ax = plt.subplots(facecolor='white')
sns.pointplot(x='trials', y='invest', data=longdata,
              ci=68, label='investment', ax=ax)
sns.pointplot(x='trials', y='returns', data=longdata,
              ci=68, color='g', label='returns', ax=ax)
ax.set_ylabel('Amount', fontsize=25)
ax.set_xlabel('Trials', fontsize=25)
ax.tick_params('both', labelsize=20)


###################################################################
# Format data for iPOMDP
###################################################################
# Single vector with 20 rows per game, in the format inv ret inv ret ...

alldata_vector = []
invest_grid, return_grid = [], []

# Get vectors of values
investments = list(longdata['prop_inv'])
returns = list(longdata['prop_ret'])

# Grid to round amounts to
inv_grid = np.array([0, 0.25, 0.5, 0.75, 1])
ret_grid = np.array([0, 1/6, 1/3, 1/2, 2/3])

for i in range(len(longdata)):
    # get values
    inv = investments[i]
    ret = returns[i]
    # Find closest value in grid
    inv_ingrid = np.argmin(np.abs(inv_grid - inv))
    ret_ingrid = np.argmin(np.abs(ret_grid - ret))
    invest_grid.append(inv_ingrid)
    return_grid.append(ret_ingrid)

    # Append to vector
    alldata_vector.append(inv_ingrid)
    alldata_vector.append(ret_ingrid)


longdata['inv_grid'] = invest_grid
longdata['ret_grid'] = return_grid

# Save longdata to csv
longdata.to_csv('data/trust_long_data.csv')

# Save to binary for the C++ script
np.asarray(alldata_vector).astype('int32').tofile('data/trust_data.bin')


