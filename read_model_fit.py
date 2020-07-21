#-*- coding: utf-8  -*-
"""
Author: michel-pierre.coll
Date: 2020-07-21 08:19:07
Description: Read .bin output from Huela's Irritability cpp scripts. Based
on "SEReadFlexible.m" from https://github.com/AndreasHula/Irritability.
"""

import numpy as np
import pandas as pd

# Participants in outputs
pairs = range(0, 128)
n_pairs = len(pairs)

# Initialise empty lists to capture results
inv_likelihood = []
inv_tom = []
inv_guilt = []
inv_aversion = []
inv_plan = []
inv_temp = []
inv_shift = []
inv_irritability = []
tru_likelihood=  []
tru_tom = []
tru_guilt = []
tru_aversion = []
tru_plan = []
tru_temp = []
tru_shift = []
tru_irritability = []

# Empty array of subjects x trials x possible actions
inv_expectations = np.empty(shape=(n_pairs, 10, 5))
tru_expectations = np.empty(shape=(n_pairs, 10, 5))

# Empty array of possible actions, guilt levels, subj, trial
inv_beliefs = np.empty(shape=(5, 3, n_pairs, 10))
tru_beliefs = np.empty(shape=(5, 3, n_pairs, 10))
# Empty array of possible actions, irritability levels, subj, trial
inv_tirrit = np.empty(shape=(5, 5, n_pairs, 10))
tru_tirrit= np.empty(shape=(5, 5, n_pairs, 10))

inv_shiftt = np.empty(shape=(6, 5, n_pairs, 10))
tru_shiftt = np.empty(shape=(6, 5, n_pairs, 10))

# Empty array of 5 actions? x 10 trials, x 3 guilt x 5 actions x 5 irr x subj
inv_act = np.empty(shape=(5, 10, 3, 5, 5, n_pairs))
tru_act = np.empty(shape=(5, 10, 3, 5, 5, n_pairs))

for sidx, s in enumerate(pairs):
    filename = 'model_fit/outputs/trust_params_pair' + str(s) + '.bin'
    with open(filename, mode='rb') as file: # b is important -> binary
        # Loop participants and read the main model parameters
            # Investor
            inv_likelihood.append(np.fromfile(file, np.double, count=1)[0])
            inv_tom.append(np.fromfile(file, np.int32, count=1)[0])
            inv_guilt.append(np.fromfile(file, np.int32, count=1)[0])
            inv_aversion.append(np.fromfile(file, np.int32, count=1)[0])
            inv_plan.append(np.fromfile(file, np.int32, count=1)[0])
            inv_temp.append(np.fromfile(file, np.double, count=1)[0])
            inv_shift.append(np.fromfile(file, np.int32, count=1)[0])
            inv_irritability.append(np.fromfile(file, np.int32, count=1)[0])
            # Trustee
            tru_likelihood.append(np.fromfile(file, np.double, count=1)[0])
            tru_tom.append(np.fromfile(file, np.int32, count=1)[0])
            tru_guilt.append(np.fromfile(file, np.int32, count=1)[0])
            tru_aversion.append(np.fromfile(file, np.int32, count=1)[0])
            tru_plan.append(np.fromfile(file, np.int32, count=1)[0])
            tru_temp.append(np.fromfile(file, np.double, count=1)[0])
            tru_shift.append(np.fromfile(file, np.int32, count=1)[0])
            tru_irritability.append(np.fromfile(file, np.int32, count=1)[0])

            # Trial wise values for investors
            for t in range(10):
                for k in range(5):
                    inv_expectations[sidx, t, k] = np.fromfile(file, np.double,
                                                            count=1)[0]
                    for g in range(3):
                        inv_beliefs[k, g, sidx, t] =  np.fromfile(file, np.double,
                                                            count=1)[0]
                    for irr in range(5):
                        inv_tirrit[k, irr, sidx, t] =  np.fromfile(file,
                                                                    np.double,
                                                                    count=1)[0]
                        inv_shiftt[k, irr, sidx, t]=  np.fromfile(file, np.double,
                                                            count=1)[0]
                        inv_shiftt[k+1, irr, sidx, t]=  np.fromfile(file, np.double,
                                                                count=1)[0]
                    for g in range(3):
                        for irr in range(5):
                            for act in range(5):
                                inv_act[act,t,g,k,irr,sidx] = np.fromfile(file,
                                                                    np.double,
                                                                    count=1)[0]
            # Trial wise values for trustee
            for t in range(10):
                for k in range(5):
                    tru_expectations[sidx, t, k] = np.fromfile(file, np.double,
                                                            count=1)[0]
                    for g in range(3):
                        tru_beliefs[k, g, sidx, t] =  np.fromfile(file, np.double,
                                                            count=1)[0]
                    for irr in range(5):
                        tru_tirrit[k, irr, sidx, t] =  np.fromfile(file,
                                                                    np.double,
                                                                    count=1)[0]
                        tru_shiftt[k, irr, sidx, t]=  np.fromfile(file, np.double,
                                                            count=1)[0]
                        tru_shiftt[k+1, irr, sidx, t]=  np.fromfile(file, np.double,
                                                                count=1)[0]
                    for g in range(3):
                        for irr in range(5):
                            for act in range(5):
                                tru_act[act,t,g,k,irr,sidx] = np.fromfile(file,
                                                                    np.double,
                                                                    count=1)[0]

# Put parameters in data frames

wide = pd.read_csv('data/trust_wide_data.csv')