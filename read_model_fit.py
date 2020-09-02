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

for codeversion in ['model_fit_v1', 'model_fit_v2']:

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
        filename = codeversion + '/outputs/trust_params_pair' + str(s) + '.bin'
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

    # Add columns
    wide['inv_likelihood'] = inv_likelihood
    wide['inv_tom'] = inv_tom
    wide['inv_guilt'] = inv_guilt
    wide['inv_aversion'] = inv_aversion
    wide['inv_plan'] = inv_plan
    wide['inv_temp'] = inv_temp
    wide['inv_shift'] = inv_shift
    wide['inv_irritability'] = inv_irritability

    wide['tru_likelihood'] = tru_likelihood
    wide['tru_tom'] = tru_tom
    wide['tru_guilt'] = tru_guilt
    wide['tru_aversion'] = tru_aversion
    wide['tru_plan'] = tru_plan
    wide['tru_temp'] = tru_temp
    wide['tru_shift'] = tru_shift
    wide['tru_irritability'] = tru_irritability

    # Save
    wide.to_csv('data/trust_wide_data' + codeversion + '.csv')

    long = pd.read_csv('data/trust_long_data.csv')

    long['inv_likelihood'] = np.repeat(inv_likelihood, 10)
    long['inv_tom'] = np.repeat(inv_tom, 10)
    long['inv_guilt'] = np.repeat(inv_guilt, 10)
    long['inv_aversion'] = np.repeat(inv_aversion, 10)
    long['inv_plan'] = np.repeat(inv_plan, 10)
    long['inv_temp'] = np.repeat(inv_temp, 10)
    long['inv_shift'] = np.repeat(inv_shift, 10)
    long['inv_irritability'] = np.repeat(inv_irritability, 10)

    long['tru_likelihood'] = np.repeat(tru_likelihood, 10)
    long['tru_tom'] = np.repeat(tru_tom, 10)
    long['tru_guilt'] = np.repeat(tru_guilt, 10)
    long['tru_aversion'] = np.repeat(tru_aversion, 10)
    long['tru_plan'] = np.repeat(tru_plan, 10)
    long['tru_temp'] = np.repeat(tru_temp, 10)
    long['tru_shift'] = np.repeat(tru_shift, 10)
    long['tru_irritability'] = np.repeat(tru_irritability, 10)

    long.to_csv('data/trust_long_data'  + codeversion +  '.csv')


v1 = pd.read_csv('data/trust_wide_data' + 'model_fit_v1' + '.csv')
v2 = pd.read_csv('data/trust_wide_data' + 'model_fit_v2' + '.csv')

import matplotlib.pyplot as plt
for toplot in ['inv_likelihood', 'inv_tom',
               'inv_guilt', 'inv_aversion', 'inv_plan', 'inv_temp',
               'inv_shift', 'inv_irritability',
               'tru_likelihood', 'tru_tom', 'tru_guilt', 'tru_aversion',
               'tru_plan', 'tru_temp', "tru_shift", 'tru_irritability']:
    plt.figure(figsize=(10, 5))
    plt.plot(v1['pairs_valid'], v1[toplot], label='old_code')
    plt.plot(v1['pairs_valid'], v2[toplot], label='new_code')
    plt.xlabel('Pairs')
    plt.ylabel(toplot)
    plt.title(toplot)
    plt.legend()

    plt.figure(figsize=(10, 5))
    plt.plot(v1['pairs_valid'], v1[toplot]-v2[toplot], label='Difference')
    plt.xlabel('Pairs')
    plt.ylabel(toplot)
    plt.title(toplot + ' Difference between old and new code')
    plt.legend()