'''


Automatic Detection of Inverse Stable Isotopic Labeling
Input: .csv alignment from MZmine 3 and at least 3 (12C, 13C, 13C+12C-precursor) raw LC-MS data in .h5 format.
Output: .csv of features showing untargeted inverse labeling. 
Authors: Tashi Liebergesell (tashi.liebergesell@utah.edu) & Ethan Murdock (Ethan.Murdock@hsc.utah.edu)
Thanks to Benjamin Bowen and Thomas Harwood (Lawrence Berkeley National Laboratory) for help with the Metabolite Atlas Python package.


'''

#%% imports
import time
start_time = time.time()
import pandas as pd
import numpy as np
import os
import copy

os.chdir("c:/users/mudoe/desktop/hetexp working/metatlas-main")  #change file path as needed
from metatlas.io import feature_tools as ft

import math
from scipy.signal import argrelextrema
import glob

# progress bar
from tqdm import tqdm

#%%files list
os.chdir("c:/users/mudoe/desktop/InverSIL/355") #Change path as needed
samples = glob.glob(os.path.join('AWP355_1213.csv'))  #Rename as needed. 
samples = pd.DataFrame(samples)
samples[0] = samples[0].apply(lambda x: os.path.basename(x).split('_')[0])

#%%calculate unique 12C and 13C features found in each condition

df = pd.read_csv("./AWP355_1213.csv")    
df = df.dropna(axis = 1)
df.columns = pd.Series(['RT','MZ_12C','13'])
df12 = df[df['MZ_12C'] == (df['MZ_12C'] + df['13'])] #unique 12C features
df13 = df[df['13'] == (df['MZ_12C'] + df['13'])] #unique 13C features
print(str(len(df12)) + " unique 12C features")
print(str(len(df13)) + " unique 13C features")

#%%combos
combos = {}
for n in range(len(samples)):
    df = pd.read_csv(samples.iat[n,0] + "_1213.csv")
    df = df.dropna(axis = 1)
    df.columns = pd.Series(['RT','MZ_12C','13'])
    df = df.replace(0, np.nan)
    align_12C_13C = df[df['MZ_12C'].notnull() & df['13'].notnull()].index 
    df.drop(align_12C_13C, inplace = True) 
    tw = df[df["MZ_12C"].notnull()].loc[:,["RT","MZ_12C"]]
    th = df[df["13"].notnull()].loc[:,["RT","13"]]
    mz_tolerance = 5    # Define m/z tolerance [m/z]. Change if needed
    rt_tolerance = 0.1  # Define RT tolerance [min]. 
    tw["_temp"] = 0
    th["_temp"] = 0
    combined = th.merge(tw, 'outer', '_temp', suffixes=['_13C', '_12C'])
    
    combined = combined[
        ((combined['13']  - combined['MZ_12C']) > mz_tolerance) &
        ((combined['13']  - combined['MZ_12C']) <= 51) &
        ((combined['RT_13C'] - combined['RT_12C']).abs() < rt_tolerance)]
    combined = combined.drop(['_temp'], axis=1)
    combined.sort_values(by=['13'], inplace=True)
    combstore = copy.deepcopy(combined)

    
    trt = combined["RT_12C"].unique() 
    print('pair filtering 1:')
    for i in tqdm(range(len(combined))):
        for j in range(i+1,len(combined)):
            if combined.iat[i,0] != 0:
                if combined.iat[i,3]-(combined.iat[i,3]/200000) < combined.iat[j,3] < combined.iat[i,3]+(combined.iat[i,3]/200000):
                    if abs(combined.iat[i,2] - combined.iat[j,2]) < .05:
                            combined.iat[i,0] = 0
    combined = combined.replace(0, np.nan)
    combined = combined.dropna(axis = 0)
    combined.sort_values(by=['MZ_12C'], inplace=True,ascending=False)
    
    print('pair filtering 2: ')
    for i in tqdm(range(len(combined))):
        for j in range(i+1,len(combined)):
            if combined.iat[i,0] != 0:
                if combined.iat[i,1]-(combined.iat[i,1]/200000) < combined.iat[j,1] < combined.iat[i,1]+(combined.iat[i,1]/200000):
                    if abs(combined.iat[i,0] - combined.iat[j,0]) < .05:
                            combined.iat[i,0] = 0
    combined = combined.replace(0, np.nan)
    combined = combined.dropna(axis = 0)
    combined.sort_values(by=['13'], inplace=True)
    combos[samples.iat[n,0]] = combined

print(str(len(combined)) + "Likely 12C/13C pairs")    
    

#%%load data
# will need to convert data to .h5 format

os.chdir("c:/users/mudoe/desktop/InverSIL/355") #Change path as needed
data = {}
data["atlas"] = 0
data["lcmsrun"] = "AWP273_P1.h5" #placeholder for the file location. Rename if needed. Make sure 13C + 12C-precursor samples are named, P1, P2, P3, etc. 
data["file_index"] = 1
data["polarity"] = "positive"  #can switch polarity if needed
ppm = 15 #change tolderance if needed

data12 = {}
data12["atlas"] = 0
data12["lcmsrun"] = "AWP273_12C.h5" #placeholder for the file location. Rename if needed. 
data12["file_index"] = 2
data12["polarity"] = "positive"  #can switch polarity if needed
ppm2 = 15 #change tolderance if needed

data13 = {}
data13["atlas"] = 0
data13["lcmsrun"] = "AWP273_13C.h5" #placeholder for the file location. Rename if needed. 
data13["file_index"] = 3
data13["polarity"] = "positive"   #can switch polarity if needed
ppm2 = 15 #change tolderance if needed


#%%loop to check pairs 
collection = {}
for n in range(len(samples)):
    key = list(combos.keys())[n]
    combined = copy.deepcopy(combos.get(key))
    data12['lcmsrun'] = key + '_12C.h5'
    data13['lcmsrun'] = key + '_13C.h5'
    twelvers = []
   
    for i in range(len(combined)):
        twelvers.append([combined.iat[i,3],combined.iat[i,2],combined.iat[i,2]+.1,combined.iat[i,2]-.1,
                    ppm2,i+1,0,i+2,i+3]) #load in 12
        twelvers.append([combined.iat[i,1],combined.iat[i,0],combined.iat[i,0]+.1,combined.iat[i,0]-.1,
                    ppm2,len(combined)+i+2,0,len(combined)+i+3,len(combined)+i+4]) #load in 13
    twelvers = pd.DataFrame(twelvers)
    twelvers = twelvers.rename(columns={0:"mz",1:"rt_peak",2:"rt_max",3:"rt_min",4:"ppm_tolerance",
                                      5:"group_index",6:"extra_time",7:"label",8:"index"})
    data12.update({"atlas":twelvers})
    data13.update({"atlas":twelvers})
    d12 = ft.get_data(data12,return_data=True,save_file=False)
    d13 = ft.get_data(data13,return_data=True,save_file=False)
    twelvers.sort_values(by=['label'], inplace=True)   
    for i in range(int(len(twelvers)/2)):  
        int13 = d13["ms1_summary"].loc[d13["ms1_summary"]["label"] == (twelvers.iat[i,7] + 1 + len(combined)),:] #13 in 13
        int12 = d12["ms1_summary"].loc[d12["ms1_summary"]["label"] == twelvers.iat[i,7],:] #12 in 12
        twinth = d13["ms1_summary"].loc[d13["ms1_summary"]["label"] == twelvers.iat[i,7],:] #12 in 13
        thintw = d12["ms1_summary"].loc[d12["ms1_summary"]["label"] == (twelvers.iat[i,7] + 1 + len(combined)),:] #13 in 12
        x = 0
        y = 0
        z = 0  
       
        if len(twinth) == 1 and len(int13) > 0:
            if (twinth.iat[0,3]/int13.iat[0,3]) > .5:      
                x = 1
        if len(thintw) == 1 and len(int12) > 0:
            if (thintw.iat[0,3]/int12.iat[0,3]) > .5:
                y = 1
        if len(int13) == 1 and len(int12) == 1:
            bigger = max(int13.iat[0,3],int12.iat[0,3])
            smaller = min(int13.iat[0,3],int12.iat[0,3])
            if bigger/smaller > 5: #5 for qtof, 8 for orbi
                z = 1
        else:
            z = 1
        if (x+y+z) > 0:
            combined.iat[i,0] = np.nan
    combined = combined.dropna(axis = 0)
    combos[key] = combined
print("--- %s seconds ---" % (time.time() - start_time))


#%%
import time
import statistics
import copy

start_time = time.time()

specs = {}
choice = {}
for n in range(len(samples)):
    key = list(combos.keys())[n]
    combined = copy.deepcopy(combos.get(key))
    specs[key] = {}
    choice[key] = {}
    data12['lcmsrun'] = key + '_12C.h5'
    data13['lcmsrun'] = key + '_13C.h5'
    for o in [1,2,3]:   # Can change number to however many 13C + 12C-precursor samples you have
        data['lcmsrun'] = key + '_P' + str(o) + ".h5"
        specs[key]['P' + str(o)] = {}
        choice[key]['P' + str(o)] = {}

        print('Getting isotopologue data for file {}:'.format(o))
        for i in tqdm(range(len(combined))):
            steps = math.ceil(combined.iat[i, 1] - combined.iat[i, 3])
            subst = []
            RT = (combined.iat[i, 0] + combined.iat[i, 2]) / 2
            for j in range(-1, steps + 1):
                subst.append([combined.iat[i, 3] + (j * 1.003586), RT, RT + .1, RT - .1, ppm, j + 2, 0, j + 3, j + 1])
            subst = pd.DataFrame(subst)
            subst = subst.rename(columns={0: "mz", 1: "rt_peak", 2: "rt_max", 3: "rt_min", 4: "ppm_tolerance",
                                          5: "group_index", 6: "extra_time", 7: "label", 8: "index"})
            data12.update({"atlas": subst})
            data13.update({"atlas": subst})
            data.update({"atlas": subst})

            try:
                d = ft.get_data(data, return_data=True, save_file=False)
                d12 = ft.get_data(data12, return_data=True, save_file=False)
                d13 = ft.get_data(data13, return_data=True, save_file=False)
            except:
                # print('Sample: {} o: {} i: {} no data returned'.format(n, o, i))
                continue

            intensity = np.array(d["ms1_summary"]["peak_height"])
            imax = intensity.argmax()
            d = d['ms1_summary']
            d12 = d12['ms1_summary']
            d13 = d13['ms1_summary']
            avgrt = statistics.mean(d['rt_peak'])

            # ### New addition
            label = list(range(2, max(data['atlas']['label'])))
            d_fill = pd.DataFrame(label)
            d_fill = d_fill.rename(columns={0: 'label'})
            d_fill = d_fill.reindex(columns=d.columns)
            for l in range(len(d_fill)):
                hd = d.loc[d['label'] == d_fill.iat[l, 0], :]
                if len(hd) > 0:
                    hd = hd.reset_index(drop=True)
                    d_fill.loc[l] = hd.loc[0]

            d = copy.deepcopy(d_fill)
            d['12'] = 0
            d['13'] = 0
            d['12'] = d['12'].astype(float)
            d['13'] = d['13'].astype(float)

            for l in range(len(d)):
                h12 = d12.loc[d12['label'] == d.iat[l, 0], :]
                if len(h12) > 0:
                    d.loc[l, '12'] = float(h12.iat[0, 3]) 
                    d.loc[l,'mz_centroid'] = float(h12.iat[0,4])
                    d.loc[l,'rt_peak'] = float(h12.iat[0,5])

            for l in range(len(d)):
                h13 = d13.loc[d13['label'] == d.iat[l, 0], :]
                if len(h13) > 0:
                    d.loc[l, '13'] = float(h13.iat[0, 3])
                    d.loc[l,'mz_centroid'] = float(h13.iat[0,4])
                    d.loc[l,'rt_peak'] = float(h13.iat[0,5])

            d['12m/z'] = combined.iat[i, 3]
            d['13m/z'] = combined.iat[i, 1]
            specs[key]['P' + str(o)][str([d.iat[0, 5], d.iat[imax, 3]])] = d[['label', 'mz_centroid', '12', 'peak_height', '13', 'rt_peak', '12m/z', '13m/z']]
            worth = 0
            for k in range(2, len(d) - 2):
                if d.iat[k, 3] > d.iat[k, 6] * 2 and d.iat[k, 3] > d.iat[k, 7] * 2:
                    worth = 1
            if worth == 1:
                choice[key]['P' + str(o)][str([avgrt, intensity[imax]])] = d[['label', 'mz_centroid', '12', 'peak_height', '13', 'rt_peak', '12m/z', '13m/z']]

print("--- %s seconds ---" % (time.time() - start_time))


 #%%save dictionary of raw data

import pickle

with open ("c:/users/mudoe/desktop/InverSIL/d355/finalchoice2.pkl", "wb") as file:
    pickle.dump(choice, file)
    file.close()
     


#%%importing raw file for detected features

import pickle
    
with open ("c:/users/mudoe/desktop/InverSIL/d355/finalchoice2.pkl", "rb") as file:
    loaded_dict = pickle.load(file)
    file.close()

#%%Filtering out false positives
    
import pandas as pd    
AWP = loaded_dict['AWP355']  # Change name as needed. Needs to be the same name chosen for the individual files. 
items_to_delete = []  

# Loop through all keys in the AWP dictionary (e.g., P1, P2, P3)
for key in AWP.keys():
    for dataframe_key, dataframe in AWP[key].items(): 
        dataframe['peak_height'] = dataframe['peak_height'].fillna(0)
        max_value12 = dataframe['12'].max()
        max_valuepeak = dataframe['peak_height'].max()
        max_value13 = dataframe['13'].max()
        
        mzcent12mz = dataframe.loc[dataframe['12'] == max_value12, 'mz_centroid'].values[0] 
        mzcentpeakmz = dataframe.loc[dataframe['peak_height'] == max_valuepeak, 'mz_centroid'].values[0]
        mzcent13mz = dataframe.loc[dataframe['13'] == max_value13, 'mz_centroid'].values[0]
        
        if mzcent12mz < mzcentpeakmz and mzcentpeakmz < mzcent13mz and mzcent12mz < mzcent13mz:
            pass
        else:
            items_to_delete.append((key, dataframe_key))

# Perform the deletion outside the loop
for key, dataframe_key in items_to_delete:
    del AWP[key][dataframe_key]  # Change name as needed. Needs to be the same name chosen for the individual files.

for key in AWP.keys():
    for dataframe_key, dataframe in AWP[key].items():
        dataframe.rename(columns={'12': '12C_intensity', '13': '13C_intensity', 'peak_height': '13C+12C_precursor_intensity'}, inplace=True)




#%%
# save reduced pkl file
with open ("/Users/tashiliebergesell/Desktop/AWP273_reduced.pkl", "wb") as file:
     pickle.dump(loaded_dict, file)
     file.close()  
     
#%% pull out and compile data

data_list = []

for key in ['P3']:
    for feature_key, feature_df in loaded_dict['AWP355'][key].items():  
       
        max_12C_intensity = feature_df['12C_intensity'].max()  # Extracts the maximum value from the 12C_intensity column in the current feature_df dataframe
        max_13C_precursor_intensity = feature_df['13C+12C_precursor_intensity'].max()
        max_13C_intensity = feature_df['13C_intensity'].max()
       
        mz_centroid_12C = feature_df.loc[feature_df['12C_intensity'] == max_12C_intensity, 'mz_centroid'].values[0]  
        mz_centroid_13C_precursor = feature_df.loc[feature_df['13C+12C_precursor_intensity'] == max_13C_precursor_intensity, 'mz_centroid'].values[0]
        mz_centroid_13C = feature_df.loc[feature_df['13C_intensity'] == max_13C_intensity, 'mz_centroid'].values[0]
        RT = feature_df.loc[feature_df['12C_intensity'] == max_12C_intensity, 'rt_peak'].values[0]
        delta_mass = mz_centroid_13C - mz_centroid_13C_precursor
        
        # Append data to the list
        data_list.append({
            'P': key,
            'RT': RT,
            'Mz_Centroid_12C': mz_centroid_12C,
            'Max_12C_intensity': max_12C_intensity,
            'Mz_Centroid_13C+12C_precursor': mz_centroid_13C_precursor,
            'Max_13C+12C_precursor_intensity': max_13C_precursor_intensity,
            'Mz_Centroid_13C': mz_centroid_13C,
            'Max_13C_intensity': max_13C_intensity,
            'Delta_mass': delta_mass
        })

# Convert the list to a DataFrame
combined_data = pd.DataFrame(data_list)

compiled_data = combined_data.sort_values(by=["P", "RT"])

compiled_data.to_csv('c:/users/mudoe/desktop/InverSIL/d355/P3_3.csv')
