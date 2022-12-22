#########calculates flamespeeds for the .cti files at different equivalence ratios. Uses initial guess from the previous model #############################


#####will be used to calculate the burning velocity of Methane-Air mix at an equivalence ratio = 1, with increasing volume 
#####fraction of 2-BTP (0-3%)

import cantera as ct
import numpy as np
import pandas as pd
import os
import csv 
import sys
import re


print("Running Cantera Version: " + str(ct.__version__))

path_ = sys.argv[1]

file_name = path_.split('/')[-1]


To = 298
Po = ct.one_atm

directory = path_


gas = ct.Solution(directory)

BTPmole_list = list(np.linspace(0.0, 0.16268, 10))#### for this run, this is not really volume fraction. This value is the moles of BTP (normalized to oxygen) volume fraction can be found on line 45 


#can use below to only test one vol_frac
#vol_frac_list = [0.095] 
#0.095 is the volume fraction of METHANE that leads to an equivalence ratio = 1 


results = {}

for i in  range(len(BTPmole_list)):
    try: 
        
        tol_ss = [1.0e-13, 1.0e-9]  #abs and rel tolerances for steady state problem
        tol_ts = [1.0e-13, 1.0e-9]  #abs and rel tie tolernces for time step function
        
        x = BTPmole_list[i]
        norm_ox = (1-x)*.21
        
        
        print(f'****************************starting new volume fraction: {x}**************************')

        vol_frac_dict = {'CH4(3)': 0.5, 'O2(4)': 1, 'N2': 3.76, '2-BTP(1)': x }
        #{'CH4': 0.4998684556695607, 'O2': 1.0, 'N2': 3.7619047619047623} is when eq ratio = 1 with no suppressant
        #volume fraction of BTP = x/(x+1+3.76+0.5). when BTP vol frac is 3%, moles of BTP (normalized to O2) is 0.16268
        print(vol_frac_dict)
        print(f"O2/CH4 ratio = {vol_frac_dict['O2(4)']/vol_frac_dict['CH4(3)']}. Complete combustion takes 2")
        gas.TPX = To, Po, vol_frac_dict
        width = 0.08
        flame = ct.FreeFlame(gas, width=width)
        flame.flame.set_steady_tolerances(default=tol_ss)   #set tolerances
        flame.flame.set_transient_tolerances(default=tol_ts)
        #flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1) 
        flame.set_refine_criteria(ratio=5, slope=0.25, curve=0.27)
        flame.max_time_step_count = 900
        loglevel = 1 
        
        ######################################################################################
#         if i!=0:
#             d = f'./data/range10pts_test_{BTPmole_list[i-1]}.csv'
#             if os.path.exists(d):  
#                 arr2 = ct.SolutionArray(gas)
#                 arr2.read_csv(d)
#                 flame.set_initial_guess(data=arr2)
#                 print(' initial guess has been set')
        #######################################################################################        
                
        #"False" stops the calculation from retrying over and over, thanks Chao 
        #flame.solve(loglevel=loglevel, auto=False)
        flame.solve(loglevel=loglevel, auto=True)
        Su = flame.velocity[0]
        vol_frac_of_btp = x/(x+1+3.76+0.5)
        results[vol_frac_of_btp] = Su
        sltn = flame.to_solution_array()
        df1 = sltn.to_pandas()
        #edited this here!! index=False
        df1.to_csv(f'{file_name}_{x}.csv', index=False)
    except Exception as e: 
        print(f'********************   passed BTP mole:{BTPmole_list[i]}, error: {e}    *************************************')
        pass


vol_fracs_of_BTP = list(results.keys())
flame_speeds = list(results.values())


print("volume fractions of BTP are:")
print(vol_fracs_of_BTP)

print("flame speeds are:")
print(flame_speeds)


with open(f'./data/speeds_{file_name}', 'w+') as g:
    g.write(directory)
    g.write('\n')
    writers = csv.writer(g)
    writers.writerow(vol_fracs_of_BTP)
    writers.writerow(flame_speeds)

        

