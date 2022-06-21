import numpy as np
import pandas as pd
from scipy.linalg import pinv
import itertools
import abc #abc is the abstract class
import HSP_reconstruct as HSP

class SolvMixer(HSP.HSPToolkit):
    def get_solv_pool(self):
        mixed_output = []
        for a in range(self.input_df['Group'].max()):
            G_df = self.input_df[self.input_df['Group']== a+1] #choose rows with same group number
            G_individual = []
            mix_G = []
            D, P, H = 0, 0, 0
            mixed_name = ''
            
            for index, row in G_df.iterrows():
                solvent = self.db_df[self.db_df['CAS'] == row['CAS']]
                print(row['CAS'])
                ratio = float(row['Ratio'])
                name = "{:.0%} {}".format(ratio, row['Solvent'])
                G_individual.append([row, solvent])
                D += ratio*float(solvent['D'])
                P += ratio*float(solvent['P'])
                H += ratio*float(solvent['H'])
                mixed_name += name + '+'
            mixed_name = mixed_name[:-1]
            mixed_G = [a+1, D, P, H, G_individual]
            mixed_output.append(mixed_G)
        
        return mixed_output