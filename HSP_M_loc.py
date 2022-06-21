import numpy as np
import pandas as pd
from scipy.linalg import pinv
import itertools
import abc #abc is the abstract class
import HSP_SolvP as HSP
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import os

@np.vectorize
def get_distance(d, p, h, ds, ps, hs):
    # x, y, z = cord  #cord is the dph tuple of the studied material.
    # sx, sy, sz = cord_s  #cord_s is the dph tuple of the ith solvent.
    distance = np.sqrt(4*(d-ds)**2+(p-ps)**2+(h-hs)**2)
    return distance

@np.vectorize
def get_weight(indicator):
    weight = indicator
    return weight

class MLoc(HSP.HSPToolkit):

    def __init__(self, input_solv, db):
        super().__init__(input_solv, db)
        self._organize_dphi()

    def get_solv_pool(self):
        mixed_output = []
        mixed_exp = []

        for a in range(self.input_df['Group'].max()):
            G_df = self.input_df[self.input_df['Group']== a+1] #choose rows with same group number
            G_individual = []

            D, P, H = 0, 0, 0
            mixed_name = ''
            
            for index, row in G_df.iterrows():
                solvent = self.db_df[self.db_df['CAS'] == row['CAS']]
                ratio = float(row['Ratio'])
                name = "{:.0%} {}".format(ratio, row['Solvent'])
                G_individual.append([row, solvent])
                D += ratio*float(solvent['D'])
                P += ratio*float(solvent['P'])
                H += ratio*float(solvent['H'])
                indicator = row['Indicator']
                
                mixed_name += name + '+'

            mixed_name = mixed_name[:-1]
            mixed_G = [a+1, D, P, H, G_individual,mixed_name]
            mixed_output.append(mixed_G)

            exp_G = [a+1, D, P, H, indicator]
            mixed_exp.append(exp_G)

        return [mixed_exp, mixed_output]


    def _organize_dphi(self):
        all_dphi = np.array(self.solv_pool[0])
        self.all_index, self.all_ds, self.all_ps, self.all_hs, self.all_ind = all_dphi.T
        
    def init_guess(self, cord = None):
        if cord is not None:
            return cord
        w2 = get_weight(self.all_ind)**2
        sum_w2 = np.sum(w2)
        d = np.sum(w2*self.all_ds)/sum_w2
        p = np.sum(w2*self.all_ps)/sum_w2
        h = np.sum(w2*self.all_hs)/sum_w2
        return d, p, h
    
    def get_delta(self, d, p, h, alpha):
        all_r = get_distance(d, p, h, self.all_ds, self.all_ps, self.all_hs)
        delta_d = -alpha * 4*np.sum(get_weight(self.all_ind) * (d-self.all_ds) / all_r)
        delta_p = -alpha * np.sum(get_weight(self.all_ind) * (p-self.all_ps) / all_r)
        delta_h = -alpha * np.sum(get_weight(self.all_ind) * (h-self.all_hs) / all_r)
        return delta_d, delta_p, delta_h

    def gradient_descent(self, cord_0, alpha = 0.01, n_max = 10000, tol = 0.00001):

        def is_converge(delta_d, delta_p, delta_h):
            if abs(delta_d) < tol and abs(delta_p) < tol and abs(delta_h) <tol:
                return True
            return False
        
        d, p, h = cord_0

        for n in range(n_max):
            delta_d, delta_p, delta_h = self.get_delta(d, p, h, alpha)
            if is_converge(delta_d, delta_p, delta_h):
                print("Converge at {} iteration".format(n))
                return [(d, p, h), n]
            d += delta_d
            p += delta_p
            h += delta_h
        
        print("Warning: Reach maximum iterations. Fail to converge!")
        return [(d, p, h), n]
    
    def get_log(self, cord_0, cord_result, alpha, n_max, tol):
        
        with open(str(self.output_pref) + "log_M_loc.txt", "w") as output:
            output.write("Initial guess = {}\n".format(cord_0))
            output.write("Maxmium iteration = {}\n".format(n_max))
            output.write("Learning rate = {}\n".format(alpha))
            output.write("Tolerance = {}\n".format(tol))
            output.write("Iteration times = {}\n".format(cord_result[1]))
            output.write("D of the material = {}\n".format(cord_result[0][0]))
            output.write("P of the material = {}\n".format(cord_result[0][1]))
            output.write("H of the material = {}\n".format(cord_result[0][2]))
            output.write("T of the material = {}\n".format(np.linalg.norm(cord_result[0])))
    
    def get_compare(self, cord_result, mixed_output):
        M_D, M_P, M_H = cord_result[0]
        M_T = np.sqrt(M_D**2+M_P**2+M_H**2)
        col_out_No = []
        col_solv_name = []
        col_S_D = []
        col_S_P = []
        col_S_H = []
        col_S_T = []
        col_ddD = []
        col_ddP = []
        col_ddH = []
        col_ddT = []
        col_ddR = []

        for in_solv in mixed_output:
            out_No = in_solv[0]
            col_out_No.append(out_No)
            solv_name = in_solv[5]
            col_solv_name.append(solv_name)
            S_D = in_solv[1]
            col_S_D.append(S_D)
            S_P = in_solv[2]
            col_S_P.append(S_P)
            S_H = in_solv[3]
            col_S_H.append(S_H)
            S_T = np.sqrt(S_D**2+S_P**2+S_H**2)
            col_S_T.append(S_T)
            ddD = np.abs(S_D-M_D)
            col_ddD.append(ddD)
            ddP = np.abs(S_P-M_P)
            col_ddP.append(ddP)
            ddH = np.abs(S_H-M_H)
            col_ddH.append(ddH)
            ddT = np.abs(S_T-M_T)
            col_ddT.append(ddT)
            ddR = np.sqrt(4*ddD**2+ddP**2+ddH**2)
            col_ddR.append(ddR)
        
        col_out_No.append('MLoc')
        col_S_D.append(M_D)
        col_S_P.append(M_P)
        col_S_H.append(M_H)
        col_S_T.append(M_T)
        col_ddD.append('')
        col_ddP.append('')
        col_ddH.append('')
        col_ddT.append('')
        min_ddR = np.min(col_ddR)
        min_index = np.argmin(col_ddR)
        col_solv_name.append(self.proc_name+' Best solvent: '+ str(col_out_No[min_index]))
        col_ddR.append(0)



        out_compare_data = {
            'No.' : col_out_No,
            'Name' : col_solv_name,
            'D' : col_S_D,
            'P' : col_S_P,
            'H' : col_S_H,
            'T' : col_S_T,
            'dD' : col_ddD,
            'dP' : col_ddP,
            'dH' : col_ddH,
            'dT' : col_ddT,
            'R' : col_ddR,
        }

        df_compare = pd.DataFrame(data = out_compare_data)
        df_compare_sort = df_compare.sort_values(by=['R'])
        df_compare_sort.to_excel(str(self.output_pref)+'HSP_compare.xlsx', index = 0)


    def get_plot(self, cord_result, c_solv = 'darkgreen', c_m = 'red'):

        mixed_output = self.solv_pool[1]
        import matplotlib as mpl
        from mpl_toolkits.mplot3d import Axes3D, proj3d
        fig = plt.figure(figsize = (12, 8))
        mpl.rcParams['legend.fontsize'] = 15
        mpl.rcParams['axes.formatter.useoffset'] = True
        mpl.rcParams.update({'font.size': 12})
        ax = fig.gca(projection = '3d')

        dm, pm, hm = cord_result[0]
        iter_time = cord_result[1]
        l_m = '{}\n\nD = {:.2f}\nP = {:.2f}\nH = {:.2f}\n'.format(self.proc_name, dm, pm, hm)
        ax.scatter(dm, pm, hm, c = c_m, label = l_m, linewidth = 5)

        for op in mixed_output:
            individual = op[4]

            D = [op[1]]
            P = [op[2]]
            H = [op[3]]

            x, y, z = D, P, H

            ax.scatter(x, y, z, c = c_solv, linewidth = 5)

        ax.set_xlabel(r'D(MPa$^{1/2}$)', fontsize = 15, labelpad = 8)
        ax.set_ylabel(r'P(MPa$^{1/2}$)', fontsize = 15, labelpad = 10)
        ax.set_zlabel(r'H(MPa$^{1/2}$)', fontsize = 15, labelpad = 8)

        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False

        ax.legend(loc = 'best', bbox_to_anchor = (-0.08, 0, 0, 0), bbox_transform=ax.transData) #If a 4-tuple or BboxBase is given, then it specifies the bbox (x, y, width, height) that the legend is placed in. To put the legend in the best location in the bottom right quadrant of the axes (or figure)

        #the next step is to plot the projection to x,y,z plane
        x_lim = ax.get_xlim3d()
        y_lim = ax.get_ylim3d()
        z_lim = ax.get_zlim3d()

        ax.scatter(dm, pm, z_lim[0], c = c_m, marker = '^', alpha = 0.2)
        ax.scatter(x_lim[0], pm, hm, c = c_m, marker = 'x', alpha = 0.2)
        ax.scatter(dm, y_lim[1], hm, c = c_m, marker = 's', alpha = 0.2)

        for op in mixed_output:
            individual = op[4]

            D = [op[1]]
            P = [op[2]]
            H = [op[3]]

            x, y, z = D, P, H
            ax.scatter(x, y, z_lim[0], c = c_solv, marker = '^', alpha = 0.2)
            ax.scatter(x_lim[0], y, z, c = c_solv, marker = 'x', alpha = 0.2)
            ax.scatter(x, y_lim[1], z, c = c_solv, marker = 's', alpha = 0.2)

        plt.tight_layout()

        plt.savefig(str(self.output_pref) + '3D_M_Loc.png')
        
    def run_all(self, cord = None, alpha = 0.01, n_max = 10000, tol = 0.00001):
        initial_guess = self.init_guess(cord = cord)
        print(initial_guess)
        cord_result = self.gradient_descent(cord_0 = initial_guess, alpha = alpha, n_max = n_max, tol = tol)
        self.get_plot(cord_result)
        self.get_log(initial_guess, cord_result, alpha = alpha, n_max = n_max, tol = tol)
        mixed_output = self.get_solv_pool()[1]
        self.get_compare(cord_result, mixed_output)
        return cord_result

# mop = MLoc(r'input_Safa_SA45_0520.xlsx', r'db.xlsx')
# cord_result = mop.run_all(alpha = 0.001, n_max = 100000, tol = 0.0001)
# print(cord_result)

# predictor = HSP.SolvPredictor(r'input_solv_sel.xlsx', r'db.xlsx')
# predictor.run_all(3, cord_result[0][0], cord_result[0][1], cord_result[0][2], rep_time = 50, std = 0.1, tol = 0.1, red_tol = 0.01)

# import shutil
# ori_path = predictor.folder_path
# target_path = mop.folder_path
# # shutil.move(ori_path, target_path)

# class MLocSolvPred():

#     def __init__(self, exp_result, input_solv, db):
#         self.mop = MLoc(exp_result, db)
#         print(self.mop.folder_path)
#         self.pred = HSP.SolvPredictor(input_solv, db, folder_name = self.mop.folder_path)
    
#     def run_all(self, alpha = 0.001, n_max = 1000000, tol_mop = 0.0005, n = 3, rep_time = 50, std = 0.1, tol_pred = 0.1, red_tol = 0.01):
#         cord_result = self.mop.run_all(alpha = alpha, n_max = n_max, tol = tol_mop)
#         self.pred.run_all(n, cord_result[0][0], cord_result[0][1], cord_result[0][2], rep_time = rep_time, std = std, tol = tol_pred, red_tol = red_tol)

# mp = MLocSolvPred(r'input_Basi_PPI2_0523.xlsx', r'input_solv_sel.xlsx', r'db.xlsx')
# mp.run_all()