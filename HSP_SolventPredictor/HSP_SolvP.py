import numpy as np
import pandas as pd
from scipy.linalg import pinv
import itertools
import abc
import os

class HSPToolkit(abc.ABC):
    """
    This is an abstract class
    """
    def __init__(self, input_solv, db, folder_name = ''):
        self.input_df = pd.read_excel(input_solv)
        self.db_df = pd.read_excel(db)
        self.solv_pool = self.get_solv_pool()
        self.input_name = str(input_solv)
        self.proc_name = self.input_name_proc()
        self.folder_path = self.new_folder_path(folder_name)
        self.check_dir()
        self.output_pref = self.path_n_pref()
    
    def check_dir(self):
        if os.path.exists(self.folder_path):
            return True
        else:
            os.mkdir(self.folder_path)


    def input_name_proc(self):
        raw_name = self.input_name
        name_no_ext = list(raw_name.split('.'))
        name_no_input = list(name_no_ext[0].split('_'))
        name_no_input.remove('input')
        proc_name = '_'.join(name_no_input)
        return proc_name
    
    def new_folder_path(self, folder_name):
        
        if not folder_name:
            folder_name = self.input_name_proc()
            current_path = os.getcwd()
            folder_name = self.input_name_proc()
            comb_path_list = [current_path, folder_name]
            print(comb_path_list)
            new_folder_path = '\\'.join(comb_path_list) + '\\'
        else:
            new_folder_path = folder_name

        return new_folder_path
    
    def path_n_pref(self):
        new_path = self.folder_path
        name_pref = self.input_name_proc()
        output_pref = new_path + name_pref +'_'
        return output_pref

    @abc.abstractmethod
    def get_solv_pool(self):
        pass


class SolvPredictor(HSPToolkit):
    """
    This is a solvent mixer.
    """

    def __init__(self, input_solv, db, folder_name = ''):
        super().__init__(input_solv, db, folder_name = folder_name)
        self.n_solvent = len(self.input_df)

    def get_solv_pool(self):

        all_input_solv = []

        for index, input_row in self.input_df.iterrows():
            db_row = self.db_df[self.db_df['CAS'] == input_row['CAS']]
            both_row = [input_row, db_row]
            D = float(db_row['D'])
            P = float(db_row['P'])
            H = float(db_row['H'])
            solv_full_info = [index, D, P, H, both_row]
            all_input_solv.append(solv_full_info)
        
        return all_input_solv
    
    def init_s_tot(self):

        s_tot = np.ones((4, self.n_solvent)) #s_tot is a 4*n_solvent matrix.

        for solv in self.solv_pool:
            solv_ind = solv[0]
            s_tot[:-1,solv_ind] = solv[1:-1]
            
        return s_tot
    
    def solv_comb(self, n):

        all_solv_comb = itertools.combinations(range(self.n_solvent),n)

        return all_solv_comb
    
    def target_mat(self, D, P, H, rep_time = 50, std = 0.1):
        
        ds_HSP = [D, P, H, 1]
        ptb_mat = np.random.randn(rep_time, 4)*(std, std, std, 0) #perturb
        ds_mat = ptb_mat + np.array(ds_HSP) * np.ones((rep_time, 4))

        return ds_mat
    
    def calc_unit(self, comb, s_tot, ds_mat, n):
        
        s_temp = np.zeros((4,n))

        out_solv = []
        out_solv_no = []
        e_list = []
        
        for i, k in enumerate(comb):
            s_temp[:,i] = s_tot[:,k]  #the ith column of s_temp is the kth solvent in s_tot.
            sel_solv = self.solv_pool[k][4][0]['Solvent']
            sel_solv_no = self.solv_pool[k][4][0]['No_db']
            out_solv.append(sel_solv)
            out_solv_no.append(sel_solv_no)
        
        c = pinv(s_temp) @ ds_mat.T
        c_mean = list(np.mean(c, axis = 1))
        c_std = list(np.std(c, axis = 1))

        real_hsp = s_temp @ c_mean

        e = s_temp @ c - ds_mat.T
        e_mean = np.mean(e, axis = 1)
        e_std = np.std(e, axis = 1)

        for j, mean in enumerate(e_mean):
            e_list += [mean, e_std[j]]
        
        return out_solv + c_mean + c_std + e_list + list(real_hsp) + out_solv_no
    
    def get_col_name(self, n):

        col_name = np.zeros((3, n), dtype = object)
        for i in range(n):
            col_name[0, i] = f'Solvent{i+1}'
            col_name[1, i] = f'c_mean{i+1}'
            col_name[2, i] = f'c_std{i+1}'


        col_name = list(col_name.reshape(3*n))
        for i in 'DPHI':
            col_name.append(f'e_mean_{i}')
            col_name.append(f'e_std_{i}')

        for i in 'DPHI':
            col_name.append(i)

        for i in range(n):
            col_name.append(f'no_db_{i+1}')
        
        return col_name

    def rough_calc(self, n, D, P, H, rep_time = 50, std = 0.1):

        s_tot = self.init_s_tot()
        all_solv_comb = self.solv_comb(n)
        ds_mat = self.target_mat(D, P, H, rep_time = rep_time, std = std)

        calc_proc = []
        for comb in all_solv_comb:
            temp_calc_proc = self.calc_unit(comb, s_tot, ds_mat, n)
            calc_proc.append(temp_calc_proc)
        
        col_name = self.get_col_name(n)
        
        df_total = pd.DataFrame(calc_proc, columns = col_name)

        return df_total
    
    def is_rough_valid(self, n, df_row):
        """
        To check if any coefficient is negative or more than 100
        """
        if any(df_row[n:2*n] < 0) or any(df_row[n:2*n] > 1):
            return False
        """
        To validate the stability of perturbance
        """
        if abs(df_row['e_mean_I']) > 0.05:
            return False
        if abs(df_row['e_mean_D']) > 0.1 and abs(df_row['e_mean_P']) > 0.1 and abs(df_row['e_mean_H']) > 0.1:
            return False
        if any(abs(df_row[2*n:3*n] > 0.1)):
            return False

        return True
    
    def get_stable_result(self, n, df_total):
        """
        This is to get the intermediate results by filtering unvalid coefficients and stability. This reslt is before 100% correction.
        """
        df_new = df_total
        for index, row in df_new.iterrows():
            if not self.is_rough_valid(n, row):
                df_new = df_new.drop(index)
        
        df_new.to_excel(str(self.output_pref) + 'stable_result.xlsx', index = 0)
        return df_new
    
    def update_c(self, n, df_row):
        rough_c_list = []
        total_rough_c = 0
        fine_c_list = []

        for rough_c in df_row[n:2*n]:
            rough_c_list.append(rough_c)
        
        total_rough_c = sum(rough_c_list)
        
        for rough_c in df_row[n:2*n]:
            fine_c = rough_c/total_rough_c
            fine_c_list.append(fine_c)
        
        df_row[n:2*n] = fine_c_list

        return df_row


    def update_hsp(self, n, df_row):
        new_c_list = []
        new_hsp_list = []
        ori_hsp_list = [] #the original hsp of individual solvents

        for new_c in df_row[n:2*n]:
            new_c_list.append(new_c)

        new_D = 0
        new_P = 0
        new_H = 0

        for i, old_solv_no in enumerate(df_row[-n:]):
            k = old_solv_no
            ori_D = float(self.db_df[self.db_df['No.']==k]['D'])
            ori_P = float(self.db_df[self.db_df['No.']==k]['P'])
            ori_H = float(self.db_df[self.db_df['No.']==k]['H'])            
            new_D += new_c_list[i]*ori_D
            new_P += new_c_list[i]*ori_P
            new_H += new_c_list[i]*ori_H

        df_row['D'] = new_D
        df_row['P'] = new_P
        df_row['H'] = new_H

        return df_row

    def update_e(self, df_row, D, P, H):
        new_e = [df_row['D']-D, df_row['P']-P, df_row['H']-H]

        df_row['e_mean_D'] = new_e[0]
        df_row['e_mean_P'] = new_e[1]
        df_row['e_mean_H'] = new_e[2]

        return df_row

    def error_valid(self, df_row, tol=1):
        if abs(df_row['e_mean_D']) > tol or abs(df_row['e_mean_P']) > tol or abs(df_row['e_mean_H']) > tol:
            return False
        return True
    
    def get_norm_result(self, n, df_new, D, P, H):

        df_new = df_new.drop(columns = ['I','e_std_D', 'e_std_P', 'e_std_H', 'e_mean_I', 'e_std_I'], axis = 1) #remove I
        df_new = df_new.drop(df_new.iloc[:, 2*n:3*n], axis = 1) #remove the c_std columns

        for index, row in df_new.iterrows():
            df_row = self.update_c(n, row) #update the coefficients of each solvent into the normalized   
            df_row = self.update_hsp(n, df_row) #update the new hsps after calculating the normalized c
            df_row = self.update_e(df_row, D, P, H)#update the new errors of each hsp
            df_new.loc[index] = df_row
        
        return df_new

    def error_filter(self, df_new, tol = 1):
        for index, row in df_new.iterrows():
            if not self.error_valid(row, tol= tol):
                df_new = df_new.drop(index)

        #print(df_new.head())
        df_new.to_excel(str(self.output_pref) + 'error_filtered.xlsx', index = 0)
        return df_new
    
    def get_red_list(self, n, df_row, red_tol = 0.01):
        """
        This function is to check if there is any redundant solvents, which means the percentage is less than the redundant tolerance red_tol.
        The default threshold of red_tol is 0.01
        """
        red_n_list = []
        for order, c in enumerate(df_row[n:2*n]):
            if c <= red_tol:
                red_n_list.append(order) #get the order of red solvent
        return red_n_list

    def set_red(self, n, df_row, red_order):
        """
        This function is to remove redundant elements from the current row
        """
        df_row[red_order] = ''
        #ori_red_c = df_row.loc[:, n + red_order] #save the original value of red solvent percentage
        df_row[n + red_order] = 0
        #df_row.loc[:, -(n-red_order)] = 0

        return df_row

    def red_filter(self, D, P, H, n, df_new, red_tol = 0.01):
        for index, row in df_new.iterrows():  
            red_solv_order_list = self.get_red_list(n, row, red_tol = red_tol)
            n_of_red_solv = len(red_solv_order_list)
            if red_solv_order_list:
                for red_solv_order in red_solv_order_list:
                    df_row = self.set_red(n, row, red_solv_order)#remove the redundant solvents info, the redundant solvents will have 0 as its coefficient
                df_row = self.update_c(n, df_row)#normalize again
                df_row = self.update_hsp(n, df_row)#update hsp
                df_row = self.update_e(df_row, D, P, H)#update e
                df_new.loc[index] = df_row
            
        for index, row in df_new.iterrows():
            if not self.error_valid(row, tol=1):
                df_new = df_new.drop(index)

        temp_solv_name_list = []
        rm_index_list = []

        for index, row in df_new.iterrows():
            solv_name_list = list(row[0:n])
            if not all(solv_name_list): #empty elements in the list
                if solv_name_list in temp_solv_name_list:
                    rm_index_list.append(index) #if the current solvent combination is already in the temp list, record the current index of row, which is going to be removed
                else:
                    temp_solv_name_list.append(solv_name_list)

        df_new = df_new.drop(rm_index_list)

        if df_new.empty:
            print("No solvent matched. Please increase tolerence.")
            
        df_new.to_excel(str(self.output_pref) + 'Final_result.xlsx', index = 0)
    
    
    def get_log(self, n, D, P, H, rep_time, std, tol, red_tol):
        
        with open(str(self.output_pref) + "log_SolvPredict.txt", "w") as output:
            output.write("Solvent amount = {}\n".format(n))
            output.write("Target HSPs:\n")
            output.write("D = {}\n".format(D))
            output.write("P = {}\n".format(P))
            output.write("H = {}\n".format(H))
            output.write("Perturbation parameters:\n")
            output.write("Repeat time = {}\n".format(rep_time))
            output.write("std = {}\n".format(std))
            output.write("Filter parameters:\n")
            output.write("Tolerance of error = {}\n".format(tol))
            output.write("Tolerance of redundant = {}\n".format(red_tol))
            output.write("See Final_result.xlsx as predicted solvents combinations.")
    
    def run_all(self, n, D, P, H, rep_time = 50, std = 0.1, tol = 1, red_tol = 0.01):
        df_total = self.rough_calc(n, D, P, H, rep_time = rep_time, std = std)
        df_new = self.get_stable_result(n, df_total)
        df_new = self.get_norm_result(n, df_new, D, P, H)
        df_new = self.error_filter(df_new, tol = tol)
        self.red_filter(D, P, H, n, df_new, red_tol = red_tol)
        self.get_log(n, D, P, H, rep_time = rep_time, std = std, tol = tol, red_tol = red_tol)

# predictor = SolvPredictor(r'input_solv_sel.xlsx', r'db.xlsx')
# predictor.run_all(3, 18.5, 13.47, 10.2)