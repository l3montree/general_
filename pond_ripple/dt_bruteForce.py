from _2d_wave import Wave
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import *
import numpy as np
import pandas as pd
import os


"""
In the instance, the CFL condition is not derived, this serves to numerically determine a relationship between 
dt,dx,dy such that the solution is stable
"""


class BruteForce_dt:
    def __init__(self, dx,dy,c, start_point):

        self.dx = dx
        self.dy = dy
        self.c = c
        self.xRange = self.yRange = [0,100]

        self._use_stored = False
        self._show_plot = False

        self.excel_file = "stability_values.xlsx"
        self._verbose = True

        self.df = None

    #either create a new excel or use the one already available
    @property
    def use_stored(self):
        return self._use_stored
    @use_stored.setter
    def use_stored(self,bool):
        self._use_stored = bool

    #show scatter plot of stable points
    @property
    def show_plot(self):
        return self._show_plot
    @show_plot.setter
    def show_plot(self,bool):
        self._show_plot = bool
    
    #shows all cases, whether pass (stable) or not
    @property
    def verbose(self):
        return self._verbose
    verbose.setter
    def verbose(self,bool):
        self._verbose= bool
    

    def _dt_case_study(self, dt,dx, dy):
        sim = Wave(dx,dy,self.xRange,self.yRange,self.c)
        sim.dt = dt
        sim._is_stability_study = True
        
        for i in range(20):
            iter(sim)

            cur_quantities_vector = sim._cur_quantities_vector

            if abs(max(cur_quantities_vector)) > 2:
                return False

        return True

    def run(self):
        h_step = self.dx/3
        dt_step = self.dx/20

        dt_end = 2

        dx_list = np.arange(self.dx*1.2,self.dx*0.6,-h_step)
        dy_list = dx_list.copy()
        dt_list = np.arange(dt_end,0,-dt_step)

        data_dict = {"dt": [], "dx": [], "dy": [],"h_step":[], "c":[]}

        if not self.use_stored:
            print(f"self.use_stored = False, re-calculating data")
            """
            for dx in dx_list:
                for dy in dy_list:
                    if dx == self.dx and dy == self.dy:
                        print("dx == self.dx and dy == self.dy")
            """
            
            for dt in dt_list:
                #print(f'check: dx = {dx}, dy = {dy}, dt = {dt}')
                dx = self.dx
                dy = self.dy

                bool_ = self._dt_case_study( dt,dx, dy)

                if bool_:

                    if self.verbose:
                        print(f'pass: dx = {dx}, dy = {dy}, dt = {dt}')

                    data_dict["dt"].append(dt)
                    data_dict["dy"].append(dy)
                    data_dict["dx"].append(dx)
                    data_dict["h_step"].append(h_step)
                    data_dict["c"].append(self.c)
                    #break
                    
                else:
                    if self.verbose:
                        print(f'fail: dx = {dx}, dy = {dy}, dt = {dt}')

            d = pd.DataFrame(data=data_dict)

            """
            #consider if you want to append dfs of multiple c together
            if self.excel_file in os.listdir(os.getcwd()):
                print(f'Old file present, reading')
                df_old = pd.read_excel(self.excel_file)
                d = df_old.append(d, ignore_index = True)
                print(f'Old file rewritten')
            else:
                print("No files present, creating new one")
            """
            
            print(f'new study complete')
            d.to_excel(self.excel_file)

        else:
            print(f'self.use_stored = True, using old data file')
            d = pd.read_excel(self.excel_file)

        self.df = d

        if self.show_plot:
            # plot figure
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            ax.scatter(data_dict["dx"], data_dict["dy"], data_dict["dt"])
            # ax.plot(data_dict["dx"][dt_range], data_dict["dy"][dt_range], data_dict["dt"][dt_range])
            plt.xlabel("dx")
            plt.ylabel("dy")
            plt.clabel("dz")

            """
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)

            ax1.plot(data_dict["dx"], data_dict["dt"])

            fig2 = plt.figure()
            ax2 = fig2.add_subplot(111)

            ax2.plot(data_dict["dy"], data_dict["dt"])
            #
            """
            plt.show()

    def show_df(self):
        print(self.df)

if __name__ == "__main__":
    dx =dy =0.5
    c = 1
    start_point = [10, 10]
    sim = BruteForce_dt(dx,dy,c, start_point)
    sim.verbose = True
    sim.run()
    sim.show_df()
