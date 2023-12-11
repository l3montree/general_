from _2d_wave import wave
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import *
import sys


"""
In the instance, the CFL condition is not derived, this serves to numerically determine a relationship between 
dt,dx,dy such that the solution is stable
"""


class CFL_investigation:
    def __init__(self):
        self.use_stored = True

    def _CFL_case(self, dt, dx, dy):
        start_point = [10, 10]
        sim = wave(start_point)
        convert_to_vector = sim.convert_matrix_to_vector
        convert_to_matrix = sim.convert_vector_to_matrix

        sim.dx = dx
        sim.dy = dy
        sim.CFL = 1  # CFL
        sim.dt = dt

        for i in range(10):
            iter(sim)

            cur_quantities_vector = sim._cur_quantities_vector

            if abs(max(cur_quantities_vector)) > 10:
                return False

        return True

    def run(self):

        dx_list = [range(5, 1), 0.01]
        dy_list = dx_list.copy()
        CFL_list = dx_list.copy()
        dt_list = dx_list.copy()

        file = "pass.xlsx"
        # sys.path.append(file)

        def dt_(CFL, dx): return dx/CFL

        data_dict = {"dt": [], "dx": [], "dy": []}

        if not self.use_stored:
            for dx in dx_list:
                for dy in dy_list:
                    for dt in CFL_list:
                        if dx == 1 and dy == 1 and dt == 1:
                            pass
                        bool_ = self._CFL_case(dt, dx, dy)

                        if bool_:
                            print(f'pass: dx = {dx}, dy = {dy}, dt = {dt}')
                            # if dx == dy:
                            data_dict["dt"].append(dt)
                            data_dict["dy"].append(dy)
                            data_dict["dx"].append(dx)
                        else:
                            print(f'fail: dx = {dx}, dy = {dy}, dt = {dt}')

            d = pd.DataFrame(data=data_dict)
            d.to_excel(file)

        else:
            data_dict = pd.read_excel(file)

        # plot figure
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        ax.scatter(data_dict["dx"], data_dict["dy"], data_dict["dt"])
        # ax.plot(data_dict["dx"][dt_range], data_dict["dy"][dt_range], data_dict["dt"][dt_range])
        plt.xlabel("dx")
        plt.ylabel("dy")
        plt.clabel("dz")

        plt.show()


if __name__ == "__main__":
    sim = CFL_investigation()
    sim.run()
