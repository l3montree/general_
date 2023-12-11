from _2d_wave import wave
import matplotlib.pyplot as plt
import pandas as pd


""""""
def _CFL_case(dt,dx,dy):
    start_point = [10, 10]
    sim = wave(start_point)
    convert_to_vector = sim.convert_matrix_to_vector
    convert_to_matrix = sim.convert_vector_to_matrix

    sim.dx = dx
    sim.dy=dy
    sim.CFL = 1#CFL
    sim.dt = dt

    for i in range(10):
        iter(sim)

        cur_quantities_vector = sim._cur_quantities_vector

        if abs(max(cur_quantities_vector))> 10:
            return False
        
    return True

dx_list = [5,3,2,1,0.01]
dy_list = dx_list.copy()
CFL_list = dx_list.copy()
dt_list = dx_list.copy()

dt_ = lambda CFL,dx: dx/CFL

data_dict = {"dt":[], "dx":[], "dy":[]}

for dx in dx_list:
    for dy in dy_list:
        for dt in CFL_list:
            if dx == 1 and dy ==1 and dt== 1:
                pass
            bool_ = _CFL_case(dt,dx,dy)
            
            if bool_:
                print(f'pass: dx = {dx}, dy = {dy}, dt = {dt}')
                if dx == dy:
                    data_dict["dt"].append(dt)
                    data_dict["dy"].append(dy)
                    data_dict["dx"].append(dx)
            else:
                print(f'fail: dx = {dx}, dy = {dy}, dt = {dt}')

d = pd.DataFrame(data = data_dict)
d.to_excel("pass.xlsx")

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
dt_max = data_dict["dt"].index(max(data_dict["dt"]))
dt_min = data_dict["dt"].index(min(data_dict["dt"]))
dt_range =(dt_min,dt_max)

ax.plot(data_dict["dx"],data_dict["dy"],data_dict["dt"])
ax.plot(data_dict["dx"][dt_range],data_dict["dy"][dt_range],data_dict["dt"][dt_range])
ax.plot
plt.xlabel("dx")
plt.ylabel("dy")
plt.clabel("dz")

plt.show()

            


        