import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib import animation
import pandas as pd

from pathlib import *
import os
import importlib

class wave:

    def __init__(self, start_point, dx,dy,c):

        self.start_point = start_point

        # wave pde params
        self._c = c  # velocity of the waves

        # spatial terms
        self._dx = dx
        self._dy = dy

        self.x_range = [0,30]
        self.y_range = self.x_range

        self.cellsize = 1

        #start point
        """
        initial field Q(x,y,t= 0) = 0 for all dimensions except Q((x,y) = (startpoint), t=0) = 1
        """
        self.x_start = start_point[0]
        self.y_start = start_point[1]

        if self.x_start < self.x_range[0] or self.x_start> self.x_range[1] or self.y_start < self.y_range[0] or self.y_start> self.y_range[1]:
            raise ValueError(
                f'({self.x_start,self.y_start}) not between xrange: ({self.x_range}) and/or yrange: ({self.y_range})')

        # x,y points
        self.x_points = np.arange(
            self.x_range[0], self.x_range[1] + self.dx, self.dx, dtype=np.float32)

        self.y_points = np.arange(
            self.y_range[0], self.y_range[1] + self.dy, self.dy, dtype=np.float32)

        self.x_num_points = len(self.x_points)
        self.y_num_points = len(self.y_points)

        #initialised field?
        self.is_initialised = False

        #plot animation
        self._animate_figure = True

        #local variables for iteration
        self._is_iteration = False
        self.iteration_iter = 0
        self._cur_quantities_vector = None
        self._prev_quantities_vector = None
        self._cur_time = None
        self._A = None
        self._initial_field = None

        #dt stability study 
        self.excel_file = "stability_values.xlsx"

        # temporal terms                                                                                                                                                                                                  
        self._dt = self._dt_stability_study() #CFL condition
        self.time_end = 2

        self.time_points = np.arange(0, self.time_end+self.dt, self.dt)
        self.time_num_points = len(self.time_points)

        
    
    
    #for viewing the plot: either animate plots or produce individual plots for each timestep
    @property
    def animate_figure(self):
        return self._clear_figure
    @animate_figure.setter
    def animate_figure(self,bool):
        if bool.lower() == "true":
            self._animate_figure = True
        if bool.lower() == "false":
            self._animate_figure = False
    
    #independent variables for CFL investigation
    @property
    def CFL(self):
        return self._CFL
    @CFL.setter
    def CFL(self,val):
        self._CFL = val
    
    @property
    def dx(self):
        return self._dx
    @dx.setter
    def dx(self,val):
        self._dx = val
    
    @property
    def dy(self):
        return self._dy
    @dy.setter
    def dy(self,val):
        self._dy = val
    
    @property
    def dt(self):
        return self._dt
    @dt.setter
    def dt(self,val):
        self._dt = val
    
    @property
    def c(self):
        return self._c
    @c.setter
    def c(self,val):
        self._c = val


    def run(self):
        field = self.initiliaze_field_()
        print(f'Field initialised. Q(x,y,t = 0)')
        print(f'FDA: {self.explicit_FTCS_h4.__name__}')
        self.explicit_FTCS_h4(field)

    def initiliaze_field_(self):

        field = np.zeros(
            [self.y_num_points, self.x_num_points], dtype=np.float32)

        midpoint = (int(self.y_num_points/2), int(self.x_num_points/2))

        field[midpoint] = 1
        self.is_initialised = True

        return field

    def convert_matrix_to_vector(self, matrix: np.ndarray):
        """
        Converts matrix of dimension (m,n) to vector of dimension (1,mxn) using row major numbering
        """
        y_num, x_num = matrix.shape
        nm_len = int(x_num*y_num)
        vector = np.zeros(nm_len)
        def ix(i, j): return self.ix(i, j)

        for i in range(x_num):
            for j in range(y_num):
                vector[ix(i, j)] = matrix[j][i]
        return vector
    
    def dt_error(self):
        print(f'dx ={self.dx}, dy ={self.dy}, c ={self.c}')

        raise ValueError("\ndt_error, needs a study!")

    def _dt_stability_study(self):
        """
        finds dt according to a stability study 
        """
        
        if self.excel_file in os.listdir(os.getcwd()): #if file already exists

            df = pd.read_excel(self.excel_file)
            print(f'dt_stability:\n{df}')

            if not "c" in df.columns():
                
            if not self.c in df["c"]:
                #if self.c is not in data file --> run another study with the self.c
                self.dt_error()

            #sorts into ascending order
            df.sort_values(by=['dx'], ascending = True, inplace=True)
            #h = df["h"][0]
            
            c_index_range = [i for i in range(len(df)) if df["c"][i]==self.c]
            df["dx_c"] = df["dx"][c_index_range]
            df["dx_diff"] = abs(df["dx"][c_index_range] - self.dx)
            df["dy_diff"] = abs(df["dy"] - self.dy)
            df["diff_sum"] = df["dx_diff"] + df["dy_diff"]
                
            dt_index = df["diff_sum"].idxmin() 
            dt = df["dt"][dt_index]

            return dt * 0.8
        else: #if file doesnt already exists
            self.dt_error()


    def convert_vector_to_matrix(self, vector: np.ndarray):
        """
        inverse of above
        """
        def xi(ix): return self.xi(ix)

        matrix = np.zeros(
            [self.y_num_points, self.x_num_points], dtype=np.float32)

        for ix in range(len(vector)):
            i, j = xi(ix)
            """
            print(f'ix = {ix}')
            print(f"vector[ix] = {vector[ix]}")
            print(f"i j = ({i,j}) ")
            """

            matrix[j][i] = vector[ix]
        return matrix

    def explicit_FTCS_h4(self, initial_field_matrix):
        """
        Explicit Finite difference approximation using:
            - Time: centred difference O(h^2)
            - space: centred difference O(h^4)
        """

        # matrix -  vector conversions
        def ix(i, j): return self.ix(i, j)
        def xi(ix): return self.xi(ix)

        nm_len = int(self.x_num_points*self.y_num_points)

        # initilise matrixs/vectors
        vectors_per_pos_matrix = np.zeros([nm_len, nm_len], dtype=np.float32)

        
        # constants in the finite difference
        tau = (self.c * self.dt)**2
        alpha = 12*self.dx**2
        beta = 12*self.dy**2

        #laplacian coefficients
        Xlap_coef_list = [-tau / alpha, 16*tau /
                            alpha, 16*tau / alpha, -tau / alpha]
        Ylap_coef_list = [-tau / beta, 16 *
                            tau / beta, 16*tau / beta, -tau / beta]
        centre_coeff = (-30/alpha - 30/beta + 2/tau)*tau

        # initial field vector
        initial_field_vector = self.convert_matrix_to_vector(
            initial_field_matrix)
        
        #loop variable init
        t = 0

        #plot variables init
        X,Y = np.meshgrid(self.x_points,self.y_points)

        #label plot function
        def _plot_surface(ax,X,Y,Z,t):
            ax.plot_surface(X, Y, Z, cmap='jet')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Q')
            ax.set_title(f"2D Wave Finite Diff Approx @ t = {round(t,2)}")
            zMax = 2
            ax.set_zlim(-2*zMax,2*zMax)
            return ax
        
        A = vectors_per_pos_matrix.copy()


        if not self._is_iteration:
            while True: 

                # time_ = self.time_points[t]

                if t == 0:  # initialises cur_quantities and prev_quanities
                    _cur_quantities_vector = initial_field_vector.copy()
                    _prev_quantities_vector = np.zeros(nm_len, dtype=np.float32)

                    # initialise plot
                    fig = plt.figure()
                    ax = plt.axes(projection='3d')
                    #plt.show()

                    Z = initial_field_matrix
                    ax = _plot_surface(ax,X,Y,Z,t)

                    plt.show()
                    
                    # plt.pause(10)
                    t += self.dt

                    for i in range(self.x_num_points):
                        for j in range(self.y_num_points):
                            cur_row = ix(i, j)  # cur_row in A

                            # BC --> dirichelet
                            if i in [0, self.x_num_points-1] or j in [0, self.y_num_points-1]:
                                A[cur_row, ix(i, j)] = 0

                            else:  # [(2:xend-2),(2:yend-2)]
                                # non BC points

                                A[cur_row, cur_row] = centre_coeff

                                for i_ in [a for a in range(i-2, i+3) if a >= 0 and a <= self.x_num_points-1]:

                                    if i_ == i:
                                        for j_ in [a for a in range(j-2, j+3) if a > 0 and not a == j and a <= self.y_num_points - 1]:

                                            if j_ > j:
                                                coef_index = j_ - (j - 1)
                                            else:
                                                coef_index = j_ - (j - 2)

                                            A[cur_row, ix(
                                                i, j_)] = Ylap_coef_list[coef_index]
                                    else:
                                        if i_ > i:
                                            coef_index = i_ - (i - 1)
                                        else:
                                            coef_index = i_ - (i - 2)
                                        A[cur_row, ix(i_, j)
                                        ] = Xlap_coef_list[coef_index]

                                        """
                                        A[cur_row,ix(i+2,j)] = -tau/ alpha
                                        A[cur_row,ix(i+1,j)] = 16*tau/ alpha
                                        A[cur_row,ix(i-1,j)] = 16*tau/ alpha
                                        A[cur_row,ix(i-2,j)] = -tau/ alpha

                                        A[cur_row,ix(i,j+2)] = -tau/ beta
                                        A[cur_row,ix(i,j+1)] = 16*tau/ beta
                                        A[cur_row,ix(i,j-1)] = 16*tau/ beta
                                        A[cur_row,ix(i,j-2)] = -tau/ beta
                                        """

                else:
                    # updates the vectors
                    _prev_quantities_vector = _cur_quantities_vector.copy()
                    _cur_quantities_vector = new_time_quantities_vector.copy()

                

                new_time_quantities_vector = A@_cur_quantities_vector - _prev_quantities_vector
                """
                print("A")
                for i in range(len(A)):
                    a = [round(a, 1) for a in A[i]]
                    print(a)

                print("cur_quantities_vector")
                A = self.convert_vector_to_matrix(cur_quantities_vector)
                for i in range(len(A)):
                    a = [round(a, 1) for a in A[i]]
                    print(a)

                print("prev_quantities_vector")
                A = self.convert_vector_to_matrix(prev_quantities_vector)
                for i in range(len(A)):
                    a = [round(a, 1) for a in A[i]]
                    print(a)

                """

                Q_matrix = self.convert_vector_to_matrix(
                    new_time_quantities_vector)
                
                Z = Q_matrix #needed to update the animation

                if self.clear_figure:
                    ax.cla()
                else:
                    plt.figure()
                    ax = plt.axes(projection='3d')

                ax = _plot_surface(ax,X,Y,Z,t)
                plt.draw()
                plt.show()

                plt.pause(2)

                t += self.dt

        else:
            self.iteration_iter+=1
            #print(f'iteration, iter {self.iteration_iter}')

            if not self._cur_time:  # initialises cur_quantities and prev_quanities
                self._A = vectors_per_pos_matrix.copy()
                self._cur_time = 0
                self._cur_quantities_vector = initial_field_vector.copy()
                self._prev_quantities_vector = np.zeros(nm_len, dtype=np.float32)

                self._cur_time += self.dt


                for i in range(self.x_num_points):
                    for j in range(self.y_num_points):
                        cur_row = ix(i, j)  # cur_row in A

                        # BC --> dirichelet
                        if i in [0, self.x_num_points-1] or j in [0, self.y_num_points-1]:
                            A[cur_row, ix(i, j)] = 0

                        else:  # [(2:xend-2),(2:yend-2)]
                            # non BC points

                            A[cur_row, cur_row] = centre_coeff

                            for i_ in [a for a in range(i-2, i+3) if a >= 0 and a <= self.x_num_points-1]:

                                if i_ == i:
                                    for j_ in [a for a in range(j-2, j+3) if a > 0 and not a == j and a <= self.y_num_points - 1]:

                                        if j_ > j:
                                            coef_index = j_ - (j - 1)
                                        else:
                                            coef_index = j_ - (j - 2)

                                        A[cur_row, ix(
                                            i, j_)] = Ylap_coef_list[coef_index]
                                else:
                                    if i_ > i:
                                        coef_index = i_ - (i - 1)
                                    else:
                                        coef_index = i_ - (i - 2)
                                    A[cur_row, ix(i_, j)
                                    ] = Xlap_coef_list[coef_index]


            new_time_quantities_vector = A@self._cur_quantities_vector - self._prev_quantities_vector

            self._prev_quantities_vector = self._cur_quantities_vector.copy()
            self._cur_quantities_vector = new_time_quantities_vector.copy()

            self._cur_time += self.dt



    def explicit_FTCS_h2(self, initial_field_matrix: np.ndarray):
        """
        solves FTCS with spatial accuracy of  O(h^2)
        """

        nm_length = self.x_num_points*self.y_num_points

        vector_per_position = np.zeros(nm_length, dtype=np.float32)

        extra_vals_per_position = vector_per_position.copy()

        def ix(x, y): return self.ix(x, y)
        def xi(ix): return self.xi(ix)

        # methods
        def fill_all_position_matrix(cur_pos_vector, i, j):
            # fills in "all position" matrix
            extra_vals_per_position[ix(
                i, j)] = -1 * (prev_quantities_vector[ix(i, j)])

            all_position_coeff_matrix[ix(i, j)][:] = cur_pos_vector

        # matrix where [ [quantities @all pos @ t = 0], [quantities @all pos @ t = 0 +dt ]..... ]
        quantities_all_time_tensor = np.zeros(
            [self.time_num_points, self.y_num_points, self.x_num_points], dtype=np.float32)
        quantities_all_time_matrix = np.zeros(
            [self.time_num_points, nm_length])
        A = np.eye(nm_length, dtype=np.float32)

        """
        Ax = b
        A --> coeffcients for unknown quantities x
        x --> unknown quantities, refer to quanities in new/next time step
        b --> known quanitiies, refer to quanities in new timestep, which is derived from prev timestep values

        A?
            eye([nm_length,nm_length])
        x?
            vector of shape (1,nm_length)
        b?
            b = B*c + extra_vals
            c --> vector containing quanities at each position at prev timestep, shape = (1,nm_length) 
            B --> matrix containing coefficients of c, which when B*c you get quantity at new timestep given previous time step
            extra_vals --> vector containing extra values from BC or laplace module (space and/or time)
        """
        X, Y = np.meshgrid(self.x_points, self.y_points)
        # Create a 3D plot
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Plot the (x, y, z) data points
        # ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
        ax.plot([self.x_start, self.x_start], [self.y_start,
                self.y_start], [-1.5, 1.5], c="k", marker="*")

        # Add labels
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        plt.show()

        for ti in range(self.time_num_points):
            time = self.time_points[ti]
            print(f'time = {time}')

            # initial field @t =0
            if time == 0:
                if not self.initialised_field:  # ICs
                    raise ValueError(f'Initial Field not passed')

                prev_quantities_vector = self.convert_matrix_to_vector(
                    initial_field_matrix)
                # needs to updated for t>0
                quantities_all_time_matrix[ti][:] = prev_quantities_vector
                quantities_all_time_tensor[ti][:][:] = self.convert_vector_to_matrix(
                    prev_quantities_vector)

                # plot vals
                Z = initial_field_matrix
                ax = fig.add_subplot(111, projection='3d')

                # Plot the (x, y, z) data points
                # ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
                # ax.plot([self.x_start, self.x_start], [self.y_start,self.y_start], [-1.5, 1.5], c="k", marker="*")
                ax.plot_surface(X, Y, Z, cmap='jet')

                # Add labels
                ax.set_xlabel('X Label')
                ax.set_ylabel('Y Label')
                ax.set_zlabel('Z Label')
                ax.set_title(f'time = {time}')
                plt.show()

                continue

            all_position_coeff_matrix = np.zeros(
                [nm_length, nm_length], dtype=np.float32)

            # t> 0
            for j in range(self.y_num_points):
                for i in range(self.x_num_points):
                    x = self.x_points[i]
                    y = self.y_points[j]

                    curr_position_quantities_vector = np.zeros(
                        nm_length, dtype=np.float32)

                    # all BCs
                    if x in self.x_range and y in self.y_range:
                        # corners
                        curr_position_quantities_vector[ix(i, j)] = (self.alpha*2*self.dy + prev_quantities_vector[ix(
                            i, abs(j-2))] + self.alpha*2*self.dx + prev_quantities_vector[ix(abs(i-2), j)]) * 0.5
                        fill_all_position_matrix(
                            curr_position_quantities_vector, i, j)
                        continue
                    elif x in self.x_range:
                        # left right
                        curr_position_quantities_vector[ix(
                            i, j)] = self.alpha*2*self.dx + prev_quantities_vector[ix(abs(i-2), j)]
                        fill_all_position_matrix(
                            curr_position_quantities_vector, i, j)
                        continue
                    elif y in self.y_range:
                        # top bott
                        curr_position_quantities_vector[ix(
                            i, j)] = self.alpha*2*self.dy + prev_quantities_vector[ix(i, abs(j-2))]
                        fill_all_position_matrix(
                            curr_position_quantities_vector, i, j)
                        continue
                    elif x in self.x_range or y in self.y_range:
                        raise ValueError(f'{x,y} in BC but not implemented')

                    # filling in matrix values for new time step
                    """
                    for a new time step, quantity at non-boundary point depends on:
                        5 prev time step quantities
                        1 prev prev time step quantity
                    """
                    alpha = (self.c*self.dt)**2
                    # right point
                    curr_position_quantities_vector[ix(
                        i+1, j)] = alpha/self.dx**2

                    # left point
                    curr_position_quantities_vector[ix(
                        i-1, j)] = alpha/self.dx**2

                    # top point
                    curr_position_quantities_vector[ix(
                        i, j+1)] = alpha/self.dy**2

                    # bott point
                    curr_position_quantities_vector[ix(
                        i, j-1)] = alpha/self.dy**2

                    # centre point
                    curr_position_quantities_vector[ix(
                        i, j)] = 2*(1 - alpha/(self.dx**2) - alpha/(self.dy**2))

                    # references previous time value for position

                    if ti == 1:
                        extra_vals_per_position[ix(
                            i, j)] = -1 * (prev_quantities_vector[ix(i, j)])

                    else:
                        extra_vals_per_position[ix(
                            i, j)] = -1 * (quantities_all_time_matrix[ti-1][ix(i, j)])

                    all_position_coeff_matrix[ix(
                        i, j)][:] = curr_position_quantities_vector

                    """
                    TODO:
                    create a vectors that will be passed into the sparse matrix: 
                    
                    """

                    # matrix_current_time[ix(x_index,y_index),:] = vector_per_position #######upto here!!!
                    # problem: ix(x,y) --> n,m values can be incorrect??

            B = all_position_coeff_matrix.copy()
            c = prev_quantities_vector.copy()
            b = B @ c

            b += extra_vals_per_position  # all known quantities derived here!!

            quantities_new_timestep = np.linalg.inv(A) @ b

            prev_quantities_vector = quantities_new_timestep.copy()
            quantities_all_time_matrix[ti][:] = prev_quantities_vector

            quantities_all_time_tensor[ti][:][:] = self.convert_vector_to_matrix(
                prev_quantities_vector)

            # plot vals
            plt.cla()
            ax = fig.add_subplot(111, projection='3d')

            # Plot the (x, y, z) data points
            # ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
            # ax.plot([self.x_start, self.x_start], [self.y_start,self.y_start], [-1.5, 1.5], c="k", marker="*")

            # Add labels
            ax.set_xlabel('X Label')
            ax.set_ylabel('Y Label')
            ax.set_zlabel('Z Label')
            ax.set_title(f'time = {time}')
            ax.plot_surface(X, Y, self.convert_vector_to_matrix(
                prev_quantities_vector), cmap='jet')

        self.wave_tensor = quantities_all_time_tensor

    def print_all_time_tensor(self, delay_time=3):
        """
        spatial_matrix = np.array([self.x_num_points, self.y_num_points], dtype = np.float32)
        temporal_spatial_matrix = np.array(
            [self.time_num_points, self.x_num_points, self.y_num_points], dtype = np.float32)

        # steps through time, which is the rows, each row represents a single time step
        for v in range(self.time_num_points):
            current_time_step = self.time_points[v]
            current_timestep_matrix = spatial_matrix.copy()
            current_quantities_vector = quantities_all_time[v, :, :]
            # each time step, now fill in quantities for each (x,y)
            for index in range(nm_points):
                current_timestep_matrix[xi(
                    index)] = current_quantities_vector[index]
            # storing current_timestep_quantities matrix into all time tensor
            temporal_spatial_matrix[index_, :, :] = current_timestep_matrix
        """
        quantities_all_time_tensor = self.wave_tensor
        # print the quantities: Q(x,y) over time
        X, Y = np.meshgrid(self.x_points, self.y_points)
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')
        Z = quantities_all_time_tensor[0, :, :]
        ax.plot(self.x_start, self.y_start, 0,
                color="k", marker="*", markersize=30)
        ax.plot_surface(X, Y, Z, cmap='viridis')
        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.set_zlabel('Z-axis')
        ax.set_title(f'3D Surface Plot (t = {self.time_points[0]})')
        ax.set_zlim(-6, 6)

        ax.view_init(elev=90, azim=0)

        if not delay_time:
            delay_time = 3

        iter = 0

        while True:
            plt.pause(delay_time)
            plt.clf()

            iter += 1
            iter = iter % self.time_num_points
            if iter == 0:  # self.time_num_points:
                break

            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(111, projection='3d')
            Z = quantities_all_time_tensor[0, :, :]
            ax.plot(self.x_start, self.y_start, 0,
                    color="k", marker="*", markersize=30)
            ax.plot_surface(X, Y, Z, cmap='viridis')
            ax.set_xlabel('X-axis')
            ax.set_ylabel('Y-axis')
            ax.set_zlabel('Z-axis')
            ax.set_title(f'3D Surface Plot (t = {self.time_points[0]})')
            # ax.plot(self.x_start,self.y_start,0,color = "k", marker= "*", markersize = 30)

            Z = np.clip(quantities_all_time_tensor[iter, :, :], -2, 2)
            # print(f'time = {self.time_points[iter]}')
            # print(Z)
            ax.plot_surface(X, Y, Z, cmap='viridis')

            ax.set_xlabel('X-axis')
            ax.set_ylabel('Y-axis')
            ax.set_zlabel('Z-axis')
            ax.set_title(f'3D Surface Plot (t = {self.time_points[iter]})')
            ax.set_zlim(-6, 6)

    def solve_using_pde_fd(self):
        # row major numbering
        if not self._is_iteration:
            Z_initial = self.initiliaze_field_()
            self.explicit_FTCS_h4(Z_initial)
        # self.print_all_time_tensor()
        if self._is_iteration:
            if not self.is_initialised:
                self._initial_field = self.initiliaze_field_()
            
            self.explicit_FTCS_h4(self._initial_field)
        # self.print_all_time_tensor()

    
    def ix( self,i, j):
        return (j*(self.x_num_points)) + i

    
    def xi( self,ix):
        # remeber its row major numbering !!
        j = ix // (self.x_num_points)
        i = ix % (self.x_num_points)
        return i, j

    def print_Z(self, Z):

        X, Y = np.meshgrid(self.x_points, self.y_points)
        # Create a 3D plot
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Plot the (x, y, z) data points
        # ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
        ax.plot([self.x_start, self.x_start], [self.y_start,
                self.y_start], [-1.5, 1.5], c="k", marker="*")
        # ax.scatter(X, Y, Z, c = Z, cmap = 'jet')

        ax.plot_surface(X, Y, Z, cmap='jet')

        # Add labels
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')

    # DEBUGGING
    def make_csv(self, input):

        csv_file = f'{self.csv_iter}_output_file.csv'

        if isinstance(input, np.ndarray):
            # Save the NumPy array as a CSV file
            np.savetxt(csv_file, input)
        else:

            # Write the matrix to the CSV file
            with open(csv_file, "w", newline="") as csvfile:
                csv_writer = csv.writer(csvfile)
                csv_writer.writerows(input)

        self.csv_iter += 1
        return self.csv_iter - 1

    
    def __iter__(self):
        self._is_iteration = True
        self.solve_using_pde_fd()
        #print(f"iter = {self.iteration_iter}: dx = {self.dx}, dy = {self.dy}, dt ={self.dt}")
        return self
    
    def __next__(self):
        if not self._is_iteration:
            raise StopIteration
        


if __name__ == "__main__":
    dx = dy = 0.5
    c = 1
    start_point = [10, 10]
    sim = wave(start_point, dx,dy,c)
    sim.run()
