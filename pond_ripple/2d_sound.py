import numpy as np
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
import csv

import inspect

import time

import scipy

from matplotlib.animation import FuncAnimation


class wave:

    def __init__(self, start_point):
    
        self.wave_origin = start_point

        # wave terms
        self.A = None
        self.c = None  # wave pde

        # spatial terms
        self.dx = 1
        self.dy = self.dx

        self.x_range = [0, 10]
        self.y_range = [0, 12]

        self.x_start = start_point[0]
        self.y_start = start_point[1]

        if not (self.x_range[1] - self.x_start)< (self.x_range[1]- self.x_range[0]) or not (self.y_range[1] - self.y_start)< (self.y_range[1]- self.y_range[0]):
            raise ValueError(f'({self.x_start,self.y_start}) not between xrange: ({self.x_range}) and/or yrange: ({self.y_range})')
        
        self.x_points = np.arange(self.x_range[0], self.x_range[1] , self.dx, dtype = np.float32)
            
        self.y_points = np.arange(
            self.y_range[0], self.y_range[1], self.dy, dtype = np.float32)

        self.x_num_points = len(self.x_points)
        self.y_num_points = len(self.y_points)

        # temporal terms
        self.dt = 0.1
        self.time_end = 2

        self.time_points = np.arange(0,self.time_end+self.dt,self.dt)
        self.time_num_points = len(self.time_points)

        # 3d sin wave
        self.wave_length = 20
        self.func = lambda A, wave_length, x, y: A * \
            (np.sin(x*np.pi/wave_length)+np.sin(y*np.pi/wave_length))

        # wave pde params
        self.c = 1

        # bc params
        self.alpha = 1  # neumann

        # deugging
        self.csv_iter = 1

        # printing solution params
        self.delay_time = self.dt

        # initialise params
        self.initialised_field = False

        # sim
        self.run_sim()

    def initialise(self):
        self.growth_rate = 1.1
        self.curr_wave_centre = 0

    def bc_check(self, position):
        [x, y] = position
        if not x in [self.x_range[0], self.y_range[0]] or not y in [self.x_range[1], self.y_range[1]]:
            return True
        else:
            raise ValueError(
                f'values passed to {inspect.currentframe().f_back.f_code.co_name} are not boundary values\nValues are: {position}')

    def dirichlet_bc(self, position, value):
        # check if value if boundary:
        if self.bc_check(position):
            return value

    def neumann_bc(self, position, alpha):
        [x, y] = position
        self.bc_check(position)

    def initiliaze_field(self):
        field = np.zeros([self.x_num_points, self.y_num_points])

        # you want a convex wave
        A = 1
        radius = 5
        wave_length = radius*2
        wave_func = lambda x, y: self.func(A, wave_length, x, y)


        x_shifted = lambda i: int(i + self.x_start // self.dx) #returns pos at which x_points[pos] = x_start
        y_shifted = lambda j: int(j + self.y_start // self.dy)



        for i in [i_ for i_ in range(len(self.x_points)) if self.x_points[i_] < radius and self.x_points[i_] > - radius]:
            for j in [j_ for j_ in range(len(self.y_points)) if self.y_points[j_] < radius and self.y_points[j_] > - radius]:

                field[x_shifted(i)][y_shifted(j)] = wave_func(self.x_points[i],
                                     self.y_points[j])
        self.initialised_field = True
        return field
    
    def convert_matrix_to_vector(self,matrix:np.ndarray):
        x_,y_ = matrix.shape
        vector = np.zeros(int(x_*y_))
        ix = lambda x,y: self.ix(x,y)
        for i in range(x_):
            for j in range(y_):
                vector[ix(i, j)] = matrix[i, j]
        return vector

    def explicit_FTCS(self, initial_field_matrix: np.ndarray):

        nm_length = self.x_num_points*self.y_num_points
        matrix_after_init = []

        vector_per_position = np.zeros(nm_length, dtype= np.float32)

        extra_vals_per_position = vector_per_position.copy()

        ix = lambda x, y: self.ix(x, y)

        #methods
        def fill_all_position_matrix(cur_pos_vector,i,j):
            #fills in "all position" matrix
            extra_vals_per_position[ix(
                    i, j)] = -1 * (prev_quantities_vector[ix(i, j)])

            all_position_coeff_matrix[ix(i,j)][:] = cur_pos_vector

        def vector_to_matrix(vector):
            xi = lambda ix: self.xi(ix)

            matrix = np.zeros([self.y_num_points,self.x_num_points])

            for ix in range(len(vector)):
                matrix[xi(ix)] = vector[ix]
    
            return matrix

        # converting initial field condition matrix into vector and storing in all_time matrix
        if self.initialised_field:
            initial_field_vector = self.convert_matrix_to_vector(initial_field_matrix)
        else:
            raise ValueError(f'Initial Field not passed')

        quantities_all_time_tensor = np.zeros([self.time_num_points,self.y_num_points,self.x_num_points],dtype=np.float32) #matrix where [ [quantities @all pos @ t = 0], [quantities @all pos @ t = 0 +dt ]..... ]
        quantities_all_time_matrix = np.zeros([self.time_num_points,nm_length])
        A = np.eye(nm_length, dtype = np.float32)

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
 
        for v in range(self.time_num_points):
            time = self.time_points[v]
            if time == 0: #ICs
                prev_quantities_vector = initial_field_vector
                quantities_all_time_matrix[v][:] = prev_quantities_vector #needs to updated for t>0
                quantities_all_time_tensor[v][:][:] = vector_to_matrix(prev_quantities_vector)
                continue
            
            
            all_position_coeff_matrix = np.zeros([nm_length,nm_length], dtype=np.float32)

            for j in range(self.y_num_points):
                for i in range(self.x_num_points):
                    x = self.x_points[i]
                    y = self.y_points[j]

                    curr_position_quantities_vector = np.zeros(nm_length, dtype = np.float32)
                    
                    #all BCs 
                    if x in self.x_range and y in self.y_range:
                        curr_position_quantities_vector[ix(i, j)] = (self.alpha*2*self.dy + prev_quantities_vector[ix(
                                        i, abs(j-2))] + self.alpha*2*self.dx + prev_quantities_vector[ix(abs(i-2), j)]) * 0.5
                        fill_all_position_matrix(curr_position_quantities_vector,i,j)
                        continue
                    elif x in self.x_range:
                        curr_position_quantities_vector[ix(
                                    i, j)] = self.alpha*2*self.dx + prev_quantities_vector[ix(abs(i-2), j)]
                        fill_all_position_matrix(curr_position_quantities_vector,i,j)
                        continue
                    elif y in self.y_range:
                        curr_position_quantities_vector[ix(
                                    i, j)] = self.alpha*2*self.dy + prev_quantities_vector[ix(i, abs(j-2))]
                        fill_all_position_matrix(curr_position_quantities_vector,i,j)
                        continue
                    elif x in self.x_range or y in self.y_range:
                        raise ValueError(f'{x,y} in BC but not implemented')

                    # filling in matrix values for new time step
                    """
                    for a new time step, quantity at non-boundary point depends on:
                        5 prev time step quantities
                        1 prev prev time step quantity
                    """
                    # right point
                    curr_position_quantities_vector[ix(i+1,j)] = 1 / \
                        np.power(self.dx, 2) * np.power(self.dt, 2) * \
                        np.power(self.c, 2)
                    
                    # left point
                    curr_position_quantities_vector[ix(i-1,j)] = 1 / \
                        np.power(self.dx, 2) * np.power(self.dt, 2) * \
                        np.power(self.c, 2)

                    # top point
                    curr_position_quantities_vector[ix(i,j+1)] = 1 / \
                        np.power(self.dy, 2) * np.power(self.dt, 2) * \
                        np.power(self.c, 2)

                    # bott point
                    curr_position_quantities_vector[ix(i,j-1)] = 1 / \
                        np.power(self.dy, 2) * np.power(self.dt, 2) * \
                        np.power(self.c, 2)

                    # centre point
                    curr_position_quantities_vector[ix(i,j)] = 2 * \
                        (- np.power(self.c*self.dt, 2) *
                         (1 / np.power(self.dx, 2)+1 / np.power(self.dy, 2) + 1))

                    # references previous time value for position
                    extra_vals_per_position[ix(
                        i, j)] = -1 * (prev_quantities_vector[ix(i, j)])

                    all_position_coeff_matrix[ix(i,j)][:] = curr_position_quantities_vector

            
                    """
                    TODO:
                    create a vectors that will be passed into the sparse matrix: 
                    
                    """

                    # matrix_current_time[ix(x_index,y_index),:] = vector_per_position #######upto here!!!
                    # problem: ix(x,y) --> n,m values can be incorrect??
         
            B = all_position_coeff_matrix.copy()
            c = prev_quantities_vector.copy()
            b = B @ c

            b += extra_vals_per_position #all known quantities derived here!!

            quantities_new_timestep = np.linalg.inv(A) @ b

            
            prev_quantities_vector = quantities_new_timestep.copy()
            quantities_all_time_matrix[v][:] = prev_quantities_vector

            quantities_all_time_tensor[v][:][:] = vector_to_matrix(prev_quantities_vector)

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

        # print the quantities: Q(x,y) over time
        X, Y = np.meshgrid(self.x_points, self.y_points)
        fig = plt.figure(figsize = (10,10))
        ax = fig.add_subplot(111, projection='3d')
        Z = quantities_all_time_tensor[0, :, :]
        ax.plot(self.x_start,self.y_start,0,color = "k", marker= "*", markersize = 30)
        ax.plot_surface(X, Y, Z, cmap='viridis')
        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.set_zlabel('Z-axis')
        ax.set_title(f'3D Surface Plot (t = {self.time_points[0]})')
        ax.set_zlim(-6,6)

        delay_time = 3

        iter = 0

        while True:
            plt.pause(delay_time)

            #plt.clf()
            iter+=1
            iter = iter % self.time_num_points

            #ax.plot(self.x_start,self.y_start,0,color = "k", marker= "*", markersize = 30)

            Z = quantities_all_time_tensor[iter, :, :]
            print(Z)
            ax.plot_surface(X, Y, Z, cmap='viridis')
            
            ax.set_xlabel('X-axis')
            ax.set_ylabel('Y-axis')
            ax.set_zlabel('Z-axis')
            ax.set_title(f'3D Surface Plot (t = {self.time_points[iter]})')
            ax.set_zlim(-6,6)

            

    def run_sim(self):
        self.solve_using_pde_fd()

        # self.solve_matrix_using_vector()

    def solve_using_pde_fd(self):
        Z_initial = self.initiliaze_field()
        self.explicit_FTCS(Z_initial)

    def ix(self, i, j):
        return (j*(self.x_num_points-1)) + i + 1 

    def xi(self, ix):
        j = ix // (self.y_num_points -1 ) - 1
        i = ix % (self.x_num_points -1 ) + 1
        return j,i

    def solve_matrix_using_vector(self):

        def x_(x): return x - self.x_start
        def y_(y): return y - self.y_start

        z = np.zeros(self.m*self.n)
        self.A = 1

        for y_index in self.y_positions:
            for x_index in self.x_positions:
                # this is basically how you set quantities to each spatial value
                index = int(self.ix(x_index, y_index))
                z[index] = self.func(self.A, x_(
                    self.x_points[x_index]), y_(self.y_points[y_index]))

        def get_solution_matrix(z_vector):
            z = np.zeros([self.m, self.n])
            for index in range(len(z_vector)):
                z[self.xi(index)] = z_vector[index]

            return z

        Z = get_solution_matrix(z)
        self.print_solution(Z)

    def print_solution(self, Z):

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

        # Show the plot
        plt.show()

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


if __name__ == "__main__":
    start_point = [5, 5]
    wave_sim = wave(start_point)
