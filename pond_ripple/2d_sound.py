import numpy as np
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
import csv

import inspect

import time

import scipy


class wave:

    def __init__(self, start_point):
        self.x_start = start_point[0]
        self.y_start = start_point[1]

        self.wave_origin = start_point

        # wave terms
        self.A = None
        self.c = None  # wave pde

        # spatial terms
        self.dx = 0.5
        self.dy = self.dx

        self.x_range = [0, 100]
        self.y_range = [0, 100]

        self.x_num_points = int((self.x_range[1]-self.x_range[0])/self.dx)+1
        self.y_num_points = int((self.y_range[1]-self.y_range[0])/self.dy)+1

        self.x_points = np.linspace(
            self.x_range[0], self.x_range[1], self.x_num_points)
        self.y_points = np.linspace(
            self.y_range[0], self.y_range[1], self.y_num_points)

        self.x_positions = range(self.x_num_points-1)
        self.y_positions = range(self.y_num_points-1)

        # temporal terms
        self.dt = 0.5
        self.time_end = 10

        self.time_points = [
            t/10 for t in range(0, self.time_end*10, int(self.dt*10))]
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
        wave_length = 4
        radius = 5
        def func(x, y): return self.func(A, wave_length, x, y)

        x_curr = self.x_start
        y_curr = self.y_start

        def x_shifted(x): return int(x + self.x_start // self.dx)
        def y_shifted(y): return int(y + self.y_start // self.dy)

        for x_index in [int(i) for i in self.x_positions if self.x_points[i] < radius and self.x_points[i] > -radius]:
            for y_index in [int(j) for j in self.y_positions if self.y_points[j] < radius and self.y_points[j] > -radius]:
                i_ = x_shifted(x_index)
                j_ = y_shifted(y_index)

                field[i_, j_] = func(self.x_points[x_index],
                                     self.y_points[y_index])

        self.initialised_field = True
        return field

    def row_matrix_pde_solve(self, initial_field_matrix: np.ndarray):

        nm_length = self.x_num_points*self.y_num_points
        matrix_all_time = np.zeros([self.time_num_points, nm_length])

        matrix_current_time = np.zeros([self.x_num_points, self.y_num_points])
        vector_per_position = np.zeros([self.x_num_points*self.y_num_points])

        extra_vals_per_position = vector_per_position.copy()

        def ix(x, y): return self.ix(x, y)
        def xi(ix_): return self.xi(ix_)

        # converting initial field condition matrix into vector and storing in all_time matrix
        initial_field_vector = np.zeros(nm_length)
        if self.initialised_field:
            for i in self.x_positions:
                for j in self.y_positions:
                    initial_field_vector[ix(i, j)] = initial_field_matrix[i, j]

            matrix_all_time[0, :] = initial_field_vector
        else:
            raise ValueError(f'Initial Field not passed')

        x_values = []
        y_values = []
        value = []

        for i in range(1, self.time_num_points):
            time = self.time_points[i]
            current_quantities_vector = matrix_all_time[i-1]
            A = scipy.sparse.eye(nm_length)
            matrix_current_time = A.copy()
            for y_index in self.y_positions:
                for x_index in self.x_positions:

                    # x,y point being used in current loop
                    x = self.x_points[x_index]
                    y = self.y_points[y_index]

                    # skips all BCs and ICs
                    if x in self.x_range:
                        continue
                    if y in self.y_range:
                        continue
                    if time == 0:
                        # time approx requires the previous time step, would not exist at t = 0
                        continue  # IC

                    # filling in matrix values for new time step
                    """
                    for a new time step, quantity at non-boundary point depends on:
                        5 prev time step quantities
                        1 prev prev time step quantity
                    """
                    # right point
                    value_at_index = 1 / \
                        np.power(self.dx, 2) * np.power(self.dt, 2) * \
                        np.power(self.c, 2)
                    value.append(value_at_index)
                    x_values.append(ix(x_index+1, y_index))
                    y_values.append(self.y_points[y_index])

                    # left point
                    value_at_index = 1 / \
                        np.power(self.dx, 2) * np.power(self.dt, 2) * \
                        np.power(self.c, 2)
                    value.append(value_at_index)
                    x_values.append(ix(x_index-1, y_index))
                    y_values.append(self.y_points[y_index])

                    # top point
                    value_at_index = 1 / \
                        np.power(self.dy, 2) * np.power(self.dt, 2) * \
                        np.power(self.c, 2)
                    value.append(value_at_index)
                    x_values.append(ix(x_index, y_index+1))
                    y_values.append(self.y_points[y_index])

                    # bott point
                    value_at_index = 1 / \
                        np.power(self.dy, 2) * np.power(self.dt, 2) * \
                        np.power(self.c, 2)
                    value.append(value_at_index)
                    x_values.append(ix(x_index, y_index-1))
                    y_values.append(self.y_points[y_index])

                    # centre point
                    value_at_index = 2 * \
                        (- np.power(self.c*self.dt, 2) *
                         (1 / np.power(self.dx, 2)+1 / np.power(self.dy, 2) + 1))
                    value.append(value_at_index)
                    x_values.append(ix(x_index, y_index))
                    y_values.append(self.y_points[y_index])

                    """
                    vector_per_position[ix(x_index+1,y_index)] = 1/ np.power(self.dx,2) * np.power(self.dt,2) * np.power(self.c,2)
                    vector_per_position[ix(x_index-1,y_index)] = 1/ np.power(self.dx,2) * np.power(self.dt,2) * np.power(self.c,2)
        
                    vector_per_position[ix(x_index, y_index + 1)] = 1/ np.power(self.dy,2) * np.power(self.dt,2) * np.power(self.c,2)
                    vector_per_position[ix(x_index, y_index - 1 )] = 1/  np.power(self.dy,2) * np.power(self.dt,2) * np.power(self.c,2)
                
                    vector_per_position[ix(x_index,y_index)] = 2 * (- np.power(self.c*self.dt,2 ) * (1/ np.power(self.dx,2)+1/ np.power(self.dy,2)  + 1) )

                    """

                    # references previous time value for position
                    extra_vals_per_position[ix(
                        x_index, y_index)] = -1 * (current_quantities_vector[ix(x_index, y_index)])

                    """
                    TODO:
                    create a vectors that will be passed into the sparse matrix: 
                    
                    """

                    # matrix_current_time[ix(x_index,y_index),:] = vector_per_position #######upto here!!!
                    # problem: ix(x,y) --> n,m values can be incorrect???

            # creating a sparce matrix:
            rows = cols = nm_length
            matrix_current_time = scipy.sparse.coo_matrix(
                (value, (x_values, y_values)), [rows, cols])

            known_quantities_vector = matrix_current_time*current_quantities_vector
            # all non-boundary quantities have been derived
            known_quantities_vector += extra_vals_per_position

            # fill in boundary quantitites!!
            # BC: Neumann
            for y_index_ in self.y_positions:
                for x_index_ in self.x_positions:
                    x_curr = self.x_points[x_index_]
                    y_curr = self.y_points[y_index_]

                    if x_curr in self.x_range or y_curr in self.y_range:  # checks if x,y are boundary points
                        # boundary values should be unset, ie = 0
                        val = known_quantities_vector[ix(x_index_, y_index_)]

                        if not val == 0:
                            raise ValueError(
                                f'position {x_index_,y_index_} has nonzero value, you need to set the bc values to a position with val = 0')
                        else:
                            if x_curr in self.x_range and y_curr in self.y_range:
                                # corners
                                print()
                                known_quantities_vector[ix(x_index_, y_index_)] = (self.alpha*2*self.dy + known_quantities_vector[ix(
                                    x_index_, abs(y_index_-2))] + self.alpha*2*self.dx + known_quantities_vector[ix(abs(x_index_-2), y_index_)]) * 0.5
                            elif x_curr in self.x_range and not y_curr in self.y_range:
                                # left right bc
                                known_quantities_vector[ix(
                                    x_index_, y_index_)] = self.alpha*2*self.dx + vector_per_position[ix(abs(x_index_-2), y_index_)]
                            elif y_curr in self.y_range and not x_curr in self.x_range:
                                # top bott bc
                                known_quantities_vector[ix(
                                    x_index, y_index)] = self.alpha*2*self.dy + vector_per_position[ix(x_index, abs(y_index-2))]
                            else:
                                raise ValueError(
                                    f'[{x,y}] outputs an error when setting BC')

            b = known_quantities_vector.copy()
            current_time_quantities = A/b

            matrix_all_time[i, :] = current_time_quantities

        # converting matrix into a tensor of form [time,Q(x,y)]
        nm_points = self.x_num_points * self.y_num_points

        spatial_matrix = np.array([self.x_num_points, self.y_num_points])
        temporal_spatial_matrix = np.array(
            [self.time_num_points, self.x_num_points, self.y_num_points])

        # steps through time, which is the rows, each row represents a single time step
        for index_ in range(self.time_num_points):
            current_time_step = self.time_points[index_]
            current_timestep_matrix = spatial_matrix.copy()
            current_quantities_vector = matrix_all_time[index_, :, :]
            # each time step, now fill in quantities for each (x,y)
            for index in range(nm_points):
                current_timestep_matrix[xi(
                    index)] = current_quantities_vector[index]
            # storing current_timestep_quantities matrix into all time tensor
            temporal_spatial_matrix[index_, :, :] = current_timestep_matrix

        # print the quantities: Q(x,y) over time
        X, Y = np.meshgrid(self.x_points, self.y_points)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        while True:
            for time_index in range(self.time_num_points):
                Z = temporal_spatial_matrix[time_index, :, :]
                surf = ax.plot_surface(X, Y, Z, cmap='viridis')
                plt.pause(self.delay_time)
                surf.remove()

    def run_sim(self):
        self.solve_using_pde_fd()

        # self.solve_matrix_using_vector()

    def solve_using_pde_fd(self):
        Z_initial = self.initiliaze_field()
        self.row_matrix_pde_solve(Z_initial)

    def ix(self, x, y):
        return (y*self.x_num_points) + x

    def xi(self, ix):
        y = ix // self.y_num_points
        x = ix % self.y_num_points
        return x, y

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
    start_point = [20, 20]
    wave_sim = wave(start_point)
