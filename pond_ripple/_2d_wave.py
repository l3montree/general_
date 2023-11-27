import numpy as np
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
import csv

import inspect

import time

import scipy

from matplotlib.animation import FuncAnimation

import math


class wave:

    def __init__(self, start_point):
    
        self.wave_origin = start_point

        # spatial terms
        self.dx = 1
        self.dy = self.dx

        self.x_range = [0,50]
        self.y_range = self.x_range

        self.x_start = start_point[0]
        self.y_start = start_point[1]

        if not (self.x_range[1] - self.x_start)< (self.x_range[1]- self.x_range[0]) or not (self.y_range[1] - self.y_start)< (self.y_range[1]- self.y_range[0]):
            raise ValueError(f'({self.x_start,self.y_start}) not between xrange: ({self.x_range}) and/or yrange: ({self.y_range})')
        
        self.x_points = np.arange(self.x_range[0], self.x_range[1] + self.dx, self.dx, dtype = int)
            
        self.y_points = np.arange(
            self.y_range[0], self.y_range[1] + self.dy, self.dy, dtype = int)

        self.x_num_points = len(self.x_points)
        self.y_num_points = len(self.y_points)

        # temporal terms
        self.dt = 1
        self.time_end = self.dt*5

        self.time_points = np.arange(0,self.time_end+self.dt,self.dt)
        self.time_num_points = len(self.time_points)

        # wave pde params
        self.c = 10 #velocity of the waves

        # bc params
        self.alpha = 1  # neumann

        # deugging
        self.csv_iter = 1

        # printing solution params
        self.delay_time = self.dt

        # initialise params
        self.initialised_field = False

        # sim
        self.solve_using_pde_fd()


    def initiliaze_field(self):
        field = np.zeros([self.y_num_points, self.x_num_points])
        #####
        """
        fix initial field
            how to create a convex wave with radius = drop radius, suppose to represent a drop of water on the pond
        """ 
        #####

        # you want a convex wave
        A = 1
        radius = 5
        x0, y0 = self.x_start, self.y_start
        wave_func = lambda x, y: 0.5*A*(np.sin( (x-x0)/(4*radius)*2*np.pi + np.pi/2 ) + np.sin((y-y0)/(4*radius)*2*np.pi + np.pi/2)) 
        """
        wave eqn
        https://www.wolframalpha.com/input?i=3d+plot+sin%28x%2F%284*r%29*2*pi%2Bpi%2F2%29+%2B+sin%28y%2F%28r*4%29*2*pi%2Bpi%2F2%29+for+r+%3D+5
        """
        
        i_centre = (self.x_start//self.dx)
        j_centre = (self.y_start//self.dy)

        i_var = radius//self.dx*1.1
        j_var = radius//self.dy*1.1

        for i in range(int(i_centre - i_var), int(i_centre + i_var)):
            for j in range(int(j_centre - j_var), int(j_centre + j_var)):
                if i<0 or j<0:
                    raise ValueError(f'i,j = {i,j}, which is wrong!!!')
                x = self.x_points[i]
                y = self.y_points[j]
                #print(f'x,y = {x,y}')

                if (np.power(x-self.x_start,2)+ np.power(y-self.y_start,2) > np.power(radius,2)):
                    continue

                if wave_func(x,y)>0:
                    #print(f'i,j = {i,j} x,y = {x,y} , wave_func = {wave_func(x,y)}')
                    field[j][i] = 1#wave_func(x,y)

        self.initialised_field = True

        #self.print_Z(field)

        """
         with open('initialised_data.csv', 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerows(field)
        """
        return field
    
    def convert_matrix_to_vector(self,matrix:np.ndarray):
        y_num,x_num = matrix.shape
        nm_len = int(x_num*y_num)
        vector = np.zeros(nm_len)
        ix = lambda i,j: self.ix(i,j)

        for i in range(x_num):
            for j in range(y_num):
                vector[ix(i, j)] = matrix[j][i]
        return vector
    
    def convert_vector_to_matrix(self,vector:np.ndarray):
        xi = lambda ix: self.xi(ix)

        matrix = np.zeros([self.y_num_points,self.x_num_points], dtype = np.float32)

        for ix in range(len(vector)):
            i,j = xi(ix)
            matrix[j][i] = vector[ix]
        return matrix

    def explicit_FTCS_h4(self,initial_field_matrix):

        #matrix -  vector conversions
        ix = lambda i,j: self.ix(i,j)
        xi = lambda ix: self.xi(ix)

        nm_len = int(self.x_num_points*self.y_num_points)

        #initilise matrixs/vectors
        all_time_tensor = np.zeros([self.time_num_points,self.y_num_points, self.x_num_points],dtype=np.float32)
        vectors_per_pos_matrix = np.zeros([nm_len,nm_len], dtype=np.float32)
        Q_matrix_per_time = np.zeros([self.y_num_points,self.x_num_points], dtype = np.float32)

        #initial field vector
        initial_field_vector = self.convert_matrix_to_vector(initial_field_matrix)

        #solving
        for t in self.time_num_points:
            time_ = self.time_points[t]

            A = vectors_per_pos_matrix.copy()

            #constants in the finite difference
            tau = (self.c * self.dt)**2
            alpha = 12*self.dx**2
            beta = 12*self.dy**2

            #initialises prev_quanities
            if t == 0:
                cur_quantities_vector = initial_field_vector.copy()
                prev_quantities_vector = np.zeros(nm_len, dtype = np.float32)
            else:
                prev_quantities_vector = cur_quantities_vector.copy()
                cur_quantities_vector = new_time_quantities_vector.copy()

            for i in range(self.x_num_points):
                for j in range(self.y_num_points):
                    cur_row = ix(i,j) #cur_row in A
                    
                    #BC --> dirichelet
                    if i in [0,self.x_num_points-1] or j in [0, self.y_num_points-1]:
                        A[cur_row,ix(i,j)] = 0
                    
                    else:
                        #non BC points

                        A[cur_row,cur_row] = (-30/alpha -30/beta + 2/tau)*tau
                        
                        for i in range(i-2,i+2):
                            if i> 0:
                                pass
                                
                            elif i>self.x_num_points -1:
                                continue
                            for j in range(j-2,j+2):
                                if j<0:
                                    continue
                                elif j>self.y_num_points - 1:
                                    continue
                                
                                #wrong logic!!!
                                A[cur_row,ix(i+2,j)] = -tau/ alpha
                                A[cur_row,ix(i+1,j)] = 16*tau/ alpha
                                A[cur_row,ix(i-1,j)] = 16*tau/ alpha
                                A[cur_row,ix(i-2,j)] = -tau/ alpha

                                A[cur_row,ix(i,j+2)] = -tau/ beta
                                A[cur_row,ix(i,j+1)] = 16*tau/ beta
                                A[cur_row,ix(i,j-1)] = 16*tau/ beta
                                A[cur_row,ix(i,j-2)] = -tau/ beta
                    
            


                

            new_time_quantities_vector = A*cur_quantities_vector - prev_quantities_vector




        

        
    def explicit_FTCS_h2(self, initial_field_matrix: np.ndarray):
        
        """
        solves FTCS with spatial accuracy of  O(h^2)
        """

        nm_length = self.x_num_points*self.y_num_points

        vector_per_position = np.zeros(nm_length, dtype= np.float32)

        extra_vals_per_position = vector_per_position.copy()

        ix = lambda x, y: self.ix(x, y)
        xi = lambda ix: self.xi(ix)

        #methods
        def fill_all_position_matrix(cur_pos_vector,i,j):
            #fills in "all position" matrix
            extra_vals_per_position[ix(
                    i, j)] = -1 * (prev_quantities_vector[ix(i, j)])

            all_position_coeff_matrix[ix(i,j)][:] = cur_pos_vector


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
                if not self.initialised_field: #ICs
                    raise ValueError(f'Initial Field not passed')

                prev_quantities_vector = self.convert_matrix_to_vector(initial_field_matrix)
                quantities_all_time_matrix[ti][:] = prev_quantities_vector #needs to updated for t>0
                quantities_all_time_tensor[ti][:][:] = self.convert_vector_to_matrix(prev_quantities_vector)

                #plot vals
                Z = initial_field_matrix
                ax = fig.add_subplot(111, projection='3d')

                # Plot the (x, y, z) data points
                # ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
                ax.plot([self.x_start, self.x_start], [self.y_start,
                        self.y_start], [-1.5, 1.5], c="k", marker="*")
                ax.plot_surface(X, Y, Z, cmap='jet')

                # Add labels
                ax.set_xlabel('X Label')
                ax.set_ylabel('Y Label')
                ax.set_zlabel('Z Label')
                ax.set_title(f'time = {time}')
                plt.show()
                
                continue
            
            all_position_coeff_matrix = np.zeros([nm_length,nm_length], dtype=np.float32)
            
            # t> 0
            for j in range(self.y_num_points):
                for i in range(self.x_num_points):
                    x = self.x_points[i]
                    y = self.y_points[j]

                    curr_position_quantities_vector = np.zeros(nm_length, dtype = np.float32)
                    
                    #all BCs 
                    if x in self.x_range and y in self.y_range:
                        #corners
                        curr_position_quantities_vector[ix(i, j)] = (self.alpha*2*self.dy + prev_quantities_vector[ix(
                                        i, abs(j-2))] + self.alpha*2*self.dx + prev_quantities_vector[ix(abs(i-2), j)]) * 0.5
                        fill_all_position_matrix(curr_position_quantities_vector,i,j)
                        continue
                    elif x in self.x_range:
                        #left right
                        curr_position_quantities_vector[ix(
                                    i, j)] = self.alpha*2*self.dx + prev_quantities_vector[ix(abs(i-2), j)]
                        fill_all_position_matrix(curr_position_quantities_vector,i,j)
                        continue
                    elif y in self.y_range:
                        #top bott
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
                    alpha = (self.c*self.dt)**2
                    # right point
                    curr_position_quantities_vector[ix(i+1,j)] = alpha/self.dx**2
                    
                    # left point
                    curr_position_quantities_vector[ix(i-1,j)] = alpha/self.dx**2

                    # top point
                    curr_position_quantities_vector[ix(i,j+1)] = alpha/self.dy**2

                    # bott point
                    curr_position_quantities_vector[ix(i,j-1)] = alpha/self.dy**2

                    # centre point
                    curr_position_quantities_vector[ix(i,j)] = 2*(1 - alpha/(self.dx**2)- alpha/(self.dy**2))

                    # references previous time value for position

                    if ti == 1:
                        extra_vals_per_position[ix(
                            i, j)] = -1 * (prev_quantities_vector[ix(i, j)])

                    
                    else:
                        extra_vals_per_position[ix(
                            i, j)] = -1 * (quantities_all_time_matrix[ti-1][ix(i, j)])
                        
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
            quantities_all_time_matrix[ti][:] = prev_quantities_vector

            quantities_all_time_tensor[ti][:][:] = self.convert_vector_to_matrix(prev_quantities_vector)

            #plot vals
            plt.cla()
            ax = fig.add_subplot(111, projection='3d')

            # Plot the (x, y, z) data points
            # ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
            ax.plot([self.x_start, self.x_start], [self.y_start,
                    self.y_start], [-1.5, 1.5], c="k", marker="*")

            # Add labels
            ax.set_xlabel('X Label')
            ax.set_ylabel('Y Label')
            ax.set_zlabel('Z Label')
            ax.set_title(f'time = {time}')
            ax.plot_surface(X, Y, self.convert_vector_to_matrix(prev_quantities_vector), cmap='jet')
    
        self.wave_tensor = quantities_all_time_tensor
    
    def print_all_time_tensor(self, delay_time = 3):
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

        ax.view_init(elev = 90, azim = 0 )

        if not delay_time:
            delay_time = 3

        iter = 0

        while True:
            plt.pause(delay_time)
            plt.clf()

            iter+=1
            iter = iter % self.time_num_points
            if iter == 0:#  self.time_num_points:
                break
            
            fig = plt.figure(figsize = (10,10))
            ax = fig.add_subplot(111, projection='3d')
            Z = quantities_all_time_tensor[0, :, :]
            ax.plot(self.x_start,self.y_start,0,color = "k", marker= "*", markersize = 30)
            ax.plot_surface(X, Y, Z, cmap='viridis')
            ax.set_xlabel('X-axis')
            ax.set_ylabel('Y-axis')
            ax.set_zlabel('Z-axis')
            ax.set_title(f'3D Surface Plot (t = {self.time_points[0]})')
            #ax.plot(self.x_start,self.y_start,0,color = "k", marker= "*", markersize = 30)

            Z = np.clip(quantities_all_time_tensor[iter, :, :],-2,2)
            #print(f'time = {self.time_points[iter]}')
            #print(Z)
            ax.plot_surface(X, Y, Z, cmap='viridis')
            
            ax.set_xlabel('X-axis')
            ax.set_ylabel('Y-axis')
            ax.set_zlabel('Z-axis')
            ax.set_title(f'3D Surface Plot (t = {self.time_points[iter]})')
            ax.set_zlim(-6,6)


            

    def solve_using_pde_fd(self):
        #row major numbering
        Z_initial = self.initiliaze_field()
        self.explicit_FTCS(Z_initial)
        self.print_all_time_tensor()

    def ix(self, i, j):
        return (j*(self.x_num_points)) + i 

    def xi(self, ix):
        #remeber its row major numbering !!
        j = ix // (self.x_num_points)
        i = ix % (self.x_num_points)
        return i,j


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


if __name__ == "__main__":
    start_point = [10, 10]
    wave_sim = wave(start_point)
