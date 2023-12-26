import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib.animation import FuncAnimation
#import pandas as pd

from datetime import datetime
import time

from pathlib import *
import os

class Wave:

    def __init__(self, dx,dy,xRange, yRange, c):

        # wave pde params
        self._c = c  # velocity of the waves

        # spatial terms
        self._dx = dx
        self._dy = dy

        self.x_range = xRange
        self.y_range = yRange

        # x,y points
        self.x_points = np.arange(
            self.x_range[0], self.x_range[1] + self.dx, self.dx, dtype=np.float32)

        self.y_points = np.arange(
            self.y_range[0], self.y_range[1] + self.dy, self.dy, dtype=np.float32)

        self.x_num_points = len(self.x_points)
        self.y_num_points = len(self.y_points)

        #matrix solution
        self.A =None

        #plotting
        self.X, self.Y = np.meshgrid(self.x_points,self.y_points)

        #initialised field?
        self.is_initialised = False
        self.initial_field = None

        #plot animation
        self._animate_figure = True

        #local variables for iteration
        self._is_iteration = False
        self.iteration_iter = 0
        self._cur_quantities_vector = None
        self._prev_quantities_vector = None
        self._cur_time = None
        self._initial_field = None

        #temporal variables                                                                                                                                                                                                  
        self.dt = None #will be set later
        self.time_end = None
        self.time_points = None
        self.time_num_points = None

        #using CFL derived time
        self._use_CFL = True #False if using brute force derived dt

        #bruteForce_dt study variables 
        self.excel_file = "stability_values.xlsx"
        self._is_stability_study = False

        #droplets
        self._droplet_counter = 0
        self.droplets = 3

    #for viewing the plot: either animate plots or produce individual plots for each timestep
    @property
    def animate_figure(self):
        return self._clear_figure
    @animate_figure.setter
    def animate_figure(self,bool):
        self._animate_figure = bool
        
    #independent variables for CFL investigation
    @property
    def is_stability_study(self):
        return self._is_stability_study
    @is_stability_study.setter
    def CFL(self,bool):
        self._is_stability_study = bool
    
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
    

    @property
    def use_CFL(self):
        return self._use_CFL
    @use_CFL.setter
    def use_CFL(self,bool):
        self._use_CFL = bool

    def run(self):

        field = self.initialise_field()

        print(f'Field initialised. Q(x,y,t = 0)')
        print(f'FDA: {self.explicit_FTCS_h4.__name__}')

        self.explicit_FTCS_h4()
        print(f'self.run() complete')

    def initialise_field(self):
        """
        initialises the field
         - three points at random: with corresponding amplitudes: 3,2,1
        """
        field = np.zeros(
            [self.y_num_points, self.x_num_points], dtype=np.float32)

        points = [(int(self.y_range[1]//self.dy*0.2), int(self.x_range[1]//self.dx*0.7)),(int(self.y_range[1]//self.dy*0.5), int(self.x_range[1]//self.dx*0.3)),(int(self.y_range[1]//self.dy*0.85), int(self.x_range[1]//self.dx*0.85))] #points list
        for i in range(3):
            #point = [random.randint(self.x_range[:]),random.randint(self.y_range[:])] #uncomment if want to randomise points
            point = points[i]
            field[point] = i+1
    
        self.is_initialised = True
        self.initial_field = field

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

        return np.float32(vector)


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
    

    """
    ix, xi used for matrix - vector conversions using row major numbering
    """
    def ix( self,i, j):
        return (j*(self.x_num_points)) + i

    
    def xi( self,ix):
        j = ix // (self.x_num_points)
        i = ix % (self.x_num_points)
        return i, j

    def explicit_FTCS_h4(self):
        """
        Explicit Finite difference approximation using:
            - Time: centred difference O(h^2)
            - space: centred difference O(h^4)
        """
        if not self.is_initialised:
            raise ValueError(" explicit_FTCS_h4 called but field not initialised")     
           
        #time
        if self.use_CFL: #uses CFL derived dt
            if not self.dx == self.dy:
                raise ValueError(f'self.use_CFL = True, cfl calcs has equal grid spacing as a assumption, CFL codition')
            self.dt = round((1/5)**0.5*self.dx/self.c,2)

            self.time_end = 10
            self.time_points = np.arange(0,self.time_end+self.dt,self.dt)
        else:
            #implement brute force dt study

            if not self.dt:
                raise ValueError(f'self.use_CFL = False but self.set_temporal_points() = Nan')

        # matrix -  vector conversions
        def ix(i, j): return self.ix(i, j)
        def xi(ix): return self.xi(ix)

        nm_len = int(self.x_num_points*self.y_num_points)
        
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
            self.initial_field)
        
        #loop variable init
        t = 0
        
        #creating A
        """
        A contains all coefficients of the field quantities

        A x = b
            A --> N x M matrix containing N - x coords, M -y coords, depending on how big spatial domain is and how you discretise it
                contains coefficients of the quantities, coefficients depend on finite difference approximation

            x --> vector of quanitities in the field of size (NM x 1) at current time step
            b --> new quantity at new time step

        """
        if not self.A:
            start = datetime.now()
            self.A = np.zeros([nm_len, nm_len], dtype=np.float32)

            for i in range(self.x_num_points):
                for j in range(self.y_num_points):
                    cur_row = ix(i, j)  # cur_row in A

                    # BC --> dirichelet
                    if i in [0, self.x_num_points-1] or j in [0, self.y_num_points-1]:
                        self.A[cur_row, ix(i, j)] = 0

                    else:  # [(2:xend-2),(2:yend-2)]
                        # non BC points

                        self.A[cur_row, cur_row] = centre_coeff

                        for i_ in [a for a in range(i-2, i+3) if a >= 0 and a <= self.x_num_points-1]:

                            if i_ == i:
                                for j_ in [a for a in range(j-2, j+3) if a > 0 and not a == j and a <= self.y_num_points - 1]:

                                    if j_ > j:
                                        coef_index = j_ - (j - 1)
                                    else:
                                        coef_index = j_ - (j - 2)

                                    self.A[cur_row, ix(
                                        i, j_)] = Ylap_coef_list[coef_index]
                            else:
                                if i_ > i:
                                    coef_index = i_ - (i - 1)
                                else:
                                    coef_index = i_ - (i - 2)
                                self.A[cur_row, ix(i_, j)] = Xlap_coef_list[coef_index]

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
            end  = datetime.now()
            print(f'A created, time elapsed: {(end - start).total_seconds()* 10**3:.03f} ms')

        if not self._is_iteration: #active when brute force dt study is NOT being conducted
            
                if self._animate_figure:

                    fig = plt.figure()
                    ax = fig.add_subplot(111, projection = '3d')

                    ax.set_xlabel('X')
                    ax.set_ylabel('Y')
                    ax.set_zlabel('Q')
                    ax.set_zlim(-4,4)

                    Z = self.initial_field
                    ax.plot_surface(self.X, self.Y, Z, cmap='jet')
                    ax.set_title(f"2D Wave Finite Diff Approx @ t = {0}")

                    self._cur_quantities_vector = initial_field_vector.copy()
                    self._prev_quantities_vector = np.zeros(nm_len, dtype=np.float32)

                    def update_animation(frame):
                        print(f'update animation {frame}')
                        ax.cla()

                        if frame == 0:
                            pass
                            #start_ = datetime.now()
                        else:
                            # updates the vectors
                            
                            self._prev_quantities_vector = self._cur_quantities_vector.copy()
                            self._cur_quantities_vector = self.new_time_quantities_vector.copy()
                        
                        start_ = datetime.now()
                        self.new_time_quantities_vector = self.A@self._cur_quantities_vector - self._prev_quantities_vector
                        end_ = datetime.now()
                        Z = self.convert_vector_to_matrix(self.new_time_quantities_vector)

                        
                        print(f'time step calc complete: time  = {frame}, time elapsed: {(end_ - start_).total_seconds():.03f} s')
                        
                        ax.plot_surface(self.X, self.Y, Z, cmap='jet')
                        ax.set_zlim(-4,4)
                        ax.set_title(f"2D Wave Finite Diff Approx @ t = {round(frame,2)}")
                    
                    ani = FuncAnimation(fig, update_animation, frames = self.time_points, interval = 200, repeat = False)
                    plt.show()

                else:
                    print(f'Plotting individual timesteps')

                    while True:

                        if t == 0:
                            _cur_quantities_vector = initial_field_vector.copy()
                            _prev_quantities_vector = np.zeros(nm_len, dtype=np.float32)
                            
                            Z = self.initial_field
                            t += self.dt
                            print(f't = {t}')
                        else:
                            # updates the vectors
                            _prev_quantities_vector = _cur_quantities_vector.copy()
                            _cur_quantities_vector = new_time_quantities_vector.copy()
                        
                            new_time_quantities_vector = self.A@_cur_quantities_vector - _prev_quantities_vector

                            Z = self.convert_vector_to_matrix(new_time_quantities_vector)
                            print(f't = {t}')
                        
                        fig = plt.figure()
                        ax = plt.axes(projection='3d')
                        ax.plot_surface(self.X, self.Y, Z, cmap='jet')
                        ax.set_xlabel('X')
                        ax.set_ylabel('Y')
                        ax.set_zlabel('Q')
                        ax.set_title(f"2D Wave Finite Diff Approx @ t = {t}")
                        ax.set_zlim(-4,4)

                        t += self.dt

        else:
            self.iteration_iter+=1
            #print(f'iteration, iter {self.iteration_iter}')

            if self._cur_time == 0:  # initialises cur_quantities and prev_quanities
                self._cur_quantities_vector = initial_field_vector.copy()
                self._prev_quantities_vector = np.zeros(nm_len, dtype=np.float32)

                self._cur_time += self.dt

            new_time_quantities_vector = self.A@self._cur_quantities_vector - self._prev_quantities_vector

            self._prev_quantities_vector = self._cur_quantities_vector.copy()
            self._cur_quantities_vector = new_time_quantities_vector.copy() #used by brute force study

            self._cur_time += self.dt

    def explicit_FTCS_h2(self, initial_field_matrix: np.ndarray):
        """
           Explicit Finite difference approximation using:
            - Time: centred difference O(h^2)
            - space: centred difference O(h^2)
        """
               #set temporal points
        self.set_temporal_points()

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

                if self._animate_figure:
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

    def __iter__(self):
        self._is_iteration = True
        if not self.is_initialised:
            self.initialise_field()
        self.explicit_FTCS_h4()
        
        #print(f"iter = {self.iteration_iter}: dx = {self.dx}, dy = {self.dy}, dt ={self.dt}")
        return self
    
    def __next__(self):
        if not self._is_iteration:
            raise StopIteration


if __name__ == "__main__":
    dx = dy = 1
    c = 10
    xRange = yRange = [0,50]

    sim = Wave(dx,dy,xRange, yRange, c)
    sim.animate_figure = True
    sim.run()
