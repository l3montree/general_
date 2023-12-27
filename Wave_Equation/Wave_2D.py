import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from datetime import datetime
import argparse
"""
Author: https://github.com/l3montree
"""

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
     
        self._cur_time = None

        self._initial_field_vector =None
        self._initial_field = None

        self._cur_quantities_vector = None
        self._prev_quantities_vector = None
        self._new_time_quantities_vector = None

        #temporal terms                                                                                                                                                                                                  
        self.dt = round((1/5)**0.5*self.dx/self.c,2) #CFL for a CTCS wave equation
        self.time_end = 10
        self.time_points = np.arange(0,self.time_end+self.dt,self.dt)

        #debugging 
        self._verbose = True

    @property
    def animate_figure(self):
        return self._clear_figure
    @animate_figure.setter
    def animate_figure(self,bool):
        self._animate_figure = bool
    
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
    def verbose(self):
        return self._verbose
    @c.setter
    def verbose(self,bool):
        self._verbose = bool


    def run(self):

        print(f'Running sim')

        field = self.initialise_field()

        print(f'Field initialised ie. Q(x,y,t = 0)')
        print(f'Calculating Finite Difference Approx')

        self.CTCS()

        print(f'sim complete')

    def initialise_field(self):
        """
        initialises the field
         - three points spread across the field at different amplitudes, everywhere else, Q = 0
        """
        field = np.zeros(
            [self.y_num_points, self.x_num_points], dtype=np.float32)

        points = [(int(self.y_range[1]//self.dy*0.2), int(self.x_range[1]//self.dx*0.7)),(int(self.y_range[1]//self.dy*0.5), int(self.x_range[1]//self.dx*0.3)),(int(self.y_range[1]//self.dy*0.85), int(self.x_range[1]//self.dx*0.85))] #points list
        amplitudes = [1,2,3]

        for i in range(3):
            point = points[i]
            field[point] = amplitudes[i]
    
        self.is_initialised = True
        self._initial_field = field

    def CTCS(self):
        """
        Explicit Finite difference approximation using:
            - Time: centred difference O(h^2)
            - space: centred difference O(h^4)
        """
        if not self.is_initialised:
            raise ValueError("CTCS called but field not initialised")     

        # matrix -  vector conversions
        def ix(i, j): return self.ix(i, j)

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
        self._initial_field_vector = self.convert_matrix_to_vector(
            self._initial_field)
        
        #loop variable init
        t = 0
      
        """
        Creating A:

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
            if self.verbose:
                print(f'A created, time elapsed: {(end - start).total_seconds():.03f} s')

       
        if self._animate_figure:
            print(f'Creating 2D wave animation')

            fig = plt.figure()
            ax = fig.add_subplot(111, projection = '3d')

            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Q')
            ax.set_zlim(-4,4)

            Z = self._initial_field
            ax.plot_surface(self.X, self.Y, Z, cmap='jet')
            ax.set_title(f"2D Wave Finite Diff Approx @ t = {0}")

            self._cur_quantities_vector = self._initial_field_vector.copy()
            self._prev_quantities_vector = np.zeros(nm_len, dtype=np.float32)

            def update_animation(frame):
                start_ = datetime.now()

                ax.cla()

                if not frame == 0:
                    self._prev_quantities_vector = self._cur_quantities_vector.copy()
                    self._cur_quantities_vector = self._new_time_quantities_vector.copy()
            
                self._new_time_quantities_vector = self.A@self._cur_quantities_vector - self._prev_quantities_vector
                
                Z = self.convert_vector_to_matrix(self._new_time_quantities_vector)

                end_ = datetime.now()
                
                ax.plot_surface(self.X, self.Y, Z, cmap='jet')
                ax.set_zlim(-4,4)
                ax.set_title(f"2D Wave Finite Diff Approx @ t = {round(frame,2)}")
                
                if self.verbose:
                    print(f'calc @time step complete: time  = {round(frame,2)}, time elapsed: {(end_ - start_).total_seconds():.03f} s')
            
            ani = FuncAnimation(fig, update_animation, frames = self.time_points, interval = 200, repeat = True)
            plt.show()

        else:
            print(f'Plotting individual timesteps')

            for t in self.time_points:

                start__ = datetime.now()

                if t == 0:
                    self._cur_quantities_vector = self._initial_field_vector.copy()
                    self._prev_quantities_vector = np.zeros(nm_len, dtype=np.float32)
                    
                    Z = self._initial_field
                    t += self.dt
                else:
                    # updates the vectors
                    self._prev_quantities_vector =self._cur_quantities_vector.copy()
                    self._cur_quantities_vector = self._new_time_quantities_vector.copy()
                
                self._new_time_quantities_vector = self.A@self._cur_quantities_vector - self._prev_quantities_vector

                Z = self.convert_vector_to_matrix(self._new_time_quantities_vector)
                
                end__ = datetime.now()

                if self.verbose:
                    print(f'calc @time step complete: time  = {t}, time elapsed: {(end__ - start__).total_seconds():.03f} s')

                fig = plt.figure()
                ax = plt.axes(projection='3d')
                ax.plot_surface(self.X, self.Y, Z, cmap='jet')
                ax.set_xlabel('X')
                ax.set_ylabel('Y')
                ax.set_zlabel('Q')
                ax.set_title(f"2D Wave Finite Diff Approx @ t = {round(t,2)}")
                ax.set_zlim(-4,4)

                t += self.dt

    
    #methods used in for matrix/vector conversion
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
            matrix[j][i] = vector[ix]
        return matrix
    
    def ix( self,i, j):
        #row major numbering used for matrix --> vector conversion
        return (j*(self.x_num_points)) + i
    
    def xi( self,ix):
        #row major numbering used for vector --> matrix conversion
        j = ix // (self.x_num_points)
        i = ix % (self.x_num_points)
        return i, j
    

    #class dunder methods
    def __iter__(self):

        self._is_iteration = True
        if not self.is_initialised:
            self.initialise_field()
        self.CTCS()
        
        return self
    
    def __next__(self):
        if not self._is_iteration:
            raise StopIteration


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = 'inputs for WAVE_2D.Wave()')
    parser.add_argument('--h', type = float, help ="spatial increments",default=1)
    parser.add_argument('--xRange', type = list[float], help  = 'x - Range of domain eg. xRange = [0,50]',default= [0,50])
    parser.add_argument('--yRange',type = list[float], help  = 'y - Range of domain eg. yRange = [0,50]', default =[0,50])
    parser.add_argument('--c', type = float,help = "wave speed", default=10)
    parser.add_argument('--verbose', type = bool, help ="to view extra info", default = True)
    parser.add_argument('--to_animate', type = bool, help ="view plots as animation?", default = True)

    arg = parser.parse_args()

    sim = Wave(arg.h,arg.h,arg.xRange,arg.yRange,arg.c)
    sim.verbose = arg.verbose
    sim.animate_figure = arg.to_animate

    print(f'Running Wave_2D.Wave() with params:\nh: {arg.h}, xRange: {arg.xRange}, yRange: {arg.yRange}, c: {arg.c}')
    sim.run()

    """
    Example usage

    python3 Wave_2D.py --h 1 --xRange [0,100] --yRange [0,100] --c 10 --verbose True --to_animate True
    """
