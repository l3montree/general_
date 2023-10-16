import numpy as np
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
import csv

class wave:

    def __init__(self,start_point):
        self.x_start = start_point[0]
        self.y_start = start_point[1]

        self.wave_origin = start_point

        #wave terms
        self.A = None
        self.c = None #wave pde

        #spatial terms
        self.dx = 0.5
        self.dy = self.dx
        

        self.x_range =[0,100]
        self.y_range =[0,100]
      
        self.x_num_points = int((self.x_range[1]-self.x_range[0])/self.dx)+1
        self.y_num_points = int((self.y_range[1]-self.y_range[0])/self.dy)+1

        self.x_points  = np.linspace(self.x_range[0],self.x_range[1],self.x_num_points)
        self.y_points = np.linspace(self.y_range[0],self.y_range[1],self.y_num_points)

        self.x_positions = range(self.x_num_points)
        self.y_positions = range(self.y_num_points)

        #temporal terms
        self.dt = 0.5
        self.time_end = 10

        self.time_points =  [t/10 for t in range(0,self.time_end*10,self.dt*10)]

        self.wave_length = 20
        self.func = lambda A,x,y: A*(np.sin(x*np.pi/self.wave_length)+np.sin(y*np.pi/self.wave_length)) 

        #deugging
        self.csv_iter = 1

        #sim 
        self.run_sim()
        
    def initialise(self):
        self.growth_rate = 1.1
        self.curr_wave_centre = 0
        pass
    
    def after_initialise(self,time):
        tensor_all_time =[]
        matrix_current_time =[]
        vector_per_position = np.zeros(self.x_num_points*self.y_num_points)

        ix = lambda x,y: self.ix(x,y)
        xi = lambda ix_: self.xi(ix_)


        for time in self.time_points:
            for y_index in self.y_positions:
                for x_index in self.x_positions:
                    x = self.x_points[x_index]
                    y = self.y_points[y_index]

                    if x in self.x_range:
                        continue #left/right BC
                    if y in self.y_range:
                        continue #top/bott BC
                    if time == 0:
                        #time approx requires the previous time step, would not exist at t = 0
                        continue #IC
                    
                    vector_per_position[ix(x+1,y)] = 1/ self.dx^2
                    vector_per_position[ix(x-1,y)] = 1/ self.dx^2
                    vector_per_position[ix(x,y+1)] = 1/ self.dy^2
                    vector_per_position[ix(x,y-1)] = 1/ self.dy^2

                    vector_per_position[ix(x,y)] = 2*(1/(self.c*self.dt)^2 - 1/self.dx^2 - 1/self.dy^2)
                    


    
    def run_sim(self):
        self.solve_matrix_using_vector()

    def ix(self, x, y):
        return (y*self.n) + x
    
    def xi(self, ix):
        y = ix // self.m
        x = ix % self.m
        return x,y
       

    def solve_matrix_using_vector(self):

        x_ = lambda x: x - self.x_start
        y_= lambda y: y - self.y_start

        z = np.zeros(self.m*self.n) 
        self.A =1

        for y_index in self.y_positions:
            for x_index in self.x_positions:
                #this is basically how you set quantities to each spatial value
                index = int(self.ix(x_index,y_index))
                z[index] =  self.func(self.A,x_(self.x_points[x_index]),y_(self.y_points[y_index]))
        
        

        def get_solution_matrix(z_vector):
            z = np.zeros([self.m,self.n])
            for index in range(len(z_vector)):
                z[self.xi(index)] = z_vector[index]
            
            return z
        
        Z = get_solution_matrix(z)
        self.print_solution(Z)


    def print_solution(self,Z):

        X, Y = np.meshgrid(self.x_points, self.y_points)
        # Create a 3D plot
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Plot the (x, y, z) data points
        #ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
        ax.plot([self.x_start, self.x_start] ,[self.y_start,self.y_start],[-1.5,1.5], c = "k", marker ="*")
        #ax.scatter(X, Y, Z, c = Z, cmap = 'jet')

        ax.plot_surface(X, Y, Z, cmap = 'jet')

        

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
    start_point = [20,20]
    wave_sim = wave(start_point)
        

                    
    

        
    
    

        











      







