import numpy as np
import matplotlib.pyplot as plt

func = lambda x: np.sin(x)
# Create example data
x = np.linspace(0, 40, 41)

start_point = (10,20)

x_ = lambda x: x - start_point[0]

plt.plot(x_(x),func(x_(x)))
plt.plot(x,0, c  = "r")
plt.plot()
plt.plot(start_point, marker = "*")

# Show the plot
plt.show()