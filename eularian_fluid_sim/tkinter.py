import tkinter as tk

# Create a Tkinter window
window = tk.Tk()
window.title("Fluid Simulation")

# Set the window size
window_width = 500
window_height = 500
window.geometry(f"{window_width}x{window_height}")

# Ball properties
ball_radius = 50
ball_x, ball_y = window_width // 2, window_height // 2

# Function to move the ball with the mouse cursor
def move_ball(event):
    global ball_x, ball_y
    ball_x, ball_y = event.x, event.y
    canvas.coords(ball, ball_x - ball_radius, ball_y - ball_radius, ball_x + ball_radius, ball_y + ball_radius)

# Bind mouse motion event to move_ball function
window.bind("<Motion>", move_ball)

# Create a canvas to draw the ball
canvas = tk.Canvas(window, width=window_width, height=window_height, bg="white")
canvas.pack()

# Draw the ball
ball = canvas.create_oval(ball_x - ball_radius, ball_y - ball_radius, ball_x + ball_radius, ball_y + ball_radius, fill="blue")

# Start the Tkinter event loop
window.mainloop()

class ball():
    def __init__():
        pass
    def moveball():
        pass
    def 