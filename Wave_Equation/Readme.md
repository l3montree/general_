# 2D Wave Equation:

Complete notes here: 

https://docs.google.com/document/d/1UdAre0QydwtkgXZ3MHh-wKCYS5iFsHr1RvFlGK6Nx_w/edit?usp=sharing

To run: 
    python3 Wave_2D.py [-h] [--h] [--xRange] [--yRange] [--c] [--verbose] [--to_animate]

    [-h]            help
    [--h]           increments in spatial domain
    [--xRange]      xRange of domain eg. xRange = [0,100]
    [--yRange]      yRange of domain eg. yRange = [0,100]
    [--c]           wave speed
    [--verbose]     view extra info?
    [--to_animate]  animate plot? or view individual plots?

Example usage:

    python3 Wave_2D.py --h 1 --xRange [0,100] --yRange [0,100] --c 10 --verbose True --to_animate True

    