![](https://komarev.com/ghpvc/?username=l3montree)
# General
General contains all current projects that I am working on.
Projects or topics that I am particularly interested in are:
- PDE approximations through Finite Element and/or Finite Diference
- CFD simulations, currently focused on Ansys Fluent, however have plans to explore openFoam

Current Projects in this folder:
## Pond Ripple: 2D wave
Finite Difference Approximation of  the 2D wave equation using:
- BC: Dirichlet
- IC: 3 random points in the domain have Q(x,y,t =0) = 1, elsewhere Q = 0
- Has overall error: O(h^2) or O(h^4) depending on spatial scheme

Extra work: 
- Finding dt through brute force, then saving value, using data you can interpolate to other points or maybe use Gauss Seidel method to determine dt for other dx,dy values

## FEM
Some work into understanding Finite Element Method

## Eularian Fluid Sim
2D Finite Difference Appromixation using collated cells of inviscid, incompressible ideal fluids

## Visa Free Fiji
Shows the countrys in which fiji citizens gain visa free entry. Achieved through web scapping a website

