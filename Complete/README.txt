This code is an attempt at a solution for solving the the Lorentz system numerically using Runge-Kutta-4, which exhibits chaotic solutions when the initial conditions are changed slightly.

In this implementation, the user is asked to input the initial (x,y,z)coordinates via doubles, as well as the key parameters sigma, beta, and rho which characterize the attractor-like system.

It is important to know that when inputing the key parameters and initial coordinates, that one must input them as doubles as follows:

input x: 1.0
input y: 1.0
input z: 1.0
input dt: .001
input time: 100
input Sigma: 10
input Rho: 28
input Beta: 2.666

It should be noted that 2.666 is not equivalent to 8/3.


The easiest way to plot this program is to use gnuplot, using the splot feature. A quick and easy way to plot "LorenzHistory.txt" would be using the following code:

gnuplot> splot "LorenzHistory.txt" using 2:3:4
