#! /usr/bin/env python
"""This program is written to solve the Schrodinger equation for a potential 
energy of V(x) = x^2, V(x) = |x|, and V(x) = 1 if |x| >1 and V(x) = 0 if |x| < 1
 and initial conditions that depend on the systems parity. We will solve this 
scenario for a symetric potential using the 4th order rung kutta method. We will
now be imposing the boundary condition that phi must aproach zero for large 
values of x. This will be done by using a root finding technique to find the 
value for E which has a zero at values of x = inf. The technique chosen for 
this is the seacant method, we can bound the energy eigenvalue to make the root 
finding easier, and use the secant method to find the energy eigenvalue to
within our desired accuracy. Finally we will be finding the expectation values 
of the x^2 quantity and also the p^2 as well. These typically require integrals
from -inf to inf, but can be approximated for sufficiently large values of x that
wave function approaches 0.

The Schrodinger Equation used has units of 1= h^2/m (as in h and 
m are equal to 1).  


Aodhan Sweeney
PHZ4151C
April 22, 2019
Schrodinger's Eqn Part B."""

#Import statements
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt

#Below is a list of functions that will be used in our program
def rung_kutta4(f, r, t, h, E, v_selector):
    """The 4th order rung kutta method for solving differential equations.
    Parameters:
    f = function (schrodinger equation in our case)
    r = vectorized function holding the values we want to solve for
    t = thing we are iterating over, typically time but position in our case
    h = step size 
    E = energy used for finding the schrodinger equation
    v_selector = dictates the potential for the schrodinger eqn
    """
    k1 = h*f(r,t, E, v_selector)
    k2 = h*f(r+0.5*k1,t+0.5*h, E, v_selector)
    k3 = h*f(r+0.5*k2,t+0.5*h, E, v_selector)
    k4 = h*f(r+k3,t+h, E, v_selector)
    return (k1 + 2*k2 + 2*k3 + k4)/6

def V(x, v_selector):
    """This function returns a postion dependent potential which is symetric and
    equivalent to the position squared.
    Parameters:
    x = position
    v_selector = dictates potential
    """
    if v_selector == 1:
        #First potential aka x**2
        potential = x**2
    elif v_selector == 2:
        #second potential aka |x|
        potential = np.abs(x)
    elif v_selector == 3:
        #third potential aka step function at |x|=1
        #V(x) = 1 if |x| > 1
        #V(x) = 0 if |x| < 1
        if np.abs(x) <= 1:
            potential = 0
        elif np.abs(x) > 1:
            potential = 1
    return(potential)

def schrodinger(vect, x, E, v_selector):
    """This function sets up the system of first order differential equations
    obtained by decomposing the second order differential eqution that is the 
    schrodinger equation. Vect is a vector which will be used to solve the 
    system of differential equations simultaneously. The x is the position, and
    both the wave function and potential are dependent on it. Finally, E is the
    energy associated with the system, and will be input by the user.
    Parameters:
    vect = vectorized function we wanted
    """
    psi = vect[0]
    phi = vect[1]
    d_psi = phi
    d_phi = 2*(V(x, v_selector) - E)*psi
    return np.array([d_psi, d_phi], float)

def appender(xpoints, E, v_selector, n):
    """appender() is the function that actually solves the differential equations.
    it does this by using the vectorization properties in python. Appender will
    compile all the data that about the solutions to the differential equation 
    into one large list. Then, the psi, phi and position values can be sifted
    out by selecting the appropriate column in the final list.
    Parameters:
    xpoints = array of positions used for finding function value
    E = energy that will be used to solve the schrodinger equation
    v_selector = specifies which potential we are solving with
    n = principle quantum number, decides parity
    """
    psi_values = []
    phi_values = []
    x_values = []
    
    #The initial conditions must be set based on the parity of the system
    if n % 2 == 0:
        vect = np.array([1.0, 0.0])
    elif n % 2 == 1:
        vect = np.array([0.0, 1.0])
    
    for x in xpoints:
        x_values += [x]
        psi_values += [vect[0]]
        phi_values += [vect[1]]
        vect += rung_kutta4(schrodinger, vect, x, xstep, E, v_selector)
    return([psi_values, phi_values, x_values])

def secant(lower_e, upper_e, xpoints, v_selector, n):
    """ secant() is the method used to find the energy eigenvalue. Secant
    root finding technique works by using the slope of the function around
    given points and incrementing the bounds to tighter and tigher regions.
    Because of this, we always check to see if the difference in the upper
    and lower energies are greater than the desired accuracy. The value 
    we are incrementing is the energy, and the function value is for
    our purposes the function value at infinity (or large x in our case).
    Parameters:
    lower_e = lower bound of the energy eigenvalue which we are looking for
    upper_e = upper bound of energy eigenvalue
    xpoints = the positions we will be calculating the psi values for, this 
    needs to be done every time we increment the energy
    v_selector = tells us which potential to calculate the psi value for
    n = principle quantum number, most importantly decided the parity
    """
    accuracy = 1e-8 #Arbitrarily decided accuracy to find the eigenvalue
    while np.abs(lower_e - upper_e) > accuracy:
        values_for_e_low = appender(xpoints, lower_e, v_selector, n)
        values_for_e_high = appender(xpoints, upper_e, v_selector, n)
        
        psi_low = values_for_e_low[0]
        psi_high = values_for_e_high[0]
        
        last_psi_low = psi_low[-1]
        last_psi_high = psi_high[-1]
        
        #The real secant method is below
        slope = (last_psi_high - last_psi_low)/(upper_e - lower_e)
        E0 = upper_e - last_psi_high/slope
        lower_e = upper_e
        upper_e = E0
    
    return(E0)


#End of functions




#BEGIN PROGRAM
if __name__ == '__main__':
    #Set the low and high boundaries for the wavefunction and variables for rung kutta solver
    steps = 100
    x_low = 0 
    x_high = 4.0 
    xstep = (x_high - x_low)/steps
    xpoints = np.arange(x_low, x_high, xstep) #array of positions to be used for finding function values
    
    #n is the principle quantum number, zero=ground state, also finding 3 excited states
    n = 0
    #Energies are going to be both set initally to zero
    energies = [0.0, 0.0]
    #Potential must be chosen and specified to solve schrodinger eqn
    v_selector = int(raw_input("Choose the potential for the schrodinger eqn: \n type 1 for x^2 \n type 2 for |x| \n type 3 for a step potential at x=1: \n Potential:  "))
    
    while n < 4: #iterating through the quantum numbers
        print('\n')
        print("EIGEN FUNCTION FOR QUANTUM NUMBER: n = {}".format(n))
        #find the range of energy eigenvalue based on quantum number    
        last_psi = [0.0, 0.0] #the final psi values will both initally be set to zero
        while last_psi[0]*last_psi[1] >= 0:
            """This while loop will find our bracketed energies for whatever quantum number we are looking for.
            This can be done by finding the two energies that produce opposite signed psi values for the last
            entry in the psi array obtained from the appender function."""
            energies[0] = energies[1]
            energies[1] += .08 #Arbitrarily choosen to increment the energy range upper limit

            values_for_e_low = appender(xpoints, energies[0], v_selector, n)
            values_for_e_high = appender(xpoints, energies[1], v_selector, n)

            psi_values_low = values_for_e_low[0]
            psi_values_high = values_for_e_high[0]

            last_psi[0] = psi_values_low[-1]
            last_psi[1] = psi_values_high[-1]

        #find real energy eigenvalue using secant method
        eigen_value = secant(energies[0], energies[1], xpoints, v_selector, n)
        print("Energy eigen value is: ", eigen_value)
        #energy eigen value is found, get psi values for this specific value
        eigen_value_data = appender(xpoints, eigen_value, v_selector, n)
        psi_values = eigen_value_data[0]
        x_values = eigen_value_data[2]


        """So far we have only calculated the values for positive values of x.
        We must now find values for negative x values based on the parity of 
        the system. Parity will be based on the principle quantum number. If 
        we have even parity, psi(-x) = psi(x). If the parity is odd then
        psi(-x) = -psi(x)."""
        if n % 2 == 0:
            psi_values = psi_values[::-1] + psi_values #Psi function is even so f(-x) = f(x)
        elif n % 2 == 1:
            negative_psi_values = []
            for i in psi_values[::-1]:
                negative_psi = -1*i
                negative_psi_values.append(negative_psi)
            psi_values = negative_psi_values + psi_values
            
        negative_x = []
        for i in x_values[::-1]: 
            #We also need to make negative all the values for the positions
            #This is done by flipping the position array and multiplying by -1
            negative_value = -1*i 
            negative_x.append(negative_value)
        x_values = negative_x + x_values


        """Next we must find the normalization factor and renormalize the wave
        function. The normalization of the wave function comes form having all 
        probabilities at any x position add up to one. To get the probability 
        amplitude we square the wave function at a given x value. Instead of 
        integrating this wavefunction and setting the integral equal to one,
        we can just multiply each probablity amplitude by the width of each
        step between the values for witch we calculate the wave function. 
        This can then be used to find the normalized wave function. """

        norm_factor = np.dot(psi_values, psi_values)*xstep
        print('Normalization factor is: ', norm_factor)
        norm_wavefunction = psi_values/(np.sqrt(norm_factor))

        """Expectation values are found by taking taking the inner product of
        both wavefunctions but sandwhiching the value for which we need in the
        middle of the inner product. This gives up the expected value (or 
        average) of the value over many trials. We will be finding
        the expectation value for both x^2 and p^2."""

        x_values = np.array(x_values)
        print('<x^2>: ', np.dot(norm_wavefunction, x_values*x_values*norm_wavefunction)*xstep)
        print('<psi|psi>', np.dot(norm_wavefunction, norm_wavefunction)*xstep)

        # using Schrodinger's Equation:
        #p**2 * Psi = -d^2Psi/dX^2 = 2 *(E - V(x))*psi
        ppPsi = np.zeros(len(norm_wavefunction) - 2)
        #Trim sizes of psi_values and x_values to match that of ppPsi
        psi_list = norm_wavefunction[1:-1]
        x_list = x_values[1:-1]
        for i in range(len(x_list)):
            ppPsi[i] = 2.0*(eigen_value - V(x_list[i], v_selector))*psi_list[i]
        print("<p^2>:", np.dot(psi_list, ppPsi)*xstep)
        


        #Finally, we must increment the principle quantum number which we are using
        #and also plot the normalized psi values
        plt.plot(x_values, norm_wavefunction, label='Energy for n = {}'.format(n))
        n +=1
    
    #The rest of the code is dedicated to creating plots for specific potentials
    if v_selector == 1:
        plt.title("Wavefunction for x^2 potential")
        plt.xlabel("Position")
        plt.ylabel("Psi(x)")
        plt.legend()
        plt.savefig("wavefunction_x_squared.png")
        plt.show()
    elif v_selector == 2:
        plt.title("Wavefunction for |x| potential")
        plt.xlabel("Position")
        plt.ylabel("Psi(x)")
        plt.legend()
        plt.savefig("wavefunction_abs_x.png")
        plt.show()
    elif v_selector == 3:
        plt.title("Wavefunction for step potential")
        plt.xlabel("Position")
        plt.ylabel("Psi(x)")
        plt.legend()
        plt.savefig("wavefunction_step_potential.png")
        plt.show()

