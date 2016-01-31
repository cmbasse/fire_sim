import numpy as np

# global temperature array
# holds the current temperature for each of the grid cells
TEMP_ARR = np.array([])

# global fuel array
# holds the current fuel amount for each of the grid cells
FUEL_ARR = np.array([])


# R - combustion function
# T - float - temperature in the cell
# V - float - amount of fuel in the cell
# returns - int - will return 1 if there is combustion, 0 if not
def R(T, V):
    if T < 572:
        return 0
    if V <= 0:
        return 0
    else:
        return 1


# dTdt - performs the differential equation for temperature change in a cell
# arr - np.array - 2D array of temperatures
# kappa - float - thermal diffusivity
# k - float - cooling constant
# delta - float - combustion constant
# h - float - height/width of the cell, set to 1 by default
def dTdt(temp_arr, fuel_arr, kappa, k, delta, h=1):
    # pad the outside of the array with
    temp_arr = np.pad(temp_arr, (1, 1), 'edge')
    # fuel_arr = np.pad(temp_arr, (1, 1), 'constant')
    new_arr = []
    rows = range(1, len(temp_arr) - 1)
    cols = range(1, len(temp_arr[0]) - 1)
    T = temp_arr[1][1]
    for i in rows:
        new_row = []
        for j in cols:
            T = (kappa / h ** 2) * (
            temp_arr[i + 1][j] - 4 * temp_arr[i][j] + temp_arr[i - 1][j] + temp_arr[i][j + 1] + temp_arr[i][
                j - 1])  # Spreading
            T -= k * temp_arr[i][j]  # Cooling
            T += delta * R(temp_arr[i][j], fuel_arr[i-1][j-1])  # Combustion
            new_row.append(T)
        new_arr.append(new_row)
    return np.array(new_arr)


# dVdt - differential equation for fuel (V) in a cell
# arr - np.array - array of fuel values
# beta - float - constant
def dVdt(temp_arr, fuel_arr, beta):
    new_arr = []
    for i in range(len(fuel_arr)):
        new_row = []
        for j in range(len(fuel_arr[0])):
            new_row.append(-beta * R(temp_arr[i][j], fuel_arr[i][j]))
        new_arr.append(new_row)
    return np.array(new_arr)


# calcTemp - calculates the current temperature using the
# previous temperature array and the derivative
# arr - np.array - 2D temperature array
# timestep - float - time interval between steps
# kappa - float - thermal diffusivity
def calcTemp(temp_arr, fuel_arr, timestep, kappa, k, delta, der_func=dTdt, h=1):
    return temp_arr + timestep * dTdt(temp_arr, fuel_arr, kappa, k, delta)


# calcFuel - calculates the current amount of fuel in a cell using
# the previous temperature array and the derivative
def calcFuel(temp_arr, fuel_arr, timestep, beta):
    return fuel_arr + timestep * dVdt(temp_arr, fuel_arr, beta)

