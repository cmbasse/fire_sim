import numpy as np

# globals needed for heatmap simulations

# global temparature array
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
    if V < 0:
        return 0
    else:
        return 1


# cooling - returns the Newton cooling term for dVdt
# temp - float - current temperature
# ambientTemp - float - ambient air temperature
# k - float - Newton cooling rate
# returns - float - current cooling amount
def cooling(temp, ambientTemp, k):
    if temp > ambientTemp:
        return -k * temp
    else:
        return 0


# dTdt - performs the differential equation for temperature change in a cell
# arr - np.array - 2D array of temperatures
# kappa - float - thermal diffusivity
# k - float - cooling constant
# delta - float - combustion constant
# h - float - height/width of the cell, set to 1 by default
# returns - np.array - 2D array of derivatives
def dTdt(temp_arr, fuel_arr, kappa, k, delta, ambientTemp, h=1):
    # pad the outside of the array with
    temp_arr = np.pad(temp_arr, (1, 1), 'edge')
    new_arr = []
    rows = range(1, len(temp_arr) - 1)
    cols = range(1, len(temp_arr[0]) - 1)
    for i in rows:
        new_row = []
        for j in cols:
            T = (kappa / h ** 2) * (
            temp_arr[i + 1][j] - 4 * temp_arr[i][j] + temp_arr[i - 1][j] + temp_arr[i][j + 1] + temp_arr[i][
                j - 1])  # Thermal diffusion
            T += cooling(temp_arr[i][j], ambientTemp, k)  # Cooling
            T += delta * R(temp_arr[i][j], fuel_arr[i - 1][j - 1])  # Combustion
            new_row.append(T)
        new_arr.append(new_row)
    return np.array(new_arr)


# dVdt - differential equation for fuel (V) in a cell
# arr - np.array - array of fuel values
# beta - float - constant
# returns - np.array - 2D array of derivatives
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
# temp_arr - np.array - 2D temperature array
# fuel_arr - np.array - 2D fuel array
# timestep - float - time interval between steps
# kappa - float - thermal diffusivity
# returns - np.array - updated array of temperatures
def calcTemp(temp_arr, fuel_arr, timestep, kappa, k, delta, ambientTemp, der_func=dTdt, h=1):
    return temp_arr + timestep * dTdt(temp_arr, fuel_arr, kappa, k, delta, ambientTemp)


# calcFuel - calculates the current amount of fuel in a cell using
# the previous temperature array and the derivative
# temp_arr - np.array - 2D temperature array
# fuel_arr - np.array - 2D fuel array
# timestep - float - time interval between steps
# beta - float - constant
# returns - np.array - updated array of temperatures
def calcFuel(temp_arr, fuel_arr, timestep, beta):
    return fuel_arr + timestep * dVdt(temp_arr, fuel_arr, beta)
