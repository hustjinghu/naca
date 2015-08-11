#!/usr/bin/python3

""" naca.py

Created: 1/1/2015
Author: Michel Robijns

This file is part of naca which is released under the MIT license.
See the file LICENSE or go to http://opensource.org/licenses/MIT for full
license details.

Description:
    A set of functions to compute the coordinates of NACA airfoils. In this
    version, only the NACA 4-digit airfoil series has been implemented.
"""

from math import *
import numpy as np


def main():
    pass


if __name__ == '__main__':
    main()


def NACA4(number, N, half_cosine_spacing=True, closed_trailing_edge=True, save_to_file=False):
    """Computes coordinates of a NACA 4-digit airfoil [1].
    
    Arguments:
        number: Name of the requested 4-digt NACA airfoil entered as a string,
                i.e. '2412'
        N: Number of desired airfoil coordinates
    
    Optional arguments:
        half_cosine_spacing: Half cosine spacing ensures that the datapoints
                             are more widely spaced around the leading edge
                             where the curvature is greatest.
        closed_trailing_edge: The trailing edge has a small but finite
                              thickness when using equations in [1] unaltered.
        save_to_file: Saves the coordinates in a data file.
    
    Returns:
        A matrix with two columns pertaining to the x and y-coordinates,
        respectively. The sequence of coordinates is clockwise, starting at the
        trailing edge.
    
    Raises:
        ValueError: Airfoil number does not have four digits or N is negative
        
    References:
    
    .. [1] Abbott, I.H. and von Doenhoff, A.E., "Theory of Wing Sections,
       Including a summary of airfoil data", Dover Publications, pp. 111-123,
       1959.
    
    """
    
    # Check whether input parameters are valid or not
    if not (0 < int(number) < 10000 and N > 0):
        raise ValueError("Invalid input.")
    
    # Scale the input parameters. We want, for example, a thickess of 0.12
    # and not 12 and the maximum camber at 0.4 * c and not 4 * c.
    
    m = float(number[0])
    p = float(number[1])
    t = float(number[2:])
    
    if m != 0:
        m = m / 100
    if p != 0:
        p = p / 10
    if t != 0:
        t = t / 100
    
    # Half cosine spacing ensures that the datapoints are more widely spaced
    # around the leading edge where the curvature is greatest.
    if half_cosine_spacing:
        x = (1 - np.cos(np.linspace(0, pi, N, dtype=float))) / 2
    else:
        x = np.linspace(0, 1, N)
    
    # The trailing edge has a small but finite thickness by default. The gap
    # can be closed by utilizing a slightly different equation. See equation
    # 6.4 of Abbot and von Doenhoff
    if closed_trailing_edge:
        thickness = t / 0.20 * (0.29690 * np.sqrt(x) - 0.12600 * x - 0.35160 *
                    np.power(x, 2) + 0.28430 * np.power(x, 3) - 0.10360 *
                    np.power(x, 4))
    else:
        thickness = t / 0.20 * (0.29690 * np.sqrt(x) - 0.12600 * x - 0.35160 *
                    np.power(x, 2) + 0.28430 * np.power(x, 3) - 0.10150 *
                    np.power(x, 4))
    
    # Compute the y-coordinates of the camber line. See equation 6.4 of Abbot
    # and von Doenhoff
    fwd_x = x[x < p]
    aft_x = x[x >= p]
    
    if 0 < p < 1 and 0 < m < 1:
        fwd_camber = m / p**2 * (2 * p * fwd_x - np.power(fwd_x, 2))
        aft_camber = m / (1 - p)**2 * ((1 - 2 * p) + 2 * p * aft_x -
                     np.power(aft_x, 2))
        camber = np.append(fwd_camber, aft_camber)
    else:
        camber = np.zeros(np.size(x))
    
    y_upper = np.add(camber, thickness)
    y_lower = np.subtract(camber, thickness)
    
    # Append the arrays that contain the x and y-coordinates
    x_upper = np.flipud(x)
    x_lower = x[1:]
    x = np.append(x_upper, x_lower)
    
    y_upper = np.flipud(y_upper)
    y_lower = y_lower[1:]
    y = np.append(y_upper, y_lower)
    
    # NumPy creates a tiny numerical error when adding and subtracting the two
    # arrays. A number of order 10^-17 is shown instead of 0. Likely due to the
    # computer's internal representation of numbers.
    if closed_trailing_edge:
        y[0] = 0
        y[-1] = 0
    
    coordinates = np.column_stack((x, y))
    
    if save_to_file:
        np.savetxt("NACA " + number + ".dat", coordinates, delimiter='\t', fmt='%f', header="NACA " + number)
    
    return coordinates
