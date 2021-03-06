import os
os.environ['SDL_VIDEO_WINDOW_POS'] = "%d,%d" % (100,25)

import pygame
from pygame.locals import *
import fire
import numpy as np
from random import randint
from numpy.random import uniform, randn
import matrixes as mtx

TILE_W = 17
TILE_H = 23


class Map:
    """stores map info, and draws tiles.
    Map is stored as an array of int's which correspond to the temperature values."""

    def __init__(self, screen, fire, temp_array, fuel_array):
        """initialize default values  """
        # load our 3 images
        self.ALIVE_TREE = self.load_image('alive_tree.jpg')
        self.FIRE_TREE = self.load_image('fire_tree.jpg')
        self.BURNT_TREE = self.load_image('burnt_tree.jpg')
        # initialize class variables
        self.difference_array = np.zeros((len(temp_array), len(temp_array)))
        self.SCREEN = screen
        self.FIRE = fire
        self.temp_array = temp_array
        self.fuel_array = fuel_array
        self.draw()

    @staticmethod
    def load_image(image):
        """load image"""
        return pygame.image.load(os.path.join('img', image))

    def update(self, temp_array, fuel_array):
        """update the map"""

        temp = []
        for i in range(len(temp_array)):
            temp_row = []
            for j in range(len(temp_array)):
                # if the location has all ready started cooling keep its difference at 1
                if self.difference_array[i][j] == 1:
                    temp_row.append(1)
                elif self.temp_array[i][j] >= self.FIRE:
                    # if the location is cooling change difference to 1
                    if temp_array[i][j] < self.FIRE:
                        temp_row.append(1)
                    else:
                        temp_row.append(0)
                else:
                    temp_row.append(0)

            temp.append(temp_row)
        # update class variables
        self.difference_array = temp
        self.temp_array = temp_array
        self.fuel_array = fuel_array
        # redraw the screen

        self.draw()

    def draw(self):
        """Turns the 2d array of numbers into a grid of pictures"""
        # loop all tiles, and draw
        size = len(self.temp_array)
        for i in range(size):
            for j in range(size):
                x = TILE_W*j
                y = TILE_H*i
                # based on temperature choose which image to display
                if self.temp_array[i][j] < self.FIRE and self.fuel_array[i][j] > 0:
                    self.SCREEN.blit(self.ALIVE_TREE, (x, y))
                if self.fuel_array[i][j] <= 0 or self.difference_array[i][j] == 1:
                    self.SCREEN.blit(self.BURNT_TREE, (x, y))
                if self.temp_array[i][j] >= self.FIRE and self.fuel_array[i][j] > 0:
                    self.SCREEN.blit(self.FIRE_TREE, (x, y))
                pygame.display.flip()


def main():
    """The function that runs it all"""
    # Initialise screen
    pygame.init()
    w = 600
    #h = 740
    h = 810
    size = (w, h)
    screen = pygame.display.set_mode(size)
    pygame.display.set_caption('Fire')
    background = pygame.Surface(screen.get_size())
    background = background.convert()
    background.fill((250, 250, 250))
    screen.blit(background, (0, 0))
    pygame.display.flip()

    temp_arr = mtx.temp_arr2
    fuel_arr = 2*np.ones((len(temp_arr), len(temp_arr)))

    # run the first iteration of the fire simulation
    sim = Map(screen, 572, temp_arr, fuel_arr)

    # Event loop
    while 1:
        # run the simulation over time
        temp_arr = fire.calcTemp(temp_arr, fuel_arr, .05, 5, .1, 5000, 110, der_func=fire.dTdt, h=1)
        fuel_arr = fire.calcFuel(temp_arr, fuel_arr, .01, 15)
        # pygame.time.wait(10)
        sim.update(temp_arr, fuel_arr)
        for event in pygame.event.get():
            if event.type == QUIT:
                return

if __name__ == '__main__':
    main()
