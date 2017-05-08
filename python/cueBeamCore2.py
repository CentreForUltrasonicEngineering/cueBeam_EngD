import math
import numpy
import math
import matplotlib.pyplot as plt
import time
import copy
import warnings
import dill
import redis
from celery import Celery


""" 
This is the core of the cueBeam.

in the version 2, i refactor the earlier code so that all information is stored in a flattish object instead of multi level objects
this is to make it easier to split the work into workers
    
**Essential fields that You need to set in a CueBeamWorld object:**

1. wavenumber: float
    The environment property, note that wavenumber = 1/wavelength in that enviroment
2. rxPlane: single item of class RxPlane()
    describes the regular sampling plane
    
3. elements: list of TxElement()
    description of the transmitters
"""


# # # # # CELERY SETTINGS
# storage for dynamically stored functions
RSTORE = redis.StrictRedis('192.168.0.43', port=6379, db=1)
# connection to the celery object, broker and broker's back-end
CELERY = Celery('tasks', broker='amqp://guest:guest@192.168.0.43:5672/',backend='redis://192.168.0.43:6379/0')

CELERY.conf.accept_content = ['json', 'msgpack']
CELERY.conf.result_serializer = 'msgpack'


class CueBeamWorld:
    """Contains all the info needed for cueBeam core code to return a field"""

    class TxElement:
        """Describes a single transmitting monopole."""
        x = 0.0
        """X-location of the monopole."""
        y = 0.0
        """Y-location of the monopole."""
        z = 0.0
        """Z-location of the monopole."""
        amplitude = 1.0
        """Amplitude of the vibration of the monopole."""
        phase = 0.0
        """Relative phase of the monopole (relative to other monopoles)"""

        def __init__(self, x: float=0.0, y: float=0.0, z: float=0.0, amplitude: float=1.0, phase: float=0.0):
            self.x = x
            self.y = y
            self.z = z
            self.amplitude = amplitude
            self.phase = phase

        def set_xyz(self, x: float=0.0, y: float=0.0, z: float=0.0):
            self.x = x
            self.y = y
            self.z = z

        def set_amp_phase(self,amplitude: float=1.0, phase: float=0.0):
            self.amplitude = amplitude
            self.phase = phase

    class RxPlane:
        """describes a rectangular XZ section of the plane - and points of sampling(calculation) of the pressure"""
        x0 = -0.01e-3  # needs to be a little bit off plane to avoid division by zer
        "x-origin of the measurement plane"
        y0 = 1.0e-3
        "y-origin of the measurement plane"
        z0 = +0.0e-3
        "z-origin of the measurement plane"
        dx = +0.5e-3
        "x-step between the points in the measurement plane"
        dy = +0.5e-3
        "y-step between the points in the measurement plane"
        dz = +1.0e-3
        "z-step between the points in the measurement plane"
        nx = 1
        "number of measurement points along x-axis"
        ny = 240
        "number of measurement points along y-axis"
        nz = 160
        "number of measurement points along z-axis"
        pressurefield = numpy.zeros((ny, nz), numpy.complex)
        "the numpy array to store the result of measurement. Note: complex numbers"

        execution_time = numpy.NaN

        def __init__(self,
                     x0=-0.01e-3, y0=1.0e-3, z0=0.0e-3,
                     dx=0.5e-3, dy=0.5e-3, dz=0.5e-3,
                     nx=1, ny=1024, nz=512):
            """Creates the new instance of RxPlane"""
            self.x0 = x0
            self.y0 = y0
            self.z0 = z0
            self.dx = dx
            self.dy = dy
            self.dz = dz
            if nx != 1:
                raise AttributeError("nx must be 1: ")
            self.nx = nx

            if (ny > 4096) | (ny < 1):
                raise AttributeError("ny must have a sensible value between 1 an 4096")
            self.ny = ny

            if (nz > 4096) | (nz < 1):
                raise AttributeError("nz must have a sensible value between 1 and 4096")
            self.nz = nz

            self.clear_pressurefield()

        def verify_plane_endpoints(self):
            """"returns x,y,z location of the corner oposite to the origin corner"""
            return [self.x0 + self.nx * self.dx, self.y0 + self.ny * self.dy, self.z0 + self.nz * self.dz]

        def clear_pressurefield(self):
            """"empties the result"""
            self.pressurefield = numpy.zeros((self.ny, self.nz), numpy.complex)

    # from this line, the list of actual items starts

    """List of monopoles. List of TxElement() - the list of monopoles emitting energy into the medium"""
    elements = [
        TxElement(0.0,  5.0e-3, -25.0e-3, 1.0, 0.0),
        TxElement(0.0,  4.0e-3, -20.0e-3, 1.0, 0.0),
        TxElement(0.0,  3.0e-3, -15.0e-3, 1.0, 0.0),
        TxElement(0.0,  2.0e-3, -10.0e-3, 1.0, 0.0),
        TxElement(0.0,  1.0e-3,  -5.0e-3, 1.0, 0.0),
        TxElement(0.0,  0.0e-3,   0.0e-3, 1.0, 0.0),
        TxElement(0.0, -1.0e-3,   5.0e-3, 1.0, 0.0),
        TxElement(0.0, -2.0e-3,  10.0e-3, 1.0, 0.0),
        TxElement(0.0, -3.0e-3,  15.0e-3, 1.0, 0.0),
        TxElement(0.0, -4.0e-3,  20.0e-3, 1.0, 0.0),
        TxElement(0.0, -5.0e-3,  25.0e-3, 1.0, 0.0)
        ]

    rxPlane = RxPlane()
    "contains the description/location of the measurement points"

    """Property of the environment: the wavenumber."""
    wavenumber = 1.0 / 1.0e-3  # wavelength = 1/wavenumber; so wavenumber = 1/wavelength

    def get_ray_count(self):
        return self.rxPlane.nx * self.rxPlane.ny * self.rxPlane.nz * len(self.elements)

    def __init__(self):
        pass


class these_are_orphans:
    # old methods, not used anymore
    # static methods

    def beamsim_timed(self, params: CueBeamWorld):
        if params is None:
            params = CueBeamWorld()

        t1 = time.clock()
        self.beamsim(params)
        t2 = time.clock()
        dt = t2 - t1
        rays_per_second = self.get_ray_count() / dt
        print("got {:0.2f} M rays per second".format(rays_per_second/1e6))
        return rays_per_second

    def beamsim_instant(self,
                         k: float = 1000.0,
                         x0: float = 0.1e-3,
                         y0: float = 1e-3,
                         z0: float = 0.0,
                         nx: int = 1,
                         ny: int = 240,
                         nz: int = 160,
                         dx: float = 1.0,
                         dy: float = 1.0e-3,
                         dz: float = 1.0e-3,
                         elements_vectorized=None):
        if elements_vectorized is None:
            elements_vectorized = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        local_parameters = CueBeamSolver2.CueBeamWorld()
        # copy all the parameters from the input to the parameter object
        local_parameters.wavenumber = k
        local_parameters.rxPlane.x0 = x0
        local_parameters.rxPlane.y0 = y0
        local_parameters.rxPlane.z0 = z0
        local_parameters.rxPlane.nx = nx
        local_parameters.rxPlane.ny = ny
        local_parameters.rxPlane.nz = nz
        local_parameters.rxPlane.dx = dx
        local_parameters.rxPlane.dy = dy
        local_parameters.rxPlane.dz = dz
        assert (math.modf(len(elements_vectorized))[0] == 0, "there must be n*6 in the vector")
        element_count = math.floor(len(elements_vectorized) / int(6))
        local_parameters.elements.clear()
        for idx_element in range(0, element_count):
            element_pointer = 6 * idx_element
            new_element = CueBeamSolver2.CueBeamWorld.TxElement(
                x=elements_vectorized[element_pointer + 0],
                y=elements_vectorized[element_pointer + 1],
                z=elements_vectorized[element_pointer + 2],
                amplitude=elements_vectorized[element_pointer + 3],
                phase=elements_vectorized[element_pointer + 4]
            )
            local_parameters.elements.append(copy.copy(new_element))
        local_parameters.rxPlane.clear_pressurefield()

        # call the beamsim as is
        local_parameters.rxPlane.pressurefield = beamsim_static(local_parameters)

        return local_parameters.rxPlane.pressurefield

    def beamsim_simpler(self,
                        k: float=1000.0,
                        x0: float=0.1e-3,
                        y0: float = 1e-3,
                        z0: float = 0.0,
                        nx: int = 1,
                        ny: int = 240,
                        nz: int = 160,
                        dx: float = 1.0,
                        dy: float = 1.0e-3,
                        dz: float = 1.0e-3,
                        elements_vectorized=None):
        """wrapper for other methods so that this method returns a good result in a single call"""
        if elements_vectorized is None:
            elements_vectorized = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]

        self.wavenumber = k
        self.rxPlane.x0 = x0
        self.rxPlane.y0 = y0
        self.rxPlane.z0 = z0
        self.rxPlane.nx = nx
        self.rxPlane.ny = ny
        self.rxPlane.nz = nz
        self.rxPlane.dx = dx
        self.rxPlane.dy = dy
        self.rxPlane.dz = dz
        # now, the format for elements_vectorized is mapped after matlab's original version:
        # it is 1D, 1x 6n list where n=number of elements,
        # and the items are float (x,y,z,amplitude, phase, reserved)
        assert(math.modf(len(elements_vectorized))[0] == 0, "there must be n*6 in the vector")
        # print("len(elements_vectorized) = {}".format(len(elements_vectorized) ))
        element_count = math.floor(len(elements_vectorized) / int(6))
        # print("calculated element count: {}".format(element_count))
        self.elements.clear()
        for idx_element in range(0,element_count):
            element_pointer = 6*idx_element
            new_element = CueBeamSolver.TxElement(
                        x=elements_vectorized[element_pointer + 0],
                        y=elements_vectorized[element_pointer + 1],
                        z=elements_vectorized[element_pointer + 2],
                        amplitude=elements_vectorized[element_pointer + 3],
                        phase=elements_vectorized[element_pointer + 4]
                         )
            self.elements.append(copy.copy(new_element))
            # print("added element {}, z={}, a={}".format(idx_element,self.elements[-1].z,self.elements[-1].amplitude))

        self.rxPlane.clear_pressurefield()
        self.beamsim(self)
        return self.rxPlane.pressurefield


    def beamsim(self):
        """Execute the pressure calculation.
        
        Executes the calculation using settings stored in elements, rxPlane and wavenumber.
        
        Parameters
        ----------
        None. Takes input from the elements, rxPlane and wavenumber.
        
        Returns
        -------
        None. The result is stored in .pressurefield field
        
        
        """

        # note: there is only a single plane implemented. This is for compatibility with the old Matlab-Cuda version.
        # Essentially, there should be no problem in extending this to 3D computation as needed.
        local_rx_plane = self.rxPlane
        local_elements = copy.copy(self.elements)
        local_wavenumber = self.wavenumber
        
        # debug code:
        # print("local elements:{}".format(len(local_elements)))
        # for element in local_elements:
        #     print("z:{}, a:{}".format(element.z, element.amplitude))

        ix = 0
        for iz in range(local_rx_plane.nz):
            for iy in range(local_rx_plane.ny):
                # initialize accumulator
                pressure_re = 0.0
                pressure_im = 0.0
                # calculate the XYZ location of the transmitter
                pixel_x = local_rx_plane.x0 + ix * local_rx_plane.dx
                pixel_y = local_rx_plane.y0 + iy * local_rx_plane.dy
                pixel_z = local_rx_plane.z0 + iz * local_rx_plane.dz
                # for each element do:
                for itx in range(len(local_elements)):
                    # calculate distance between the transmitter element and receiver pixel
                    ddx = (pixel_x - local_elements[itx].x)
                    ddy = (pixel_y - local_elements[itx].y)
                    ddz = (pixel_z - local_elements[itx].z)
                    distance = math.sqrt(ddx * ddx + ddy * ddy + ddz * ddz)
                    # based on distance and wavenumber, calculate the amplitude and phase
                    k_phase = - local_wavenumber * distance + local_elements[itx].phase
                    k_amplitude = local_elements[itx].amplitude / (2 * 3.14159 * distance)
                    # increment the accumulator
                    pressure_re += math.cos(k_phase) * k_amplitude
                    pressure_im += math.sin(k_phase) * k_amplitude
                # memory write
                local_rx_plane.pressurefield[(iy, iz)] = pressure_re + 1.0j * pressure_im

                # end for iy
            # end for iz
        # store result to the object
        self.rxPlane = copy.copy(local_rx_plane)
        # end beamsim


    def get_ray_count(self):
        return self.rxPlane.nx * self.rxPlane.ny * self.rxPlane.nz * len(self.elements)

    def do_plot_real(self):
        """Makes an example plot.
        
         It is configured for the default data. Compatible with Jupyter"""

        hfig = plt.figure(num=1, figsize=(8, 6), dpi=90, facecolor='white', edgecolor='black')
        imgplot = plt.imshow(
            X=numpy.real(self.rxPlane.pressurefield),
            extent=(
                self.rxPlane.z0, self.rxPlane.z0 + self.rxPlane.nz * self.rxPlane.dz,
                self.rxPlane.y0 + self.rxPlane.ny * self.rxPlane.dy, self.rxPlane.y0
            ),
            interpolation="spline36",
            clim=(-8.0, 8.0),
            origin="lower")
        # end imshow
        plt.set_cmap("PRGn")  # purple-white-green color map
        plt.xlabel("z-axis[m]")
        plt.ylabel("y-axis[m]")
        plt.show()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
@CELERY.task()
def beamsim(world: CueBeamWorld):
    if world is None:
        warnings.warn("creating new world with default settings")
        world = CueBeamWorld()

    ix = 0
    for iz in range(world.rxPlane.nz):
        for iy in range(world.rxPlane.ny):
            # initialize accumulator
            pressure_re = 0.0
            pressure_im = 0.0
            # calculate the XYZ location of the transmitter
            pixel_x = world.rxPlane.x0 + ix * world.rxPlane.dx
            pixel_y = world.rxPlane.y0 + iy * world.rxPlane.dy
            pixel_z = world.rxPlane.z0 + iz * world.rxPlane.dz
            # for each element do:
            for itx in range(len(world.elements)):
                # calculate distance between the transmitter element and receiver pixel
                ddx = (pixel_x - world.elements[itx].x)
                ddy = (pixel_y - world.elements[itx].y)
                ddz = (pixel_z - world.elements[itx].z)
                distance = math.sqrt(ddx * ddx + ddy * ddy + ddz * ddz)
                # based on distance and wavenumber, calculate the amplitude and phase
                k_phase = - world.wavenumber * distance + world.elements[itx].phase
                k_amplitude = world.elements[itx].amplitude / (2 * 3.14159 * distance)
                # increment the accumulator
                pressure_re += math.cos(k_phase) * k_amplitude
                pressure_im += math.sin(k_phase) * k_amplitude
            # memory write
            world.rxPlane.pressurefield[(iy, iz)] = pressure_re + 1.0j * pressure_im
    return world


@CELERY.task()
def beamsim_instant(self,
                    k: float = 1000.0,
                    x0: float = 0.1e-3,
                    y0: float = 1e-3,
                    z0: float = 0.0,
                    nx: int = 1,
                    ny: int = 240,
                    nz: int = 160,
                    dx: float = 1.0,
                    dy: float = 1.0e-3,
                    dz: float = 1.0e-3,
                    elements_vectorized=None):
        if elements_vectorized is None:
            elements_vectorized = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        local_world = CueBeamWorld()
        # copy all the parameters from the input to the parameter object
        local_world.wavenumber = k
        local_world.rxPlane.x0 = x0
        local_world.rxPlane.y0 = y0
        local_world.rxPlane.z0 = z0
        local_world.rxPlane.nx = nx
        local_world.rxPlane.ny = ny
        local_world.rxPlane.nz = nz
        local_world.rxPlane.dx = dx
        local_world.rxPlane.dy = dy
        local_world.rxPlane.dz = dz
        assert (math.modf(len(elements_vectorized))[0] == 0, "there must be n*6 in the vector")
        element_count = math.floor(len(elements_vectorized) / int(6))
        local_world.elements.clear()
        for idx_element in range(0, element_count):
            element_pointer = 6 * idx_element
            new_element = CueBeamSolver2.CueBeamWorld.TxElement(
                x=elements_vectorized[element_pointer + 0],
                y=elements_vectorized[element_pointer + 1],
                z=elements_vectorized[element_pointer + 2],
                amplitude=elements_vectorized[element_pointer + 3],
                phase=elements_vectorized[element_pointer + 4]
            )
            local_world.elements.append(copy.copy(new_element))
        local_world.rxPlane.clear_pressurefield()

        # call the beamsim as is
        local_world.rxPlane.pressurefield = beamsim(local_world)

        return local_world.rxPlane.pressurefield


def do_plot_abs(the_input):
        """Makes an example plot.

         It is configured for the default data. Compatible with Jupyter"""
        pressurefield = None

        if the_input is None:
            raise Exception("You must supply a pressurefield or world:cueBeamCore2.CueBeamWorld")

        if (type(the_input) is CueBeamWorld):
            world = the_input
            pressurefield = world.rxPlane.pressurefield

        if (type(the_input) is numpy.ndarray):
            world = CueBeamWorld() # create new, default world
            pressurefield = the_input

        if pressurefield is None:
            raise Exception("Something wrong: pressurefield is still None")

        hfig = plt.figure(num=1, figsize=(8, 6), dpi=90, facecolor='white', edgecolor='black')

        imgplot = plt.imshow(
            X=numpy.abs(pressurefield),
            extent=(
                world.rxPlane.z0, world.rxPlane.z0 + world.rxPlane.nz * world.rxPlane.dz,
                world.rxPlane.y0, world.rxPlane.y0 + world.rxPlane.ny * world.rxPlane.dy
            ),
            interpolation="spline36",
            clim=(0, 8.0),
            origin="lower")
        # end imshow
        plt.set_cmap("plasma") # black-to-yellow color map
        plt.xlabel("z-axis[m]")
        plt.ylabel("y-axis[m]")
        plt.show()

# to test, do :
# import cueBeamCore2
# w = cueBeamCore2.CueBeamWorld()
# cueBeamCore2.do_plot_abs(cueBeamCore2.beamsim(w).rxPlane.pressurefield)