import math
import numpy
import math
import matplotlib.pyplot as plt
import time
import copy
import warnings
import dill
import pickle
dill.settings['recurse'] = True
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

        def set_x0y0n0(self, x0: float=0.1e-3, y0: float = 1.0e-3, z0: float = 0.0e-3):
            self.x0 = x0
            self.y0 = y0
            self.z0 = z0

        def set_nxnynz(self, nx:int = 1, ny: int = 512, nz: int = 160):
            self.nx = nx
            self.ny = ny
            self.nz = nz

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

    """just in case if i want a benchmark of the given node"""
    last_performance_rays_per_second = numpy.NaN

    def get_ray_count(self):
        return self.rxPlane.nx * self.rxPlane.ny * self.rxPlane.nz * len(self.elements)

    def __init__(self):
        pass


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# @CELERY.task()
def beamsim(world: CueBeamWorld):
    t_before = time.clock()
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
    t_after = time.clock()
    execution_time = t_after - t_before
    world.last_performance_rays_per_second = world.get_ray_count() / execution_time
    return world.rxPlane.pressurefield


@CELERY.task()
def beamsim_instant(k: float = 1000.0,
                    x0: float = 0.1e-3,
                    y0: float = 1e-3,
                    z0: float = 0.0,
                    nx: int = 1,
                    ny: int = 240,
                    nz: int = 160,
                    dx: float = 1.0,
                    dy: float = 1.0e-3,
                    dz: float = 1.0e-3,
                    elements_vectorized = None):
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
            new_element = CueBeamWorld.TxElement(
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

def send_and_receive_many(world: CueBeamWorld):
    """ uses the beamsim_instant call method"""
    elements_vectorized1 = []
    for idxElement in range(0, len(world.elements) - 1):
        elements_vectorized1.extend(
            [world.elements[idxElement].x, world.elements[idxElement].y, world.elements[idxElement].z,
             world.elements[idxElement].amplitude, world.elements[idxElement].phase, 0.0])
    time_start = time.clock()
    current_ray_count = world.get_ray_count()
    estimated_worker_performance = 300000.0
    need_workers = math.ceil(current_ray_count / estimated_worker_performance)
    each_worker_does_ylines = math.ceil(world.rxPlane.ny / need_workers )
    # update
    handles = []
    for idx_worker in range(need_workers):
        yline0 = idx_worker*each_worker_does_ylines # starts at zero
        yline_y = world.rxPlane.y0 + world.rxPlane.dy * yline0
        handles.append({
                        'yline_y': yline_y,
                        'async_handle': beamsim_instant.delay(
                                                            k=world.wavenumber,
                                                            x0=world.rxPlane.x0,
                                                            y0=yline_y,
                                                            z0=world.rxPlane.z0,
                                                            nx=world.rxPlane.nx,
                                                            ny=each_worker_does_ylines,
                                                            nz=world.rxPlane.nz,
                                                            dx=world.rxPlane.dx,
                                                            dy=world.rxPlane.dy,
                                                            dz=world.rxPlane.dz,
                                                            elements_vectorized=elements_vectorized1)
                        })
    # TODO: FRONTIER HERE ===================

    # TODO: Wait for first worker, and load the result,
    #while not (async_handle.ready()):
    #    time.sleep(0.02)

    world.rxPlane.pressurefield = pickle.loads(async_handle.result)
    time_end = time.clock()
    world.last_performance_rays_per_second = world.get_ray_count() / (time_end - time_start)
    print('performance = {} kRays/sec'.format(world.last_performance_rays_per_second / 1e3))
    return world


def send_and_receive(world: CueBeamWorld):
    """ uses the beamsim_instant call method"""
    elements_vectorized1 = []
    for idxElement in range(0, len(world.elements)-1):
        elements_vectorized1.extend([world.elements[idxElement].x,world.elements[idxElement].y,world.elements[idxElement].z,world.elements[idxElement].amplitude,world.elements[idxElement].phase,0.0])
    time_start = time.clock()
    async_handle = beamsim_instant.delay(
        k=world.wavenumber,
        x0=world.rxPlane.x0,
        y0=world.rxPlane.y0,
        z0=world.rxPlane.z0,
        nx=world.rxPlane.nx,
        ny=world.rxPlane.ny,
        nz=world.rxPlane.nz,
        dx=world.rxPlane.dx,
        dy=world.rxPlane.dy,
        dz=world.rxPlane.dz,
        elements_vectorized = elements_vectorized1)
    while not(async_handle.ready()):
        time.sleep(0.02)

    world.rxPlane.pressurefield = pickle.loads(async_handle.result)
    time_end = time.clock()
    world.last_performance_rays_per_second = world.get_ray_count() / (time_end-time_start)
    print('performance = {} kRays/sec'.format(world.last_performance_rays_per_second/1e3))
    return world


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
            X=numpy.real(pressurefield),
            extent=(
                world.rxPlane.z0, world.rxPlane.z0 + world.rxPlane.nz * world.rxPlane.dz,
                world.rxPlane.y0, world.rxPlane.y0 + world.rxPlane.ny * world.rxPlane.dy
            ),
            # interpolation="spline36",
            interpolation="nearest",
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