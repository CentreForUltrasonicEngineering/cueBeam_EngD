# import pycuda.autoinit # do not autoinit - I will initialize manually
import math
import pycuda.driver as drv
import numpy
import matplotlib.pyplot as plt
import time
from pycuda.compiler import SourceModule

import redis
from celery import Celery

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # CELERY SETTINGS # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# storage for dynamically stored functions
RSTORE = redis.StrictRedis('192.168.0.43', port=6379, db=1)
# connection to the celery object, broker and broker's back-end
CELERY = Celery('tasks', broker='amqp://guest:guest@192.168.0.43:5672/',backend='redis://192.168.0.43:6379/0')

CELERY.conf.accept_content = ['json', 'msgpack']
CELERY.conf.result_serializer = 'msgpack'


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # CueBeamWorld class
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class CueBeamWorld:
    """Contains all the info needed for cueBeam core code to return a pressure field
    it's actually not very compatible with the core cueBeam solver
    but it's soo incredibly convinient for storing world state , and it's default values are so usefull 
    that I decided to leave it in for now    
    """

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


def do_plot(the_input):
    """Makes an example plot.

     It is configured for the default data. Compatible with Jupyter"""
    pressurefield = None

    if the_input is None:
        raise Exception("You must supply a pressurefield or world:CueBeamWorld")

    if type(the_input) is CueBeamWorld:
        world = the_input
        pressurefield = world.rxPlane.pressurefield

    if type(the_input) is numpy.ndarray:
        world = CueBeamWorld()  # create new, default world
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
        # interpolation="spline36",
        interpolation="nearest",
        clim=(0, 8.0),
        origin="lower")
    # end imshow
    plt.set_cmap("plasma")  # black-to-yellow color map
    # plt.set_cmap("PRGn")
    plt.xlabel("z-axis[m]")
    plt.ylabel("y-axis[m]")
    plt.show()


@CELERY.task()
def beamsim_instant_gpu(k: float = 1000.0,
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

    # initialize manually
    drv.init()
    # local_device_count = drv.Device.count()
    # print('found {} GPUs'.format(local_device_count))

    # take device 0 for now
    gpu_context = drv.Device(0).make_context()  # note, could randomize here so that different calls get a different GPU
    gpu_context.push()
    code_text = SourceModule("""

#include <stdio.h>

#define pi2 2.0f*3.141592653589793f

__global__ void BeamSimKernel ( const float *tx, unsigned int tx_length, 
                                float *out,
                                float x0, float y0, float z0,
                                unsigned int nx, unsigned int ny, unsigned int nz, 
                                float dx, float dy, float dz,
                                float k )
    {
    unsigned int offset=0;
    unsigned int ix,iy,iz=0;
    unsigned int itx = 0; // used as a transducer iterator
    float amplitude = 0.0; 
    // float pressure; // no need for it anymore
    float distance,kd,pressure_re,pressure_im=0;
    // float directivity_cos; // directivity_cos not used in this version
    float pixel_x,pixel_y,pixel_z,dix,diy,diz =0;        
    // di* - delta distances as optimisation    
    // iz=0 for now in CUDA, use 2D calculation only
    
    // calculate iz,ix from thread built-in variables
    
    
    ix = blockIdx.x * blockDim.x + threadIdx.x; // use cuda.x-grid as world.x
    iy = blockIdx.y * blockDim.y + threadIdx.y; // use cuda.y-grid as world.y
    iz = blockIdx.z * blockDim.z + threadIdx.z; // use cuda.z-grid as world.z
    
    // make sure that this thread won't try to calculate non-existent receiver
    
    if (iy >= ny) return;        
    if (ix >= nx) return;    
    if (iz >= nz) return;
    
    // start actual calculation
    pressure_re = 0;
    pressure_im = 0;
    pixel_x = (float)ix * dx + x0;
    pixel_y = (float)iy * dy + y0; 
    pixel_z = (float)iz * dz + z0;
    
    // printf("block %d.%d: pixel  %d.%d.%d, at %0.3f|%0.3f|%0.3f\\n",blockIdx.y,blockIdx.z, ix, iy, iz, pixel_x,pixel_y,pixel_z); // note that enabling this makes this a long call, time outs the driver
    
    for (itx=0; itx<tx_length*6; itx=itx+6) // this hopefully acesses the same memory location for each thread and therefore will be cached
     {            
        dix = (pixel_x-tx[0+itx]);                            // tx.x
        diy = (pixel_y-tx[1+itx]);                            // tx.y
        diz = (pixel_z-tx[2+itx]);                            // tx.z
        distance = sqrtf( dix*dix + diy*diy + diz*diz ); 
        amplitude = tx[3+itx] / ( pi2 * distance );                   // amplitude is at itx+3  
        kd = -k * distance + tx[4 + itx];                           // phase is at itx+4            
                                   
        pressure_re = pressure_re + __cosf(kd) * amplitude;                    
        pressure_im = pressure_im + __sinf(kd) * amplitude; 
     }        
    // mem write
    // pressure=sqrtf(pressure_re*pressure_re+pressure_im*pressure_im);        
    // in CUDA, i need to calculate rx array memory offset manually for each thread:
    // offset=ix+iy*nx+iz*nx*ny; // that's a version for xyz version
    offset=2*(iz+iy*nz+ix*nz*ny); // that's a version for xyz version
    //out[offset]=(float)pressure;           
    out[offset]=(float)pressure_re;
    offset++; // go to the imaginary value pointer
    out[offset]=(float)pressure_im;
    //out[offset+nz*ny*nx]=(float)pressure_im;
}
""")

    BeamSimKernel = code_text.get_function("BeamSimKernel")
    ck = numpy.float32(k)
    cx0 = numpy.float32(x0)
    cy0 = numpy.float32(y0)
    cz0 = numpy.float32(z0)
    cnx = numpy.int32(nx)
    cny = numpy.int32(ny)
    cnz = numpy.int32(nz)
    cdx = numpy.float32(dx)
    cdy = numpy.float32(dy)
    cdz = numpy.float32(dz)
    ctx = numpy.asarray(elements_vectorized).astype(numpy.float32)
    ctx_count = numpy.int32(len(ctx) / 6)
    cout = numpy.zeros((int(cny), int(cnz))).astype(numpy.complex64)
    threads_x = 1
    threads_y = 16
    threads_z = 64
    blocks_x = 1
    blocks_y = int((int(cny) / threads_y) + 1)
    blocks_z = int((int(cnz) / threads_z) + 1)
    time_1 = time.clock()
    BeamSimKernel(
        drv.In(ctx),
        ctx_count,
        drv.Out(cout),
        cx0, cy0, cz0,
        cnx, cny, cnz,
        cdx, cdy, cdz,
        ck,
        block=(threads_x, threads_y, threads_z), grid=(blocks_x, blocks_y, blocks_z))
    # release the GPU from this thread
    gpu_context.detach()
    time_2 = time.clock()
    performance = (cnx*cny*cnz*ctx_count) / (time_2 - time_1)
    return cout


def demo():
    world = CueBeamWorld()
    world.rxPlane.ny = 2048
    world.rxPlane.nz = 512
    world.rxPlane.dy = 0.5e-3
    world.rxPlane.dz = 0.5e-3
    world.rxPlane.x0 = 0.0
    world.rxPlane.y0 = 1e-3
    world.rxPlane.z0 = 0e-3
    elements_vectorized_local = []
    for idxElement in range(0, len(world.elements) - 1):
        elements_vectorized_local.extend(
            [world.elements[idxElement].x,  # idx+0
             world.elements[idxElement].y,  # idx+1
             world.elements[idxElement].z,  # idx+2
             world.elements[idxElement].amplitude,  # idx+3
             world.elements[idxElement].phase,  # idx+4
             0.0]  # idx+5 - reserved for the element size/directivity for later
        )
    time_start = time.clock()
    pressurefield = beamsim_instant_gpu( world.wavenumber,
                                         world.rxPlane.x0,
                                         world.rxPlane.y0,
                                         world.rxPlane.z0,
                                         world.rxPlane.nx,
                                         world.rxPlane.ny,
                                         world.rxPlane.nz,
                                         world.rxPlane.dx,
                                         world.rxPlane.dy,
                                         world.rxPlane.dz,
                                         elements_vectorized_local)
    time_end = time.clock()
    local_performance = world.get_ray_count() / (time_end-time_start)
    print("got {:7.3f} MRays/s".format(local_performance*1e-6))
    world.rxPlane.pressurefield = pressurefield
    do_plot(world)


