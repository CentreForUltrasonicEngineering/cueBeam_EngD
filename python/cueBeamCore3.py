import pycuda.autoinit
import math
import pycuda.driver as drv
import numpy
import cueBeamCore2
import matplotlib.pyplot as plt
import time

from pycuda.compiler import SourceModule


def do_plot(the_input):
    """Makes an example plot.

     It is configured for the default data. Compatible with Jupyter"""
    pressurefield = None

    if the_input is None:
        raise Exception("You must supply a pressurefield or world:cueBeamCore2.CueBeamWorld")

    if (type(the_input) is cueBeamCore2.CueBeamWorld):
        world = the_input
        pressurefield = world.rxPlane.pressurefield

    if (type(the_input) is numpy.ndarray):
        world = cueBeamCore2.CueBeamWorld()  # create new, default world
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

def beamsim_instant2(k: float = 1000.0,
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
    pass

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
    float pressure,distance,kd,pressure_re,pressure_im=0;
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
    //pressure=sqrtf(pressure_re*pressure_re+pressure_im*pressure_im);        
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
# try to run the kernel with default world
import cueBeamCore2
world = cueBeamCore2.CueBeamWorld()
world.rxPlane.ny = 512
world.rxPlane.nz = 256
world.rxPlane.dy = 0.5e-3
world.rxPlane.dz = 0.5e-3
world.rxPlane.x0 = 0.0
world.rxPlane.y0 = 1e-3
world.rxPlane.z0 = 0e-3


# build the call
elements_vectorized1 = []
for idxElement in range(0, len(world.elements) - 1):
    elements_vectorized1.extend(
                                [world.elements[idxElement].x,          # idx+0
                                 world.elements[idxElement].y,          # idx+1
                                 world.elements[idxElement].z,          # idx+2
                                 world.elements[idxElement].amplitude,  # idx+3
                                 world.elements[idxElement].phase,      # idx+4
                                 0.0]                                   # idx+5 - reserved for the element size/directivity for later
                                )

ck = numpy.float32(world.wavenumber)
cx0 = numpy.float32(world.rxPlane.x0)
cy0 = numpy.float32(world.rxPlane.y0)
cz0 = numpy.float32(world.rxPlane.z0)
cnx = numpy.int32(world.rxPlane.nx)
cny = numpy.int32(world.rxPlane.ny)
cnz = numpy.int32(world.rxPlane.nz)
cdx = numpy.float32(world.rxPlane.dx)
cdy = numpy.float32(world.rxPlane.dy)
cdz = numpy.float32(world.rxPlane.dz)
ctx = numpy.asarray(elements_vectorized1).astype(numpy.float32)
ctx_length = numpy.int32(len(ctx)/6)
cout = numpy.zeros((int(cny), int(cnz))).astype(numpy.complex64)

threads_x = 1
threads_y = 1
threads_z = 64
blocks_x = 1
blocks_y = int((int(cny) / threads_y) + 1)
blocks_z = int((int(cnz) / threads_z) + 1)

# start actual call
benchmark_repeat = 1
time_1 = time.clock()
for idxBench in range(benchmark_repeat):
    kernel_time = BeamSimKernel(
            drv.In(ctx),
            ctx_length,
            drv.Out(cout),
            cx0, cy0, cz0,
            cnx, cny, cnz,
            cdx, cdy, cdz,
            ck,
            block=(threads_x, threads_y, threads_z), grid=(blocks_x, blocks_y,blocks_z),
            time_kernel=True)
time_2 = time.clock()
performance = world.get_ray_count() / (time_2-time_1) * benchmark_repeat

print("got {:8.3f} M rays/s".format(performance*1.0e-6))

world.rxPlane.pressurefield = cout

do_plot(world)
# const float *tx,
# unsigned int tx_length,
# float *out,
# float x0,
# float y0,
# float z0,
# unsigned int nx,
# unsigned int ny,
# unsigned int nz,
# float dx,
# float dy,
# float dz,
# float k