# import pycuda.autoinit # do not autoinit - I will initialize manually
import bugcatcher

import random
import time
import math

import numpy

try:
    import pycuda.driver as drv
    from pycuda.compiler import SourceModule
except:
    pass

import redis
from celery import Celery
from cuebeam import cueBeamWorld
from cuebeam.cueBeamWorld import example_plot


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # CELERY SETTINGS # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# storage for dynamically stored functions
# RSTORE = redis.StrictRedis('EEE-mimir.ds.strath.ac.uk', port=6379, db=1)
RSTORE = redis.StrictRedis('eee-ultra26.ds.strath.ac.uk', port=6379, db=1)
# connection to the celery object, broker and broker's back-end
# CELERY = Celery('tasks', broker='amqp://guest:guest@192.168.0.43:5672/',backend='redis://192.168.0.43:6379/0')
#CELERY = Celery('tasks', broker='amqp://guest:guest@EEE-mimir.ds.strath.ac.uk/',backend='redis://EEE-mimir.ds.strath.ac.uk/0')
CELERY = Celery('tasks', broker='amqp://guest:guest@eee-ultra26.ds.strath.ac.uk/',backend='redis://eee-ultra26.ds.strath.ac.uk/0')

CELERY.conf.accept_content = ['json', 'msgpack']
CELERY.conf.result_serializer = 'msgpack'

# for debugging:
cueBeamCore3Verbosity = False


def enable_verbosity():
    global cueBeamCore3Verbosity
    cueBeamCore3Verbosity = True


def disable_verbosity():
    global cueBeamCore3Verbosity
    cueBeamCore3Verbosity = False


def output3values():
    a = 3
    b = 4
    c = 5
    return a, b, c

# ############################################
# The wrapper for calling the beamsim remotely and returning the field only
# this way matlab is completely relieved from dealing with python objects
# ############################################
def beamsim_remote(k: float = 1000.0,
                   x0: float = 0.1e-3,
                   y0: float = 1e-3,
                   z0: float = 0.0,
                   nx: int = 1.0,
                   ny: int = 240,
                   nz: int = 160,
                   dx: float = 1.0e-3,
                   dy: float = 1.0e-3,
                   dz: float = 1.0e-3,
                   elements_vectorized=None):
    if cueBeamCore3Verbosity:
        print("calling remote worker and waiting")
    async_result = beamsim_through_celery.delay(k, x0, y0, z0, nx, ny, nz, dx, dy, dz, elements_vectorized)
    return async_result.get(interval=0.1)
    #
    #WaitCounter=0
    #while not(async_result.ready()):
    #    time.sleep(0.01)  # check up to 100x per second
    #    WaitCounter += 1

    #print('waited {} times.'.format(WaitCounter))
    #return async_result.result


# ############################################
# The wrapper for calling the beamsim remotely and returning the field only
# this way matlab is completely relieved from dealing with python objects
# ############################################
def beamsim_lambert_remote(
        elements_vectorized=None,
        k: float = 1000.0,
        lambert_radius: float = 100e-3,
        lambert_map_density: float = 1e-3):
    if cueBeamCore3Verbosity:
        print("calling remote worker and waiting")
    async_result = beamsim_lambert_through_celery.delay(elements_vectorized,k,lambert_radius,lambert_map_density)
    #while not(async_result.ready()):
    #    time.sleep(0.01)  # check up to 100x per second
    #return async_result.result
    return async_result.get(interval=0.1)


@CELERY.task()
def beamsim_through_celery(k: float = 1000.0,
                           x0: float = 0.1e-3,
                           y0: float = 1e-3,
                           z0: float = 0.0,
                           nx: int = 1,
                           ny: int = 240,
                           nz: int = 160,
                           dx: float = 1.0e-3,
                           dy: float = 1.0e-3,
                           dz: float = 1.0e-3,
                           elements_vectorized=None):
    if cueBeamCore3Verbosity:
        print('got work, {} beams.'.format(nx*ny*nz*(len(elements_vectorized)/6)))

    return beamsim_gpu(k, x0, y0, z0, nx, ny, nz, dx, dy, dz, elements_vectorized)


@CELERY.task()
def beamsim_lambert_through_celery(
        elements_vectorized = None,
        k: float = 1000.0,
        lambert_radius: float = 100e-3,
        lambert_map_density: float = 1e-3):
    return cuebeamlambert(elements_vectorized, k, lambert_radius,lambert_map_density)





# same as previous, but this time, nx>1 is allowed
def beamsim_3d(k:  float = 1000.0,
               x0: float = 0.1e-3,
               y0: float = 1e-3,
               z0: float = 0.0,
               nx: int = 160,
               ny: int = 240,
               nz: int = 160,
               dx: float = 1.0e-3,
               dy: float = 1.0e-3,
               dz: float = 1.0e-3,
               elements_vectorized=None):

    if elements_vectorized is None:
        elements_vectorized1 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    else:
        elements_vectorized1 = elements_vectorized

    # create the space for the result
    cout3d = numpy.zeros((int(nx), int(ny), int(nz))).astype(numpy.complex64)
    async_handles=[]  # placeholder for the async results
    # send out the work
    for idxX in range(nx):
        local_x0 = x0 + idxX * dx
        async_handles.append(
            {
                'handle': beamsim_through_celery.delay(k, local_x0, y0, z0, int(1), ny, nz, dx, dy, dz, elements_vectorized1),
                'idxX': idxX,
                'done': False
            })
    # read out the results
    results_to_read = nx
    while results_to_read>0:
        for idxX in range(nx):
            if not(async_handles[idxX]['done']):
                if async_handles[idxX]['handle'].ready():
                    cout3d[async_handles[idxX]['idxX'], :, :] = async_handles[idxX]['handle'].result
                    async_handles[idxX]['done'] = True
                    results_to_read -= 1
    return cout3d





# ############################################
# The actual beamsim function
# ############################################
def beamsim_gpu(k: float = 1000.0,
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

    # note, that on a remote worker, there is never a saying where this thread wakes up
    # therefore i must be extra carefull to always initialize all the resources needed

    # initialize manually
    drv.init()
    # local_device_count = drv.Device.count()
    # print('found {} GPUs'.format(local_device_count))
    # choose one of the GPUs at random
    gpu_to_take = random.choice(range(drv.Device.count()))
    # take device 0 for now
    gpu_context = drv.Device(gpu_to_take).make_context()
    gpu_context.push()  # make the context active
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
    
    // start actual calculation : zero the accumulators
    pressure_re = 0;
    pressure_im = 0;
    
    // where am I in space?
    pixel_x = (float)ix * dx + x0;
    pixel_y = (float)iy * dy + y0; 
    pixel_z = (float)iz * dz + z0;
    
    // debugging code only:
    // printf("block %d.%d: pixel  %d.%d.%d, at %0.3f|%0.3f|%0.3f\\n",blockIdx.y,blockIdx.z, ix, iy, iz, pixel_x,pixel_y,pixel_z); // note that enabling this makes this a long call, time outs the driver
    
    // for each transmitter-element, do this:
    for (itx=0; itx<tx_length*6; itx=itx+6) // this hopefully acesses the same memory location for each thread and therefore will be cached
     {            
        // calculate distance in cartesian space:
        dix = (pixel_x-tx[0+itx]);                            // tx.x
        diy = (pixel_y-tx[1+itx]);                            // tx.y
        diz = (pixel_z-tx[2+itx]);                            // tx.z
        distance = sqrtf( dix*dix + diy*diy + diz*diz );
        // amplitude decays with distance as the energy gets distributed on a ring around the transmit point
        // note that ring is for 2D space, and a sphere surface would be more appropriate for 3D space
        
        amplitude = tx[3+itx] / ( pi2 * distance );                   // amplitude is at itx+3  
        kd = -k * distance + tx[4 + itx];                           // phase is at itx+4            
        
        // accumulate the energy                           
        pressure_re = pressure_re + __cosf(kd) * amplitude;                    
        pressure_im = pressure_im + __sinf(kd) * amplitude; 
     }        
     
    // write the result:
        
    // in case if I want the absolute pressure only (discards the phase)
    // pressure=sqrtf(pressure_re*pressure_re+pressure_im*pressure_im); 
           
    // in CUDA, i need to calculate rx array memory offset manually for each thread:
    // offset=ix+iy*nx+iz*nx*ny; // that's a version for xyz, real-only numbers (e.g. pycuda.float32 version
    // out[offset]=(float)pressure;
               
    // this is a version for complex numbers: pycuda.complex64
    offset=2*(iz+iy*nz+ix*nz*ny); // that's a version for xyz version    
    out[offset]=(float)pressure_re;
    offset++; // go to the imaginary value pointer
    out[offset]=(float)pressure_im;    
}
""")

    # instantiate the code into the compiler
    beamsim_kernel = code_text.get_function("BeamSimKernel")

    # convert the values from pythonic to Cudific
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
    # note: must reserve the output memory right here
    # note: for 2D values, the x must be == 1
    assert(nx == 1), 'this version supports nx=1 only'
    # prevent from the image size being too large
    assert(ny < 4096+1), ' ny too large: reduce calculated field pixel count'
    assert(nz < 4096+1), ' nz too large: reduce calculated field pixel count'
    # prevent from the transducer description to be too large - you can remove this limitation later on
    assert(ctx_count < 81920+1),'too many radiator points. Reduce simulation complexity'
    cuda_out = numpy.zeros((int(cny), int(cnz))).astype(numpy.complex64)

    # prepare the GPU call : thread wave shape:
    threads_x = 1
    threads_y = 16
    threads_z = 64
    blocks_x = 1
    blocks_y = int((int(cny) / threads_y) + 1)
    blocks_z = int((int(cnz) / threads_z) + 1)

    # start the timer!
    # time_1 = time.clock()

    beamsim_kernel(
        drv.In(ctx),
        ctx_count,
        drv.Out(cuda_out),
        cx0, cy0, cz0,
        cnx, cny, cnz,
        cdx, cdy, cdz,
        ck,
        block=(threads_x, threads_y, threads_z), grid=(blocks_x, blocks_y, blocks_z))

    # time_2 = time.clock()

    # release the GPU from this thread
    # release the context, otherwise memory leak might occur
    gpu_context.pop()
    gpu_context.detach()

    # performance = numpy.float64(numpy.int128(cnx)*numpy.int128(cny)*numpy.int128(cnz)*numpy.int128(ctx_count)) / (time_2 - time_1)

    return cuda_out


def cuebeamlambert(
    elements_vectorized=None,
    k: float = 1000.0,
    lambert_radius: float = 100e-3,
    lambert_map_density: float = 100):
    # print("lambert_radius = {}, lambert_map_density={}".format(lambert_radius,lambert_map_density))
    # matlab:: [img_lambert lambert_x lambert_y lambert_z] = cueBeam.cueBeam_lambert(tx',enviroment.wavenumber,lambert_radius,lambert_map_density);
    t0 = time.time()

    if elements_vectorized is None:
        elements_vectorized = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]

    # note, that on a remote worker, there is never a saying where this thread wakes up
    # therefore i must be extra carefull to always initialize all the resources needed

    drv_init_time_before = time.time()-t0
    # initialize manually

    drv.init()

    # local_device_count = drv.Device.count()
    # print('found {} GPUs'.format(local_device_count))

    # choose one of the GPUs at random
    gpu_to_take = random.choice(range(drv.Device.count()))
    gpu_context = drv.Device(gpu_to_take).make_context()
    gpu_context.push()  # make the context active

    drv_init_time_after = time.time() - t0

    code_text_lambert=SourceModule("""
    
    #include <stdio.h>

#define pi2 2.0f*3.141592653589793f

    __global__ void BeamsimLambertKernel ( float *tx, unsigned int tx_length, float *out,
        unsigned int n, float d, float r,
        float k,float *xp, float *yp, float *zp)
{
    unsigned int offset=0;
    unsigned int ix,iy,itx=0;
    float pressure,distance,kd,pressure_re,pressure_im=0;
    float dist2=0;
    float dix,diy,diz,lambert_x,lambert_y,lambert_z=0;
    float xbase,ybase,rho2,rhoi,cosphi,sinphi,cosl,sinl=0;
    float xbase0=-sqrtf((float)2)+(float)1e-8;
    // calculate ix,iy from thread built-in variables
    ix = blockIdx.x * blockDim.x + threadIdx.x;
    iy = blockIdx.y * blockDim.y + threadIdx.y;
    //ix=0; // debug //
    //C// for (iy=0; iy<ny; iy++)
    //C//  for (ix=0; ix<nx; ix++)
    
    // make sure that this thread won't try to calculate non-existent receiver
    if (iy>n) return;
    
    if (ix>n) return;
    
    // start actual calculation
    pressure_re=0;
    pressure_im=0;
    
    xbase=(float)ix*d+xbase0;
    ybase=(float)iy*d+xbase0; // it would be an optimisation not to recalculate it for each pixel, this has to stay here be due to future port to CUDA where each pixel has it's own thread
    rho2=xbase*xbase+ybase*ybase;
    offset=ix+n*iy;
    if (rho2>(float)2)
    {
        out[2*offset]=0;
        out[2*offset+1]=0; // note, complex valued output
        
        xp[offset]=0;
        yp[offset]=0;
        zp[offset]=0;
        return;
    }
    rhoi=rsqrtf(rho2);
    cosl=-ybase*rhoi;
    cosphi=sqrtf(rho2-rho2*rho2/(float)4);
    lambert_x=r*cosl*cosphi;
    
    sinl=xbase*rhoi;
    lambert_y=r*sinl*cosphi;
    sinphi=(float)1-rho2/(float)2;
    lambert_z=r*sinphi;
    
    xp[offset]=lambert_x;
    yp[offset]=lambert_y;
    zp[offset]=lambert_z;
    
    for (itx=0; itx<tx_length*6; itx=itx+6) // this hopefully acesses the same memory location for each thread and therefore will be cached
    {
        //    distance=single(sqrt( (ix*dx+x0-tx(1,itx)).^2 + (iy*dy+y0-tx(2,itx)).^2 + (iz*dz+z0-tx(3,itx)).^2 ));
        dix=(lambert_x-tx[0+itx]);
        diy=(lambert_y-tx[1+itx]);
        diz=(lambert_z-tx[2+itx]);
        distance=sqrtf( dix*dix + diy*diy + diz*diz );
        kd=-k*distance+tx[5+itx];
        dist2=tx[4+itx]/(6.283185307179586f*distance); //equals 2*pi                
        pressure_re=pressure_re+__cosf(kd)*dist2;
        pressure_im=pressure_im+__sinf(kd)*dist2;
        
        // note: __sinf is an simlpified sin function that yields less accurate result. May need to switch to full sin for final product, ok for testing for now
        
        // note 2: function sincosf(...) may be faster in this case - calculates both sin and cos. but since it requires additional accumulators,
        // a detailed test will be required to find out what's faster.
        
    }
    // mem write    
    // pressure=sqrtf(pressure_re*pressure_re+pressure_im*pressure_im);
    // in CUDA, i need to calculate rx array memory offset manually for each thread:
    //offset=ix+nx*iy+(ny*nx)*iz;    
    // out[offset]=(float)pressure; //left for debug    
    // offset=2*(ix+n*iy);  // note, complex-valued version A 
    out[2*offset+0]=(float)pressure_re; // use real-to-complex offset conversion
    offset++; // go to the imaginary value pointer
    out[2*offset+1]=(float)pressure_im;
}

    """)

    # instantiate the code into the compiler
    beamsim_lambert_kernel = code_text_lambert.get_function("BeamsimLambertKernel")
    kernel_init_time = time.time() - t0

    # calc basic space properties
    npts = float(math.ceil(6.283185307179586 * lambert_radius / lambert_map_density))
    d = 2.0 * math.sqrt(2)/npts
    n = 1 + math.ceil(2*math.sqrt(2)/d)

    # convert the values from pythonic to Cudific
    ctx = numpy.asarray(elements_vectorized).astype(numpy.float32)
    ctx_count = numpy.int32(len(ctx) / 6)
    cn = numpy.int32(n)
    cd = numpy.float32(d)
    cr = numpy.float32(lambert_radius)
    ck = numpy.float32(k)
    cuda_out_xp = numpy.zeros((int(n), int(n))).astype(numpy.float32)
    cuda_out_yp = numpy.zeros((int(n), int(n))).astype(numpy.float32)
    cuda_out_zp = numpy.zeros((int(n), int(n))).astype(numpy.float32)
    cuda_out = numpy.zeros((int(n), int(n))).astype(numpy.complex64) # note: must reserve the output memory right here

    # prevent from the transducer description to be too large - you can remove this limitation later on
    assert (ctx_count < 300001), "transducer definition too large"

    # prepare the GPU call : thread wave shape:
    threads_x = 16
    threads_y = 64
    threads_z = 1
    blocks_x = int((int(n) / threads_x) + 1)
    blocks_y = int((int(n) / threads_y) + 1)
    blocks_z = 1

    kernel_prepare_time = time.time() - t0
    beamsim_lambert_kernel(
        drv.In(ctx),
        ctx_count,
        drv.Out(cuda_out),
        cn, cd, cr, ck,
        drv.Out(cuda_out_xp),
        drv.Out(cuda_out_yp),
        drv.Out(cuda_out_zp),
        block=(threads_x, threads_y, threads_z),
        grid=(blocks_x, blocks_y, blocks_z))

    kernel_run_time=time.time() - t0


    # release the GPU from this thread
    # release the context, otherwise memory leak might occur
    gpu_context.pop()
    gpu_context.detach()
    detach_time = time.time() - t0

    # calculate performance metrics. This is usefull for debugging the distributed computation system
    dinittime = drv_init_time_after - drv_init_time_before
    dkernel_init_time = kernel_init_time - drv_init_time_after
    dkernel_prepare_time = kernel_prepare_time - kernel_init_time
    dkernel_run_time = kernel_run_time - kernel_prepare_time
    dkernel_detach_time = detach_time - kernel_run_time
    print("lambert: init:{:06.4f}, kernel_init:{:06.4f}, kernel_prepare:{:06.4f}, kernel_run:{:06.4f}, detach:{:06.4f}".format(dinittime,dkernel_init_time,dkernel_prepare_time,dkernel_run_time,dkernel_detach_time))

    # finally...
    return cuda_out, cuda_out_xp, cuda_out_yp, cuda_out_zp


def cueBeamDemo():
    world = cueBeamWorld.CueBeamWorld()
    world.rxPlane.ny = 1024
    world.rxPlane.nz = 256
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

    # TODO: Select local or remote execution path depending on which resource is available

    pressurefield = beamsim_remote(world.wavenumber,
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
    if cueBeamCore3Verbosity:
        time_total = time_end - time_start
        data_transfer_performance = pressurefield.shape[0]*pressurefield.shape[1]*8 / time_total
        print("got {:7.3f} MRays/s, {:7.1f} MB/sec".format(local_performance*1e-6,data_transfer_performance/1024/1024))
    world.rxPlane.pressurefield = pressurefield
    example_plot(world)


