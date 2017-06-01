import time

from cueBeamWorld import CueBeamWorld
from cueBeamWorld import example_plot

from cueBeamCore3 import beamsim_gpu
from cueBeamCore3 import beamsim_through_celery

world = CueBeamWorld()
world.rxPlane.ny = 4096
world.rxPlane.nz = 2048
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

async_result = beamsim_through_celery.delay(world.wavenumber,
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
while not(async_result.ready()):
    time.sleep(0.01)

pressurefield = async_result.result

time_end = time.clock()
time_total = time_end - time_start
remote_performance = world.get_ray_count() / time_total
data_transfer_performance = pressurefield.shape[0]*pressurefield.shape[1]*8 / time_total

print("got {:7.3f} MRays/s, {:7.1f} MB/sec".format(remote_performance*1e-6, data_transfer_performance/1024/1024))

world.rxPlane.pressurefield = pressurefield
example_plot(world)

# now try the local speed

time_start = time.clock()

# TODO: Select local or remote execution path depending on which resource is available

pressurefield = beamsim_gpu(world.wavenumber,
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
time_total = time_end - time_start
local_performance = world.get_ray_count() / time_total
data_transfer_performance = pressurefield.shape[0]*pressurefield.shape[1]*8 / time_total

print("got {:7.3f} MRays/s, {:7.1f} MB/sec".format(local_performance*1e-6, data_transfer_performance/1024/1024 ))

world.rxPlane.pressurefield = pressurefield
example_plot(world)

