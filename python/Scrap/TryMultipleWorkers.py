import time
import cueBeamCore2
import copy
import math
import pickle

estimated_worker_performance = 300000.0  # rays/second
intended_worker_time = 0.3  # seconds

world = cueBeamCore2.CueBeamWorld()
world.rxPlane.set_nxnynz(1, 1024, 1024)
world.rxPlane.clear_pressurefield()
current_ray_count = world.get_ray_count()

need_workers = int(math.ceil(current_ray_count / ( estimated_worker_performance * intended_worker_time)))
each_worker_does_y_lines = math.ceil(world.rxPlane.ny / need_workers)
print("going for {} workers doing {} lines each".format(need_workers, each_worker_does_y_lines))

# update the world ny size based on the lines done
world.rxPlane.ny = int(need_workers * each_worker_does_y_lines)
world.rxPlane.clear_pressurefield()

elements_vectorized_local = []
for idxElement in range(0, len(world.elements) - 1):
    elements_vectorized_local.extend(
                                    [world.elements[idxElement].x,
                                     world.elements[idxElement].y,
                                     world.elements[idxElement].z,
                                     world.elements[idxElement].amplitude,
                                     world.elements[idxElement].phase,
                                     0.0
                                     ])

handles = []
time_start = time.clock()
for idx_worker in range(need_workers):
    y_line_n = idx_worker * each_worker_does_y_lines  # starts correctly at zero - unlike in matlab!
    y_line_y = world.rxPlane.y0 + world.rxPlane.dy * y_line_n
    handles.append({
        'received': False,
        'y_line_y': y_line_y,
        'y_line_n': y_line_n,
        'async_handle': cueBeamCore2.beamsim_instant.delay(  # note: generates and issues the remote call right here.
            k=world.wavenumber,
            x0=world.rxPlane.x0,
            y0=y_line_y,
            z0=world.rxPlane.z0,
            nx=world.rxPlane.nx,
            ny=each_worker_does_y_lines,
            nz=world.rxPlane.nz,
            dx=world.rxPlane.dx,
            dy=world.rxPlane.dy,
            dz=world.rxPlane.dz,
            elements_vectorized=elements_vectorized_local)
    })

items_to_receive = need_workers
while items_to_receive > 0:
    time.sleep(0.05)
    for idx_worker in range(need_workers):
        if not(handles[idx_worker]['received']):
            if handles[idx_worker]['async_handle'].ready():
                idx1 = handles[idx_worker]['y_line_n']
                idx2 = idx1+each_worker_does_y_lines
                data_in = pickle.loads(handles[idx_worker]['async_handle'].result)
                world.rxPlane.pressurefield[idx1:idx2, :] = data_in
                handles[idx_worker]['received'] = True
                items_to_receive -= 1
time_end = time.clock()
cueBeamCore2.do_plot_abs(world)
overall_performance = world.get_ray_count() / (time_end - time_start)
print("made {:06.2f} k rays/sec".format(overall_performance/1.0e3))