import time
import sys
import cueBeamCore2
import dill
import pickle
dill.settings['recurse'] = True

w = cueBeamCore2.CueBeamWorld()
w.rxPlane.set_nxnynz(1,32,160)
w.rxPlane.set_nxnynz(1,33,160)

local_result = cueBeamCore2.beamsim.delay(pickle.dumps(w))
#
# async_handle = cueBeamCore2.beamsim.delay(dill.dumps(w))
#
# while not(async_handle.ready()):
#     sys.stdout.write('.')
#     time.sleep(0.5)
#     pass
# sys.stdout.write('\n')
# time.sleep(0.5)
# remote_result = dill.loads(async_handle.result)
# cueBeamCore2.do_plot_abs(remote_result)