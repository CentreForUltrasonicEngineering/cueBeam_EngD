import time
import sys
import cueBeamCore2
w = cueBeamCore2.CueBeamWorld()
# local_result = cueBeamCore2.beamsim(w)
async_handle = cueBeamCore2.beamsim.delay(w)

while not(async_handle.ready()):
    time.sleep(0.1)
    sys.stdout.write('.')
    pass
remote_result = async_handle.result

cueBeamCore2.do_plot_abs(remote_result)