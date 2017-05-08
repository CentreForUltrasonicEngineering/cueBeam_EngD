import time
import sys
import cueBeamCore2
w = cueBeamCore2.CueBeamWorld()
w.rxPlane.set_nxnynz(1,64,160)
q=cueBeamCore2.send_and_receive(w)

cueBeamCore2.do_plot_abs(q)
# now, the idea is to split the workload into units that would execute in approximately given time




