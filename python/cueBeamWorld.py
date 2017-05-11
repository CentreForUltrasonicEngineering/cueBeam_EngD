
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # CueBeamWorld class
# # Utility class. Usefull for python-side works and debugging, but do not actually use for production
# # because it is difficult to modify it correctly from Matlab side.
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


def example_plot(the_input):
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
