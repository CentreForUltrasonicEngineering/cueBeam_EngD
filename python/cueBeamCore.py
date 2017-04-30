import math
import numpy
import matplotlib.pyplot as plt


class CueBeamSolver:
    """ 
    This is the core of the cueBeam.
    
    **Essential fields that You need to set:**
    
    wavenumber: float
        The enviroment property, note that wavenumber = 1/wavelength in that enviroment
    elements: list of TxElement
        description of the transmitters
    rxPlane: single item of class RxPlane()
        describes the regular sampling plane
    
    
    
    
    """

    wavenumber = 1.0 / 1.0e-3  # wavelength = 1/wavenumber; so wavenumber = 1/wavelength
    """Property of the environment: the wavenumber."""

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

        def __init__(self, x:float=0.0, y:float=0.0, z:float=0.0, amplitude:float=1.0, phase:float=0.0):
            self.x = x
            self.y = y
            self.z = z
            self.amplitude = amplitude
            self.phase = phase

    elements = [TxElement(0.0, 0.0, -15.0e-3, 1.0, 0.0),
                TxElement(0.0, 0.0, -10.0e-3, 1.0, 0.0),
                TxElement(0.0, 0.0,  -5.0e-3, 1.0, 0.0),
                TxElement(0.0, 0.0,   0.0e-3, 1.0, 0.0),
                TxElement(0.0, 0.0,   5.0e-3, 1.0, 0.0),
                TxElement(0.0, 0.0,  10.0e-3, 1.0, 0.0),
                TxElement(0.0, 0.0,  15.0e-3, 1.0, 0.0)
                ]
    """List of monopoles.
    
    List of TxElement() - the list of monopoles emitting energy into the medium"""

    class RxPlane:
        """describes a rectangular XZ section of the plane - and points of sampling(calculation) of the pressure"""
        x0 = -0.01e-3  # needs to be a little bit off plane to avoid division by zer
        "x-origin of the measurement plane"
        y0 = -1.0e-3
        "y-origin of the measurement plane"
        z0 = +0.0e-3
        "z-origin of the measurement plane"
        dx = +1.0e-3
        "x-step between the points in the measurement plane"
        dy = +1.0e-3
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

        def __init__(self, x0=-5.0e-3, y0=-5.0e-3, z0=0.0, dx=0.1e-3, dy=0.1e-3, dz=0.1e-3, nx=1, ny=10, nz=10):
            "Creates the new instance of RxPlane"
            self.x0 = x0
            self.y0 = y0
            self.z0 = z0
            self.dx = dx
            self.dy = dy
            self.dz = dz
            if nx != 1:
                raise AttributeError("nx must be 1: ")
            if (ny > 4096) | (ny < 1):
                raise AttributeError("ny must have a sensible value between 1 an 4096")
            self.ny = ny
            if (nz > 4096) | (nz < 1):
                raise AttributeError("nz must have a sensible value between 1 and 4096")
            self.nz = nz
            self.clear_pressurefield()

        def verify_plane_endpoints(self):
            "returns x,y,z location of the corner oposite to the origin corner"
            return [self.x0 + self.nx * self.dx, self.y0 + self.ny * self.dy, self.z0 + self.nz * self.dz]

        def clear_pressurefield(self):
            "empties the result"
            self.pressurefield = numpy.zeros((self.ny, self.nz), numpy.complex)

    rxPlane = RxPlane()
    "contains the description/location of the measurement points"

    def beamsim(self):
        """Execute the pressure calculation.
        
        Executes the calculation using settings stored in elements, rxPlane and wavenumber.
        
        Parameters
        ----------
        None. Takes input from the elements, rxPlane and wavenumber.
        
        Returns
        -------
        None. The result is stored in .pressurefield field
        
        
        """

        # note: there is only a single plane implemented. This is for compatibility with the old Matlab-Cuda version.
        # Essentially, there should be no problem in extending this to 3D computation as needed.
        ix = 0
        for iz in range(self.rxPlane.nz):
            for iy in range(self.rxPlane.ny):
                pressure_re = 0.0
                pressure_im = 0.0
                pixel_x = self.rxPlane.x0 + ix * self.rxPlane.dx
                pixel_y = self.rxPlane.y0 + iy * self.rxPlane.dy
                pixel_z = self.rxPlane.z0 + iz * self.rxPlane.dz
                for itx in range(len(self.elements)):
                    ddx = (pixel_x - self.elements[itx].x)
                    ddy = (pixel_y - self.elements[itx].y)
                    ddz = (pixel_z - self.elements[itx].z)
                    distance = math.sqrt(ddx * ddx + ddy * ddy + ddz * ddz)
                    kphase = - self.wavenumber * distance + self.elements[itx].phase
                    kamplitude = self.elements[itx].amplitude / (2 * 3.14159 * distance)
                    pressure_re += math.cos(kphase) * kamplitude
                    pressure_im += math.sin(kphase) * kamplitude
                # mem write
                self.rxPlane.pressurefield[(iy, iz)] = pressure_re + 1.0j * pressure_im

                # end for iy
                # end for iz
                # end beamsim

    def do_plot(self):
        """Makes an example plot.
        
         It is configured for the default data. Compatible with Jupyter"""

        hfig = plt.figure(num=1, figsize=(8, 6), dpi=90, facecolor='white', edgecolor='black')
        imgplot = plt.imshow(
            X=numpy.real(self.rxPlane.pressurefield),
            extent=(
                self.rxPlane.z0, self.rxPlane.z0 + self.rxPlane.nz * self.rxPlane.dz,
                self.rxPlane.y0 + self.rxPlane.ny * self.rxPlane.dy, self.rxPlane.y0
            ),
            interpolation="spline36",
            clim=(-7.0, 7.0),
            origin="lower")
        # end imshow
        plt.set_cmap("PRGn")
        plt.xlabel("z-axis[m]")
        plt.ylabel("y-axis[m]")
        plt.show()

