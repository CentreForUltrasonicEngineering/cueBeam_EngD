import math
import numpy
import matplotlib.pyplot as plt
import time
import copy


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

        def __init__(self, x: float=0.0, y: float=0.0, z: float=0.0, amplitude: float=1.0, phase: float=0.0):
            self.x = x
            self.y = y
            self.z = z
            self.amplitude = amplitude
            self.phase = phase

        def setXYZ(self, x: float=0.0, y: float=0.0, z: float=0.0):
            self.x = x
            self.y = y
            self.z = z

        def setAmpPhase(self,amplitude: float=1.0, phase: float=0.0):
            self.amplitude = amplitude
            self.phase = phase

    elements = [TxElement(0.0,  3.0e-3, -15.0e-3, 1.0, 0.0),
                TxElement(0.0,  2.0e-3, -10.0e-3, 1.0, 0.0),
                TxElement(0.0,  1.0e-3,  -5.0e-3, 1.0, 0.0),
                TxElement(0.0,  0.0e-3,   0.0e-3, 1.0, 0.0),
                TxElement(0.0, -1.0e-3,   5.0e-3, 1.0, 0.0),
                TxElement(0.0, -2.0e-3,  10.0e-3, 1.0, 0.0),
                TxElement(0.0, -3.0e-3,  15.0e-3, 1.0, 0.0)
                ]
    """List of monopoles.
    
    List of TxElement() - the list of monopoles emitting energy into the medium"""

    class RxPlane:
        """describes a rectangular XZ section of the plane - and points of sampling(calculation) of the pressure"""
        x0 = -0.01e-3  # needs to be a little bit off plane to avoid division by zer
        "x-origin of the measurement plane"
        y0 =  1.0e-3
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

        def __init__(self, x0=-0.01e-3, y0=1.0e-3, z0=0.0e-3, dx=0.5e-3, dy=0.5e-3, dz=0.5e-3, nx=1, ny=240, nz=160):
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

    def beamsim_timed(self):
        t1 = time.clock()
        self.beamsim()
        t2 = time.clock()
        dt = t2 - t1
        rays_per_second = self.get_ray_count() / dt
        print("got {:0.2f} M rays per second".format(rays_per_second/1e6))
        return rays_per_second


    def beamsim_simpler(self,
                        k: float=1000.0,
                        x0: float=0.0,
                        y0: float = 0.0,
                        z0: float = 0.0,
                        nx: int = 1,
                        ny: int = 240,
                        nz: int = 160,
                        dx: float = 1.0,
                        dy: float = 1.0e-3,
                        dz: float = 1.0e-3,
                        elements_vectorized=None):
        """wrapper for other methods so that this method returns a good result in a single call"""
        if elements_vectorized is None:
            elements_vectorized = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]

        self.wavenumber = k
        self.rxPlane.x0 = x0
        self.rxPlane.y0 = y0
        self.rxPlane.z0 = z0
        self.rxPlane.nx = nx
        self.rxPlane.ny = ny
        self.rxPlane.nz = nz
        self.rxPlane.dx = dx
        self.rxPlane.dy = dy
        self.rxPlane.dz = dz
        # now, the format for elements_vectorized is mapped after matlab's original version:
        # it is 1D, 1x 6n list where n=number of elements,
        # and the items are float (x,y,z,amplitude, phase, reserved)
        assert(math.modf(len(elements_vectorized))[0] == 0, "there must be n*6 in the vector")
        element_count = len(elements_vectorized) / 6
        self.elements.clear()
        for idx_element in range(element_count):
            element_pointer = 6*idx_element
            new_element = CueBeamSolver.TxElement(
                        x=elements_vectorized[element_pointer + 0],
                        y=elements_vectorized[element_pointer + 1],
                        z=elements_vectorized[element_pointer + 2],
                        amplitude=elements_vectorized[element_pointer + 3],
                        phase=elements_vectorized[element_pointer + 4]
                         )
            self.elements.append(new_element)

        return self.elements



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
        local_rx_plane = self.rxPlane
        local_elements = self.elements
        local_wavenumber = self.wavenumber

        ix = 0
        for iz in range(local_rx_plane.nz):
            for iy in range(local_rx_plane.ny):
                # initialize accumulator
                pressure_re = 0.0
                pressure_im = 0.0
                # calculate the XYZ location of the transmitter
                pixel_x = local_rx_plane.x0 + ix * local_rx_plane.dx
                pixel_y = local_rx_plane.y0 + iy * local_rx_plane.dy
                pixel_z = local_rx_plane.z0 + iz * local_rx_plane.dz
                # for each element do:
                for itx in range(len(local_elements)):
                    # calculate distance between the transmitter element and receiver pixel
                    ddx = (pixel_x - local_elements[itx].x)
                    ddy = (pixel_y - local_elements[itx].y)
                    ddz = (pixel_z - local_elements[itx].z)
                    distance = math.sqrt(ddx * ddx + ddy * ddy + ddz * ddz)
                    # based on distance and wavenumber, calculate the amplitude and phase
                    k_phase = - local_wavenumber * distance + local_elements[itx].phase
                    k_amplitude = local_elements[itx].amplitude / (2 * 3.14159 * distance)
                    # increment the accumulator
                    pressure_re += math.cos(k_phase) * k_amplitude
                    pressure_im += math.sin(k_phase) * k_amplitude
                # memory write
                local_rx_plane.pressurefield[(iy, iz)] = pressure_re + 1.0j * pressure_im

                # end for iy
            # end for iz
        # store result to the object
        self.rxPlane = copy.copy(local_rx_plane)
        # end beamsim


    def get_ray_count(self):
        return self.rxPlane.nx * self.rxPlane.ny * self.rxPlane.nz * len(self.elements)

    def do_plot_real(self):
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
            clim=(-8.0, 8.0),
            origin="lower")
        # end imshow
        plt.set_cmap("PRGn")  # purple-white-green color map
        plt.xlabel("z-axis[m]")
        plt.ylabel("y-axis[m]")
        plt.show()

    def do_plot_abs(self):
        """Makes an example plot.

         It is configured for the default data. Compatible with Jupyter"""

        hfig = plt.figure(num=1, figsize=(8, 6), dpi=90, facecolor='white', edgecolor='black')
        imgplot = plt.imshow(
            X=numpy.abs(self.rxPlane.pressurefield),
            extent=(
                self.rxPlane.z0, self.rxPlane.z0 + self.rxPlane.nz * self.rxPlane.dz,
                self.rxPlane.y0 + self.rxPlane.ny * self.rxPlane.dy, self.rxPlane.y0
            ),
            interpolation="spline36",
            clim=(0, 8.0),
            origin="lower")
        # end imshow
        plt.set_cmap("plasma") # black-to-yellow color map
        plt.xlabel("z-axis[m]")
        plt.ylabel("y-axis[m]")
        plt.show()
