# from distutils.core import setup
import setuptools
from setuptools import setup

setup(
    name='cueBeam',
    version='0.0.1',
    packages=setuptools.find_packages(exclude=['contrib', 'docs', 'tests*']),
    url='http://www.strath.ac.uk',
    license='License :: Free for non-commercial use', # license='Creative commons BY-NC-SA',
    author='Jerzy Dziewierz',
    author_email='jerzy.dziewierz@strath.ac.uk',
    description='Simple acoustic interference and beam forming simulation based on Huygen''s principle',
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        #     # Indicate who your project is intended for
        'Intended Audience :: Manufacturing',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Physics',
        'Operating System :: OS Independent',
        'Topic :: Software Development :: Libraries',
        'Topic :: Utilities',
        'License :: Free for non-commercial use',
        'Programming Language :: Python :: 3 :: Only'
    ],
     keywords='cueBeam cueArt acoustics ultrasound NDE beamforming beam forming'
)

# from distutils.core import setup
# from setuptools import setup
# import setuptools
#
# setup(
#     name='LaudaDriver',
#     version='1.0.0a7',
#     packages=setuptools.find_packages(exclude=['contrib', 'docs', 'tests*']),  # ['LaudaDriver'],
#     url='https://uni.lardner.io/Jurek/Lauda_driver',
#     license='License :: OSI Approved :: Universal Permissive License (UPL)',
#     author='Jerzy Dziewierz',
#     author_email='jerzy.dziewierz@strath.ac.uk',
#     description='Simple helper to control the Lauda thermostat',
#     long_description="""
# Lauda Eco Silver thermostat helper library
# ==========================================
#
# A small helper module to control the Lauda-Brinkman ECO Silver thermostat,
# commonly found in the University of Strathclyde's CMAC Labs
#
# Usage is fairly obvious:
#
#
#     import lauda
#
#     lauda.set_pumping_speed(lauda.PumpingSpeed.a_Typical); lauda.set_temp(25.0); lauda.start()
#
#     current_temp_readout = lauda.read_current_temp()
#
#
# Note that the COM port used is hard-coded to "COM9".
# At this time this can only be changed by editing the source code.
# """,
#
#     classifiers=[
#     # How mature is this project? Common values are
#     #   3 - Alpha
#     #   4 - Beta
#     #   5 - Production/Stable
#     'Development Status :: 3 - Alpha',
#     # Indicate who your project is intended for
#     'Intended Audience :: Manufacturing',
#     'Intended Audience :: Science/Research',
#     'Intended Audience :: Developers',
#     'Natural Language :: English',
#     'Topic :: Scientific/Engineering :: Chemistry',
#     'Operating System :: OS Independent',
#     'Topic :: Software Development :: Libraries',
#     'Topic :: Utilities',
#     'License :: OSI Approved :: Universal Permissive License (UPL)',
#
#     # Specify the Python versions you support here. In particular, ensure
#     # that you indicate whether you support Python 2, Python 3 or both.
#     'Programming Language :: Python :: 3 :: Only'],
#     keywords='Lauda thermostat driver helper',
#     install_requires=['pyserial']
# )
