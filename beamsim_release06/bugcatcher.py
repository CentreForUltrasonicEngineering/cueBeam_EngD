# notify bugsnag that i am imported
def install_and_import(package):
    import importlib
    print('starting to import {}'.format(package))
    try:
        print('importing')
        importlib.import_module(package)
        print('import succeeded');
    except ImportError:
        print('in ImportError');
        # import pip
        # print('trying to install the package')
        # pip.main(['install', package])
        import os
        os.system('pip install bugsnag')
        print('')
        print('apparently installed')
    finally:
        print('refreshing site')
        import site
        importlib.import_module('site')
        importlib.reload(site)
        print('importing...')
        globals()[package] = importlib.import_module(package)
        print('done.')


install_and_import('bugsnag')

from bugsnag.handlers import BugsnagHandler

import logging
import os
import platform
import socket

bugsnag.configure(
    api_key="870c14af39685068c999a54d623454bf",
    project_root="T:\git\cueBeam\python",
    app_version="0.0.1",
    ignore_classes=["django.http.Http404"],
    release_stage="development",
    asynchronous=False) # trying to disable async - possibly matlab disbands the instance before the error is reported

logger = logging.getLogger("basic")
logger.setLevel(logging.INFO)
logger.addHandler(BugsnagHandler())

bugsnag.notify(Exception("I have been imported"), context="bugcatcher.start", meta_data={"startup_data":{"login": os.getlogin(), "hostname": socket.gethostname(), "node":platform.node()}})

# raise ValueError('trying to catch something from Matlab')

# zepta = 3/0

# print("this is the end - should not be executed {}".format(zepta))