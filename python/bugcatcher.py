# notify bugsnag that i am imported
import bugsnag
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