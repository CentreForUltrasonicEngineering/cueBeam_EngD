syspy=py.importlib.import_module('sys');
syspy.path.append([pwd '//..//python']);
cueBeamPy=py.importlib.import_module('cueBeamCore3');
py.importlib.reload(cueBeamPy);