from distutils.core import setup, Extension
import subprocess
module=Extension('wlc', ['wlc_python.c'])
module.libraries = ['wlc']
# module.extra_compile_args=[subprocess.check_output(["wlc-config", "--libs"]),
                           # subprocess.check_output(["wlc-config", "--cflags"])];
module.extra_compile_args=[subprocess.check_output(["wlc-config", "--cflags"])]
module.extra_link_args=[subprocess.check_output(["wlc-config", "--libs"])]
setup(name='wlc', version='0.0', ext_modules=[module])
