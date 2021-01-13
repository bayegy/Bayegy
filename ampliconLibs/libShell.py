from .osEnv import OSEnv
from .settings import path
import os


class PathNotConfigedException(Exception):
    pass


def solve_conda_env(config_key=False, env_path=False):
    home = env_path or path.get(config_key, "")
    if not home:
        raise PathNotConfigedException("{} not configed".format(config_key))
    envs = dict(path=os.path.join(home, 'bin'))
    lib_path = os.path.join(home, 'lib')
    # pythonpath = ""
    for d in os.listdir(lib_path):
        if d.startswith('python'):
            envs['pythonpath'] = os.path.join(lib_path, d, 'site-packages')
            break
    perl_lib_path = os.path.join(lib_path, 'site_perl')
    if os.path.exists(perl_lib_path):
        perl_v = os.listdir(perl_lib_path)
        if perl_v:
            envs['perl5lib'] = os.path.join(perl_lib_path, perl_v[0])
    r_libs = os.path.join(lib_path, 'R/library')
    if os.path.exists(r_libs):
        envs['r_libs'] = r_libs
    return envs


class qiime2(OSEnv):

    def __init__(self):
        super().__init__(**solve_conda_env('qiime2_home'))

class lefse(OSEnv):

    def __init__(self):
        paths = dict(
            pythonpath=path.get("lefse_pylib_home"),
            path=path.get("lefse_py_home"),
            r_libs=path.get("lefse_rlib_home")
        )
        for k, v in paths.items():
            if not v:
                raise PathNotConfigedException("{} not configed".format(k))
        super().__init__(**paths)

class picrust2(OSEnv):

    def __init__(self):
        super().__init__(**solve_conda_env('picrust2_home'))


class qiime1(OSEnv):

    def __init__(self):
        super().__init__(**solve_conda_env('qiime1_home'))
