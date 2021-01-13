import os
import sys


class OSEnv:
    """docstring for OSEnv"""

    def __init__(self, **kwargs):
        self.kwargs = {k.upper(): v for k, v in kwargs.items()}
        # backup and recovery
        self.paths_backup = {}
        self.syspath_backup = []

    def __enter__(self):
        for env, val in self.kwargs.items():
            origin_env = os.environ.get(env, "")
            # backup setattr(self, env, origin_env)
            self.paths_backup[env] = origin_env
            origin_env = origin_env.strip(":")
            val = val.strip(":")
            if env == "PYTHONPATH":
                pp = val.split(':')
                self.syspath_backup = sys.path.copy()
                sys.path = pp + sys.path
            os.environ[env] = ':'.join([val, origin_env]).strip(":")
            print("Current {} is {}".format(env, os.environ.get(env, "")))
        return self

    def __exit__(self, type, value, trace):
        for env, val in self.kwargs.items():
            if env == "PYTHONPATH":
                sys.path = self.syspath_backup
            origin_env = self.paths_backup.get(env)
            if origin_env:
                os.environ[env] = origin_env
            else:
                del os.environ[env]
        if trace:
            print(trace)
