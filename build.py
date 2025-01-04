import subprocess
from hatchling.builders.hooks.plugin.interface import BuildHookInterface


class CustomBuildHook(BuildHookInterface):
    def initialize(self, version, build_data):
        # Run setup.py build_ext before main build
        subprocess.check_call(["python", "setup.py", "build_ext", "--inplace"])
        return super().initialize(version, build_data)
