import subprocess
import sysconfig

from hatchling.builders.hooks.plugin.interface import BuildHookInterface


class CustomBuildHook(BuildHookInterface):
    def initialize(self, version, build_data):
        # Set the custom wheel tag
        build_data["tag"] = self._get_wheel_tag()
        print(f"Using custom wheel tag: {build_data['tag']}")
        # Run setup.py build_ext before main build
        subprocess.check_call(["python", "setup.py", "build_ext", "--inplace"])
        return super().initialize(version, build_data)

    def _get_wheel_tag(self):
        # Without the tag, the wheel will be named jcvi-0.0.0-py3-none-any.whl
        platform_tag = sysconfig.get_platform().replace("-", "_").replace(".", "_")
        python_version = sysconfig.get_python_version().replace(".", "")  # e.g., "310"
        python_impl = "cp"  # Assuming CPython. Modify if using PyPy or others.
        abi_tag = f"{python_impl}{python_version}"
        return f"{python_impl}{python_version}-{abi_tag}-{platform_tag}"
