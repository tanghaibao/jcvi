import subprocess

from hatchling.builders.hooks.plugin.interface import BuildHookInterface
from setuptools.command.bdist_wheel import get_abi_tag, get_platform, tags


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
        impl_name = tags.interpreter_name()
        impl_ver = tags.interpreter_version()
        abi_tag = get_abi_tag()
        plat_tag = get_platform(None)
        plat_tag = (  # macosx_11.0 => macosx_11_0
            plat_tag.lower().replace("-", "_").replace(".", "_").replace(" ", "_")
        )
        return f"{impl_name}{impl_ver}-{abi_tag}-{plat_tag}"
