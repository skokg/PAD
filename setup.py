"""Build the platform-specific ctypes library with setuptools."""

import sys

from setuptools import Extension, setup
from setuptools.command.bdist_wheel import bdist_wheel
from setuptools.command.build_ext import build_ext
from setuptools.errors import PlatformError


class BuildNativeLibrary(build_ext):
    """Keep the shared-library filename independent of the Python ABI."""

    def run(self):
        if sys.platform not in ("linux", "win32", "darwin"):
            raise PlatformError("pad_python supports Linux, Windows and macOS.")
        super().run()

    def get_ext_filename(self, fullname):
        return {"linux": "libpad_python.so", "win32": "pad_python.dll", "darwin": "libpad_python.dylib"}[sys.platform]

    def get_export_symbols(self, ext):
        # The C++ source exports the C ABI; there is no Python PyInit function.
        return []

    def get_libraries(self, ext):
        # Never link this ctypes library to a version-specific Python runtime.
        return ext.libraries

    def build_extensions(self):
        if self.compiler.compiler_type == "msvc":
            compile_args = ["/O2", "/std:c++14", "/EHsc", "/W4"]
            link_args = []
        else:
            compile_args = ["-O3", "-std=c++11", "-Wall", "-Wextra", "-Wpedantic", "-pthread", "-fvisibility=hidden"]
            link_args = ["-pthread", "-Wl,--no-undefined"] if sys.platform == "linux" else ["-pthread"]
        if sys.platform == "darwin":
            # Build a dylib instead of Python's usual macOS extension bundle.
            for name in ("linker_so", "linker_so_cxx"):
                command = getattr(self.compiler, name, None)
                if command:
                    setattr(self.compiler, name, ["-dynamiclib" if part == "-bundle" else part for part in command if part not in ("-undefined", "dynamic_lookup")])
        for extension in self.extensions:
            extension.extra_compile_args = compile_args
            extension.extra_link_args = link_args
        super().build_extensions()


class BuildPlatformWheel(bdist_wheel):
    """Publish a separately built wheel for each selected CPython version."""

    def get_tag(self):
        _, _, platform = super().get_tag()
        # The Python selector is version-specific; the native ABI remains ctypes.
        return f"cp{sys.version_info.major}{sys.version_info.minor}", "none", platform


# Metadata and runtime dependencies are defined in pyproject.toml.
# This is a plain shared library loaded by ctypes, not imported as an extension.
setup(
    ext_modules=[
        Extension(
            "libpad_python",
            sources=["src/cpp/CC_PAD_python_lib.cc"],
            include_dirs=["src/cpp"],
            depends=["src/cpp/CU_PAD_code.cc", "src/cpp/CU_kdtree_with_index.cc"],
            language="c++",
        ),
    ],
    cmdclass={"build_ext": BuildNativeLibrary, "bdist_wheel": BuildPlatformWheel},
)
