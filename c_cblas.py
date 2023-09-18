r"""Wrapper for cblas.h

Generated with:
/home/james.bland/.conda/envs/py39/bin/ctypesgen -o c_cblas.py -llibcblas.so /usr/include/cblas.h

Do not modify this file.
"""

__docformat__ = "restructuredtext"

# Begin preamble for Python

import ctypes
import sys
from ctypes import *  # noqa: F401, F403

_int_types = (ctypes.c_int16, ctypes.c_int32)
if hasattr(ctypes, "c_int64"):
    # Some builds of ctypes apparently do not have ctypes.c_int64
    # defined; it's a pretty good bet that these builds do not
    # have 64-bit pointers.
    _int_types += (ctypes.c_int64,)
for t in _int_types:
    if ctypes.sizeof(t) == ctypes.sizeof(ctypes.c_size_t):
        c_ptrdiff_t = t
del t
del _int_types



class UserString:
    def __init__(self, seq):
        if isinstance(seq, bytes):
            self.data = seq
        elif isinstance(seq, UserString):
            self.data = seq.data[:]
        else:
            self.data = str(seq).encode()

    def __bytes__(self):
        return self.data

    def __str__(self):
        return self.data.decode()

    def __repr__(self):
        return repr(self.data)

    def __int__(self):
        return int(self.data.decode())

    def __long__(self):
        return int(self.data.decode())

    def __float__(self):
        return float(self.data.decode())

    def __complex__(self):
        return complex(self.data.decode())

    def __hash__(self):
        return hash(self.data)

    def __le__(self, string):
        if isinstance(string, UserString):
            return self.data <= string.data
        else:
            return self.data <= string

    def __lt__(self, string):
        if isinstance(string, UserString):
            return self.data < string.data
        else:
            return self.data < string

    def __ge__(self, string):
        if isinstance(string, UserString):
            return self.data >= string.data
        else:
            return self.data >= string

    def __gt__(self, string):
        if isinstance(string, UserString):
            return self.data > string.data
        else:
            return self.data > string

    def __eq__(self, string):
        if isinstance(string, UserString):
            return self.data == string.data
        else:
            return self.data == string

    def __ne__(self, string):
        if isinstance(string, UserString):
            return self.data != string.data
        else:
            return self.data != string

    def __contains__(self, char):
        return char in self.data

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        return self.__class__(self.data[index])

    def __getslice__(self, start, end):
        start = max(start, 0)
        end = max(end, 0)
        return self.__class__(self.data[start:end])

    def __add__(self, other):
        if isinstance(other, UserString):
            return self.__class__(self.data + other.data)
        elif isinstance(other, bytes):
            return self.__class__(self.data + other)
        else:
            return self.__class__(self.data + str(other).encode())

    def __radd__(self, other):
        if isinstance(other, bytes):
            return self.__class__(other + self.data)
        else:
            return self.__class__(str(other).encode() + self.data)

    def __mul__(self, n):
        return self.__class__(self.data * n)

    __rmul__ = __mul__

    def __mod__(self, args):
        return self.__class__(self.data % args)

    # the following methods are defined in alphabetical order:
    def capitalize(self):
        return self.__class__(self.data.capitalize())

    def center(self, width, *args):
        return self.__class__(self.data.center(width, *args))

    def count(self, sub, start=0, end=sys.maxsize):
        return self.data.count(sub, start, end)

    def decode(self, encoding=None, errors=None):  # XXX improve this?
        if encoding:
            if errors:
                return self.__class__(self.data.decode(encoding, errors))
            else:
                return self.__class__(self.data.decode(encoding))
        else:
            return self.__class__(self.data.decode())

    def encode(self, encoding=None, errors=None):  # XXX improve this?
        if encoding:
            if errors:
                return self.__class__(self.data.encode(encoding, errors))
            else:
                return self.__class__(self.data.encode(encoding))
        else:
            return self.__class__(self.data.encode())

    def endswith(self, suffix, start=0, end=sys.maxsize):
        return self.data.endswith(suffix, start, end)

    def expandtabs(self, tabsize=8):
        return self.__class__(self.data.expandtabs(tabsize))

    def find(self, sub, start=0, end=sys.maxsize):
        return self.data.find(sub, start, end)

    def index(self, sub, start=0, end=sys.maxsize):
        return self.data.index(sub, start, end)

    def isalpha(self):
        return self.data.isalpha()

    def isalnum(self):
        return self.data.isalnum()

    def isdecimal(self):
        return self.data.isdecimal()

    def isdigit(self):
        return self.data.isdigit()

    def islower(self):
        return self.data.islower()

    def isnumeric(self):
        return self.data.isnumeric()

    def isspace(self):
        return self.data.isspace()

    def istitle(self):
        return self.data.istitle()

    def isupper(self):
        return self.data.isupper()

    def join(self, seq):
        return self.data.join(seq)

    def ljust(self, width, *args):
        return self.__class__(self.data.ljust(width, *args))

    def lower(self):
        return self.__class__(self.data.lower())

    def lstrip(self, chars=None):
        return self.__class__(self.data.lstrip(chars))

    def partition(self, sep):
        return self.data.partition(sep)

    def replace(self, old, new, maxsplit=-1):
        return self.__class__(self.data.replace(old, new, maxsplit))

    def rfind(self, sub, start=0, end=sys.maxsize):
        return self.data.rfind(sub, start, end)

    def rindex(self, sub, start=0, end=sys.maxsize):
        return self.data.rindex(sub, start, end)

    def rjust(self, width, *args):
        return self.__class__(self.data.rjust(width, *args))

    def rpartition(self, sep):
        return self.data.rpartition(sep)

    def rstrip(self, chars=None):
        return self.__class__(self.data.rstrip(chars))

    def split(self, sep=None, maxsplit=-1):
        return self.data.split(sep, maxsplit)

    def rsplit(self, sep=None, maxsplit=-1):
        return self.data.rsplit(sep, maxsplit)

    def splitlines(self, keepends=0):
        return self.data.splitlines(keepends)

    def startswith(self, prefix, start=0, end=sys.maxsize):
        return self.data.startswith(prefix, start, end)

    def strip(self, chars=None):
        return self.__class__(self.data.strip(chars))

    def swapcase(self):
        return self.__class__(self.data.swapcase())

    def title(self):
        return self.__class__(self.data.title())

    def translate(self, *args):
        return self.__class__(self.data.translate(*args))

    def upper(self):
        return self.__class__(self.data.upper())

    def zfill(self, width):
        return self.__class__(self.data.zfill(width))


class MutableString(UserString):
    """mutable string objects

    Python strings are immutable objects.  This has the advantage, that
    strings may be used as dictionary keys.  If this property isn't needed
    and you insist on changing string values in place instead, you may cheat
    and use MutableString.

    But the purpose of this class is an educational one: to prevent
    people from inventing their own mutable string class derived
    from UserString and than forget thereby to remove (override) the
    __hash__ method inherited from UserString.  This would lead to
    errors that would be very hard to track down.

    A faster and better solution is to rewrite your program using lists."""

    def __init__(self, string=""):
        self.data = string

    def __hash__(self):
        raise TypeError("unhashable type (it is mutable)")

    def __setitem__(self, index, sub):
        if index < 0:
            index += len(self.data)
        if index < 0 or index >= len(self.data):
            raise IndexError
        self.data = self.data[:index] + sub + self.data[index + 1 :]

    def __delitem__(self, index):
        if index < 0:
            index += len(self.data)
        if index < 0 or index >= len(self.data):
            raise IndexError
        self.data = self.data[:index] + self.data[index + 1 :]

    def __setslice__(self, start, end, sub):
        start = max(start, 0)
        end = max(end, 0)
        if isinstance(sub, UserString):
            self.data = self.data[:start] + sub.data + self.data[end:]
        elif isinstance(sub, bytes):
            self.data = self.data[:start] + sub + self.data[end:]
        else:
            self.data = self.data[:start] + str(sub).encode() + self.data[end:]

    def __delslice__(self, start, end):
        start = max(start, 0)
        end = max(end, 0)
        self.data = self.data[:start] + self.data[end:]

    def immutable(self):
        return UserString(self.data)

    def __iadd__(self, other):
        if isinstance(other, UserString):
            self.data += other.data
        elif isinstance(other, bytes):
            self.data += other
        else:
            self.data += str(other).encode()
        return self

    def __imul__(self, n):
        self.data *= n
        return self


class String(MutableString, ctypes.Union):

    _fields_ = [("raw", ctypes.POINTER(ctypes.c_char)), ("data", ctypes.c_char_p)]

    def __init__(self, obj=b""):
        if isinstance(obj, (bytes, UserString)):
            self.data = bytes(obj)
        else:
            self.raw = obj

    def __len__(self):
        return self.data and len(self.data) or 0

    def from_param(cls, obj):
        # Convert None or 0
        if obj is None or obj == 0:
            return cls(ctypes.POINTER(ctypes.c_char)())

        # Convert from String
        elif isinstance(obj, String):
            return obj

        # Convert from bytes
        elif isinstance(obj, bytes):
            return cls(obj)

        # Convert from str
        elif isinstance(obj, str):
            return cls(obj.encode())

        # Convert from c_char_p
        elif isinstance(obj, ctypes.c_char_p):
            return obj

        # Convert from POINTER(ctypes.c_char)
        elif isinstance(obj, ctypes.POINTER(ctypes.c_char)):
            return obj

        # Convert from raw pointer
        elif isinstance(obj, int):
            return cls(ctypes.cast(obj, ctypes.POINTER(ctypes.c_char)))

        # Convert from ctypes.c_char array
        elif isinstance(obj, ctypes.c_char * len(obj)):
            return obj

        # Convert from object
        else:
            return String.from_param(obj._as_parameter_)

    from_param = classmethod(from_param)


def ReturnString(obj, func=None, arguments=None):
    return String.from_param(obj)


# As of ctypes 1.0, ctypes does not support custom error-checking
# functions on callbacks, nor does it support custom datatypes on
# callbacks, so we must ensure that all callbacks return
# primitive datatypes.
#
# Non-primitive return values wrapped with UNCHECKED won't be
# typechecked, and will be converted to ctypes.c_void_p.
def UNCHECKED(type):
    if hasattr(type, "_type_") and isinstance(type._type_, str) and type._type_ != "P":
        return type
    else:
        return ctypes.c_void_p


# ctypes doesn't have direct support for variadic functions, so we have to write
# our own wrapper class
class _variadic_function(object):
    def __init__(self, func, restype, argtypes, errcheck):
        self.func = func
        self.func.restype = restype
        self.argtypes = argtypes
        if errcheck:
            self.func.errcheck = errcheck

    def _as_parameter_(self):
        # So we can pass this variadic function as a function pointer
        return self.func

    def __call__(self, *args):
        fixed_args = []
        i = 0
        for argtype in self.argtypes:
            # Typecheck what we can
            fixed_args.append(argtype.from_param(args[i]))
            i += 1
        return self.func(*fixed_args + list(args[i:]))


def ord_if_char(value):
    """
    Simple helper used for casts to simple builtin types:  if the argument is a
    string type, it will be converted to it's ordinal value.

    This function will raise an exception if the argument is string with more
    than one characters.
    """
    return ord(value) if (isinstance(value, bytes) or isinstance(value, str)) else value

# End preamble

_libs = {}
_libdirs = []

# Begin loader

"""
Load libraries - appropriately for all our supported platforms
"""
# ----------------------------------------------------------------------------
# Copyright (c) 2008 David James
# Copyright (c) 2006-2008 Alex Holkner
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.
#  * Neither the name of pyglet nor the names of its
#    contributors may be used to endorse or promote products
#    derived from this software without specific prior written
#    permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
# ----------------------------------------------------------------------------

import ctypes
import ctypes.util
import glob
import os.path
import platform
import re
import sys


def _environ_path(name):
    """Split an environment variable into a path-like list elements"""
    if name in os.environ:
        return os.environ[name].split(":")
    return []


class LibraryLoader:
    """
    A base class For loading of libraries ;-)
    Subclasses load libraries for specific platforms.
    """

    # library names formatted specifically for platforms
    name_formats = ["%s"]

    class Lookup:
        """Looking up calling conventions for a platform"""

        mode = ctypes.DEFAULT_MODE

        def __init__(self, path):
            super(LibraryLoader.Lookup, self).__init__()
            self.access = dict(cdecl=ctypes.CDLL(path, self.mode))

        def get(self, name, calling_convention="cdecl"):
            """Return the given name according to the selected calling convention"""
            if calling_convention not in self.access:
                raise LookupError(
                    "Unknown calling convention '{}' for function '{}'".format(
                        calling_convention, name
                    )
                )
            return getattr(self.access[calling_convention], name)

        def has(self, name, calling_convention="cdecl"):
            """Return True if this given calling convention finds the given 'name'"""
            if calling_convention not in self.access:
                return False
            return hasattr(self.access[calling_convention], name)

        def __getattr__(self, name):
            return getattr(self.access["cdecl"], name)

    def __init__(self):
        self.other_dirs = []

    def __call__(self, libname):
        """Given the name of a library, load it."""
        paths = self.getpaths(libname)

        for path in paths:
            # noinspection PyBroadException
            try:
                return self.Lookup(path)
            except Exception:  # pylint: disable=broad-except
                pass

        raise ImportError("Could not load %s." % libname)

    def getpaths(self, libname):
        """Return a list of paths where the library might be found."""
        if os.path.isabs(libname):
            yield libname
        else:
            # search through a prioritized series of locations for the library

            # we first search any specific directories identified by user
            for dir_i in self.other_dirs:
                for fmt in self.name_formats:
                    # dir_i should be absolute already
                    yield os.path.join(dir_i, fmt % libname)

            # check if this code is even stored in a physical file
            try:
                this_file = __file__
            except NameError:
                this_file = None

            # then we search the directory where the generated python interface is stored
            if this_file is not None:
                for fmt in self.name_formats:
                    yield os.path.abspath(os.path.join(os.path.dirname(__file__), fmt % libname))

            # now, use the ctypes tools to try to find the library
            for fmt in self.name_formats:
                path = ctypes.util.find_library(fmt % libname)
                if path:
                    yield path

            # then we search all paths identified as platform-specific lib paths
            for path in self.getplatformpaths(libname):
                yield path

            # Finally, we'll try the users current working directory
            for fmt in self.name_formats:
                yield os.path.abspath(os.path.join(os.path.curdir, fmt % libname))

    def getplatformpaths(self, _libname):  # pylint: disable=no-self-use
        """Return all the library paths available in this platform"""
        return []


# Darwin (Mac OS X)


class DarwinLibraryLoader(LibraryLoader):
    """Library loader for MacOS"""

    name_formats = [
        "lib%s.dylib",
        "lib%s.so",
        "lib%s.bundle",
        "%s.dylib",
        "%s.so",
        "%s.bundle",
        "%s",
    ]

    class Lookup(LibraryLoader.Lookup):
        """
        Looking up library files for this platform (Darwin aka MacOS)
        """

        # Darwin requires dlopen to be called with mode RTLD_GLOBAL instead
        # of the default RTLD_LOCAL.  Without this, you end up with
        # libraries not being loadable, resulting in "Symbol not found"
        # errors
        mode = ctypes.RTLD_GLOBAL

    def getplatformpaths(self, libname):
        if os.path.pathsep in libname:
            names = [libname]
        else:
            names = [fmt % libname for fmt in self.name_formats]

        for directory in self.getdirs(libname):
            for name in names:
                yield os.path.join(directory, name)

    @staticmethod
    def getdirs(libname):
        """Implements the dylib search as specified in Apple documentation:

        http://developer.apple.com/documentation/DeveloperTools/Conceptual/
            DynamicLibraries/Articles/DynamicLibraryUsageGuidelines.html

        Before commencing the standard search, the method first checks
        the bundle's ``Frameworks`` directory if the application is running
        within a bundle (OS X .app).
        """

        dyld_fallback_library_path = _environ_path("DYLD_FALLBACK_LIBRARY_PATH")
        if not dyld_fallback_library_path:
            dyld_fallback_library_path = [
                os.path.expanduser("~/lib"),
                "/usr/local/lib",
                "/usr/lib",
            ]

        dirs = []

        if "/" in libname:
            dirs.extend(_environ_path("DYLD_LIBRARY_PATH"))
        else:
            dirs.extend(_environ_path("LD_LIBRARY_PATH"))
            dirs.extend(_environ_path("DYLD_LIBRARY_PATH"))
            dirs.extend(_environ_path("LD_RUN_PATH"))

        if hasattr(sys, "frozen") and getattr(sys, "frozen") == "macosx_app":
            dirs.append(os.path.join(os.environ["RESOURCEPATH"], "..", "Frameworks"))

        dirs.extend(dyld_fallback_library_path)

        return dirs


# Posix


class PosixLibraryLoader(LibraryLoader):
    """Library loader for POSIX-like systems (including Linux)"""

    _ld_so_cache = None

    _include = re.compile(r"^\s*include\s+(?P<pattern>.*)")

    name_formats = ["lib%s.so", "%s.so", "%s"]

    class _Directories(dict):
        """Deal with directories"""

        def __init__(self):
            dict.__init__(self)
            self.order = 0

        def add(self, directory):
            """Add a directory to our current set of directories"""
            if len(directory) > 1:
                directory = directory.rstrip(os.path.sep)
            # only adds and updates order if exists and not already in set
            if not os.path.exists(directory):
                return
            order = self.setdefault(directory, self.order)
            if order == self.order:
                self.order += 1

        def extend(self, directories):
            """Add a list of directories to our set"""
            for a_dir in directories:
                self.add(a_dir)

        def ordered(self):
            """Sort the list of directories"""
            return (i[0] for i in sorted(self.items(), key=lambda d: d[1]))

    def _get_ld_so_conf_dirs(self, conf, dirs):
        """
        Recursive function to help parse all ld.so.conf files, including proper
        handling of the `include` directive.
        """

        try:
            with open(conf) as fileobj:
                for dirname in fileobj:
                    dirname = dirname.strip()
                    if not dirname:
                        continue

                    match = self._include.match(dirname)
                    if not match:
                        dirs.add(dirname)
                    else:
                        for dir2 in glob.glob(match.group("pattern")):
                            self._get_ld_so_conf_dirs(dir2, dirs)
        except IOError:
            pass

    def _create_ld_so_cache(self):
        # Recreate search path followed by ld.so.  This is going to be
        # slow to build, and incorrect (ld.so uses ld.so.cache, which may
        # not be up-to-date).  Used only as fallback for distros without
        # /sbin/ldconfig.
        #
        # We assume the DT_RPATH and DT_RUNPATH binary sections are omitted.

        directories = self._Directories()
        for name in (
            "LD_LIBRARY_PATH",
            "SHLIB_PATH",  # HP-UX
            "LIBPATH",  # OS/2, AIX
            "LIBRARY_PATH",  # BE/OS
        ):
            if name in os.environ:
                directories.extend(os.environ[name].split(os.pathsep))

        self._get_ld_so_conf_dirs("/etc/ld.so.conf", directories)

        bitage = platform.architecture()[0]

        unix_lib_dirs_list = []
        if bitage.startswith("64"):
            # prefer 64 bit if that is our arch
            unix_lib_dirs_list += ["/lib64", "/usr/lib64"]

        # must include standard libs, since those paths are also used by 64 bit
        # installs
        unix_lib_dirs_list += ["/lib", "/usr/lib"]
        if sys.platform.startswith("linux"):
            # Try and support multiarch work in Ubuntu
            # https://wiki.ubuntu.com/MultiarchSpec
            if bitage.startswith("32"):
                # Assume Intel/AMD x86 compat
                unix_lib_dirs_list += ["/lib/i386-linux-gnu", "/usr/lib/i386-linux-gnu"]
            elif bitage.startswith("64"):
                # Assume Intel/AMD x86 compatible
                unix_lib_dirs_list += [
                    "/lib/x86_64-linux-gnu",
                    "/usr/lib/x86_64-linux-gnu",
                ]
            else:
                # guess...
                unix_lib_dirs_list += glob.glob("/lib/*linux-gnu")
        directories.extend(unix_lib_dirs_list)

        cache = {}
        lib_re = re.compile(r"lib(.*)\.s[ol]")
        # ext_re = re.compile(r"\.s[ol]$")
        for our_dir in directories.ordered():
            try:
                for path in glob.glob("%s/*.s[ol]*" % our_dir):
                    file = os.path.basename(path)

                    # Index by filename
                    cache_i = cache.setdefault(file, set())
                    cache_i.add(path)

                    # Index by library name
                    match = lib_re.match(file)
                    if match:
                        library = match.group(1)
                        cache_i = cache.setdefault(library, set())
                        cache_i.add(path)
            except OSError:
                pass

        self._ld_so_cache = cache

    def getplatformpaths(self, libname):
        if self._ld_so_cache is None:
            self._create_ld_so_cache()

        result = self._ld_so_cache.get(libname, set())
        for i in result:
            # we iterate through all found paths for library, since we may have
            # actually found multiple architectures or other library types that
            # may not load
            yield i


# Windows


class WindowsLibraryLoader(LibraryLoader):
    """Library loader for Microsoft Windows"""

    name_formats = ["%s.dll", "lib%s.dll", "%slib.dll", "%s"]

    class Lookup(LibraryLoader.Lookup):
        """Lookup class for Windows libraries..."""

        def __init__(self, path):
            super(WindowsLibraryLoader.Lookup, self).__init__(path)
            self.access["stdcall"] = ctypes.windll.LoadLibrary(path)


# Platform switching

# If your value of sys.platform does not appear in this dict, please contact
# the Ctypesgen maintainers.

loaderclass = {
    "darwin": DarwinLibraryLoader,
    "cygwin": WindowsLibraryLoader,
    "win32": WindowsLibraryLoader,
    "msys": WindowsLibraryLoader,
}

load_library = loaderclass.get(sys.platform, PosixLibraryLoader)()


def add_library_search_dirs(other_dirs):
    """
    Add libraries to search paths.
    If library paths are relative, convert them to absolute with respect to this
    file's directory
    """
    for path in other_dirs:
        if not os.path.isabs(path):
            path = os.path.abspath(path)
        load_library.other_dirs.append(path)


del loaderclass

# End loader

add_library_search_dirs([])

# Begin libraries
_libs["libcblas.so"] = load_library("libcblas.so")

# 1 libraries
# End libraries

# No modules

enum_CBLAS_ORDER = c_int# /usr/include/cblas.h: 5

CblasRowMajor = 101# /usr/include/cblas.h: 5

CblasColMajor = 102# /usr/include/cblas.h: 5

enum_CBLAS_TRANSPOSE = c_int# /usr/include/cblas.h: 6

CblasNoTrans = 111# /usr/include/cblas.h: 6

CblasTrans = 112# /usr/include/cblas.h: 6

CblasConjTrans = 113# /usr/include/cblas.h: 6

AtlasConj = 114# /usr/include/cblas.h: 6

enum_CBLAS_UPLO = c_int# /usr/include/cblas.h: 8

CblasUpper = 121# /usr/include/cblas.h: 8

CblasLower = 122# /usr/include/cblas.h: 8

enum_CBLAS_DIAG = c_int# /usr/include/cblas.h: 9

CblasNonUnit = 131# /usr/include/cblas.h: 9

CblasUnit = 132# /usr/include/cblas.h: 9

enum_CBLAS_SIDE = c_int# /usr/include/cblas.h: 10

CblasLeft = 141# /usr/include/cblas.h: 10

CblasRight = 142# /usr/include/cblas.h: 10

# /usr/include/cblas.h: 17
for _lib in _libs.values():
    if _lib.has("cblas_errprn", "cdecl"):
        _func = _lib.get("cblas_errprn", "cdecl")
        _restype = c_int
        _errcheck = None
        _argtypes = [c_int, c_int, String]
        cblas_errprn = _variadic_function(_func,_restype,_argtypes,_errcheck)

# /usr/include/cblas.h: 18
if _libs["libcblas.so"].has("cblas_xerbla", "cdecl"):
    _func = _libs["libcblas.so"].get("cblas_xerbla", "cdecl")
    _restype = None
    _errcheck = None
    _argtypes = [c_int, String, String]
    cblas_xerbla = _variadic_function(_func,_restype,_argtypes,_errcheck)

# /usr/include/cblas.h: 25
if _libs["libcblas.so"].has("cblas_sdsdot", "cdecl"):
    cblas_sdsdot = _libs["libcblas.so"].get("cblas_sdsdot", "cdecl")
    cblas_sdsdot.argtypes = [c_int, c_float, POINTER(c_float), c_int, POINTER(c_float), c_int]
    cblas_sdsdot.restype = c_float

# /usr/include/cblas.h: 27
if _libs["libcblas.so"].has("cblas_dsdot", "cdecl"):
    cblas_dsdot = _libs["libcblas.so"].get("cblas_dsdot", "cdecl")
    cblas_dsdot.argtypes = [c_int, POINTER(c_float), c_int, POINTER(c_float), c_int]
    cblas_dsdot.restype = c_double

# /usr/include/cblas.h: 29
if _libs["libcblas.so"].has("cblas_sdot", "cdecl"):
    cblas_sdot = _libs["libcblas.so"].get("cblas_sdot", "cdecl")
    cblas_sdot.argtypes = [c_int, POINTER(c_float), c_int, POINTER(c_float), c_int]
    cblas_sdot.restype = c_float

# /usr/include/cblas.h: 31
if _libs["libcblas.so"].has("cblas_ddot", "cdecl"):
    cblas_ddot = _libs["libcblas.so"].get("cblas_ddot", "cdecl")
    cblas_ddot.argtypes = [c_int, POINTER(c_double), c_int, POINTER(c_double), c_int]
    cblas_ddot.restype = c_double

# /usr/include/cblas.h: 36
if _libs["libcblas.so"].has("cblas_cdotu_sub", "cdecl"):
    cblas_cdotu_sub = _libs["libcblas.so"].get("cblas_cdotu_sub", "cdecl")
    cblas_cdotu_sub.argtypes = [c_int, POINTER(None), c_int, POINTER(None), c_int, POINTER(None)]
    cblas_cdotu_sub.restype = None

# /usr/include/cblas.h: 38
if _libs["libcblas.so"].has("cblas_cdotc_sub", "cdecl"):
    cblas_cdotc_sub = _libs["libcblas.so"].get("cblas_cdotc_sub", "cdecl")
    cblas_cdotc_sub.argtypes = [c_int, POINTER(None), c_int, POINTER(None), c_int, POINTER(None)]
    cblas_cdotc_sub.restype = None

# /usr/include/cblas.h: 41
if _libs["libcblas.so"].has("cblas_zdotu_sub", "cdecl"):
    cblas_zdotu_sub = _libs["libcblas.so"].get("cblas_zdotu_sub", "cdecl")
    cblas_zdotu_sub.argtypes = [c_int, POINTER(None), c_int, POINTER(None), c_int, POINTER(None)]
    cblas_zdotu_sub.restype = None

# /usr/include/cblas.h: 43
if _libs["libcblas.so"].has("cblas_zdotc_sub", "cdecl"):
    cblas_zdotc_sub = _libs["libcblas.so"].get("cblas_zdotc_sub", "cdecl")
    cblas_zdotc_sub.argtypes = [c_int, POINTER(None), c_int, POINTER(None), c_int, POINTER(None)]
    cblas_zdotc_sub.restype = None

# /usr/include/cblas.h: 50
if _libs["libcblas.so"].has("cblas_snrm2", "cdecl"):
    cblas_snrm2 = _libs["libcblas.so"].get("cblas_snrm2", "cdecl")
    cblas_snrm2.argtypes = [c_int, POINTER(c_float), c_int]
    cblas_snrm2.restype = c_float

# /usr/include/cblas.h: 51
if _libs["libcblas.so"].has("cblas_sasum", "cdecl"):
    cblas_sasum = _libs["libcblas.so"].get("cblas_sasum", "cdecl")
    cblas_sasum.argtypes = [c_int, POINTER(c_float), c_int]
    cblas_sasum.restype = c_float

# /usr/include/cblas.h: 53
if _libs["libcblas.so"].has("cblas_dnrm2", "cdecl"):
    cblas_dnrm2 = _libs["libcblas.so"].get("cblas_dnrm2", "cdecl")
    cblas_dnrm2.argtypes = [c_int, POINTER(c_double), c_int]
    cblas_dnrm2.restype = c_double

# /usr/include/cblas.h: 54
if _libs["libcblas.so"].has("cblas_dasum", "cdecl"):
    cblas_dasum = _libs["libcblas.so"].get("cblas_dasum", "cdecl")
    cblas_dasum.argtypes = [c_int, POINTER(c_double), c_int]
    cblas_dasum.restype = c_double

# /usr/include/cblas.h: 56
if _libs["libcblas.so"].has("cblas_scnrm2", "cdecl"):
    cblas_scnrm2 = _libs["libcblas.so"].get("cblas_scnrm2", "cdecl")
    cblas_scnrm2.argtypes = [c_int, POINTER(None), c_int]
    cblas_scnrm2.restype = c_float

# /usr/include/cblas.h: 57
if _libs["libcblas.so"].has("cblas_scasum", "cdecl"):
    cblas_scasum = _libs["libcblas.so"].get("cblas_scasum", "cdecl")
    cblas_scasum.argtypes = [c_int, POINTER(None), c_int]
    cblas_scasum.restype = c_float

# /usr/include/cblas.h: 59
if _libs["libcblas.so"].has("cblas_dznrm2", "cdecl"):
    cblas_dznrm2 = _libs["libcblas.so"].get("cblas_dznrm2", "cdecl")
    cblas_dznrm2.argtypes = [c_int, POINTER(None), c_int]
    cblas_dznrm2.restype = c_double

# /usr/include/cblas.h: 60
if _libs["libcblas.so"].has("cblas_dzasum", "cdecl"):
    cblas_dzasum = _libs["libcblas.so"].get("cblas_dzasum", "cdecl")
    cblas_dzasum.argtypes = [c_int, POINTER(None), c_int]
    cblas_dzasum.restype = c_double

# /usr/include/cblas.h: 66
if _libs["libcblas.so"].has("cblas_isamax", "cdecl"):
    cblas_isamax = _libs["libcblas.so"].get("cblas_isamax", "cdecl")
    cblas_isamax.argtypes = [c_int, POINTER(c_float), c_int]
    cblas_isamax.restype = c_int

# /usr/include/cblas.h: 67
if _libs["libcblas.so"].has("cblas_idamax", "cdecl"):
    cblas_idamax = _libs["libcblas.so"].get("cblas_idamax", "cdecl")
    cblas_idamax.argtypes = [c_int, POINTER(c_double), c_int]
    cblas_idamax.restype = c_int

# /usr/include/cblas.h: 68
if _libs["libcblas.so"].has("cblas_icamax", "cdecl"):
    cblas_icamax = _libs["libcblas.so"].get("cblas_icamax", "cdecl")
    cblas_icamax.argtypes = [c_int, POINTER(None), c_int]
    cblas_icamax.restype = c_int

# /usr/include/cblas.h: 69
if _libs["libcblas.so"].has("cblas_izamax", "cdecl"):
    cblas_izamax = _libs["libcblas.so"].get("cblas_izamax", "cdecl")
    cblas_izamax.argtypes = [c_int, POINTER(None), c_int]
    cblas_izamax.restype = c_int

# /usr/include/cblas.h: 80
if _libs["libcblas.so"].has("cblas_sswap", "cdecl"):
    cblas_sswap = _libs["libcblas.so"].get("cblas_sswap", "cdecl")
    cblas_sswap.argtypes = [c_int, POINTER(c_float), c_int, POINTER(c_float), c_int]
    cblas_sswap.restype = None

# /usr/include/cblas.h: 82
if _libs["libcblas.so"].has("cblas_scopy", "cdecl"):
    cblas_scopy = _libs["libcblas.so"].get("cblas_scopy", "cdecl")
    cblas_scopy.argtypes = [c_int, POINTER(c_float), c_int, POINTER(c_float), c_int]
    cblas_scopy.restype = None

# /usr/include/cblas.h: 84
if _libs["libcblas.so"].has("cblas_saxpy", "cdecl"):
    cblas_saxpy = _libs["libcblas.so"].get("cblas_saxpy", "cdecl")
    cblas_saxpy.argtypes = [c_int, c_float, POINTER(c_float), c_int, POINTER(c_float), c_int]
    cblas_saxpy.restype = None

# /usr/include/cblas.h: 86
for _lib in _libs.values():
    if not _lib.has("catlas_saxpby", "cdecl"):
        continue
    catlas_saxpby = _lib.get("catlas_saxpby", "cdecl")
    catlas_saxpby.argtypes = [c_int, c_float, POINTER(c_float), c_int, c_float, POINTER(c_float), c_int]
    catlas_saxpby.restype = None
    break

# /usr/include/cblas.h: 88
for _lib in _libs.values():
    if not _lib.has("catlas_sset", "cdecl"):
        continue
    catlas_sset = _lib.get("catlas_sset", "cdecl")
    catlas_sset.argtypes = [c_int, c_float, POINTER(c_float), c_int]
    catlas_sset.restype = None
    break

# /usr/include/cblas.h: 91
if _libs["libcblas.so"].has("cblas_dswap", "cdecl"):
    cblas_dswap = _libs["libcblas.so"].get("cblas_dswap", "cdecl")
    cblas_dswap.argtypes = [c_int, POINTER(c_double), c_int, POINTER(c_double), c_int]
    cblas_dswap.restype = None

# /usr/include/cblas.h: 93
if _libs["libcblas.so"].has("cblas_dcopy", "cdecl"):
    cblas_dcopy = _libs["libcblas.so"].get("cblas_dcopy", "cdecl")
    cblas_dcopy.argtypes = [c_int, POINTER(c_double), c_int, POINTER(c_double), c_int]
    cblas_dcopy.restype = None

# /usr/include/cblas.h: 95
if _libs["libcblas.so"].has("cblas_daxpy", "cdecl"):
    cblas_daxpy = _libs["libcblas.so"].get("cblas_daxpy", "cdecl")
    cblas_daxpy.argtypes = [c_int, c_double, POINTER(c_double), c_int, POINTER(c_double), c_int]
    cblas_daxpy.restype = None

# /usr/include/cblas.h: 97
for _lib in _libs.values():
    if not _lib.has("catlas_daxpby", "cdecl"):
        continue
    catlas_daxpby = _lib.get("catlas_daxpby", "cdecl")
    catlas_daxpby.argtypes = [c_int, c_double, POINTER(c_double), c_int, c_double, POINTER(c_double), c_int]
    catlas_daxpby.restype = None
    break

# /usr/include/cblas.h: 99
for _lib in _libs.values():
    if not _lib.has("catlas_dset", "cdecl"):
        continue
    catlas_dset = _lib.get("catlas_dset", "cdecl")
    catlas_dset.argtypes = [c_int, c_double, POINTER(c_double), c_int]
    catlas_dset.restype = None
    break

# /usr/include/cblas.h: 102
if _libs["libcblas.so"].has("cblas_cswap", "cdecl"):
    cblas_cswap = _libs["libcblas.so"].get("cblas_cswap", "cdecl")
    cblas_cswap.argtypes = [c_int, POINTER(None), c_int, POINTER(None), c_int]
    cblas_cswap.restype = None

# /usr/include/cblas.h: 104
if _libs["libcblas.so"].has("cblas_ccopy", "cdecl"):
    cblas_ccopy = _libs["libcblas.so"].get("cblas_ccopy", "cdecl")
    cblas_ccopy.argtypes = [c_int, POINTER(None), c_int, POINTER(None), c_int]
    cblas_ccopy.restype = None

# /usr/include/cblas.h: 106
if _libs["libcblas.so"].has("cblas_caxpy", "cdecl"):
    cblas_caxpy = _libs["libcblas.so"].get("cblas_caxpy", "cdecl")
    cblas_caxpy.argtypes = [c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int]
    cblas_caxpy.restype = None

# /usr/include/cblas.h: 108
for _lib in _libs.values():
    if not _lib.has("catlas_caxpby", "cdecl"):
        continue
    catlas_caxpby = _lib.get("catlas_caxpby", "cdecl")
    catlas_caxpby.argtypes = [c_int, POINTER(None), POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    catlas_caxpby.restype = None
    break

# /usr/include/cblas.h: 110
for _lib in _libs.values():
    if not _lib.has("catlas_cset", "cdecl"):
        continue
    catlas_cset = _lib.get("catlas_cset", "cdecl")
    catlas_cset.argtypes = [c_int, POINTER(None), POINTER(None), c_int]
    catlas_cset.restype = None
    break

# /usr/include/cblas.h: 113
if _libs["libcblas.so"].has("cblas_zswap", "cdecl"):
    cblas_zswap = _libs["libcblas.so"].get("cblas_zswap", "cdecl")
    cblas_zswap.argtypes = [c_int, POINTER(None), c_int, POINTER(None), c_int]
    cblas_zswap.restype = None

# /usr/include/cblas.h: 115
if _libs["libcblas.so"].has("cblas_zcopy", "cdecl"):
    cblas_zcopy = _libs["libcblas.so"].get("cblas_zcopy", "cdecl")
    cblas_zcopy.argtypes = [c_int, POINTER(None), c_int, POINTER(None), c_int]
    cblas_zcopy.restype = None

# /usr/include/cblas.h: 117
if _libs["libcblas.so"].has("cblas_zaxpy", "cdecl"):
    cblas_zaxpy = _libs["libcblas.so"].get("cblas_zaxpy", "cdecl")
    cblas_zaxpy.argtypes = [c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int]
    cblas_zaxpy.restype = None

# /usr/include/cblas.h: 119
for _lib in _libs.values():
    if not _lib.has("catlas_zaxpby", "cdecl"):
        continue
    catlas_zaxpby = _lib.get("catlas_zaxpby", "cdecl")
    catlas_zaxpby.argtypes = [c_int, POINTER(None), POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    catlas_zaxpby.restype = None
    break

# /usr/include/cblas.h: 121
for _lib in _libs.values():
    if not _lib.has("catlas_zset", "cdecl"):
        continue
    catlas_zset = _lib.get("catlas_zset", "cdecl")
    catlas_zset.argtypes = [c_int, POINTER(None), POINTER(None), c_int]
    catlas_zset.restype = None
    break

# /usr/include/cblas.h: 128
if _libs["libcblas.so"].has("cblas_srotg", "cdecl"):
    cblas_srotg = _libs["libcblas.so"].get("cblas_srotg", "cdecl")
    cblas_srotg.argtypes = [POINTER(c_float), POINTER(c_float), POINTER(c_float), POINTER(c_float)]
    cblas_srotg.restype = None

# /usr/include/cblas.h: 129
if _libs["libcblas.so"].has("cblas_srotmg", "cdecl"):
    cblas_srotmg = _libs["libcblas.so"].get("cblas_srotmg", "cdecl")
    cblas_srotmg.argtypes = [POINTER(c_float), POINTER(c_float), POINTER(c_float), c_float, POINTER(c_float)]
    cblas_srotmg.restype = None

# /usr/include/cblas.h: 130
if _libs["libcblas.so"].has("cblas_srot", "cdecl"):
    cblas_srot = _libs["libcblas.so"].get("cblas_srot", "cdecl")
    cblas_srot.argtypes = [c_int, POINTER(c_float), c_int, POINTER(c_float), c_int, c_float, c_float]
    cblas_srot.restype = None

# /usr/include/cblas.h: 132
if _libs["libcblas.so"].has("cblas_srotm", "cdecl"):
    cblas_srotm = _libs["libcblas.so"].get("cblas_srotm", "cdecl")
    cblas_srotm.argtypes = [c_int, POINTER(c_float), c_int, POINTER(c_float), c_int, POINTER(c_float)]
    cblas_srotm.restype = None

# /usr/include/cblas.h: 135
if _libs["libcblas.so"].has("cblas_drotg", "cdecl"):
    cblas_drotg = _libs["libcblas.so"].get("cblas_drotg", "cdecl")
    cblas_drotg.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double)]
    cblas_drotg.restype = None

# /usr/include/cblas.h: 136
if _libs["libcblas.so"].has("cblas_drotmg", "cdecl"):
    cblas_drotmg = _libs["libcblas.so"].get("cblas_drotmg", "cdecl")
    cblas_drotmg.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), c_double, POINTER(c_double)]
    cblas_drotmg.restype = None

# /usr/include/cblas.h: 137
if _libs["libcblas.so"].has("cblas_drot", "cdecl"):
    cblas_drot = _libs["libcblas.so"].get("cblas_drot", "cdecl")
    cblas_drot.argtypes = [c_int, POINTER(c_double), c_int, POINTER(c_double), c_int, c_double, c_double]
    cblas_drot.restype = None

# /usr/include/cblas.h: 139
if _libs["libcblas.so"].has("cblas_drotm", "cdecl"):
    cblas_drotm = _libs["libcblas.so"].get("cblas_drotm", "cdecl")
    cblas_drotm.argtypes = [c_int, POINTER(c_double), c_int, POINTER(c_double), c_int, POINTER(c_double)]
    cblas_drotm.restype = None

# /usr/include/cblas.h: 146
if _libs["libcblas.so"].has("cblas_sscal", "cdecl"):
    cblas_sscal = _libs["libcblas.so"].get("cblas_sscal", "cdecl")
    cblas_sscal.argtypes = [c_int, c_float, POINTER(c_float), c_int]
    cblas_sscal.restype = None

# /usr/include/cblas.h: 147
if _libs["libcblas.so"].has("cblas_dscal", "cdecl"):
    cblas_dscal = _libs["libcblas.so"].get("cblas_dscal", "cdecl")
    cblas_dscal.argtypes = [c_int, c_double, POINTER(c_double), c_int]
    cblas_dscal.restype = None

# /usr/include/cblas.h: 148
if _libs["libcblas.so"].has("cblas_cscal", "cdecl"):
    cblas_cscal = _libs["libcblas.so"].get("cblas_cscal", "cdecl")
    cblas_cscal.argtypes = [c_int, POINTER(None), POINTER(None), c_int]
    cblas_cscal.restype = None

# /usr/include/cblas.h: 149
if _libs["libcblas.so"].has("cblas_zscal", "cdecl"):
    cblas_zscal = _libs["libcblas.so"].get("cblas_zscal", "cdecl")
    cblas_zscal.argtypes = [c_int, POINTER(None), POINTER(None), c_int]
    cblas_zscal.restype = None

# /usr/include/cblas.h: 150
if _libs["libcblas.so"].has("cblas_csscal", "cdecl"):
    cblas_csscal = _libs["libcblas.so"].get("cblas_csscal", "cdecl")
    cblas_csscal.argtypes = [c_int, c_float, POINTER(None), c_int]
    cblas_csscal.restype = None

# /usr/include/cblas.h: 151
if _libs["libcblas.so"].has("cblas_zdscal", "cdecl"):
    cblas_zdscal = _libs["libcblas.so"].get("cblas_zdscal", "cdecl")
    cblas_zdscal.argtypes = [c_int, c_double, POINTER(None), c_int]
    cblas_zdscal.restype = None

# /usr/include/cblas.h: 156
for _lib in _libs.values():
    if not _lib.has("cblas_crotg", "cdecl"):
        continue
    cblas_crotg = _lib.get("cblas_crotg", "cdecl")
    cblas_crotg.argtypes = [POINTER(None), POINTER(None), POINTER(None), POINTER(None)]
    cblas_crotg.restype = None
    break

# /usr/include/cblas.h: 157
for _lib in _libs.values():
    if not _lib.has("cblas_zrotg", "cdecl"):
        continue
    cblas_zrotg = _lib.get("cblas_zrotg", "cdecl")
    cblas_zrotg.argtypes = [POINTER(None), POINTER(None), POINTER(None), POINTER(None)]
    cblas_zrotg.restype = None
    break

# /usr/include/cblas.h: 158
for _lib in _libs.values():
    if not _lib.has("cblas_csrot", "cdecl"):
        continue
    cblas_csrot = _lib.get("cblas_csrot", "cdecl")
    cblas_csrot.argtypes = [c_int, POINTER(None), c_int, POINTER(None), c_int, c_float, c_float]
    cblas_csrot.restype = None
    break

# /usr/include/cblas.h: 160
for _lib in _libs.values():
    if not _lib.has("cblas_zdrot", "cdecl"):
        continue
    cblas_zdrot = _lib.get("cblas_zdrot", "cdecl")
    cblas_zdrot.argtypes = [c_int, POINTER(None), c_int, POINTER(None), c_int, c_double, c_double]
    cblas_zdrot.restype = None
    break

# /usr/include/cblas.h: 172
if _libs["libcblas.so"].has("cblas_sgemv", "cdecl"):
    cblas_sgemv = _libs["libcblas.so"].get("cblas_sgemv", "cdecl")
    cblas_sgemv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_TRANSPOSE, c_int, c_int, c_float, POINTER(c_float), c_int, POINTER(c_float), c_int, c_float, POINTER(c_float), c_int]
    cblas_sgemv.restype = None

# /usr/include/cblas.h: 177
if _libs["libcblas.so"].has("cblas_sgbmv", "cdecl"):
    cblas_sgbmv = _libs["libcblas.so"].get("cblas_sgbmv", "cdecl")
    cblas_sgbmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_TRANSPOSE, c_int, c_int, c_int, c_int, c_float, POINTER(c_float), c_int, POINTER(c_float), c_int, c_float, POINTER(c_float), c_int]
    cblas_sgbmv.restype = None

# /usr/include/cblas.h: 182
if _libs["libcblas.so"].has("cblas_strmv", "cdecl"):
    cblas_strmv = _libs["libcblas.so"].get("cblas_strmv", "cdecl")
    cblas_strmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, POINTER(c_float), c_int, POINTER(c_float), c_int]
    cblas_strmv.restype = None

# /usr/include/cblas.h: 186
if _libs["libcblas.so"].has("cblas_stbmv", "cdecl"):
    cblas_stbmv = _libs["libcblas.so"].get("cblas_stbmv", "cdecl")
    cblas_stbmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, c_int, POINTER(c_float), c_int, POINTER(c_float), c_int]
    cblas_stbmv.restype = None

# /usr/include/cblas.h: 190
if _libs["libcblas.so"].has("cblas_stpmv", "cdecl"):
    cblas_stpmv = _libs["libcblas.so"].get("cblas_stpmv", "cdecl")
    cblas_stpmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, POINTER(c_float), POINTER(c_float), c_int]
    cblas_stpmv.restype = None

# /usr/include/cblas.h: 193
if _libs["libcblas.so"].has("cblas_strsv", "cdecl"):
    cblas_strsv = _libs["libcblas.so"].get("cblas_strsv", "cdecl")
    cblas_strsv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, POINTER(c_float), c_int, POINTER(c_float), c_int]
    cblas_strsv.restype = None

# /usr/include/cblas.h: 197
if _libs["libcblas.so"].has("cblas_stbsv", "cdecl"):
    cblas_stbsv = _libs["libcblas.so"].get("cblas_stbsv", "cdecl")
    cblas_stbsv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, c_int, POINTER(c_float), c_int, POINTER(c_float), c_int]
    cblas_stbsv.restype = None

# /usr/include/cblas.h: 201
if _libs["libcblas.so"].has("cblas_stpsv", "cdecl"):
    cblas_stpsv = _libs["libcblas.so"].get("cblas_stpsv", "cdecl")
    cblas_stpsv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, POINTER(c_float), POINTER(c_float), c_int]
    cblas_stpsv.restype = None

# /usr/include/cblas.h: 205
if _libs["libcblas.so"].has("cblas_dgemv", "cdecl"):
    cblas_dgemv = _libs["libcblas.so"].get("cblas_dgemv", "cdecl")
    cblas_dgemv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_TRANSPOSE, c_int, c_int, c_double, POINTER(c_double), c_int, POINTER(c_double), c_int, c_double, POINTER(c_double), c_int]
    cblas_dgemv.restype = None

# /usr/include/cblas.h: 210
if _libs["libcblas.so"].has("cblas_dgbmv", "cdecl"):
    cblas_dgbmv = _libs["libcblas.so"].get("cblas_dgbmv", "cdecl")
    cblas_dgbmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_TRANSPOSE, c_int, c_int, c_int, c_int, c_double, POINTER(c_double), c_int, POINTER(c_double), c_int, c_double, POINTER(c_double), c_int]
    cblas_dgbmv.restype = None

# /usr/include/cblas.h: 215
if _libs["libcblas.so"].has("cblas_dtrmv", "cdecl"):
    cblas_dtrmv = _libs["libcblas.so"].get("cblas_dtrmv", "cdecl")
    cblas_dtrmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, POINTER(c_double), c_int, POINTER(c_double), c_int]
    cblas_dtrmv.restype = None

# /usr/include/cblas.h: 219
if _libs["libcblas.so"].has("cblas_dtbmv", "cdecl"):
    cblas_dtbmv = _libs["libcblas.so"].get("cblas_dtbmv", "cdecl")
    cblas_dtbmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, c_int, POINTER(c_double), c_int, POINTER(c_double), c_int]
    cblas_dtbmv.restype = None

# /usr/include/cblas.h: 223
if _libs["libcblas.so"].has("cblas_dtpmv", "cdecl"):
    cblas_dtpmv = _libs["libcblas.so"].get("cblas_dtpmv", "cdecl")
    cblas_dtpmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, POINTER(c_double), POINTER(c_double), c_int]
    cblas_dtpmv.restype = None

# /usr/include/cblas.h: 226
if _libs["libcblas.so"].has("cblas_dtrsv", "cdecl"):
    cblas_dtrsv = _libs["libcblas.so"].get("cblas_dtrsv", "cdecl")
    cblas_dtrsv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, POINTER(c_double), c_int, POINTER(c_double), c_int]
    cblas_dtrsv.restype = None

# /usr/include/cblas.h: 230
if _libs["libcblas.so"].has("cblas_dtbsv", "cdecl"):
    cblas_dtbsv = _libs["libcblas.so"].get("cblas_dtbsv", "cdecl")
    cblas_dtbsv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, c_int, POINTER(c_double), c_int, POINTER(c_double), c_int]
    cblas_dtbsv.restype = None

# /usr/include/cblas.h: 234
if _libs["libcblas.so"].has("cblas_dtpsv", "cdecl"):
    cblas_dtpsv = _libs["libcblas.so"].get("cblas_dtpsv", "cdecl")
    cblas_dtpsv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, POINTER(c_double), POINTER(c_double), c_int]
    cblas_dtpsv.restype = None

# /usr/include/cblas.h: 238
if _libs["libcblas.so"].has("cblas_cgemv", "cdecl"):
    cblas_cgemv = _libs["libcblas.so"].get("cblas_cgemv", "cdecl")
    cblas_cgemv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_TRANSPOSE, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_cgemv.restype = None

# /usr/include/cblas.h: 243
if _libs["libcblas.so"].has("cblas_cgbmv", "cdecl"):
    cblas_cgbmv = _libs["libcblas.so"].get("cblas_cgbmv", "cdecl")
    cblas_cgbmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_TRANSPOSE, c_int, c_int, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_cgbmv.restype = None

# /usr/include/cblas.h: 248
if _libs["libcblas.so"].has("cblas_ctrmv", "cdecl"):
    cblas_ctrmv = _libs["libcblas.so"].get("cblas_ctrmv", "cdecl")
    cblas_ctrmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, POINTER(None), c_int, POINTER(None), c_int]
    cblas_ctrmv.restype = None

# /usr/include/cblas.h: 252
if _libs["libcblas.so"].has("cblas_ctbmv", "cdecl"):
    cblas_ctbmv = _libs["libcblas.so"].get("cblas_ctbmv", "cdecl")
    cblas_ctbmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, c_int, POINTER(None), c_int, POINTER(None), c_int]
    cblas_ctbmv.restype = None

# /usr/include/cblas.h: 256
if _libs["libcblas.so"].has("cblas_ctpmv", "cdecl"):
    cblas_ctpmv = _libs["libcblas.so"].get("cblas_ctpmv", "cdecl")
    cblas_ctpmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, POINTER(None), POINTER(None), c_int]
    cblas_ctpmv.restype = None

# /usr/include/cblas.h: 259
if _libs["libcblas.so"].has("cblas_ctrsv", "cdecl"):
    cblas_ctrsv = _libs["libcblas.so"].get("cblas_ctrsv", "cdecl")
    cblas_ctrsv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, POINTER(None), c_int, POINTER(None), c_int]
    cblas_ctrsv.restype = None

# /usr/include/cblas.h: 263
if _libs["libcblas.so"].has("cblas_ctbsv", "cdecl"):
    cblas_ctbsv = _libs["libcblas.so"].get("cblas_ctbsv", "cdecl")
    cblas_ctbsv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, c_int, POINTER(None), c_int, POINTER(None), c_int]
    cblas_ctbsv.restype = None

# /usr/include/cblas.h: 267
if _libs["libcblas.so"].has("cblas_ctpsv", "cdecl"):
    cblas_ctpsv = _libs["libcblas.so"].get("cblas_ctpsv", "cdecl")
    cblas_ctpsv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, POINTER(None), POINTER(None), c_int]
    cblas_ctpsv.restype = None

# /usr/include/cblas.h: 271
if _libs["libcblas.so"].has("cblas_zgemv", "cdecl"):
    cblas_zgemv = _libs["libcblas.so"].get("cblas_zgemv", "cdecl")
    cblas_zgemv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_TRANSPOSE, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_zgemv.restype = None

# /usr/include/cblas.h: 276
if _libs["libcblas.so"].has("cblas_zgbmv", "cdecl"):
    cblas_zgbmv = _libs["libcblas.so"].get("cblas_zgbmv", "cdecl")
    cblas_zgbmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_TRANSPOSE, c_int, c_int, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_zgbmv.restype = None

# /usr/include/cblas.h: 281
if _libs["libcblas.so"].has("cblas_ztrmv", "cdecl"):
    cblas_ztrmv = _libs["libcblas.so"].get("cblas_ztrmv", "cdecl")
    cblas_ztrmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, POINTER(None), c_int, POINTER(None), c_int]
    cblas_ztrmv.restype = None

# /usr/include/cblas.h: 285
if _libs["libcblas.so"].has("cblas_ztbmv", "cdecl"):
    cblas_ztbmv = _libs["libcblas.so"].get("cblas_ztbmv", "cdecl")
    cblas_ztbmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, c_int, POINTER(None), c_int, POINTER(None), c_int]
    cblas_ztbmv.restype = None

# /usr/include/cblas.h: 289
if _libs["libcblas.so"].has("cblas_ztpmv", "cdecl"):
    cblas_ztpmv = _libs["libcblas.so"].get("cblas_ztpmv", "cdecl")
    cblas_ztpmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, POINTER(None), POINTER(None), c_int]
    cblas_ztpmv.restype = None

# /usr/include/cblas.h: 292
if _libs["libcblas.so"].has("cblas_ztrsv", "cdecl"):
    cblas_ztrsv = _libs["libcblas.so"].get("cblas_ztrsv", "cdecl")
    cblas_ztrsv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, POINTER(None), c_int, POINTER(None), c_int]
    cblas_ztrsv.restype = None

# /usr/include/cblas.h: 296
if _libs["libcblas.so"].has("cblas_ztbsv", "cdecl"):
    cblas_ztbsv = _libs["libcblas.so"].get("cblas_ztbsv", "cdecl")
    cblas_ztbsv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, c_int, POINTER(None), c_int, POINTER(None), c_int]
    cblas_ztbsv.restype = None

# /usr/include/cblas.h: 300
if _libs["libcblas.so"].has("cblas_ztpsv", "cdecl"):
    cblas_ztpsv = _libs["libcblas.so"].get("cblas_ztpsv", "cdecl")
    cblas_ztpsv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, POINTER(None), POINTER(None), c_int]
    cblas_ztpsv.restype = None

# /usr/include/cblas.h: 308
if _libs["libcblas.so"].has("cblas_ssymv", "cdecl"):
    cblas_ssymv = _libs["libcblas.so"].get("cblas_ssymv", "cdecl")
    cblas_ssymv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_float, POINTER(c_float), c_int, POINTER(c_float), c_int, c_float, POINTER(c_float), c_int]
    cblas_ssymv.restype = None

# /usr/include/cblas.h: 312
if _libs["libcblas.so"].has("cblas_ssbmv", "cdecl"):
    cblas_ssbmv = _libs["libcblas.so"].get("cblas_ssbmv", "cdecl")
    cblas_ssbmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_int, c_float, POINTER(c_float), c_int, POINTER(c_float), c_int, c_float, POINTER(c_float), c_int]
    cblas_ssbmv.restype = None

# /usr/include/cblas.h: 316
if _libs["libcblas.so"].has("cblas_sspmv", "cdecl"):
    cblas_sspmv = _libs["libcblas.so"].get("cblas_sspmv", "cdecl")
    cblas_sspmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_float, POINTER(c_float), POINTER(c_float), c_int, c_float, POINTER(c_float), c_int]
    cblas_sspmv.restype = None

# /usr/include/cblas.h: 320
if _libs["libcblas.so"].has("cblas_sger", "cdecl"):
    cblas_sger = _libs["libcblas.so"].get("cblas_sger", "cdecl")
    cblas_sger.argtypes = [enum_CBLAS_ORDER, c_int, c_int, c_float, POINTER(c_float), c_int, POINTER(c_float), c_int, POINTER(c_float), c_int]
    cblas_sger.restype = None

# /usr/include/cblas.h: 323
if _libs["libcblas.so"].has("cblas_ssyr", "cdecl"):
    cblas_ssyr = _libs["libcblas.so"].get("cblas_ssyr", "cdecl")
    cblas_ssyr.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_float, POINTER(c_float), c_int, POINTER(c_float), c_int]
    cblas_ssyr.restype = None

# /usr/include/cblas.h: 326
if _libs["libcblas.so"].has("cblas_sspr", "cdecl"):
    cblas_sspr = _libs["libcblas.so"].get("cblas_sspr", "cdecl")
    cblas_sspr.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_float, POINTER(c_float), c_int, POINTER(c_float)]
    cblas_sspr.restype = None

# /usr/include/cblas.h: 329
if _libs["libcblas.so"].has("cblas_ssyr2", "cdecl"):
    cblas_ssyr2 = _libs["libcblas.so"].get("cblas_ssyr2", "cdecl")
    cblas_ssyr2.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_float, POINTER(c_float), c_int, POINTER(c_float), c_int, POINTER(c_float), c_int]
    cblas_ssyr2.restype = None

# /usr/include/cblas.h: 333
if _libs["libcblas.so"].has("cblas_sspr2", "cdecl"):
    cblas_sspr2 = _libs["libcblas.so"].get("cblas_sspr2", "cdecl")
    cblas_sspr2.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_float, POINTER(c_float), c_int, POINTER(c_float), c_int, POINTER(c_float)]
    cblas_sspr2.restype = None

# /usr/include/cblas.h: 337
if _libs["libcblas.so"].has("cblas_dsymv", "cdecl"):
    cblas_dsymv = _libs["libcblas.so"].get("cblas_dsymv", "cdecl")
    cblas_dsymv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_double, POINTER(c_double), c_int, POINTER(c_double), c_int, c_double, POINTER(c_double), c_int]
    cblas_dsymv.restype = None

# /usr/include/cblas.h: 341
if _libs["libcblas.so"].has("cblas_dsbmv", "cdecl"):
    cblas_dsbmv = _libs["libcblas.so"].get("cblas_dsbmv", "cdecl")
    cblas_dsbmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_int, c_double, POINTER(c_double), c_int, POINTER(c_double), c_int, c_double, POINTER(c_double), c_int]
    cblas_dsbmv.restype = None

# /usr/include/cblas.h: 345
if _libs["libcblas.so"].has("cblas_dspmv", "cdecl"):
    cblas_dspmv = _libs["libcblas.so"].get("cblas_dspmv", "cdecl")
    cblas_dspmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_double, POINTER(c_double), POINTER(c_double), c_int, c_double, POINTER(c_double), c_int]
    cblas_dspmv.restype = None

# /usr/include/cblas.h: 349
if _libs["libcblas.so"].has("cblas_dger", "cdecl"):
    cblas_dger = _libs["libcblas.so"].get("cblas_dger", "cdecl")
    cblas_dger.argtypes = [enum_CBLAS_ORDER, c_int, c_int, c_double, POINTER(c_double), c_int, POINTER(c_double), c_int, POINTER(c_double), c_int]
    cblas_dger.restype = None

# /usr/include/cblas.h: 352
if _libs["libcblas.so"].has("cblas_dsyr", "cdecl"):
    cblas_dsyr = _libs["libcblas.so"].get("cblas_dsyr", "cdecl")
    cblas_dsyr.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_double, POINTER(c_double), c_int, POINTER(c_double), c_int]
    cblas_dsyr.restype = None

# /usr/include/cblas.h: 355
if _libs["libcblas.so"].has("cblas_dspr", "cdecl"):
    cblas_dspr = _libs["libcblas.so"].get("cblas_dspr", "cdecl")
    cblas_dspr.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_double, POINTER(c_double), c_int, POINTER(c_double)]
    cblas_dspr.restype = None

# /usr/include/cblas.h: 358
if _libs["libcblas.so"].has("cblas_dsyr2", "cdecl"):
    cblas_dsyr2 = _libs["libcblas.so"].get("cblas_dsyr2", "cdecl")
    cblas_dsyr2.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_double, POINTER(c_double), c_int, POINTER(c_double), c_int, POINTER(c_double), c_int]
    cblas_dsyr2.restype = None

# /usr/include/cblas.h: 362
if _libs["libcblas.so"].has("cblas_dspr2", "cdecl"):
    cblas_dspr2 = _libs["libcblas.so"].get("cblas_dspr2", "cdecl")
    cblas_dspr2.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_double, POINTER(c_double), c_int, POINTER(c_double), c_int, POINTER(c_double)]
    cblas_dspr2.restype = None

# /usr/include/cblas.h: 370
if _libs["libcblas.so"].has("cblas_chemv", "cdecl"):
    cblas_chemv = _libs["libcblas.so"].get("cblas_chemv", "cdecl")
    cblas_chemv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_chemv.restype = None

# /usr/include/cblas.h: 374
if _libs["libcblas.so"].has("cblas_chbmv", "cdecl"):
    cblas_chbmv = _libs["libcblas.so"].get("cblas_chbmv", "cdecl")
    cblas_chbmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_chbmv.restype = None

# /usr/include/cblas.h: 378
if _libs["libcblas.so"].has("cblas_chpmv", "cdecl"):
    cblas_chpmv = _libs["libcblas.so"].get("cblas_chpmv", "cdecl")
    cblas_chpmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, POINTER(None), POINTER(None), POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_chpmv.restype = None

# /usr/include/cblas.h: 382
if _libs["libcblas.so"].has("cblas_cgeru", "cdecl"):
    cblas_cgeru = _libs["libcblas.so"].get("cblas_cgeru", "cdecl")
    cblas_cgeru.argtypes = [enum_CBLAS_ORDER, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), c_int]
    cblas_cgeru.restype = None

# /usr/include/cblas.h: 385
if _libs["libcblas.so"].has("cblas_cgerc", "cdecl"):
    cblas_cgerc = _libs["libcblas.so"].get("cblas_cgerc", "cdecl")
    cblas_cgerc.argtypes = [enum_CBLAS_ORDER, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), c_int]
    cblas_cgerc.restype = None

# /usr/include/cblas.h: 388
if _libs["libcblas.so"].has("cblas_cher", "cdecl"):
    cblas_cher = _libs["libcblas.so"].get("cblas_cher", "cdecl")
    cblas_cher.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_float, POINTER(None), c_int, POINTER(None), c_int]
    cblas_cher.restype = None

# /usr/include/cblas.h: 391
if _libs["libcblas.so"].has("cblas_chpr", "cdecl"):
    cblas_chpr = _libs["libcblas.so"].get("cblas_chpr", "cdecl")
    cblas_chpr.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_float, POINTER(None), c_int, POINTER(None)]
    cblas_chpr.restype = None

# /usr/include/cblas.h: 394
if _libs["libcblas.so"].has("cblas_cher2", "cdecl"):
    cblas_cher2 = _libs["libcblas.so"].get("cblas_cher2", "cdecl")
    cblas_cher2.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), c_int]
    cblas_cher2.restype = None

# /usr/include/cblas.h: 397
if _libs["libcblas.so"].has("cblas_chpr2", "cdecl"):
    cblas_chpr2 = _libs["libcblas.so"].get("cblas_chpr2", "cdecl")
    cblas_chpr2.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None)]
    cblas_chpr2.restype = None

# /usr/include/cblas.h: 401
if _libs["libcblas.so"].has("cblas_zhemv", "cdecl"):
    cblas_zhemv = _libs["libcblas.so"].get("cblas_zhemv", "cdecl")
    cblas_zhemv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_zhemv.restype = None

# /usr/include/cblas.h: 405
if _libs["libcblas.so"].has("cblas_zhbmv", "cdecl"):
    cblas_zhbmv = _libs["libcblas.so"].get("cblas_zhbmv", "cdecl")
    cblas_zhbmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_zhbmv.restype = None

# /usr/include/cblas.h: 409
if _libs["libcblas.so"].has("cblas_zhpmv", "cdecl"):
    cblas_zhpmv = _libs["libcblas.so"].get("cblas_zhpmv", "cdecl")
    cblas_zhpmv.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, POINTER(None), POINTER(None), POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_zhpmv.restype = None

# /usr/include/cblas.h: 413
if _libs["libcblas.so"].has("cblas_zgeru", "cdecl"):
    cblas_zgeru = _libs["libcblas.so"].get("cblas_zgeru", "cdecl")
    cblas_zgeru.argtypes = [enum_CBLAS_ORDER, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), c_int]
    cblas_zgeru.restype = None

# /usr/include/cblas.h: 416
if _libs["libcblas.so"].has("cblas_zgerc", "cdecl"):
    cblas_zgerc = _libs["libcblas.so"].get("cblas_zgerc", "cdecl")
    cblas_zgerc.argtypes = [enum_CBLAS_ORDER, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), c_int]
    cblas_zgerc.restype = None

# /usr/include/cblas.h: 419
if _libs["libcblas.so"].has("cblas_zher", "cdecl"):
    cblas_zher = _libs["libcblas.so"].get("cblas_zher", "cdecl")
    cblas_zher.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_double, POINTER(None), c_int, POINTER(None), c_int]
    cblas_zher.restype = None

# /usr/include/cblas.h: 422
if _libs["libcblas.so"].has("cblas_zhpr", "cdecl"):
    cblas_zhpr = _libs["libcblas.so"].get("cblas_zhpr", "cdecl")
    cblas_zhpr.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, c_double, POINTER(None), c_int, POINTER(None)]
    cblas_zhpr.restype = None

# /usr/include/cblas.h: 425
if _libs["libcblas.so"].has("cblas_zher2", "cdecl"):
    cblas_zher2 = _libs["libcblas.so"].get("cblas_zher2", "cdecl")
    cblas_zher2.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), c_int]
    cblas_zher2.restype = None

# /usr/include/cblas.h: 428
if _libs["libcblas.so"].has("cblas_zhpr2", "cdecl"):
    cblas_zhpr2 = _libs["libcblas.so"].get("cblas_zhpr2", "cdecl")
    cblas_zhpr2.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None)]
    cblas_zhpr2.restype = None

# /usr/include/cblas.h: 441
if _libs["libcblas.so"].has("cblas_sgemm", "cdecl"):
    cblas_sgemm = _libs["libcblas.so"].get("cblas_sgemm", "cdecl")
    cblas_sgemm.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_TRANSPOSE, enum_CBLAS_TRANSPOSE, c_int, c_int, c_int, c_float, POINTER(c_float), c_int, POINTER(c_float), c_int, c_float, POINTER(c_float), c_int]
    cblas_sgemm.restype = None

# /usr/include/cblas.h: 446
if _libs["libcblas.so"].has("cblas_ssymm", "cdecl"):
    cblas_ssymm = _libs["libcblas.so"].get("cblas_ssymm", "cdecl")
    cblas_ssymm.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_SIDE, enum_CBLAS_UPLO, c_int, c_int, c_float, POINTER(c_float), c_int, POINTER(c_float), c_int, c_float, POINTER(c_float), c_int]
    cblas_ssymm.restype = None

# /usr/include/cblas.h: 451
if _libs["libcblas.so"].has("cblas_ssyrk", "cdecl"):
    cblas_ssyrk = _libs["libcblas.so"].get("cblas_ssyrk", "cdecl")
    cblas_ssyrk.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, c_int, c_int, c_float, POINTER(c_float), c_int, c_float, POINTER(c_float), c_int]
    cblas_ssyrk.restype = None

# /usr/include/cblas.h: 455
if _libs["libcblas.so"].has("cblas_ssyr2k", "cdecl"):
    cblas_ssyr2k = _libs["libcblas.so"].get("cblas_ssyr2k", "cdecl")
    cblas_ssyr2k.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, c_int, c_int, c_float, POINTER(c_float), c_int, POINTER(c_float), c_int, c_float, POINTER(c_float), c_int]
    cblas_ssyr2k.restype = None

# /usr/include/cblas.h: 460
if _libs["libcblas.so"].has("cblas_strmm", "cdecl"):
    cblas_strmm = _libs["libcblas.so"].get("cblas_strmm", "cdecl")
    cblas_strmm.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_SIDE, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, c_int, c_float, POINTER(c_float), c_int, POINTER(c_float), c_int]
    cblas_strmm.restype = None

# /usr/include/cblas.h: 465
if _libs["libcblas.so"].has("cblas_strsm", "cdecl"):
    cblas_strsm = _libs["libcblas.so"].get("cblas_strsm", "cdecl")
    cblas_strsm.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_SIDE, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, c_int, c_float, POINTER(c_float), c_int, POINTER(c_float), c_int]
    cblas_strsm.restype = None

# /usr/include/cblas.h: 471
if _libs["libcblas.so"].has("cblas_dgemm", "cdecl"):
    cblas_dgemm = _libs["libcblas.so"].get("cblas_dgemm", "cdecl")
    cblas_dgemm.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_TRANSPOSE, enum_CBLAS_TRANSPOSE, c_int, c_int, c_int, c_double, POINTER(c_double), c_int, POINTER(c_double), c_int, c_double, POINTER(c_double), c_int]
    cblas_dgemm.restype = None

# /usr/include/cblas.h: 476
if _libs["libcblas.so"].has("cblas_dsymm", "cdecl"):
    cblas_dsymm = _libs["libcblas.so"].get("cblas_dsymm", "cdecl")
    cblas_dsymm.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_SIDE, enum_CBLAS_UPLO, c_int, c_int, c_double, POINTER(c_double), c_int, POINTER(c_double), c_int, c_double, POINTER(c_double), c_int]
    cblas_dsymm.restype = None

# /usr/include/cblas.h: 481
if _libs["libcblas.so"].has("cblas_dsyrk", "cdecl"):
    cblas_dsyrk = _libs["libcblas.so"].get("cblas_dsyrk", "cdecl")
    cblas_dsyrk.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, c_int, c_int, c_double, POINTER(c_double), c_int, c_double, POINTER(c_double), c_int]
    cblas_dsyrk.restype = None

# /usr/include/cblas.h: 485
if _libs["libcblas.so"].has("cblas_dsyr2k", "cdecl"):
    cblas_dsyr2k = _libs["libcblas.so"].get("cblas_dsyr2k", "cdecl")
    cblas_dsyr2k.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, c_int, c_int, c_double, POINTER(c_double), c_int, POINTER(c_double), c_int, c_double, POINTER(c_double), c_int]
    cblas_dsyr2k.restype = None

# /usr/include/cblas.h: 490
if _libs["libcblas.so"].has("cblas_dtrmm", "cdecl"):
    cblas_dtrmm = _libs["libcblas.so"].get("cblas_dtrmm", "cdecl")
    cblas_dtrmm.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_SIDE, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, c_int, c_double, POINTER(c_double), c_int, POINTER(c_double), c_int]
    cblas_dtrmm.restype = None

# /usr/include/cblas.h: 495
if _libs["libcblas.so"].has("cblas_dtrsm", "cdecl"):
    cblas_dtrsm = _libs["libcblas.so"].get("cblas_dtrsm", "cdecl")
    cblas_dtrsm.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_SIDE, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, c_int, c_double, POINTER(c_double), c_int, POINTER(c_double), c_int]
    cblas_dtrsm.restype = None

# /usr/include/cblas.h: 501
if _libs["libcblas.so"].has("cblas_cgemm", "cdecl"):
    cblas_cgemm = _libs["libcblas.so"].get("cblas_cgemm", "cdecl")
    cblas_cgemm.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_TRANSPOSE, enum_CBLAS_TRANSPOSE, c_int, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_cgemm.restype = None

# /usr/include/cblas.h: 506
if _libs["libcblas.so"].has("cblas_csymm", "cdecl"):
    cblas_csymm = _libs["libcblas.so"].get("cblas_csymm", "cdecl")
    cblas_csymm.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_SIDE, enum_CBLAS_UPLO, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_csymm.restype = None

# /usr/include/cblas.h: 511
if _libs["libcblas.so"].has("cblas_csyrk", "cdecl"):
    cblas_csyrk = _libs["libcblas.so"].get("cblas_csyrk", "cdecl")
    cblas_csyrk.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_csyrk.restype = None

# /usr/include/cblas.h: 515
if _libs["libcblas.so"].has("cblas_csyr2k", "cdecl"):
    cblas_csyr2k = _libs["libcblas.so"].get("cblas_csyr2k", "cdecl")
    cblas_csyr2k.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_csyr2k.restype = None

# /usr/include/cblas.h: 520
if _libs["libcblas.so"].has("cblas_ctrmm", "cdecl"):
    cblas_ctrmm = _libs["libcblas.so"].get("cblas_ctrmm", "cdecl")
    cblas_ctrmm.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_SIDE, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int]
    cblas_ctrmm.restype = None

# /usr/include/cblas.h: 525
if _libs["libcblas.so"].has("cblas_ctrsm", "cdecl"):
    cblas_ctrsm = _libs["libcblas.so"].get("cblas_ctrsm", "cdecl")
    cblas_ctrsm.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_SIDE, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int]
    cblas_ctrsm.restype = None

# /usr/include/cblas.h: 531
if _libs["libcblas.so"].has("cblas_zgemm", "cdecl"):
    cblas_zgemm = _libs["libcblas.so"].get("cblas_zgemm", "cdecl")
    cblas_zgemm.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_TRANSPOSE, enum_CBLAS_TRANSPOSE, c_int, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_zgemm.restype = None

# /usr/include/cblas.h: 536
if _libs["libcblas.so"].has("cblas_zsymm", "cdecl"):
    cblas_zsymm = _libs["libcblas.so"].get("cblas_zsymm", "cdecl")
    cblas_zsymm.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_SIDE, enum_CBLAS_UPLO, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_zsymm.restype = None

# /usr/include/cblas.h: 541
if _libs["libcblas.so"].has("cblas_zsyrk", "cdecl"):
    cblas_zsyrk = _libs["libcblas.so"].get("cblas_zsyrk", "cdecl")
    cblas_zsyrk.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_zsyrk.restype = None

# /usr/include/cblas.h: 545
if _libs["libcblas.so"].has("cblas_zsyr2k", "cdecl"):
    cblas_zsyr2k = _libs["libcblas.so"].get("cblas_zsyr2k", "cdecl")
    cblas_zsyr2k.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_zsyr2k.restype = None

# /usr/include/cblas.h: 550
if _libs["libcblas.so"].has("cblas_ztrmm", "cdecl"):
    cblas_ztrmm = _libs["libcblas.so"].get("cblas_ztrmm", "cdecl")
    cblas_ztrmm.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_SIDE, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int]
    cblas_ztrmm.restype = None

# /usr/include/cblas.h: 555
if _libs["libcblas.so"].has("cblas_ztrsm", "cdecl"):
    cblas_ztrsm = _libs["libcblas.so"].get("cblas_ztrsm", "cdecl")
    cblas_ztrsm.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_SIDE, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, enum_CBLAS_DIAG, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int]
    cblas_ztrsm.restype = None

# /usr/include/cblas.h: 565
if _libs["libcblas.so"].has("cblas_chemm", "cdecl"):
    cblas_chemm = _libs["libcblas.so"].get("cblas_chemm", "cdecl")
    cblas_chemm.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_SIDE, enum_CBLAS_UPLO, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_chemm.restype = None

# /usr/include/cblas.h: 570
if _libs["libcblas.so"].has("cblas_cherk", "cdecl"):
    cblas_cherk = _libs["libcblas.so"].get("cblas_cherk", "cdecl")
    cblas_cherk.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, c_int, c_int, c_float, POINTER(None), c_int, c_float, POINTER(None), c_int]
    cblas_cherk.restype = None

# /usr/include/cblas.h: 574
if _libs["libcblas.so"].has("cblas_cher2k", "cdecl"):
    cblas_cher2k = _libs["libcblas.so"].get("cblas_cher2k", "cdecl")
    cblas_cher2k.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, c_float, POINTER(None), c_int]
    cblas_cher2k.restype = None

# /usr/include/cblas.h: 579
if _libs["libcblas.so"].has("cblas_zhemm", "cdecl"):
    cblas_zhemm = _libs["libcblas.so"].get("cblas_zhemm", "cdecl")
    cblas_zhemm.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_SIDE, enum_CBLAS_UPLO, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, POINTER(None), POINTER(None), c_int]
    cblas_zhemm.restype = None

# /usr/include/cblas.h: 584
if _libs["libcblas.so"].has("cblas_zherk", "cdecl"):
    cblas_zherk = _libs["libcblas.so"].get("cblas_zherk", "cdecl")
    cblas_zherk.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, c_int, c_int, c_double, POINTER(None), c_int, c_double, POINTER(None), c_int]
    cblas_zherk.restype = None

# /usr/include/cblas.h: 588
if _libs["libcblas.so"].has("cblas_zher2k", "cdecl"):
    cblas_zher2k = _libs["libcblas.so"].get("cblas_zher2k", "cdecl")
    cblas_zher2k.argtypes = [enum_CBLAS_ORDER, enum_CBLAS_UPLO, enum_CBLAS_TRANSPOSE, c_int, c_int, POINTER(None), POINTER(None), c_int, POINTER(None), c_int, c_double, POINTER(None), c_int]
    cblas_zher2k.restype = None

# /usr/include/cblas.h: 594
for _lib in _libs.values():
    if _lib.has("cblas_errprn", "cdecl"):
        _func = _lib.get("cblas_errprn", "cdecl")
        _restype = c_int
        _errcheck = None
        _argtypes = [c_int, c_int, String]
        cblas_errprn = _variadic_function(_func,_restype,_argtypes,_errcheck)

CBLAS_INDEX = c_int# /usr/include/cblas.h: 15

# No inserted files

# No prefix-stripping

