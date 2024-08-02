
from distutils.sysconfig import get_python_inc
import sys
import os
import subprocess

# 将 YCM 安装路径替换为实际路径
ycm_path = os.path.expanduser("~/.vim/bundle/YouCompleteMe/third_party/ycmd")

# 将路径添加到 sys.path
if ycm_path not in sys.path:
    sys.path.append(ycm_path)

import ycm_core

version = subprocess.getoutput("gcc --version | head -n1 | awk '{print $NF}'")
major_version = version.split('.')[0]

home_dir = os.path.expanduser("~")
local = os.path.join(home_dir, '.local')

result = subprocess.run(["mpicc", "--showme:incdirs"], stdout=subprocess.PIPE)
mpi_include = result.stdout.decode().strip().split(' ')[-1]

library_include_paths = {
    'vtk': os.path.join(local, 'vtk', 'include', 'vtk-9.2'),
    'gmsh': os.path.join(local, 'include'),
    'mfem': os.path.join(local, 'mfem', 'include'),
    'cgal': os.path.join(local, 'cgal', 'include'),
    'googletest': os.path.join(local, 'googletest', 'include')
}

include_flags = [item for path in library_include_paths.values() for item in ('-isystem', path)]
flags = [
    '-Wall',
    '-Wextra',
    '-Werror',
    '-std=c++20',
    '-isystem',
    '-x', 'c++',
    '-isystem', get_python_inc(),
    '-isystem', '/usr/include/x86_64-linux-gnu/c++/' + major_version,
    '-isystem', '/usr/include/c++/' + major_version,
    '-isystem', ' /usr/include/c++/' + major_version + '/backward',
    '-isystem', '/usr/lib/gcc/x86_64-linux-gnu/' + major_version + '/include',
    '-isystem', '/usr/include',
    '-isystem', '/usr/include/eigen3',
    '-isystem', '/usr/local/include',
    '-isystem', mpi_include,
    '-isystem', '/usr/lib/x86_64-linux-gnu/openmpi/include/',
    '-isystem', './include',
] + include_flags

def Settings( **kwargs ):
  language = kwargs[ 'language' ]

  if language == 'cfamily':
    # If the file is a header, try to find the corresponding source file and
    # retrieve its flags from the compilation database if using one. This is
    # necessary since compilation databases don't have entries for header files.
    # In addition, use this source file as the translation unit. This makes it
    # possible to jump from a declaration in the header file to its definition
    # in the corresponding source file.
    return { 'flags': flags}

  return {}



