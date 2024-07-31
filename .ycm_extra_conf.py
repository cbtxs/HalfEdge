
from distutils.sysconfig import get_python_inc
import sys
import platform
import os
import os.path as p
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

DIR_OF_THIS_SCRIPT = p.abspath(p.dirname(__file__))
DIR_OF_THIRD_PARTY = p.join(DIR_OF_THIS_SCRIPT, 'third_party')
SOURCE_EXTENSIONS = ['.cpp', '.cxx', '.cc', '.c', '.m', '.mm', '.inl']

library_include_paths = {
    'vtk': os.path.join(local, 'vtk', 'include', 'vtk-9.2'),
    'gmsh': os.path.join(local, 'include'),
    'mfem': os.path.join(local, 'mfem', 'include'),
    'cgal': os.path.join(local, 'cgal', 'include'),
    'googletest': os.path.join(local, 'googletest', 'include')
}

include_flags = [item for path in library_include_paths.values() for item in ('-I', path)]
flags = [
    '-Wall',
    '-Wextra',
    '-Werror',
    '-std=c++17',
    '-isystem',
    '-x', 'c++',
    '-isystem', get_python_inc(),
    '-isystem', '/usr/include/x86_64-linux-gnu/c++/' + major_version,
    '-isystem', '/usr/include/c++/' + major_version,
    '-isystem', ' /usr/include/c++/' + major_version + '/backward',
    '-isystem', '/usr/lib/gcc/x86_64-linux-gnu/' + major_version + '/include',
    '-isystem', '/usr/include',
    '-isystem', '/usr/include/cairomm-1.0',
    '-isystem', '/usr/include/cairo',
    '-isystem', '/usr/include/glib-2.0',
    '-isystem', '/usr/include/harfbuzz',
    '-isystem', '/usr/include/gtk-3.0',
    '-isystem', '/usr/include/atk-1.0',
    '-isystem', '/usr/include/freetype2',
    '-isystem', '/usr/include/gdk-pixbuf-2.0',
    '-isystem', '/usr/include/pango-1.0',
    '-isystem', '/usr/include/sigc++-2.0',
    '-isystem', '/usr/lib/x86_64-linux-gnu/glib-2.0/include',
    '-isystem', '/usr/lib/x86_64-linux-gnu/cairomm-1.0/include',
    '-isystem', '/usr/lib/x86_64-linux-gnu/sigc++-2.0/include',
    '-isystem', '/usr/lib/gcc/x86_64-linux-gnu/11/include',
    '-isystem', '/usr/include/eigen3',
    '-isystem', '/usr/local/include',
    '-isystem', mpi_include,
    '-I', '/usr/lib/x86_64-linux-gnu/openmpi/include/',
    '-I', './include',
    '-I', './thirdparty/include',
] + include_flags

def FlagsForFile(filename, **kwargs):
    extension = os.path.splitext(filename)[1]
    if extension in SOURCE_EXTENSIONS:
        final_flags = flags
    else:
        final_flags = []

    return {
        'flags': final_flags,
        'do_cache': True
    }


compilation_database_folder = ''

if p.exists( compilation_database_folder ):
  database = ycm_core.CompilationDatabase( compilation_database_folder )
else:
  database = None


def IsHeaderFile( filename ):
  extension = p.splitext( filename )[ 1 ]
  return extension in [ '.h', '.hxx', '.hpp', '.hh' ]


def FindCorrespondingSourceFile( filename ):
  if IsHeaderFile( filename ):
    basename = p.splitext( filename )[ 0 ]
    for extension in SOURCE_EXTENSIONS:
      replacement_file = basename + extension
      if p.exists( replacement_file ):
        return replacement_file
  return filename


def PathToPythonUsedDuringBuild():
  try:
    filepath = p.join( DIR_OF_THIS_SCRIPT, 'PYTHON_USED_DURING_BUILDING' )
    with open( filepath ) as f:
      return f.read().strip()
  # We need to check for IOError for Python 2 and OSError for Python 3.
  except ( IOError, OSError ):
    return None


def Settings( **kwargs ):
  language = kwargs[ 'language' ]

  if language == 'cfamily':
    # If the file is a header, try to find the corresponding source file and
    # retrieve its flags from the compilation database if using one. This is
    # necessary since compilation databases don't have entries for header files.
    # In addition, use this source file as the translation unit. This makes it
    # possible to jump from a declaration in the header file to its definition
    # in the corresponding source file.
    filename = FindCorrespondingSourceFile( kwargs[ 'filename' ] )

    if not database:
      return {
        'flags': flags,
        'include_paths_relative_to_dir': DIR_OF_THIS_SCRIPT,
        'override_filename': filename
      }

    compilation_info = database.GetCompilationInfoForFile( filename )
    if not compilation_info.compiler_flags_:
      return {}

    # Bear in mind that compilation_info.compiler_flags_ does NOT return a
    # python list, but a "list-like" StringVec object.
    final_flags = list( compilation_info.compiler_flags_ )

    # NOTE: This is just for YouCompleteMe; it's highly likely that your project
    # does NOT need to remove the stdlib flag. DO NOT USE THIS IN YOUR
    # ycm_extra_conf IF YOU'RE NOT 100% SURE YOU NEED IT.
    try:
      final_flags.remove( '-stdlib=libc++' )
    except ValueError:
      pass

    return {
      'flags': final_flags,
      'include_paths_relative_to_dir': compilation_info.compiler_working_dir_,
      'override_filename': filename
    }

  if language == 'python':
    return {
      'interpreter_path': PathToPythonUsedDuringBuild()
    }

  return {}


def GetStandardLibraryIndexInSysPath( sys_path ):
  for index, path in enumerate( sys_path ):
    if p.isfile( p.join( path, 'os.py' ) ):
      return index
  raise RuntimeError( 'Could not find standard library path in Python path.' )


def PythonSysPath( **kwargs ):
  sys_path = kwargs[ 'sys_path' ]

  interpreter_path = kwargs[ 'interpreter_path' ]
  major_version = subprocess.check_output( [
    interpreter_path, '-c', 'import sys; print( sys.version_info[ 0 ] )' ]
  ).rstrip().decode( 'utf8' )

  sys_path.insert( GetStandardLibraryIndexInSysPath( sys_path ) + 1,
                   p.join( DIR_OF_THIRD_PARTY, 'python-future', 'src' ) )
  sys_path[ 0:0 ] = [ p.join( DIR_OF_THIS_SCRIPT ),
                      p.join( DIR_OF_THIRD_PARTY, 'bottle' ),
                      p.join( DIR_OF_THIRD_PARTY, 'cregex',
                              'regex_{}'.format( major_version ) ),
                      p.join( DIR_OF_THIRD_PARTY, 'frozendict' ),
                      p.join( DIR_OF_THIRD_PARTY, 'jedi_deps', 'jedi' ),
                      p.join( DIR_OF_THIRD_PARTY, 'jedi_deps', 'numpydoc' ),
                      p.join( DIR_OF_THIRD_PARTY, 'jedi_deps', 'parso' ),
                      p.join( DIR_OF_THIRD_PARTY, 'requests_deps', 'requests' ),
                      p.join( DIR_OF_THIRD_PARTY, 'requests_deps',
                                                  'urllib3',
                                                  'src' ),
                      p.join( DIR_OF_THIRD_PARTY, 'requests_deps',
                                                  'chardet' ),
                      p.join( DIR_OF_THIRD_PARTY, 'requests_deps',
                                                  'certifi' ),
                      p.join( DIR_OF_THIRD_PARTY, 'requests_deps',
                                                  'idna' ),
                      p.join( DIR_OF_THIRD_PARTY, 'waitress' ) ]

  return sys_path
