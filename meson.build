project('zeldovich-PLT',
        'cpp',
        default_options : ['cpp_std=c++11', 'debug=true', 'optimization=3'],
        # default_options : ['cpp_std=c++11', 'debug=true', 'b_sanitize=address'],
        )
add_global_arguments('-march=native', language : 'cpp')

disk = get_option('DISK')
if disk
  add_project_arguments(['-DDISK'], language: 'cpp')
endif

ph_proj = subproject('ParseHeader')
libparseheader_dep = ph_proj.get_variable('libparseheader_dep')

src = [
    'src/block_array.cpp',
    'src/output.cpp',
    'src/parameters.cpp',
    'src/power_spectrum.cpp',
    'src/STimer.cc',
    'src/zeldovich.cpp',
]

inc = include_directories('include')

dep = [dependency('gsl'),
       dependency('fftw3'),
       dependency('openmp'),
       libparseheader_dep,
       ]

exe = executable('zeldovich',
                 src,
                 include_directories : inc,
                 dependencies : dep,
                 install : true,
                 install_rpath : '$ORIGIN/../lib64',
                 )
