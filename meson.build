project(
    'zeldovich-PLT',
    'cpp',
    default_options: ['cpp_std=gnu++17', 'debug=true', 'optimization=3', 'warning_level=3'],
    # default_options : ['cpp_std=gnu++17', 'debug=true', 'b_sanitize=address,undefined'],
    meson_version: '>=1.1.0',
)
add_project_arguments('-march=native', language: 'cpp')

if get_option('optimization') == '3'
    add_project_arguments(['-Ofast'], language: 'cpp')
endif

disk = get_option('DISK')
if disk
    add_project_arguments(['-DDISK'], language: 'cpp')
endif

ph_proj = subproject('ParseHeader', default_options: ['b_sanitize=none'])
parseheader_dep = ph_proj.get_variable('parseheader_dep')

src = [
    'src/block_array.cpp',
    'src/output.cpp',
    'src/parameters.cpp',
    'src/power_spectrum.cpp',
    'src/STimer.cc',
    'src/zeldovich.cpp',
]

inc = include_directories('include')

dep = [
    dependency('gsl'),
    dependency('fftw3'),
    dependency('fmt'),
    dependency('openmp'),
    parseheader_dep,
]

exe = executable(
    'zeldovich',
    src,
    include_directories: inc,
    dependencies: dep,
    install: true,
    install_rpath: '$ORIGIN/../lib64',
)
