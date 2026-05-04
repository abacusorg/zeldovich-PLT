project(
    'zeldovich-PLT',
    'cpp',
    default_options: ['cpp_std=gnu++17', 'debug=true', 'optimization=3', 'warning_level=3'],
    # default_options : ['cpp_std=gnu++17', 'debug=true', 'b_sanitize=address,undefined'],
    meson_version: '>=1.1.0',
)

add_project_arguments('-march=' + get_option('march'), language: 'cpp')

cxx = meson.get_compiler('cpp')

if get_option('optimization') == '3'
    # Promote optimization=3 to -Ofast by default
    if cxx.get_id() == 'intel-llvm'
        opt_flag_ = '-fp-model=fast=1'
    elif cxx.get_id() == 'clang'
        opt_flag_ = '-ffast-math'
    else
        opt_flag_ = '-Ofast'
    endif
    add_project_arguments(opt_flag_, language: 'cpp')
endif

disk = get_option('DISK')
if disk
    add_project_arguments(['-DDISK'], language: 'cpp')
endif

fs = import('fs')

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
cpp_args = []

gsl_dep = dependency('gsl', required: get_option('allow_v1'))

if gsl_dep.found()
    cpp_args += ['-DHAVE_GSL']
endif

dep = [
    dependency('fftw3'),
    dependency('fmt'),
    dependency('openmp'),
    gsl_dep,
    parseheader_dep,
]

exe = executable(
    'zeldovich',
    src,
    include_directories: inc,
    dependencies: dep,
    install: true,
    install_rpath: '$ORIGIN/../lib64:$ORIGIN/../lib',
    cpp_args: cpp_args,
)

fs.copyfile(
  'eigmodes128',
)
