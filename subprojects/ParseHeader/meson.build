project(
    'ParseHeader',
    'cpp',
    default_options: [
        'cpp_std=c++17',
        'debug=true',
        'optimization=3',
        'b_ndebug=false',
    ],
)

# for std::filesystem on older compilers
add_project_link_arguments(['-lstdc++fs'], language: 'cpp')

fmt_dep = dependency('fmt', required: true)

flex = find_program('flex')
bison = find_program('bison')

phParser = custom_target(
    input: ['src/phParser.yy'],
    output: [
        'phParser.tab.hh',
        'phParser.tab.cc',
        'stack.hh',
        'location.hh',
        'position.hh',
    ],
    command: [bison, '-b@OUTDIR@/@BASENAME@', '@INPUT@'],
)

phParser_header = phParser[0]

phScanner = custom_target(
    input: ['src/phScanner.ll'],
    output: ['phScanner.cc', 'phScanner.hh'],
    command: [flex, '-o@OUTPUT0@', '--header-file=@OUTPUT1@', '@INPUT@'],
)

src = [
    'src/HeaderStream.cc',
    'src/ParseHeader.cc',
    'src/phDriver.cc',
    phParser,
    phScanner,
    'src/stringutil.cc',
]
public_inc = include_directories('include')
# private_inc = include_directories('src')
libparseheader = library(
    'parseheader',
    src,
    include_directories: [public_inc],
    install: true,
    dependencies: [fmt_dep],
    # TODO: ParseHeader currently unsafe under sanitizers
    override_options: ['b_sanitize=none'],
)

parseheader_dep = declare_dependency(
    sources: [phParser_header],
    include_directories: public_inc,
    link_with: libparseheader,
    dependencies: [fmt_dep],
)
