project('ParseHeader',
        'cpp',
        default_options : ['cpp_std=c++11', 'debug=true', 'optimization=3'],
        )

flex = find_program('flex')
bison = find_program('bison')

phParser = custom_target(
    input : ['src/phParser.yy'],
    output : ['phParser.tab.cc', 'phParser.tab.hh', 'stack.hh', 'location.hh', 'position.hh'],
    command : [bison, '-b@OUTDIR@/@BASENAME@', '@INPUT@'],
)

phScanner = custom_target(
    input : ['src/phScanner.ll'],
    output : ['phScanner.cc', 'phScanner.hh'],
    command : [flex, '-o@OUTPUT0@', '--header-file=@OUTPUT1@', '@INPUT@'],
)

src = ['src/HeaderStream.cc',
       'src/ParseHeader.cc',
       'src/phDriver.cc',
       phParser,
       phScanner,
       'src/stringutil.cc',
       ]
inc = include_directories('include')
libparseheader = library('parseheader',
                         src,
                         include_directories : inc,
                         install : true,
                         )

libparseheader_dep = declare_dependency(
    include_directories : inc,
    link_with : libparseheader,
    )
