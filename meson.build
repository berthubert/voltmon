project('voltmon', 'cpp', default_options : ['cpp_std=c++17'])

sqlitedep = dependency('sqlite3', version : '>3')
thread_dep = dependency('threads')
json_dep = dependency('nlohmann_json')
fmt_dep = dependency('fmt', static: true)

#cpphttplib = dependency('cpp-httplib')
sqlitewriter_dep = dependency('sqlitewriter', static: true)
doctest_dep=dependency('doctest')
#simplesockets_dep = dependency('simplesockets', static: true)
argparse_dep = dependency('argparse', version: '>=3')

executable('streamint', 'streamint.cc', 
	dependencies: [sqlitedep, json_dep, fmt_dep, sqlitewriter_dep,argparse_dep])

executable('metachan', 'metachan.cc',
	dependencies: [argparse_dep, json_dep])
       
