#
# Copyright (c) 2025      High Performance Computing Center Stuttgart,
#                         University of Stuttgart.  All rights reserved.
#
# Authors: Kingshuk Haldar <kingshuk.haldar@hlrs.de>
#

import os
import subprocess
import re

AddOption('--compdb', dest='compdb', nargs=1, type='int', action='store', default=0,
          help='produce compilation database [0/1]')
AddOption('--build-type', '--bt', dest='buildtype', nargs=1, type='choice', choices=['release', 'debug'], action='store', default='release',
          help='specify build type [release/debug]')
AddOption('--compile-link-opts', '--clo', dest='clopts', nargs=1, type='string', action='store',
          help='provide other compile-cum-link flags')

env= Environment(LIBS= 'm', ENV= os.environ)
if 'CC' in env['ENV']:
    env['CC']= env['ENV']['CC']

env.Append(CPPDEFINES= [('_POSIX_C_SOURCE', 200809),#strdup, clock_gettime
                        ('CLOCKTALK_VERSION_MAJOR', 0),
                        ('CLOCKTALK_VERSION_MINOR', 0),
                        ('CLOCKTALK_VERSION_PATCH', 1)])

env.Append(CFLAGS= ['-std=c17', '-Wall', '-Werror'])
env.Append(CFLAGS= ['-g3', '-O0'] if GetOption('buildtype')== 'debug' else '-O3')

env.Append(LINKFLAGS= '-Wl,--no-undefined')

env.Append(CPATH= [os.getcwd()+ '/src'])

if None!= GetOption('clopts'):
    opts= re.split(r'[,]', GetOption('clopts'))
    env.Append(CFLAGS= opts)
    env.Append(LINKFLAGS= opts)

if(GetOption('compdb')!= 0):
    env.Tool('compilation_db')
    env.CompilationDatabase('compile_commands.json')

clocktalk= env.SConscript('src/SConscript',
                          exports= 'env',
                          variant_dir= 'build',
                          duplicate= 0
)
env.Install(target='install/bin', source=clocktalk, duplicate=0)
