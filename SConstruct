#
# Copyright (c) 2025      High Performance Computing Center Stuttgart,
#                         University of Stuttgart.  All rights reserved.
#
# Authors: Kingshuk Haldar <kingshuk.haldar@hlrs.de>
#

import os
import subprocess
import re

Decider('MD5-timestamp')

AddOption('--compdb', action='store_true',
          help='produce compilation database')

AddOption('--build-type', '--bt', dest='buildtype', nargs=1,
          type='choice', choices=['release', 'debug'],
          action='store', default='release',
          help='specify build type [release/debug]')

AddOption('--asan', action='store_true',
          help='employ address sanitizer')

AddOption('--ubsan', action='store_true',
          help='employ undefined behaviour sanitizer')

AddOption('--cflags', action='store', type= 'string', default= '',
          help= 'provide comma-separated compile flags')

env= Environment(LIBS= 'm', ENV= os.environ)
if 'CC' in env['ENV']:
    env['CC']= env['ENV']['CC']

env.Append(CPPDEFINES= ['CLOCKTALK_NORELEASE'])  # between-release development
env.Append(CPPDEFINES= [('_POSIX_C_SOURCE', 200809),#strdup, clock_gettime
                        ('CLOCKTALK_VERSION_MAJOR', 0),
                        ('CLOCKTALK_VERSION_MINOR', 0),
                        ('CLOCKTALK_VERSION_PATCH', 1)])

env.Append(CFLAGS= ['-std=c17', '-Wall', '-Werror'])
env.Append(CFLAGS= ['-g3', '-O0'] if 'debug'== GetOption('buildtype') else '-O3')

if GetOption('asan'):
  env.Append(CFLAGS= ['-fsanitize=address', '-fno-omit-frame-pointer'])
  env.Append(LINKFLAGS= '-fsanitize=address')

if GetOption('ubsan'):
  env.Append(CFLAGS= '-fsanitize=undefined')
  env.Append(LINKFLAGS= '-fsanitize=undefined')

if GetOption('cflags'):
    env.Append(CFAGS= re.split(r'[,]', GetOption('cflags')))

env.Append(LINKFLAGS= '-Wl,--no-undefined')

env.Append(CPATH= [os.getcwd()+ '/src'])

if GetOption('compdb'):
    env.Tool('compilation_db')
    env.CompilationDatabase('compile_commands.json')

clocktalk= env.SConscript('src/SConscript',
                          exports= 'env',
                          variant_dir= 'build',
                          duplicate= 0
)
env.Install(target='install/bin', source=clocktalk, duplicate=0)
