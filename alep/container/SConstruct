# -*-python-*-

import os
from openalea.sconsx import config, environ

ALEASolution = config.ALEASolution
pj = os.path.join

name = 'container'

options = Variables('options.py', ARGUMENTS)
tools = ['boost_python']

env = ALEASolution(options, tools)
prefix = env['build_prefix']

# Build stage
SConscript(pj(prefix,"src/cpp/SConscript"), exports={"env": env})
SConscript(pj(prefix,"src/wrapper/SConscript"), exports={"env":env})
Default("build")
