import os

platform = ARGUMENTS.get('OS', Platform())

# Expose TPL dependencies to all aspects of the build

tpl_inc = ['#packages/fftw/include', '#packages/teuchos/src']
tpl_lib = ['#packages/fftw/lib', '#packages/teuchos/src']

bas_env = Environment( PLATFORM = platform,
                       CXXFLAGS = ' -DWJB_SCONS',
                       INCDIR  = tpl_inc,
                       LIBDIR  = tpl_lib,
                       CPPPATH = [tpl_inc],
                       LIBPATH = [tpl_lib],
                       LIBS = ['teuchos', 'fftw3'] )

debug = ARGUMENTS.get('debug', 0)
if int(debug):
  bas_env.Prepend(CCFLAGS = '-g')
else:
  bas_env.Prepend(CCFLAGS = '-O2')


bas_env.SConscript(['src/SConscript',
                    'test/SConscript',
                    'packages/teuchos/SConscript',
                    'packages/fftw/SConscript'], exports='bas_env')


# WJB: consider moving inflater into packages/SConscript

def buildTeuchosInflater(target, source, env):
  # Whatever it takes to inflate teuchos (cvs checkout vs. extract tarball)

  source_name = str(source[0])
  print "source = ", source_name
  target_name = str(target[0])
  print "target = ", target_name, "does not exist;  INFLATING... "

  #os.system('tar tvzf ' + source_name + ' | tee ' + target_name)
  os.system('tar xvzf ' + source_name + ' > ' + target_name)
  return

#bas_env.Command('#packages/teuchos/fileList.log', '#packages/teuchos.tgz',
                #buildTeuchosInflater)


dict = bas_env.Dictionary()
####for key in ['OBJSUFFIX', 'LIBSUFFIX', 'PLATFORM']:
for key in ['CC', 'CPPPATH', 'CCFLAGS']:
  print "keyTOP = %s, value = %s" % (key, dict[key])

