from distutils.core import setup, Extension

# define the extension module
find_seq_module = Extension('find_seq_module', sources=['find_seq_module.c'])

# run the setup
setup(ext_modules=[find_seq_module])
