# This script generates polyamri.h at the top level of the build directory.
import sys, os, os.path
source_dir = sys.argv[1]
target_dir = sys.argv[2]

def visit_c_files(headers, dirname, files):
    if 'tests' in files:
        files.remove('tests')
    h_files = [f for f in files if f[-2:] == '.h' and f[0] != '.']
    headers.extend(h_files)

def find_c_headers(dirname):
    headers = []
    os.path.walk(dirname, visit_c_files, headers)
    return headers

# Generate the top-level header.
c_headers = find_c_headers('%s/polyamri'%source_dir)
header_name = '%s/polyamri.h'%target_dir
if not os.path.exists(target_dir):
    os.makedirs(target_dir)
header = open(header_name, 'w')
header.write('// polyamri.h -- automatically generated.\n')
header.write('// This file is part of the polyamri HPC library. See the license\n')
header.write('// in the actual source files for details of distribution.\n\n')
header.write('#ifndef POLYAMRI_LIBRARY_H\n')
header.write('#define POLYAMRI_LIBRARY_H\n\n')
header.write('#ifdef __cplusplus\n')
header.write('extern "C" {\n')
header.write('#endif\n\n')
header.write('#include "polyamri/polyamri.h"\n')
for h in c_headers:
    if h != 'polyamri.h':
        header.write('#include "polyamri/%s"\n'%h)
header.write('\n')
header.write('#ifdef __cplusplus\n')
header.write('}\n')
header.write('#endif\n\n')
header.write('#endif\n\n')
header.close()

