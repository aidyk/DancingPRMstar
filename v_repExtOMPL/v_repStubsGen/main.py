import argparse
import tempita
import os
import re
import sys

from parse import parse
import model

parser = argparse.ArgumentParser(description='Generate stubs for V-REP plugin.')
parser.add_argument('xml_file', type=str, default='',
                   help='the XML file with the callback definitions')
parser.add_argument('--cpp', '-C', type=str, default='',
                   help='the C++ source file to generate')
parser.add_argument('--hpp', '-H', type=str, default='',
                   help='the C++ header file to generate')
parser.add_argument('--include', '-I', type=str, default='',
                   help='extra header file to include')
args = parser.parse_args()

plugin = parse(args.xml_file)

template_data = dict()
template_data['args'] = args
template_data['plugin'] = plugin
template_data['model'] = model
template_data['header_file'] = args.hpp if args.hpp else re.sub(r'\.c(|xx|pp|c)$','.h',args.cpp)

mkpath = lambda f: os.path.join(os.path.dirname(os.path.realpath(__file__)), f)

template = dict()
template[mkpath('stubs.h.in')] = args.hpp
template[mkpath('stubs.cpp.in')] = args.cpp

for template_file, out_file in template.items():
    if out_file:
        template_data['out_file'] = out_file
        with open(template_file, 'r') as f:
            template_text = f.read()
        tmpl = tempita.Template(template_text)
        f = sys.stdout if out_file == '-' else open(out_file, 'w')
        f.write(tmpl.substitute(**template_data))
        f.close()

