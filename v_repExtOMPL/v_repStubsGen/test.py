from function import *
from struct import *
from variable import *
from enum import *

a1 = Variable('s','const char*','NULL')
f = Function('foo', args=[a1], body=['if(s)',['print(s);'],'else',['print("Hello world!");']])
s = Struct('Bar', fields=[a1])
e = Enum('Baz', ['a', 'b', 'c'], 23)

objs = [f, s, e]

for obj in objs:
    print(obj.declaration())

for obj in objs:
    print(obj.definition())
