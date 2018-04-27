modules = ['tools','utils']

for module in modules:
	exec('from '+module+' import *')

