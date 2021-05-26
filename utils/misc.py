import sys
import os
import functools
import time
import datetime

from numpy.random import choice

def get_indices(l):
	"""
	[Takes a list and returns dict giving indices matching each possible 
	list member]
	Example:
		Input [0, 1, 1, 0, 0]
		Output {0:[0,3,4], 1:[1,2]}
	"""
	ret=dict()
	for member in set(l):
		ret[member] = list()
	i=0
	for el in l:
		ret[el].append(i)
		i+=1
	return(ret)

def all_zero(l):
	"""
	[Returns TRUE if supplied list contains all zeros
	Returns FALSE if list contains ANY non-zero values
	Returns FALSE if list is empty]
	"""
	values=set(l)
	if len(values) > 1:
		return(False)
	elif len(values)==1 and l[0] in [0, 0.0, "0", "0.0"]:
		return(True)
	else:
		return(False)

def weighted_draw(d, num_samples=1):
	choices = list(d.keys())
	weights = list(d.values())
	return(choice(choices, num_samples, p=weights))

def timer(func):
	"""[print the runtime of the decorated function in the format HH:MM:SS]"""
	@functools.wraps(func)
	def wrapper_timer(*args, **kwargs):
		start_time = time.perf_counter()
		value = func(*args, **kwargs)
		end_time = time.perf_counter()
		run_time = end_time - start_time
		final_runtime = str(datetime.timedelta(seconds=run_time))
		print(f"Finshed {func.__name__!r} in {final_runtime}\n")
		return value
	return wrapper_timer

def progressbar(it, prefix="", size=60, file=sys.stdout):
	count = len(it)
	def show(j):
		x = int(size*j/count)
		file.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), j, count))
		file.flush()
	show(0)
	for i, item in enumerate(it):
		yield item
		show(i+1)
	file.write("\n")
	file.flush()

def isnotebook():
	try:
		shell = get_ipython().__class__.__name__
		if shell == 'ZMQInteractiveShell':
			# Jupyter notebook or qtconsole
			return True
		elif shell == 'TerminalInteractiveShell':
			# Terminal running IPython
			return False
		else:
			# Other type (?)
			return False  
	except NameError:
		# Probably standard Python interpreter
		return False      
