from time import strftime
import os

def printTime(x, fileName = None):
	if fileName:
		fileCon = open(fileName, "a")
		print(strftime("%c"), x, file = fileCon)
		fileCon.close()
	else:
		print(strftime("%c"), x)
	
def dictSubset(x, keys):
	
	assert isinstance(x, dict)
	assert isinstance(keys, (tuple, list))
	xSubset = filter(lambda i:i[0] in keys, x.items())
	return xSubset

def sortDictByValue(x, reverse = False):
	
	assert isinstance(x, dict)
	xSort = sorted(x.items(), key = lambda i: i[1], reverse = reverse)
	return xSort

def sortDictByKey(x, reverse = False):
	
	assert isinstance(x, dict)
	xSort = sorted(x.items(), key = lambda i: i[0], reverse = reverse)
	return  xSort

def basename2(x):
	
	return os.path.splitext(os.path.basename(x))[0]
