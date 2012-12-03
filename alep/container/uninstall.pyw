import os
import shutil

def d (name) :
	try :
		os.remove(name)
	except Exception :
		print "file not found %s" % name

def dd (name) :
	try :
		shutil.rmtree(name)
	except Exception :
		print "dir not found %s" % name

d("C:/Python25/Lib/site-packages/shared_libs/container.dll")
d("C:/Python25/Lib/site-packages/shared_libs/container.dll.egm")
d("C:/Python25/Lib/site-packages/shared_libs/libcontainer.a.egm")
d("C:/Python25/Lib/site-packages/shared_libs/libcontainer.a")

dd("C:/Python25/Lib/site-packages/openalea.container-2.0.0-py2.5-win32.egg")

