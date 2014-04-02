#!/project/soft/linux64/epd_free-7.3-2-rh5-x86_64/bin/python
import os, sys
sys.path.append('/project/home/wo213/dev/snvdb')
os.environ['DJANGO_SETTINGS_MODULE'] = 'snvproj.settings'
import django.core.handlers.wsgi
application = django.core.handlers.wsgi.WSGIHandler()

