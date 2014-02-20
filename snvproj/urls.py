#NEW APP
from django.conf.urls import patterns, include, url
import snv.views

# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
		url(r'^$', snv.views.home, name='home',),
        url(r'^info/$', snv.views.info, name='info',),
        url(r'^about/$', snv.views.about, name='about',),
        url(r'^getting-started/$', snv.views.gettingstarted, name='getting-started',),
        #url(r'^$', snv.views.home.as_view(), name='home',),
        url(r'^listuniprot$', snv.views.UniprotList.as_view(), name='uniprot-list',),
        url(r'^uniprot/(?P<pk>(\D+|\d+)+)/$', snv.views.UniprotView.as_view(), name='uniprot-view',),
        url(r'^disease/(?P<pk>\d+)/$', snv.views.DiseaseView.as_view(), name='disease-view',),
		url(r'^snv/(?P<pk>(\D+|\d+)+)/$', snv.views.SnvView.as_view(), name='snv-view',),
        url(r'^interaction/(?P<pk>(\D+|\d+)+)/$',snv.views.InteractionView.as_view(), name='interaction-view'),
        url(r'^login/$', 'django.contrib.auth.views.login'),
        url(r'^logout/$', 'django.contrib.auth.views.logout'),
        url(r'^search-form/$', snv.views.search_form),
        url(r'^search/$', snv.views.search, name='search-results',),
        url(r'^viewer/(?P<uniprot_acc>(\D+|\d+)+)/$',snv.views.SuperpositionView.as_view(),name='superposition-view'),
        #url(r'^media/(?P<path>.*)$', 'django.views.static.serve',{'document_root': settings.MEDIA_ROOT}),

    # url(r'^snvproj/', include('snvproj.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    # url(r'^admin/', include(admin.site.urls)),

)
