from django.contrib import admin
from django.urls import path, include
from django.conf import settings
from processMotif import views
from django.conf.urls.static import static

urlpatterns = [
    path('upload/', views.upload, name="upload"),
    path('find/', views.find, name="find"),
]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)