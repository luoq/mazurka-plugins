PLUGIN_LIBRARY_NAME := mazurka-plugins
VAMP_SDK_INCLUDE := /usr/include/vamp-sdk
VAMP_SDK_LIB := /usr/lib
CXXFLAGS := -I $(VAMP_SDK_INCLUDE) -I../include/ -Wall -fPIC -O2
LDFLAGS := -shared -Wl,-soname=$(PLUGIN) $(VAMP_SDK_LIB)/libvamp-sdk.a -lfftw3
PLUGIN_EXT := .so

PLUGIN := $(PLUGIN_LIBRARY_NAME)$(PLUGIN_EXT)
PLUGIN_SOURCES := MazurkaPlugin.cpp MzChronogram.cpp MzPitchPower.cpp MzSilenceFinder.cpp MzSpectrogramClient.cpp MzSummation.cpp \
	MazurkaTransformer.cpp  MzHarmonicSpectrum.cpp  MzPowerCurve.cpp  MzSpectralFlatness.cpp  MzSpectrogramFFTW.cpp \
	MazurkaWindower.cpp     MzNevermore.cpp         MzPowerscape.cpp  MzSpectralFlux.cpp      MzSpectrogramHost.cpp \
	plugins.cpp

PLUGIN_OBJECTS := $(PLUGIN_SOURCES:.cpp=.o)
PLUGIN_OBJECTS := $(PLUGIN_OBJECTS:.c=.o)

$(PLUGIN): $(PLUGIN_OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS)
plugins.cpp:
	./makepluginscpp > plugins.cpp
clean:
	rm plugins.cpp
	rm -f *.o *.so
