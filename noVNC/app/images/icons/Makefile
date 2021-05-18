ICONS := \
	bwb-16x16.png \
	bwb-24x24.png \
	bwb-32x32.png \
	bwb-48x48.png \
	bwb-64x64.png

ANDROID_LAUNCHER := \
	bwb-48x48.png \
	bwb-72x72.png \
	bwb-96x96.png \
	bwb-144x144.png \
	bwb-192x192.png

IPHONE_LAUNCHER := \
	bwb-60x60.png \
	bwb-120x120.png

IPAD_LAUNCHER := \
	bwb-76x76.png \
	bwb-152x152.png

ALL_ICONS := $(ICONS) $(ANDROID_LAUNCHER) $(IPHONE_LAUNCHER) $(IPAD_LAUNCHER)

all: $(ALL_ICONS)

bwb-%.png: bwb-icon.svg
	convert -density $(shell expr `echo $* | cut -d x -f 1` \* 90 / 48) -resize $* \
		-background transparent "$<" "$@"

clean:
	rm -f *.png
