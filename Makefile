all: removePrimer qualTrim stitch

removePrimer: removePrimer.c removePrimer.h
	gcc -g -Wall -O3 -std=c99 -o removePrimer removePrimer.c -lz

qualTrim: qualTrim.c qualTrim.h
	gcc -g -Wall -O3 -std=c99 -o qualTrim qualTrim.c -lz

stitch: stitch.c stitch.h
	gcc -g -Wall -O3 -std=c99 -o stitch stitch.c -lz
