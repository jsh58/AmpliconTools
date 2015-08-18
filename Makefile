all: removePrimer qualTrim stitch

removePrimer: removePrimer.c removePrimer.h
	gcc -g -Wall -O3 -std=c99 -o removePrimer removePrimer.c

qualTrim: qualTrim.c qualTrim.h
	gcc -g -Wall -O3 -std=c99 -o qualTrim qualTrim.c

stitch: stitch.c stitch.h
	gcc -g -Wall -O3 -std=c99 -o stitch stitch.c
