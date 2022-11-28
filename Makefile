CXXFLAGS = -Wall -O3
OBJS = Textures.o Image.o lib/lodepng/lodepng.o
GROUP = d
ASS = 10
.PHONY: clean zip

asdf: main
	-rm -f result.ppm 
	./main

main: $(OBJS)

clean:
	-rm -f *.o
	-rm main

zip: 
	zip -r group_$(GROUP)_assignment_$(ASS).zip main.cpp Textures.h Textures.cpp Image.h Image.cpp
