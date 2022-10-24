CXXFLAGS = -Wall -O3
OBJS = Textures.o Image.o
GROUP = f
ASS = 5

main: $(OBJS)

.PHONY: clean zip
clean:
	rm -f *.o
	rm main

zip: 
	zip -r group_$(GROUP)_assignment_$(ASS).zip main.cpp Textures.h Textures.cpp Image.h Image.cpp
