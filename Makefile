CXXFLAGS = -Wall -O3
OBJS = Textures.o Image.o

main: $(OBJS)

.PHONY: clean
clean:
	rm -f *.o
	rm main
