CXXFLAGS = -Wall -O3
# CXXFLAGS = -Wall -O3 -DANTI_ALIASING
SOURCE_FILES = Textures.cpp Image.cpp lib/lodepng/lodepng.cpp
OBJS = $(SOURCE_FILES:%.cpp=objects/%.o)

GROUP = d
ASS = 10
.PHONY: clean zip

asdf: main
	-rm -f result.ppm 
	./main

main: $(OBJS)

$(OBJS): objects/%.o: %.h

$(OBJS): objects/%.o: %.cpp
	@echo "Compiling" $<
	@mkdir -p $(@D)
	@gcc $(CFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	-rm -f $(OBJS)
	-rm main

zip: 
	zip -r group_$(GROUP)_assignment_$(ASS).zip main.cpp Textures.h Textures.cpp Image.h Image.cpp lib objects result.ppm textures Makefile Material.h glm
