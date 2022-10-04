CXXFLAGS = -O3 -Wall

# main: main.cpp

video: main
	./main
	ffmpeg -y -framerate 30 -i ./frames/result_%d.ppm -c:v libx264 -crf 25 -vf "scale=1024:768,format=yuv420p" -movflags +faststart output.mp4                                                                                   
	mpv output.mp4 --loop
