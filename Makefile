BD = ./build
SD = ./source
ID = ./include
OBJS = $(BD)/main.o $(BD)/image.o $(BD)/capsule.o $(BD)/interface.o $(BD)/objects.o $(BD)/config.o $(BD)/script.o
CC = g++ -O3
LINKING_GSL = -lgsl -lgslcblas -lm
LINKING_OPENCV = -lopencv_core -lopencv_highgui
LINKING_BOOST = -lboost_system -lboost_filesystem -lboost_program_options
FLAG_OPENMP = -fopenmp

OpenCapsule:			$(OBJS)
				$(CC) -Wall -o OpenCapsule $(OBJS) $(FLAG_OPENMP) $(LINKING_OPENCV) $(LINKING_GSL) $(LINKING_BOOST) 

$(BD)/main.o:			$(SD)/main.cpp	
				$(CC) -c -I$(ID) $(SD)/main.cpp -o $(BD)/main.o

$(BD)/image.o:			$(SD)/image.cpp $(ID)/image.h 
				$(CC) -c -I$(ID) $(SD)/image.cpp -o $(BD)/image.o

$(BD)/capsule.o:	    	$(SD)/capsule.cpp $(ID)/capsule.h 
				$(CC) -c -I$(ID) $(FLAG_OPENMP) $(SD)/capsule.cpp -o $(BD)/capsule.o

$(BD)/interface.o: 		$(SD)/interface.cpp $(ID)/interface.h $
				$(CC) -c -I$(ID) $(FLAG_OPENMP) $(SD)/interface.cpp -o $(BD)/interface.o

$(BD)/objects.o:    		$(SD)/objects.cpp $(ID)/objects.h
				$(CC) -c -I$(ID) $(SD)/objects.cpp -o $(BD)/objects.o

$(BD)/config.o:			$(SD)/config.cpp $(ID)/config.h
				$(CC) -c -I$(ID) $(SD)/config.cpp -o $(BD)/config.o

$(BD)/script.o:			$(SD)/script.cpp $(ID)/script.h
				$(CC) -c -I$(ID) $(SD)/script.cpp -o $(BD)/script.o

clean:
				rm -f *.o *~ ./out/* ./tmp/* ./global_out/* $(BD)/* $(SD)/*~ $(ID)/*~		

install:
				cp OpenCapsule /usr/local/bin/

