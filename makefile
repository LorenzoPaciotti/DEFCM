all: defc_v2
	gcc -ggdb -g -Wall -lm -o defc_v2 defc_v2.c
	gcc -ggdb -g -Wall -lm -o fcm fcm_standard_paciotti_test.c
	gcc -ggdb -g -Wall -lm -o data_generator data_generator.c
clean:
	rm defc_v2 fcm
