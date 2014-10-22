all:
	gcc -ggdb -g -Wall -o defc_v2 defc_v2.c -lm
	gcc -ggdb -g -Wall -o fcm fcm_standard_paciotti_test.c -lm
	gcc -ggdb -g -Wall -o data_generator data_generator.c -lm
clean:
	rm defc_v2 fcm
