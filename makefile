all: defc_v2
	gcc -ggdb -g -Wall -o defc_v2 defc_v2.c -lm
	gcc -ggdb -g -Wall -o fcm fcm_standard_paciotti_test.c -lm
clean:
	rm ./*
