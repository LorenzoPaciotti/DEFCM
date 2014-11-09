all:
	gcc -ggdb -g -Wall -o defc2 defc_v2.c -lm
	gcc -ggdb -g -Wall -o fcm fcm_standard_paciotti_test.c -lm
	gcc -ggdb -g -Wall -o data_gen1 data_generator.c -lm
	gcc -ggdb -g -Wall -o data_gen2 data_generator_2.c -lm
	gcc -ggdb -g -Wall -o data_gen3 data_generator_3.c -lm
	gcc -ggdb -g -Wall -o data_gen4 data_generator_4.c -lm
	gcc -ggdb -g -Wall -o defc3 defc_v3.c -lm
	gcc -ggdb -g -Wall -o defc5 defc_v5.c -lm
clean:
	rm defc_v2 fcm
