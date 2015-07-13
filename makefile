all:
	gcc -ggdb -g -Wall -o defc2.x defc_v2.c -lm
	gcc -ggdb -g -Wall -o fcm.x fcm_standard_paciotti_test.c -lm
	gcc -ggdb -g -Wall -o data_gen1.x data_generator.c -lm
	gcc -ggdb -g -Wall -o data_gen2.x data_generator_2.c -lm
	gcc -ggdb -g -Wall -o data_gen3.x data_generator_3.c -lm
	gcc -ggdb -g -Wall -o data_gen4.x data_generator_4.c -lm
	gcc -ggdb -g -Wall -o defc3.x defc_v3.c -lm
	gcc -ggdb -g -Wall -o defc5.x defc_v5.c -lm
	gcc -ggdb -g -Wall -o defc7.x defc_v7.c -lm
	gcc -ggdb -g -Wall -o defc9.x defc_v9.c -lm
	gcc -ggdb -g -Wall -o defc9b.x defc_v9b.c -lm
