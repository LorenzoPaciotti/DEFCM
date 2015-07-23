all:
	gcc -ggdb -g -Wall -o fcm.x fcm.c -lm
	gcc -ggdb -g -Wall -o data_gen1.x data_generator.c -lm
	gcc -ggdb -g -Wall -o data_gen4.x data_generator_4.c -lm
	gcc -ggdb -g -Wall -o data_gen5.x data_generator_5.c -lm
	gcc -ggdb -g -Wall -o defc9b.x defc_v9b.c -lm
	gcc -ggdb -g -Wall -o defc10.x defc_v10.c -lm
