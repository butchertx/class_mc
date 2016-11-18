#ifndef RANDOM_H_
#define RANDOM_H_

void rand_init(int *);
void rand_init_(int *);
double drand1();
double drand1_();
void write_seeds(char seedfilename[]);
void write_seeds_(char seedfilename[]);
void read_seeds(char seedfilename[]);
void read_seeds_(char seedfilename[]);

#endif
