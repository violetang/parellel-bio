/**
 * @file cline.h
 * @author Karin Dorman, kdorman@iastate.edu
 * @date Mon Dec 10 16:20:59 CST 2012
 *
 * Header file for cline.c.
 */

#ifndef __H_CLINE__
#define __H_CLINE__

#include <errno.h>
#include <stdio.h>
#include <string.h>

int usage_error(char **argv, int i, void *obj);
int read_int(int argc, char **argv, int i, void *obj);
unsigned int read_uint(int argc, char **argv, int i, void *obj);
long read_long(int argc, char **argv, int i, void *obj);
double read_double(int argc, char **argv, int i, void *obj);
char read_char(int argc,  char **argv, int i, void *obj);
void print_usage(const char *);
void fprint_usage(FILE *, const char *, void *obj);

#endif
