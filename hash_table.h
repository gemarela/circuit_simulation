//
//  hash_table.h
//  part4
//
//  Created by MARELAS GEORGE on 11/16/14.
//  Copyright (c) 2014 MARELAS GEORGE. All rights reserved.
//

#ifndef __part4__hash_table__
#define __part4__hash_table__

#include "utilities.h"
#define TABLESIZE 128	//Max size of the hash table


struct hashTable {
    char *nodeName;
    int key;
    struct hashTable *next;
};

struct hashTable;

struct hashTable *hashEntry;

struct hashTable *hashT[1024];


void init_hash();

int searchHashTable(char *string);

int addEntry(char *nodeName, int key);

void freeHashTable();


#endif /* defined(__part4__hash_table__) */
