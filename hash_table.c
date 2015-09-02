/* HASH TABLE */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hash_table.h"

/****************************************************************
 *Purpose: arxikopoiei thn hash table
 *Parametrs: none
 *Preconditions:none
 *Postconditions:tupwnei mhnuma
 *****************************************************************/
void init_hash() {
    int i;
    
    
    for(i=0; i<1024; i++) {
        hashT[i] = (struct hashTable*)malloc(sizeof(struct hashTable));
        hashT[i]->next = NULL;
    }
    printf("initialized\n");
}



/* sdbm algorithm
 * this algorithm was created for sdbm (a public-domain reimplementation of ndbm)
 * database library. it was found to do well in scrambling bits, causing
 * better distribution of the keys and fewer splits. it also happens to be
 * a good general hashing function with good distribution.
 * the actual function is hash(i) = hash(i - 1) * 65599 + str[i];
 * what is included below is the faster version used in gawk.
 * [there is even a faster, duff-device version] the magic constant 65599
 * was picked out of thin air while experimenting with different constants,
 * and turns out to be a prime. this is one of
 * the algorithms used in berkeley db (see sleepycat) and elsewhere.
 */

static unsigned long
sdbm(str)
 char *str;
{
    unsigned long hash = 0;
    int c;
    
    while ((c = *str++))
        hash = c + (hash << 6) + (hash << 16) - hash;
        
        return hash;
}


/**********************************************************************************************
 *Purpose: Anathetei kainourgies  times stous komvous.Ean to node pou dinei einai to 0 to agnoei
 den to prosthetei sto hash table kai epistrefei thn palia timh pou eixe to key,
 elegxei ean yparxei idi o komvos sto hash table kai epistrefei  thn timh
 pou tou antistixei(key) diaforetika kanei antistixish tou komvou se kainourgia timh.
 *Parametrs:nodeName:to onoma tou komvou
 key:o ari8mos pou pairnei o kombos
 *Preconditions:
 *Postconditions:epistrefei thn antistoixh timh pou dinete ston komvo
 ***********************************************************************************************/
int addEntry(char *nodeName, int key) {
    
    int i;
    
    unsigned long index;
    struct hashTable *hashNode,*tempHash;
    
    if( !strcmp(nodeName,"0") )
    {
	
        return key;
    }
    
    index = sdbm(nodeName) % 1024;
    
    
    for(tempHash=hashT[index]; tempHash->next != NULL; tempHash=tempHash->next) {
        if( !strcmp(tempHash->nodeName, nodeName) ) return key;
    }
    
    hashNode = (struct hashTable*)malloc(sizeof(struct hashTable));
    
    hashNode->nodeName = malloc((strlen(nodeName) + 1) * sizeof(char));
    strncpy(hashNode->nodeName, nodeName, strlen(nodeName) + 1);
    hashNode->key = key;
    
    hashNode->next = hashT[index];
    hashT[index] = hashNode;
    
    
    //hashArray(tempHash);
    
    ++key;
    return key;
}


int searchHashTable(char *string)
{
    int index;
    int result=-1;
    struct hashTable *hsearch;
    
    index = sdbm(string) % 1024;
    //printf("to index einai %d\n", index);
    for(hsearch=hashT[index]; hsearch->next != NULL; hsearch=hsearch->next)
    {
        if( !strcmp(hsearch->nodeName, string) )
        {
            result = hsearch->key;
            //printf("mesa sto hashTable einai %s to result %d\n", hsearch->nodeName, result);
            break;
        }
    }
    
    return result;
}


void freeHashTable()
{
    int index;
    struct hashTable *runner;
    runner=hashEntry;
    
    for(index=0; index<1024; index++)
    {
        
        runner=hashT[index];
        if(runner->nodeName!=NULL)
        {
            free(runner->nodeName);
            free(runner);
        }
        runner=runner->next;
    }
    
    free(hashEntry);
    
    return;
}