#include "ValueAndRowIndexAndColumnIndex.h"

int ValueAndRowIndexAndColumnIndex::compare
    (const void *ptr1, const void *ptr2) {
     	ValueAndRowIndexAndColumnIndex *ptrValue1 =
     	  (ValueAndRowIndexAndColumnIndex*)ptr1;
     	ValueAndRowIndexAndColumnIndex *ptrValue2 =
     	  (ValueAndRowIndexAndColumnIndex*)ptr2;
     	  
	if (ptrValue1->value == ptrValue2->value) return(0);
	if (ptrValue1->value < ptrValue2->value) return(-1);
	return(1);
}

ValueAndRowIndexAndColumnIndex::ValueAndRowIndexAndColumnIndex
    (double value, int indexRow, int indexColumn) {
	    this->value = value;
    	this->indexRow = indexRow;
    	this->indexColumn = indexColumn;
}

ValueAndRowIndexAndColumnIndex::ValueAndRowIndexAndColumnIndex() {
		this->value = -1;
    	this->indexRow = -1;
    	this->indexColumn = -1;
	
}


ValueAndRowIndexAndColumnIndex::~ValueAndRowIndexAndColumnIndex() {
}
