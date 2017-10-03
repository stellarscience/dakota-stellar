#ifndef VALUEANDROWINDEXANDCOLUMNINDEX_H_
#define VALUEANDROWINDEXANDCOLUMNINDEX_H_

#include <cstdio>
#include <cstdlib>

class ValueAndRowIndexAndColumnIndex {
	public:
	    double value;
	    int indexRow;
	    int indexColumn;
	    static int compare(const void *ptr1, const void *ptr2);
	    ValueAndRowIndexAndColumnIndex
	        (double value, 
	         int indexRow, 
	         int indexColum);
	    ValueAndRowIndexAndColumnIndex();
	    ~ValueAndRowIndexAndColumnIndex();
};

#endif /*VALUEANDINDEX_H_*/
