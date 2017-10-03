#ifndef MAINEFFECTSANALYZER3_H
#define MAINEFFECTSANALYZER3_H

#include <cstdio>
#include <cstdlib>

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <stdexcept>

#include "DataValue.h"
#include "ColumnHeader.h"


//It is hoped that, by caching computed values,
//we can dramatically shorten the execution time.
#define CACHE_COMPUTED_VALUES


class MainEffectsAnalyzer3 {
    public:

	MainEffectsAnalyzer3();
        /**
         * Construct a MainEffectsAnalyzer3.
         */
        MainEffectsAnalyzer3
            (std::vector< ColumnHeader > columnHeaders,
             std::vector< std::vector < DataValue > > data);

        /**
         * Construct a MainEffectsAnalyzer3.
         * @param columnHeader Each element in the header contains
         * information on one factor or one response.  That 
         * information includes 1) the title, or name, of the
         * factor and response and 2) the units of measure.  
         * There should one, and only one, header for each factor
         * (even if the raw input data file contains multiple
         * columns for the same factor).
         * There should one, and only one, header for each response
         * (even if the raw input data file contains multiple
         * columns for the same response).
         * @param lineNumberOfFirstDataRow In a typical data
         * file, the first few rows of the file contain header
         * information.  How many rows of header information does
         * this file have?  EXAMPLE:  Suppose, in our data file,
         * the first row contains the names of the columns, 
         * the second row contains the units of measure,
         * the third row is blank,
         * and the fourth row contains the first line of data 
         * then we would set lineNumberOfFirstDataRow to 3.          
         * param mapFilenameColumnsToColumnHeaders An index array
         * that maps columns in the data file to column header
         * entries.  EXAMPLE;  The array {3,0,1,2} 
         * maps the 1st column in the data file to columnHeaders[3],  
         * maps the 2nd column in the data file to columnHeaders[0],
         * maps the 3rd column in the data file to columnHeaders[1], 
         * maps the 4th column in the data file to columnHeaders[2].
         */
         void importCvsFile
            (std::vector< ColumnHeader > columnHeaders,
             std::string filename,
             int lineNumberOfFirstDataRow,
             int mapFilenameColumnsToColumnHeaders[]);

        /* Destructor */
        virtual ~MainEffectsAnalyzer3();


        /**
         * Given an indentifier of a data column in the data array,
         * return the index of the column.
         * <p>
         * The input can be either: <br>
         * &nbsp;&nbsp;&nbsp;The name of the column <br>
         * &nbsp;&nbsp;&nbsp;The capital letters of the col name <br>
         * &nbsp;&nbsp;&nbsp;"a" or "b" or "c" or ... <br>
         * &nbsp;&nbsp;&nbsp;"A" or "B" or "C" or ... <br>
         * <p>
         * EXAMPLES: <br>
         * if the data array has the following 3 columns: <br>
         * &nbsp;&nbsp;&nbsp;&nbsp;OneTwoThree<br>
         * &nbsp;&nbsp;&nbsp;&nbsp;FourFiveSix<br>
         * &nbsp;&nbsp;&nbsp;&nbsp;SevenEightNine<br>
         * then <br>
         * "OneTwoThree" returns the index value 0 <br>
         * "OTT" returns the index value 0 <br>
         * "a" returns the index value 0 <br>
         * "A" returns the index value 0 <br>
         * "DoesNotExist" throws an exception <br>
         * @param columnIdentifier 
         * A string that identifies a column.
         * @return The index to the column in the data array.
         * @throws NoSuchColumnException Thrown if the column
         * can not be found in the data array.
         */
        int toIndexInputColumn(std::string columnId);



        /**
         * Get all of the unique factors that are in
         * the targeted factor column.
         * <p>
         * The unique factors are returned as a set
         * of DataValues.
         * <p>
         * EXAMPLE <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 22 23 <br>
         * &nbsp;&nbsp; 32 33 <br>
         * &nbsp;&nbsp; 42 43 <br>
         * 51 52 &nbsp;&nbsp; <br>
         * 51 62 &nbsp;&nbsp; <br>
         * and if the target column is the first column <br>
         * then we would return the set {DataValue(11), DataValue(51)}.
         */
        virtual std::vector< DataValue > getNonEmptyUniqueFactors
            (int indexFactorColumn);


        /**
         * Is the DataValue already in the set?
         * @param dataValue Is this DataValue in set?
         * @param set Does this set contain the dataValue?
         */
        virtual bool isDataValueInSet
               (DataValue dataValue,
                std::vector< DataValue > set);


        /**
         * Count the number of rows in the data array
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would return the value 1 <br>
         * because only one row has entries in both the 
         * first and last column.
         * @param indexFactorColumn The data array index of a column.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual int getNumberOfObservations
            (int indexFactorColumn,
             int indexResponseColumn);



        /**
         * Count the number of rows in the data array
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would return the value 1 <br>
         * because only one row has entries in both the 
         * first and last column.
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError

         */          
        virtual int getNumberOfObservations
            (std::string factorColumn,
             std::string responseColumn);

        /**
         * Count the number of rows in the data array
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would return the value 1 <br>
         * because only one row has entries in both the 
         * first and last column.
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual int getNumberOfObservations
            (std::string factorColumn,
             int indexResponseColumn);


        /**
         * Count the number of rows in the data array
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would return the value 1 <br>
         * because only one row has entries in both the 
         * first and last column.
         * @param indexFactorColumn The data array index of a column.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual int getNumberOfObservations
            (int indexFactorColumn,
             std::string responseColumn);



        /**
         * Count the number of rows in the data array
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 2 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 1) have a non-empty value
         * in the response column.
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual int getNumberOfObservations
            (std::string factorColumn,
             DataValue factorValue,
             std::string responseColumn);


        /**
         * Count the number of rows in the data array
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 2 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual int getNumberOfObservations
            (std::string factorColumn,
             DataValue factorValue,
             int indexResponseColumn);




        /**
         * Count the number of rows in the data array
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 2 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual int getNumberOfObservations
            (int indexFactorColumn,
             DataValue factorValue,
             std::string responseColumn);



        /**
         * Count the number of rows in the data array
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 2 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual int getNumberOfObservations
            (int indexFactorColumn,
             DataValue factorValue,
             int indexResponseColumn);



        /**
         * Count the number of rows in the data array
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 2 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */        
        virtual int getNumberOfObservations
            (std::string factorColumn,
             std::string factorValue,
             std::string responseColumn);



        /**
         * Count the number of rows in the data array
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 2 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual int getNumberOfObservations
            (int indexFactorColumn,
             std::string factorValue,
             int indexResponseColumn);

        virtual int getNumberOfObservations
            (std::string factorColumn,
             std::string factorValue,
             int indexResponseColumn);

        /**
         * Count the number of rows in the data array
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 2 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual int getNumberOfObservations
            (int indexFactorColumn,
             std::string factorValue,
             std::string responseColumn);


        /**
         * Count the number of rows in the data array
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 2 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */        
        virtual int getNumberOfObservations
            (std::string factorColumn,
             double factorValue,
             std::string responseColumn);


        /**
         * Count the number of rows in the data array
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 2 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual int getNumberOfObservations
            (std::string factorColumn,
             double factorValue,
             int indexResponseColumn);

        /**
         * Count the number of rows in the data array
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 2 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual int getNumberOfObservations
            (int indexFactorColumn,
             double factorValue,
             std::string responseColumn);


        /**
         * Count the number of rows in the data array
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 2 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual int getNumberOfObservations
            (int indexFactorColumn,
             double factorValue,
             int indexResponseColumn);

        /**
         * Count the number of rows in the data array
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 2 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */        
        virtual int getNumberOfObservations
            (std::string factorColumn,
             int factorValue,
             std::string responseColumn);


        /**
         * Count the number of rows in the data array
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 2 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual int getNumberOfObservations
            (std::string factorColumn,
             int factorValue,
             int indexResponseColumn);

        /**
         * Count the number of rows in the data array
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 2 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual int getNumberOfObservations
            (int indexFactorColumn,
             int factorValue,
             std::string responseColumn);

        /**
         * Count the number of rows in the data array
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 2 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual int getNumberOfObservations
            (int indexFactorColumn,
             int factorValue,
             int indexResponseColumn);



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/






        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would return the value 382 <br>
         * (which is the sum of 13+113+123+133). <br>
         * @param indexFactorColumn The data array index of a column.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (int indexFactorColumn,
             int indexResponseColumn);




        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would return the value 382 <br>
         * (which is the sum of 13+113+123+133). <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (std::string factorColumn,
             std::string responseColumn);



        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would return the value 382 <br>
         * (which is the sum of 13+113+123+133). <br>
         * @param indexFactorColumn The data array index of a column.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (int indexFactorColumn,
             std::string responseColumn);


        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would return the value 382 <br>
         * (which is the sum of 13+113+123+133). <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (std::string factorColumn,
             int indexResponseColumn);


        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor entry
         * matches the desired value, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 126 <br>
         * (which is the sum of 13+113). <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (int indexFactorColumn,
             DataValue factorValue,
             int indexResponseColumn);


        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor entry
         * matches the desired value, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 126 <br>
         * (which is the sum of 13+113). <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (std::string factorColumn,
             DataValue factorValue,
             std::string responseColumn);



        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor entry
         * matches the desired value, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 126 <br>
         * (which is the sum of 13+113). <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (int indexFactorColumn,
             DataValue factorValue,
             std::string responseColumn);


        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor entry
         * matches the desired value, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 126 <br>
         * (which is the sum of 13+113). <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (std::string factorColumn,
             DataValue factorValue,
             int indexResponseColumn);


        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor entry
         * matches the desired value, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 126 <br>
         * (which is the sum of 13+113). <br>
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (int indexFactorColumn,
             double factorValue,
             int indexResponseColumn);


        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor entry
         * matches the desired value, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 126 <br>
         * (which is the sum of 13+113). <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (std::string factorColumn,
             double factorValue,
             std::string responseColumn);



        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor entry
         * matches the desired value, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 126 <br>
         * (which is the sum of 13+113). <br>
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (int indexFactorColumn,
             double factorValue,
             std::string responseColumn);


        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor entry
         * matches the desired value, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 126 <br>
         * (which is the sum of 13+113). <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (std::string factorColumn,
             double factorValue,
             int indexResponseColumn);





        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor entry
         * matches the desired value, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 126 <br>
         * (which is the sum of 13+113). <br>
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (int indexFactorColumn,
             int factorValue,
             int indexResponseColumn);


        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor entry
         * matches the desired value, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 126 <br>
         * (which is the sum of 13+113). <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (std::string factorColumn,
             int factorValue,
             std::string responseColumn);



        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor entry
         * matches the desired value, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 126 <br>
         * (which is the sum of 13+113). <br>
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (int indexFactorColumn,
             int factorValue,
             std::string responseColumn);


        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor entry
         * matches the desired value, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 126 <br>
         * (which is the sum of 13+113). <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (std::string factorColumn,
             int factorValue,
             int indexResponseColumn);



        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor entry
         * matches the desired value, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 126 <br>
         * (which is the sum of 13+113). <br>
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (int indexFactorColumn,
             std::string factorValue,
             int indexResponseColumn);


        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor entry
         * matches the desired value, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 126 <br>
         * (which is the sum of 13+113). <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (std::string factorColumn,
             std::string factorValue,
             std::string responseColumn);



        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor entry
         * matches the desired value, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 126 <br>
         * (which is the sum of 13+113). <br>
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (int indexFactorColumn,
             std::string factorValue,
             std::string responseColumn);


        /**
         * Compute the sum of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor entry
         * matches the desired value, then the value of 
         * the response is added to the sum.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 126 <br>
         * (which is the sum of 13+113). <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual double getSumOfObservations
            (std::string factorColumn,
             std::string factorValue,
             int indexResponseColumn);

    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/






        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would return the value 127.3333 <br>
         * (which is the average of 13,113,123,133). <br>
         * @param indexFactorColumn The data array index of a column.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (int indexFactorColumn,
             int indexResponseColumn);







        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would return the value 127.3333 <br>
         * (which is the average of 13,113,123,133). <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (std::string factorColumn,
             std::string responseColumn);







        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would return the value 127.3333 <br>
         * (which is the average of 13,113,123,133). <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (std::string factorColumn,
             int indexResponseColumn);







        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would return the value 127.3333 <br>
         * (which is the average of 13,113,123,133). <br>
         * @param indexFactorColumn The data array index of a column.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (int indexFactorColumn,
             std::string responseColumn);







        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor
         * entry equals the desired value, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 63 <br>
         * (which is the average of 13,113). <br>
         * @param indexFactoryColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (int indexFactorColumn,
             DataValue factorValue,
             int indexResponseColumn);


        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor
         * entry equals the desired value, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 63 <br>
         * (which is the average of 13,113). <br>
         * @param indexFactorColumn The data array index of a column.
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).         
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (std::string factorColumn,
             DataValue factorValue,
             std::string responseColumn);



        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor
         * entry equals the desired value, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 63 <br>
         * (which is the average of 13,113). <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (std::string factorColumn,
             DataValue factorValue,
             int indexResponseColumn);


        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor
         * entry equals the desired value, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 63 <br>
         * (which is the average of 13,113). <br>
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).         
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (int indexFactorColumn,
             DataValue factorValue,
             std::string responseColumn);



        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor
         * entry equals the desired value, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 63 <br>
         * (which is the average of 13,113). <br>
         * @param indexFactoryColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (int indexFactorColumn,
             double factorValue,
             int indexResponseColumn);


        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor
         * entry equals the desired value, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 63 <br>
         * (which is the average of 13,113). <br>
         * @param indexFactorColumn The data array index of a column.
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).         
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (std::string factorColumn,
             double factorValue,
             std::string responseColumn);



        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor
         * entry equals the desired value, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 63 <br>
         * (which is the average of 13,113). <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (std::string factorColumn,
             double factorValue,
             int indexResponseColumn);


        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor
         * entry equals the desired value, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 63 <br>
         * (which is the average of 13,113). <br>
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).         
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (int indexFactorColumn,
             double factorValue,
             std::string responseColumn);



        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor
         * entry equals the desired value, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 63 <br>
         * (which is the average of 13,113). <br>
         * @param indexFactoryColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (int indexFactorColumn,
             int factorValue,
             int indexResponseColumn);


        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor
         * entry equals the desired value, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 63 <br>
         * (which is the average of 13,113). <br>
         * @param indexFactorColumn The data array index of a column.
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).         
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (std::string factorColumn,
             int factorValue,
             std::string responseColumn);



        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor
         * entry equals the desired value, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 63 <br>
         * (which is the average of 13,113). <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (std::string factorColumn,
             int factorValue,
             int indexResponseColumn);


        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor
         * entry equals the desired value, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 63 <br>
         * (which is the average of 13,113). <br>
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).         
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (int indexFactorColumn,
             int factorValue,
             std::string responseColumn);






        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor
         * entry equals the desired value, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 63 <br>
         * (which is the average of 13,113). <br>
         * @param indexFactoryColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (int indexFactorColumn,
             std::string factorValue,
             int indexResponseColumn);


        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor
         * entry equals the desired value, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 63 <br>
         * (which is the average of 13,113). <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).         
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (std::string factorColumn,
             std::string factorValue,
             std::string responseColumn);



        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor
         * entry equals the desired value, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 63 <br>
         * (which is the average of 13,113). <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (std::string factorColumn,
             std::string factorValue,
             int indexResponseColumn);


        /**
         * Compute the average of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor
         * entry equals the desired value, then the value of 
         * the response is included in the average.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 63 <br>
         * (which is the average of 13,113). <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).         
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getAverageObservation
            (int indexFactorColumn,
             std::string factorValue,
             std::string responseColumn);





    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/






        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * sum = 13 + 113 + 123 + 133 = 382 <br>
         * numberOfValues = 4 <br>
         * avg = sum/numberOfValues =  95.5 <br>
         * sumOfSquares = (13-95.5)^2 + (113-95.5)^2 <br>
         * + (123-95.5)^2 + (133-95.5)^2 = 9275 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (int indexFactorColumn,
             int indexResponseColumn);



        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * sum = 13 + 113 + 123 + 133 = 382 <br>
         * numberOfValues = 4 <br>
         * avg = sum/numberOfValues =  95.5 <br>
         * sumOfSquares = (13-95.5)^2 + (113-95.5)^2 <br>
         * + (123-95.5)^2 + (133-95.5)^2 = 9275 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).      
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (std::string factorColumn,
             std::string responseColumn);



        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * sum = 13 + 113 + 123 + 133 = 382 <br>
         * numberOfValues = 4 <br>
         * avg = sum/numberOfValues =  95.5 <br>
         * sumOfSquares = (13-95.5)^2 + (113-95.5)^2 <br>
         * + (123-95.5)^2 + (133-95.5)^2 = 9275 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (std::string factorColumn,
             int indexResponseColumn);


        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * sum = 13 + 113 + 123 + 133 = 382 <br>
         * numberOfValues = 4 <br>
         * avg = sum/numberOfValues =  95.5 <br>
         * sumOfSquares = (13-95.5)^2 + (113-95.5)^2 <br>
         * + (123-95.5)^2 + (133-95.5)^2 = 9275 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).      
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (int indexFactorColumn,
             std::string responseColumn);




        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * matches the desired value, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6344 <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (int indexFactorColumn,
             DataValue factorValue,
             int indexResponseColumn);



        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * matches the desired value, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6344 <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).       
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (std::string factorColumn,
             DataValue factorValue,
             std::string responseColumn);



        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * matches the desired value, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6344 <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (std::string factorColumn,
             DataValue factorValue,
             int indexResponseColumn);



        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * matches the desired value, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6344 <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).       
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (int indexFactorColumn,
             DataValue factorValue,
             std::string responseColumn);




        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * matches the desired value, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6344 <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (int indexFactorColumn,
             double factorValue,
             int indexResponseColumn);



        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * matches the desired value, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6344 <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).       
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (std::string factorColumn,
             double factorValue,
             std::string responseColumn);



        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * matches the desired value, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6344 <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (std::string factorColumn,
             double factorValue,
             int indexResponseColumn);



        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * matches the desired value, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6344 <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).       
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (int indexFactorColumn,
             double factorValue,
             std::string responseColumn);






        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * matches the desired value, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6344 <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (int indexFactorColumn,
             int factorValue,
             int indexResponseColumn);



        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * matches the desired value, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6344 <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).       
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (std::string factorColumn,
             int factorValue,
             std::string responseColumn);



        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * matches the desired value, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6344 <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (std::string factorColumn,
             int factorValue,
             int indexResponseColumn);



        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * matches the desired value, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6344 <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).       
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (int indexFactorColumn,
             int factorValue,
             std::string responseColumn);






        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * matches the desired value, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6344 <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (int indexFactorColumn,
             std::string factorValue,
             int indexResponseColumn);



        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * matches the desired value, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6344 <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).       
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (std::string factorColumn,
             std::string factorValue,
             std::string responseColumn);



        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * matches the desired value, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6344 <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (std::string factorColumn,
             std::string factorValue,
             int indexResponseColumn);



        /**
         * Compute the sum-of-squares of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * matches the desired value, then the value of 
         * the response is added to the sum-of-squares.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6344 <br>
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).       
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getSumOfSquares
            (int indexFactorColumn,
             std::string factorValue,
             std::string responseColumn);





    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/






        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * sum = 13 + 113 + 123 + 133 = 382 <br>
         * numberOfValues = 4 <br>
         * avg = sum/numberOfValues =  95.5 <br>
         * sumOfSquares = (13-95.5)^2 + (113-95.5)^2 <br>
         * + (123-95.5)^2 + (133-95.5)^2 = 9275 <br>
         * variance = 9275/3 = 3092 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param indexFactorColumn The data array index of a column.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (int indexFactorColumn,
             int indexResponseColumn);







        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * sum = 13 + 113 + 123 + 133 = 382 <br>
         * numberOfValues = 4 <br>
         * avg = sum/numberOfValues =  95.5 <br>
         * sumOfSquares = (13-95.5)^2 + (113-95.5)^2 <br>
         * + (123-95.5)^2 + (133-95.5)^2 = 9275 <br>
         * variance = 9275/3 = 3092 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).    
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (std::string factorColumn,
             std::string responseColumn);






        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * sum = 13 + 113 + 123 + 133 = 382 <br>
         * numberOfValues = 4 <br>
         * avg = sum/numberOfValues =  95.5 <br>
         * sumOfSquares = (13-95.5)^2 + (113-95.5)^2 <br>
         * + (123-95.5)^2 + (133-95.5)^2 = 9275 <br>
         * variance = 9275/3 = 3092 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (std::string factorColumn,
             int indexResponseColumn);








        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * sum = 13 + 113 + 123 + 133 = 382 <br>
         * numberOfValues = 4 <br>
         * avg = sum/numberOfValues =  95.5 <br>
         * sumOfSquares = (13-95.5)^2 + (113-95.5)^2 <br>
         * + (123-95.5)^2 + (133-95.5)^2 = 9275 <br>
         * variance = 9275/3 = 3092 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).    
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (int indexFactorColumn,
             std::string responseColumn);





        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * equals the desired value, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6100 <br>
         * variance = 6100/1 = 6100 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (int indexFactorColumn,
             DataValue factorValue,
             int indexResponseColumn);


        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * equals the desired value, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6100 <br>
         * variance = 6100/1 = 6100 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).    
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (std::string factorColumn,
             DataValue factorValue,
             std::string responseColumn);



        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * equals the desired value, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6100 <br>
         * variance = 6100/1 = 6100 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).    
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (int indexFactorColumn,
             DataValue factorValue,
             std::string responseColumn);



        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * equals the desired value, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6100 <br>
         * variance = 6100/1 = 6100 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (std::string factorColumn,
             DataValue factorValue,
             int indexResponseColumn);



        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * equals the desired value, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6100 <br>
         * variance = 6100/1 = 6100 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (int indexFactorColumn,
             double factorValue,
             int indexResponseColumn);


        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * equals the desired value, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6100 <br>
         * variance = 6100/1 = 6100 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).    
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (std::string factorColumn,
             double factorValue,
             std::string responseColumn);



        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * equals the desired value, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6100 <br>
         * variance = 6100/1 = 6100 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).    
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (int indexFactorColumn,
             double factorValue,
             std::string responseColumn);



        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * equals the desired value, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6100 <br>
         * variance = 6100/1 = 6100 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (std::string factorColumn,
             double factorValue,
             int indexResponseColumn);




        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * equals the desired value, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6100 <br>
         * variance = 6100/1 = 6100 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (int indexFactorColumn,
             int factorValue,
             int indexResponseColumn);


        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * equals the desired value, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6100 <br>
         * variance = 6100/1 = 6100 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).    
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (std::string factorColumn,
             int factorValue,
             std::string responseColumn);



        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * equals the desired value, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6100 <br>
         * variance = 6100/1 = 6100 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).    
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (int indexFactorColumn,
             int factorValue,
             std::string responseColumn);



        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * equals the desired value, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6100 <br>
         * variance = 6100/1 = 6100 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (std::string factorColumn,
             int factorValue,
             int indexResponseColumn);






        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * equals the desired value, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6100 <br>
         * variance = 6100/1 = 6100 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (int indexFactorColumn,
             std::string factorValue,
             int indexResponseColumn);


        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * equals the desired value, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6100 <br>
         * variance = 6100/1 = 6100 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).    
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (std::string factorColumn,
             std::string factorValue,
             std::string responseColumn);



        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * equals the desired value, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6100 <br>
         * variance = 6100/1 = 6100 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).    
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (int indexFactorColumn,
             std::string factorValue,
             std::string responseColumn);



        /**
         * Compute the variance of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column and if the factor value
         * equals the desired value, then the value of 
         * the response is included in the variance.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would use the first and second
         * rows. <br>
         * sum = 13 + 113 = 126 <br>
         * numberOfValues = 2 <br>
         * avg = sum/numberOfValues =  63 <br>
         * sumOfSquares = (13-63)^2 + (113-63)^2 = 6100 <br>
         * variance = 6100/1 = 6100 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */          
        virtual double getVariance
            (std::string factorColumn,
             std::string factorValue,
             int indexResponseColumn);




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



        /**
         * Compute d.  d is numberOfObservations - 1.
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would return the value 0 <br>
         * because only one row has entries in both the 
         * first and last column.
         * @param indexFactorColumn The data array index of a column.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual int getD
            (int indexFactorColumn,
             int indexResponseColumn);



        /**
         * Compute d.  d is numberOfObservations - 1.
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would return the value 0 <br>
         * because only one row has entries in both the 
         * first and last column.
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError

         */          
        virtual int getD
            (std::string factorColumn,
             std::string responseColumn);

        /**
         * Compute d.  d is numberOfObservations - 1.
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would return the value 0 <br>
         * because only one row has entries in both the 
         * first and last column.
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual int getD
            (std::string factorColumn,
             int indexResponseColumn);


        /**
         * Compute d.  d is numberOfObservations - 1.
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would return the value 0 <br>
         * because only one row has entries in both the 
         * first and last column.
         * @param indexFactorColumn The data array index of a column.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual int getD
            (int indexFactorColumn,
             std::string responseColumn);



        /**
         * Compute d.  d is numberOfObservations - 1.
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 1 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 1) have a non-empty value
         * in the response column.
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual int getD
            (std::string factorColumn,
             DataValue factorValue,
             std::string responseColumn);







        /**
         * Compute d.  d is numberOfObservations - 1.
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 1 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual int getD
            (std::string factorColumn,
             DataValue factorValue,
             int indexResponseColumn);




        /**
         * Compute d.  d is numberOfObservations - 1.
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 1 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual int getD
            (int indexFactorColumn,
             DataValue factorValue,
             std::string responseColumn);



        /**
         * Compute d.  d is numberOfObservations - 1.
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 1 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * <p>
         * How do we determine if the data array's factor 
         * value matches the desired value?  If the data array's
         * factor value is stored as a double, then the value
         * of that double is compared with the double value that
         * is stored inside of the factorValue parameter (the
         * dataType of the factorValue does not have to be
         * set to DOUBLE).  If the data array's factor
         * value is stored as an int, then the value of that int
         * is compared with the int value that is stored inside
         * of the factorValue parameter (the dataType of the
         * factorValue does not have to be set to INTEGER).  
         * If the data array's factor value is stored as a string,
         * then the value that string is compared with the string
         * value that is stored inside of the factorValue
         * parameter (the dataType of the factorValue does not 
         * have to be set to STRING).
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual int getD
            (int indexFactorColumn,
             DataValue factorValue,
             int indexResponseColumn);



        /**
         * Compute d.  d is numberOfObservations - 1.
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 1 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */        
        virtual int getD
            (std::string factorColumn,
             std::string factorValue,
             std::string responseColumn);



        /**
         * Compute d.  d is numberOfObservations - 1.
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 1 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual int getD
            (int indexFactorColumn,
             std::string factorValue,
             int indexResponseColumn);

        virtual int getD
            (std::string factorColumn,
             std::string factorValue,
             int indexResponseColumn);

        /**
         * Compute d.  d is numberOfObservations - 1.
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 1 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual int getD
            (int indexFactorColumn,
             std::string factorValue,
             std::string responseColumn);


        /**
         * Compute d.  d is numberOfObservations - 1.
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 1 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */        
        virtual int getD
            (std::string factorColumn,
             double factorValue,
             std::string responseColumn);


        /**
         * Compute d.  d is numberOfObservations - 1.
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 1 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual int getD
            (std::string factorColumn,
             double factorValue,
             int indexResponseColumn);

        /**
         * Compute d.  d is numberOfObservations - 1.
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 1 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual int getD
            (int indexFactorColumn,
             double factorValue,
             std::string responseColumn);


        /**
         * Compute d.  d is numberOfObservations - 1.
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 1 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual int getD
            (int indexFactorColumn,
             double factorValue,
             int indexResponseColumn);

        /**
         * Compute d.  d is numberOfObservations - 1.
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 1 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */        
        virtual int getD
            (std::string factorColumn,
             int factorValue,
             std::string responseColumn);


        /**
         * Compute d.  d is numberOfObservations - 1.
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 1 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual int getD
            (std::string factorColumn,
             int factorValue,
             int indexResponseColumn);

        /**
         * Compute d.  d is numberOfObservations - 1.
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 1 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         */          
        virtual int getD
            (int indexFactorColumn,
             int factorValue,
             std::string responseColumn);

        /**
         * Compute d.  d is numberOfObservations - 1.
         * where 1) the factor column has an entry in the row
         * AND 2) the response column has an entry in the row
         * AND 3) the value of the factor matches the desired
         * value.
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12  13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 31 32 &nbsp;&nbsp; <br>
         * 31 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * and if the desired factor value is 11 <br>
         * then we would return the value 1 <br>
         * because only two rows 1) have the desired value
         * in the factor column and 2) have a non-empty value
         * in the response column.
         * @param indexFactorColumn The data array index of a column.
         * @param factorValue the disired factor value.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         */          
        virtual int getD
            (int indexFactorColumn,
             int factorValue,
             int indexResponseColumn);


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



        /**
         * Compute the sum-of-squares-between-groups 
         * of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * sum-of-squares-between-groups computation.
         * <p>
         * sum-of-squares-between-groups = <br>
         *   SUM(numberOfObservationsInGroup*(avgGroup-avgTotal)^2) <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * totalSum = 13 + 113 + 123 + 133 = 382 <br>
         * totalNumberOfValues = 4 <br>
         * totalAverage = totalSum/totalNumberOfValues = 382/4 = 95.5 <br>
         * group 1 contains all rows with a factor value of 11. <br>
         * &nbsp;&nbsp;&nbsp;group 1 contains rows 0 and 1 <br>
         * &nbsp;&nbsp;&nbsp;group1NumberOfObserations = 2 <br>
         * &nbsp;&nbsp;&nbsp;group1Sum = 13+113 = 126; <br>
         * &nbsp;&nbsp;&nbsp;group1Average = sum/count = 126/2 = 63 <br>
         * group 2 contains all rows with a factor value of 21. <br>
         * &nbsp;&nbsp;&nbsp;group 2 contains rows 3 and 5 <br>
         * &nbsp;&nbsp;&nbsp;group2NumberOfObservations = 2 <br>
         * &nbsp;&nbsp;&nbsp;group2Sum = 123+133 = 256 <br>
         * &nbsp;&nbsp;&nbsp;group2Average = sum/count = 256/2 = 128 <br>
         * sumOfSquaresBetweenGroups = <br>
         * &nbsp;&nbsp;&nbsp;group1Count(group1Avg-totalAvg)^2 +<br>
         * &nbsp;&nbsp;&nbsp;group2count(group2Avg-totalAvg)^2 <br>
         * &nbsp;&nbsp;&nbsp;= 2*(63-95.5)^2 + 2*(128-95.5)^2 <br>
         * &nbsp;&nbsp;&nbsp;= 4225 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */ 
        virtual double getSumOfSquaresBetweenGroups
                (int indexFactorColumn,
                 int indexResponseColumn);



        /**
         * Compute the sum-of-squares-between-groups 
         * of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * sum-of-squares-between-groups computation.
         * <p>
         * sum-of-squares-between-groups = <br>
         *   SUM(numberOfObservationsInGroup*(avgGroup-avgTotal)^2) <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * totalSum = 13 + 113 + 123 + 133 = 382 <br>
         * totalNumberOfValues = 4 <br>
         * totalAverage = totalSum/totalNumberOfValues = 382/4 = 95.5 <br>
         * group 1 contains all rows with a factor value of 11. <br>
         * &nbsp;&nbsp;&nbsp;group 1 contains rows 0 and 1 <br>
         * &nbsp;&nbsp;&nbsp;group1NumberOfObserations = 2 <br>
         * &nbsp;&nbsp;&nbsp;group1Sum = 13+113 = 126; <br>
         * &nbsp;&nbsp;&nbsp;group1Average = sum/count = 126/2 = 63 <br>
         * group 2 contains all rows with a factor value of 21. <br>
         * &nbsp;&nbsp;&nbsp;group 2 contains rows 3 and 5 <br>
         * &nbsp;&nbsp;&nbsp;group2NumberOfObservations = 2 <br>
         * &nbsp;&nbsp;&nbsp;group2Sum = 123+133 = 256 <br>
         * &nbsp;&nbsp;&nbsp;group2Average = sum/count = 256/2 = 128 <br>
         * sumOfSquaresBetweenGroups = <br>
         * &nbsp;&nbsp;&nbsp;group1Count(group1Avg-totalAvg)^2 +<br>
         * &nbsp;&nbsp;&nbsp;group2count(group2Avg-totalAvg)^2 <br>
         * &nbsp;&nbsp;&nbsp;= 2*(63-95.5)^2 + 2*(128-95.5)^2 <br>
         * &nbsp;&nbsp;&nbsp;= 4225 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */ 
        virtual double getSumOfSquaresBetweenGroups
                (std::string factorColumn,
                 std::string responseColumn);



        /**
         * Compute the sum-of-squares-between-groups 
         * of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * sum-of-squares-between-groups computation.
         * <p>
         * sum-of-squares-between-groups = <br>
         *   SUM(numberOfObservationsInGroup*(avgGroup-avgTotal)^2) <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * totalSum = 13 + 113 + 123 + 133 = 382 <br>
         * totalNumberOfValues = 4 <br>
         * totalAverage = totalSum/totalNumberOfValues = 382/4 = 95.5 <br>
         * group 1 contains all rows with a factor value of 11. <br>
         * &nbsp;&nbsp;&nbsp;group 1 contains rows 0 and 1 <br>
         * &nbsp;&nbsp;&nbsp;group1NumberOfObserations = 2 <br>
         * &nbsp;&nbsp;&nbsp;group1Sum = 13+113 = 126; <br>
         * &nbsp;&nbsp;&nbsp;group1Average = sum/count = 126/2 = 63 <br>
         * group 2 contains all rows with a factor value of 21. <br>
         * &nbsp;&nbsp;&nbsp;group 2 contains rows 3 and 5 <br>
         * &nbsp;&nbsp;&nbsp;group2NumberOfObservations = 2 <br>
         * &nbsp;&nbsp;&nbsp;group2Sum = 123+133 = 256 <br>
         * &nbsp;&nbsp;&nbsp;group2Average = sum/count = 256/2 = 128 <br>
         * sumOfSquaresBetweenGroups = <br>
         * &nbsp;&nbsp;&nbsp;group1Count(group1Avg-totalAvg)^2 +<br>
         * &nbsp;&nbsp;&nbsp;group2count(group2Avg-totalAvg)^2 <br>
         * &nbsp;&nbsp;&nbsp;= 2*(63-95.5)^2 + 2*(128-95.5)^2 <br>
         * &nbsp;&nbsp;&nbsp;= 4225 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */ 
        virtual double getSumOfSquaresBetweenGroups
                (std::string factorColumn,
                 int indexResponseColumn);



        /**
         * Compute the sum-of-squares-between-groups 
         * of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * sum-of-squares-between-groups computation.
         * <p>
         * sum-of-squares-between-groups = <br>
         *   SUM(numberOfObservationsInGroup*(avgGroup-avgTotal)^2) <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * totalSum = 13 + 113 + 123 + 133 = 382 <br>
         * totalNumberOfValues = 4 <br>
         * totalAverage = totalSum/totalNumberOfValues = 382/4 = 95.5 <br>
         * group 1 contains all rows with a factor value of 11. <br>
         * &nbsp;&nbsp;&nbsp;group 1 contains rows 0 and 1 <br>
         * &nbsp;&nbsp;&nbsp;group1NumberOfObserations = 2 <br>
         * &nbsp;&nbsp;&nbsp;group1Sum = 13+113 = 126; <br>
         * &nbsp;&nbsp;&nbsp;group1Average = sum/count = 126/2 = 63 <br>
         * group 2 contains all rows with a factor value of 21. <br>
         * &nbsp;&nbsp;&nbsp;group 2 contains rows 3 and 5 <br>
         * &nbsp;&nbsp;&nbsp;group2NumberOfObservations = 2 <br>
         * &nbsp;&nbsp;&nbsp;group2Sum = 123+133 = 256 <br>
         * &nbsp;&nbsp;&nbsp;group2Average = sum/count = 256/2 = 128 <br>
         * sumOfSquaresBetweenGroups = <br>
         * &nbsp;&nbsp;&nbsp;group1Count(group1Avg-totalAvg)^2 +<br>
         * &nbsp;&nbsp;&nbsp;group2count(group2Avg-totalAvg)^2 <br>
         * &nbsp;&nbsp;&nbsp;= 2*(63-95.5)^2 + 2*(128-95.5)^2 <br>
         * &nbsp;&nbsp;&nbsp;= 4225 <br>
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */ 
        virtual double getSumOfSquaresBetweenGroups
                (int indexFactorColumn,
                 std::string responseColumn);


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



        /**
         * Compute the sum-of-squares-within-groups 
         * of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * sum-of-squares-within-groups computation.
         * <p>
         * sum-of-squares-within-groups = <br>
         *   SUM(GroupSumOfSquares) <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * totalSum = 13 + 113 + 123 + 133 = 382 <br>
         * totalNumberOfValues = 4 <br>
         * totalAverage = totalSum/totalNumberOfValues = 382/4 = 95.5 <br>
         * group 1 contains all rows with a factor value of 11. <br>
         * &nbsp;&nbsp;&nbsp;group 1 contains rows 0 and 1 <br>
         * &nbsp;&nbsp;&nbsp;group1NumberOfObserations = 2 <br>
         * &nbsp;&nbsp;&nbsp;group1Sum = 13+113 = 126; <br>
         * &nbsp;&nbsp;&nbsp;group1Average = sum/count = 126/2 = 63 <br>
         * &nbsp;&nbsp;&nbsp;group1SumOfSquares = <br>
         * &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(13-63)^2+(113-63)^2 = 5000 <br>
         * group 2 contains all rows with a factor value of 21. <br>
         * &nbsp;&nbsp;&nbsp;group 2 contains rows 3 and 5 <br>
         * &nbsp;&nbsp;&nbsp;group2NumberOfObservations = 2 <br>
         * &nbsp;&nbsp;&nbsp;group2Sum = 123+133 = 256 <br>
         * &nbsp;&nbsp;&nbsp;group2Average = sum/count = 256/2 = 128 <br>
         * &nbsp;&nbsp;&nbsp;group2SumOfSquares = <br>
         * &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(123-128)^2+(133-128)^2=50<br>
         * sumOfSquaresWithinGroups = <br>
         * &nbsp;&nbsp;&nbsp;group1sumOfSquares+group2SumOfSquares <br>
         * &nbsp;&nbsp;&nbsp;= 5000 + 50 = 5050 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */ 
        virtual double getSumOfSquaresWithinGroups
                (int indexFactorColumn,
                 int indexResponseColumn);




        /**
         * Compute the sum-of-squares-within-groups 
         * of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * sum-of-squares-within-groups computation.
         * <p>
         * sum-of-squares-within-groups = <br>
         *   SUM(GroupSumOfSquares) <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * totalSum = 13 + 113 + 123 + 133 = 382 <br>
         * totalNumberOfValues = 4 <br>
         * totalAverage = totalSum/totalNumberOfValues = 382/4 = 95.5 <br>
         * group 1 contains all rows with a factor value of 11. <br>
         * &nbsp;&nbsp;&nbsp;group 1 contains rows 0 and 1 <br>
         * &nbsp;&nbsp;&nbsp;group1NumberOfObserations = 2 <br>
         * &nbsp;&nbsp;&nbsp;group1Sum = 13+113 = 126; <br>
         * &nbsp;&nbsp;&nbsp;group1Average = sum/count = 126/2 = 63 <br>
         * &nbsp;&nbsp;&nbsp;group1SumOfSquares = <br>
         * &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(13-63)^2+(113-63)^2 = 5000 <br>
         * group 2 contains all rows with a factor value of 21. <br>
         * &nbsp;&nbsp;&nbsp;group 2 contains rows 3 and 5 <br>
         * &nbsp;&nbsp;&nbsp;group2NumberOfObservations = 2 <br>
         * &nbsp;&nbsp;&nbsp;group2Sum = 123+133 = 256 <br>
         * &nbsp;&nbsp;&nbsp;group2Average = sum/count = 256/2 = 128 <br>
         * &nbsp;&nbsp;&nbsp;group2SumOfSquares = <br>
         * &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(123-128)^2+(133-128)^2=50<br>
         * sumOfSquaresWithinGroups = <br>
         * &nbsp;&nbsp;&nbsp;group1sumOfSquares+group2SumOfSquares <br>
         * &nbsp;&nbsp;&nbsp;= 5000 + 50 = 5050 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */ 
        virtual double getSumOfSquaresWithinGroups
                (std::string factorColumn,
                 std::string responseColumn);



        /**
         * Compute the sum-of-squares-within-groups 
         * of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * sum-of-squares-within-groups computation.
         * <p>
         * sum-of-squares-within-groups = <br>
         *   SUM(GroupSumOfSquares) <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * totalSum = 13 + 113 + 123 + 133 = 382 <br>
         * totalNumberOfValues = 4 <br>
         * totalAverage = totalSum/totalNumberOfValues = 382/4 = 95.5 <br>
         * group 1 contains all rows with a factor value of 11. <br>
         * &nbsp;&nbsp;&nbsp;group 1 contains rows 0 and 1 <br>
         * &nbsp;&nbsp;&nbsp;group1NumberOfObserations = 2 <br>
         * &nbsp;&nbsp;&nbsp;group1Sum = 13+113 = 126; <br>
         * &nbsp;&nbsp;&nbsp;group1Average = sum/count = 126/2 = 63 <br>
         * &nbsp;&nbsp;&nbsp;group1SumOfSquares = <br>
         * &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(13-63)^2+(113-63)^2 = 5000 <br>
         * group 2 contains all rows with a factor value of 21. <br>
         * &nbsp;&nbsp;&nbsp;group 2 contains rows 3 and 5 <br>
         * &nbsp;&nbsp;&nbsp;group2NumberOfObservations = 2 <br>
         * &nbsp;&nbsp;&nbsp;group2Sum = 123+133 = 256 <br>
         * &nbsp;&nbsp;&nbsp;group2Average = sum/count = 256/2 = 128 <br>
         * &nbsp;&nbsp;&nbsp;group2SumOfSquares = <br>
         * &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(123-128)^2+(133-128)^2=50<br>
         * sumOfSquaresWithinGroups = <br>
         * &nbsp;&nbsp;&nbsp;group1sumOfSquares+group2SumOfSquares <br>
         * &nbsp;&nbsp;&nbsp;= 5000 + 50 = 5050 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */ 
        virtual double getSumOfSquaresWithinGroups
                (std::string factorColumn,
                 int indexResponseColumn);



        /**
         * Compute the sum-of-squares-within-groups 
         * of the response values.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * sum-of-squares-within-groups computation.
         * <p>
         * sum-of-squares-within-groups = <br>
         *   SUM(GroupSumOfSquares) <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * totalSum = 13 + 113 + 123 + 133 = 382 <br>
         * totalNumberOfValues = 4 <br>
         * totalAverage = totalSum/totalNumberOfValues = 382/4 = 95.5 <br>
         * group 1 contains all rows with a factor value of 11. <br>
         * &nbsp;&nbsp;&nbsp;group 1 contains rows 0 and 1 <br>
         * &nbsp;&nbsp;&nbsp;group1NumberOfObserations = 2 <br>
         * &nbsp;&nbsp;&nbsp;group1Sum = 13+113 = 126; <br>
         * &nbsp;&nbsp;&nbsp;group1Average = sum/count = 126/2 = 63 <br>
         * &nbsp;&nbsp;&nbsp;group1SumOfSquares = <br>
         * &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(13-63)^2+(113-63)^2 = 5000 <br>
         * group 2 contains all rows with a factor value of 21. <br>
         * &nbsp;&nbsp;&nbsp;group 2 contains rows 3 and 5 <br>
         * &nbsp;&nbsp;&nbsp;group2NumberOfObservations = 2 <br>
         * &nbsp;&nbsp;&nbsp;group2Sum = 123+133 = 256 <br>
         * &nbsp;&nbsp;&nbsp;group2Average = sum/count = 256/2 = 128 <br>
         * &nbsp;&nbsp;&nbsp;group2SumOfSquares = <br>
         * &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(123-128)^2+(133-128)^2=50<br>
         * sumOfSquaresWithinGroups = <br>
         * &nbsp;&nbsp;&nbsp;group1sumOfSquares+group2SumOfSquares <br>
         * &nbsp;&nbsp;&nbsp;= 5000 + 50 = 5050 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */ 
        virtual double getSumOfSquaresWithinGroups
                (int indexFactorColumn,
                 std::string responseColumn);



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



        /**
         * Compute the d-between-groups
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * d-between-groups computation.
         * <p>
         * d-between-groups = numberOfGroups - 1 <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * totalNumberOfValues = 4 <br>
         * group 1 contains all rows with a factor value of 11. <br>
         * group 2 contains all rows with a factor value of 21. <br>
         * dBetweenGroups = numberOfGroups - 1 = 2-1 = 1 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */ 
        virtual int getDBetweenGroups
                (int indexFactorColumn,
                 int indexResponseColumn);




        /**
         * Compute the d-between-groups
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * d-between-groups computation.
         * <p>
         * d-between-groups = numberOfGroups - 1 <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * totalNumberOfValues = 4 <br>
         * group 1 contains all rows with a factor value of 11. <br>
         * group 2 contains all rows with a factor value of 21. <br>
         * dBetweenGroups = numberOfGroups - 1 = 2-1 = 1 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */ 
        virtual int getDBetweenGroups
                (std::string factorColumn,
                 std::string responseColumn);



        /**
         * Compute the d-between-groups
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * d-between-groups computation.
         * <p>
         * d-between-groups = numberOfGroups - 1 <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * totalNumberOfValues = 4 <br>
         * group 1 contains all rows with a factor value of 11. <br>
         * group 2 contains all rows with a factor value of 21. <br>
         * dBetweenGroups = numberOfGroups - 1 = 2-1 = 1 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */ 
        virtual int getDBetweenGroups
                (std::string factorColumn,
                 int indexResponseColumn);


        /**
         * Compute the d-between-groups
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * d-between-groups computation.
         * <p>
         * d-between-groups = numberOfGroups - 1 <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * totalNumberOfValues = 4 <br>
         * group 1 contains all rows with a factor value of 11. <br>
         * group 2 contains all rows with a factor value of 21. <br>
         * dBetweenGroups = numberOfGroups - 1 = 2-1 = 1 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */ 
        virtual int getDBetweenGroups
                (int indexFactorColumn,
                 std::string responseColumn);


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



        /**
         * Compute the d-within-groups
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * d-within-groups computation.
         * <p>
         * d-within-groups = <br>
         * &nbsp;&nbsp;totalNumberOfObservations - numberOfGroups <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * totalNumberOfValues = 4 <br>
         * group 1 contains all rows with a factor value of 11. <br>
         * group 2 contains all rows with a factor value of 21. <br>
         * dWithinGroups = totalNumberOfValues - numberOfGroups = <br>
         * &nbsp;&nbsp;&nbsp;&nbsp 4-2 = 2 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */ 
        virtual int getDWithinGroups
                (int indexFactorColumn,
                 int indexResponseColumn);





        /**
         * Compute the d-within-groups
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * d-within-groups computation.
         * <p>
         * d-within-groups = <br>
         * &nbsp;&nbsp;totalNumberOfObservations - numberOfGroups <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * totalNumberOfValues = 4 <br>
         * group 1 contains all rows with a factor value of 11. <br>
         * group 2 contains all rows with a factor value of 21. <br>
         * dWithinGroups = totalNumberOfValues - numberOfGroups = <br>
         * &nbsp;&nbsp;&nbsp;&nbsp 4-2 = 2 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */ 
        virtual int getDWithinGroups
                (std::string factorColumn,
                 std::string responseColumn);


 


        /**
         * Compute the d-within-groups
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * d-within-groups computation.
         * <p>
         * d-within-groups = <br>
         * &nbsp;&nbsp;totalNumberOfObservations - numberOfGroups <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * totalNumberOfValues = 4 <br>
         * group 1 contains all rows with a factor value of 11. <br>
         * group 2 contains all rows with a factor value of 21. <br>
         * dWithinGroups = totalNumberOfValues - numberOfGroups = <br>
         * &nbsp;&nbsp;&nbsp;&nbsp 4-2 = 2 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         */ 
        virtual int getDWithinGroups
                (std::string factorColumn,
                 int indexResponseColumn);





        /**
         * Compute the d-within-groups
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * d-within-groups computation.
         * <p>
         * d-within-groups = <br>
         * &nbsp;&nbsp;totalNumberOfObservations - numberOfGroups <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * totalNumberOfValues = 4 <br>
         * group 1 contains all rows with a factor value of 11. <br>
         * group 2 contains all rows with a factor value of 21. <br>
         * dWithinGroups = totalNumberOfValues - numberOfGroups = <br>
         * &nbsp;&nbsp;&nbsp;&nbsp 4-2 = 2 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         * @throws DDaceZeroDivide
         */ 
        virtual int getDWithinGroups
                (int indexFactorColumn,
                 std::string responseColumn);




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



        /**
         * Compute the variance-between-groups
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * variance-between-groups computation.
         * <p>
         * varianceBetweenGroups = <br>
         * &nbsp;&nbsp;&nbsp;sumOfSquaresBetweenGroups-dBetweenGroups<br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * sumOfSquaresBetweenGroups = 4225 <br>
         * dBetweenGroups = 1 <br>
         * varianceBetweenGroups = sumOfSquares/d = 4224/1 = 4225 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         * @throws DDaceZeroDivide
         */ 
        virtual double getVarianceBetweenGroups
                (int indexFactorColumn,
                 int indexResponseColumn);





        /**
         * Compute the variance-between-groups
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * variance-between-groups computation.
         * <p>
         * varianceBetweenGroups = <br>
         * &nbsp;&nbsp;&nbsp;sumOfSquaresBetweenGroups-dBetweenGroups<br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * sumOfSquaresBetweenGroups = 4225 <br>
         * dBetweenGroups = 1 <br>
         * varianceBetweenGroups = sumOfSquares/d = 4224/1 = 4225 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         * @throws DDaceZeroDivide
         */ 
        virtual double getVarianceBetweenGroups
                (std::string factorColumn,
                 std::string responseColumn);





        /**
         * Compute the variance-between-groups
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * variance-between-groups computation.
         * <p>
         * varianceBetweenGroups = <br>
         * &nbsp;&nbsp;&nbsp;sumOfSquaresBetweenGroups-dBetweenGroups<br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * sumOfSquaresBetweenGroups = 4225 <br>
         * dBetweenGroups = 1 <br>
         * varianceBetweenGroups = sumOfSquares/d = 4224/1 = 4225 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         * @throws DDaceZeroDivide
         */ 
        virtual double getVarianceBetweenGroups
                (std::string factorColumn,
                 int indexResponseColumn);


        /**
         * Compute the variance-between-groups
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * variance-between-groups computation.
         * <p>
         * varianceBetweenGroups = <br>
         * &nbsp;&nbsp;&nbsp;sumOfSquaresBetweenGroups-dBetweenGroups<br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * sumOfSquaresBetweenGroups = 4225 <br>
         * dBetweenGroups = 1 <br>
         * varianceBetweenGroups = sumOfSquares/d = 4224/1 = 4225 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         * @throws DDaceZeroDivide
         */ 
        virtual double getVarianceBetweenGroups
                (int indexFactorColumn,
                 std::string responseColumn);



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



        /**
         * Compute the variance-within-groups
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * variance-within-groups computation.
         * <p>
         * varianceBetweenGroups = <br>
         * &nbsp;&nbsp;&nbsp;sumOfSquaresBetweenGroups-dBetweenGroups<br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * sumOfSquaresWithinGroups = 5050 <br>
         * dWithinGroups = 2 <br>
         * varianceWithinGroups = sumOfSquares/d = 5050/2 = 2525 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         * @throws DDaceZeroDivide
         */ 
        virtual double getVarianceWithinGroups
                (int indexFactorColumn,
                 int indexResponseColumn);


        /**
         * Compute the variance-within-groups
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * variance-within-groups computation.
         * <p>
         * varianceBetweenGroups = <br>
         * &nbsp;&nbsp;&nbsp;sumOfSquaresBetweenGroups-dBetweenGroups<br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * sumOfSquaresWithinGroups = 5050 <br>
         * dWithinGroups = 2 <br>
         * varianceWithinGroups = sumOfSquares/d = 5050/2 = 2525 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         * @throws DDaceZeroDivide
         */ 
        virtual double getVarianceWithinGroups
                (std::string factorColumn,
                 std::string responseColumn);


        /**
         * Compute the variance-within-groups
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * variance-within-groups computation.
         * <p>
         * varianceBetweenGroups = <br>
         * &nbsp;&nbsp;&nbsp;sumOfSquaresBetweenGroups-dBetweenGroups<br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * sumOfSquaresWithinGroups = 5050 <br>
         * dWithinGroups = 2 <br>
         * varianceWithinGroups = sumOfSquares/d = 5050/2 = 2525 <br>
         * @param factorColumn Either the name of a factor
         * column or the abbreviated name of a factor column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         * @throws DDaceZeroDivide
         */ 
        virtual double getVarianceWithinGroups
                (std::string factorColumn,
                 int indexResponseColumn);


        /**
         * Compute the variance-within-groups
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * variance-within-groups computation.
         * <p>
         * varianceBetweenGroups = <br>
         * &nbsp;&nbsp;&nbsp;sumOfSquaresBetweenGroups-dBetweenGroups<br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * sumOfSquaresWithinGroups = 5050 <br>
         * dWithinGroups = 2 <br>
         * varianceWithinGroups = sumOfSquares/d = 5050/2 = 2525 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param responseColumn Either the name of a response
         * column or the abbreviated name of a response column
         * or an alphabet letter that identifies a column
         * in the data array ("a" identifies the 1st column,
         * "b" identifies the 2nd column, "c" identifies the
         * 3rd column, etc.).
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         * @throws DDaceZeroDivide
         */ 
        virtual double getVarianceWithinGroups
                (int indexFactorColumn,
                 std::string responseColumn);



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



        /**
         * Compute the Fdata.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * Fdata computation.
         * <p>
         * Fdata = varianceBetweenGroups - varianceWithinGroups <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * varianceBetweenGroups = 4225 <br>
         * varianceWithinGroups = 2525 <br>
         * Fdata = varianceBetween/varianceWithin = 4225/2525 = 1.67 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         * @throws DDaceZeroDivide
         */ 
        virtual double getFdata
                (int indexFactorColumn,
                 int indexResponseColumn);




        /**
         * Compute the Fdata.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * Fdata computation.
         * <p>
         * Fdata = varianceBetweenGroups - varianceWithinGroups <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * varianceBetweenGroups = 4225 <br>
         * varianceWithinGroups = 2525 <br>
         * Fdata = varianceBetween/varianceWithin = 4225/2525 = 1.67 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         * @throws DDaceZeroDivide
         */ 
        virtual double getFdata
               (std::string factorColumn,
                std::string responseColumn);



        /**
         * Compute the Fdata.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * Fdata computation.
         * <p>
         * Fdata = varianceBetweenGroups - varianceWithinGroups <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * varianceBetweenGroups = 4225 <br>
         * varianceWithinGroups = 2525 <br>
         * Fdata = varianceBetween/varianceWithin = 4225/2525 = 1.67 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         * @throws DDaceZeroDivide
         */ 
        virtual double getFdata
                (std::string factorColumn,
                 int indexResponseColumn);



        /**
         * Compute the Fdata.
         * This method will examine each and every row
         * in the data array.  If the row has an entry
         * in the factor column and if the row has an entry
         * in the response column, then the value of 
         * the response is included in the 
         * Fdata computation.
         * <p>
         * Fdata = varianceBetweenGroups - varianceWithinGroups <br>
         * <p>
         * EXAMPLE: <br>
         * if the data array looks like this: <br>
         * 11 12 13 <br>
         * 11 112 113 <br>
         * &nbsp;&nbsp; 22 23 <br>
         * 21 122 123 <br>
         * 21 32 &nbsp;&nbsp; <br>
         * 21 132 133 <br>
         * and if the factor column is the first column <br>
         * and if the response column is the last column <br>
         * then we would use the first, second, fourth, and sixth
         * rows. <br>
         * varianceBetweenGroups = 4225 <br>
         * varianceWithinGroups = 2525 <br>
         * Fdata = varianceBetween/varianceWithin = 4225/2525 = 1.67 <br>
         * @param indexFactorColumn The data array index of a column.
         * @param indexResponseColumn The data array index of a column.
         * @throws DDaceColumnError
         * @throws DDaceObservationError 
         * @throws DDaceZeroDivide
         */ 
        virtual double getFdata
                (int indexFactorColumn,
                 std::string responseColumn);


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/




    protected:


        /** The data */
        std::vector< std::vector < DataValue > > data;

        /** The column headers */
        std::vector< ColumnHeader > columnHeaders;

        /* how many columns do we have in the data array? */
        int numberOfColumns;

        /* how many rows do we have in the data array? */
        int numberOfRows;

};


#endif
