#include "MainEffectsAnalyzer3.h"
#include <cstring> // for strlen

using std::runtime_error;
    
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    MainEffectsAnalyzer3::MainEffectsAnalyzer3()
{
}

    MainEffectsAnalyzer3::MainEffectsAnalyzer3
            (std::vector< ColumnHeader > columnHeaders,
             std::vector< std::vector < DataValue > > data) {

        /* save the column headers & data */
        this->columnHeaders = columnHeaders;
        this->data = data;        

        this->numberOfRows = data.size();
        this->numberOfColumns = 0;
        if (this->numberOfRows>0)
            this->numberOfColumns = data[0].size();

    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    MainEffectsAnalyzer3::~MainEffectsAnalyzer3(){
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    int MainEffectsAnalyzer3::toIndexInputColumn(std::string columnId){

        /* does columnId match a column title?           */
        /* does columnId match an abbreviated col title? */
        for (int i=0; i<this->numberOfColumns; i++) {
            if (columnId == this->columnHeaders[i].getTitle())
                return(i);
            if (columnId == this->columnHeaders[i].getAbbreviatedTitle())
                return(i);
        }

        /* is columnId a single alphabet character? */
        char *ptr = (char *)(columnId.c_str());
        if (strlen(ptr)==1)
            if (isalpha(*ptr)) {
                if (islower(*ptr))
                   return(*ptr - 'a');
                if (isupper(*ptr))
                   return(*ptr - 'A');
            }

        throw runtime_error
            (columnId + " is not a column.");
     
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    std::vector< DataValue> MainEffectsAnalyzer3::getNonEmptyUniqueFactors
        (int indexFactorColumn) {


        /* convert the number of columns in the data array to a string */
        char sNumberOfColumns[128];
        sprintf(sNumberOfColumns,"%d",this->numberOfColumns);
            
        /* error check the input */
        if (indexFactorColumn<0)
            throw runtime_error
                ("Index value of column can not be negative.");
        if (indexFactorColumn >= this->numberOfColumns)
            throw runtime_error
                ("Index value must be smaler than "
                 + std::string(sNumberOfColumns));
        if (columnHeaders[indexFactorColumn].getFactorOrResponse()!=
            ColumnHeader::FACTOR)
            throw runtime_error
                (std::string("factor index must point ")
                +std::string("to a column containing factors"));

        /* store our unique factors here */
        std::vector< DataValue > uniqueFactors;
     
        /* process every row in the data array */
        for (int row=0; row<this->numberOfRows; row++) {

           /* retreive he factor value that is in the row */
           DataValue value = this->data[row][indexFactorColumn];
 
           /* if the retrieved DataValue is empty then do nothing */
           if (value.getDataType() == DataValue::EMPTY) {
               continue;
           }

           /* if the factor value is already in our set */
           /* then do nothing                           */
           if (isDataValueInSet(value, uniqueFactors)) {
                continue;
           }
           /* add the factor value to our set */
           uniqueFactors.push_back(value);
        }//for

        /* return the unique factors */
        return(uniqueFactors);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    bool MainEffectsAnalyzer3::isDataValueInSet
               (DataValue dataValue,
                std::vector< DataValue > set){

        std::vector< DataValue >::iterator iterator;
        for (iterator=set.begin(); iterator!=set.end(); iterator++) {
            if (dataValue.equals((DataValue)(*iterator)))
                return(true);
        }
        return(false);
    }





    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


        
    int MainEffectsAnalyzer3::getNumberOfObservations
        (int indexFactorColumn,
         int indexResponseColumn){

        

        /* convert the number of columns in the data array to a string */
        char sNumberOfColumns[128];
        sprintf(sNumberOfColumns,"%d",this->numberOfColumns);
            
        /* error check the inputs */
        if (indexFactorColumn<0)
            throw runtime_error
                ("Index value of column can not be negative.");
        if (indexResponseColumn<0)
            throw runtime_error
                ("Index value of column can not be negative.");
        if (indexFactorColumn >= this->numberOfColumns)
            throw runtime_error
                ("Index value must be smaler than "
                 + std::string(sNumberOfColumns));
        if (indexResponseColumn >= this->numberOfColumns)
            throw runtime_error
                ("Index value must be smaler than "
                 + std::string(sNumberOfColumns));
        if (this->numberOfRows==0)
            return(0);
        if (columnHeaders[indexFactorColumn].getFactorOrResponse()!=
            ColumnHeader::FACTOR)
            throw runtime_error
                (std::string("factor index must point ")
                +std::string("to a column containing factors"));
        if (columnHeaders[indexResponseColumn].getFactorOrResponse()!=
            ColumnHeader::RESPONSE)
            throw runtime_error
                (std::string("response index must point ")
                +std::string("to a column containing responses"));

        

        int numberOfValues = 0;

        for (int row=0; row<this->numberOfRows; row++) {

            /* get the factor value & and the response value */
            DataValue dataValueFactor = data[row][indexFactorColumn];
            DataValue dataValueResponse = data[row][indexResponseColumn];

            /* if either value is empty then do nothing */
            if (dataValueFactor.getDataType()==DataValue::EMPTY)
               continue;
            if (dataValueResponse.getDataType()==DataValue::EMPTY)
               continue;

            /* We have a factor value and a response value. */
            /* Increment our counter.                      */
            numberOfValues++;

        }//for

        /* return the number of values */
        return(numberOfValues);

    }//method


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    int MainEffectsAnalyzer3::getNumberOfObservations
         (std::string factorColumn,
          std::string responseColumn){
         
        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        int numberOfObservations = 
            this->getNumberOfObservations
               (indexFactorColumn,
                indexResponseColumn);
        return(numberOfObservations);
            
     }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/




    int MainEffectsAnalyzer3::getNumberOfObservations
        (std::string factorColumn,
         int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int numberOfObservations = 
            this->getNumberOfObservations
               (indexFactorColumn,
                indexResponseColumn);
        return(numberOfObservations);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getNumberOfObservations
        (int indexFactorColumn,
         std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        int numberOfObservations = 
            this->getNumberOfObservations
               (indexFactorColumn,
                indexResponseColumn);
        return(numberOfObservations);    
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/




    int MainEffectsAnalyzer3::getNumberOfObservations
        (int indexFactorColumn,
         DataValue factorValue,
         int indexResponseColumn){

        /* convert the number of columns in the data array to a string */
        char sNumberOfColumns[128];
        sprintf(sNumberOfColumns,"%d",this->numberOfColumns);
            
        /* error check the inputs */
        if (indexFactorColumn<0)
            throw runtime_error
                ("Index value of column can not be negative.");
        if (indexResponseColumn<0)
            throw runtime_error
                ("Index value of column can not be negative.");
        if (indexFactorColumn >= this->numberOfColumns)
            throw runtime_error
                ("Index value must be smaler than "
                 + std::string(sNumberOfColumns));
        if (indexResponseColumn >= this->numberOfColumns)
            throw runtime_error
                ("Index value must be smaler than "
                 + std::string(sNumberOfColumns));
        if (this->numberOfRows==0)
            return(0);
        if (columnHeaders[indexFactorColumn].getFactorOrResponse()!=
            ColumnHeader::FACTOR)
            throw runtime_error
                (std::string("factor index must point ")
                +std::string("to a column containing factors"));
        if (columnHeaders[indexResponseColumn].getFactorOrResponse()!=
            ColumnHeader::RESPONSE)
            throw runtime_error
                (std::string("response index must point ")
                +std::string("to a column containing responses"));

        

        int numberOfValues = 0;

        for (int row=0; row<this->numberOfRows; row++) {

            /* get the factor value & and the response value */
            DataValue dataValueFactor = data[row][indexFactorColumn];
            DataValue dataValueResponse = data[row][indexResponseColumn];

            /* if either value is empty then do nothing */
            if (dataValueFactor.getDataType()==DataValue::EMPTY)
               continue;
            if (dataValueResponse.getDataType()==DataValue::EMPTY)
               continue;



            /* if the dataValueFactor does NOT match the factorValue */
            /* then do nothing.                                      */
            if (dataValueFactor.getDataType()==DataValue::DOUBLE){
                if (dataValueFactor.getDoubleValue()
                    !=factorValue.getDoubleValue())
                    continue;
            }
            else if (dataValueFactor.getDataType()==DataValue::STRING){
                if (dataValueFactor.getStringValue()
                    !=factorValue.getStringValue())
                    continue;
            }
            else if (dataValueFactor.getDataType()==DataValue::INTEGER){
                if (dataValueFactor.getIntValue()
                    !=factorValue.getIntValue()) 
                    continue;                 
            }
            else {
               continue;
            }

            /* We have a factor value and a response value. */
            /* The input value matches the desired value.  */
            /* Increment our counter.                      */
            numberOfValues++;

        }//for

        /* return the number of values */
        return(numberOfValues);


    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/





    int MainEffectsAnalyzer3::getNumberOfObservations
         (std::string factorColumn,
          DataValue factorValue,
          std::string responseColumn){
         
        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        int numberOfObservations = 
            this->getNumberOfObservations
               (indexFactorColumn,
                factorValue,
                indexResponseColumn);
        return(numberOfObservations);
            
     }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/




    int MainEffectsAnalyzer3::getNumberOfObservations
        (std::string factorColumn,
         DataValue factorValue, 
         int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int numberOfObservations = 
            this->getNumberOfObservations
               (indexFactorColumn,
                factorValue,
                indexResponseColumn);
        return(numberOfObservations);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    int MainEffectsAnalyzer3::getNumberOfObservations
        (int indexFactorColumn,
         DataValue factorValue,
         std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        int numberOfObservations = 
            this->getNumberOfObservations
               (indexFactorColumn,
                factorValue,
                indexResponseColumn);
        return(numberOfObservations);    
    }





    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/





    int MainEffectsAnalyzer3::getNumberOfObservations
        (std::string factorColumn,
         std::string factorValue,
         std::string responseColumn){          

         int indexFactorColumn = toIndexInputColumn(factorColumn);
         int indexResponseColumn = toIndexInputColumn(responseColumn);
         DataValue dataValue(factorValue);
         int numberOfObservations =
             getNumberOfObservations
                (indexFactorColumn,
                 dataValue,
                 indexResponseColumn);
         return(numberOfObservations);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    int MainEffectsAnalyzer3::getNumberOfObservations
        (int indexFactorColumn,
         std::string factorValue,
         int indexResponseColumn){


         DataValue dataValue(factorValue);
         int numberOfObservations =
             getNumberOfObservations
                (indexFactorColumn,
                 dataValue,
                 indexResponseColumn);
         return(numberOfObservations);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    int MainEffectsAnalyzer3::getNumberOfObservations
        (std::string factorColumn,
         std::string factorValue,
         int indexResponseColumn){

         int indexFactorColumn = toIndexInputColumn(factorColumn);
         DataValue dataValue(factorValue);
         int numberOfObservations =
             getNumberOfObservations
                (indexFactorColumn,
                 dataValue,
                 indexResponseColumn);
         return(numberOfObservations);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    int MainEffectsAnalyzer3::getNumberOfObservations
        (int indexFactorColumn,
         std::string factorValue,
         std::string responseColumn){

         int indexResponseColumn = toIndexInputColumn(responseColumn);
         DataValue dataValue(factorValue);
         int numberOfObservations =
             getNumberOfObservations
                (indexFactorColumn,
                 dataValue,
                 indexResponseColumn);
         return(numberOfObservations);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    int MainEffectsAnalyzer3::getNumberOfObservations
        (std::string factorColumn,
         double factorValue,
         std::string responseColumn){

         int indexFactorColumn = toIndexInputColumn(factorColumn);
         int indexResponseColumn = toIndexInputColumn(responseColumn);
         DataValue dataValue(factorValue);
         int numberOfObservations =
             getNumberOfObservations
                (indexFactorColumn,
                 dataValue,
                 indexResponseColumn);
         return(numberOfObservations);

    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    int MainEffectsAnalyzer3::getNumberOfObservations
        (int indexFactorColumn,
         double factorValue,
         int indexResponseColumn){

         DataValue dataValue(factorValue);
         int numberOfObservations =
             getNumberOfObservations
                (indexFactorColumn,
                 dataValue,
                 indexResponseColumn);
         return(numberOfObservations);

    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    int MainEffectsAnalyzer3::getNumberOfObservations
        (std::string factorColumn,
         double factorValue,
         int indexResponseColumn){

         int indexFactorColumn = toIndexInputColumn(factorColumn);
         DataValue dataValue(factorValue);
         int numberOfObservations =
             getNumberOfObservations
                (indexFactorColumn,
                 dataValue,
                 indexResponseColumn);
         return(numberOfObservations);

    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    int MainEffectsAnalyzer3::getNumberOfObservations
        (int indexFactorColumn,
         double factorValue,
         std::string responseColumn){

         int indexResponseColumn = toIndexInputColumn(responseColumn);
         DataValue dataValue(factorValue);
         int numberOfObservations =
             getNumberOfObservations
                (indexFactorColumn,
                 dataValue,
                 indexResponseColumn);
         return(numberOfObservations);


    }






    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    int MainEffectsAnalyzer3::getNumberOfObservations
        (int indexFactorColumn,
         int factorValue,
         int indexResponseColumn){

         DataValue dataValue(factorValue);
         int numberOfObservations =
             getNumberOfObservations
                (indexFactorColumn,
                 dataValue,
                 indexResponseColumn);
         return(numberOfObservations);

    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    int MainEffectsAnalyzer3::getNumberOfObservations
        (std::string factorColumn,
         int factorValue,
         int indexResponseColumn){

         int indexFactorColumn = toIndexInputColumn(factorColumn);
         DataValue dataValue(factorValue);
         int numberOfObservations =
             getNumberOfObservations
                (indexFactorColumn,
                 dataValue,
                 indexResponseColumn);
         return(numberOfObservations);

    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    int MainEffectsAnalyzer3::getNumberOfObservations
        (int indexFactorColumn,
         int factorValue,
         std::string responseColumn){

         int indexResponseColumn = toIndexInputColumn(responseColumn);
         DataValue dataValue(factorValue);
         int numberOfObservations =
             getNumberOfObservations
                (indexFactorColumn,
                 dataValue,
                 indexResponseColumn);
         return(numberOfObservations);


    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    int MainEffectsAnalyzer3::getNumberOfObservations
        (std::string factorColumn,
         int factorValue,
         std::string responseColumn){

         int indexFactorColumn = toIndexInputColumn(factorColumn);
         int indexResponseColumn = toIndexInputColumn(responseColumn);
         DataValue dataValue(factorValue);
         int numberOfObservations =
             getNumberOfObservations
                (indexFactorColumn,
                 dataValue,
                 indexResponseColumn);
         return(numberOfObservations);


    }

    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfObservations
        (int indexFactorColumn,
         int indexResponseColumn){

        /* convert the number of columns in the data array to a string */
        char sNumberOfColumns[128];
        sprintf(sNumberOfColumns,"%d",this->numberOfColumns);
            
        /* error check the inputs */
        if (indexFactorColumn<0)
            throw runtime_error
                ("Index value of column can not be negative.");
        if (indexResponseColumn<0)
            throw runtime_error
                ("Index value of column can not be negative.");
        if (indexFactorColumn >= this->numberOfColumns)
            throw runtime_error
                ("Index value must be smaler than "
                 + std::string(sNumberOfColumns));
        if (indexResponseColumn >= this->numberOfColumns)
            throw runtime_error
                ("Index value must be smaler than "
                 + std::string(sNumberOfColumns));
        if (this->numberOfRows==0)
            return(0);
        if (columnHeaders[indexFactorColumn].getFactorOrResponse()!=
            ColumnHeader::FACTOR)
            throw runtime_error
                (std::string("factor index must point ")
                +std::string("to a column containing factors"));
        if (columnHeaders[indexResponseColumn].getFactorOrResponse()!=
            ColumnHeader::RESPONSE)
            throw runtime_error
                (std::string("response index must point ")
                +std::string("to a column containing responses"));

        
        double sumOfValues = 0.0;

        for (int row=0; row<this->numberOfRows; row++) {

            /* get the factor value & and the response value */
            DataValue dataValueFactor = data[row][indexFactorColumn];
            DataValue dataValueResponse = data[row][indexResponseColumn];

            /* if either value is empty then do nothing */
            if (dataValueFactor.getDataType()==DataValue::EMPTY)
               continue;
            if (dataValueResponse.getDataType()==DataValue::EMPTY)
               continue;

            /* We have a factor value and a response value. */
            /* add the value to our sum.                    */
            std::string dataTypeResponse = 
                 dataValueResponse.getDataType();
            if (dataTypeResponse==DataValue::DOUBLE){
                sumOfValues += dataValueResponse.getDoubleValue();
            }else if (dataTypeResponse==DataValue::INTEGER) {
                sumOfValues += (double)(dataValueResponse.getIntValue());
            }else if (dataTypeResponse==DataValue::STRING) {
                char *ptr = 
                 (char *)(dataValueResponse.getStringValue().c_str());
                sumOfValues += (double)(atof(ptr));
            } else {
               continue;
            }         
            

        }//for

        /* return the sum */
        return(sumOfValues);

    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfObservations
        (std::string factorColumn,
         std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double sumOfObservations = 
            this->getSumOfObservations
               (indexFactorColumn,
                indexResponseColumn);
        return(sumOfObservations);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfObservations
        (int indexFactorColumn,
         std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double sumOfObservations = 
            this->getSumOfObservations
               (indexFactorColumn,
                indexResponseColumn);
        return(sumOfObservations);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfObservations
        (std::string factorColumn,
         int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        double sumOfObservations = 
            this->getSumOfObservations
               (indexFactorColumn,
                indexResponseColumn);
        return(sumOfObservations);
    }

    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfObservations
        (int indexFactorColumn,
         DataValue factorValue,
         int indexResponseColumn){

        /* convert the number of columns in the data array to a string */
        char sNumberOfColumns[128];
        sprintf(sNumberOfColumns,"%d",this->numberOfColumns);
            
        /* error check the inputs */
        if (indexFactorColumn<0)
            throw runtime_error
                ("Index value of column can not be negative.");
        if (indexResponseColumn<0)
            throw runtime_error
                ("Index value of column can not be negative.");
        if (indexFactorColumn >= this->numberOfColumns)
            throw runtime_error
                ("Index value must be smaler than "
                 + std::string(sNumberOfColumns));
        if (indexResponseColumn >= this->numberOfColumns)
            throw runtime_error
                ("Index value must be smaler than "
                 + std::string(sNumberOfColumns));
        if (this->numberOfRows==0)
            return(0);
        if (columnHeaders[indexFactorColumn].getFactorOrResponse()!=
            ColumnHeader::FACTOR)
            throw runtime_error
                (std::string("factor index must point ")
                +std::string("to a column containing factors"));
        if (columnHeaders[indexResponseColumn].getFactorOrResponse()!=
            ColumnHeader::RESPONSE)
            throw runtime_error
                (std::string("response index must point ")
                +std::string("to a column containing responses"));

        

        double sumOfValues = 0.0;

        for (int row=0; row<this->numberOfRows; row++) {

            /* get the factor value & and the response value */
            DataValue dataValueFactor = data[row][indexFactorColumn];
            DataValue dataValueResponse = data[row][indexResponseColumn];

            /* if either value is empty then do nothing */
            if (dataValueFactor.getDataType()==DataValue::EMPTY)
               continue;
            if (dataValueResponse.getDataType()==DataValue::EMPTY)
               continue;



            /* if the dataValueFactor does NOT match the factorValue */
            /* then do nothing.                                      */
            if (dataValueFactor.getDataType()==DataValue::DOUBLE){
                if (dataValueFactor.getDoubleValue()
                    !=factorValue.getDoubleValue())
                    continue;
            }
            else if (dataValueFactor.getDataType()==DataValue::STRING){
                if (dataValueFactor.getStringValue()
                    !=factorValue.getStringValue())
                    continue;
            }
            else if (dataValueFactor.getDataType()==DataValue::INTEGER){
                if (dataValueFactor.getIntValue()
                    !=factorValue.getIntValue()) 
                    continue;                 
            }
            else {
               continue;
            }

            /* We have a factor value and a response value. */
            /* The input value matches the desired value.   */
            /* add the response value to our sum            */
            std::string dataTypeResponse = 
                 dataValueResponse.getDataType();
            if (dataTypeResponse==DataValue::DOUBLE){
                sumOfValues += dataValueResponse.getDoubleValue();
            }else if (dataTypeResponse==DataValue::INTEGER) {
                sumOfValues += (double)(dataValueResponse.getIntValue());
            }else if (dataTypeResponse==DataValue::STRING) {
                char *ptr = 
                 (char *)(dataValueResponse.getStringValue().c_str());
                sumOfValues += (double)(atof(ptr));
            } else {
               continue;
            }         
            


        }//for

        /* return the sum of the values */
        return(sumOfValues);
    }

    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfObservations
        (std::string factorColumn,
         DataValue factorValue,
         std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double sumOfObservations = 
            this->getSumOfObservations
               (indexFactorColumn,
                factorValue,
                indexResponseColumn);
        return(sumOfObservations);
    }

    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfObservations
        (int indexFactorColumn,
         DataValue factorValue,
         std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double sumOfObservations = 
            this->getSumOfObservations
               (indexFactorColumn,
                factorValue,
                indexResponseColumn);
        return(sumOfObservations);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfObservations
        (std::string factorColumn,
         DataValue factorValue,
         int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        double sumOfObservations = 
            this->getSumOfObservations
               (indexFactorColumn,
                factorValue,
                indexResponseColumn);
        return(sumOfObservations);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfObservations
        (int indexFactorColumn,
         double factorValue,
         int indexResponseColumn){

        DataValue dataValue(factorValue);
        double sumOfObservations = 
            this->getSumOfObservations
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(sumOfObservations);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfObservations
        (int indexFactorColumn,
         double factorValue,
         std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        DataValue dataValue(factorValue);
        double sumOfObservations = 
            this->getSumOfObservations
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(sumOfObservations);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfObservations
        (std::string factorColumn,
         double factorValue,
         int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        DataValue dataValue(factorValue);
        double sumOfObservations = 
            this->getSumOfObservations
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(sumOfObservations);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfObservations
        (std::string factorColumn,
         double factorValue,
         std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        DataValue dataValue(factorValue);
        double sumOfObservations = 
            this->getSumOfObservations
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(sumOfObservations);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfObservations
        (int indexFactorColumn,
         int factorValue,
         int indexResponseColumn){

        DataValue dataValue(factorValue);
        double sumOfObservations = 
            this->getSumOfObservations
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(sumOfObservations);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfObservations
        (int indexFactorColumn,
         int factorValue,
         std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        DataValue dataValue(factorValue);
        double sumOfObservations = 
            this->getSumOfObservations
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(sumOfObservations);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfObservations
        (std::string factorColumn,
         int factorValue,
         int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        DataValue dataValue(factorValue);
        double sumOfObservations = 
            this->getSumOfObservations
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(sumOfObservations);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
 

    double MainEffectsAnalyzer3::getSumOfObservations
        (std::string factorColumn,
         int factorValue,
         std::string responseColumn) {

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        DataValue dataValue(factorValue);
        double sumOfObservations = 
            this->getSumOfObservations
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(sumOfObservations);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfObservations
        (int indexFactorColumn,
         std::string factorValue,
         int indexResponseColumn){

        DataValue dataValue(factorValue);
        double sumOfObservations = 
            this->getSumOfObservations
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(sumOfObservations);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfObservations
        (int indexFactorColumn,
         std::string factorValue,
         std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        DataValue dataValue(factorValue);
        double sumOfObservations = 
            this->getSumOfObservations
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(sumOfObservations);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfObservations
        (std::string factorColumn,
         std::string factorValue,
         int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        DataValue dataValue(factorValue);
        double sumOfObservations = 
            this->getSumOfObservations
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(sumOfObservations);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfObservations
        (std::string factorColumn,
         std::string factorValue,
         std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        DataValue dataValue(factorValue);
        double sumOfObservations = 
            this->getSumOfObservations
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(sumOfObservations);
    }

    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (int indexFactorColumn,
             int indexResponseColumn){

        int numberOfObservations = 
            getNumberOfObservations
                (indexFactorColumn,
                 indexResponseColumn);
        double sumOfObservations =
            getSumOfObservations
                (indexFactorColumn,
                 indexResponseColumn);
       if (numberOfObservations==0)
           throw runtime_error
              ("Need at least one observation to compute an average");
       double average = sumOfObservations/(double)numberOfObservations;
       return(average);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (std::string factorColumn,
             std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double averageObservation = 
            this->getAverageObservation
               (indexFactorColumn,
                indexResponseColumn);
        return(averageObservation);

    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (std::string factorColumn,
             int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        double averageObservation = 
            this->getAverageObservation
               (indexFactorColumn,
                indexResponseColumn);
        return(averageObservation);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (int indexFactorColumn,
             std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double averageObservation = 
            this->getAverageObservation
               (indexFactorColumn,
                indexResponseColumn);
        return(averageObservation);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (int indexFactorColumn,
             DataValue factorValue,
             int indexResponseColumn){

        int numberOfObservations = 
            getNumberOfObservations
                (indexFactorColumn,
                 factorValue,
                 indexResponseColumn);
        double sumOfObservations =
            getSumOfObservations
                (indexFactorColumn,
                 factorValue,
                 indexResponseColumn);
       if (numberOfObservations==0)
           throw runtime_error
              ("Need at least one observation to compute an average");
       double average = sumOfObservations/(double)numberOfObservations;
       return(average);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (std::string factorColumn,
             DataValue factorValue,
             std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double averageObservation = 
            this->getAverageObservation
               (indexFactorColumn,
                factorValue,
                indexResponseColumn);
        return(averageObservation);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (std::string factorColumn,
             DataValue factorValue,
             int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        double averageObservation = 
            this->getAverageObservation
               (indexFactorColumn,
                factorValue,
                indexResponseColumn);
        return(averageObservation);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (int indexFactorColumn,
             DataValue factorValue,
             std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double averageObservation = 
            this->getAverageObservation
               (indexFactorColumn,
                factorValue,
                indexResponseColumn);
        return(averageObservation);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (int indexFactorColumn,
             double factorValue,
             int indexResponseColumn){

        DataValue dataValue(factorValue);
        double averageObservation = 
            this->getAverageObservation
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(averageObservation);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (std::string factorColumn,
             double factorValue,
             std::string responseColumn){

        DataValue dataValue(factorValue);
        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double averageObservation = 
            this->getAverageObservation
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(averageObservation);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (std::string factorColumn,
             double factorValue,
             int indexResponseColumn){

        DataValue dataValue(factorValue);
        int indexFactorColumn = toIndexInputColumn(factorColumn);
        double averageObservation = 
            this->getAverageObservation
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(averageObservation);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (int indexFactorColumn,
             double factorValue,
             std::string responseColumn){

        DataValue dataValue(factorValue);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double averageObservation = 
            this->getAverageObservation
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(averageObservation);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (int indexFactorColumn,
             int factorValue,
             int indexResponseColumn){

        DataValue dataValue(factorValue);
        double averageObservation = 
            this->getAverageObservation
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(averageObservation);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (std::string factorColumn,
             int factorValue,
             std::string responseColumn){

        DataValue dataValue(factorValue);
        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double averageObservation = 
            this->getAverageObservation
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(averageObservation);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (std::string factorColumn,
             int factorValue,
             int indexResponseColumn){

        DataValue dataValue(factorValue);
        int indexFactorColumn = toIndexInputColumn(factorColumn);
        double averageObservation = 
            this->getAverageObservation
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(averageObservation);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (int indexFactorColumn,
             int factorValue,
             std::string responseColumn){

        DataValue dataValue(factorValue);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double averageObservation = 
            this->getAverageObservation
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(averageObservation);
    }





    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (int indexFactorColumn,
             std::string factorValue,
             int indexResponseColumn){

        DataValue dataValue(factorValue);
        double averageObservation = 
            this->getAverageObservation
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(averageObservation);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (std::string factorColumn,
             std::string factorValue,
             std::string responseColumn){

        DataValue dataValue(factorValue);
        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double averageObservation = 
            this->getAverageObservation
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(averageObservation);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (std::string factorColumn,
             std::string factorValue,
             int indexResponseColumn){

        DataValue dataValue(factorValue);
        int indexFactorColumn = toIndexInputColumn(factorColumn);
        double averageObservation = 
            this->getAverageObservation
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(averageObservation);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getAverageObservation
            (int indexFactorColumn,
             std::string factorValue,
             std::string responseColumn){

        DataValue dataValue(factorValue);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double averageObservation = 
            this->getAverageObservation
               (indexFactorColumn,
                dataValue,
                indexResponseColumn);
        return(averageObservation);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfSquares
            (int indexFactorColumn,
             int indexResponseColumn){

        /* convert the number of columns in the data array to a string */
        char sNumberOfColumns[128];
        sprintf(sNumberOfColumns,"%d",this->numberOfColumns);
            
        /* error check the inputs */
        if (indexFactorColumn<0)
            throw runtime_error
                ("Index value of column can not be negative.");
        if (indexResponseColumn<0)
            throw runtime_error
                ("Index value of column can not be negative.");
        if (indexFactorColumn >= this->numberOfColumns)
            throw runtime_error
                ("Index value must be smaler than "
                 + std::string(sNumberOfColumns));
        if (indexResponseColumn >= this->numberOfColumns)
            throw runtime_error
                ("Index value must be smaler than "
                 + std::string(sNumberOfColumns));
        if (this->numberOfRows==0)
            return(0);
        if (columnHeaders[indexFactorColumn].getFactorOrResponse()!=
            ColumnHeader::FACTOR)
            throw runtime_error
                (std::string("factor index must point ")
                +std::string("to a column containing factors"));
        if (columnHeaders[indexResponseColumn].getFactorOrResponse()!=
            ColumnHeader::RESPONSE)
            throw runtime_error
                (std::string("response index must point ")
                +std::string("to a column containing responses"));

        
        /* get the average value */
        double averageValue = getAverageObservation
               (indexFactorColumn,
                indexResponseColumn);
        double sumOfSquares = 0.0;

        for (int row=0; row<this->numberOfRows; row++) {

            /* get the factor value & and the response value */
            DataValue dataValueFactor = data[row][indexFactorColumn];
            DataValue dataValueResponse = data[row][indexResponseColumn];

            /* if either value is empty then do nothing */
            if (dataValueFactor.getDataType()==DataValue::EMPTY)
               continue;
            if (dataValueResponse.getDataType()==DataValue::EMPTY)
               continue;

            /* We have a factor value and a response value. */
            /* get the response value                       */
            double value = averageValue;
            std::string dataTypeResponse = 
                 dataValueResponse.getDataType();
            if (dataTypeResponse==DataValue::DOUBLE){
                value = dataValueResponse.getDoubleValue();
            }else if (dataTypeResponse==DataValue::INTEGER) {
                value = (double)(dataValueResponse.getIntValue());
            }else if (dataTypeResponse==DataValue::STRING) {
                char *ptr = 
                 (char *)(dataValueResponse.getStringValue().c_str());
                value = (double)(atof(ptr));
            } else {
               continue;
            }         
            
            /* add (value-average)^2 to sumOfSquares */
            double difference = (value-averageValue);
            sumOfSquares += difference*difference;

        }//for

        /* return the sumOfSquares */
        return(sumOfSquares);

        
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    double MainEffectsAnalyzer3::getSumOfSquares
            (std::string factorColumn,
             std::string responseColumn){

         
        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double sumOfSquares = 
            this->getSumOfSquares
               (indexFactorColumn,
                indexResponseColumn);
        return(sumOfSquares);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfSquares
            (std::string factorColumn,
             int indexResponseColumn){

         
        int indexFactorColumn = toIndexInputColumn(factorColumn);
        double sumOfSquares = 
            this->getSumOfSquares
               (indexFactorColumn,
                indexResponseColumn);
        return(sumOfSquares);
    }


    double MainEffectsAnalyzer3::getSumOfSquares
            (int indexFactorColumn,
             std::string responseColumn){
         
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double sumOfSquares = 
            this->getSumOfSquares
               (indexFactorColumn,
                indexResponseColumn);
        return(sumOfSquares);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfSquares
            (int indexFactorColumn,
             DataValue factorValue,
             int indexResponseColumn){

        
        /* convert the number of columns in the data array to a string */
        char sNumberOfColumns[128];
        sprintf(sNumberOfColumns,"%d",this->numberOfColumns);
            
        /* error check the inputs */
        if (indexFactorColumn<0)
            throw runtime_error
                ("Index value of column can not be negative.");
        if (indexResponseColumn<0)
            throw runtime_error
                ("Index value of column can not be negative.");
        if (indexFactorColumn >= this->numberOfColumns)
            throw runtime_error
                ("Index value must be smaler than "
                 + std::string(sNumberOfColumns));
        if (indexResponseColumn >= this->numberOfColumns)
            throw runtime_error
                ("Index value must be smaler than "
                 + std::string(sNumberOfColumns));
        if (this->numberOfRows==0)
            return(0);
        if (columnHeaders[indexFactorColumn].getFactorOrResponse()!=
            ColumnHeader::FACTOR)
            throw runtime_error
                (std::string("factor index must point ")
                +std::string("to a column containing factors"));
        if (columnHeaders[indexResponseColumn].getFactorOrResponse()!=
            ColumnHeader::RESPONSE)
            throw runtime_error
                (std::string("response index must point ")
                +std::string("to a column containing responses"));

        
        /* get the average value */
        double averageValue = getAverageObservation
               (indexFactorColumn,
                factorValue,
                indexResponseColumn);

        double sumOfSquares = 0.0;

        for (int row=0; row<this->numberOfRows; row++) {

            /* get the factor value & and the response value */
            DataValue dataValueFactor = data[row][indexFactorColumn];
            DataValue dataValueResponse = data[row][indexResponseColumn];

            /* if either value is empty then do nothing */
            if (dataValueFactor.getDataType()==DataValue::EMPTY)
               continue;
            if (dataValueResponse.getDataType()==DataValue::EMPTY)
               continue;



            /* if the dataValueFactor does NOT match the factorValue */
            /* then do nothing.                                      */
            if (dataValueFactor.getDataType()==DataValue::DOUBLE){
                if (dataValueFactor.getDoubleValue()
                    !=factorValue.getDoubleValue())
                    continue;
            }
            else if (dataValueFactor.getDataType()==DataValue::STRING){
                if (dataValueFactor.getStringValue()
                    !=factorValue.getStringValue())
                    continue;
            }
            else if (dataValueFactor.getDataType()==DataValue::INTEGER){
                if (dataValueFactor.getIntValue()
                    !=factorValue.getIntValue()) 
                    continue;                 
            }
            else {
               continue;
            }


            /* We have a factor value and a response value.    */
            /* The input value matches the desired value.      */
            /* get the response value                          */
            double value = 0.0;
            std::string dataTypeResponse = 
                 dataValueResponse.getDataType();
            if (dataTypeResponse==DataValue::DOUBLE){
                value = dataValueResponse.getDoubleValue();
            }else if (dataTypeResponse==DataValue::INTEGER) {
                value = (double)(dataValueResponse.getIntValue());
            }else if (dataTypeResponse==DataValue::STRING) {
                char *ptr = 
                 (char *)(dataValueResponse.getStringValue().c_str());
                value = (double)(atof(ptr));
            } else {
               continue;
            }         
            

            double difference = value - averageValue;
            sumOfSquares += difference * difference;


        }//for

        /* return the sum of the squares */
        return(sumOfSquares);


    }

    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    double MainEffectsAnalyzer3::getSumOfSquares
            (std::string factorColumn,
             DataValue factorValue,
             std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double sumOfSquares = 
            this->getSumOfSquares
               (indexFactorColumn,
                factorValue,
                indexResponseColumn);
        return(sumOfSquares);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    double MainEffectsAnalyzer3::getSumOfSquares
            (std::string factorColumn,
             DataValue factorValue,
             int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        double sumOfSquares = 
            this->getSumOfSquares
               (indexFactorColumn,
                factorValue,
                indexResponseColumn);
        return(sumOfSquares);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/



    double MainEffectsAnalyzer3::getSumOfSquares
            (int indexFactorColumn,
             DataValue factorValue,
             std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double sumOfSquares = 
            this->getSumOfSquares
               (indexFactorColumn,
                factorValue,
                indexResponseColumn);
        return(sumOfSquares);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfSquares
            (int indexFactorColumn,
             double factorValue,
             int indexResponseColumn){

        DataValue dataValue = DataValue(factorValue);
        double sumOfSquares = getSumOfSquares
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(sumOfSquares);        
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfSquares
            (std::string factorColumn,
             double factorValue,
             std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        DataValue dataValue = DataValue(factorValue);
        double sumOfSquares = getSumOfSquares
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(sumOfSquares);        
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfSquares
            (std::string factorColumn,
             double factorValue,
             int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        DataValue dataValue = DataValue(factorValue);
        double sumOfSquares = getSumOfSquares
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(sumOfSquares);        
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfSquares
            (int indexFactorColumn,
             double factorValue,
             std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        DataValue dataValue = DataValue(factorValue);
        double sumOfSquares = getSumOfSquares
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(sumOfSquares);        
    }






    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfSquares
            (int indexFactorColumn,
             int factorValue,
             int indexResponseColumn){

        DataValue dataValue = DataValue(factorValue);
        double sumOfSquares = getSumOfSquares
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(sumOfSquares);        
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfSquares
            (std::string factorColumn,
             int factorValue,
             std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        DataValue dataValue = DataValue(factorValue);
        double sumOfSquares = getSumOfSquares
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(sumOfSquares);        
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfSquares
            (std::string factorColumn,
             int factorValue,
             int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        DataValue dataValue = DataValue(factorValue);
        double sumOfSquares = getSumOfSquares
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(sumOfSquares);        
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfSquares
            (int indexFactorColumn,
             int factorValue,
             std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        DataValue dataValue = DataValue(factorValue);
        double sumOfSquares = getSumOfSquares
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(sumOfSquares);        
    }





    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfSquares
            (int indexFactorColumn,
             std::string factorValue,
             int indexResponseColumn){

        DataValue dataValue = DataValue(factorValue);
        double sumOfSquares = getSumOfSquares
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(sumOfSquares);        
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfSquares
            (std::string factorColumn,
             std::string factorValue,
             std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        DataValue dataValue = DataValue(factorValue);
        double sumOfSquares = getSumOfSquares
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(sumOfSquares);        
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfSquares
            (std::string factorColumn,
             std::string factorValue,
             int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        DataValue dataValue = DataValue(factorValue);
        double sumOfSquares = getSumOfSquares
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(sumOfSquares);        
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfSquares
            (int indexFactorColumn,
             std::string factorValue,
             std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        DataValue dataValue = DataValue(factorValue);
        double sumOfSquares = getSumOfSquares
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(sumOfSquares);        
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (int indexFactorColumn,
             int indexResponseColumn){

        
        /* get the average value */
        double sumOfSquares = getSumOfSquares
               (indexFactorColumn,
                indexResponseColumn);
        
       /* get the number of values */
       int numberOfValues = getNumberOfObservations
               (indexFactorColumn,
                indexResponseColumn);

       if (numberOfValues<=1)
           throw runtime_error
                 ("Need 2 or more observations to compute variance");

       /* compute the variance */
       double variance = sumOfSquares / (double)(numberOfValues-1);


        /* return the variance */
        return(variance);
        
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (std::string factorColumn,
             std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double variance = getVariance
             (indexFactorColumn,
              indexResponseColumn);
        return(variance);
    }

    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (int indexFactorColumn,
             std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double variance = getVariance
             (indexFactorColumn,
              indexResponseColumn);
        return(variance);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (std::string factorColumn,
             int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        double variance = getVariance
             (indexFactorColumn,
              indexResponseColumn);
        return(variance);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (int indexFactorColumn,
             DataValue factorValue,
             int indexResponseColumn){

        /* get the average value */
        double sumOfSquares = getSumOfSquares
               (indexFactorColumn,
                factorValue,
                indexResponseColumn);
        
       /* get the number of values */
       int numberOfValues = getNumberOfObservations
               (indexFactorColumn,
                factorValue,
                indexResponseColumn);

       if (numberOfValues<=1)
           throw runtime_error
                 ("Need 2 or more observations to compute variance");

       /* compute the variance */
       double variance = sumOfSquares / (double)(numberOfValues-1);


        /* return the variance */
        return(variance);
    
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (std::string factorColumn,
             DataValue factorValue,
             std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double variance = getVariance
             (indexFactorColumn,
              factorValue,
              indexResponseColumn);
        return(variance);

    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (int indexFactorColumn,
             DataValue factorValue,
             std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double variance = getVariance
             (indexFactorColumn,
              factorValue,
              indexResponseColumn);
        return(variance);

    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (std::string factorColumn,
             DataValue factorValue,
             int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        double variance = getVariance
             (indexFactorColumn,
              factorValue,
              indexResponseColumn);
        return(variance);

    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (int indexFactorColumn,
             double factorValue,
             int indexResponseColumn){

        DataValue dataValue(factorValue);
        double variance = getVariance
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(variance);

    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (std::string factorColumn,
             double factorValue,
             std::string responseColumn){

        DataValue dataValue(factorValue);
        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double variance = getVariance
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(variance);

    }






    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (int indexFactorColumn,
             double factorValue,
             std::string responseColumn){

        DataValue dataValue(factorValue);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double variance = getVariance
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(variance);

    }





    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (std::string factorColumn,
             double factorValue,
             int indexResponseColumn){

        DataValue dataValue(factorValue);
        int indexFactorColumn = toIndexInputColumn(factorColumn);
        double variance = getVariance
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(variance);

    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (int indexFactorColumn,
             int factorValue,
             int indexResponseColumn){

        DataValue dataValue(factorValue);
        double variance = getVariance
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(variance);
    
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (std::string factorColumn,
             int factorValue,
             std::string responseColumn){

        DataValue dataValue(factorValue);
        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double variance = getVariance
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(variance);

    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (int indexFactorColumn,
             int factorValue,
             std::string responseColumn){

        DataValue dataValue(factorValue);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double variance = getVariance
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(variance);

    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (std::string factorColumn,
             int factorValue,
             int indexResponseColumn){

        DataValue dataValue(factorValue);
        int indexFactorColumn = toIndexInputColumn(factorColumn);
        double variance = getVariance
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(variance);

    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (int indexFactorColumn,
             std::string factorValue,
             int indexResponseColumn){

        DataValue dataValue(factorValue);
        double variance = getVariance
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(variance);
    
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (std::string factorColumn,
             std::string factorValue,
             std::string responseColumn){

        DataValue dataValue(factorValue);
        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double variance = getVariance
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(variance);

    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (int indexFactorColumn,
             std::string factorValue,
             std::string responseColumn){

        DataValue dataValue(factorValue);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double variance = getVariance
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(variance);

    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getVariance
            (std::string factorColumn,
             std::string factorValue,
             int indexResponseColumn){

        DataValue dataValue(factorValue);
        int indexFactorColumn = toIndexInputColumn(factorColumn);
        double variance = getVariance
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(variance);

    }

    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (int indexFactorColumn,
             int indexResponseColumn){

        int numberOfObservations = getNumberOfObservations
            (indexFactorColumn,
             indexResponseColumn);
        if (numberOfObservations==0)
           throw runtime_error
             ("You must have at least one observation to compute d.");
        int d = numberOfObservations - 1;
        return(d);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (std::string factorColumn,
             std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        int d = getD
             (indexFactorColumn,
              indexResponseColumn);
        return(d);

    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (int indexFactorColumn,
             std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        int d = getD
             (indexFactorColumn,
              indexResponseColumn);
        return(d);

    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (std::string factorColumn,
             int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int d = getD
             (indexFactorColumn,
              indexResponseColumn);
        return(d);

    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (int indexFactorColumn,
             DataValue factorValue,
             int indexResponseColumn){

        int numberOfObservations = getNumberOfObservations
            (indexFactorColumn,
             factorValue,
             indexResponseColumn);
        if (numberOfObservations==0)
           throw runtime_error
             ("You must have at least one observation to compute d.");
        int d = numberOfObservations - 1;
        return(d);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (std::string factorColumn,
             DataValue factorValue,
             std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        int d = getD
             (indexFactorColumn,
              indexResponseColumn);
        return(d);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (int indexFactorColumn,
             DataValue factorValue,
             std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        int d = getD
             (indexFactorColumn,
              indexResponseColumn);
        return(d);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (std::string factorColumn,
             DataValue factorValue,
             int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int d = getD
             (indexFactorColumn,
              indexResponseColumn);
        return(d);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (int indexFactorColumn,
             double factorValue,
             int indexResponseColumn){

        DataValue dataValue(factorValue);
        int d = getD
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(d);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (std::string factorColumn,
             double factorValue,
             std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        DataValue dataValue(factorValue);
        int d = getD
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(d);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (int indexFactorColumn,
             double factorValue,
             std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        DataValue dataValue(factorValue);
        int d = getD
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(d);
    }





    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (std::string factorColumn,
             double factorValue,
             int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        DataValue dataValue(factorValue);
        int d = getD
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(d);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (int indexFactorColumn,
             int factorValue,
             int indexResponseColumn){

        DataValue dataValue(factorValue);
        int d = getD
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(d);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (std::string factorColumn,
             int factorValue,
             std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        DataValue dataValue(factorValue);
        int d = getD
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(d);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (int indexFactorColumn,
             int factorValue,
             std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        DataValue dataValue(factorValue);
        int d = getD
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(d);
    }





    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (std::string factorColumn,
             int factorValue,
             int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        DataValue dataValue(factorValue);
        int d = getD
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(d);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (int indexFactorColumn,
             std::string factorValue,
             int indexResponseColumn){

        DataValue dataValue(factorValue);
        int d = getD
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(d);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (std::string factorColumn,
             std::string factorValue,
             std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        DataValue dataValue(factorValue);
        int d = getD
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(d);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (int indexFactorColumn,
             std::string factorValue,
             std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        DataValue dataValue(factorValue);
        int d = getD
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(d);
    }





    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getD
            (std::string factorColumn,
             std::string factorValue,
             int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        DataValue dataValue(factorValue);
        int d = getD
             (indexFactorColumn,
              dataValue,
              indexResponseColumn);
        return(d);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/




    double MainEffectsAnalyzer3::getSumOfSquaresBetweenGroups
                (int indexFactorColumn,
                 int indexResponseColumn){

        double sumOfSquaresBetweenGroups = 0.0;


        /* Get the total average */
        double totalAverage = getAverageObservation 
                (indexFactorColumn,
                 indexResponseColumn);


        /* get unique factors */
        std::vector< DataValue > uniqueFactors = 
            getNonEmptyUniqueFactors(indexFactorColumn);
 
       
        /* process each and every unique factor */
        std::vector< DataValue >::iterator iterator;
        for (iterator = uniqueFactors.begin();
             iterator != uniqueFactors.end();
             iterator++) {

            /* get one unique factor */
            DataValue uniqueFactor = (DataValue)(*iterator);

            /* a group contains rows which have this unique factor */
            /* collect some stats on this group                    */
            int groupNumberOfObservations = getNumberOfObservations
                (indexFactorColumn,
                 uniqueFactor,
                 indexResponseColumn);
                
            double groupAverage = getAverageObservation 
                (indexFactorColumn,
                 uniqueFactor,
                 indexResponseColumn);

            double groupAvgMinusTotalAvg =
                groupAverage - totalAverage;
            double groupSumOfSquares =
                (double)(groupNumberOfObservations)
              * groupAvgMinusTotalAvg * groupAvgMinusTotalAvg;
            
            sumOfSquaresBetweenGroups += groupSumOfSquares;
        }

        return(sumOfSquaresBetweenGroups);

    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/




    double MainEffectsAnalyzer3::getSumOfSquaresBetweenGroups
                (std::string factorColumn,
                 std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double sumOfSquaresBetweenGroups = 
            getSumOfSquaresBetweenGroups
               (indexFactorColumn,
                indexResponseColumn);
        return(sumOfSquaresBetweenGroups);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/




    double MainEffectsAnalyzer3::getSumOfSquaresBetweenGroups
                (std::string factorColumn,
                 int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        double sumOfSquaresBetweenGroups = 
            getSumOfSquaresBetweenGroups
               (indexFactorColumn,
                indexResponseColumn);
        return(sumOfSquaresBetweenGroups);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/




    double MainEffectsAnalyzer3::getSumOfSquaresBetweenGroups
                (int indexFactorColumn,
                 std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double sumOfSquaresBetweenGroups = 
            getSumOfSquaresBetweenGroups
               (indexFactorColumn,
                indexResponseColumn);
        return(sumOfSquaresBetweenGroups);
    }






    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfSquaresWithinGroups
                (int indexFactorColumn,
                 int indexResponseColumn){

        
        double sumOfSquaresWithinGroups = 0.0;


        /* get unique factors */
        std::vector< DataValue > uniqueFactors = 
            getNonEmptyUniqueFactors(indexFactorColumn);
 
       
        /* process each and every unique factor */
        std::vector< DataValue >::iterator iterator;
        for (iterator = uniqueFactors.begin();
             iterator != uniqueFactors.end();
             iterator++) {

            /* get one unique factor */
            DataValue uniqueFactor = (DataValue)(*iterator);

            /* a group contains rows which have this unique factor */
            /* collect some stats on this group                    */
            double sumOfSquares = getSumOfSquares
                (indexFactorColumn,
                 uniqueFactor,
                 indexResponseColumn);           
                            
            sumOfSquaresWithinGroups += sumOfSquares;
        }

        return(sumOfSquaresWithinGroups);

    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfSquaresWithinGroups
                (std::string factorColumn,
                 std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double sumOfSquaresWithinGroups = 
            getSumOfSquaresWithinGroups
               (indexFactorColumn,
                indexResponseColumn);
        return(sumOfSquaresWithinGroups);

    } 



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfSquaresWithinGroups
                (std::string factorColumn,
                 int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        double sumOfSquaresWithinGroups = 
            getSumOfSquaresWithinGroups
               (indexFactorColumn,
                indexResponseColumn);
        return(sumOfSquaresWithinGroups);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getSumOfSquaresWithinGroups
                (int indexFactorColumn,
                 std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double sumOfSquaresWithinGroups = 
            getSumOfSquaresWithinGroups
               (indexFactorColumn,
                indexResponseColumn);
        return(sumOfSquaresWithinGroups);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    int MainEffectsAnalyzer3::getDBetweenGroups
                (int indexFactorColumn,
                 int indexResponseColumn){

        /* get unique factors */
        std::vector< DataValue > uniqueFactors = 
            getNonEmptyUniqueFactors(indexFactorColumn);

        /* how many groups are there? */
        int numberOfGroups = uniqueFactors.size();

        /* compute d-between-groups */
        int d = numberOfGroups - 1;

        /* return d-between-groups */
        return(d);

    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    int MainEffectsAnalyzer3::getDBetweenGroups
                (std::string factorColumn,
                 std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        int dBetweenGroups = 
            getDBetweenGroups
               (indexFactorColumn,
                indexResponseColumn);
        return(dBetweenGroups);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    int MainEffectsAnalyzer3::getDBetweenGroups
                (std::string factorColumn,
                 int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int dBetweenGroups = 
            getDBetweenGroups
               (indexFactorColumn,
                indexResponseColumn);
        return(dBetweenGroups);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    int MainEffectsAnalyzer3::getDBetweenGroups
                (int indexFactorColumn,
                 std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        int dBetweenGroups = 
            getDBetweenGroups
               (indexFactorColumn,
                indexResponseColumn);
        return(dBetweenGroups);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getDWithinGroups
                (int indexFactorColumn,
                 int indexResponseColumn){

        /* get total number of observations */
        int totalNumberOfObservations = getNumberOfObservations
            (indexFactorColumn,
             indexResponseColumn);

        /* get unique factors */
        std::vector< DataValue > uniqueFactors = 
            getNonEmptyUniqueFactors(indexFactorColumn);

        /* how many groups are there? */
        int numberOfGroups = uniqueFactors.size();

        /* compute d-within-groups */
        int d = totalNumberOfObservations - numberOfGroups;

        /* return d-within-groups */
        return(d);        

    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getDWithinGroups
                (std::string factorColumn,
                 std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        int dWithinGroups = 
            getDWithinGroups
               (indexFactorColumn,
                indexResponseColumn);
        return(dWithinGroups);

    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getDWithinGroups
                (std::string factorColumn,
                 int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int dWithinGroups = 
            getDWithinGroups
               (indexFactorColumn,
                indexResponseColumn);
        return(dWithinGroups);

    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    int MainEffectsAnalyzer3::getDWithinGroups
                (int indexFactorColumn,
                 std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        int dWithinGroups = 
            getDWithinGroups
               (indexFactorColumn,
                indexResponseColumn);
        return(dWithinGroups);

    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getVarianceBetweenGroups
                (int indexFactorColumn,
                 int indexResponseColumn){

        /* get the sumOfSquaresBetweenGroups */
        double sumOfSquaresBetweenGroups =
               getSumOfSquaresBetweenGroups
                    (indexFactorColumn,
                     indexResponseColumn);

        /* get the dBetweenGroups */
        double dBetweenGroups =
               getDBetweenGroups
                    (indexFactorColumn,
                     indexResponseColumn);

        /* do some error checking */
        if (dBetweenGroups == 0)
            throw runtime_error
                ("Need at least 1 observation to compute variance.");

        /* compute the varianceBetweenGroups */
        double varianceBetweenGroups =
               sumOfSquaresBetweenGroups / dBetweenGroups; 

        /* return the varianceBetweenGroups */
        return(varianceBetweenGroups);

    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

   double MainEffectsAnalyzer3::getVarianceBetweenGroups
                (std::string factorColumn,
                 std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double varianceBetweenGroups = 
            getVarianceBetweenGroups
               (indexFactorColumn,
                indexResponseColumn);
        return(varianceBetweenGroups);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


   double MainEffectsAnalyzer3::getVarianceBetweenGroups
                (std::string factorColumn,
                 int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        double varianceBetweenGroups = 
            getVarianceBetweenGroups
               (indexFactorColumn,
                indexResponseColumn);
        return(varianceBetweenGroups);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

   double MainEffectsAnalyzer3::getVarianceBetweenGroups
                (int indexFactorColumn,
                 std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double varianceBetweenGroups = 
            getVarianceBetweenGroups
               (indexFactorColumn,
                indexResponseColumn);
        return(varianceBetweenGroups);

    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

   double MainEffectsAnalyzer3::getVarianceWithinGroups
                (std::string factorColumn,
                 std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double varianceWithinGroups = 
            getVarianceWithinGroups
               (indexFactorColumn,
                indexResponseColumn);
        return(varianceWithinGroups);
    }





    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

   double MainEffectsAnalyzer3::getVarianceWithinGroups
                (std::string factorColumn,
                int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        double varianceWithinGroups = 
            getVarianceWithinGroups
               (indexFactorColumn,
                indexResponseColumn);
        return(varianceWithinGroups);
    }




    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

   double MainEffectsAnalyzer3::getVarianceWithinGroups
                (int indexFactorColumn,
                 std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double varianceWithinGroups = 
            getVarianceWithinGroups
               (indexFactorColumn,
                indexResponseColumn);
        return(varianceWithinGroups);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/

    double MainEffectsAnalyzer3::getVarianceWithinGroups
                (int indexFactorColumn,
                 int indexResponseColumn){

        /* get the sumOfSquaresWithinGroups */
        double sumOfSquaresWithinGroups =
               getSumOfSquaresWithinGroups
                    (indexFactorColumn,
                     indexResponseColumn);

        /* get the dWithinGroups */
        double dWithinGroups =
               getDWithinGroups
                    (indexFactorColumn,
                     indexResponseColumn);

        /* do some error checking */
        if (dWithinGroups == 0)
            throw runtime_error
                ("d Within Groups must be greater than 0.");

        /* compute the varianceWithinGroups */
        double varianceWithinGroups =
               sumOfSquaresWithinGroups / dWithinGroups; 

        /* return the varianceWithinGroups */
        return(varianceWithinGroups);

    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getFdata
                (int indexFactorColumn,
                 int indexResponseColumn){

       /* compute Variance Between Groups */
       double varianceBetweenGroups = getVarianceBetweenGroups
              (indexFactorColumn,
               indexResponseColumn);

       /* compute Variance Within Groups */
       double varianceWithinGroups = getVarianceWithinGroups
              (indexFactorColumn,
               indexResponseColumn);

       /* do some error checking */
       if (varianceWithinGroups==0)
           throw runtime_error
                ("Variance Within Groups must be greater than zero.");

       /* compute fdata */
       double fdata = varianceBetweenGroups / varianceWithinGroups;

       /* return fdata */
       return(fdata);

    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getFdata
                (std::string factorColumn,
                 std::string responseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double fdata = getFdata
               (indexFactorColumn,
                indexResponseColumn);
        return(fdata);
    }


    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getFdata
                (std::string factorColumn,
                 int indexResponseColumn){

        int indexFactorColumn = toIndexInputColumn(factorColumn);
        double fdata = getFdata
               (indexFactorColumn,
                indexResponseColumn);
        return(fdata);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/


    double MainEffectsAnalyzer3::getFdata
                (int indexFactorColumn,
                 std::string responseColumn){

        int indexResponseColumn = toIndexInputColumn(responseColumn);
        double fdata = getFdata
               (indexFactorColumn,
                indexResponseColumn);
        return(fdata);
    }



    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*/
