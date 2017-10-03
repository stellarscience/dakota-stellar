/*******************************************************/
/* Copyright (c) 2009-2010 Sandia Corporation          */
/* All Rights Reserved                                 */
/*******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <time.h>


/*---- TO COMPILE, LINK WITH THIS OBJECT FILE. */
extern void  genDownloadPage (const char * const  szAddress);


/*---- NEED A LIMIT TO PREVENT BUFFER OVERFLOW. */
static const int  MAX_FORM_SIZE = 1000;
static const int  MAX_PROP_SIZE = 100;

static const char  REGFILE[] = "hopspack_files/hopspack_reg.txt";
static const char  HOMEURL[] = "http://www.sandia.gov/hopspack";


/** Return the decimal value corresponding to a hex digit.
 */
static int  toHexValue (const char  c)
{
    if ((c >= '0') && (c <= '9'))
        return( ((int) c) - ((int) '0') );
    if ((c >= 'A') && (c <= 'F'))
        return( 10 + ((int) c) - ((int) 'A') );
    if ((c >= 'a') && (c <= 'f'))
        return( 10 + ((int) c) - ((int) 'f') );

    return( -1 );
}


/** Unencode from URL encoding (very simple implementation).
 */
static void  urlUnencode (const char * const  szEncoded,
                                char * const  szUnencoded,
                          const int           nMaxChars)
{
    int  k;

    if (szEncoded == NULL)
    {
        sprintf (szUnencoded, "(null)");
        return;
    }

    int  i = 0;
    int  j = 0;
    while (szEncoded[i] != 0)
    {
        if (szEncoded[i] == '%')
        {
            if ((szEncoded[i+1] != 0) && (szEncoded[i+2] != 0))
            {
                k = 16 * toHexValue (szEncoded[i+1])
                       + toHexValue (szEncoded[i+2]);
                szUnencoded[j] = (char) k;
                i += 3;
                j++;
            }
            else    /*-- ERRONEOUS ENCODING, BUT KEEP GOING */
            {
                szUnencoded[j] = szEncoded[i];
                i++;
                j++;
            }
        }
        else if (szEncoded[i] == '+')
        {
            szUnencoded[j] = ' ';
            i++;
            j++;
        }
        else
        {
            szUnencoded[j] = szEncoded[i];
            i++;
            j++;
        }

        if (j >= nMaxChars)
            break;
    }

    szUnencoded[j] = 0;
    return;
}


/** Send an HTML error message if there is no form; for instance, if
 *  someone tries to download through a bookmark.
 */
static void  replyNoForm (void)
{
    printf ("Content-type: text/html%c%c", 10, 10);
    printf ("<html>%c<head><title>HOPSPACK Illegal Request</title></head>%c",
            10, 10);
    printf ("<body text='#000000'>%c", 10);
    printf ("<img src='%s/HOPSPACK_SNL_NNSA_100x734.png' alt='Logo'"
            "     border='0' height='100' width='734'>%c", HOMEURL, 10);
    printf ("</br></br>%c", 10);
    printf ("<hr>%c", 10);
    printf ("<h3>Illegal Request</h3>%c", 10);
    printf ("</br></br>%c", 10);
    printf ("Your options are:%c", 10);
    printf ("<ul>%c", 10);
    printf ("  <li><a href='%s/index.html'>Sign in</a> with an email address%c",
            HOMEURL, 10);
    printf ("  <li><a href='%s/registration.shtml'>Register</a> your email address%c",
            HOMEURL, 10);
    printf ("</ul>%c", 10);
    printf ("</body>%c</html>%c", 10, 10);
    fflush (stdout);
    return;
}


/** Send an HTML error message that the form is incorrect.
 */
static void  replyFormError (const char * const  szMessage)
{
    printf ("Content-type: text/html%c%c", 10, 10);
    printf ("<html>%c<head><title>Registration Failed</title></head>%c", 10, 10);
    printf ("<body text='#000000'>%c", 10);
    printf ("<img src='%s/HOPSPACK_SNL_NNSA_100x734.png' alt='Logo'"
            "     border='0' height='100' width='734'>%c", HOMEURL, 10);
    printf ("</br></br>%c", 10);
    printf ("<hr>%c", 10);
    printf ("<h3>Registration Failed</h3>%c", 10);
    printf ("%s%c", szMessage, 10);
    printf ("</br></br>%c", 10);
    printf ("Please use the 'back' button on your browser and correct the"
            " problem with your registration.%c", 10);
    printf ("</body>%c</html>%c", 10, 10);
    fflush (stdout);
    return;
}


/** Find the next form field, delimited by '&', and replace the delimiter with 0.
 *  Return the new current position to resume parsing.
 */
static int  delimitNextFormField (      char * const  szBuffer,
                                  const int           nPos)
{
    int  k;

    k = nPos;
    while ((szBuffer[k] != 0) && (szBuffer[k] != '&'))
        k++;

    if (szBuffer[k] != 0)
        szBuffer[k] = 0;

    return( k + 1 );
}


/** Return 1 if the desired property is found before the first '='.
 *  Return 0 if not.
 */
static int  parseFormProperty (      char * const  szBuffer,
                               const char * const  szPropertyName,
                                     char * const  szValue,
                               const int           nValueSize)
{
    int  k = 0;
    while ((szBuffer[k] != 0) && (szBuffer[k] != '='))
        k++;

    if (szBuffer[k] == 0)
        return( 0 );

    if (strncmp (szBuffer, szPropertyName, k) != 0)
        return( 0 );

    urlUnencode (&szBuffer[k + 1], szValue, nValueSize);
    szValue[nValueSize] = 0;

    return( 1 );
}


/** Return 1 if the input is valid, 0 if not.
 */
static int  parseAndValidate (      char * const  szUrlEncBuffer,
                                    char * const  szEmailAddr,
                                    char * const  szFullName,
                                    char * const  szCompany,
                                    char * const  szSector,
                                    char * const  szCountry,
                                    char * const  szRegion,
                                    char * const  szApp,
                                    char * const  szHow,
                                    char * const  szUserList,
                              const int           nMaxPropSize)
{
    /*---- PARSE EACH PROPERTY BEFORE UNENCODING, IN CASE THE VALUE
     *---- CONTAINS '&' OR '=' CHARACTERS. */

    int  nBufferLen = strlen (szUrlEncBuffer);
    int  nPos = 0;
    int  nNext;

    nNext = delimitNextFormField (szUrlEncBuffer, nPos);
    if (parseFormProperty (&szUrlEncBuffer[nPos], "emailaddr",
                           szEmailAddr, nMaxPropSize) == 0)
        return( 0 );

    nPos = nNext;
    nNext = delimitNextFormField (szUrlEncBuffer, nPos);
    if (parseFormProperty (&szUrlEncBuffer[nPos], "fullname",
                           szFullName, nMaxPropSize) == 0)
        return( 0 );

    nPos = nNext;
    nNext = delimitNextFormField (szUrlEncBuffer, nPos);
    if (parseFormProperty (&szUrlEncBuffer[nPos], "company",
                           szCompany, nMaxPropSize) == 0)
        return( 0 );

    nPos = nNext;
    nNext = delimitNextFormField (szUrlEncBuffer, nPos);
    if (parseFormProperty (&szUrlEncBuffer[nPos], "sector",
                           szSector, nMaxPropSize) == 0)
        return( 0 );

    nPos = nNext;
    nNext = delimitNextFormField (szUrlEncBuffer, nPos);
    if (parseFormProperty (&szUrlEncBuffer[nPos], "country",
                           szCountry, nMaxPropSize) == 0)
        return( 0 );

    nPos = nNext;
    nNext = delimitNextFormField (szUrlEncBuffer, nPos);
    if (parseFormProperty (&szUrlEncBuffer[nPos], "region",
                           szRegion, nMaxPropSize) == 0)
        return( 0 );

    nPos = nNext;
    nNext = delimitNextFormField (szUrlEncBuffer, nPos);
    if (parseFormProperty (&szUrlEncBuffer[nPos], "app",
                           szApp, nMaxPropSize) == 0)
        return( 0 );

    nPos = nNext;
    nNext = delimitNextFormField (szUrlEncBuffer, nPos);
    if (parseFormProperty (&szUrlEncBuffer[nPos], "how",
                           szHow, nMaxPropSize) == 0)
        return( 0 );

    if (nBufferLen <= (nNext - 1))
    {
        szUserList[0] = 'N';
        szUserList[1] = 0;
    }
    else
    {
        nPos = nNext;
        nNext = delimitNextFormField (szUrlEncBuffer, nPos);
        if (parseFormProperty (&szUrlEncBuffer[nPos], "userlist",
                               szUserList, nMaxPropSize) == 0)
            return( 0 );
    }

    return( 1 );
}


/** This program is invoked as a CGI module from the registration.shtml page,
 *  which submits a form of information.  It reads the information, parses
 *  into fields, saves to a file, and returns the download page.
 *
 *  Information is not checked for validity.  The web page ensures all
 *  required fields are nonempty, but little more than that.  If the same
 *  user registers twice, the second record is appended to the file so that
 *  both records are available.
 */
int  main (int argc, char *argv[])
{
    char    szUrlEncRegForm[MAX_FORM_SIZE + 1];
    char    szCopyEncRegForm[MAX_FORM_SIZE + 1];

    char    szRegForm[MAX_FORM_SIZE + 1];
    char    szEmailAddr[MAX_PROP_SIZE + 1];
    char    szFullName[MAX_PROP_SIZE + 1];
    char    szCompany[MAX_PROP_SIZE + 1];
    char    szSector[MAX_PROP_SIZE + 1];
    char    szCountry[MAX_PROP_SIZE + 1];
    char    szRegion[MAX_PROP_SIZE + 1];
    char    szApp[MAX_PROP_SIZE + 1];
    char    szHow[MAX_PROP_SIZE + 1];
    char    szUserList[MAX_PROP_SIZE + 1];

    int     i;
    char    c;
    FILE *  fp;
    time_t  rawTime;
    char    szTime[30];


    /*---- READ THE FORM PARAMETERS. */
    szUrlEncRegForm[0] = 0;
    fgets (szUrlEncRegForm, MAX_FORM_SIZE, stdin);
    if (strlen (szUrlEncRegForm) < 10)
    {
        replyNoForm();
        return( 0 );
    }
    szRegForm[MAX_FORM_SIZE] = 0;

    strncpy (szCopyEncRegForm, szUrlEncRegForm, MAX_FORM_SIZE);

    /*---- UNCOMMENT THIS FOR DEBUGGING, AND CREATE THE FILE AS ugo+w. */
    /*
    FILE * fpd = fopen ("hopspack_files/tbdreginput.txt", "a");
    fprintf (fpd, "!%s!\n", szUrlEncRegForm);
    fclose (fpd);
    */

    if (parseAndValidate (szUrlEncRegForm,
                          szEmailAddr,
                          szFullName,
                          szCompany,
                          szSector,
                          szCountry,
                          szRegion,
                          szApp,
                          szHow,
                          szUserList,
                          MAX_PROP_SIZE) == 0)
    {
        replyFormError ("Invalid form entry detected.");
        return( 0 );
    }

    /*---- EMAIL CANNOT HAVE ANY WHITESPACE. */
    i = 0;
    while (szEmailAddr[i] != 0)
    {
        c = szEmailAddr[i];
        if ((c == ' ') || (c < 32) || (c >= 127))
        {
            replyFormError ("No whitespace allowed in an email address.");
            return( 0 );
        }
        i++;
    }

    /*---- TIME STAMP HAS FORMAT ddd mmm DD HH:MM:SS YYYY. */
    time (&rawTime);
    asctime_r (localtime (&rawTime), szTime);
    szTime[strlen (szTime) - 1] = 0;

    /*---- ALTHOUGH ALL FIELDS ARE PARSED, SAVE IT IN THE ORIGINAL ENCODED FORM.
     *---- PARSING DOES VERIFY CORRECTNESS. */
    fp = fopen (REGFILE, "a");
    fprintf (fp, "%s %s\n", szEmailAddr, szTime);
    fprintf (fp, "%s\n", szCopyEncRegForm);
    fclose (fp);

    /*---- FINISH BY RETURNING THE DOWNLOAD PAGE. */
    genDownloadPage (szEmailAddr);

    return( 0 );
}
