/*******************************************************/
/* Copyright (c) 2009-2010 Sandia Corporation          */
/* All Rights Reserved                                 */
/*******************************************************/

#include <stdio.h>
#include <string.h>


/*---- TO COMPILE, LINK WITH THIS OBJECT FILE. */
extern void  genDownloadPage (const char * const  szAddress);


/*---- NEED A LIMIT TO PREVENT BUFFER OVERFLOW. */
static const int  MAX_EMAIL_ADDR_SIZE =  200;
static const int  MAX_REG_DATA_SIZE   = 1000;

static const char  CONTACT[] = "tplante@sandia.gov";
static const char  REGFILE[] = "hopspack_files/hopspack_reg.txt";
static const char  HOMEURL[] = "http://www.sandia.gov/hopspack";


/** Send a generic HTML error message for unlikely internal errors.
 */
static void  replyInternalError (const int  nErrNum)
{
    printf ("Content-type: text/html%c%c", 10, 10);
    printf ("<html>%c<head><title>Error</title></head>%c", 10, 10);
    printf ("<body>%c", 10);
    printf ("HOPSPACK download script internal error %d%c", nErrNum, 10);
    printf ("<br></br>%c", 10);
    printf ("This should not be happening.%c", 10);
    printf ("Please contact <a href='mailto:%s'>%s</a>%c", CONTACT, CONTACT, 10);
    printf ("</body>%c</html>%c", 10, 10);
    fflush (stdout);
    return;
}


/** Send an HTML error message that the address is not registered.
 */
static void  replyNotRegistered (const char * const  szAddress)
{
    genSandiaCompliantHeader ("HOPSPACK - Not Registered", NULL);

    printf ("%c", 10);
    printf ("<div style='font:normal 1.25em Verdana, Arial, Helvetica, sans-serif;%c", 10);
    printf ("            padding:20px; border:0;'>%c", 10);

    printf ("<table>%c", 10);
    printf ("  <tr>%c%c", 10, 10);
    printf ("    <td>%c", 10);
    printf ("      <!-- Left side content -->%c", 10);
    printf ("      <h2 style='color:rgb(0%,0%,0%)'>Not Registered</h2>%c", 10);

    printf ("      <p><br>%c", 10);
    printf ("      The email address '%s' has not been registered.%c",
            szAddress, 10);
    printf ("      </p>%c", 10);
    printf ("      <p><br>%c", 10);
    printf ("      Your options are:%c", 10);
    printf ("      <div style='padding-left:20px'>%c", 10);
    printf ("      <ul>%c", 10);
    printf ("        <li><a href='%s/index.html'>Go back</a> and try again%c",
            HOMEURL, 10);
    printf ("        <li><a href='%s/registration.shtml'>Register</a> your email address%c",
            HOMEURL, 10);
    printf ("      </ul>%c", 10);
    printf ("      </div>%c", 10);
    printf ("      </p>%c", 10);
    printf ("    </td>%c%c", 10, 10);

    printf ("    <td width='200'>&nbsp;</td>%c%c", 10, 10);

    printf ("    <td style='padding-left:10px;vertical-align:top'>%c", 10);
    printf ("      <!-- Right side logo and nav -->%c", 10);
    printf ("      <img src='%s/HOPSPACK_Logo_120x118.png'%c", HOMEURL, 10);
    printf ("           alt='HOPSPACK logo' border='0' height='120' width='118'>%c", 10);
    printf ("      <p><br>%c", 10);
    printf ("         <strong>Page Contact</strong><br>%c", 10);
    printf ("         Principle Member Technical Staff<br>%c", 10);
    printf ("         <a href='mailto:%s'>Todd Plantenga</a><br>%c", CONTACT, 10);
    printf ("         (925) 294-3486%c", 10);
    printf ("      </p>%c", 10);
    printf ("      <p><br>%c", 10);
    printf ("        <strong>Related Links</strong><br>%c", 10);
    printf ("        <a href='https://software.sandia.gov/trac/hopspack/'>HOPSPACK Wiki</a>%c", 10);
    printf ("      </p>%c", 10);
    printf ("    </td>%c%c", 10, 10);
    printf ("  </tr>%c", 10);
    printf ("</table>%c%c", 10, 10);
    printf ("<p><br></p>%c%c", 10, 10);
    printf ("</div>%c%c", 10, 10);

    genSandiaCompliantFooter ();
    printf ("</body>%c</html>%c", 10, 10);
    fflush (stdout);
    return;
}


/** Send an HTML error message if there is no form; for instance, if
 *  someone tries to download through a bookmark.
 */
static void  replyNoForm (void)
{
    genSandiaCompliantHeader ("HOPSPACK - Illegal Request", NULL);

    printf ("%c", 10);
    printf ("<div style='font:normal 1.25em Verdana, Arial, Helvetica, sans-serif;%c", 10);
    printf ("            padding:20px; border:0;'>%c", 10);

    printf ("<table>%c", 10);
    printf ("  <tr>%c%c", 10, 10);
    printf ("    <td>%c", 10);
    printf ("      <!-- Left side content -->%c", 10);
    printf ("      <h2 style='color:rgb(0%,0%,0%)'>Illegal Request</h2>%c", 10);
    
    printf ("      <p><br>%c", 10);
    printf ("      Your options are:%c", 10);
    printf ("      <div style='padding-left:20px'>%c", 10);
    printf ("      <ul>%c", 10);
    printf ("        <li><a href='%s/index.html'>Sign in</a> with an email address%c",
            HOMEURL, 10);
    printf ("        <li><a href='%s/registration.shtml'>Register</a> your email address%c",
            HOMEURL, 10);
    printf ("      </ul>%c", 10);
    printf ("      </div>%c", 10);
    printf ("      </p>%c", 10);
    printf ("    </td>%c%c", 10, 10);

    printf ("    <td width='200'>&nbsp;</td>%c%c", 10, 10);

    printf ("    <td style='padding-left:10px;vertical-align:top'>%c", 10);
    printf ("      <!-- Right side logo and nav -->%c", 10);
    printf ("      <img src='%s/HOPSPACK_Logo_120x118.png'%c", HOMEURL, 10);
    printf ("           alt='HOPSPACK logo' border='0' height='120' width='118'>%c", 10);
    printf ("      <p><br>%c", 10);
    printf ("         <strong>Page Contact</strong><br>%c", 10);
    printf ("         Principle Member Technical Staff<br>%c", 10);
    printf ("         <a href='mailto:%s'>Todd Plantenga</a><br>%c", CONTACT, 10);
    printf ("         (925) 294-3486%c", 10);
    printf ("      </p>%c", 10);
    printf ("      <p><br>%c", 10);
    printf ("        <strong>Related Links</strong><br>%c", 10);
    printf ("        <a href='https://software.sandia.gov/trac/hopspack/'>HOPSPACK Wiki</a>%c", 10);
    printf ("      </p>%c", 10);
    printf ("    </td>%c%c", 10, 10);
    printf ("  </tr>%c", 10);
    printf ("</table>%c%c", 10, 10);
    printf ("<p><br></p>%c%c", 10, 10);
    printf ("</div>%c%c", 10, 10);

    genSandiaCompliantFooter ();
    printf ("</body>%c</html>%c", 10, 10);
    fflush (stdout);
    return;
}


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
                                char * const  szUnencoded)
{
    int  k;

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
    }

    szUnencoded[j] = 0;
    return;
}


/** This program is invoked as a CGI module from the index.html download page.
 *  It reads the email address, checks for registration, and returns
 *  either an error page (no registration) or the download page.
 *
 *  Since the server usually runs this as a user with limited permissions,
 *  the registration file must be read/writable for all.
 */
int  main (int argc, char *argv[])
{
    char    szRequest[MAX_EMAIL_ADDR_SIZE + 1];
    char    szReqAddr[MAX_EMAIL_ADDR_SIZE + 1];
    char    szNextAddr[MAX_EMAIL_ADDR_SIZE + 1];
    char    szNextLine[MAX_REG_DATA_SIZE + 1];
    int     nAddrLength;
    int     k;
    FILE *  fp;


    /*---- READ THE FORM PARAMETERS. */
    szRequest[0] = 0;
    fgets (szRequest, MAX_EMAIL_ADDR_SIZE, stdin);
    if (strlen (szRequest) < 10)
    {
        replyNoForm();
        return( 0 );
    }

    /*---- FIND THE EMAIL ADDRESS. */
    if (strncmp (szRequest, "emailaddr=", 10) != 0)
    {
        replyInternalError (2);
        return( 0 );
    }
    urlUnencode (&szRequest[10], szReqAddr);
    nAddrLength = strlen (szReqAddr);

    if (nAddrLength <= 0)
    {
        replyInternalError (3);
        return( 0 );
    }

    /*---- CHECK IF THE ADDRESS IS REGISTERED. */
    fp = fopen (REGFILE, "r");
    if (fp == NULL)
    {
        replyInternalError (4);
        return( 0 );
    }
    while (feof (fp) == 0)
    {
        fgets (szNextAddr, MAX_EMAIL_ADDR_SIZE, fp);
        fgets (szNextLine, MAX_REG_DATA_SIZE, fp);

        /*---- REMOVE TRAILING EOL CHARACTERS FROM THE STORED ADDRESS. */
        k = strlen (szNextAddr);
        while (k > 0)
        {
            if ((szNextAddr[k-1] != 10) && (szNextAddr[k-1] != 13))
                break;
            k--;
        }
        szNextAddr[k] = 0;

        /*---- FIND THE FIRST BLANK; THE REST IS A TIMESTAMP. */
        k = 0;
        while ((szNextAddr[k] != 0) && (szNextAddr[k] != ' '))
            k++;
        szNextAddr[k] = 0;

        /*---- TEST IF THIS IS THE ADDRESS. */
        if (   (nAddrLength == k)
            && (strncmp (szReqAddr, szNextAddr, nAddrLength) == 0))
        {
            /*---- ADDRESS IS REGISTERED. */
            fclose (fp);
            genDownloadPage (szReqAddr);
            fflush (stdout);
            return( 0 );
        }
    }
    fclose (fp);

    /*---- NOT REGISTERED.  COULD LOG THESE FAILED ATTEMPTS. */
    replyNotRegistered (szReqAddr);

    fflush (stdout);
    return( 0 );
}
