/*******************************************************/
/* Copyright (c) 2009-2010 Sandia Corporation          */
/* All Rights Reserved                                 */
/*******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


static const char  CONTACT[] = "tplante@sandia.gov";
static const char  LICFILE[] = "hopspack_license_text.txt";
static const char  HDRFILE[] = "hopspack_snl_header.txt";
static const char  HOMEURL[] = "http://www.sandia.gov/hopspack";
static const char  BANNERSSI[] = "/web/hopspack/sandia_banner_ssi.txt";
static const char  FOOTERSSI[] = "/web/hopspack/sandia_footer_ssi.txt";


/** Dynamically generate Sandia compliant header for an HTML page.
 *  Return -1 if there was an error, +1 if successful.
 */
int  genSandiaCompliantHeader (const char * const  szWebPageTitle,
                               const char * const  szJavaScript)
{
    FILE *  fp;
    char    szNextLine[200 + 1];
    char *  pDel;
    int     k;

    fp = fopen (HDRFILE, "r");
    if (fp == NULL)
    {
        printf ("Content-type: text/html%c%c", 10, 10);
        printf ("<html>%c<head><title>Error</title></head>%c", 10, 10);
        printf ("<body>%c", 10);
        printf ("HOPSPACK download script internal error 50%c", 10);
        printf ("<br></br>%c", 10);
        printf ("This should not be happening.%c", 10);
        printf ("Please contact <a href='mailto:%s'>%s</a>%c",
                CONTACT, CONTACT, 10);
        printf ("</body>%c</html>%c", 10, 10);
        fflush (stdout);

        return( -1 );
    }

    /*---- PRINT THE HEADER. */
    printf ("Content-type: text/html%c%c", 10, 10);
    while (feof (fp) == 0)
    {
        if (fgets (szNextLine, 200, fp) != NULL)
        {
            /*---- SUBSTITUTE FOR THE TITLE, IF FOUND. */
            char *  pDel = strchr (szNextLine, '%');
            if (   (pDel != NULL)
                && (*(pDel + 1) == 'T')
                && (*(pDel + 2) == 'I')
                && (*(pDel + 3) == 'T')
                && (*(pDel + 4) == 'L')
                && (*(pDel + 5) == 'E')
                && (*(pDel + 6) == '%') )
            {
                k = pDel - szNextLine;
                szNextLine[k] = 0;
                printf ("%s%s%s",
                        szNextLine, szWebPageTitle, &szNextLine[k + 7]);
            }
            else
            {
                printf ("%s", szNextLine);
            }
        }
    }
    fflush (stdout);
    fclose (fp);

    /*---- INSERT ANY JAVASCRIPT. */
    if (szJavaScript != NULL)
    {
        printf ("%c%s%c", 10, szJavaScript, 10);
        fflush (stdout);
    }

    /*---- CLOSE THE HEADER AND PRINT THE BANNER. */
    printf ("</head>%c", 10);
    printf ("%c", 10);
    printf ("<body>%c", 10);
    printf ("<!-- sandia blue banner server-side include -->%c", 10);
    fflush (stdout);

    fp = fopen (BANNERSSI, "r");
    if (fp == NULL)
    {
        printf ("Content-type: text/html%c%c", 10, 10);
        printf ("<html>%c<head><title>Error</title></head>%c", 10, 10);
        printf ("<body>%c", 10);
        printf ("HOPSPACK download script internal error 51%c", 10);
        printf ("<br></br>%c", 10);
        printf ("This should not be happening.%c", 10);
        printf ("Please contact <a href='mailto:%s'>%s</a>%c",
                CONTACT, CONTACT, 10);
        printf ("</body>%c</html>%c", 10, 10);
        fflush (stdout);

        return( -1 );
    }
    while (feof (fp) == 0)
    {
        if (fgets (szNextLine, 200, fp) != NULL)
            printf ("%s", szNextLine);
    }
    fflush (stdout);
    fclose (fp);

    return( 1 );
}


/** Dynamically generate Sandia compliant footer for an HTML page.
 *  Return -1 if there was an error, +1 if successful.
 */
int  genSandiaCompliantFooter (void)
{
    FILE *  fp;
    char    szNextLine[200 + 1];

    fp = fopen (FOOTERSSI, "r");
    if (fp == NULL)
    {
        printf ("Content-type: text/html%c%c", 10, 10);
        printf ("<html>%c<head><title>Error</title></head>%c", 10, 10);
        printf ("<body>%c", 10);
        printf ("HOPSPACK download script internal error 52%c", 10);
        printf ("<br></br>%c", 10);
        printf ("This should not be happening.%c", 10);
        printf ("Please contact <a href='mailto:%s'>%s</a>%c",
                CONTACT, CONTACT, 10);
        printf ("</body>%c</html>%c", 10, 10);
        fflush (stdout);

        return( -1 );
    }

    printf ("<!-- sandia footer server-side include -->%c", 10);

    while (feof (fp) == 0)
    {
        if (fgets (szNextLine, 200, fp) != NULL)
            printf ("%s", szNextLine);
    }
    fflush (stdout);
    fclose (fp);

    return( 1 );
}


/** Dynamically generate and send the HTML download page.
 */
void  genDownloadPage (const char * const  szAddress)
{
    FILE *  fp;
    char  szQuery[1000 + 1];
    char  szJScript[2000 + 1];
    char  szNextLine[100 + 1];


    /*---- EMBED szAddress IN THE JAVASCRIPT FUNCTION.
     *---- GET THE REFERRER FROM THIS REQUEST BECAUSE THE LOGGER WILL ALWAYS
     *---- SEE THIS PAGE AS ITS REFERRER.  IF THIS FIELD IS EMPTY IT COULD
     *---- MEAN THE USER DOWNLOADED FROM A BOOKMARK INSTEAD OF THRU THE ADDRESS
     *---- CHECK, OR MAYBE THEIR BROWSER HAS REFERRER DISABLED. */
    sprintf (szQuery, "addr=%s&referer=%s",
             szAddress,
             getenv ("HTTP_REFERER"));

    /*---- APPENDING THE DOCUMENT ELEMENT WILL GENERATE A CALL TO THE
     *---- SERVER.  SERVER REPLY IS IGNORED, BUT SOME REPLY IS EXPECTED OR
     *---- APACHE WILL LOG THE ERROR "Premature end of script headers". */
    sprintf (szJScript, "%s%c %s%c%s%c%s%s%s%c %s%c%s%c%s%c%s%c %s%c",
             "<script type='text/javascript'> <!--", 10,
             "function checkdnld (szDownload) {", 10,
             "  var  sc = document.createElement('script');", 10,
             "  sc.src = 'hopspack_logdnld?",
             szQuery,
             "' + '&download=' + szDownload;", 10,
             "  sc.type = 'text/javascript';", 10,
             "  document.getElementsByTagName('head')[0].appendChild (sc);", 10,
             "  document.getElementsByTagName('head')[0].removeChild (sc);", 10,
             " }", 10,
             "// --> </script>", 10);

    if (genSandiaCompliantHeader ("HOPSPACK - Download Page", szJScript) != 1)
        return;


    printf ("%c", 10);
    printf ("<div style='font:normal 1.25em Verdana, Arial, Helvetica, sans-serif;%c", 10);
    printf ("            padding:20px; border:0;'>%c", 10);

    printf ("<table>%c", 10);
    printf ("  <tr>%c%c", 10, 10);
    printf ("    <td>%c", 10);
    printf ("      <!-- Left side content -->%c", 10);
    printf ("      <h2 style='color:rgb(0%,0%,0%)'>HOPSPACK Registration</h2>%c", 10);

    fp = fopen (LICFILE, "r");
    if (fp == NULL)
    {
        printf ("HOPSPACK download script internal error 55%c", 10);
        printf ("</br>%c", 10);
        printf ("Please contact <a href='mailto:%s'>%s</a>%c",
                CONTACT, CONTACT, 10);
        printf ("</body>%c</html>%c", 10, 10);
        fflush (stdout);
        return;
    }
    while (feof (fp) == 0)
    {
        if (fgets (szNextLine, 100, fp) != NULL)
            printf ("%s", szNextLine);
    }
    fclose (fp);

    printf ("    </td>%c%c", 10, 10);

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
    printf ("</div>%c", 10);


    genSandiaCompliantFooter();

    printf ("</body>%c</html>%c", 10, 10);
    fflush (stdout);
    return;
}
