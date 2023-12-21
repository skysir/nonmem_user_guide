/*
**      nmhelp_win.c
**      help display utility for Windows
**      Written by JM Gries for the NONMEM Project Group - 1996
**      Modified by RJ Bauer for WIndows XP, August 2009
*/



#include <stdio.h>

#include <stdlib.h>

#include <string.h>

#include <ctype.h>



#define antislash '\\'

#define SCCS "@#@ Written by JM Gries for the NONMEM Project Group"

#define DEBUG( a )  { fprintf(stderr, "debug is %s\n", a ); }



        #define bin ""

        #define pager "more"

        #define printer "prn:"

        #define gotsome "~htmp"

        #define gotsome2 "~htmp2"

        #define defdir "c:\\nmv\\help\\"

        void Stolower( char * );

        int CountTokens(char * );

        char ** CreateArray(int ,int );

        int FreeArray(char **, int );

        int ReadTokens( char *, int);

        void PrintMenu( int , char **, char *, int  );

        int GetChoice( int , int );

        void undosify( char *);

        void dosify( char *);

        void CleanString( char * );

        void FindIt( char *, char * );

        int IsFileExist( char *, char * );

        int GetEnvString( char * );
        void system2(char *);





typedef struct {

        char **table1, **table2;

        } TWOTABLES;



static TWOTABLES general;



int main( argc, argv )

int argc;

char ** argv;

{

        char tmp[1024+1],lf[128+1];

        char tmp2[1024+1],sccs[128+1];

        int lname,count,oldans, ans,i;

        char ** tablenames, ** tablekw;




        char word[15]="";

        char space[15]=" ";




        strcpy(sccs, SCCS);

        FindIt(argv[0],lf);



        if(argc == 1)

        {

                sprintf( tmp, "%s%s < %slisthelp",bin, pager, lf);

                system2( tmp );

                exit(0);

        }



        argv++;

        argc--;



        while( argc-- )

        {

                if( strcasecmp( *argv, "HELP")==0 || strcasecmp( *argv, "NMHELP")==0 )

                {

                fprintf(stderr,"Point2\n");

                        sprintf( tmp, "%s < %shelphelp", pager, lf);

                        system2( tmp );

                        break;

                }



                if( strcasecmp( *argv, "-W" )==0 )

                {

                        strcpy( word, "-w" );

                        strcpy( space, "" );

                        argv++;

                        continue;

                }



                if( strcasecmp( *argv, "-A" )==0 )

                {

                        strcpy( space, "" );

                        argv++;

                        continue;

                }



                if( strcasecmp( *argv, "-R" )==0 )

                {

                        strcpy( space, " " );

                        strcpy( word, "" );

                        argv++;

                        continue;

                }



                strcpy( tmp2, *argv );

                Stolower( tmp2 );

                argv++;

                if( strcmp( tmp2, "\"")==0 || strcmp( tmp2, ":")==0 || strcmp( tmp2, "@" )==0 || strcmp( tmp2, "#" )==0 || strcmp( tmp2, ";" )==0 || strcmp( tmp2, "*" )==0 )

                {

                fprintf(stderr,"Point3\n");

                        sprintf( tmp, "%sgrep \"\\%s\" %sindex | %sgawk \"BEGIN {FS=\\\"~\\\"};{print $1,\\\" ~ \\\",$2}\" | sort > %s", bin, tmp2, lf, bin, gotsome2 );

                        system2( tmp );

                        sprintf( tmp, "%sgawk \"{if(last!=$0)print;last=$0}\" %s > %s", bin, gotsome2, gotsome );

                        system2( tmp );


                        sprintf( tmp, "del %s > nul", gotsome2 );


                        system2(tmp);

                }

                else

                {

                        while( strchr(tmp2,'-') != NULL )

                                *(strchr(tmp2,'-'))='.';

                        lname=strlen(tmp2);

                        if( lname < 3 )

                        {

                                sprintf( tmp, "%sgrep -i -w %s %sindex | %sgawk \"BEGIN {FS=\\\"~\\\"};{print $1,\\\" ~ \\\",$2}\" | sort > %s", bin, tmp2, lf, bin, gotsome2 );

                                system2( tmp );

                                sprintf( tmp, "%sgawk \"{if(last!=$0)print;last=$0}\" %s > %s", bin, gotsome2, gotsome );

                                system2( tmp );


                                sprintf( tmp, "del %s > nul", gotsome2 );


                                system2(tmp);

                        }

                        else

                        {

                                sprintf( tmp, "%sgrep -i %s \"%s%s\" %sindex | %sgawk \"BEGIN {FS=\\\"~\\\"};{print $1,\\\" ~ \\\",$2}\" | sort > %s", bin, word, space, tmp2, lf, bin, gotsome2 );

                                system2( tmp );

                                sprintf( tmp, "%sgawk \"{if(last!=$0)print;last=$0}\" %s > %s", bin, gotsome2, gotsome );

                                system2( tmp );


                                sprintf( tmp, "del %s > nul", gotsome2 );


                                system2(tmp);

                        }

                }

                count = CountTokens( gotsome );

                if( count == 0 )

                {

                        fprintf( stdout, "No help available for %s\n", tmp2 );

                        if( tmp2[strlen(tmp2)-1] == 's' )

                        {

                                tmp2[strlen(tmp2)-1]=0;

                                fprintf( stdout, "Try nmhelp %s\n",tmp2 );

                        }

                        continue;

                }


                i = ReadTokens ( gotsome, count);

                if( i!= count )

                {

                        fprintf(stderr,"Error: tokens# different from line#\n");

                        exit(-1);

                }

                tablenames = general.table1;

                tablekw = general.table2;



                oldans=0;

                do {

                        PrintMenu(count, tablekw, tmp2,(oldans<count) ? (oldans+1) : 0);

                        ans=GetChoice(count,oldans);

                        if(ans == 0) break;

                        if(ans == 255) continue;

                        if(ans<0) {


                                sprintf( tmp, "copy %s%s %s", lf, tablenames[(-ans)-1],printer );


                                system2( tmp );

                                oldans=-ans;

                        } else {


                                sprintf( tmp, "%s < %s%s", pager, lf, tablenames[ans-1] );


                                system2( tmp );

                                oldans=ans;

                        }

                }while( ans != 0);



                FreeArray(tablenames,count);

                FreeArray(tablekw,count);

                general.table1=NULL;

                general.table2=NULL;


                sprintf( tmp, "del %s > nul", gotsome );


                system2(tmp);

        }

}



void Stolower( string )

char *string;

{

        char *ptr=string;

        while( *ptr )

        {

                *ptr=tolower( *ptr );

                ptr++;

        }

}



int CountTokens( filename )

char *filename;

{

        FILE *infile;

        char tmp[1024+1];

        int line = 0;

        char *got;



        if( (infile=fopen( filename, "r" )) == NULL )

        {

                fprintf(stderr," Unable to open %s. ABORT\n", filename );

                exit( -1 );

        }



        do {

                got=fgets( tmp, 1024, infile);

                if( got!= NULL )

                        line++;

        } while( !feof( infile ) );

        fclose( infile );

        return( line );

}



char ** CreateArray( nbtokens, sizetoken)

int nbtokens,sizetoken;

{

        char ** tmp;

        int i;

        tmp = (char **)calloc( (unsigned)nbtokens, sizeof( char *) );

        for( i=0;i<nbtokens;i++)

                tmp[i] = (char *)malloc( sizetoken);

        return(tmp);

}



int FreeArray( table, nbtokens)

char ** table;

int nbtokens;

{

        int i;

        for( i=0;i<nbtokens;i++)

                free(table[i]);

        free( table);

        return(0);

}



int ReadTokens ( filename, nbtokens)

char *filename;

int nbtokens;

{

        FILE *infile;

        char tmp[1024+1],filen[128];

        int line = 0,i;

        char * got,* got2;

        char ** tablenames;

        char ** tablekw;


        tablenames = CreateArray(nbtokens,30);

        if( (infile=fopen( filename, "r" )) == NULL )

        {

                fprintf(stderr," Unable to open %s. ABORT\n", filename );

                FreeArray( tablenames, nbtokens);

                exit( -1 );

        }

        tablekw = CreateArray( nbtokens, 128);

        while( !feof( infile ) )

        {

                got=fgets( tmp,128,infile);

                if( got != NULL)

                {

                        i = sscanf(tmp,"%12s  ~", filen);

                        if(i!=1) {

/* prbm */

                                FreeArray( tablenames, nbtokens);

                                FreeArray( tablekw, nbtokens);

                                perror("ERROR: file name missing");

                        }

                        strcpy( tablenames[line], filen);

                        got2 = strchr(tmp,'~');

                        if( got2==NULL) {

/* error */

                                strncpy(tablekw[line],"ERROR MISSING KW",127);

                                perror("KW error");

                        }

                        else {

                                got2 += 2;

                                strncpy(tablekw[line], got2, 31);

                                (tablekw[line])[32]=0;

                                got2 = strchr(tablekw[line],'\n');

                                if(got2 !=NULL) *got2 = 0; /* RJB, August 2009 */

                        }

                        line++;

                        if( line == nbtokens ) break;

                }

        }

        general.table1 = tablenames;

        general.table2 = tablekw;



        fclose( infile );

        return( line );

}



void PrintMenu( nbtokens, tablekw, topic, nextone )

int nbtokens, nextone;

char ** tablekw, *topic;

{

        int i,half=(nbtokens+1)/2,odd=nbtokens%2;



        fprintf(stdout, "Help for %s is available in the following %d document(s):\n", topic, nbtokens);



        if( half != 1) {

                for(i=1;i<=(half-1);i++)

                        fprintf(stdout,"%d\t%-30s\t%d\t%-30s\n",i,tablekw[i-1],i+half,tablekw[i+half-1]);

        }

        if(odd)

                fprintf(stdout,"%d\t%-30s\n",half,tablekw[half-1]);

        else

                fprintf(stdout,"%d\t%-30s\t%d\t%-30s\n",half,tablekw[half-1],half+half,tablekw[half+half-1]);



        fprintf(stdout, "Your choices:\n");

        fprintf(stdout, "  enter p number to print a document (e.g., p1 or p 1):\n");

        fprintf(stdout, "  enter the number of the document you wish to see\n");

        fprintf(stdout, "  enter 0 or q to quit, return for next choice (%d)\n",nextone);

}



int GetChoice(count,oldans)

int count,oldans;

{

        char tmp[80];

        int print=0,quit=0,nb=0,val=0;



        fgets(tmp,79,stdin);



        Stolower(tmp);



        if( strchr(tmp,'q')!=NULL) quit=1;

        if( strchr(tmp,'p')!=NULL) print=1;



        CleanString( tmp);



        nb=sscanf(tmp,"%d",&val);

        if(nb==0 && quit==0 ) return(255);

        if(nb==-1) val=(oldans+1)>count ? 0: oldans+1;



        if( val==0 || quit) return(0);

        if( val < 0) return(255);

        if( val > count) {

                return( 255);

        } else {

                if(print) return( -val);

                else return(val);

        }

}



void dosify( string )

char *string;

{

        char *got;



        while( (got=strchr(string,'/'))!= NULL )

                *got = '\\';

}



void undosify( string )

char *string;

{

        char *got;



        while( (got=strchr(string,antislash))!= NULL )

                *got = '/';

}



void CleanString( string )

char *string;

{

        char *pchar=string;



        while( *(pchar)!= 0 ) {

                if( isdigit(*pchar) == 0 )

                        *pchar = ' ';

                pchar ++;

        }

}



void FindIt( prgmname, dir )

char *prgmname, *dir;

{

        char where[128+1], *got;



        strcpy(where,prgmname);


        undosify(where);


        got=strrchr(where,'/');

        if(got!=NULL) {

                *(++got)=0;

        }

        else {

                strcpy(where,"./");

        }

        

        if( IsFileExist( "index", "./" ) == 1 ) {

                 strcpy(dir,"./");

        }

        else if( IsFileExist( "index", where ) == 1 ) {

                strcpy( dir, where );

        }

        else {

                if( GetEnvString( where ) == 1)

                {

                        if( IsFileExist( "index", where ) == 1 ) {

                                strcpy( dir, where );

                        }

                }

                else if( IsFileExist( "index", defdir ) == 0 ) {

                        fprintf(stderr, "ERROR cannot find the index file. Aborting\n");

                        exit(-1);

                } else strcpy( dir, defdir );

        }


        dosify( dir );


}



int IsFileExist( filename, where )

char *filename, *where;

{

        char tmp[1024];

        FILE *in;



        strcpy( tmp, where );

        strcat( tmp, filename );




        dosify( tmp );


        if( (in = fopen( tmp, "r")) == NULL )

                return( 0 );

        else {

                fclose( in );

                return(1);

        }

}



int GetEnvString( where )

char *where;

{

        char *got;



        if( (got = getenv("NMHELP")) == NULL )

        {

                strcpy( where, "");

                return( 0);

        }

        strcpy( where, got);

        if( where[strlen(where)-1] != '/' && where[strlen(where)-1] != '\\' )

                strcat( where, "/");




        dosify( where );


        return(1);

}

void system2(sline) 
     char *sline;
{
/*     
     printf(sline);
     printf("\n");
*/
     system(sline);
     return;
}
   
      