// This is a class to provide a standard interface for loading 
// parameter files and writing headers.

// However, you must extend this class to provide the core functionality
// for the primary code.  The derived class should contain all of the
// keyword variables that one wants, as well as a constructor to set 
// the defaults on those variables, a function assign_key(key,value)
// to set those variables, and functions to set up any derived values and
// to print.

#define MAXCOMMENT 65536
#define MAXLINE 1000

class Header {
  protected:
    char comment[MAXCOMMENT];  // A general text field
    virtual int assign_key(char *key, char *value) { return 0; }
        // This will decide what the keys mean
        // Must overload this function

  private:
    int commentsize;
    int append_comment_line(char *line,char *endcomment) {
        if (strncmp(line,endcomment,strlen(endcomment))==0) {
            // We've found the end of the comment block.  Do not append.
            // Return 0 because we're not in a comment anymore
            return 0;
        }
        // Otherwise, we're ok to try to add
        if (strlen(line)+commentsize>MAXCOMMENT) {
            // Yikes: too big!  We will generate a warning, not append, 
	    // but continue reading.
            fprintf(stderr,"Warning: Comment block size exceeded.  Truncating comment.  Increase MAXCOMMENT\n");
        } else {
            strcpy(comment+commentsize,line);
            commentsize+=strlen(line);
        }
        return 1;  // Still in comment
    }

    void strip_whitespace(char **s) {
        // This moves *s up to point to the first non-whitespace character.
        // Also replaces any trailing whitespace with \0
        char *p;
        while (isspace(**s)) (*s)++;
        for (p=*s;*p!='\0';p++);   // Point to the end
        while (p>*s&&isspace(p[-1])) *(--p)='\0';
    }

    void convert_to_lowercase(char *s) {
        char *p;
        for (p=s;*p!='\0';p++) *p = tolower(*p);
    }

    void parseline(char *line,char **_key,char **_value) {
        // Strip comments and white space, then try to parse into key = value;
        // Return key and/or value = NULL if this parsing fails.
        // This routine will alter 'line'.
        int j;
        char *key,*value;
        // Strip out comments that begin with # (but not \#)
        for (j=0;line[j]!='\0'&&j<MAXLINE;j++)
            if (line[j]=='#'&&(j==0||line[j-1]!='\\')) line[j] = '\0';
        // Now look for the first =
        for (value=line;*value!='\0';value++)
            if (*value=='=') { *value='\0'; value++; break; }
        // value either points to the beginning of a second string
        // or it points to the \0 of the first one.
        key = line;        strip_whitespace(&key);
        convert_to_lowercase(key);    // We force lower case key comparison
        strip_whitespace(&value);
        if (*value=='\0') { *_value = NULL; return; } else { *_value=value; }
        if (*key=='\0') { *_key = NULL; return; } else { *_key=key; }
        return;
    }

  public:
    Header() { commentsize = 0; }
    ~Header() { return; }

    int append_file_to_comments(char filename[]) {
        // Copy this whole file into the comments.
        // Return 0 if ok.  Add a \n just to be sure.
        FILE *fp;
        char buffer[MAXCOMMENT];
        char no_match[50] = "nothing to see here";
        int size;
        fp = fopen(filename,"r");        if (fp==NULL) {
            fprintf(stderr,"Warning: File to be included in comments not found.\n"); return 1;
        }
        // Write the file name
        sprintf(buffer,"\nIncluding file: %s\n",filename);
        append_comment_line(buffer,no_match);
        // Copy the whole file into the comment string
        size = fread(buffer,sizeof(char),MAXCOMMENT-2,fp);
        if (size>MAXCOMMENT-commentsize) {
            size = MAXCOMMENT-commentsize-2;
        }
        buffer[size]='\n'; buffer[size+1]='\0';
        strcpy(comment+commentsize,buffer);
        commentsize+=(size+1);
        fclose(fp);
        return 0;
    }

    int load(char *input_filename) {
        // Read the input parameter file and stuff it into the Parameter instance.
        // Return 0 if all is well, 1 if this failed.
        FILE *fp;
        char line[1000],*key=NULL,*value=NULL;
        char endcomment[100];
        int incomment;
        //
        fp = fopen(input_filename,"r");
        if (fp==NULL) {
            fprintf(stderr,"File %s not found\n", input_filename); return 1;
        }
        incomment = 0; endcomment[0] = '\0';
        //
        while (fgets(line,MAXLINE,fp)!=NULL) {
            line[MAXLINE-1]='\0'; line[MAXLINE-2]='\n';  // Truncate long lines
            if (incomment) {
                // We're in a comment.  Keep appending until we reach the end.
                incomment = append_comment_line(line,endcomment);
            } else {
                // We're not in a comment.  Try to parse the line.
                parseline(line,&key,&value);
                // printf("%s = %s\n", key,value);
                if (key==NULL||value==NULL) continue;  // Failed to parse. Skip.
                if (strcmp(key,"comment")==0) {
                    strncpy(endcomment,value,100);
                    incomment = 1;
                    continue;
                }
                // Otherwise, assign this key to a value
                if (assign_key(key,value)) {
                    // Failed to find a match.
                    // For now, we throw an error for this
                    // but one might not want to be this dramatic.
                    fprintf(stderr,"Couldn't parse line: %s\n",line);
                    return 1;
                }
            }
        }
        return 0;
    }
};

