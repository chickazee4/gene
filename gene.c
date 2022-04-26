#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

int aflag=0;
int eflag=0;
int fflag=0;
int Fflag=0;
int iflag=0;
int lflag=0;
int nflag=0; // to be used later. borked right now.
int Rflag=0;
int autofastaid=0;
int usefile=1;

char *legalchars="ACGTacgt";
long start = 0;
long end = 0;
long length = 0;
char **seqs;
char **fastaid;
char *filename;
int indct;
int *indices;
int rowlen=50;

static char bases[256] = {0}; 

int
indexok(int ix){
	for(int i=0; i<indct; i++){
		if(indices[i]==ix)
			return 1;
	}
	return 0;
}

void
loadfromfile(FILE *fp){
	// we need to do a lot more to deal with FASTA input than raw input since there can be multiple sequences.
	if(fflag>0){
		char line[256]; // under essentially no circumstances should a line in a FASTA file exceed 256 bytes so let's avoid wasting mem
		int seqi=-1, ok=0, allow_read=1, seqlen, applen, idlen;
		while(fgets(line, 256, fp)){
			if(line[0]=='>'){
				if(seqi==-1||aflag>0||(nflag>0 && (ok=indexok(seqi))==1)){ // make sure it's reading only sequences it was told to
					fastaid=realloc(fastaid, (seqi+2)*sizeof(char*));
					idlen=strlen(line);
					if(autofastaid==1){
						fastaid[seqi+1]=malloc(idlen-1);
						strncpy(fastaid[seqi+1], line, idlen-1); // idlen-1 necessary bc otherwise it will include trailing \n from the line
					} else {
						fastaid[seqi+1]=fastaid[0];
					}
					seqs=realloc(seqs, (seqi+2)*sizeof(char*));
					seqi++;
					allow_read=1;
				} else {
					allow_read=0;
				}
			} else { // if line is part of the sequence.
				if(seqi>-1 && allow_read==1){
					if(seqs[seqi]!=NULL)
						seqlen=strlen(seqs[seqi]);
					else seqlen=0;
					applen=strlen(line);
					seqs[seqi]=realloc(seqs[seqi], (seqlen+applen)*sizeof(char));
					strcat(seqs[seqi], line);
				}
			}
		}
	} else {
		seqs=malloc(sizeof(char*));
		seqs[0]=(char*)calloc(end - start, sizeof(char));
		fseek(fp, end - start, SEEK_SET);
		fread(seqs[0], sizeof(char), end - start, fp);
	}
	fclose(fp);
}

char *rev(char *inp, int len) // string reverse
{
    char *p1=inp;
    char *p2=inp+len-1;
    while (p1<p2) {
        char c=*p1;
        *p1++=*p2;
        *p2--=c;
    }
    return inp;
}


char
invert(char base){
	return bases[(unsigned char)base]==0 ? base : bases[(unsigned char)base];
}

void
output(char *seq){
	int qlen=strlen(seq);
	if (end==0 || eflag==0){
		end=qlen;
	}
	if(Rflag>0){
		for(int i=0; i<qlen; i++){
			seq[i]=invert(seq[i]);
		}
		seq=rev(seq, qlen);
	}
	for(int oi=0, fi=0; oi<(end-1); oi++){
		if(fflag>0){
			if(strchr(legalchars, seq[oi])!=NULL){
				if(fi>=(start-1) && fi<=end) {
					putc(seq[oi], stdout);
					if(start!=0){
						if(Fflag>0 && (fi - (start-2)) % rowlen == 0)
							putc('\n', stdout);
					} else {
						if(Fflag>0 && (fi+1) % rowlen == 0)
							putc('\n', stdout);
					}
				}
				fi++;
			}
			else end++;
		} else {
			if(strchr(legalchars, seq[oi])!=NULL){
				putc(seq[oi], stdout);
				fi++;
			}
		}
	}
	putc('\n', stdout);
}

int
main(int argc, char **argv)
{
	bases['A']='T';
	bases['C']='G';
	bases['G']='C';
	bases['T']='A';
	bases['a']='t';
	bases['c']='g';
	bases['g']='c';
	bases['t']='a';
	int alen;
	if(argc>1) {
		for(int i=1, j=1; i<argc; i++) {
			if(argv[i][0]=='-') {
				alen=strlen(argv[i]);
				for(int j=1; j<alen; j++){
					switch(argv[i][j]){
						case 'a':
							aflag++;
							break;
						case 'e':
							eflag++;
							end=atol(argv[i+1]);
							i++;
							break;
						case 'f':
							fflag++;
							break;
						case 'F':
							Fflag++;
							if(i<argc-1 && argv[i+1][0]=='>'){
								fastaid=malloc(sizeof(char*));
								fastaid[0]=argv[i+1];
								i++;
							} else {
								autofastaid=1;
							}
							break;
						case 'i':
							iflag++;
							if(i<argc-1)
								filename=argv[i+1];
							else
								fprintf(stderr, "W: No filename specified, ignoring -i flag.");
							i++;
							break;
						case 'l':
							lflag++;
							if(i<argc-1)
								length=atol(argv[i+1]);
							else
								fprintf(stderr, "W: No length specified, ignoring -l flag.");
							i++;
							break;
						case 'm':
							if(i<argc-1){
								switch(argv[i+1][0]){
									case 'd':
										legalchars="ACGTacgt";
										break;
									case 'r':
										legalchars="ACGUacgu";
										bases['A']='U';
										bases['C']='G';
										bases['G']='C';
										bases['U']='A';
										bases['a']='u';
										bases['c']='g';
										bases['g']='c';
										bases['u']='a';
										break;
									case 'n':
										legalchars="ACGTURYKMSWBDHVNacgturykmswbdhvn";
										bases['A']='T';
										bases['C']='G';
										bases['G']='C';
										bases['T']='A';
										bases['U']='A';
										bases['a']='t';
										bases['c']='g';
										bases['g']='c';
										bases['t']='a';
										bases['u']='a';
										break;
									case 'a':
										legalchars="ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv";
										break;
									case 'x':
										legalchars="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
										break;
									default:
										fprintf(stderr, "W: Invalid parse mode specified. Will default to treating input as a DNA sequence.\n");
										break;
								}
								i++;
							} else fprintf(stderr, "W: Invalid parse mode specified. Will default to treating input as a DNA sequence.\n");
							break;
						case 'r':
							if(i<argc-1){
								rowlen=atoi(argv[i+1]);
								i++;
							} else fprintf(stderr, "W: No row length specified. Will default to 50.");
							break;
						case 'R':
							Rflag++;
							break;
						case 's':
							start=atol(argv[i+1]);
							i++;
							break;
						default:
							break;
					}
				}
			}
			else {
				seqs=malloc(sizeof(char*));
				seqs[0]=argv[i];
				usefile=0;
			}
		}
	}
	if(lflag>0){
		if(start!=0)
			end=start+length;
		else if(end!=0)
			start=end-length;
		else end=length;
	}
	if(usefile>0){
		if(iflag>0){
			FILE *fp;
			if((fp=fopen(filename, "r"))==NULL){
				fprintf(stderr, "E: Input file not found.");
				return 1;
			}
			loadfromfile(fp);
		} else loadfromfile(stdin);
	}
	int si=0;
	while(seqs[si]!=NULL){
		if(Fflag>0){
			if(fastaid==NULL || fastaid[si]==NULL){
				fprintf(stderr, "E: Need to specify an ID or use a FASTA formatted file if outputting as FASTA.\n");
				exit(1);
			} else {
				printf("%s\n", fastaid[si]);
			}
		}
		output(seqs[si]);
		si++;
	}
}
