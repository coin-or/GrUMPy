#include<stdio.h>
#include<stdlib.h>

int main() {
  FILE * pFile;
  int iSize;
  char * buffer;
  size_t result;

  pFile = fopen("p0201_GLPK.in", "r");
  if (pFile==NULL) {
    fputs("File Error", stderr);
    exit(1);
  }

  // obtain file size
  fseek(pFile, 0, SEEK_END);
  iSize = ftell(pFile);
  rewind(pFile);

  // allocate memory to contain the file
  buffer = (char*) malloc (sizeof(char)*iSize);
  if (buffer==NULL) {
    fputs("Memory Error", stderr);
    exit(2);
  }

  // copy the file into the buffer
  result = fread(buffer,1,iSize,pFile);
  if (result!=iSize) {
    fputs("Reading Error", stderr);
    exit(3);
  }

  fclose(pFile);

  int length = 0;
  int i = 0;
  int j = 0;
  int k;
  char temp[1000];
  for(i=0;i<iSize;i++) {
    if (buffer[i] !='\n') {
      temp[j]=buffer[i];
      j++;
      continue;
    }
    temp[j++]='\n';
    printf("%.*s", j, temp);
    j=0;
  }
  return 0;
}
