// PrTool.cpp
//

/** PrTool.cpp PDB2fasta
*
*		@Name: 		a small pdb/fasta tool
*		@Author:	Li Zhixin
*		@Version: 	1.0
*		@Tool: 		gcc -std=c99
*		@Start: 	----/04/24
*		@Done:		----/05/09
*/
#include<stdio.h>
#include "stdafx.h"
#include<string.h>
#include<stdlib.h>
#define STANDARD_LINE_WIDTH (80)
#define READ_BUFF (160)
const char THREE_LIST[21][4] = { "???","ALA","CYS","ASP","GLU","PHE" ,"GLY","HIS",
"ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR" };

const char ONE_LIST[21] = { '?','A','C','D','E','F','G','H','I','K','L','M',
'N','P','Q','R','S','T','V','W','Y' };

const float MO_WEIGHT[21] = { -1,89.079,121.145,133.089,147.116,165.177,75.052,
155.141,131.16,146.17,131.16,149.199,132.104,115.117,146.131,174.188,105.078,
119.105,117.133,204.213,181.176 };

int showHelp()// show cmd help
{
	//printf("    prTool 1.0 Li Zhixin\n");//normal disply
	system("echo \"\033[44;37;5m prTool 1.0 Li Zhixin \033[0m\"");//color disply
	printf("prTool <option> <input file> <output file> \n");

	printf("<option>:\n");
	printf("-A: read sequence from ATOM. (input: PDB, output: FASTA)\n");
	printf("-S: read sequence from SEQRES. (input: PDB, output: FASTA)\n");
	printf("-W: read molecular weight from FASTA. (input: FASTA, output: TXT)\n");
	printf("-h: show this help.\n");
	printf("-a: show about.\n");
	return 0;
}
char getShortAmino(const char three[4])// three to one amino abbr.
{
	for (int i = 1; i <= 20; i++)
	{
		if (strcmp(three, THREE_LIST[i]) == 0)
			return ONE_LIST[i];
	}
	return ONE_LIST[0];
}
float getWeight(const char a)// get one weight
{
	for (int i = 1; i <= 20; i++)
	{
		if (a == ONE_LIST[i])
			return MO_WEIGHT[i];
	}
	return MO_WEIGHT[0];
}
/* return 0: not atom line, 1: new amino, 2:success, 3:ter ( better to use enum)*/
int splitAtom(char*line, int NUM, char name[4], int *num)
{
	const char*splitter = " ";//divided by space
	char*p = strtok(line, splitter);// in lib string, token
	int paranum = 0;
	if (strcmp("ATOM", p) == 0)//find line atom
	{
		while (p != NULL)
		{
			p = strtok(NULL, splitter);
			paranum++;
			if (paranum == 3)strcpy(name, p);
			else if (paranum == 5)
			{
				*num = atoi(p);
				break;
			}
		}
		if (*num > NUM)return 1;
		else return 2;
	}
	else if (strcmp("TER", p) == 0)
	{
		return 3;
	}
	return 0;
}
/*CHAIN:A,B,...  nlCnt:line cnt, to fit 80 width
return 2:new chain 0:success*/
int splitSeq(char*line, char CHAIN, FILE* outputFile, int*nlCnt)
{
	const char*splitter = " ";
	char name[4];
	char*p = strtok(line, splitter);
	int paranum = 0;
	if (strcmp("SEQRES", p) == 0)
	{
		while (p != NULL)
		{
			p = strtok(NULL, splitter);
			if (p == NULL || p[0] == '\n')break;
			paranum++;
			if (paranum == 2)
			{
				if (p[0] != CHAIN)return 2;
			}
			else if (paranum >= 4)
			{
				strcpy(name, p);
				if (*nlCnt == STANDARD_LINE_WIDTH)
				{
					fprintf(outputFile, "\n");
					*nlCnt = 0;
				}
				fprintf(outputFile, "%c", getShortAmino(name));
				*nlCnt += 1;
			}
		}

	}

	return 0;
}
int atom(FILE*fp, FILE*outputFile)
{
	char line[READ_BUFF];
	char starLetter = 'A';
	char name[10];
	int num = 0;
	int NUM = -1;
	int newChain = 0;//flag to show new chain
	int newlineCnt = 0;
	while (!feof(fp))
	{
		fgets(line, READ_BUFF, fp);  //read one line
		int rt = splitAtom(line, NUM, name, &num);
		if (rt == 1)
		{
			if (newChain == 0)
			{
				fprintf(outputFile, ">AA00:%c|PDBID|CHAIN|SEQUENCE\n", starLetter);
				newChain = 1;
			}

			if (newlineCnt == STANDARD_LINE_WIDTH)
			{
				fprintf(outputFile, "\n");
				newlineCnt = 0;
			}
			fprintf(outputFile, "%c", getShortAmino(name)); //output
			newlineCnt++;
			NUM = num;
		}
		else if (rt == 3)
		{
			fprintf(outputFile, "\n");
			NUM = -1;
			starLetter++;
			newChain = 0;
			newlineCnt = 0;
		}
	}
	//printf("\n");
	return 0;
}
int seqres(FILE*fp, FILE*outputFile)//seq mode, call splitSeq
{
	char line[READ_BUFF];
	int newlineCnt = 0;
	char starLetter = 'A';
	char name[4];
	int num = 0;
	int NUM = -1;
	int newChain = 0;//flag to show new chain
	fprintf(outputFile, ">AA00:%c|PDBID|CHAIN|SEQUENCE\n", starLetter);
	while (!feof(fp))
	{
		fgets(line, READ_BUFF, fp);  //read one line
		int rt = splitSeq(line, starLetter, outputFile, &newlineCnt);
		if (rt == 2)
		{
			starLetter++;
			newlineCnt = 0;
			fprintf(outputFile, "\n>AA00:%c|PDBID|CHAIN|SEQUENCE\n", starLetter);
			fseek(fp, - STANDARD_LINE_WIDTH, SEEK_CUR);//back to last line

		}
	}
	return 0;
}
int weight(FILE*fp, FILE*outputFile)// weight mode
{
	char line[READ_BUFF];
	char chainCnt = 'A';
	int length = 0;
	double wt = .0;
	while (!feof(fp))
	{
		if (fgets(line, READ_BUFF, fp) == NULL)break;  //read one line
		if (line[0] == '>')
		{
			if (length > 0)
			{
				fprintf(outputFile, "\tAmino number:\t\t%d\n\tMolecular weight:\t%.3f\n", length, wt - 18 * (length - 1));
				length = 0;
				wt = 0;
			}
			fprintf(outputFile, "====================\nChain-%c:\n", chainCnt);
			chainCnt++;

		}
		/*else if (feof(fp))
		{
		fprintf(outputFile, "\tAmino number:\t\t%d\n\tMolecular weight:\t%.3f\n", length, wt);
		break;
		}*/
		else
		{
			for (int i = 0; line[i] != '\n'; i++)
			{
				float rt = getWeight(line[i]);
				if (rt > 0)
				{
					length++;
					wt += rt;
				}
				else break;
			}
		}
	} // while end 
	if (length > 0)
		fprintf(outputFile, "\tAmino number:\t\t%d\n\tMolecular weight:\t%.3f\n", length, wt - 18 * (length - 1));

	return 0;
}
int readArg(int argc, char* argv[])// cmd
{
	if (argc == 1)
	{
		system("echo \"PrTool v1.0 [Written by \033[1;32;5m Li Zhixin \033[0m\"]");//cool~
																				   //printf("Email: waynefn@live.com\n");
	}
	else if (argc>1 && argv[1][0] == '-')
	{
		if (argc == 2)
		{
			if (argv[1][1] == 'h')
				showHelp();
			else if (argv[1][1] == 'a')
			{

				system("echo \"Written by \033[1;32;5m Li Zhixin \033[0m\"");//cool~
				printf("Email: waynefn@live.com\n");
			}
			else
				printf("-%c is undefined now...\n", argv[1][1]);
		}
		else if (argc == 3)
		{
			printf("Commend error, do you forget output file?\n");
		}
		else if (argc == 4)
		{
			FILE*fp1 = fopen(argv[2], "rt");
			if (fp1 == NULL)
			{
				printf("Open input file error.\n");
				return -1;
			}
			FILE*fp2 = fopen(argv[3], "wt");
			if (fp2 == NULL)
			{
				printf("Open output file error.\n");
				return -2;
			}

			if (argv[1][1] == 'A')
			{
				printf("Input PDB file: %s (mode=ATOM)\n", argv[2]);
				atom(fp1, fp2);
			}
			else if (argv[1][1] == 'S')
			{
				printf("Input PDB file: %s (mode=SEQRES)\n", argv[2]);
				seqres(fp1, fp2);
			}
			else if (argv[1][1] == 'W')
			{
				printf("Input FASTA file: %s (mode=molecular weight)\n", argv[2]);
				weight(fp1, fp2);
			}
			else
			{
				printf("Unknown commend -%c, type -h for help.\n", argv[1][1]);
				return -1;
			}
			fclose(fp1);
			fclose(fp2);
			printf("Output file: %s\n", argv[3]);
		}
		else
		{
			printf("Unknown commend, type -h for help.\n");
		}

	}
	else
	{
		printf("Input error, type -h for help.\n");
	}
	return 0;
}

int main(int argc, char* argv[])
{
	readArg(argc, argv);

	return 0;
}

