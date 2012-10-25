/*
	Muhammad Abdul-Rahim, Shelly Bagchi, Adam McCormack, Johnathan Ross
	GWU HCPL REU, Fall 2012
	Sequence Alignment using Needleman-Wunsch & Smith-Waterman Algorithms
	Non-Parallel C Implementation
*/



#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>



unsigned long int sx, tx, fy, fx, a;
int g, m, mism, optimal;
unsigned int algorithm;

/*
	Global Variables:
	
	sx		size of sequence 1 (s)
	tx		size of sequence 2 (t)
	fy		y-size of similarity matrix f
	fx		x-size of similarity matrix f
	a		maximum size of aligned sequences (sx + tx)
	g		gap penalty
	m		match bonus
	mism		mismatch penalty
	optimal		score of optimal alignment, i.e. maximum possible score
	algorithm	Needleman-Wunsch = 0, Smith-Waterman = 1
*/




inline void
print1dint(const unsigned long int n, int arr[n])
{
	unsigned long int i;

	for(i=0; i<n; i++)
		printf("%d\t",arr[i]);

	putchar('\n');
}




inline void
print2dint(const unsigned long int n, const unsigned long int m, int arr2d[n][m])
{
	unsigned long int i,j;

	printf("2D int array size: y-size = %lu, x-size = %lu\n",n,m);

	for(i=0; i<n; i++)
	{
		for(j=0; j<m; j++)
			printf("%d\t",arr2d[i][j]);
		putchar('\n');
	}

	putchar('\n');
}




inline void
print3dint(const long unsigned int n, const long unsigned int m, const long unsigned int l, int arr3d[n][m][l])
{
	unsigned long int i,j,k;

	printf("3D int array size: y-size = %lu, x-size = %lu\n",n,m);

	for(i=0; i<n; i++)
	{
		for(j=0; j<m; j++)
		{
			putchar('{');
			for(k=0; k<l; k++)
			{
				printf("%d",arr3d[i][j][k]);
				if( k<(l-1) )
					putchar(',');
			}
			printf("}\t");
		}
		putchar('\n');
	}

	putchar('\n');
}




int
max(const unsigned long int i, const unsigned long int j, int paths[fy][fx][3],
	const int left, const int diag, const int up)
{
	if(diag >= left && diag >= up)	 paths[i][j][1] = 1;
	if(left >= diag && left >= up) 	 paths[i][j][0] = 1;
	if(up	>= left && up	>= diag) paths[i][j][2] = 1;

	if(paths[i][j][1]) 	return diag;
	else if(paths[i][j][0]) return left;
	else if(paths[i][j][2]) return up;
	
	return 0;
}




/*
	Compute 2D similarity matrix
*/


void
findSim(unsigned long int i, unsigned long int j, char s[sx], char t[tx], int f[fy][fx], int paths[fy][fx][3])
{
	if(i==0 && j==0)
	{
		f[0][0] = 0;
		for(i=0; i<fy; i++)
		{
			for(j=0; j<fx; j++)
			{
				if(i==0 && j==0);
				else findSim(i,j,s,t,f,paths);
			}
			
		}
		
		optimal = f[i-1][j-1];
	}

	else if(i==0)
	{
		f[i][j] = algorithm?0:j*g;
		paths[i][j][0] = 1;
	}

	else if(j==0)
	{
		f[i][j] = algorithm?0:i*g;
		paths[i][j][2] = 1;
	}

	else
	{
		f[i][j] = max(i,j,paths, (f[i][j-1]+g), ( f[i-1][j-1]+( s[j-1]==t[i-1] ? m : mism ) ), (f[i-1][j]+g) );
		
		if(algorithm && f[i][j]<0)
			f[i][j] = 0;
	}
}







void bumpRow(int f1d[2][fx])
{
	unsigned long int i;
	
	print1dint(fx,f1d[0]);
	
	for(i=0; i<fx; i++)
	{
		f1d[0][i] = f1d[1][i];
		f1d[1][i] = 0;
	}
}




/*
	Find similarity matrix, but only one row at a time (space optimization)
*/


void
findSim1d(unsigned long int i, unsigned long int j, char s[sx], char t[tx], int f1d[2][fx], int pathsFrom1d[fy][fx][3])
{
	if(i==0 && j==0)
	{
		f1d[0][0] = 0;
		for(i=0; i<fy; i++)
		{
			for(j=0; j<fx; j++)
			{
				if(i==0 && j==0);
				else findSim1d(i,j,s,t,f1d,pathsFrom1d);
			}
			
			if(i!=0) bumpRow(f1d);
			
		}
		
		print1dint(fx,f1d[0]);
		optimal = f1d[1][j-1];
	}

	else if(i==0)
	{
		f1d[0][j] = algorithm?0:j*g;
		pathsFrom1d[i][j][0] = 1;
	}

	else if(j==0)
	{
		f1d[1][j] = algorithm?0:i*g;
		pathsFrom1d[i][j][2] = 1;
	}

	else
	{
		f1d[1][j] = max(i,j,pathsFrom1d, (f1d[1][j-1]+g), ( f1d[0][j-1]+( s[j-1]==t[i-1] ? m : mism ) ), (f1d[0][j]+g) );
		
		if(algorithm && f1d[1][j]<0)
			f1d[1][j] = 0;
	}
}




/*
	Used in actual sequence alignment
*/


void
insertGap(int index, const unsigned long int n, char arr[n])
{
	unsigned long int i;
	char temp = '\0', temp1 = '\0';
	
	if(index<0) index=0;
	
	for(i=0; i<n; i++)
	{
		if(i==(unsigned long)index)
		{
			temp = arr[i];
			arr[i] = '-';
		}
		
		else if(i>(unsigned long)index)
		{
			temp1 = arr[i];
			arr[i] = temp;
			temp = temp1;
		}
	}
}



/*
	Aligns actual sequences
*/

void
align(char sAligned[a], char tAligned[a], int paths[fy][fx][3])
{
	long int i=(fy-1),j=(fx-1);
	
	while( i>=0 && j>=0 )
	{
		if(paths[i][j][1])		//diag => match
		{
			i--;
			j--;
		}
		else if(paths[i][j][0])		//left => insert in t
		{
			insertGap( i,a,tAligned);
			j--;
		}
		else if(paths[i][j][2])		//up => insert in s
		{
			insertGap( j,a,sAligned);
			i--;
		}
		else break;
		
	}

}




int
score(char sAligned[a], char tAligned[a])
{
	unsigned long int i, sum=0;
	
	for(i=0; i<a; i++)
	{	
		if(sAligned[i]=='\0' || tAligned[i]=='\0') 	break;
		else if(sAligned[i] == tAligned[i]) 		sum += m;	// match
		else if(sAligned[i]=='-' || tAligned[i]=='-') 	sum += g;	// gap
		else  						sum+=mism;	// mismatch
		
	}
	
	return sum;
}



// Main Function

int
main(int argc, char* argv[])
{
	unsigned long int i;
	time_t start, end;
	
	if(argc==3)
	{
		// defaults
				
		g = -1;			// gap penalty
		m = 2;			// match bonus
		mism = -1;		// mismatch penalty
		algorithm = 0;		// Needleman-Wunsch
	}
	
	else if(argc==7)
	{
		g = atoi(argv[3]);
		m = atoi(argv[4]);
		mism = atoi(argv[5]);
		algorithm = atoi(argv[6]);
	}
	
	else
	{
		printf("\nIncorrect usage.\nSyntax: %s Sequence1(s) Sequence2(t)\n\nOR\n\n%s Sequence1(s) Sequence2(t) gapPenalty(default=-1) matchBonus(default=2) mismatchPenalty(default=-1) whichAlgorithm(default NW=0;SW=1)\n\n", argv[0], argv[0]);
		return -1;
	}
	

	sx = strlen(argv[1]);
	tx = strlen(argv[2]);
	
	char s[sx], t[tx];
	
	strcpy(s,argv[1]);
	strcpy(t,argv[2]);

	fy = tx+1;
	fx = sx+1;

	int f[fy][fx];  // similarity matrix
	memset(f,0,sizeof(f));

	int paths[fy][fx][3];
	memset(paths,0,sizeof(paths));
	
	
	// Variables for 1D implementation test

	int f1d[2][fx];  // similarity matrix (1D)
	memset(f1d,0,sizeof(f1d));

	int pathsFrom1d[fy][fx][3];
	memset(pathsFrom1d,0,sizeof(pathsFrom1d));
	
	// end
	
	
	a = sx + tx;

	char sAligned[a];
	char tAligned[a];

	for(i=0; i<a; i++)
	{
		if(i<sx) sAligned[i] = s[i];
		else sAligned[i] = '\0';
		if(i<tx) tAligned[i] = t[i];
		else tAligned[i] = '\0';
	}
	
	// End initializations
	

	// Print current config
	
	printf("\nSequence s:\n\t");
	puts(s);
	printf("Sequence t:\n\t");
	puts(t);
	printf("Gap Penalty: %d\n",g);
	printf("Match Bonus: %d\n",m);
	printf("   Mismatch: %d\n",mism);
	printf("      Using: ");
	if( algorithm==0)
		puts("Needleman-Wunsch Algorithm");
	else
		puts("Smith-Waterman");
	puts("\n");
	
	
	
	
	// Start timer
	start = time(NULL);
	
	
	
	// Start computations
	
	printf("Computing similarity...\n");
	findSim(0,0,s,t,f,paths);
	printf("Similarity matrix computed.\n");
	printf("The optimal score is %d.\n\n", optimal);

	print2dint(fy,fx,f);
	print3dint(fy,fx,3,paths);
	
	
	// 1D implementation test
	
	printf("Starting 1D computation...\n");
	findSim1d(0,0,s,t,f1d,pathsFrom1d);
	putchar('\n');

	print3dint(fy,fx,3,pathsFrom1d);
	printf("1D computation finished.\n");
	printf("The optimal score is %d.\n\n", optimal);
	
	// end test
	

	printf("Computing optimal alignment...\n");
	align(sAligned,tAligned,pathsFrom1d);
	printf("Optimal global alignment computed.\n\n");

	printf("Aligned sequence s:\n");
	puts(sAligned);
	printf("Aligned sequence t:\n");
	puts(tAligned);
	
	printf("The alignment score is:  %d\n", score(sAligned, tAligned) );
	
	// End computations
	
	
	
	
	// End timer
	end = time(NULL);
	printf("Time taken for computation: %f sec\n", difftime(end,start));


	puts("\n");
	
	
	
	return 0;
}
