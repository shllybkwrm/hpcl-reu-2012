/*
	Shelly Bagchi, HCPL REU, Fall 2012
*/



#include <stdio.h>


int sx,tx,fy,fx,a,g,m,mism,optimal;

void print1dint(int n, int arr[n]) {
	int i;

	//printf("1D int array size: n = %d\n",n);
	for(i=0; i<n; i++) {

		printf("%d ",arr[i]);

	}

	printf("\n\n");

	return;
}



void print2dint(int n, int m, int arr2d[n][m]) {
	int i,j;

	printf("2D int array size: n = %d, m = %d\n",n,m);

	for(i=0; i<n; i++) {

		for(j=0; j<m; j++) {

			printf("%d\t",arr2d[i][j]);

		}

		printf("\n");
	}

	printf("\n");


	return;
}



void print3dint(int n, int m, int l, int arr3d[n][m][l]) {
	int i,j,k;

	printf("3D int array size: n = %d, m = %d\n",n,m);

	for(i=0; i<n; i++) {

		for(j=0; j<m; j++) {

			printf("{");

			for(k=0; k<l; k++) {

				printf("%d",arr3d[i][j][k]);
				if( k<(l-1) ) printf(",");

			}

			printf("}\t");

		}

		printf("\n");
	}

	printf("\n");


	return;
}




int max(int i, int j, int paths[fy][fx][3], int left, int diag, int up) {

	if(left >= diag  &&  left >= up) paths[i][j][0] = 1;
	if(diag >= left  &&  diag >= up) paths[i][j][1] = 1;
	if(up >= left  &&  up >= diag) paths[i][j][2] = 1;


	if(paths[i][j][0]) return left;
	else if(paths[i][j][1]) return diag;
	else if(paths[i][j][2]) return up;
	
	return 999;
}



void sim(int i, int j, char s[sx], char t[tx], int f[fy][fx], int paths[fy][fx][3]) {
	
	if(i==0 && j==0) {
		f[0][0] = 0;

		for(i=0; i<fy; i++) {
			for(j=0; j<fx; j++) {
				if(i==0 && j==0);
				else sim(i,j,s,t,f,paths);
			}
			
		}
		
		optimal = f[i-1][j-1];

	}

	else if(i==0) {
		f[i][j] = j*g;
		paths[i][j][0] = 1;
	}

	else if(j==0) {
		f[i][j] = i*g;
		paths[i][j][2] = 1;
	}

	else {
		f[i][j] = max(i,j,paths, (f[i][j-1]+g), ( f[i-1][j-1]+( s[j-1]==t[i-1] ? m : mism ) ), (f[i-1][j]+g) );
	}



	return;
}



void insertGap(int index, int n, char arr[n]) {
	int i;
	char temp, temp1;
	
	if(index<0) index=0;
	
	for(i=0; i<n; i++) {
	
		if(i==index) {
			temp = arr[i];
			arr[i] = '-';
		}
		
		else if(i>index) {
			temp1 = arr[i];
			arr[i] = temp;
			temp = temp1;
		}

	}


	return;
}



void align(char sAligned[a], char tAligned[a], int paths[fy][fx][3]) {
	int i=(fy-1),j=(fx-1);

	
	while( i>=0 && j>=0 ) {
		//printf("i = %d\n", i);
		//printf("j = %d\n", j);

		if(paths[i][j][0]) {  //left => insert in t
			insertGap( (j-1),a,tAligned);
			j--;
		}
		else if(paths[i][j][1]) {  //diag => match
			i--;
			j--;
		}
		else if(paths[i][j][2]) {  //up => insert in s
			insertGap( (i-1),a,sAligned);
			i--;
		}
		else break;


		/*printf("partially aligned sequence s:\n");
		puts(sAligned);
		printf("partially aligned sequence t:\n");
		puts(tAligned);*/
		
	}


	return;
}



int score(char sAligned[a], char tAligned[a]) {
	int i, sigma=0;
	
	for(i=0; i<a; i++) {
		
		if(sAligned[i]=='\0' || tAligned[i]=='\0') break;
		else if(sAligned[i] == tAligned[i]) sigma += m;			// match
		else if(sAligned[i]=='-' || tAligned[i]=='-') sigma += g;	// gap
		else  sigma += mism;						// mismatch
		
	}
	
	return sigma;
}



int main() {
	int i,j,k;
	g = -1;		// gap penalty
	m = 2;		// match
	mism = -1;	// mismatch penalty

	char s[] = "ACACACTA";
	char t[] = "AGCACACA";

	sx = sizeof(s)/sizeof(char) -1;
	tx = sizeof(t)/sizeof(char) -1;

	fy = tx+1;
	fx = sx+1;

	int f[fy][fx];  // similarity matrix
	for(i=0; i<fy; i++) {
		for(j=0; j<fx; j++) {
			f[i][j] = -1;
		}
	}


	int paths[fy][fx][3];

	for(i=0; i<fy; i++) {
		for(j=0; j<fx; j++) {
			for(k=0; k<3; k++) {
				paths[i][j][k] = 0;
			}
		}
	}


	a = sx + tx;

	char sAligned[a];
	char tAligned[a];

	for(i=0; i<a; i++) {
		if(i<sx) sAligned[i] = s[i];
		else sAligned[i] = '\0';
	}
	for(j=0; j<a; j++) {
		if(j<tx) tAligned[j] = t[j];
		else tAligned[j] = '\0';
	}


	printf("Sequence s:\n");
	puts(s);
	printf("Sequence t:\n");
	puts(t);
	printf("Gap Penalty: %d\n\n",g);


	sim(0,0,s,t,f,paths);
	printf("Similarity matrix computed.\n");
	printf("The optimal score is %d.\n\n", optimal);

	print2dint(fy,fx,f);
	print3dint(fy,fx,3,paths);


	align(sAligned,tAligned,paths);
	printf("Optimal global alignment computed.\n\n");

	printf("Aligned sequence s:\n");
	puts(sAligned);
	printf("Aligned sequence t:\n");
	puts(tAligned);
	
	
	printf("The alignment score is:  %d\n", score(sAligned, tAligned) );
	puts("\n\n");

	return 0;
}
