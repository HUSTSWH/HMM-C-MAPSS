#include <cstdio>
#include <cstring>
#include <algorithm>

using namespace std;

#define prev oasdfasdfpoqert
#define index aserhpasdjifpaohi

double dp[400][5];
int prev[400][5];
int machnum[100000];
int index[100000];
double pool[100000][11];
int main ()
{
    freopen("train_data.txt","r",stdin);
    char s[300];
    scanf("%[^\n]s",s);
    int n=0;
    while(scanf("%d %d",machnum+n,index+n)!=EOF){
	for(int j=0; j<11; j++)
	    scanf("%lf",&pool[n][j]);
	n++;
    }
    double SX, SX2;
    int l[60],r[60];
    int pos=0;
    for(int i=0; i<60; i++) {
	l[i]=pos+8;
	while(pos<n && machnum[pos]==machnum[pos+1])pos++;
	r[i]=++pos;
    }
    for(int t=0; t<60; t++){
        SX = SX2 = 0.0;
	for(int i=0; i<r[t]-l[t]; i++){
            for(int k=0; k<11; k++){
		SX += pool[i+l[t]][k];
		SX2 += pool[i+l[t]][k]*pool[i+l[t]][k];
            }
            dp[i][1] = SX2/(i+1)/11-SX*SX/(i+1)/(i+1)/121;
	    dp[i][2] = dp[i][3] = dp[i][4] = 3e100;
	}
	for(int i=4; i<r[t]-l[t]; i++){
	    SX = SX2 = 0;
            for(int j=0; j<i; j++){
		for(int k=0; k<11; k++){
		    SX += pool[i-j+l[t]][k];
	            SX2 += pool[i-j+l[t]][k]*pool[i-j+l[t]][k];
		}
		double DX = SX2/(j+1)/11 - SX/(j+1)/11*SX/(j+1)/11;
		for(int k=2; k<5; k++){
	            if(dp[i][k] > dp[i-j-1][k-1] + DX){
			dp[i][k] = dp[i-j-1][k-1] + DX;
			prev[i][k] = i-j-1;
	            }
		}
	    }
	}
	int d = r[t]-l[t]-1;
	int c = prev[d][4];
	int b = prev[c][3];
	int a = prev[b][2];
	printf("%d %d %d %d\n",a+9,b+9,c+9,d+9);
    }
    return 0;
}

