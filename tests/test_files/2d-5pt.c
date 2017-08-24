double a[M][N];
double b[M][N];
double c0;
double c1;

for(int j=1; j < M-1; j++){
for(int i=1; i < N-1; i++){
b[j][i] = c0*a[j][i]
+ c1 * (a[j][i-1] + a[j][i+1])
+ c1 * (a[j-1][i] + a[j+1][i])
;
}
}
