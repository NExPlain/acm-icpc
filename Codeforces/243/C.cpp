#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<vector>
#include<algorithm>
using namespace std;

typedef long long LL;
#define CLR(a,b) memset(a,b,sizeof(a))

const int N = 1e5 + 20;
int n,m,s,e;
vector<int> vec[N];
int a[N],b[N];
int f[N][302];
int ans;
void check(int k,int x)
{
    ans = max(ans,min(k, (s-x)/e));
}
void solve()
{
    ans = 0;
    CLR(f, -1);
    f[0][0] = 0;
    for(int i = 1 ;i <= m ; i ++){
        int x;
        cin >> x;
        for(int j = 1 ; j <= 300 ; j ++){
            if(f[i-1][j-1] != -1){
                int t = lower_bound(vec[x].begin(),vec[x].end(),f[i-1][j-1] + 1) - vec[x].begin();
                if(t != vec[x].size()) f[i][j] = vec[x][t];
            }
        }
        for(int j = 0 ;j <= 300; j ++){
            if(f[i][j] == -1) f[i][j] = f[i-1][j];
            else if(f[i-1][j] != -1)f[i][j] = min(f[i][j], f[i-1][j]);
        }
        for(int j = 0 ; j <= 300; j ++){
            if(f[i][j] != -1) check(j, i + f[i][j]);
        }
    }
    printf("%d\n", ans);
}
int main()
{
    cin >> n >> m >> s >> e;
    for(int i = 1 ; i <= n ; i ++){
        cin >> a[i];
        vec[a[i]].push_back(i);
    }
    solve();
    return 0;
}
