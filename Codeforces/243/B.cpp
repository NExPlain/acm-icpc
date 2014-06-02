#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<string>
#include<algorithm>
#include<vector>
#include<map>
#include<set>
#include<queue>
#include<cmath>
using namespace std;

typedef long long LL;
#define CLR(a,b) memset(a,b,sizeof(a))
#define INF 0x3f3f3f3f
#define FI first
#define SE second
#define PB push_back
#define MP make_pair

const int N = 100+20;
int n,m,K;
int a[N][N];

void solve()
{
    int ans = INF;
    for(int sta = 0 ; sta < (1 << n) ; sta ++){
        int lef = K;
        bool ok = 1;
        for(int i = 0 ; i < m ; i ++){
            int same = n;
            for(int j = 0 ; j < n ; j ++){
                if((a[j][i] ^ ((sta >> j)&1)) == 1)same--;
            }
            int use = 0;
            if(same > n - same)use = n-same;
            else use = same;
            if(use <= lef)lef -= use;
            else ok = 0;
        }
        if(ok){
            ans = min(ans,K-lef);
        }
    }
    if(ans != INF)printf("%d\n",ans);
    else printf("-1\n");
}
int main()
{
    scanf("%d%d%d",&n,&m,&K);
    for(int i = 0; i < n ; i ++){
        for(int j = 0 ; j < m ; j ++){
            scanf("%d",&a[i][j]);
        }
    }
    if(n <= K){
        solve();
        return 0;
    }
    int ans = INF;
    for(int i = 0 ; i < n; i ++){
        int lef = K;
        bool ok = 1;
        for(int k = 0 ; k < n ; k ++){
            if(k == i) continue;
            int same = 0;
            for(int j = 0 ; j < m ; j ++){
                if(a[i][j] == a[k][j]){
                    same ++;
                }
            }
            int use = 0;
            if(same > m-same){
                use = m-same;
            }else{
                use = same;
            }
            if(use <= lef){
                lef -= use;
            }else{
                ok = 0;
            }
        }
        if(ok)ans = min(ans,K-lef);
    }
    if(ans != INF)
        printf("%d\n",ans);
    else
        printf("-1\n");
    return 0;
}

/*
3 7 10
1 0 0 1 1 0 1
0 0 1 1 1 1 1
1 1 1 1 1 1 1
*/
