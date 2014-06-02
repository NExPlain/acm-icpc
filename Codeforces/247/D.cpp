#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<set>
#include<vector>
#include<algorithm>
using namespace std;

typedef unsigned long long LL;
typedef pair<LL,int> P;

const int N = 60;
LL n;
int k;
LL C[N][N];
set<P> vec[N];
set<P>::iterator ite;
bool vis[N];

void solve()
{
    int kk = k-1;
    LL ans = 0;
    if(n == 0 && k != 1){
        ans = 1;
    }
    while(n){
        ite = vec[kk].upper_bound(make_pair(n,100));
        --ite;
        ite = vec[kk].lower_bound(make_pair((*ite).first, 0));
        vis[(*ite).second] = 1;
        n -= (*ite).first;
        vec[kk].erase(ite);
        kk--;
    }
    for(int i = 0 ;i < N ; i++){
        if(vis[i]){
            ans += ((LL)1 << i);
        }
    }
    cout << ans << endl;
}
int main()
{
    memset(C,0,sizeof(C));
    C[0][0] = 1;
    vec[0].insert(make_pair(C[0][0],0));
    for(int i = 1 ; i < N ; i ++){
        C[i][0] = 1;
        vec[0].insert(make_pair(C[i][0],i));
        for(int j = 1 ;j <= i ; j ++){
            C[i][j] = C[i][j-1] * (i-j+1) / j;
            vec[j].insert(make_pair(C[i][j],i));
        }
    }
    cin >> n >> k;
    solve();
    return 0;
}
