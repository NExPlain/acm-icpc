#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<algorithm>
#include<cstring>
#include<string>
#include<vector>
#include<cmath>
#include<set>
#include<stack>
#include<map>
using namespace std;

typedef pair<int,int> PII;
typedef long long LL;
#define CLR(a,b) memset(a,b,sizeof(a))
#define PB push_back
#define MP make_pair
#define INF 0x3f3f3f3f
const int N = 100000+20;

int lowbit(int x){return x & (-x);}
/*
void add (int idx , int v) {
    for (int i = idx ; i <= n ; i += lowbit (i))
        s[i] += v;
}
int sum (int idx) {
    int ret = 0;
    for (int i = idx ; i > 0 ; i -= lowbit (i))
        ret += s[i];
    return ret;
}
 */

int n,m;
set<pair<LL,int> > g[N];
LL v[N];
PII sv[N];
bool vis[N];

int uset[N];
LL miset[N];
int cnt[N];
int ranks[N];
void makeSet(int size) {
    for(int i = 0;i < size;i++)  uset[i] = i;
    for(int i = 0;i < size;i++)  {miset[i] = v[i];cnt[i] = 1;}
    for(int i = 0;i < size;i++)  ranks[i] = 0;
}
int find(int x) {
    if (x != uset[x]) uset[x] = find(uset[x]);
    return uset[x];
}
void unionSet(int x, int y,LL &ans) {
    if ((x = find(x)) == (y = find(y))) return;
    ans += miset[x] * 2 * cnt[x] * cnt[y];
    miset[x] = miset[y] = min(miset[x],miset[y]);
    cnt[x] = cnt[y] = cnt[x] + cnt[y];
    if (ranks[x] > ranks[y]) uset[y] = x;
    else {
        uset[x] = y;
        if (ranks[x] == ranks[y]) ranks[y]++;
    }
}
set<pair<LL,int> >::iterator ite;
void solve()
{
    LL ans = 0;
    makeSet(n);
    for(int i = n - 1 ; i >= 0 ;i --){
        int x = sv[i].second;
        ite = g[x].upper_bound(MP(v[x],x));
        for(;ite!=g[x].end();ite++){
            int nxt = (*ite).second;
            unionSet(x,nxt,ans);
        }
      //  cout << x << " " << ans << endl;
    }
    double res = ans*1.0 / ((LL)n * (n-1));
    printf("%.6f\n",res);
}
int main()
{
    cin >> n >> m;
    for(int i = 0 ; i < n ;i ++)cin >> v[i];
    cout << endl;
    for(int i = 0 ; i < m ;i ++){
        int x,y;
        cin >> x >> y;
        x--; y--;
        g[x].insert(MP(v[y],y));
        g[y].insert(MP(v[x],x));
    }
    for(int i = 0 ; i< n ;i ++)sv[i] = MP(v[i],i);
    sort(sv,sv+n);
    solve();
    return 0;
}
