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
const int N = 1000+20;

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
vector<int> g[N];
LL v[N];

int main()
{
    cin >> n >> m;
    for(int i = 0 ; i < n ;i ++)cin >> v[i];
    LL ans = 0;
    for(int i = 0 ; i < m ;i ++){
        int x,y;
        cin >> x >> y;
        x--; y--;
        ans += min(v[x],v[y]);
    }
    cout << ans << endl;
    return 0;
}
