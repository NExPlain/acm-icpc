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

const int N = 200+20;
int n,K;
int a[N];
void solve()
{
    int ans = -INF;
    priority_queue<int> big;
    priority_queue<int,vector<int>,greater<int> > smal;
    for(int l = 0 ; l < n ; l ++){
        for(int r = l ; r < n ; r ++){
            while(smal.size())smal.pop();
            while(big.size())big.pop();
            int sum = 0;
            for(int k = 0 ; k < n ; k ++){
                if(l <= k && k <= r){
                    smal.push(a[k]);
                    sum += a[k];
                }else{
                    big.push(a[k]);
                }
            }
            for(int k = 0 ; k < K ; k++){
                if(smal.size() && big.size()){
                    int x = big.top();
                    big.pop();
                    int y = smal.top();
                    smal.pop();
                    if(x > y){
                        sum += x - y;
                    }
                }
            }
            ans = max(ans, sum);
        }
    }
    printf("%d\n",ans);
}
int main()
{
    scanf("%d%d",&n,&K);
    for(int i = 0; i < n ; i++){
        scanf("%d",&a[i]);
    }
    solve();
    return 0;
}
