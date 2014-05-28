#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<algorithm>
#include<map>
using namespace std;

typedef pair<int,int> PII;
typedef long long LL;
#define CLR(a,b) memset(a,b,sizoef(a))
#define MP make_pair
#define INF 0x7f7f7f7f
#define ABS(x) ((x)>0?(x):(-x))

int a,b,c,d,x,y;
LL ext_gcd(LL a,LL b)
{
    LL t,d;
    if(b==0)
    {
        x=1;
        y=0;
        return a;
    }
    d=ext_gcd(b,a%b);
    t=x;
    x=y;
    y=t-(a/b)*y;
    return d;
}

bool change;
void solve()
{
    int k = y * d / a;
    PII ans = MP(INF,INF);
    int use = 0;
    for(int i = k-2 ; i <= k+2; i ++){
        int nx = x + b/d * i;
        int ny = y - a/d * i;
        int tot1 = ABS(nx) + ABS(ny);
        int tot2 = a * ABS(nx) + b * ABS(ny);
        if(MP(tot1,tot2) < ans){
            ans = MP(tot1, tot2);
            use = i;
        }
    }
    x = x + b/d * use;
    y = y - a/d * use;
    if(change)swap(x,y);
    //  printf("%d %d\n",x,y);
    if(x < 0)x = -x;
    if(y < 0)y = -y;
    printf("%d %d\n",x,y);
}
int main()
{
    while(~scanf("%d%d%d",&a,&b,&c)){
        if(a == 0 && b == 0 && c == 0)break;
        if(a < b){
            swap(a,b);
            change = 1;
        }else{
            change = 0;
        }
        d = ext_gcd(a,b);
        x *= c/d;
        y *= c/d;
        solve();
    }
    return 0;
}