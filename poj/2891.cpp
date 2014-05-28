#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<algorithm>
#include<map>
#include<climits>
using namespace std;

typedef pair<int,int> PII;
typedef long long LL;
#define CLR(a,b) memset(a,b,sizoef(a))
#define MP make_pair
#define INF 0x7f7f7f7f
#define ABS(x) ((x)>0?(x):(-x))

int n;
LL x,y;
LL gcd(LL a,LL b)
{
    if(b==0)return a;
    return gcd(b, a%b);
}
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

int main()
{
    while(~scanf("%d",&n)){
        LL a1,b1;
        bool ok = 1;
        scanf("%lld%lld",&a1,&b1);
        if(n == 1){
            printf("%lld\n",a1+b1);
            continue;
        }
        for(int i = 1 ; i < n ; i ++){
            LL a2,b2;
            scanf("%lld%lld",&a2,&b2);
            LL gg = gcd(a1,a2);
            if((b2 - b1) % gg != 0)ok = 0;
            if(ok){
                LL d = ext_gcd(a1,a2);
                x *= (b2-b1)/d;
                x = x - (x*d/a2)*(a2/d);
                if(x < 0) x += a2/d;
                LL c = a1 * x + b1;
                b1 = c;
                a1 = a1 / d * a2;
            }
        }
        if(ok)
            printf("%lld\n",b1);
        else
            printf("-1\n");
    }
    return 0;
}