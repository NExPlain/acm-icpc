#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<algorithm>
#include<cstring>
using namespace std;

typedef long long LL;
const int MOD = 1000000000+7;
int n,k,d;
LL dp[101][101];

LL get_dp(int sum, int maximum)
{
    if(sum == maximum)return dp[sum][maximum] = 1;
    if(sum < maximum) return dp[sum][maximum] = 0;
    if(dp[sum][maximum] != -1) return dp[sum][maximum];
    int ret = 0;
    for(int i = 1; i <= k ; i ++){
        if(sum >= i){
            if(i < maximum)
                ret = (ret + get_dp(sum - i, maximum)) % MOD;
            else if(i == maximum){
                for(int j = 1 ; j <= maximum; j ++){
                    ret = (ret + get_dp(sum - i , j)) % MOD;
                }
            }
        }
    }
    return dp[sum][maximum] = ret % MOD;
}
int main()
{
    cin >> n >> k >> d;
    memset(dp, -1, sizeof(dp));
    int ans = 0;
    for(int i = d ; i <= k ; i++){
        ans = (ans + get_dp(n, i)) % MOD;
    }
    printf("%d\n",ans);
    return 0;
}
